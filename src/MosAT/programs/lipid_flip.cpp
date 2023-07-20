#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h> 
#include <mpi.h>
#include <fstream>
#include <vector>
using namespace std;

#include "xdr/include/xdrfile_xtc.h"                         //used to read xtc files 
#include "xdr/include/xdr_seek.h"                            //used to get and set the file position in xtc and trr files
#include "xdr/include/xdrfile_trr.h"                         //used to read trr files
#include "xdr/include/xdrfile.h"                             //used to read xtc and trr files
#include "xdr/include/trr_header.h"                          //used to read the header info of trr files
#include "headers/multi_dim_vec.h"                           //This defines multidimensional vectors
#include "headers/switch.h"                                  //This defines a switch (on, off)
#include "headers/file_reader.h"                             //This has basic routines for reading text files
#include "headers/vector_mpi.h"                              //This has routines for collecting vector data
#include "headers/mosat_routines.h"                         //This is where most of the functions called in main are located
#include "headers/file_naming.h"                             //This has routines for added tags to an existing file name    
#include "headers/file_naming_mpi.h"                         //This has routines for added tags to an existing file name (mpi)
#include "headers/command_line_args_mpi.h"                   //This has routines for adding command line arguments
#include "MosAT/program_variables/pv_lipid_flip.h"          //This has the variables specific to the analysis program
#include "headers/array.h"                                   //This has routines used for working with arrays
#include "headers/performance.h"                             //This has a class for logging performance data
#include "headers/index.h"                                   //This has a class for working with index files
#include "headers/traj.h"                                    //This has a class for working with the trajectory
#include "headers/leaflet_finder.h"                          //This has routines used to find leaflets in membrane simulations
#include "headers/protein_finder.h"                          //This has routines used to find protein atoms
#include "headers/sol_finder.h"                              //This has routines used to find the solvent
#include "headers/grid.h"                                    //This has routines used for working with a grid
#include "headers/protein.h"                                 //This has routines used for working with protein data
#include "headers/force_serial.h"                            //This has routines used for forcing the code to run on a single mpi process
#include "headers/param.h"                                   //This has routines used for reading complex parameter data

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function checks for lipid flips                                                                      //
//                                                                                                           //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void flip(Trajectory &traj,system_variables &s,program_variables &p,Param &param,dv2d &zcoord,iv1d &resid)
{
    int i           = 0;                      //standard variable used in loops
    int j           = 0;                      //standard variable used in loops

    //allocate memory to hold current frame lipids and resid
    dv1d current_frame_z(0,0);
    iv1d current_resid(0,0);

    //add the time as first entry
    current_frame_z.push_back((double)traj.get_frame_global());

    for(i=0; i<traj.target_leaflet.size(); i++) //loop over target leaflet atoms
    {
        //jump to the next lipid
        i = traj.next_target_lipid(i);

        //get the first and last atom of the current lipid
        int min = traj.t_lip_start(i);
        int max = traj.t_lip_end(i);

        for(j=0; j<param.main_size_y(); j++) //loop over lipid types 
        {
            if(strcmp(traj.res_name[min].c_str(), param.param_main_s[j][0].c_str()) == 0) //lipid type is correct
            {
                //extract the atoms used for the center
                sv1d target_atoms = param.get_column_sec_s(j,0);

                //compute the center for the lipid
                dv1d r_center = traj.center(target_atoms,min,max);
   
                //add the current lipid z coord.  
                current_frame_z.push_back(r_center[2]);

                //store resid
                current_resid.push_back(traj.res_nr[min]);
            }
        }
    }
    //add the current frame z-coordinates
    zcoord.push_back(current_frame_z);

    if(s.world_rank == 0 && traj.current_frame == 0)
    {
        resid = current_resid; 
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Collect z-coords and write out data.                                                                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,dv2d &zcoord,iv1d &resid)
{
    int i = 0;
    int j = 0;

    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();
    
    if(s.world_rank == 0)
    {
        printf("\nfinalizing analysis. This requires communication of the z-coordinates and could take some time. \n");
    }

    //here we collect the zcoord measurements from each core
    collect_dv2d(s.world_size,s.world_rank,zcoord);

    //write the z-coord data to output file
    if(s.world_rank == 0)
    {
        FILE *flip_file = fopen(p.flip_file_name.c_str(), "w");
        if(flip_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",p.flip_file_name.c_str());
        }
        else
        {
            fprintf(flip_file,"# %10s %5s: %s\n","column","1","trajectory frame");
            fprintf(flip_file,"# %10s %5d-%5d: %s\n","column",2,zcoord[0].size(),"lipid z-coordinates");
            for(i=0; i<zcoord.size(); i++) 
            {
                for(j=0; j<zcoord[i].size(); j++)
                {
                    fprintf(flip_file," %10f ",zcoord[i][j]);
                }
                fprintf(flip_file,"\n");
            }
            fclose(flip_file);
        }

        string resid_file_name = add_tag(p.flip_file_name,"_resid");
        FILE *resid_file = fopen(resid_file_name.c_str(), "w");
        if(resid_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",resid_file_name.c_str());
        }
        else
        {
            fprintf(resid_file,"# %15s %10s \n","column number","res_id");
            for(i=0; i<resid.size(); i++)
            {
                fprintf(resid_file,"  %15d %10d \n",i+1,resid[i]);
            }
            fclose(resid_file);
        }
    }

    //compute and return time spent in function
    return (clock() - s.t)/CLOCKS_PER_SEC;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This is the main function which executes the other functions.                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, const char * argv[]) 
{
    //Here we set up the MPI environment
    MPI_Init(NULL, NULL);
    
    //nowe we define the system and program variables and initialize them
    system_variables s;
    program_variables p;
    initialize_system_variables(&s);
    initialize_program_variables(&p);

    //create object for logging performance data
    Performance perf; 

    //set the name of your analysis program here
    s.program_name = "Lipid Flip";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                          s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                          s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                         s.world_rank, s.cl_tags, &p.b_print,   0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                           s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                                      s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                         s.world_rank, s.cl_tags, &p.b_lsq,     0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                           s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",           s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-crd",    p.param_file_name,            "Selection card with target lipids (crd)",                             s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-flip",   p.flip_file_name,             "Output file with z-coordinates data for each trajectory frame (dat)", s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,           "PDB file with sorted leaflets (pdb)",                                 s.world_rank, s.cl_tags, &p.b_lf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,  "File with additional leaflet finder parameters (prm)",                s.world_rank, s.cl_tags, &p.b_lf_param,0);
    add_argument_mpi_i(argc,argv,"-leaf",   &p.leaflet,                   "Which leaflet? (0:both 1:upper 2:lower)",                             s.world_rank, s.cl_tags, nullptr,      1);
    conclude_input_arguments_mpi(argc,argv,s.world_rank,s.program_name,s.cl_tags);

    //create a trajectory
    Trajectory traj; 

    //set trajectory parameters
    traj.set_block_parallel(on);
    traj.set_traj(p.in_file_name);
    traj.set_ref(p.ref_file_name);
    traj.set_traj_w(p.out_file_name,p.b_print);
    traj.set_lsq(p.lsq_index_file_name,p.b_lsq,p.lsq_dim,p.lsq_ref);
    traj.set_res(p.stride,p.start_frame,p.end_frame);

    //analyze the trajectory (log time spent) 
    perf.log_time(traj.build(),"Analyze Trajectory");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //check file extensions                                                                                     
    check_extension_mpi(s.world_rank,"-crd",p.param_file_name,".crd");
    check_extension_mpi(s.world_rank,"-flip",p.flip_file_name,".dat");

    if(p.b_lf_pdb == 1)
    {
        check_extension_mpi(s.world_rank,"-lf_pdb",p.lf_pdb_file_name,".pdb");
    }
    if(p.b_lf_param == 1)
    {
        check_extension_mpi(s.world_rank,"-lf_prm",p.leaflet_finder_param_name,".prm");
    }

    //create parameter files
    Param param;

    //read parameter files
    param.get_param(p.param_file_name,2,1,1);

    //check the integrity of the parameter files
    if(param.check_file() == 0) //bad files. kill program 
    {
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }

    //run leaflet finder
    traj.get_leaflets(p.leaflet,p.leaflet_finder_param_name,p.b_lf_param);

    //print a pdb with distinguished leaflets
    traj.write_leaflets(p.lf_pdb_file_name,p.b_lf_pdb);

    //print info about the target leaflet
    traj.get_leaflet_stats();

    //print info about the lipid selection
    traj.get_lipid_selection_stats(param.get_column_s(0),"-crd");

    //count lipid in the selected leaflet
    p.num_lip_t = traj.count_target_lipids_type(param.get_column_s(0));

    //allocate memory to hold z-coords
    dv2d zcoord(0,dv1d(p.num_lip_t+1,0.0));

    //allocate memory to hold resids
    iv1d resid(0,0);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //print info about the worlk load distribution
    traj.workload();

    //print that analysis is beginning
    traj.report_progress();

    s.t = clock();
    //read read frames of the trajectory and perform analysis
    for(traj.current_frame=0; traj.current_frame<traj.get_num_frames(); traj.current_frame++)
    {
        traj.read_traj_frame();

        traj.do_fit();

        flip(traj,s,p,param,zcoord,resid);

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect the contacts from each mpi process and compute the average
    perf.log_time(finalize_analysis(traj,s,p,zcoord,resid),"Fin Ana");

    //splice temporary traj file together (log time spent)
    perf.log_time(traj.finalize_trajectory(),"Finalize Trajectory");

    //print the performance stats
    perf.print_stats();

    //print closing statements
    print_closing(s.world_rank);

    //relinquish the mpi environment
    MPI_Finalize();

    return 0;
}

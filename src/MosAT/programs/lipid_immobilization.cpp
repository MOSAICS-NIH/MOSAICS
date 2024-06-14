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

#include "xdr/include/xdrfile_xtc.h"                               //used to read xtc files 
#include "xdr/include/xdr_seek.h"                                  //used to get and set the file position in xtc and trr files
#include "xdr/include/xdrfile_trr.h"                               //used to read trr files
#include "xdr/include/xdrfile.h"                                   //used to read xtc and trr files
#include "xdr/include/trr_header.h"                                //used to read the header info of trr files
#include "headers/multi_dim_vec.h"                                 //This defines multidimensional vectors
#include "headers/switch.h"                                        //This defines a switch (on, off)
#include "headers/file_reader.h"                                   //This has basic routines for reading text files
#include "headers/vector_mpi.h"                                    //This has routines for collecting vector data
#include "headers/mosat_routines.h"                                //This is where most of the functions called in main are located
#include "headers/file_naming.h"                                   //This has routines for added tags to an existing file name    
#include "headers/file_naming_mpi.h"                               //This has routines for added tags to an existing file name (mpi)
#include "headers/command_line_args_mpi.h"                         //This has routines for adding command line arguments
#include "MosAT/program_variables/pv_lipid_immobilization.h"       //This has the variables specific to the analysis program
#include "headers/array.h"                                         //This has routines used for working with arrays
#include "headers/performance.h"                                   //This has a class for logging performance data
#include "headers/index.h"                                         //This has a class for working with index files
#include "headers/traj.h"                                          //This has a class for working with the trajectory
#include "headers/leaflet_finder.h"                                //This has routines used to find leaflets in membrane simulations
#include "headers/protein_finder.h"                                //This has routines used to find protein atoms
#include "headers/sol_finder.h"                                    //This has routines used to find the solvent
#include "headers/grid.h"                                          //This has routines used for working with a grid
#include "headers/protein.h"                                       //This has routines used for working with protein data
#include "headers/force_serial.h"                                  //This has routines used for forcing the code to run on a single mpi process
#include "headers/param.h"                                         //This has routines used for reading complex parameter data
#include "headers/atom_select.h"                                   //This has routines used for making atom selections using a selection text

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes the resid for target lipids to file                                                  //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void lookup_lipids(Trajectory &traj,system_variables &s,program_variables &p,Index &crd,string lookup_file_name)
{
    int    i        = 0;                                //standard variable used in loops
    int    j        = 0;                                //standard variable used in loops
    int    count    = 0;                                //count the target lipids as they are encountered

    if(s.world_rank == 0 && traj.current_frame == 0)
    {
        FILE *lookup_file = fopen(lookup_file_name.c_str(), "w");
        if(lookup_file == NULL)
        {
            printf("failure opening %s for writing. \n",lookup_file_name.c_str());
            fflush(stdout);
        }
        else
        {
            //print header info
            fprintf(lookup_file," %10s %10s \n","lip_nr","res_id");
            fprintf(lookup_file,"%22s \n","______________________");

            for(i=0; i<traj.target_leaflet.size(); i++) //loop over target leaflet
            {
                //jump to the next lipid
                i = traj.next_target_lipid(i);

                int min = traj.t_lip_start(i);         //first atom of the current lipid
                int max = traj.t_lip_end(i);           //last atom of the current lipid

                for(j=0; j<crd.index_s.size(); j++) //loop over lipid types
                {
                    if(strcmp(traj.res_name[min].c_str(), crd.index_s[j].c_str()) == 0) //lipid type is correct
                    {
                        fprintf(lookup_file," %10d %10d \n",count,traj.res_nr[min]);
                        count++;
                    }
                }
            }
           
            //close file
            fclose(lookup_file);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function stamps each lipid's position onto a grid.                                                   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void track_lipids(Trajectory &traj,system_variables &s,program_variables &p,Index &crd,vector <Grid_i> &lip_pos)
{
    int    i        = 0;                                //standard variable used in loops
    int    j        = 0;                                //standard variable used in loops
    int    k        = 0;                                //standard variable used in loops
    int    count    = 0;                                //count the target lipids as they are encountered

    for(i=0; i<traj.target_leaflet.size(); i++) //loop over target leaflet
    {
        //jump to the next lipid
        i = traj.next_target_lipid(i);

        int min = traj.t_lip_start(i);         //first atom of the current lipid
        int max = traj.t_lip_end(i);           //last atom of the current lipid

        for(j=0; j<crd.index_s.size(); j++) //loop over lipid types
        {
            if(strcmp(traj.res_name[min].c_str(), crd.index_s[j].c_str()) == 0) //lipid type is correct
            {
                sv1d target_atoms(0);

                for(k=min; k<=max; k++) //loop over current residue 
                {
                    target_atoms.push_back(traj.atom_name[k]);
                }

                dv1d center = traj.center(target_atoms,min,max);

                lip_pos[count].stamp(center[0],center[1],p.radius,1);

                count++;
	    }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Collect grids for each lipid and write data to file                                                       //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,vector <Grid_i> &lip_pos)
{
    int    i        = 0;                                //standard variable used in loops
    int    j        = 0;                                //standard variable used in loops
    int    k        = 0;                                //standard variable used in loops

    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nfinalizing analysis. This requires communication of the grid and could take some time. \n");
    }

    for(i=0; i<lip_pos.size(); i++) //loop over lipids
    {	
        //collect lip_pos from all ranks    
        lip_pos[i].collect_grid();

        //write lip_pos to file
        if(s.world_rank == 0)
        {
            //create grid to handle insignificant data
            Grid_i nan;
            nan.set_dim(p.APS,p.num_g_x,p.num_g_y);

            //set nan
            for(j=0; j<lip_pos[i].num_x(); j++) //loop over x
            {
                for(k=0; k<lip_pos[i].num_y(); k++) //loop over y
                {
            	if(lip_pos[i].grid[k][j] > 0) //some samples were collected
                    {  
                        nan.grid[k][j] = 0;
                    }
                    else //no samples were collected
                    {
                        nan.grid[k][j] = 1;
                    }
                }
            }
            lip_pos[i].write_grid(nan);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

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
    s.program_name = "Lipid Immobilization";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                  s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (tpr, gro)",                                  s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                 s.world_rank, s.cl_tags, &p.b_print,   0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                   s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                              s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                               s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting",                                       s.world_rank, s.cl_tags, &p.b_lsq,     0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                   s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",   s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lim",    p.lim_file_name,              "Base file name for lipid immobilization data (dat)",          s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-crd",    p.crd_file_name,              "Name of the selection card with lipid types (crd)",           s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,           "PDB file with sorted leaflets ",                              s.world_rank, s.cl_tags, &p.b_lf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,  "File with additional leaflet finder parameters ",             s.world_rank, s.cl_tags, &p.b_lf_param,0);
    add_argument_mpi_d(argc,argv,"-APS",    &p.APS,                       "Area per grid box (nm^2)",                                    s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-r",      &p.radius,                    "Radius of the target atom (nm)",                              s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-bx",     &p.box_x,                     "Grid x dimension (nm)",                                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-by",     &p.box_y,                     "Grid y dimension (nm)",                                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-leaf",   &p.leaflet,                   "Which leaflet? (0:both 1:upper 2:lower)",                     s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-odf",    &p.odf,                       "What is the output data format? (0:matrix 1:vector)",         s.world_rank, s.cl_tags, nullptr,      0);
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
    check_extension_mpi(s.world_rank,"-crd",p.crd_file_name,".crd");
    check_extension_mpi(s.world_rank,"-lim",p.lim_file_name,".dat");
    if(p.b_lf_pdb == 1)
    {
        check_extension_mpi(s.world_rank,"-lf_pdb",p.lf_pdb_file_name,".pdb");
    }
    if(p.b_lf_param == 1)
    {
        check_extension_mpi(s.world_rank,"-lf_prm",p.leaflet_finder_param_name,".prm");
    }

    //create index objects
    Index crd;

    //read the index files
    crd.get_index(p.crd_file_name);

    //run leaflet finder
    traj.get_leaflets(p.leaflet,p.leaflet_finder_param_name,p.b_lf_param);

    //print a pdb with distinguished leaflets
    traj.write_leaflets(p.lf_pdb_file_name,p.b_lf_pdb);

    //print info about the target leaflet
    traj.get_leaflet_stats();

    //print info about the lipid selection
    traj.get_lipid_selection_stats(crd.get_column_s(1,0),"-crd");

    //set the grid cell size 
    p.cell_size = sqrt(p.APS);

    //get the grid dimensions
    get_grid_size(p.box_x,p.box_y,traj.ibox,&p.num_g_x,&p.num_g_y,p.cell_size);

    //count target lipids 
    int num_target_lipids = traj.count_target_lipids_type(crd.get_column_s(1,0));

    //make a grid for each lipid to track its position
    vector <Grid_i>  lip_pos(num_target_lipids);

    double mem = (double)p.num_g_x*(double)p.num_g_y*4.0*(double)num_target_lipids/1000000.0;
    if(s.world_rank == 0)
    {
        printf("Estimated memory requirements: %f MB. \n\n",mem);
    }

    int i = 0;
    for(i=0; i<lip_pos.size(); i++) //loop over lipids
    {
        //set the grid dimensions
        lip_pos[i].set_dim(p.APS,p.num_g_x,p.num_g_y);

        //create a filename for lip pos grids
        string this_file_name = add_tag(p.lim_file_name,"_" + to_string(i));

        //set the output file name for grid
        lip_pos[i].set_output(this_file_name,p.odf);
    }

    //print info about the grid
    print_grid_stats(p.box_x,p.box_y,traj.ibox,p.num_g_x,p.num_g_y,p.cell_size,s.world_rank);

    //create file name for lipid lookup data
    string lookup_file_name = add_tag(p.lim_file_name,"_lookup");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //print info about the worlk load distribution
    traj.workload();

    //print that analysis is beginning
    traj.report_progress();

    s.t = clock();
    //read read frames of the trajector and perform analysis
    for(traj.current_frame=0; traj.current_frame<traj.get_num_frames(); traj.current_frame++)
    {
        traj.read_traj_frame();

        traj.do_fit();

        lookup_lipids(traj,s,p,crd,lookup_file_name);

        track_lipids(traj,s,p,crd,lip_pos);

        traj.set_beta_lf();

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect dist from mpi processes and compute average
    perf.log_time(finalize_analysis(traj,s,p,lip_pos),"Fin Ana");

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

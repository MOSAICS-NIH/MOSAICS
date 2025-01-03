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
#include "MosAT/program_variables/pv_lipid_msd.h"           //This has the variables specific to the analysis program
#include "headers/array.h"                                   //This has routines used for working with arrays
#include "headers/performance.h"                             //This has a class for logging performance data
#include "headers/index.h"                                   //This has a class for working with index files
#include "headers/traj.h"                                    //This has a class for working with the trajectory
#include "headers/leaflet_finder.h"                          //This has routines used to find leaflets in membrane simulations
#include "headers/protein_finder.h"                          //This has routines used to find protein atoms
#include "headers/sol_finder.h"                              //This has routines used to find the solvent
#include "headers/grid.h"                                    //This has routines used for working with a grid
#include "headers/protein.h"                                 //This has routines used for working with protein data
#include "headers/parallel.h"                                //This has routines for different parallelization schemes
#include "headers/force_serial.h"                            //This has routines used for forcing the code to run on a single mpi process
#include "headers/param.h"                                   //This has routines used for reading complex parameter data

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the center for each lipid and stores their coordinates                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void store_coords(Trajectory &traj,system_variables &s,program_variables &p,Param &param,dv3d &my_coords,iv1d &resid)
{
    int    i           = 0;                          //standard variable used in loops
    int    j           = 0;                          //standard variable used in loops
    int    prev_lip    = -1;                         //last lipid encountered
    int    lip_count   = -1;                         //count lipids as they are encountered 
    int    ef_lip      = -1;                         //the lipid index for the mpi process

    for(i=0; i<traj.target_leaflet.size(); i++) //loop over target leaflet 
    {
        //get the first and last atom of the current lipid
        int min = traj.t_lip_start(i);
        int max = traj.t_lip_end(i);

        //jump to the next lipid
        i = traj.next_target_lipid(i);

        for(j=0; j<param.main_size_y(); j++) //loop over lipid types 
        {
            if(strcmp(traj.res_name[min].c_str(), param.param_main_s[j][0].c_str()) == 0) //lipid type is correct
            {
                //check for a new lipid
                if(traj.res_nr[min] != prev_lip)
                {
                   lip_count++;
                   prev_lip = traj.res_nr[min];
                }

                if(lip_count >= traj.lipid_start && lip_count <= traj.lipid_end) //each core does their contributing lipids
                {
                    ef_lip++;

                    //extract atoms involved in geometric center 
                    sv1d target_atoms = param.get_column_sec_s(j,0);

                    //compute the center for the lipid
                    dv1d r_center = traj.center(target_atoms,min,max);

                    my_coords[ef_lip][traj.current_frame][0] = r_center[0]; 
                    my_coords[ef_lip][traj.current_frame][1] = r_center[1];

                    if(traj.current_frame == 0)
                    {
                        resid.push_back(traj.res_nr[min]);
                    }

                    //free memory
                    vector<string>().swap(target_atoms);
                    vector<double>().swap(r_center);
		}
            }
        }
    }
    //printf("step %d world rank %d resid.capacity %d my_coords.capacity() %d \n",traj.current_frame,s.world_rank,resid.capacity(),my_coords.capacity());
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function analyzes the coordinates and computes the msd.                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analyisis(Trajectory &traj,system_variables &s,program_variables &p,dv3d &my_coords,iv1d &resid)
{
    int    i  = 0;                          //standard variable used in loops
    int    j  = 0;                          //standard variable used in loops
    int    k  = 0;                          //standard variable used in loops
    double x  = 0.0;                        //the current x value
    double y  = 0.0;                        //the current y value
    double x0 = 0.0;                        //the initial x value
    double y0 = 0.0;                        //the initial y value

    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock(); 

    dv2d msd(traj.get_num_frames()-1,dv1d(traj.my_lipids,0.0));    //holds the msd data for each cores lipids

    for(i=0; i<traj.my_lipids; i++) //loop over lipids
    {
        if(s.world_rank == 0)
        {
            printf("working on lipid %d \n",i);
        }

        for(j=1; j<traj.get_num_frames(); j++) //loop over deltas
        {
            double this_msd = 0.0;

            for(k=0; k<traj.get_num_frames()-j; k++) //slide current delta
            {
                x0 = my_coords[i][k][0];
                y0 = my_coords[i][k][1];
                x  = my_coords[i][k+j][0];
                y  = my_coords[i][k+j][1];

                this_msd = this_msd +      (x-x0)*(x-x0) + (y-y0)*(y-y0); 
            }
            this_msd = this_msd/(double)(traj.get_num_frames()-j);

            msd[j-1][i] = this_msd;
        }
    }

    //collect the msd from mpi processes
    FILE *msd_file;
    if(s.world_rank == 0)
    {
        msd_file = fopen(p.msd_file_name.c_str(), "w");
    }

    //print the header information
    collect_iv1d(s.world_size,s.world_rank,resid);

    if(s.world_rank == 0)
    {
        fprintf(msd_file," #Column 1: %10s \n","time (ps)");
        fprintf(msd_file," #Column 2-%d: %10s \n",resid.size()+1,"MSD for the individual lipids (nm^2)");
        fprintf(msd_file," #Column %d: %10s \n",resid.size()+2,"MSD averaged over the lipids (nm^2)");
        fprintf(msd_file," %10s ","#time"); 
        for(j=0; j<p.num_lipids; j++) //loop over lipids
        {
            fprintf(msd_file," %10d ",resid[j]);
        }
        fprintf(msd_file," %10s \n","Average");
        fflush (msd_file);
    }


    for(i=0; i<traj.get_num_frames()-1; i++) //loop over deltas
    {
        dv1d msd_lipids = collect_and_clone_dv1d(s.world_size,s.world_rank,msd[i]);

        int count      = 0;
        double avg_msd = 0;

        if(s.world_rank == 0)
        {
            fprintf(msd_file," %10f ",(double)(i+1)*p.ef_dt);

            for(j=0; j<p.num_lipids; j++) //loop over lipids
            {
                avg_msd = avg_msd + msd_lipids[j];
                count++;
                fprintf(msd_file," %10f ",msd_lipids[j]);
            }   
            avg_msd = avg_msd/(double)count;

            fprintf(msd_file," %10f \n",avg_msd);
            fflush (msd_file);
        }

        //free memory
        vector<double>().swap(msd_lipids);
    }

    if(s.world_rank == 0)
    {
        fclose(msd_file);
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
    s.program_name = "lipid MSD";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                  s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                  s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                 s.world_rank, s.cl_tags, &p.b_print,   0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                   s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                              s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                               s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                 s.world_rank, s.cl_tags, &p.b_lsq,     0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                   s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",   s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-crd",    p.param_file_name,            "Selection card (crd)",                                        s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-msd",    p.msd_file_name,              "Output data file with MSD data (dat)",                        s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-dt",     &p.ef_dt,                     "Time between analyzed frames (accounting for stride, ps)",    s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,           "PDB file with sorted leaflets (pdb)",                         s.world_rank, s.cl_tags, &p.b_lf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,  "File with additional leaflet finder parameters (prm) ",       s.world_rank, s.cl_tags, &p.b_lf_param,0);
    conclude_input_arguments_mpi(argc,argv,s.world_rank,s.program_name,s.cl_tags);

    //create a trajectory
    Trajectory traj; 

    //set trajectory parameters
    traj.set_block_parallel(off);
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
    check_extension_mpi(s.world_rank,"-msd",p.msd_file_name,".dat");
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

    //count the number of lipids in the target leaflet(s) of the target type
    p.num_lipids = traj.count_target_lipids_type(param.get_column_s(0));

    //parallelize by lipids
    traj.parallelize_by_lipid(p.num_lipids);

    //create vector to hold lipid coordinates
    dv3d my_coords(traj.my_lipids,dv2d(traj.get_num_frames(),dv1d(2,0.0)));

    //create vector to hold resid for lipids
    iv1d resid(0,0);

    //print memory estimate
    double mem_1 = (double)traj.my_lipids*(double)traj.get_num_frames()*2.0*8.0/1000000.0; //my_coords
    double mem_2 = (double)(traj.get_num_frames()-1)*(double)traj.my_lipids*8.0/1000000.0; //msd
    double mem_3 = (double)(traj.get_num_frames()-1)*(double)p.num_lipids*8.0/1000000.0;   //msd_lipids
    double mem_4 = (double)traj.my_lipids*4.0/1000000.0;                                   //resid
    double mem   = mem_1 + mem_2 + mem_3 + mem_4;

    if(s.world_rank == 0)
    {
        printf("Estimated memory: (MB). \n");
	printf("%s\n","------------------------");
        printf("  %10s: %10.1f \n","my_coords",mem_1);
        printf("  %10s: %10.1f \n","msd",mem_2);
        printf("  %10s: %10.1f \n","msd_lipids",mem_3);
        printf("  %10s: %10.1f \n","resid",mem_4);
        printf("  %10s: %10.1f \n\n","total",mem);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //print info about the worlk load distribution
    traj.workload_lipid();

    //print that analysis is beginning
    traj.report_progress();

    s.t = clock();
    //read read frames of the trajectory and perform analysis
    for(traj.current_frame=0; traj.current_frame<traj.get_num_frames(); traj.current_frame++)
    {
        traj.read_traj_frame();

        traj.do_fit();

        store_coords(traj,s,p,param,my_coords,resid);

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //compute msd values
    perf.log_time(finalize_analyisis(traj,s,p,my_coords,resid),"Fin Ana");

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

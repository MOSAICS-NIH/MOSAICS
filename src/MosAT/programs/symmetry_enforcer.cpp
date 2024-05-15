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

#include "xdr/include/xdrfile_xtc.h"                       //used to read xtc files 
#include "xdr/include/xdr_seek.h"                          //used to get and set the file position in xtc and trr files
#include "xdr/include/xdrfile_trr.h"                       //used to read trr files
#include "xdr/include/xdrfile.h"                           //used to read xtc and trr files
#include "xdr/include/trr_header.h"                        //used to read the header info of trr files
#include "headers/multi_dim_vec.h"                         //This defines multidimensional vectors
#include "headers/switch.h"                                //This defines a switch (on, off)
#include "headers/file_reader.h"                           //This has basic routines for reading text files
#include "headers/vector_mpi.h"                            //This has routines for collecting vector data
#include "headers/mosat_routines.h"                        //This is where most of the functions called in main are located
#include "headers/file_naming.h"                           //This has routines for added tags to an existing file name    
#include "headers/file_naming_mpi.h"                       //This has routines for added tags to an existing file name (mpi)
#include "headers/command_line_args_mpi.h"                 //This has routines for adding command line arguments
#include "MosAT/program_variables/pv_symmetry_enforcer.h"  //This has the variables specific to the analysis program
#include "headers/array.h"                                 //This has routines used for working with arrays
#include "headers/performance.h"                           //This has a class for logging performance data
#include "headers/index.h"                                 //This has a class for working with index files
#include "headers/traj.h"                                  //This has a class for working with the trajectory
#include "headers/leaflet_finder.h"                        //This has routines used to find leaflets in membrane simulations
#include "headers/protein_finder.h"                        //This has routines used to find protein atoms
#include "headers/sol_finder.h"                            //This has routines used to find the solvent
#include "headers/grid.h"                                  //This has routines used for working with a grid
#include "headers/protein.h"                               //This has routines used for working with protein data
#include "headers/force_serial.h"                          //This has routines used for forcing the code to run on a single mpi process
#include "headers/param.h"                                 //This has routines used for reading complex parameter data

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function swaps coordinates for each cycle                                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void swap_coords(Trajectory &traj,system_variables &s,program_variables &p,Param &param,int cycle)
{
    int i = 0;        //standard variable used in loops
    int j = 0;        //standard variable used in loops
    int k = 0;        //standard variable used in loops

    rvec *x;                                        //holds a copy of the atomic coordinates
    x = (rvec *)calloc(traj.atoms() , sizeof(*x));  //allocate memory to hold coords

    //make a copy of the atomic coordinates
    for(i=0; i<traj.atoms(); i++) //loop over atoms
    {
        x[i][0] = traj.r[i][0];
        x[i][1] = traj.r[i][1];
        x[i][2] = traj.r[i][2];
    }

    for(i=0; i<param.main_size_y(); i++) //loop over the protomers  
    {
        int partner = (i + cycle)%param.main_size_y();    //the exchange partner

        for(j=0; j<param.sec_size_y(i); j++) //loop over current protomer atoms
        {
            traj.r[param.param_sec_i[i][j][0]-1][0] = x[param.param_sec_i[partner][j][0]-1][0];
            traj.r[param.param_sec_i[i][j][0]-1][1] = x[param.param_sec_i[partner][j][0]-1][1];
            traj.r[param.param_sec_i[i][j][0]-1][2] = x[param.param_sec_i[partner][j][0]-1][2];
        }
    }

    free(x);
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
    s.program_name = "Symmetry Enforcer";

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
    add_argument_mpi_s(argc,argv,"-crd",    p.param_file_name,            "Selection card with indices for each protomer (crd)",         s.world_rank, s.cl_tags, nullptr,      1);
    conclude_input_arguments_mpi(argc,argv,s.world_rank,s.program_name,s.cl_tags);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //check file extensions                                                                                     
    check_extension_mpi(s.world_rank,"-crd",p.param_file_name,".crd");

    //create parameter files
    Param param;

    //read parameter files
    param.get_param(p.param_file_name,1,0,1);

    //check the integrity of the parameter files
    if(param.check_file() == 0) //bad files. kill program 
    {
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }

    //check that protomers are the same size
    for(p.cycle=0; p.cycle<param.main_size_y(); p.cycle++) //loop over the protomers  
    {
        if(param.sec_size_y(p.cycle) != param.sec_size_y(0))
        {
            if(s.world_rank == 0)
            {
                printf("Protomer size difference found. Protomer %d size %d. Protomer %d size %d. \n",p.cycle,param.sec_size_y(p.cycle),0,param.sec_size_y(0));
            }

	    MPI_Finalize();
            exit(EXIT_SUCCESS);
	}
    }   
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //loop over cycles
    for(p.cycle=0; p.cycle<param.main_size_y(); p.cycle++) //loop over protomers
    {
        if(s.world_rank == 0)
        {
            printf("Working on cycle %d \n",p.cycle);
        }

	//reset counter for printing time stats
	s.counter = 0;

        //create a trajectory
        Trajectory traj; 

	//get filename for current cycle 
        string tag = "_cycle_" + to_string(p.cycle);
        string this_out_file_name = add_tag(p.out_file_name,tag);

	//set trajectory parameters
        traj.set_block_parallel(on);
        traj.set_traj(p.in_file_name);
        traj.set_ref(p.ref_file_name);
        traj.set_traj_w(this_out_file_name,p.b_print);
        traj.set_lsq(p.lsq_index_file_name,p.b_lsq,p.lsq_dim,p.lsq_ref);
        traj.set_res(p.stride,p.start_frame,p.end_frame);

        //analyze the trajectory (log time spent) 
        perf.log_time(traj.build(),"Analyze Trajectory");

        //print info about the worlk load distribution
        traj.workload();

        //print that analysis is beginning
        traj.report_progress();

        s.t = clock();
        //read read frames of the trajectory and perform analysis
        for(traj.current_frame=0; traj.current_frame<traj.get_num_frames(); traj.current_frame++)
        {
            traj.read_traj_frame();

            swap_coords(traj,s,p,param,p.cycle);

            traj.do_fit();

            traj.write_traj_frame();

            time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
        }

        //log time spent in main loop
        perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

        //splice temporary traj file together (log time spent)
        perf.log_time(traj.finalize_trajectory(),"Finalize Trajectory");

	MPI_Barrier(MPI_COMM_WORLD);
    }

    //print the performance stats
    perf.print_stats();

    //print closing statements
    print_closing(s.world_rank);

    //relinquish the mpi environment
    MPI_Finalize();

    return 0;
}

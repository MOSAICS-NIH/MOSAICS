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
#include "headers/mosat_routines.h"                          //This is where most of the functions called in main are located
#include "headers/file_naming.h"                             //This has routines for added tags to an existing file name    
#include "headers/file_naming_mpi.h"                         //This has routines for added tags to an existing file name (mpi)
#include "headers/command_line_args_mpi.h"                   //This has routines for adding command line arguments
#include "MosAT/program_variables/pv_pbc_gen.h"              //This has the variables specific to the analysis program
#include "headers/array.h"                                   //This has routines used for working with arrays
#include "headers/performance.h"                             //This has a class for logging performance data
#include "headers/index.h"                                   //This has a class for working with index files
#include "headers/traj.h"                                    //This has a class for working with the trajectory
#include "headers/leaflet_finder.h"                          //This has routines used to find leaflets in membrane simulations
#include "headers/protein_finder.h"                          //This has routines used to find protein atoms
#include "headers/sol_finder.h"                              //This has routines used to find the solvent
#include "headers/grid.h"                                    //This has routines used for working with a grid
#include "headers/protein.h"                                 //This has routines used for working with protein data
#include "headers/fit.h"                                     //This has routines used for fitting data
#include "headers/parallel.h"                                //This has routines for different parallelization schemes
#include "headers/force_serial.h"                            //This has routines used for forcing the code to run on a single mpi process
#include "headers/param.h"                                   //This has routines used for reading complex parameter data
#include "headers/grid_lt.h"                                 //This has routines used for reading in grid data
#include "headers/voronoi.h"                                 //This has routines used for computing voronoi diagrams
#include "headers/atom_select.h"                             //This has routines used for making atom selections using a selection text

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function places the molecules back inside the box                                                    //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void expand(Trajectory &traj,Trajectory &traj_expanded,system_variables &s,program_variables &p)
{
    int i          = 0;     //standard variable used in loops
    int x          = 0;     //standard variable used in loops
    int y          = 0;     //standard variable used in loops
    int counter    = 0;     //keeps track of the expanded atom index
    int itteration = 0;     //keeps track of the current copy index

    for(x=-1; x<=1; x++) //loop over x shifts
    {
        for(y=-1; y<=1; y++) //loop over y shifts
        { 
            for(i=0; i<traj.atoms(); i++) //loop over system atoms
            {
                traj_expanded.r[counter][0] = traj.r[i][0] + x*traj.box[XX][XX];
                traj_expanded.r[counter][1] = traj.r[i][1] + y*traj.box[YY][YY];
                traj_expanded.r[counter][2] = traj.r[i][2];

                counter++;
            }
        }
    }
    counter    = 0;
    itteration = 0;

    //expand the reference structure
    if(traj.current_frame == 0)
    {
        for(x=-1; x<=1; x++) //loop over x shifts
        {
            for(y=-1; y<=1; y++) //loop over y shifts
            {
                for(i=0; i<traj.atoms(); i++) //loop over system atoms
                {
                    traj_expanded.atom_name[counter] = traj.atom_name[i];
                    traj_expanded.res_name[counter]  = traj.res_name[i];
                    traj_expanded.atom_nr[counter]   = counter + 1;
                    traj_expanded.res_nr[counter]    = traj.res_nr[i] + itteration*traj.res_nr[traj.atoms()-1];

                    traj_expanded.beta[counter]      = traj.beta[i];                
                    traj_expanded.weight[counter]    = traj.weight[i];                
                    traj_expanded.element[counter]   = traj.element[i];                
                    traj_expanded.chain_id[counter]  = traj.chain_id[i];                

                    //printf("traj_expanded.atom_name[%d] %s \n",counter,traj_expanded.atom_name[counter].c_str());
                    counter++;
                }
                itteration++;
            }
        }
        counter    = 0;
        itteration = 0;

        for(x=-1; x<=1; x++) //loop over x shifts
        {
            for(y=-1; y<=1; y++) //loop over y shifts
            {
                for(i=0; i<traj.get_num_residues(); i++)
                {
                    traj_expanded.res_start[counter] = traj.res_start[i] + itteration*traj.atoms();
                    traj_expanded.res_end[counter]   = traj.res_end[i] + itteration*traj.atoms(); 
                }
                itteration++;
            }
        }
    }

    traj_expanded.box[XX][XX] = traj.box[XX][XX]*3.0;
    traj_expanded.box[YY][YY] = traj.box[YY][YY]*3.0;
    traj_expanded.box[ZZ][ZZ] = traj.box[ZZ][ZZ];
    traj_expanded.box[XX][YY] = traj.box[XX][YY];
    traj_expanded.box[XX][ZZ] = traj.box[XX][ZZ];
    traj_expanded.box[YY][XX] = traj.box[YY][XX];
    traj_expanded.box[YY][ZZ] = traj.box[YY][ZZ];
    traj_expanded.box[ZZ][XX] = traj.box[ZZ][XX];
    traj_expanded.box[ZZ][YY] = traj.box[ZZ][YY];

    for(i=0; i<200; i++)
    {
        traj_expanded.title[i] = traj.title[i];
    }

    traj_expanded.step = traj.step;
    traj_expanded.time = traj.time;
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
    s.program_name = "PBC Gen";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                                   s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                                   s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                                  s.world_rank, s.cl_tags, &p.b_print,   1);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                                    s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                                               s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                                s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                                  s.world_rank, s.cl_tags, &p.b_lsq,     0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                                    s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",                    s.world_rank, s.cl_tags, nullptr,      0);
    conclude_input_arguments_mpi(argc,argv,s.world_rank,s.program_name,s.cl_tags);

    //create a trajectory
    Trajectory traj; 

    //set trajectory parameters
    traj.set_block_parallel(on);
    traj.set_traj(p.in_file_name);
    traj.set_traj_w(p.out_file_name,p.b_print);
    traj.set_ref(p.ref_file_name);
    traj.set_lsq(p.lsq_index_file_name,p.b_lsq,p.lsq_dim,p.lsq_ref);
    traj.set_res(p.stride,p.start_frame,p.end_frame);

    //analyze the trajectory (log time spent) 
    perf.log_time(traj.build(),"Analyze Trajectory");

    //create a second trajectory for exanded system
    Trajectory traj_expanded;

    //set expanded trajectory parameters
    traj_expanded.set_block_parallel(on);
    traj_expanded.set_traj(p.in_file_name);
    traj_expanded.set_traj_w(p.out_file_name,p.b_print);
    traj_expanded.set_ref(p.ref_file_name);
    traj_expanded.set_lsq(p.lsq_index_file_name,p.b_lsq,p.lsq_dim,p.lsq_ref);
    traj_expanded.set_res(p.stride,p.start_frame,p.end_frame);

    //expand the trajectory
    traj_expanded.expand_traj(traj);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

        expand(traj,traj_expanded,s,p);

        traj_expanded.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
        fflush(stdin);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //splice temporary traj file together (log time spent)
    perf.log_time(traj_expanded.finalize_trajectory(),"Finalize Trajectory");

    //print the performance stats
    perf.print_stats();

    //print closing statements
    print_closing(s.world_rank);

    //relinquish the mpi environment
    MPI_Finalize();

    return 0;
}

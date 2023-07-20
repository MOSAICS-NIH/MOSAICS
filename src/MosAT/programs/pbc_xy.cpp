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
#include "MosAT/program_variables/pv_pbc_xy.h"              //This has the variables specific to the analysis program
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function removes jumps in the x or y directions                                                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void fix_pbc_xy(Trajectory &traj,system_variables &s,program_variables &p,double jump_x[],double jump_y[],double prev_x[],double prev_y[])
{
    int i = 0;  //standard variable used in loops
    int j = 0;  //standard variable used in loops

    if(traj.current_frame == 0) //store the coordinate
    {
        for(i=0; i<traj.atoms(); i++)
        {
            prev_x[i] = traj.r[i][0];
            prev_y[i] = traj.r[i][1];
        }
    }
    else if(traj.current_frame >= 1)
    {
        for(i=0; i<traj.atoms(); i++) //loop over system atoms
        {
            int min = traj.get_res_start(i);      //get the first atom of the current residue                                                         
            int max = traj.get_res_end(i);        //get the last atom of the current residue                                                                                                     
            i       = traj.next_residue(i);       //jump to the next residue                                                                                                    

            //check that the jump value is consistent over all atoms of the reside
            //if not this could inidacte that molecules are not whole throughout.
            for(j=min; j<=max; j++) //loop over current residue
            {
                if(jump_x[j] != jump_x[min])
                {
                    printf("broken molecule detected %s %d on frame %d jump_x[%d] %f jump_x[%d] %f \n",traj.res_name[j].c_str(),traj.res_nr[j],traj.current_frame,j,jump_x[j],min,jump_x[min]);
                }
                if(jump_y[j] != jump_y[min])
                {
                    printf("broken molecule detected %s %d on frame %d jump_y[%d] %f jump_y[%d] %f \n",traj.res_name[j].c_str(),traj.res_nr[j],traj.current_frame,j,jump_y[j],min,jump_y[min]);
                }
            }

            //unwrap the current frame
            for(j=min; j<=max; j++) //loop over current residue
            {
                //check jumping in x direction
                if(traj.r[j][0] - prev_x[j] > p.cutoff*traj.box[XX][XX]) //jumped accross box left side   
                {
                    jump_x[j] = jump_x[j] - traj.box[XX][XX];
                }
                else if(traj.r[j][0] - prev_x[j] < -1*p.cutoff*traj.box[XX][XX]) //jumped accross box right side
                {
                    jump_x[j] = jump_x[j] + traj.box[XX][XX];
                }

                //check jumping in y direction
                if(traj.r[j][1] - prev_y[j] > p.cutoff*traj.box[YY][YY]) //jumped accross box bottom side   
                {
                    jump_y[j] = jump_y[j] - traj.box[YY][YY];
                }
                else if(traj.r[j][1] - prev_y[j] < -1*p.cutoff*traj.box[YY][YY]) //jumped accross box top side
                {
                    jump_y[j] = jump_y[j] + traj.box[YY][YY];
                }

                //update prev coord
                prev_x[j] = traj.r[j][0];
                prev_y[j] = traj.r[j][1];

                //fix jumping   
                traj.r[j][0] = traj.r[j][0] + jump_x[j];
                traj.r[j][1] = traj.r[j][1] + jump_y[j];
            }
        }
    }
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
    s.program_name = "PBC-XY";

    //force program to run in serial?
    enum Switch serial         = on;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                                                          s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                                                          s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                                                         s.world_rank, s.cl_tags, &p.b_print,   1);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                                                           s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                                                                      s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                                                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                                                         s.world_rank, s.cl_tags, &p.b_lsq,     0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                                                           s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",                                           s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-cutoff", &p.cutoff,                    "Percentage of box dimension that an atom must moved between frames to be counted as a jump (0 to 1)", s.world_rank, s.cl_tags, nullptr,      1);
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
    //array telling if lipid jumped
    double jump_x[traj.atoms()];              //keeps track of how far jumps have moved the atoms
    double jump_y[traj.atoms()];              //keeps track of how far jumps have moved the atoms
    init_darray(jump_x,traj.atoms());
    init_darray(jump_y,traj.atoms());

    //array holding previous z 
    double prev_x[traj.atoms()];         //holds the x-coord from the previous step
    double prev_y[traj.atoms()];         //holds the y-coord from the previous step
    init_darray(prev_x,traj.atoms());
    init_darray(prev_y,traj.atoms());
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

        fix_pbc_xy(traj,s,p,jump_x,jump_y,prev_x,prev_y);

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

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

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
#include "MosAT/program_variables/pv_pbc_z.h"               //This has the variables specific to the analysis program
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
// This function removes jumps in the z direction                                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void fix_pbc_z(Trajectory &traj,system_variables &s,program_variables &p,int jump[],double prev_z[],iv2d &jump_record)
{
    int i = 0;      //standard variable used in loops
    int count = 0;  //count the non sol atoms

    if(traj.current_frame == 0) //store the coordinate
    {
        for(i=0; i<traj.atoms(); i++)
        {
            prev_z[i] = traj.r[i][2];
        }
    }
    else if(traj.current_frame >= 1)
    {
        for(i=0; i<traj.atoms(); i++) //loop over system atoms
        {
            if(traj.b_sol[i] == 0) //atom is not the solvent
            {
                if(traj.r[i][2] - prev_z[i] > p.cutoff*traj.box[ZZ][ZZ]) //jumped accross bottom of box   
                {
                    jump[i] = jump[i] - 1;

                    //report jump
                    p.jump_count = p.jump_count + 1;
                    if(p.jump_count < 100)
                    {
                        printf("Jump %3d detected between frames %9d and %9d for atom %10d. Jump %d \n",p.jump_count,traj.get_frame_global()-1,traj.get_frame_global(),i,jump[i]);
                    }
                    if(p.jump_count == 99)
                    {
                        printf("100 jumps detected. Will no longer report jumps. \n");
                    }
                }
                else if(traj.r[i][2] - prev_z[i] < -1*p.cutoff*traj.box[ZZ][ZZ]) //jumped accross top of box
                {
                    jump[i] = jump[i] + 1;

                    //report jump
                    p.jump_count = p.jump_count + 1;
                    if(p.jump_count < 100)
                    {
                        printf("Jump %3d detected between frames %9d and %9d for atom %10d. Jump %d \n",p.jump_count,traj.get_frame_global()-1,traj.get_frame_global(),i,jump[i]);
                    }
                    if(p.jump_count == 99)
                    {
                        printf("100 jumps detected. Will no longer report jumps. \n");
                    }
                }
               
                //record the jump state 
                if(p.b_record == 1)
                {
                    jump_record[traj.get_frame_global()][count] = jump[i];
                }

                //update prev coord
                prev_z[i] = traj.r[i][2];

                //fix jumping   
                traj.r[i][2] = traj.r[i][2] + jump[i]*traj.box[ZZ][ZZ];
 
                count++;
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes the jump record to file                                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double print_record(Trajectory &traj,system_variables &s,program_variables &p,iv2d &jump_record)
{
    int i = 0;
    int j = 0;

    s.t = clock();

    printf("\nReporting jump records. \n"); 

    if(p.b_record == 1)
    {
        FILE *record_file = fopen(p.record_file_name.c_str(), "w");
        if(record_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",p.record_file_name.c_str());
        }
        else
        {
            for(i=0; i<jump_record[0].size(); i++) //loop over atoms
            {
                for(j=0; j<jump_record.size(); j++) //loop over frames
                {
                    fprintf(record_file," %10d ",jump_record[j][i]);
                }
                fprintf(record_file,"\n");
            }
            fclose(record_file);
        }
    }

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
    s.program_name = "PBC-Z";

    //force program to run in serial?
    enum Switch serial         = on;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                                                       s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                                                       s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                                                      s.world_rank, &p.b_print,   1);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                                                        s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                                                                   s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                                                    s.world_rank, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                                                      s.world_rank, &p.b_lsq,     0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                                                        s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",                                        s.world_rank, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-cutoff", &p.cutoff,                    "Percentage of box height that an atom must moved between frames to be counted as a jump (0 to 1)", s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-sf_pdb", p.sf_pdb_file_name,           "PDB file with selected sol (pdb)",                                                                 s.world_rank, &p.b_sf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-sf_prm", p.solvent_finder_param_name,  "File with additional solvent finder parameters (prm)",                                             s.world_rank, &p.b_sf_param,0);
    add_argument_mpi_s(argc,argv,"-record", p.record_file_name,           "Output data file with the jump records (dat)",                                                     s.world_rank, &p.b_record,  0);
    conclude_input_arguments_mpi(argc,argv,s.world_rank,s.program_name);

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
    if(p.b_sf_pdb == 1)
    {
        check_extension_mpi(s.world_rank,"-sf_pdb",p.sf_pdb_file_name,".pdb");
    }
    if(p.b_sf_param == 1)
    {
        check_extension_mpi(s.world_rank,"-sf_prm",p.solvent_finder_param_name,".prm");
    }
    if(p.b_record == 1)
    {
        check_extension_mpi(s.world_rank,"-record",p.record_file_name,".dat");
    }
    //array telling if lipid jumped
    int jump[traj.atoms()];              //tells if an atom has jumped (pbc)
    init_iarray(jump,traj.atoms());

    //array holding previous z 
    double prev_z[traj.atoms()];         //holds the z-coord from the previous step
    init_darray(prev_z,traj.atoms());

    //run solvent finder
    traj.get_solvent(p.solvent_finder_param_name,p.b_sf_param);

    //print pdb with distinguished solven 
    traj.write_sol(p.sf_pdb_file_name,p.b_sf_pdb);

    //print info about the water
    traj.get_sol_stats();

    //allocate memory for jump record
    iv2d jump_record;
    if(p.b_record == 1)
    {
        //print memory estimate for jump record
        int   size = traj.get_ef_frames()*(traj.atoms()-traj.sol.size());
        double mem = 4.0*(double)size/1000000.0;
        printf("Estimated memory for jump record: %f MB \n\n",mem);

        jump_record.resize(traj.get_ef_frames());
        for(p.i=0; p.i<traj.get_ef_frames(); p.i++)
        {
            jump_record[p.i].resize(traj.atoms()-traj.sol.size(),0);
        }
    }

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

        fix_pbc_z(traj,s,p,jump,prev_z,jump_record);

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //report the jump record 
    perf.log_time(print_record(traj,s,p,jump_record),"Jump Record");

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

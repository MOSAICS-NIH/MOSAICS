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
#include "headers/file_naming_mpi.h"                         //This has routines for added tags to an existing file name mpi   
#include "headers/file_naming.h"                             //This has routines for added tags to an existing file name    
#include "headers/command_line_args_mpi.h"                   //This has routines for adding command line arguments
#include "MosAT/program_variables/pv_atomic_rmsf.h"          //This has the variables specific to the analysis program
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
// This function computes the diviation from ref                                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_dev(Trajectory &traj,system_variables &s,program_variables &p,dv2d &dev)
{
    int     i = 0;       //standard variable used in loops

    for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
    {
        dv1d a; //hold the coordinates of the atom
        dv1d b; //hold the coordinates of the ref atom

        //get coordinates for the atom
        a.push_back(traj.r[traj.prot[i]-1][0]);
        a.push_back(traj.r[traj.prot[i]-1][1]);
        a.push_back(traj.r[traj.prot[i]-1][2]);

        //get coordinates for the ref atom
        b.push_back(traj.r_ref[traj.prot[i]-1][0]);
        b.push_back(traj.r_ref[traj.prot[i]-1][1]);
        b.push_back(traj.r_ref[traj.prot[i]-1][2]);

        //add deviation to long term sum
        dev[i][traj.current_frame] = traj.get_dist(a,b);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects deviation data and writes the RMSD to the B-factor of a PDB file.                  //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,dv2d &dev)
{
    int i             = 0;         //standard variable used in loops  
    int j             = 0;         //standard variable used in loops 

    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    //initialize b-factor
    for(i=0; i<traj.atoms(); i++)
    {
        traj.beta[i] = 0.0;
    }

    //collect deviations and compute rmsf
    for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
    {
        //collect deviations for current atom
        collect_dv1d(s.world_size,s.world_rank,dev[i]);

        if(s.world_rank == 0)
        {
            double rmsf = 0;
            for(j=0; j<traj.get_ef_frames(); j++) //loop over trajectory frames
            {
                rmsf = rmsf + pow(dev[i][j],2.0);
            }
            rmsf = sqrt(rmsf/(double)traj.get_ef_frames());
 
            //set b-factor
            traj.beta[traj.prot[i]-1] = rmsf;
        }
    }

    //write b-factor data to PDB file
    if(s.world_rank == 0)
    {
        //open file for writing
        FILE *pdb_file = fopen(p.rmsf_file_name.c_str(), "w");
        if(pdb_file == NULL)
        {
            printf("failure opening %s (pdb single). Make sure the file exists. \n",p.rmsf_file_name.c_str());
        }
        write_frame_pdb(traj.ibox,traj.atoms(),traj.atom_nr,traj.res_nr,traj.res_name,traj.atom_name,traj.r_ref,traj.title,s.world_rank,&pdb_file,traj.beta,traj.weight,traj.element,traj.chain_id,1);
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
    s.program_name = "Atomic RMSF";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                               s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (tpr, gro)",                                               s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                              s.world_rank, s.cl_tags, &p.b_print,   0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                                s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                                           s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                            s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting",                                                    s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                                s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",                s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-pf_pdb", p.pf_pdb_file_name,           "PDB file with selected protein ",                                          s.world_rank, s.cl_tags, &p.b_pf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-pf_prm", p.protein_finder_param_name,  "File with additional protein finder parameters ",                          s.world_rank, s.cl_tags, &p.b_pf_param,0);
    add_argument_mpi_s(argc,argv,"-rmsf",   p.rmsf_file_name,              "Output file with RMSF data (PDB)",                                        s.world_rank, s.cl_tags, nullptr,      1);
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
    check_extension_mpi(s.world_rank,"-rmsf",p.rmsf_file_name,".pdb");
    if(p.b_pf_pdb == 1)
    {
        check_extension_mpi(s.world_rank,"-pf_pdb",p.pf_pdb_file_name,".pdb");
    }
    if(p.b_pf_param == 1)
    {
        check_extension_mpi(s.world_rank,"-pf_prm",p.protein_finder_param_name,".prm");
    }

    //create index for reading lsq index 
    Index lsq_index;

    //read the index files
    lsq_index.get_index(p.lsq_index_file_name);

    //run leaflet/proten/solvent finder
    traj.get_protein(p.protein_finder_param_name,p.b_pf_param);

    //print a pdb with distinguished leaflets/protein/solvent
    traj.write_protein(p.pf_pdb_file_name,p.b_pf_pdb);

    //create structure to hold deviation from ref values
    dv2d dev(traj.prot.size(),dv1d(traj.get_num_frames(),0.0));

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

        //move refernce structure to origin
        if(traj.current_frame == 0)
        {
            traj.place_ref_at_origin(lsq_index,p.lsq_dim,p.lsq_ref);
        }
	
        //do fit but leave at the origin
        traj.do_this_fit(lsq_index,p.lsq_dim,p.lsq_ref);

        get_dev(traj,s,p,dev);

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect deviations from mpi processes and write output
    perf.log_time(finalize_analysis(traj,s,p,dev),"Fin Ana");

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

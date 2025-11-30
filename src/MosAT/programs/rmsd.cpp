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
#include "MosAT/program_variables/pv_rmsd.h"                 //This has the variables specific to the analysis program
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
#include "headers/atom_select.h"                             //This has routines used for making atom selections using a selection text

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function calculates the center for the reference structure                                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
dv1d get_ref_center(Trajectory &traj,system_variables &s,program_variables &p,Index &rmsd_index)
{
    int    i         = 0;     //standard variable used in loops
    double cx        = 0;     //used to add up the x component from each atom
    double cy        = 0;     //used to add up the y component from each atom
    double cz        = 0;     //used to add up the z component from each atom

    dv1d   r_center(3,0.0); //the geometric center

    for(i=0; i<rmsd_index.index_i.size(); i++) //loop over group of atoms
    {
        cx = cx + traj.r_ref[rmsd_index.index_i[i]-1][0];
        cy = cy + traj.r_ref[rmsd_index.index_i[i]-1][1];
        cz = cz + traj.r_ref[rmsd_index.index_i[i]-1][2];
    }
    r_center[0] = cx/(double)rmsd_index.index_i.size();
    r_center[1] = cy/(double)rmsd_index.index_i.size();
    r_center[2] = cz/(double)rmsd_index.index_i.size();

    return r_center;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function calculates the RMSD for the current frame                                                   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_rmsd(Trajectory &traj,system_variables &s,program_variables &p,dv1d &rmsd,Index &rmsd_index,Index &lsq_index)
{
    int i = 0;     //standard variable used in loops

    //copy coordinates of ref structure
    rvec *r_ref_tmp;                                          //holds ref coords copy
    r_ref_tmp = (rvec*)calloc(traj.atoms() , sizeof(rvec));   //allocate memory
    for(i=0; i<traj.atoms(); i++)                             //loop over atoms
    {
        r_ref_tmp[i][0] = traj.r_ref[i][0];
        r_ref_tmp[i][1] = traj.r_ref[i][1];
        r_ref_tmp[i][2] = traj.r_ref[i][2];
    }

    //compute centers for moving ref and frame to origin
    dv1d center     = traj.center_i(lsq_index.index_i);
    dv1d center_ref = get_ref_center(traj,s,p,lsq_index);

    //center frame at origin 
    for(i=0; i<traj.atoms(); i++) //loop over atoms
    {
        traj.r[i][0] = traj.r[i][0] - center[0];
        traj.r[i][1] = traj.r[i][1] - center[1];
        traj.r[i][2] = traj.r[i][2] - center[2];

        traj.r_ref[i][0] = traj.r_ref[i][0] - center_ref[0];
        traj.r_ref[i][1] = traj.r_ref[i][1] - center_ref[1];
        traj.r_ref[i][2] = traj.r_ref[i][2] - center_ref[2];
    }

    //tag atoms for fitting
    real   dummy_mass[traj.atoms()];
    for(i=0; i<traj.atoms(); i++)
    {
        dummy_mass[i] = 0.0;
    }
    for(i=0; i<lsq_index.index_i.size(); i++)
    { 
        dummy_mass[lsq_index.index_i[i]-1] = 1.0;
    }

    //do lsq fitting
    do_fit_ndim(3, traj.atoms(), dummy_mass, traj.r_ref, traj.r);

    //tag atoms for rmsd computation 
    for(i=0; i<traj.atoms(); i++)
    {
        dummy_mass[i] = 0.0;
    }
    for(i=0; i<rmsd_index.index_i.size(); i++)
    {
        dummy_mass[rmsd_index.index_i[i]-1] = 1.0;
    }

    //compute rmsd
    rmsd.push_back(rmsdev(traj.atoms(), dummy_mass, traj.r, traj.r_ref));

    //restore coords for ref structure
    for(i=0; i<traj.atoms(); i++)
    {
        traj.r_ref[i][0] = r_ref_tmp[i][0];
        traj.r_ref[i][1] = r_ref_tmp[i][1];
        traj.r_ref[i][2] = r_ref_tmp[i][2];
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Collect RMSD data                                                                                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,dv1d &rmsd)
{
    int i = 0; //standard variable used in loops
    int j = 0; //standard variable used in loops 

    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nfinalizing analysis. This requires communication of the RMSD data and could take some time. \n");
    }

    collect_dv1d(s.world_size,s.world_rank,rmsd);

    if(s.world_rank == 0)
    {
        FILE *rmsd_file = fopen(p.rmsd_file_name.c_str(),"w");
        fprintf(rmsd_file," %9s   %10s \n","#step","#rmsd_(nm)");
        for(i=0; i<rmsd.size(); i++)
        {
            fprintf(rmsd_file," %9d   %10f \n",i,rmsd[i]);          
        }
	fclose(rmsd_file);
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
    s.program_name = "RMSD";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                      s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                      s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                     s.world_rank, s.cl_tags, &p.b_print,   0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                                  s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                   s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                     s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-rmsd",   p.rmsd_file_name,             "Output data file (dat) with RMSD data",                           s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-n",      p.index_file_name,            "Index file (ndx) for calculating RMSD",                           s.world_rank, s.cl_tags, nullptr,      1);
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
    check_extension_mpi(s.world_rank,"-rmsd",p.rmsd_file_name,".dat");
    check_extension_mpi(s.world_rank,"-n",p.index_file_name,".ndx");
    check_extension_mpi(s.world_rank,"-lsq",p.lsq_index_file_name,".ndx");

    //create index for lipid types and atom types
    Index rmsd_index;
    Index lsq_index;

    //read the index files
    rmsd_index.get_index(p.index_file_name);
    lsq_index.get_index(p.lsq_index_file_name);

    //holds the rmsd for each frame
    dv1d rmsd(0);

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

        get_rmsd(traj,s,p,rmsd,rmsd_index,lsq_index);

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //splice temporary traj file together (log time spent)
    perf.log_time(traj.finalize_trajectory(),"Finalize Trajectory");

    //finalize analysis 
    perf.log_time(finalize_analysis(traj,s,p,rmsd),"Finilize Ana.");

    //print the performance stats
    perf.print_stats();

    //print closing statements
    print_closing(s.world_rank);

    //relinquish the mpi environment
    MPI_Finalize();

    return 0;
}

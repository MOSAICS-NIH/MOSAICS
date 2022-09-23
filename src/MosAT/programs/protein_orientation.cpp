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
#include "headers/linear_algebra.h"                          //This has linear algebra stuff
#include "headers/mosat_routines.h"                         //This is where most of the functions called in main are located
#include "headers/file_naming.h"                             //This has routines for added tags to an existing file name    
#include "headers/file_naming_mpi.h"                         //This has routines for added tags to an existing file name (mpi)
#include "headers/command_line_args_mpi.h"                   //This has routines for adding command line arguments
#include "MosAT/program_variables/pv_protein_orientation.h" //This has the variables specific to the analysis program
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
// This function writes a frame to a pdb file adding dummy atoms to represent a vector                       //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_frame_pdb_vector(Trajectory &traj,system_variables &s,program_variables &p,int global_frame,double com_u[],double com_l[],double v[],FILE **pdb_vec_file)
{
    int i = 0;    //standard variable used in loops

    //print the title
    fprintf(*pdb_vec_file,"%-6s    %s","TITLE",traj.title);

    //print the box
    gmx_write_pdb_box(*pdb_vec_file,-1,traj.box);

    //print the model number
    fprintf(*pdb_vec_file,"MODEL%9d\n",global_frame);

    //print the atom lines
    for(i=0; i<traj.atoms(); i++)
    {
        fprintf_atomline_pdb(pdb_vec_file,"ATOM",traj.atom_nr[i],traj.atom_name[i].c_str(),' ',traj.res_name[i].c_str(),traj.chain_id[i],traj.res_nr[i],' ',10.0*traj.r[i][0],10.0*traj.r[i][1],10.0*traj.r[i][2],traj.weight[i],traj.beta[i],traj.element[i].c_str());
    }

    //print dummy atoms
    fprintf(*pdb_vec_file,"ATOM%7d  %-3s %3s%6d    %8.3f%8.3f%8.3f  %4.2f  %4.2f\n",traj.atoms()+1,"DUM","DUM",traj.res_nr[traj.atoms()-1]+1,10*com_u[0],10*com_u[1],10*com_u[2],1.00,0.00);
    fprintf(*pdb_vec_file,"ATOM%7d  %-3s %3s%6d    %8.3f%8.3f%8.3f  %4.2f  %4.2f\n",traj.atoms()+2,"DUM","DUM",traj.res_nr[traj.atoms()-1]+2,10*com_l[0],10*com_l[1],10*com_l[2],1.00,0.00);
    fprintf(*pdb_vec_file,"ATOM%7d  %-3s %3s%6d    %8.3f%8.3f%8.3f  %4.2f  %4.2f\n",traj.atoms()+3,"DUM","DUM",traj.res_nr[traj.atoms()-1]+3,0.0,0.0,0.0,1.00,0.00);
    fprintf(*pdb_vec_file,"ATOM%7d  %-3s %3s%6d    %8.3f%8.3f%8.3f  %4.2f  %4.2f\n",traj.atoms()+4,"DUM","DUM",traj.res_nr[traj.atoms()-1]+4,10*v[0],10*v[1],10*v[2],1.00,0.00);

    //print the TER and ENDMDL tafs
    fprintf(*pdb_vec_file,"TER\n");
    fprintf(*pdb_vec_file,"ENDMDL\n");
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the orientation angles theta and phi for the protein.                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_orientation(Trajectory &traj,system_variables &s,program_variables &p,Index &upper,Index &lower,double theta[],double phi[])
{
    int i = 0;               //standard variable used in loops
    int j = 0;               //standard variable used in loops
    double com_u[3];         //center of mass for upper selection
    double com_l[3];         //center of mass for lower selection
    double v[3];             //the orientation vector
    init_darray(com_u,3);
    init_darray(com_l,3);

    //reset com
    com_u[0] = 0;
    com_u[1] = 0;
    com_u[2] = 0;
    com_l[0] = 0;
    com_l[1] = 0;
    com_l[2] = 0;

    //now we compute the upper plane center of mass
    for(i=0; i<upper.index_i.size(); i++) //loop over upper plane atoms
    {
        for(j=0; j<3; j++) //loop over 3 dimensions
        {
            com_u[j] = com_u[j] + traj.r[upper.index_i[i]-1][j]/upper.index_i.size();
        }
    }
    //now we compute the lower plane center of mass
    for(i=0; i<lower.index_i.size(); i++) //loop over lower plane atoms
    {
        for(j=0; j<3; j++) //loop over 3 dimensions
        {
            com_l[j] = com_l[j] + traj.r[lower.index_i[i]-1][j]/lower.index_i.size();
        }
    }

    //compute the orientation vector
    for(i=0; i<3; i++)
    {
        v[i] = com_u[i] - com_l[i];
    }

    //print the orientatioin vector 
    if(p.print_vec == 1)
    {
        printf("rank %3d frame %5d orientation_x %10.3f orientation_y %10.3f orientation_z %10.3f \n",s.world_rank,traj.current_frame,v[0],v[1],v[2]);
    }

    //compute theta
    double z[3]          = {0, 0, 1};
    double mag_v         = sqrt( dot(v,v) );
    double mag_z         = sqrt( dot(z,z) );
    double v_dot_z       = dot(v,z);
    double cos_theta     = v_dot_z/(mag_v*mag_z);
    double pi            = 3.141592654;
    theta[traj.current_frame] = acos(cos_theta)*(180/pi);

    //compute phi
    double proj_v[3]    = {v[0],v[1],0};
    double mag_proj_v   = sqrt( dot(proj_v,proj_v) );
    double y[3]         = {0,1,0};
    double mag_y        = sqrt( dot(y,y) );
    double proj_v_dot_y = dot(proj_v,y);
    double cos_phi      = proj_v_dot_y/(mag_proj_v*mag_y);
    phi[traj.current_frame]  = acos(cos_phi)*(180/pi);

    if(v[0] < 0)
    {
        phi[traj.current_frame] = 2*(180 - phi[traj.current_frame]) + phi[traj.current_frame];
    }

    //now we print the pdb with the orientation vector
    if(p.print_vec_pdb == 1 && s.world_rank == 0)
    {
        int global_frame = 0;
        for(i=0; i<s.world_rank; i++)
        {
            global_frame = global_frame + traj.get_num_frames_world()[i];
        }
        global_frame = global_frame + traj.current_frame;

        string out_ori_vec_file_name = p.orientation_file_name;
        out_ori_vec_file_name.pop_back();
        out_ori_vec_file_name.pop_back();
        out_ori_vec_file_name.pop_back();
        out_ori_vec_file_name.pop_back();
        out_ori_vec_file_name = out_ori_vec_file_name + "_" + "ori_vec" + ".pdb";

        FILE *pdb_vec_file = fopen(out_ori_vec_file_name.c_str(), "w");
        if(pdb_vec_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",out_ori_vec_file_name.c_str());
        }

        write_frame_pdb_vector(traj,s,p,0,com_u,com_l,v,&pdb_vec_file);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects the orientation angles theta and phi from the nodes and prints them to the output  //
// file.                                                                                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,double theta[],double phi[])
{
    int i                = 0;                  //standard variable used in loops
    int j                = 0;                  //standard variable used in loops
    int k                = 0;                  //standard variable used in loops
    double average_theta = 0;                  //the average theta
    double world_theta[traj.get_ef_frames()];  //used for collecting theta
    double world_phi[traj.get_ef_frames()];    //used for collecting phi
    init_darray(world_theta,traj.get_ef_frames());
    init_darray(world_phi,traj.get_ef_frames());

    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nFinalizing analysis. This could take some time. \n");
    }

    //first we collect the anlges from nodes and print them to file
    collect_1d_double_array(s.world_size,traj.get_num_frames_world(),s.world_rank,world_theta,theta);
    collect_1d_double_array(s.world_size,traj.get_num_frames_world(),s.world_rank,world_phi,phi);

    if(s.world_rank == 0)
    {
        fprintf(p.orientation_file,"%10s %10s %10s \n","traj frame","theta(deg)","phi(deg)");
        for(i=0; i<traj.get_ef_frames(); i++)
        {
            average_theta = average_theta + world_theta[i];
            fprintf(p.orientation_file,"%10d %10.3f %10.3f \n",i,world_theta[i],world_phi[i]);
        }
        average_theta = average_theta/traj.get_ef_frames();
        printf("average theta %10.4f \n",average_theta);
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
    s.program_name = "Protein Orientation";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                                 s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                                 s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                                s.world_rank, &p.b_print,   0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                                  s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                                             s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                              s.world_rank, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                                s.world_rank, &p.b_lsq,     0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                                  s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",                  s.world_rank, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-upper",  p.upper_file_name,            "Index file with upper atoms (ndx)",                                          s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lower",  p.lower_file_name,            "Index file with lower atoms (ndx)",                                          s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ori",    p.orientation_file_name,      "Output data file with theta and phi angles for each trajectory frame (dat)", s.world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-p_vec",  &p.print_vec,                 "Print orientation vec (0:no 1:yes)",                                         s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-o_pdb",  &p.print_vec_pdb,             "print orientation vec to pdb (0:no 1:yes)",                                  s.world_rank, nullptr,      0);
    conclude_input_arguments_mpi(argc,argv,s.world_rank,s.program_name);

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
    check_extension_mpi(s.world_rank,"-upper",p.upper_file_name,".ndx");
    check_extension_mpi(s.world_rank,"-lower",p.lower_file_name,".ndx");
    check_extension_mpi(s.world_rank,"-ori",p.orientation_file_name,".dat");

    //open additional files here
    p.orientation_file = fopen(p.orientation_file_name.c_str(), "w");
    if(p.orientation_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",p.orientation_file_name.c_str());
    }

    //create index objects
    Index upper;
    Index lower;

    //read the index files
    upper.get_index(p.upper_file_name);
    lower.get_index(p.lower_file_name);

    double theta[traj.get_num_frames()];  //holds theta
    double phi[traj.get_num_frames()];    //holds phi
    init_darray(theta,traj.get_num_frames());
    init_darray(phi,traj.get_num_frames());
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

        get_orientation(traj,s,p,upper,lower,theta,phi);

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect angles from mpi processes and print output data
    perf.log_time(finalize_analysis(traj,s,p,theta,phi),"Fin Ana");

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

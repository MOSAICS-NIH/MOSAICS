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
#include "MosAT/program_variables/pv_p2_3d.h"               //This has the variables specific to the analysis program
#include "headers/array.h"                                   //This has routines used for working with arrays
#include "headers/performance.h"                             //This has a class for logging performance data
#include "headers/index.h"                                   //This has a class for working with index files
#include "headers/traj.h"                                    //This has a class for working with the trajectory
#include "headers/leaflet_finder.h"                          //This has routines used to find leaflets in membrane simulations
#include "headers/protein_finder.h"                          //This has routines used to find protein atoms
#include "headers/sol_finder.h"                              //This has routines used to find the solvent
#include "headers/grid.h"                                    //This has routines used for working with a grid
#include "headers/grid_3d.h"                                 //This has routines used for working with a grid
#include "headers/protein.h"                                 //This has routines used for working with protein data
#include "headers/force_serial.h"                            //This has routines used for forcing the code to run on a single mpi process
#include "headers/param.h"                                   //This has routines used for reading complex parameter data

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the second rank order parameter P2 and adds it to the grid                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void p2(Trajectory &traj,system_variables &s,program_variables &p,Param &param,Grid_3d &order,Grid_3d &theta)
{
    //use cos_theta = u_dot_v/(mag_u*mag_v) to compute the angle between 2 vecors u and v
    int i = 0;                                    //standard variable used in loops
    int j = 0;                                    //standard variable used in loops
    int k = 0;                                    //standard variable used in loops
    int l = 0;                                    //standard variable used in loops
    int m = 0;                                    //standard variable used in loops
    int n = 0;                                    //standard variable used in loops
    double pi = 3.141592654;                      //pi

    //clear the current frame grids
    order.clean_frame();
    theta.clean_frame();

    for(i=0; i<traj.target_leaflet.size(); i++) //loop over target_leaflet atoms
    {
        //get the first and last atoms of the current lipid
        int min = traj.t_lip_start(i);
        int max = traj.t_lip_end(i);

        //jump to the next lipid
        i = traj.next_target_lipid(i);

        for(j=0; j<param.main_size_y(); j++) //loop over lipid types 
        {
            if(strcmp(traj.res_name[min].c_str(), param.param_main_s[j][0].c_str() ) == 0) //lipid type is correct
            {
                double avg_p2    = 0.0;
                double avg_theta = 0.0;
                int count        = 0;

                for(k=min; k<=max; k++) //loop over current residue atoms
                {
                    for(l=0; l<param.sec_size_y(j); l++) //loop over target atoms
                    {
                        if(strcmp(traj.atom_name[k].c_str(), param.param_sec_s[j][l][0].c_str() ) == 0) //first atom of bond
                        {
                            for(m=min; m<=max; m++) //loop over current residue atoms
                            {
                                if(strcmp(traj.atom_name[m].c_str(), param.param_sec_s[j][l][1].c_str() ) == 0) //second atom of bond
                                {
                                    double delta[3] = {0, 0, 0};
                                    double z[3]     = {0, 0, 1};

                                    for(n=0; n<3; n++) //loop over 3 dimensoins
                                    {
                                        delta[n] = traj.r[m][n] - traj.r[k][n];
                                    }

                                    double mag_delta   = sqrt(dot(delta,delta));
                                    double delta_dot_z = dot(delta,z);
                                    double cos_theta   = delta_dot_z/mag_delta;
                                    double theta       = acos(cos_theta)*(180/pi); //angle in degrees
                                    double p2          = 0.5*( (3*(cos_theta)*(cos_theta)) - 1);

                                    avg_p2    = avg_p2    + p2;
                                    avg_theta = avg_theta + theta;
                                    count++;
                                }
                            }
                        }
                    }
                }
                avg_p2    = avg_p2/(double)count;
                avg_theta = avg_theta/(double)count;

                //add data to mapping atoms
                for(k=min; k<=max; k++) //loop over current residue atoms
                {
                    if(strcmp(traj.atom_name[k].c_str(), param.param_main_s[j][1].c_str() ) == 0) //mapping atom 1
                    {
                        double hx = traj.r[k][0];
                        double hy = traj.r[k][1];
                        double hz = traj.r[k][2];
                        order.stamp(hx,hy,hz,p.radius,avg_p2);
                        theta.stamp(hx,hy,hz,p.radius,avg_theta);
                    }
                    else if(strcmp(traj.atom_name[k].c_str(), param.param_main_s[j][2].c_str() ) == 0) //mapping atom 2
                    {
                        double hx = traj.r[k][0];
                        double hy = traj.r[k][1];
                        double hz = traj.r[k][2];
                        order.stamp(hx,hy,hz,p.radius,avg_p2);
                        theta.stamp(hx,hy,hz,p.radius,avg_theta);
                    }
                }
            }
        }
    }

    //get the average for the current frame
    order.norm_frame();
    theta.norm_frame();

    //add the average single frame p2 parameter to long term sum
    order.add_frame();
    theta.add_frame();

    //now we print the single frame p2 data
    if(p.b_stdev == 1)
    {
        order.write_frame(traj.get_frame_global(),p.ex_val);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// collect p2 and compute average                                                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,Grid_3d &order,Grid_3d &theta)
{
    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nFinalizing analysis. This requires communicating data over the grid and could take some time depending on the resolution. \n");
    }

    //collect p2, theta and rho from all ranks
    order.collect_grid();
    theta.collect_grid();

    //normalize p2 and write the <p2> and rho to file
    if(s.world_rank == 0)
    {
        order.normalize();
        theta.normalize();

        order.exclude_data(p.cutoff,1);
        theta.exclude_data(p.cutoff,0);

        order.write_grid(p.ex_val);
        order.write_rho(p.ex_val);
        theta.write_grid(p.ex_val);
    }

    //compute standard deviation
    order.get_stdev(p.b_stdev,p.b_clean,traj,p.ex_val);

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
    s.program_name = "P2 3D";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                       s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                       s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                      s.world_rank, s.cl_tags, &p.b_print,   0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                        s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                                   s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                    s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                      s.world_rank, s.cl_tags, &p.b_lsq,     0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                        s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",        s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-bond",   p.bond_file_name,             "Selection card with lipid types and bond definitions, etc. (crd)", s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-p2",     p.p2_file_name,               "Output grid data with the spatially resolved time average (dx)",   s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,           "PDB file with sorted leaflets (pdb)",                              s.world_rank, s.cl_tags, &p.b_lf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,  "File with additional leaflet finder parameters (prm)",             s.world_rank, s.cl_tags, &p.b_lf_param,0);
    add_argument_mpi_d(argc,argv,"-APS",    &p.APS,                       "Area per grid box (nm^2)",                                         s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-r",      &p.radius,                    "Radius of the target atom (nm)",                                   s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-cutoff", &p.cutoff,                    "Cutoff for excluding data (chi)",                                  s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-bx",     &p.box_x,                     "Grid x dimension (nm)",                                            s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-by",     &p.box_y,                     "Grid y dimension (nm)",                                            s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-bz",     &p.box_z,                     "Grid z dimension (nm)",                                            s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-leaf",   &p.leaflet,                   "Which leaflet? (0:both 1:upper 2:lower)",                          s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-stdev",  &p.b_stdev,                   "Compute the STDEV and STEM? (0:no 1:yes)",                         s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-clean",  &p.b_clean,                   "Remove single frame files? (0:no 1:yes)",                          s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-ex_val", &p.ex_val,                    "Set excluded lattice points to this value",                        s.world_rank, s.cl_tags, nullptr,      0);
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
    check_extension_mpi(s.world_rank,"-bond",p.bond_file_name,".crd");
    check_extension_mpi(s.world_rank,"-p2",p.p2_file_name,".dx");

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
    param.get_param(p.bond_file_name,4,3,2);

    //run leaflet finder
    traj.get_leaflets(p.leaflet,p.leaflet_finder_param_name,p.b_lf_param);

    //print a pdb with distinguished leaflets
    traj.write_leaflets(p.lf_pdb_file_name,p.b_lf_pdb);

    //print info about the target leaflet
    traj.get_leaflet_stats();

    //print info about the lipid selection
    traj.get_lipid_selection_stats(param.get_column_s(0),"-bond");

    //create a grid to hold p2 and theta
    Grid_3d order;
    Grid_3d theta;

    //get the grid dimensions
    order.get_dim(p.box_x,p.box_y,p.box_z,traj.ibox,p.APS);
    if(s.world_rank == 0) //attempt to estimate requirements befor a crash
    {
        float mem_d = 4.0*(float)order.num_x()*(float)order.num_y()*(float)order.num_z()*8.0;
        float mem_i = 2.0*(float)order.num_x()*(float)order.num_y()*(float)order.num_z()*4.0;
        float mem_t = 2.0*(mem_d + mem_i);
        printf("Estimated memory to hold the grid: %f (MB) \n\n",mem_t/1000000.0);
    }
    theta.get_dim(p.box_x,p.box_y,p.box_z,traj.ibox,p.APS);

    //create file name for theta
    p.theta_file_name = add_tag(p.p2_file_name,"_theta");

    //set the output file name for grid
    order.set_output(p.p2_file_name);
    theta.set_output(p.theta_file_name);

    //print info about the grid
    order.print_dim();
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

        p2(traj,s,p,param,order,theta);

        traj.set_beta_lf();

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect the p2 values from mpi processes and compute the final average
    perf.log_time(finalize_analysis(traj,s,p,order,theta),"Fin Ana");

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

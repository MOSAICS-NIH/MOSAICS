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
#include "MosAT/program_variables/pv_3d_enrichment.h"       //This has the variables specific to the analysis program
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
// This function adds the lipid density to the grid for current frame                                        //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void enrichment(Trajectory &traj,system_variables &s,program_variables &p,Param &param_1,Param &param_2,Grid_3d &rho1,Grid_3d &rho2,Grid_3d &rho_t)
{
    int i = 0;                                   //standard variable used in loops
    int j = 0;                                   //standard variable used in loops
    int k = 0;                                   //standard variable used in loops
    int l = 0;                                   //standard variable used in loops

    //clear the current frame grids
    rho1.clean_frame();
    rho2.clean_frame();
    rho_t.clean_frame();

    for(i=0; i<traj.target_leaflet.size(); i++) //loop over the target leaflet atoms
    {
        //get the first and last atom of the current lipid
        int min = traj.t_lip_start(i);
        int max = traj.t_lip_end(i);

        //jump to the next lipid
        i = traj.next_target_lipid(i);

        for(j=0; j<param_1.main_size_y(); j++) //loop over lipid types 1
        {
            if(strcmp(traj.res_name[min].c_str(), param_1.param_main_s[j][0].c_str() ) == 0) //lipid 1 type is correct
            {
                for(k=min; k<=max; k++) //loop over current residue atoms
                {
                    for(l=0; l<param_1.sec_size_y(j); l++) //loop over target atoms
                    {
                        if(strcmp(traj.atom_name[k].c_str(), param_1.param_sec_s[j][l][0].c_str() ) == 0) //atom is a target atom
                        {
                            double rx = traj.r[k][0];
                            double ry = traj.r[k][1];
                            double rz = traj.r[k][2];

                            rho1.stamp(rx,ry,rz,p.radius,1.0);
                            rho_t.stamp(rx,ry,rz,p.radius,1.0);
                        }
                    }
                }
            }
        }

        for(j=0; j<param_2.main_size_y(); j++) //loop over lipid types 2
        {
            if(strcmp(traj.res_name[min].c_str(), param_2.param_main_s[j][0].c_str() ) == 0) //lipid 2 type is correct
            {
                for(k=min; k<=max; k++) //loop over current residue atoms
                {
                    for(l=0; l<param_2.sec_size_y(j); l++) //loop over target atoms
                    {
                        if(strcmp(traj.atom_name[k].c_str(), param_2.param_sec_s[j][l][0].c_str() ) == 0) //atom is a target atom
                        {
                            double rx = traj.r[k][0];
                            double ry = traj.r[k][1];
                            double rz = traj.r[k][2];

                            rho2.stamp(rx,ry,rz,p.radius,1.0);
                            rho_t.stamp(rx,ry,rz,p.radius,1.0);

                        }
                    }
                }
            }
        }
    }

    //get the average for the current frame
    rho1.norm_frame();
    rho2.norm_frame();
    rho_t.norm_frame();

    //now we add the single frame data to the long term sums 
    rho1.add_frame();
    rho2.add_frame();
    rho_t.add_frame();

    //now we print the single frame rho data
    if(p.print_stride > 0)
    {
        if(traj.get_frame_global()%p.print_stride == 0)
        {
            rho1.write_frame(traj.get_frame_global(),p.ex_val);
            rho2.write_frame(traj.get_frame_global(),p.ex_val);
            rho_t.write_frame(traj.get_frame_global(),p.ex_val);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects rho for the 2 lipid types and computes the percent enrichmnet                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,Grid_3d &rho1,Grid_3d &rho2,Grid_3d &rho_t,Grid_3d &enrich)
{
    int i = 0;                              //standard variable used in loops
    int j = 0;                              //standard variable used in loops
    int k = 0;                              //standard variable used in loops

    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nfinalizing analysis. This requires communication of the grid and could take some time. \n");
    }

    //collect rho1, rho2 and rho_t from all ranks
    rho1.collect_grid();
    rho2.collect_grid();
    rho_t.collect_grid();

    //write the enrichment and rho1/2 to file
    if(s.world_rank == 0)
    {
        //copy rho_t for excluding data
        rho1.copy_rho(rho_t);
        rho2.copy_rho(rho_t);

        //exclude data
        rho1.exclude_data(p.cutoff,0);
        rho2.exclude_data(p.cutoff,0);
        rho_t.exclude_data(p.cutoff,1);

        //now write the lipid density to output file
        rho1.write_grid(p.ex_val);
        rho2.write_grid(p.ex_val);
        rho_t.write_grid(p.ex_val);

        //compute enrichment
        for(i=0; i<enrich.num_z(); i++) //loop over z-dimension
        {
            for(j=0; j<enrich.num_y(); j++) //loop over y-dimension
            {
                for(k=0; k<enrich.num_x(); k++) //loop over x-dimension 
                {
                    if(rho2.grid[i][j][k] > 0) //not division by zero
                    {
                        double ratio         = rho1.grid[i][j][k]/rho2.grid[i][j][k];
                        enrich.grid[i][j][k] = ((ratio - p.bulk)/p.bulk)*100;
                    }
                    else //division by zero
                    {
                        enrich.grid[i][j][k] = 999999999999;
                    }
                }
            }
        }
        enrich.copy_rho(rho_t);
        enrich.exclude_data(p.cutoff,0);
        enrich.write_grid(p.ex_val);
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
    s.program_name = "3D Enrichment";

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
    add_argument_mpi_s(argc,argv,"-crd_1",  p.param_1_file_name,          "Selection card for lipids A (crd)",                           s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-crd_2",  p.param_2_file_name,          "Selection card for lipids B (crd)",                           s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-enrich", p.enrich_file_name,           "Output file with spatially resolved enrichment factor (dx)",  s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,           "PDB file with sorted leaflets (pdb)",                         s.world_rank, s.cl_tags, &p.b_lf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,  "File with additional leaflet finder parameters (prm)",        s.world_rank, s.cl_tags, &p.b_lf_param,0);
    add_argument_mpi_d(argc,argv,"-APS",    &p.APS,                       "Area per grid box (nm^2)",                                    s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-r",      &p.radius,                    "Radius of the target atom (nm)",                              s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-cutoff", &p.cutoff,                    "Cutoff for excluding data (chi)",                             s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-bx",     &p.box_x,                     "Grid x dimension (nm)",                                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-by",     &p.box_y,                     "Grid y dimension (nm)",                                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-bz",     &p.box_z,                     "Grid z dimension (nm)",                                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-leaf",   &p.leaflet,                   "Which leaflet? (0:both 1:upper 2:lower)",                     s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-ps",     &p.print_stride,              "Print single frame stamping data every nth frame (0:off)",    s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-ex_val", &p.ex_val,                    "Set excluded lattice points to this value",                   s.world_rank, s.cl_tags, nullptr,      0);
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
    check_extension_mpi(s.world_rank,"-crd_1",p.param_1_file_name,".crd");
    check_extension_mpi(s.world_rank,"-crd_2",p.param_2_file_name,".crd");
    check_extension_mpi(s.world_rank,"-enrich",p.enrich_file_name,".dx");

    if(p.b_lf_pdb == 1)
    {
        check_extension_mpi(s.world_rank,"-lf_pdb",p.lf_pdb_file_name,".pdb");
    }
    if(p.b_lf_param == 1)
    {
        check_extension_mpi(s.world_rank,"-lf_prm",p.leaflet_finder_param_name,".prm");
    }

    //create parameter files
    Param param_1;
    Param param_2;

    //read parameter files
    param_1.get_param(p.param_1_file_name,2,1,1);
    param_2.get_param(p.param_2_file_name,2,1,1);

    //check the integrity of the parameter files
    if(param_1.check_file() == 0 || param_2.check_file() == 0) //bad files. kill program 
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
    traj.get_lipid_selection_stats(param_1.get_column_s(0),"-crd_1");
    traj.get_lipid_selection_stats(param_2.get_column_s(0),"-crd_2");

    //create a grid to hold enrichment etc.
    Grid_3d rho1;
    Grid_3d rho2;
    Grid_3d rho_t;
    Grid_3d enrich;

    //get the grid dimensions
    rho1.get_dim(p.box_x,p.box_y,p.box_z,traj.ibox,p.APS);
    if(s.world_rank == 0) //attempt to estimate requirements befor a crash
    {
        float mem_d = 4.0*(float)rho1.num_x()*(float)rho1.num_y()*(float)rho1.num_z()*8.0;
        float mem_i = 2.0*(float)rho1.num_x()*(float)rho1.num_y()*(float)rho1.num_z()*4.0;
        float mem_t = 4.0*(mem_d + mem_i);
        printf("Estimated memory to hold the grid: %f (MB) \n\n",mem_t/1000000.0);
    }
    rho2.get_dim(p.box_x,p.box_y,p.box_z,traj.ibox,p.APS);
    rho_t.get_dim(p.box_x,p.box_y,p.box_z,traj.ibox,p.APS);
    enrich.get_dim(p.box_x,p.box_y,p.box_z,traj.ibox,p.APS);

    //create file names for rho data
    string rho1_file_name  = add_tag(p.enrich_file_name,"_rho_A");
    string rho2_file_name  = add_tag(p.enrich_file_name,"_rho_B");
    string rho_t_file_name = add_tag(p.enrich_file_name,"_rho_t");

    //set the output file name for grid
    rho1.set_output(rho1_file_name);
    rho2.set_output(rho2_file_name);
    rho_t.set_output(rho_t_file_name);
    enrich.set_output(p.enrich_file_name);

    //print info about the grid
    enrich.print_dim();

    //count lipid in the selected leaflet
    p.num_lip_t1 = traj.count_target_lipids_type(param_1.get_column_s(0));
    p.num_lip_t2 = traj.count_target_lipids_type(param_2.get_column_s(0));
    p.bulk = (double)p.num_lip_t1/(double)p.num_lip_t2;

    if(s.world_rank == 0)
    {
        printf("num_lip_t1 %d num_lip_t2 %d bulk %f \n\n",p.num_lip_t1,p.num_lip_t2,p.bulk);
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

        enrichment(traj,s,p,param_1,param_2,rho1,rho2,rho_t);

        traj.set_beta_lf();

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect lipid density from mpi processes and compute enrichment
    perf.log_time(finalize_analysis(traj,s,p,rho1,rho2,rho_t,enrich),"Fin Ana");

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

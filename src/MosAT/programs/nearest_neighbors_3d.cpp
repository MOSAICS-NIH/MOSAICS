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
#include "MosAT/program_variables/pv_nearest_neighbors_3d.h"//This has the variables specific to the analysis program
#include "headers/array.h"                                   //This has routines used for working with arrays
#include "headers/performance.h"                             //This has a class for logging performance data
#include "headers/index.h"                                   //This has a class for working with index files
#include "headers/traj.h"                                    //This has a class for working with the trajectory
#include "headers/leaflet_finder.h"                          //This has routines used to find leaflets in membrane simulations
#include "headers/protein_finder.h"                          //This has routines used to find protein atoms
#include "headers/sol_finder.h"                              //This has routines used to find the solvent
#include "headers/grid.h"                                    //This has routines used for working with a grid
#include "headers/grid_3d.h"                                 //This has routines used for working with a 3d grid
#include "headers/protein.h"                                 //This has routines used for working with protein data
#include "headers/force_serial.h"                            //This has routines used for forcing the code to run on a single mpi process
#include "headers/param.h"                                   //This has routines used for reading complex parameter data

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the number of local lipids and adds it to the grid.                                //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void neighbors(Trajectory &traj,system_variables &s,program_variables &p,Param &param_1,Param &param_2,Grid_3d &nbr)
{
    int i           = 0;                      //standard variable used in loops
    int j           = 0;                      //standard variable used in loops
    int k           = 0;                      //standard variable used in loops
    int l           = 0;                      //standard variable used in loops
    double distance = 0;                      //distance between centers

    //clear the current frame grids
    nbr.clean_frame();

    for(i=0; i<traj.target_leaflet.size(); i++) //loop over target leaflet atoms
    {
        //jump to the next lipid
        i = traj.next_target_lipid(i);

        //get the first and last atom of the current lipid
        int min_1 = traj.t_lip_start(i);                      
        int max_1 = traj.t_lip_end(i);

        for(j=0; j<param_1.main_size_y(); j++) //loop over lipid types 
        {
            if(strcmp(traj.res_name[min_1].c_str(), param_1.param_main_s[j][0].c_str()) == 0) //lipid type is correct
            {
                double neighbors = 0;                                                               //how many centers are close to the lipid

                //extract the atoms used for the center
                sv1d target_atoms_1 = param_1.get_column_sec_s(j,0);

                //compute the center for the lipid
                dv1d r_center_1 = traj.center(target_atoms_1,min_1,max_1);

                for(k=0; k<traj.target_leaflet.size(); k++) //loop over target leaflet atoms
                {
                    //jump to the next lipid
                    k = traj.next_target_lipid(k);

                    //get the first and last atom of the current lipid
                    int min_2 = traj.t_lip_start(k);                                      //first atom of the current lipid
                    int max_2 = traj.t_lip_end(k);                                        //last atom of the current lipid

                    if(traj.res_nr[min_1] != traj.res_nr[min_2]) //dont count atoms on same lipid
                    {
                        for(l=0; l<param_2.main_size_y(); l++) //loop over lipid types 2 
                        {
                            if(strcmp(traj.res_name[min_2].c_str(), param_2.param_main_s[l][0].c_str()) == 0) //lipid type is correct
                            {
                                //extract the atoms used for the center
                                sv1d target_atoms_2 = param_2.get_column_sec_s(l,0);

                                //compute the center for the lipid
                                dv1d r_center_2 = traj.center(target_atoms_2,min_2,max_2);

                                //compute distance between centers
                                double dif_x = r_center_2[0] - r_center_1[0];
                                double dif_y = r_center_2[1] - r_center_1[1];
                                double dif_z = r_center_2[2] - r_center_1[2];

                                distance = sqrt(dif_x*dif_x + dif_y*dif_y + dif_z*dif_z);

                                if(distance < p.local_rad)
                                {
                                    neighbors++;
                                }
                            }
                        }
                    }
                }

                //find the mapping atoms and add neighbors to the grid
                for(l=min_1; l<=max_1; l++) //loop over current residue
                {
                    if(strcmp(traj.atom_name[l].c_str(), param_1.param_main_s[j][1].c_str()) == 0) //atom is mapping atom 1
                    {
                        double rx = traj.r[l][0];   //mapping atom x-coord
                        double ry = traj.r[l][1];   //mapping atom y-coord
                        double rz = traj.r[l][2];   //mapping atom y-coord

                        nbr.stamp(rx,ry,rz,p.radius,neighbors);
                    }
                    if(strcmp(traj.atom_name[l].c_str(), param_1.param_main_s[j][2].c_str()) == 0) //atom is mapping atom 2
                    {
                        double rx = traj.r[l][0];   //mapping atom x-coord
                        double ry = traj.r[l][1];   //mapping atom y-coord
                        double rz = traj.r[l][2];   //mapping atom y-coord

                        nbr.stamp(rx,ry,rz,p.radius,neighbors);
                    }
                }
            }
        }
    }

    //get the average for the current frame
    nbr.norm_frame();

    //add the current frame grid to long term sum
    nbr.add_frame();

    //now we print the single frame local density data for stdev computation
    if(p.b_stdev == 1)
    {
        nbr.write_frame(traj.get_frame_global(),p.ex_val);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Collect densities and compute average                                                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,Grid_3d &nbr)
{
    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nfinalizing analysis. This requires communication of the grid and could take some time. \n");
    }

    //now we collect nbr and rho from all nodes
    nbr.collect_grid();

    //normalize nbr and write the <nbr> and rho to file
    if(s.world_rank == 0)
    {
        nbr.normalize();

        nbr.exclude_data(p.cutoff,1);

        nbr.write_grid(p.ex_val);
        nbr.write_rho(p.ex_val);
    }

    //compute standard deviation
    nbr.get_stdev(p.b_stdev,p.b_clean,traj,p.ex_val);

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
    s.program_name = "Nearest Neighbors 3d";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                        s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                        s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                       s.world_rank, s.cl_tags, &p.b_print,   0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                         s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                                    s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                     s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                       s.world_rank, s.cl_tags, &p.b_lsq,     0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                         s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",         s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-crd_1",  p.param_1_file_name,          "Selection card 1 (crd)",                                            s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-crd_2",  p.param_2_file_name,          "Selection card 2 (crd)",                                            s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-nbrs",   p.nbrs_file_name,             "Output grid with spatially resolved nearest neighbors count (dx)",  s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,           "PDB file with sorted leaflets (pdb)",                               s.world_rank, s.cl_tags, &p.b_lf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,  "File with additional leaflet finder parameters (prm)",              s.world_rank, s.cl_tags, &p.b_lf_param,0);
    add_argument_mpi_d(argc,argv,"-l_rad",  &p.local_rad,                 "Radius of sphere for counting neighbors (nm)",                      s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-APS",    &p.APS,                       "Area per grid box (nm^2)",                                          s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-r",      &p.radius,                    "Radius of the target atom (nm)",                                    s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-cutoff", &p.cutoff,                    "Cutoff for excluding data (chi)",                                   s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-bx",     &p.box_x,                     "Grid x dimension (nm)",                                             s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-by",     &p.box_y,                     "Grid y dimension (nm)",                                             s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-bz",     &p.box_z,                     "Grid z dimension (nm)",                                             s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-leaf",   &p.leaflet,                   "Which leaflet? (0:both 1:upper 2:lower)",                           s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-ex_val", &p.ex_val,                    "Set excluded lattice points to this value",                         s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-stdev",  &p.b_stdev,                   "Compute the STDEV and STEM? (0:no 1:yes)",                          s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-clean",  &p.b_clean,                   "Remove single frame files? (0:no 1:yes)",                           s.world_rank, s.cl_tags, nullptr,      0);
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
    check_extension_mpi(s.world_rank,"-nbrs",p.nbrs_file_name,".dx");

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
    param_1.get_param(p.param_1_file_name,4,3,1);
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

    //create a grid to hold neighbors count
    Grid_3d nbr;

    //get the grid dimensions
    nbr.get_dim(p.box_x,p.box_y,p.box_z,traj.ibox,p.APS);
    if(s.world_rank == 0) //attempt to estimate requirements befor a crash
    {
        float mem_d = 4.0*(float)nbr.num_x()*(float)nbr.num_y()*(float)nbr.num_z()*8.0;
        float mem_i = 2.0*(float)nbr.num_x()*(float)nbr.num_y()*(float)nbr.num_z()*4.0;
        float mem_t = 1.0*(mem_d + mem_i);
        printf("Estimated memory to hold the grid: %f (MB) \n\n",mem_t/1000000.0);
    }

    //set the output file name for grid
    nbr.set_output(p.nbrs_file_name);

    //print info about the grid
    nbr.print_dim();
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

        neighbors(traj,s,p,param_1,param_2,nbr);

        traj.set_beta_lf();

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect densities from mpi processes and compute average
    perf.log_time(finalize_analysis(traj,s,p,nbr),"Fin Ana");

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

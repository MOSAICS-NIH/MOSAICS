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

#include "xdr/include/xdrfile_xtc.h"                               //used to read xtc files 
#include "xdr/include/xdr_seek.h"                                  //used to get and set the file position in xtc and trr files
#include "xdr/include/xdrfile_trr.h"                               //used to read trr files
#include "xdr/include/xdrfile.h"                                   //used to read xtc and trr files
#include "xdr/include/trr_header.h"                                //used to read the header info of trr files
#include "headers/multi_dim_vec.h"                                 //This defines multidimensional vectors
#include "headers/switch.h"                                        //This defines a switch (on, off)
#include "headers/file_reader.h"                                   //This has basic routines for reading text files
#include "headers/vector_mpi.h"                                    //This has routines for collecting vector data
#include "headers/mosat_routines.h"                               //This is where most of the functions called in main are located
#include "headers/file_naming.h"                                   //This has routines for added tags to an existing file name    
#include "headers/file_naming_mpi.h"                               //This has routines for added tags to an existing file name (mpi)
#include "headers/command_line_args_mpi.h"                         //This has routines for adding command line arguments
#include "MosAT/program_variables/pv_lipid_density_3d.h"          //This has the variables specific to the analysis program
#include "headers/array.h"                                         //This has routines used for working with arrays
#include "headers/performance.h"                                   //This has a class for logging performance data
#include "headers/index.h"                                         //This has a class for working with index files
#include "headers/traj.h"                                          //This has a class for working with the trajectory
#include "headers/leaflet_finder.h"                                //This has routines used to find leaflets in membrane simulations
#include "headers/protein_finder.h"                                //This has routines used to find protein atoms
#include "headers/sol_finder.h"                                    //This has routines used to find the solvent
#include "headers/grid.h"                                          //This has routines used for working with a grid
#include "headers/grid_3d.h"                                       //This has routines used for working with a grid
#include "headers/protein.h"                                       //This has routines used for working with protein data
#include "headers/force_serial.h"                                  //This has routines used for forcing the code to run on a single mpi process
#include "headers/param.h"                                         //This has routines used for reading complex parameter data

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function deposits density to the grid around target atoms                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void deposit_density(Trajectory &traj,system_variables &s,program_variables &p,Param &param,Grid_3d_i &rho)
{
    int    i        = 0;                                //standard variable used in loops
    int    j        = 0;                                //standard variable used in loops
    int    k        = 0;                                //standard variable used in loops
    int    l        = 0;                                //standard variable used in loops
    int    m        = 0;                                //standard variable used in loops

    for(i=0; i<traj.target_leaflet.size(); i++) //loop over target leaflet atoms
    {
        //jump to the next lipid
        i = traj.next_target_lipid(i);
        
        //get the first and last atom of the current lipid
        int min = traj.t_lip_start(i);
        int max = traj.t_lip_end(i);
        for(j=0; j<param.main_size_y(); j++) //loop over lipid types  
        {
            if(strcmp(traj.res_name[min].c_str(), param.param_main_s[j][0].c_str()) == 0) //lipid type is correct
            {
                for(k=min; k<=max; k++)//loop over lipid atoms
                {
                    for(l=0; l<param.sec_size_y(j); l++) //loop over target atoms
                    {
                        if(strcmp(traj.atom_name[k].c_str(), param.param_sec_s[j][l][0].c_str()) == 0) //mapping atom 
                        {
                            double hx = traj.r[k][0];
                            double hy = traj.r[k][1];
                            double hz = traj.r[k][2];

                            if(p.b_dist == 1)
                            {
                                for(m=0; m<traj.prot.size(); m++) //loop over protein atoms
                                {
                                    double dif_x = hx - traj.r[traj.prot[m]-1][0];
                                    double dif_y = hy - traj.r[traj.prot[m]-1][1];
                                    double dif_z = hz - traj.r[traj.prot[m]-1][2];

                                    double distance = sqrt(dif_x*dif_x + dif_y*dif_y + dif_z*dif_z);

                                    if(distance <= p.dist_cutoff)
                                    {
                                        rho.stamp(hx,hy,hz,p.radius,1);
                                        break;
                                    }
                                }
                            }
                            else 
                            {
                                rho.stamp(hx,hy,hz,p.radius,1);
                            }
                        }
                    }
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Collect distances and compute average.                                                                    //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,Grid_3d_i &rho,Grid_3d_i &nan)
{
    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nfinalizing analysis. This requires communication of the grid and could take some time. \n");
    }

    //collect rho and rho from all ranks
    rho.collect_grid();

    //normalize rho and write the density to file
    if(s.world_rank == 0)
    {
        rho.write_grid(nan,p.ex_val);
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
    s.program_name = "Lipid Density 3D";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                  s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                  s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                 s.world_rank, &p.b_print,   0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                   s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                              s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                               s.world_rank, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                 s.world_rank, &p.b_lsq,     0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                   s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",   s.world_rank, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-rho",    p.rho_file_name,              "Output file with the lipid density map (dx)",                 s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-crd",    p.param_file_name,            "Name of the selection card (crd)",                            s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,           "PDB file with sorted leaflets (pdb)",                         s.world_rank, &p.b_lf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-pf_pdb", p.pf_pdb_file_name,           "PDB file with selected protein (pdb)",                        s.world_rank, &p.b_pf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,  "File with additional leaflet finder parameters (prm)",        s.world_rank, &p.b_lf_param,0);
    add_argument_mpi_s(argc,argv,"-pf_prm", p.protein_finder_param_name,  "File with additional protein finder parameters (prm)",        s.world_rank, &p.b_pf_param,0);
    add_argument_mpi_d(argc,argv,"-APS",    &p.APS,                       "Area per grid box (nm^2)",                                    s.world_rank, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-r",      &p.radius,                    "Radius of the target atom (nm)",                              s.world_rank, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-bx",     &p.box_x,                     "Grid x dimension (nm)",                                       s.world_rank, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-by",     &p.box_y,                     "Grid y dimension (nm)",                                       s.world_rank, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-bz",     &p.box_z,                     "Grid z dimension (nm)",                                       s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-leaf",   &p.leaflet,                   "Which leaflet? (0:both 1:upper 2:lower)",                     s.world_rank, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-dist",   &p.dist_cutoff,               "How far from protein should lipids be counted? (nm)",         s.world_rank, &p.b_dist,    0);
    add_argument_mpi_i(argc,argv,"-ex_val", &p.ex_val,                    "Set excluded lattice points to this value",                   s.world_rank, nullptr,      0);
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
    check_extension_mpi(s.world_rank,"-crd",p.param_file_name,".crd");
    check_extension_mpi(s.world_rank,"-rho",p.rho_file_name,".dx");

    if(p.b_lf_pdb == 1)
    {
        check_extension_mpi(s.world_rank,"-lf_pdb",p.lf_pdb_file_name,".pdb");
    }
    if(p.b_lf_param == 1)
    {
        check_extension_mpi(s.world_rank,"-lf_prm",p.leaflet_finder_param_name,".prm");
    }
    if(p.b_pf_pdb == 1)
    {
        check_extension_mpi(s.world_rank,"-pf_pdb",p.pf_pdb_file_name,".pdb");
    }
    if(p.b_pf_param == 1)
    {
        check_extension_mpi(s.world_rank,"-pf_prm",p.protein_finder_param_name,".prm");
    }

    //create parameter files
    Param param;

    //read parameter files
    param.get_param(p.param_file_name,2,1,1);

    //check the integrity of the parameter files
    if(param.check_file() == 0) //bad files. kill program 
    {
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }

    //run leaflet/protein finder
    traj.get_leaflets(p.leaflet,p.leaflet_finder_param_name,p.b_lf_param);
    traj.get_protein(p.protein_finder_param_name,p.b_pf_param);

    //print a pdb with distinguished leaflets/protein
    traj.write_leaflets(p.lf_pdb_file_name,p.b_lf_pdb);
    traj.write_protein(p.pf_pdb_file_name,p.b_pf_pdb);

    //print info about the target leaflet
    traj.get_leaflet_stats();

    //print info about the lipid selection
    traj.get_lipid_selection_stats(param.get_column_s(0),"-crd");

    //print info about the protein
    traj.get_prot_stats();

    //create a grid to hold nearest neighbor distances
    Grid_3d_i rho;
    Grid_3d_i nan;

    //get the grid dimensions
    rho.get_dim(p.box_x,p.box_y,p.box_z,traj.ibox,p.APS);
    if(s.world_rank == 0) //attempt to estimate requirements befor a crash
    {
        float mem_t = 2.0*(float)rho.num_x()*(float)rho.num_y()*(float)rho.num_z()*4.0;
        printf("Estimated memory to hold the grid: %f (MB) \n\n",mem_t/1000000.0);
    }
    nan.get_dim(p.box_x,p.box_y,p.box_z,traj.ibox,p.APS);

    //set the output file name for grid
    rho.set_output(p.rho_file_name);

    //print info about the grid
    rho.print_dim();
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

        deposit_density(traj,s,p,param,rho);

        traj.set_beta_lf();

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect rho from mpi processes and compute average
    perf.log_time(finalize_analysis(traj,s,p,rho,nan),"Fin Ana");

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

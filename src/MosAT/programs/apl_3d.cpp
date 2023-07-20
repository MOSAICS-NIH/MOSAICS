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
#include "MosAT/program_variables/pv_apl_3d.h"              //This has the variables specific to the analysis program
#include "headers/performance.h"                             //This has a class for logging performance data
#include "headers/array.h"                                   //This has routines used for working with arrays
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
#include "headers/voronoi.h"                                 //This has routines used for computing voronoi diagrams
#include "headers/histo.h"                                   //This has routines used for making histograms

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the distances and adds them to the grid                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void area(Trajectory &traj,Trajectory &traj_v,system_variables &s,program_variables &p,Index &param1,Param &param2,Grid_3d &apl,dv1d &area_measurements)
{
    int    i         = 0;                      //standard variable used in loops
    int    j         = 0;                      //standard variable used in loops
    int    k         = 0;                      //standard variable used in loops
    int    l         = 0;                      //standard variable used in loops
    int dynamic      = 1;                      //use a dynamic box for voronoi diagrams?
    int pbc          = 1;                      //account for periodic boundary conditions in voronoi diagrams?
    double hx        = 0;                      //x coord of mapping atom
    double hy        = 0;                      //y coord of mapping atom
    
    //get the voronoi diagram
    Grid_i voronoi = voronoi_diagram(traj_v,apl.get_aps(),apl.num_x(),apl.num_y(),param2,p.c_dist,p.voro_stamp_rad,p.v_prot,dynamic,pbc); 

    //write the voronoi diagram to file
    if(p.b_voronoi == 1)
    {
        //set output filename
        string voronoi_file_name = chop_and_add_tag(p.apl_file_name,"_voro_" + to_string(traj.get_frame_global()) + ".dat"); 
        voronoi.set_output(voronoi_file_name,p.out_data_format);

        //create grid for excluding data
        Grid_i nan;
        nan.set_dim(voronoi.get_aps(),voronoi.num_x(),voronoi.num_y());

        //write voronoi diagram to file
        voronoi.write_grid(nan); 
    }

    //clear the current frame grids
    apl.clean_frame();

    for(i=0; i<traj.target_leaflet.size(); i++) //loop over the target leaflet atoms
    {
        //get the first and last atom of the current lipid
        int min = traj.t_lip_start(i);
        int max = traj.t_lip_end(i);

        //jump to the next lipid
        i = traj.next_target_lipid(i);

        for(j=0; j<param1.index_s.size(); j+=3) //loop over lipid types
        {
            if(strcmp(traj.res_name[min].c_str(), param1.index_s[j].c_str()) == 0) //residue is a target lipid
            {
                int count = 0;

                //compute the area of the lipid
                for(k=0; k<voronoi.num_x(); k++) //loop over x
                {
                    for(l=0; l<voronoi.num_y(); l++) //loop over y
                    {
                        if(voronoi.grid[l][k] == traj.res_nr[min])
                        {
                            count++;
                        }
                    }
                }
                double area = (double)count*voronoi.get_aps();

                //store measurement for histogram
                area_measurements.push_back(area);

                for(k=min; k<=max; k++) //loop over current residue atoms
                {
                    if(strcmp(traj.atom_name[k].c_str(), param1.index_s[j+1].c_str()) == 0 || strcmp(traj.atom_name[k].c_str(), param1.index_s[j+2].c_str()) == 0) //mapping atom 
                    {
                        double hx = traj.r[k][0];
                        double hy = traj.r[k][1];
                        double hz = traj.r[k][2];
 
                        apl.stamp(hx,hy,hz,p.radius,area);
                    }
                }
            }
        }
    }

    //get the average for the current frame
    apl.norm_frame();

    //add the current frame grid to long term sum
    apl.add_frame();

    //now we print the single frame apl data
    if(p.b_stdev == 1)
    {
        apl.write_frame(traj.get_frame_global(),p.ex_val);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Collects distances and computes the average                                                               //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,Grid_3d &apl,dv1d &area_measurements)
{
    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nfinalizing analysis. This requires communication of the grid and could take some time. \n");
    }

    //collect apl and rho from all ranks
    apl.collect_grid();

    //normalize apl and write the <apl> and rho to file
    if(s.world_rank == 0)
    {
        apl.normalize();
        
        apl.exclude_data(p.cutoff,1);
       
        apl.write_grid(p.ex_val);
        apl.write_rho(p.ex_val);
    }

    //compute standard deviation
    apl.get_stdev(p.b_stdev,p.b_clean,traj,p.ex_val);

    //collect area measurements
    collect_dv1d(s.world_size,s.world_rank,area_measurements);

    if(s.world_rank == 0)
    {
        printf("\nWorking on probability histogram.\n");
    
        //get filename for histogram
        string histo_file_name = chop_and_add_tag(p.apl_file_name,"_histo") + ".dat";
        
        //bin area measurements
        Histogram_d histo;
        histo.bin_data(area_measurements,p.bin_width);
        histo.write_histo(histo_file_name,"area (nm^2)");

        //print average and stdev
        printf("average: %f nm^2 \n",histo.get_avg(area_measurements));
        printf("  stdev: %f nm^2 \n",histo.get_stdev(area_measurements));
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
    s.program_name = "APL 3d";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name);

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                                   s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-traj_v", p.traj_voro_file_name,        "Input trajectory file for constructing voronoi diagrams (xtc, trr, pdb, gro)", s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                                   s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                                  s.world_rank, s.cl_tags, &p.b_print,   0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                                    s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                                               s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                                s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                                  s.world_rank, s.cl_tags, &p.b_lsq,     0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                                    s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",                    s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-apl",    p.apl_file_name,              "Output file with spatially resolved area per lipid (dx)",                      s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-crd_1",  p.param1_file_name,           "Selection card with target lipid types and mapping atoms (crd)",               s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-crd_2",  p.param2_file_name,           "Selection card with voronoi diagram lipid types (crd)",                        s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,           "PDB file with sorted leaflets (pdb)",                                          s.world_rank, s.cl_tags, &p.b_lf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,  "File with additional leaflet finder parameters (prm)",                         s.world_rank, s.cl_tags, &p.b_lf_param,0);
    add_argument_mpi_s(argc,argv,"-pf_pdb", p.pf_pdb_file_name,           "PDB file with protein marked (pdb)",                                           s.world_rank, s.cl_tags, &p.b_pf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-pf_prm", p.protein_finder_param_name,  "File with additional protein finder parameters (prm)",                         s.world_rank, s.cl_tags, &p.b_pf_param,0);
    add_argument_mpi_d(argc,argv,"-APS",    &p.APS,                       "Area per grid box (nm^2)",                                                     s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-r",      &p.radius,                    "Radius of the target atom (nm)",                                               s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-cutoff", &p.cutoff,                    "Cutoff for excluding data (chi)",                                              s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-bx",     &p.box_x,                     "Grid x dimension (nm)",                                                        s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-by",     &p.box_y,                     "Grid y dimension (nm)",                                                        s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-bz",     &p.box_z,                     "Grid z dimension (nm)",                                                        s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-leaf",   &p.leaflet,                   "Which leaflet? (0:both 1:upper 2:lower)",                                      s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-stdev",  &p.b_stdev,                   "Compute the STDEV and STEM? (0:no 1:yes)",                                     s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-clean",  &p.b_clean,                   "Remove single frame files? (0:no 1:yes)",                                      s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-c_dist", &p.c_dist,                    "Distance cutoff for counting protein atoms in voronoi diagram (nm)",           s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-voro",   &p.b_voronoi,                 "Write single frame voronoi files? (0:no 1:yes)",                               s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-voro_r", &p.voro_stamp_rad,            "Stamping radius used when computing Voronoi diagrams (nm)",                    s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-bin",    &p.bin_width,                 "Bin width for histogram (nm^2)",                                               s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-ex_val", &p.ex_val,                    "Set excluded lattice points to this value",                                    s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-odf",    &p.out_data_format,           "What is the output data format? (0:matrix 1:vector)",                          s.world_rank, s.cl_tags, nullptr,      0);
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

    //create a trajectory for voronoi diagrams
    Trajectory traj_v;

    //set trajectory parameters
    traj_v.set_block_parallel(on);
    traj_v.set_traj(p.traj_voro_file_name);
    traj_v.set_ref(p.ref_file_name);
    traj_v.set_traj_w(p.out_file_name,p.b_print);
    traj_v.set_lsq(p.lsq_index_file_name,p.b_lsq,p.lsq_dim,p.lsq_ref);
    traj_v.set_res(p.stride,p.start_frame,p.end_frame);

    //analyze the trajectory (log time spent) 
    perf.log_time(traj_v.build(),"Analyze Traj Voro");

    //compare trajectories 
    if(traj.get_frames() != traj_v.get_frames() || traj.atoms() != traj_v.atoms())
    {
        if(s.world_rank == 0)
        {
            printf("Mismatch between trajectories detected. \n");
            printf(" %40s %-10d \n","Number of trajectory frames: (-traj):",traj.get_frames());
            printf(" %40s %-10d \n","Number of trajectory frames: (-traj_v):",traj_v.get_frames());
            printf(" %40s %-10d \n","Number of atoms (-traj):",traj.atoms());
            printf(" %40s %-10d \n","Number of atoms (-traj_v):",traj_v.atoms());
        }
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //check file extensions                                                                                     
    check_extension_mpi(s.world_rank,"-crd_1",p.param1_file_name,".crd");
    check_extension_mpi(s.world_rank,"-crd_2",p.param2_file_name,".crd");
    check_extension_mpi(s.world_rank,"-apl",p.apl_file_name,".dx");

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
    Index param1;
    Param param2;

    //read parameter files
    param1.get_index(p.param1_file_name);
    param2.get_param(p.param2_file_name,2,1,1);

    //check the integrity of the parameter files
    if(param2.check_file() == 0) //bad files. kill program 
    {
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }

    //run leaflet/protein finder
    traj.get_leaflets(p.leaflet,p.leaflet_finder_param_name,p.b_lf_param);
    traj.get_protein(p.protein_finder_param_name,p.b_pf_param);
    traj_v.get_leaflets(p.leaflet,p.leaflet_finder_param_name,p.b_lf_param);
    traj_v.get_protein(p.protein_finder_param_name,p.b_pf_param);

    //print a pdb with distinguished leaflets/protein
    traj.write_leaflets(p.lf_pdb_file_name,p.b_lf_pdb);
    traj.write_protein(p.pf_pdb_file_name,p.b_pf_pdb);

    //print info about the target leaflet
    traj.get_leaflet_stats();

    //print info about the lipid selection
    traj.get_lipid_selection_stats(param1.get_column_s(3,0),"-crd_1");
    traj.get_lipid_selection_stats(param2.get_column_s(0),"-crd_2");

    //create a grid to hold the apl etc.
    Grid_3d   apl;

    //get the grid dimensions
    apl.get_dim(p.box_x,p.box_y,p.box_z,traj.ibox,p.APS);
    if(s.world_rank == 0) //attempt to estimate requirements befor a crash
    {
        float mem_d = 4.0*(float)apl.num_x()*(float)apl.num_y()*(float)apl.num_z()*8.0;
        float mem_i = 2.0*(float)apl.num_x()*(float)apl.num_y()*(float)apl.num_z()*4.0;
        float mem_t = 1.0*(mem_d + mem_i);
        printf("Estimated memory to hold the grid: %f (MB) \n\n",mem_t/1000000.0);
    }

    //set the output file name for grid
    apl.set_output(p.apl_file_name);

    //print info about the grid
    apl.print_dim();

    //create vector for storing area measurements
    dv1d area_measurements(0,0.0);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //print info about the worlk load distribution
    traj.workload();

    //print that analysis is beginning
    traj.report_progress();

    s.t = clock();
    traj_v.current_frame=0;
    //read read frames of the trajector and perform analysis
    for(traj.current_frame=0; traj.current_frame<traj.get_num_frames(); traj.current_frame++,traj_v.current_frame++)
    {
        traj.read_traj_frame();
        traj_v.read_traj_frame();

        traj.do_fit();

        area(traj,traj_v,s,p,param1,param2,apl,area_measurements);

        traj.set_beta_lf();

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //combine distances from mpi processes and get the average
    perf.log_time(finalize_analysis(traj,s,p,apl,area_measurements),"Fin Ana");

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

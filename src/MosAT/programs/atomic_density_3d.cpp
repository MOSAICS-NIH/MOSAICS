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
#include "MosAT/program_variables/pv_atomic_density_3d.h"          //This has the variables specific to the analysis program
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
// This function checks if the atom is in the shape                                                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int check_in_shape(Trajectory &traj,system_variables &s,program_variables &p,Index &shapes,double hx,double hy,double hz)
{
    int i          = 0;        //standard variable used in loops
    int b_in_shape = 1;        //set to 1 initially so if a shape is not provided it automatically passes

    if(p.b_shapes == 1)
    {
        b_in_shape = 0;

        for(i=0; i<shapes.index_s.size(); i+=5)
        {
            double cent_x = shapes.index_d[i];
            double cent_y = shapes.index_d[i+1];
            double cent_z = shapes.index_d[i+2];
            double rad    = shapes.index_d[i+3];
            double height = shapes.index_d[i+4];

            double dx = hx - cent_x;
            double dy = hy - cent_y;
            double dz = hz - cent_z;

            double distance = sqrt(dx*dx + dy*dy);

            if(distance <= rad)
            {
                if(hz > cent_z - 0.5*height && hz < cent_z + 0.5*height)
                {
                    b_in_shape = 1;
                }
            }
        }
    }
    return b_in_shape;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function deposits density to the grid around target atoms                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void deposit_density(Trajectory &traj,system_variables &s,program_variables &p,Index &target_atoms,Grid_3d_d &rho,
		     Index &shapes,iv1d &atom_count)
{
    int    i  = 0;                                //standard variable used in loops
    int    j  = 0;                                //standard variable used in loops
    int count = 0;                                //how many atoms stamped data in the target region

    for(i=0; i<target_atoms.index_s.size(); i++) //loop over target atoms
    {
        double hx = traj.r[target_atoms.index_i[i]-1][0];
        double hy = traj.r[target_atoms.index_i[i]-1][1];
        double hz = traj.r[target_atoms.index_i[i]-1][2];

        int b_in_shape = check_in_shape(traj,s,p,shapes,hx,hy,hz);

        if(b_in_shape == 1)
        {
            if(p.b_dist == 1)
            {
                for(j=0; j<traj.prot.size(); j++) //loop over protein atoms
                {
                    double dif_x = hx - traj.r[traj.prot[j]-1][0];
                    double dif_y = hy - traj.r[traj.prot[j]-1][1];
                    double dif_z = hz - traj.r[traj.prot[j]-1][2];
            
                    double distance = sqrt(dif_x*dif_x + dif_y*dif_y + dif_z*dif_z);
            
                    if(distance <= p.dist_cutoff)
                    {
                        count++;
                        rho.stamp(hx,hy,hz,p.radius,1.0);
                        break;
                    }
                }
            }
            else
            {
                count++;
                rho.stamp(hx,hy,hz,p.radius,1.0);
            }
        }
    }
    atom_count[traj.current_frame] = count;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Collect distances and compute average.                                                                    //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,Grid_3d_d &rho,Grid_3d_i &nan,iv1d &atom_count)
{
    int i  = 0;                                //standard variable used in loops

    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nfinalizing analysis. This requires communication of the grid and could take some time. \n");
    }

    //here we collect the atom counts from each core
    collect_iv1d(s.world_size,s.world_rank,atom_count);

    //collect rho and rho from all ranks
    rho.collect_grid();

    //normalize rho and write the density to file
    if(s.world_rank == 0)
    {
        rho.normalize_constant((double)traj.get_ef_frames()); 
        rho.write_grid(nan,p.ex_val);

        //write atom count data to file
        string atom_count_file_name = chop_and_add_tag(p.rho_file_name,"_atom_count.dat");
	FILE *atom_count_file = fopen(atom_count_file_name.c_str(),"w");
        for(i=0; i<atom_count.size(); i++) //loop over frames
        {
            fprintf(atom_count_file," %10d %10d \n",i,atom_count[i]);
        }
        fclose(atom_count_file);
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
    s.program_name = "Atomic Density 3D";

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
    add_argument_mpi_s(argc,argv,"-rho",    p.rho_file_name,              "Output file with the atomic density map (dx)",                s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-atoms",  p.target_atoms_file_name,     "Name of the index file with target atoms (ndx)",              s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-APS",    &p.APS,                       "Area per grid box (nm^2)",                                    s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-r",      &p.radius,                    "Radius of the target atom (nm)",                              s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-bx",     &p.box_x,                     "Grid x dimension (nm)",                                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-by",     &p.box_y,                     "Grid y dimension (nm)",                                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-bz",     &p.box_z,                     "Grid z dimension (nm)",                                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-dist",   &p.dist_cutoff,               "How far from protein should lipids be counted? (nm)",         s.world_rank, s.cl_tags, &p.b_dist,    0);
    add_argument_mpi_d(argc,argv,"-ex_val", &p.ex_val,                    "Set excluded lattice points to this value",                   s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-pf_pdb", p.pf_pdb_file_name,           "PDB file with selected protein (pdb)",                        s.world_rank, s.cl_tags, &p.b_pf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-pf_prm", p.protein_finder_param_name,  "File with additional protein finder parameters (prm)",        s.world_rank, s.cl_tags, &p.b_pf_param,0);
    add_argument_mpi_s(argc,argv,"-shape",  p.shape_file_name,            "Selection card (ndx) with shapes limiting the analysis",      s.world_rank, s.cl_tags, &p.b_shapes,  0);
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
    check_extension_mpi(s.world_rank,"-atoms",p.target_atoms_file_name,".ndx");
    check_extension_mpi(s.world_rank,"-rho",p.rho_file_name,".dx");

    if(p.b_pf_pdb == 1)
    {
        check_extension_mpi(s.world_rank,"-pf_pdb",p.pf_pdb_file_name,".pdb");
    }
    if(p.b_pf_param == 1)
    {
        check_extension_mpi(s.world_rank,"-pf_prm",p.protein_finder_param_name,".prm");
    }
    if(p.b_shapes == 1)
    {
        check_extension_mpi(s.world_rank,"-shape",p.shape_file_name,".crd");
    }

    //create index for atoms to include in the analysis
    Index target_atoms;
    Index shapes; 

    //read the index files
    target_atoms.get_index(p.target_atoms_file_name);
    if(p.b_shapes == 1)
    { 
        shapes.get_index(p.shape_file_name);
    }

    //run protein finder
    traj.get_protein(p.protein_finder_param_name,p.b_pf_param);

    //print a pdb with distinguished protein
    traj.write_protein(p.pf_pdb_file_name,p.b_pf_pdb);

    //print info about the protein
    traj.get_prot_stats();

    //create a grid to hold density data
    Grid_3d_d rho;
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

    //create vector to hold atoms in target region
    iv1d atom_count(traj.get_num_frames(),0);
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

        deposit_density(traj,s,p,target_atoms,rho,shapes,atom_count);

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect rho from mpi processes and compute average
    perf.log_time(finalize_analysis(traj,s,p,rho,nan,atom_count),"Fin Ana");

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

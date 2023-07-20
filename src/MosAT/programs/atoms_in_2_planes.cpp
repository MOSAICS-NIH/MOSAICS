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
#include "MosAT/program_variables/pv_atoms_in_2_planes.h"   //This has the variables specific to the analysis program
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
// This function searches for pairs of atoms sitting on top of each other                                    //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void atoms_in_a_plane(Trajectory &traj,system_variables &s,program_variables &p)
{
    
    int i = 0;                 //standard variable used in loops
    int j = 0;                 //standard variable used in loops
    int used[traj.atoms()];    //records if the atom is already included in a pair
    init_iarray(used,traj.atoms());
    
    for(i=0; i<traj.prot.size(); i++) //search for atoms in upper set
    {   
        for(j=0; j<traj.prot.size(); j++) //search for atoms in lower set
        {   
            double x1 = traj.r[traj.prot[i]-1][0];
            double y1 = traj.r[traj.prot[i]-1][1];
            double z1 = traj.r[traj.prot[i]-1][2];
            
            double x2 = traj.r[traj.prot[j]-1][0];
            double y2 = traj.r[traj.prot[j]-1][1];
            double z2 = traj.r[traj.prot[j]-1][2];
            
            double dif_x = x1 - x2;
            double dif_y = y1 - y2;
            double dif_z = z1 - z2;
            
            double delta_xy = sqrt(dif_x*dif_x + dif_y*dif_y);
            double delta_z  = sqrt(dif_z*dif_z);
           
            //check if atoms are compatable to make a pair 
            if(delta_z > p.cutoff_z && delta_xy < p.cutoff_xy && strcmp(traj.atom_name[traj.prot[i]-1].c_str(), p.atom_type.c_str()) == 0 && strcmp(traj.atom_name[traj.prot[j]-1].c_str(), p.atom_type.c_str()) == 0 && used[traj.prot[i]-1] == 0 && used[traj.prot[j]-1] == 0)
            {   
                used[traj.prot[i]-1] = 1;
                used[traj.prot[j]-1] = 1;
                if(traj.r[i][2] > traj.r[j][2])
                {   
                    printf("Upper Residue: %6d atom: %6d Lower Residue: %6d atom: %6d \n",traj.res_nr[traj.prot[i]-1],traj.atom_nr[traj.prot[i]-1],traj.res_nr[traj.prot[j]-1],traj.atom_nr[traj.prot[j]-1]);
                }
                else
                {   
                    printf("Upper Residue: %6d atom: %6d Lower Residue: %6d atom: %6d \n",traj.res_nr[traj.prot[j]-1],traj.atom_nr[traj.prot[j]-1],traj.res_nr[traj.prot[i]-1],traj.atom_nr[traj.prot[i]-1]);
                }
            }
        }
    }
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
    s.program_name = "Atoms in 2 PLanes";

    //force program to run in serial?
    enum Switch serial         = on;

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
    add_argument_mpi_d(argc,argv,"-z",      &p.cutoff_z,                  "How far away in z to still count atoms (nm)",                 s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-xy",     &p.cutoff_xy,                 "How far away in xy to still count atoms (nm)",                s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-pf_pdb", p.pf_pdb_file_name,           "PDB file with selected protein (B-factor) (pdb)",             s.world_rank, s.cl_tags, &p.b_pf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-pf_prm", p.protein_finder_param_name,  "File with additional protein finder parameters (prm)",        s.world_rank, s.cl_tags, &p.b_pf_param,0);
    add_argument_mpi_s(argc,argv,"-type",   p.atom_type,                  "Atom type when finding pairs (string)",                       s.world_rank, s.cl_tags, nullptr,      1);
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
    if(p.b_pf_pdb == 1)
    {
        check_extension_mpi(s.world_rank,"-pf_pdb",p.pf_pdb_file_name,".pdb");
    }
    if(p.b_pf_param == 1)
    {
        check_extension_mpi(s.world_rank,"-pf_prm",p.protein_finder_param_name,".prm");
    }

    //run protein finder
    traj.get_protein(p.protein_finder_param_name,p.b_pf_param);

    //write pdb with distinguished protein
    traj.write_protein(p.pf_pdb_file_name,p.b_pf_pdb);

    //print info about the protein
    traj.get_prot_stats();
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

        atoms_in_a_plane(traj,s,p);

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

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

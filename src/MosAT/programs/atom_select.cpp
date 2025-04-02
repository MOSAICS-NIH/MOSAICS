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

#include "xdr/include/xdrfile_xtc.h"                   //used to read xtc files 
#include "xdr/include/xdr_seek.h"                      //used to get and set the file position in xtc and trr files
#include "xdr/include/xdrfile_trr.h"                   //used to read trr files
#include "xdr/include/xdrfile.h"                       //used to read xtc and trr files
#include "xdr/include/trr_header.h"                    //used to read the header info of trr files
#include "headers/multi_dim_vec.h"                     //This defines multidimensional vectors
#include "headers/switch.h"                            //This defines a switch (on, off)
#include "headers/file_reader.h"                       //This has basic routines for reading text files
#include "headers/vector_mpi.h"                        //This has routines for collecting vector data
#include "headers/mosat_routines.h"                    //This is where most of the functions called in main are located
#include "headers/file_naming.h"                       //This has routines for added tags to an existing file name    
#include "headers/file_naming_mpi.h"                   //This has routines for added tags to an existing file name (mpi)
#include "headers/command_line_args_mpi.h"             //This has routines for adding command line arguments
#include "MosAT/program_variables/pv_atom_select.h"    //This has the variables specific to the analysis program
#include "headers/array.h"                             //This has routines used for working with arrays
#include "headers/performance.h"                       //This has a class for logging performance data
#include "headers/index.h"                             //This has a class for working with index files
#include "headers/traj.h"                              //This has a class for working with the trajectory
#include "headers/leaflet_finder.h"                    //This has routines used to find leaflets in membrane simulations
#include "headers/protein_finder.h"                    //This has routines used to find protein atoms
#include "headers/sol_finder.h"                        //This has routines used to find the solvent
#include "headers/grid.h"                              //This has routines used for working with a grid
#include "headers/protein.h"                           //This has routines used for working with protein data
#include "headers/force_serial.h"                      //This has routines used for forcing the code to run on a single mpi process
#include "headers/atom_select.h"                       //This has routines used for making atom selections using a selection text

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
    s.program_name = "Atom Select";

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
    add_argument_mpi_s(argc,argv,"-sel",    p.selection_text_file_name,   "Selection card with the selection text (sel)",                s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,           "PDB file with sorted leaflets (pdb)",                         s.world_rank, s.cl_tags, &p.b_lf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,  "File with additional leaflet finder parameters (prm)",        s.world_rank, s.cl_tags, &p.b_lf_param,0);
    add_argument_mpi_s(argc,argv,"-pf_pdb", p.pf_pdb_file_name,           "PDB file with selected protein (pdb)",                        s.world_rank, s.cl_tags, &p.b_pf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-pf_prm", p.protein_finder_param_name,  "File with additional protein finder parameters (prm)",        s.world_rank, s.cl_tags, &p.b_pf_param,0);
    add_argument_mpi_s(argc,argv,"-sf_pdb", p.sf_pdb_file_name,           "PDB file with selected sol (pdb)",                            s.world_rank, s.cl_tags, &p.b_sf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-sf_prm", p.solvent_finder_param_name,  "File with additional solvent finder parameters (prm)",        s.world_rank, s.cl_tags, &p.b_sf_param,0);
    add_argument_mpi_i(argc,argv,"-leaf",   &p.leaflet,                   "Which leaflet? (0:both 1:upper 2:lower)",                     s.world_rank, s.cl_tags, nullptr,      1);
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
    check_extension_mpi(s.world_rank,"-sel",p.selection_text_file_name,".sel");
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
    if(p.b_sf_pdb == 1)
    {
        check_extension_mpi(s.world_rank,"-sf_pdb",p.sf_pdb_file_name,".pdb");
    }
    if(p.b_sf_param == 1)
    {
        check_extension_mpi(s.world_rank,"-sf_prm",p.solvent_finder_param_name,".prm");
    }

    //run leaflet/proten/solvent finder
    traj.get_leaflets(p.leaflet,p.leaflet_finder_param_name,p.b_lf_param);
    traj.get_protein(p.protein_finder_param_name,p.b_pf_param);
    traj.get_solvent(p.solvent_finder_param_name,p.b_sf_param);

    //print a pdb with distinguished leaflets/protein/solvent
    traj.write_leaflets(p.lf_pdb_file_name,p.b_lf_pdb);
    traj.write_protein(p.pf_pdb_file_name,p.b_pf_pdb);
    traj.write_sol(p.sf_pdb_file_name,p.b_sf_pdb);

    //print info about the target leaflet
    traj.get_leaflet_stats();

    //print info about the protein
    traj.get_prot_stats();

    //print info about the water
    traj.get_sol_stats();

    //create object to tag the selected atoms
    iv1d custom_sel(traj.atoms(),0);

    //create a object to hold an atom selection
    Selection sel; 
  
    //select the atoms
    sel.get_selection(traj,p.selection_text_file_name);

    //generate pdb file name for highlighting the selection 
    string pdb_filename = chop_and_add_tag(p.selection_text_file_name,".pdb");

    //highlight the selected atoms
    sel.highlight_sel(traj,pdb_filename);

    //refine the selection
    custom_sel = sel.tag_atoms(traj);

    //generate a file name for index file with selected atoms
    string index_filename = chop_and_add_tag(p.selection_text_file_name,".ndx");

    //write an index with selected atoms
    int i = 0;  //standard variable used for loops
    FILE *index_file = fopen(index_filename.c_str(), "w");
    if(index_file == NULL)
    {
        printf("failure opening %s for writing. \n",index_filename.c_str());
        fflush(stdin);
    }
    else
    {
        for(i=0; i<traj.atoms(); i++) //loop over system atoms
        {
            if(custom_sel[i] == 1)
            {
                fprintf(index_file," %d \n",i+1);
            }
        }
        fclose(index_file);
    }

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

        traj.do_fit();

        //highlight selected atoms in traj for using -o option
        int i = 0;
        int j = 0;
        for(i=0; i<traj.atoms(); i++) //loop over system atoms
        {
            traj.beta[i] = 0.0;

            for(j=0; j<sel.sel.size(); j++) //loop over selection
            {
                if(traj.atom_nr[i] == sel.sel[j])
                {
                    traj.beta[i] = 1.0;
                }
            }
        }

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

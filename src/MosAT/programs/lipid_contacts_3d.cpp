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
#include "MosAT/program_variables/pv_lipid_contacts_3d.h"    //This has the variables specific to the analysis program
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
#include "headers/atom_select.h"                             //This has routines used for making atom selections using a selection text

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the number of lipid contacts and adds it to the grid                               //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void lip_cont(Trajectory &traj,system_variables &s,program_variables &p,Param &param,Grid_3d &lcont,iv1d &custom_sel)
{
    int    i        = 0;                      //standard variable used in loops
    int    j        = 0;                      //standard variable used in loops
    int    k        = 0;                      //standard variable used in loops
    int    l        = 0;                      //standard variable used in loops
    int    m        = 0;                      //standard variable used in loops
    int    n        = 0;                      //standard variable used in loops
    double distance = 0;                      //how far between atoms
    double hx       = 0;                      //head atom x-component 
    double hy       = 0;                      //head atom y-component
    double hz       = 0;                      //head atom z-component

    iv1d b_counted(traj.atoms(),0);           //keep track of whether a residue has been added or not

    //clear the current frame grids
    lcont.clean_frame();

    for(i=0; i<traj.target_leaflet.size(); i++) //loop over target leaflet atoms
    {
        //get the first and last atom of the current lipid
        int min = traj.t_lip_start(i);
        int max = traj.t_lip_end(i);

        //jump to the next lipid
        i = traj.next_target_lipid(i);

        for(j=0; j<param.main_size_y(); j++) //loop over lipid types
        {
            if(strcmp(traj.res_name[min].c_str(), param.param_main_s[j][0].c_str()) == 0) //lipid type is correct
            {
                int contacts = 0;

                //reset b_counted for each lipid
                if(p.b_mol == 1) 
                {
                    for(k=0; k<traj.atoms(); k++) //loop over system atoms
                    {
                        b_counted[k] = 0;
                    }
                }

                for(k=min; k<=max; k++) //loop over current residue atoms
                {
                    for(l=0; l<param.sec_size_y(j); l++) //loop over target atoms
                    {
                        if(strcmp(traj.atom_name[k].c_str(), param.param_sec_s[j][l][0].c_str()) == 0) //atom is a target atom
                        {
                            //get lip-prot contacts
                            if(p.count_prot == 1)
                            {
                                for(m=0; m<traj.prot.size(); m++) //loop over protein atoms
                                {
                                    if(b_counted[traj.prot[m]-1] == 0 || p.b_mol == 0) //check that molecule has not been added already 
                                    {
                                        double dif_x = traj.r[k][0] - traj.r[traj.prot[m]-1][0];
                                        double dif_y = traj.r[k][1] - traj.r[traj.prot[m]-1][1];
                                        double dif_z = traj.r[k][2] - traj.r[traj.prot[m]-1][2];

                                        distance = sqrt(dif_x*dif_x + dif_y*dif_y + dif_z*dif_z);

                                        if(distance < p.contact_cutoff)
                                        {
                                            contacts++;

                                            //flag that residue has been counted
                                            if(p.b_mol == 1) 
                                            {
                                                int min_1 = traj.res_start[traj.res_nr[traj.prot[m]-1]-1];
                                                int max_1 = traj.res_end[traj.res_nr[traj.prot[m]-1]-1];

                                                for(n=min_1; n<=max_1; n++) //loop over current residue
                                                {
                                                    b_counted[n] = 1;
                                                }
                                            }
                                        }
                                    }
                                }
                            }

                            //get lip-lip contacts
                            if(p.count_lip == 1)
                            {
                                for(m=0; m<traj.full_membrane.size(); m++) //loop over membrane atoms
                                {
                                    if(b_counted[traj.full_membrane[m]-1] == 0 || p.b_mol == 0) //check that molecule has not been added already 
                                    {
                                        if(traj.res_nr[traj.full_membrane[m]-1] != traj.res_nr[k]) //dont count contacts between the lipid and itself
                                        {
                                            double dif_x = traj.r[k][0] - traj.r[traj.full_membrane[m]-1][0];
                                            double dif_y = traj.r[k][1] - traj.r[traj.full_membrane[m]-1][1];
                                            double dif_z = traj.r[k][2] - traj.r[traj.full_membrane[m]-1][2];

                                            distance = sqrt(dif_x*dif_x + dif_y*dif_y + dif_z*dif_z);

                                            if(distance < p.contact_cutoff)
                                            {
                                                contacts++;

                                                //flag that residue has been counted
                                                if(p.b_mol == 1)
                                                {
                                                    int min_1 = traj.res_start[traj.res_nr[traj.full_membrane[m]-1]-1];
                                                    int max_1 = traj.res_end[traj.res_nr[traj.full_membrane[m]-1]-1];

                                                    for(n=min_1; n<=max_1; n++) //loop over current residue
                                                    {
                                                        b_counted[n] = 1;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }

                            //get lip_sol contacts
                            if(p.count_sol == 1)
                            {
                                for(m=0; m<traj.sol.size(); m++) //loop over sol atoms
                                {
                                    if(b_counted[traj.sol[m]-1] == 0 || p.b_mol == 0) //check that molecule has not been added already 
                                    {
                                        double dif_x = traj.r[k][0] - traj.r[traj.sol[m]-1][0];
                                        double dif_y = traj.r[k][1] - traj.r[traj.sol[m]-1][1];
                                        double dif_z = traj.r[k][2] - traj.r[traj.sol[m]-1][2];

                                        distance = sqrt(dif_x*dif_x + dif_y*dif_y + dif_z*dif_z);

                                        if(distance < p.contact_cutoff)
                                        {
                                            contacts++;

                                            //flag that residue has been counted
                                            if(p.b_mol == 1)
                                            {
                                                int min_1 = traj.res_start[traj.res_nr[traj.sol[m]-1]-1];
                                                int max_1 = traj.res_end[traj.res_nr[traj.sol[m]-1]-1];

                                                for(n=min_1; n<=max_1; n++) //loop over current residue
                                                {
                                                    b_counted[n] = 1;
                                                }
                                            }
                                        }
                                    }
                                }
                            }

                            //count contacts with a custom selection
                            if(p.b_sel_text == 1)
                            {
                                for(m=0; m<traj.atoms(); m++) //loop over system atoms
                                {
                                    if(custom_sel[m] == 1)
                                    {
                                        if(b_counted[m] == 0 || p.b_mol == 0) //check that molecule has not been added already 
                                        {
                                            double dif_x = traj.r[k][0] - traj.r[m][0];
                                            double dif_y = traj.r[k][1] - traj.r[m][1];
                                            double dif_z = traj.r[k][2] - traj.r[m][2];

                                            distance = sqrt(dif_x*dif_x + dif_y*dif_y + dif_z*dif_z);

                                            if(distance < p.contact_cutoff)
                                            {
                                                contacts++;

                                                //flag that residue has been counted
                                                if(p.b_mol == 1)
                                                {
                                                    int min_1 = traj.res_start[traj.res_nr[m]-1];
                                                    int max_1 = traj.res_end[traj.res_nr[m]-1];

                                                    for(n=min_1; n<=max_1; n++) //loop over current residue
                                                    {
                                                        b_counted[n] = 1;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                //add the number of contacts to the grid
                for(k=min; k<=max; k++) //loop over current residue atoms
                {
                    if(strcmp(traj.atom_name[k].c_str(), param.param_main_s[j][1].c_str()) == 0) //mapping atom 1
                    {
                        hx = traj.r[k][0]; 
                        hy = traj.r[k][1];
                        hz = traj.r[k][2];

                        lcont.stamp(hx,hy,hz,p.radius,contacts);
                    }
                    else if(strcmp(traj.atom_name[k].c_str(), param.param_main_s[j][2].c_str()) == 0) //mapping atom 2
                    {
                        hx = traj.r[k][0];
                        hy = traj.r[k][1];
                        hz = traj.r[k][2];

                        lcont.stamp(hx,hy,hz,p.radius,contacts);
                    }
                }
            }
        }
    }

    //get the average for the current frame
    lcont.norm_frame();

    //add the current frame grid to long term sum
    lcont.add_frame();

    //now we print the single frame lc data
    if(p.b_stdev == 1)
    {
        lcont.write_frame(traj.get_frame_global(),p.ex_val);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Collect the contacts and compute the average                                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,Grid_3d &lcont)
{
    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nFinalizing analysis. This requires communicating data over the grid and could take some time depending on the resolution. \n");
    }

    //collect lcont and rho from all ranks
    lcont.collect_grid();

    //normalize lcont and write the <lcont> and rho to file
    if(s.world_rank == 0)
    {
        lcont.normalize();

        lcont.exclude_data(p.cutoff,1);

        lcont.write_grid(p.ex_val);
        lcont.write_rho(p.ex_val);
    }

    //compute standard deviation
    lcont.get_stdev(p.b_stdev,p.b_clean,traj,p.ex_val);

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
    s.program_name = "Lipid Contacts 3D";

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
    add_argument_mpi_s(argc,argv,"-crd",    p.param_file_name,            "Selection card (crd)",                                        s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lc",     p.lc_file_name,               "Output file with spatially resolved contacts count (dx)",     s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,           "PDB file with sorted leaflets (pdb)",                         s.world_rank, s.cl_tags, &p.b_lf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,  "File with additional leaflet finder parameters (prm)",        s.world_rank, s.cl_tags, &p.b_lf_param,0);
    add_argument_mpi_s(argc,argv,"-pf_pdb", p.pf_pdb_file_name,           "PDB file with selected protein (pdb)",                        s.world_rank, s.cl_tags, &p.b_pf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-pf_prm", p.protein_finder_param_name,  "File with additional protein finder parameters (prm)",        s.world_rank, s.cl_tags, &p.b_pf_param,0);
    add_argument_mpi_s(argc,argv,"-sf_pdb", p.sf_pdb_file_name,           "PDB file with selected sol (pdb)",                            s.world_rank, s.cl_tags, &p.b_sf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-sf_prm", p.solvent_finder_param_name,  "File with additional solvent finder parameters (prm)",        s.world_rank, s.cl_tags, &p.b_sf_param,0);
    add_argument_mpi_d(argc,argv,"-APS",    &p.APS,                       "Area per grid box (nm^2)",                                    s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-r",      &p.radius,                    "Radius of the target atom (nm)",                              s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-cutoff", &p.cutoff,                    "Cutoff for excluding data (chi)",                             s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-bx",     &p.box_x,                     "Grid x dimension (nm)",                                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-by",     &p.box_y,                     "Grid y dimension (nm)",                                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-bz",     &p.box_z,                     "Grid z dimension (nm)",                                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-leaf",   &p.leaflet,                   "Which leaflet? (0:both 1:upper 2:lower)",                     s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-cdist",  &p.contact_cutoff,            "The contact cutoff distance (nm)",                            s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-lip",    &p.count_lip,                 "Count lipids? (0:no 1:yes)",                                  s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-prot",   &p.count_prot,                "Count protein? (0:no 1:yes)",                                 s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-sol",    &p.count_sol,                 "Count sol? (0:no 1:yes)",                                     s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-stdev",  &p.b_stdev,                   "Compute the STDEV and STEM? (0:no 1:yes)",                    s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-clean",  &p.b_clean,                   "Remove single frame files? (0:no 1:yes)",                     s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-mol",    &p.b_mol,                     "Count residues instead of contacts? (0:no 1:yes)",            s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-ex_val", &p.ex_val,                    "Set excluded lattice points to this value",                   s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-sel",    p.selection_text_file_name,   "Input file with the atom selection text (sel)",               s.world_rank, s.cl_tags, &p.b_sel_text,0);
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
    check_extension_mpi(s.world_rank,"-crd",p.param_file_name,".crd");
    check_extension_mpi(s.world_rank,"-lc",p.lc_file_name,".dx");

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
    if(p.b_sel_text == 1)
    {
        check_extension_mpi(s.world_rank,"-sel",p.selection_text_file_name,".sel");
    }

    //create parameter files
    Param param;

    //read parameter files
    param.get_param(p.param_file_name,4,3,1);

    //check the integrity of the parameter files
    if(param.check_file() == 0) //bad files. kill program 
    {
        MPI_Finalize();
        exit(EXIT_SUCCESS);
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

    //print info about the lipid selection
    traj.get_lipid_selection_stats(param.get_column_s(0),"-crd");

    //print info about the protein
    traj.get_prot_stats();

    //print info about the water
    traj.get_sol_stats();

    //create a grid to hold lipid contacts
    Grid_3d lcont;

    //get the grid dimensions
    lcont.get_dim(p.box_x,p.box_y,p.box_z,traj.ibox,p.APS);
    if(s.world_rank == 0) //attempt to estimate requirements befor a crash
    {
        float mem_d = 4.0*(float)lcont.num_x()*(float)lcont.num_y()*(float)lcont.num_z()*8.0;
        float mem_i = 2.0*(float)lcont.num_x()*(float)lcont.num_y()*(float)lcont.num_z()*4.0;
        float mem_t = 1.0*(mem_d + mem_i);
        printf("Estimated memory to hold the grid: %f (MB) \n\n",mem_t/1000000.0);
    }

    //set the output file name for grid
    lcont.set_output(p.lc_file_name);

    //print info about the grid
    lcont.print_dim();

    //create object to hold protein atom refinement
    iv1d custom_sel(traj.atoms(),0);

    if(p.b_sel_text == 1)
    {
        //create a object to hold an atom selection
        Selection this_sel;

        //select the atoms
        this_sel.get_selection(traj,p.selection_text_file_name);

        //generate pdb file name for highlighting the selection 
        string pdb_filename = chop_and_add_tag(p.selection_text_file_name,".pdb");

        //highlight the selected atoms
        this_sel.highlight_sel(traj,pdb_filename);

        //refine the selection
        custom_sel = this_sel.tag_atoms(traj);
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

        lip_cont(traj,s,p,param,lcont,custom_sel);

        traj.set_beta_lf();

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect contacts from mpi processes and compute the average
    perf.log_time(finalize_analysis(traj,s,p,lcont),"Fin Ana");

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

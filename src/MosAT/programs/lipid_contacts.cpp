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
#include "MosAT/program_variables/pv_lipid_contacts.h"       //This has the variables specific to the analysis program
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
#include "headers/param.h"                                   //This has routines used for reading complex parameter data
#include "headers/atom_select.h"                             //This has routines used for making atom selections using a selection text

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the number of lipid contacts and adds it to the grid                               //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void lip_cont(Trajectory &traj,system_variables &s,program_variables &p,Param &param,Grid &lcont,Grid &lmol,iv1d &custom_sel)
{
    int    i        = 0;                      //standard variable used in loops
    int    j        = 0;                      //standard variable used in loops
    int    k        = 0;                      //standard variable used in loops
    int    l        = 0;                      //standard variable used in loops
    int    m        = 0;                      //standard variable used in loops
    int    n        = 0;                      //standard variable used in loops
    double hx       = 0;                      //head atom x-component 
    double hy       = 0;                      //head atom y-component

    dv2d lipid_centers(0,dv1d(3,0.0));        //hold the center of each lipid molecule in the target leaflet
    dv2d protein_centers(0,dv1d(3,0.0));      //hold the center of each protein residue
    dv2d sol_centers(0,dv1d(3,0.0));          //hold the center of each solvent molecule
    dv2d full_mem_centers(0,dv1d(3,0.0));     //hold the center of each lipid molecule in the full mem
    dv2d system_centers(0,dv1d(3,0.0));       //hold the center of each residue in the system

    //clear the current frame grids
    lcont.clean_frame();
    lmol.clean_frame();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Compute the centers used for screening distances                                                          //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    lipid_centers = traj.get_centers_target_lf(); 
    if(p.count_lip == 1)
    {
        full_mem_centers = traj.get_centers_full_mem();     
    }
    if(p.count_prot == 1)
    {
        protein_centers = traj.get_centers_prot();
    }
    if(p.count_sol == 1)
    { 
        sol_centers = traj.get_centers_sol();
    }
    if(p.b_sel_text == 1)
    {
        system_centers = traj.get_centers_system();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Compute the number of contacts with the lipids                                                            //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int counter_a = 0;
    for(i=0; i<traj.target_leaflet.size(); i++) //loop over target leaflet atoms
    {
        //get the first and last atom of the current lipid
        int min_l = traj.t_lip_start(i);
        int max_l = traj.t_lip_end(i);

        //jump to the next lipid
        i = traj.next_target_lipid(i);

        for(j=0; j<param.main_size_y(); j++) //loop over lipid types
        {
            if(strcmp(traj.res_name[min_l].c_str(), param.param_main_s[j][0].c_str()) == 0) //lipid type is correct
            {
                int counter_b = 0;    //count residues as they are encountered
                int contacts  = 0;    //count contacts between current residue pair        
                int molecules = 0;    //count the number of molecules in contact with the lipids

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // Get contacts with protein                                                                                 //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                if(p.count_prot == 1)
                {
                    for(k=0; k<traj.prot.size(); k++) //loop over protein atoms
                    {
                        //get the first and last atom of the current lipid
                        int min_p = traj.p_res_start(k);
                        int max_p = traj.p_res_end(k);

                        //jump to the next lipid
                        k = traj.next_prot_res(k);

                        double center_dx = lipid_centers[counter_a][0] - protein_centers[counter_b][0];
                        double center_dy = lipid_centers[counter_a][1] - protein_centers[counter_b][1];
                        double center_dz = lipid_centers[counter_a][2] - protein_centers[counter_b][2];

                        double center_dist = sqrt(center_dx*center_dx + center_dy*center_dy + center_dz*center_dz);
     
                        if(center_dist < p.screen_dist) //count contacts
                        {
                            int these_contacts = 0;

                            for(l=min_l; l<=max_l; l++) //loop over current lipid atoms
                            {
                                for(m=min_p; m<=max_p; m++) //loop over current protein residue atoms
                                {
                                    double dx = traj.r[l][0] - traj.r[m][0]; 
                                    double dy = traj.r[l][1] - traj.r[m][1];
                                    double dz = traj.r[l][2] - traj.r[m][2];

                                    double dist = sqrt(dx*dx + dy*dy + dz*dz);
 
                                    if(dist < p.contact_cutoff)
                                    {
                                        these_contacts++;
                                    }
                                }
                            }
                            if(these_contacts > 0)
                            {
                                molecules++;
                            }
                            contacts = contacts + these_contacts; 
                        }
                        counter_b++;
                    }
                }
                counter_b = 0;

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // Get contacts with lipids                                                                                  //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                if(p.count_lip == 1)
                {
                    for(k=0; k<traj.full_membrane.size(); k++) //loop over full membrane lipid atoms
                    {
                        //get the first and last atom of the current lipid
                        int min_fm = traj.fm_lip_start(k);
                        int max_fm = traj.fm_lip_end(k);

                        //jump to the next lipid
                        k = traj.next_full_mem_lipid(k);
                      
                        if(traj.res_nr[min_fm] != traj.res_nr[min_l]) //dont count contacts with self
                        {
                            double center_dx = lipid_centers[counter_a][0] - full_mem_centers[counter_b][0];
                            double center_dy = lipid_centers[counter_a][1] - full_mem_centers[counter_b][1];
                            double center_dz = lipid_centers[counter_a][2] - full_mem_centers[counter_b][2];

                            double center_dist = sqrt(center_dx*center_dx + center_dy*center_dy + center_dz*center_dz);

                            if(center_dist < p.screen_dist) //count contacts
                            {
                                int these_contacts = 0;

                                for(l=min_l; l<=max_l; l++) //loop over current lipid atoms
                                {
                                    for(m=min_fm; m<=max_fm; m++) //loop over current lipid atoms
                                    {
                                        double dx = traj.r[l][0] - traj.r[m][0];
                                        double dy = traj.r[l][1] - traj.r[m][1];
                                        double dz = traj.r[l][2] - traj.r[m][2];

                                        double dist = sqrt(dx*dx + dy*dy + dz*dz);

                                        if(dist < p.contact_cutoff)
                                        {
                                            these_contacts++;
                                        }
                                    }
                                }
                                if(these_contacts > 0)
                                {
                                    molecules++;
                                }
                                contacts = contacts + these_contacts;
                            }
                        }
                        counter_b++;
                    }
                }
		counter_b = 0;

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // Get contacts with solvent                                                                                 //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                if(p.count_sol == 1)
                {
                    for(k=0; k<traj.sol.size(); k++) //loop over solvent atoms
                    {
                        //get the first and last atom of the current water
                        int min_sol = traj.sol_start(k);
                        int max_sol = traj.sol_end(k);

                        //jump to the next water
                        k = traj.next_water(k);

                        double center_dx = lipid_centers[counter_a][0] - sol_centers[counter_b][0];
                        double center_dy = lipid_centers[counter_a][1] - sol_centers[counter_b][1];
                        double center_dz = lipid_centers[counter_a][2] - sol_centers[counter_b][2];

                        double center_dist = sqrt(center_dx*center_dx + center_dy*center_dy + center_dz*center_dz);

                        if(center_dist < p.screen_dist) //count contacts
                        {
                            int these_contacts = 0;

                            for(l=min_l; l<=max_l; l++) //loop over current lipid atoms
                            {
                                for(m=min_sol; m<=max_sol; m++) //loop over current lipid atoms
                                {
                                    double dx = traj.r[l][0] - traj.r[m][0];
                                    double dy = traj.r[l][1] - traj.r[m][1];
                                    double dz = traj.r[l][2] - traj.r[m][2];

                                    double dist = sqrt(dx*dx + dy*dy + dz*dz);

                                    if(dist < p.contact_cutoff)
                                    {
                                        these_contacts++;
                                    }
                                }
                            }
                            if(these_contacts > 0)
                            {
                                molecules++;
                            }
                            contacts = contacts + these_contacts;
                        }
                        counter_b++;
                    }
                }
                counter_b = 0;


                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // Get contacts with custom atom selection                                                                   //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                if(p.b_sel_text == 1)
                {
                    for(k=0; k<traj.atoms(); k++) //loop over system atoms
                    {
                        //get the first and last atom of the current residue
                        int min_sys = traj.get_res_start(k);
                        int max_sys = traj.get_res_end(k);

                        //jump to the next residue
                        k = traj.next_residue(k);

                        int check_this_res = 0;
                        for(l=min_sys; l<=max_sys; l++) //loop over current residue atoms
                        {
                            if(custom_sel[l] == 1)
                            {
                                check_this_res = 1; 
                            }
                        }

                        if(check_this_res == 1)
                        {
                            double center_dx = lipid_centers[counter_a][0] - system_centers[counter_b][0];
                            double center_dy = lipid_centers[counter_a][1] - system_centers[counter_b][1];
                            double center_dz = lipid_centers[counter_a][2] - system_centers[counter_b][2];

                            double center_dist = sqrt(center_dx*center_dx + center_dy*center_dy + center_dz*center_dz);

                            if(center_dist < p.screen_dist) //count contacts
                            {
                                int these_contacts = 0;

                                for(l=min_l; l<=max_l; l++) //loop over current lipid atoms
                                {
                                    for(m=min_sys; m<=max_sys; m++) //loop over current residue atoms
                                    {
                                        if(custom_sel[m] == 1)
                                        {
                                            double dx = traj.r[l][0] - traj.r[m][0];
                                            double dy = traj.r[l][1] - traj.r[m][1];
                                            double dz = traj.r[l][2] - traj.r[m][2];

                                            double dist = sqrt(dx*dx + dy*dy + dz*dz);

                                            if(dist < p.contact_cutoff)
                                            {
                                                these_contacts++;
                                            }
                                        }
                                    }
                                }
                                if(these_contacts > 0)
                                {
                                    molecules++;
                                }
                                contacts = contacts + these_contacts;
                            }
                        }
                        counter_b++;
                    }
                }
                counter_b = 0;

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // Find mapping atoms and stamp data to the grid                                                             //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                for(k=min_l; k<=max_l; k++) //loop over current lipid atoms
                {
                    if(strcmp(traj.atom_name[k].c_str(), param.param_main_s[j][1].c_str()) == 0) //mapping atom 1
                    {
                        hx = traj.r[k][0];
                        hy = traj.r[k][1];

                        lcont.stamp(hx,hy,p.radius,contacts);
                        lmol.stamp(hx,hy,p.radius,molecules);
                    }
                    else if(strcmp(traj.atom_name[k].c_str(), param.param_main_s[j][2].c_str()) == 0) //mapping atom 2
                    {
                        hx = traj.r[k][0];
                        hy = traj.r[k][1];

                        lcont.stamp(hx,hy,p.radius,contacts);
                        lmol.stamp(hx,hy,p.radius,molecules);
                    }
                }
            }
        }
        counter_a++;
    }

    //get the average for the current frame
    lcont.norm_frame();
    lmol.norm_frame();

    //add the current frame grid to long term sum
    lcont.add_frame();
    lmol.add_frame();

    //now we print the single frame lc data
    if(p.b_stdev == 1)
    {
        lcont.write_frame(traj.get_frame_global());
        lmol.write_frame(traj.get_frame_global());
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Collect the contacts and compute the average                                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,Grid &lcont,Grid &lmol)
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
    lmol.collect_grid();

    //normalize lcont and write the <lcont> and rho to file
    if(s.world_rank == 0)
    {
        lcont.normalize();
        lmol.normalize();

        lcont.exclude_data(p.cutoff,1);
        lmol.exclude_data(p.cutoff,0);

        lcont.write_grid();
        lcont.write_rho();
	lmol.write_grid();
        lmol.write_rho();
    }

    //compute standard deviation
    lcont.get_stdev(p.b_stdev,p.b_clean,traj);
    lmol.get_stdev(p.b_stdev,p.b_clean,traj);

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
    s.program_name = "Lipid Contacts";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                                 s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                                 s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                                s.world_rank, s.cl_tags, &p.b_print,   0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                                  s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                                             s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                              s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                                s.world_rank, s.cl_tags, &p.b_lsq,     0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                                  s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",                  s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-crd",    p.param_file_name,            "Selection card (crd)",                                                       s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lc",     p.lc_file_name,               "Output file with spatially resolved contacts count (dat)",                   s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,           "PDB file with sorted leaflets (pdb)",                                        s.world_rank, s.cl_tags, &p.b_lf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,  "File with additional leaflet finder parameters (prm)",                       s.world_rank, s.cl_tags, &p.b_lf_param,0);
    add_argument_mpi_s(argc,argv,"-pf_pdb", p.pf_pdb_file_name,           "PDB file with selected protein (pdb) ",                                      s.world_rank, s.cl_tags, &p.b_pf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-pf_prm", p.protein_finder_param_name,  "File with additional protein finder parameters (prm)",                       s.world_rank, s.cl_tags, &p.b_pf_param,0);
    add_argument_mpi_s(argc,argv,"-sf_pdb", p.sf_pdb_file_name,           "PDB file with selected sol (pdb)",                                           s.world_rank, s.cl_tags, &p.b_sf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-sf_prm", p.solvent_finder_param_name,  "File with additional solvent finder parameters (prm)",                       s.world_rank, s.cl_tags, &p.b_sf_param,0);
    add_argument_mpi_d(argc,argv,"-APS",    &p.APS,                       "Area per grid box (nm^2)",                                                   s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-r",      &p.radius,                    "Radius of the target atom (nm)",                                             s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-cutoff", &p.cutoff,                    "Cutoff for excluding data (chi)",                                            s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-bx",     &p.box_x,                     "Grid x dimension (nm)",                                                      s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-by",     &p.box_y,                     "Grid y dimension (nm)",                                                      s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-leaf",   &p.leaflet,                   "Which leaflet? (0:both 1:upper 2:lower)",                                    s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-odf",    &p.out_data_format,           "What is the output data format? (0:matrix 1:vector)",                        s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-cdist",  &p.contact_cutoff,            "The contact cutoff distance (nm)",                                           s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-lip",    &p.count_lip,                 "Count lipids? (0:no 1:yes)",                                                 s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-prot",   &p.count_prot,                "Count protein? (0:no 1:yes)",                                                s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-sol",    &p.count_sol,                 "Count sol? (0:no 1:yes)",                                                    s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-stdev",  &p.b_stdev,                   "Compute the STDEV and STEM? (0:no 1:yes)",                                   s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-clean",  &p.b_clean,                   "Remove single frame files? (0:no 1:yes)",                                    s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-sel",    p.selection_text_file_name,   "Input file with the atom selection text (sel)",                              s.world_rank, s.cl_tags, &p.b_sel_text,0);
    add_argument_mpi_d(argc,argv,"-screen", &p.screen_dist,               "Screen residues whose centers are within this disatnce (nm) of each other",  s.world_rank, s.cl_tags, nullptr,      0);
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
    check_extension_mpi(s.world_rank,"-lc",p.lc_file_name,".dat");

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

    //create a grid to hold lipid contacts or the mols in contact with the lipid
    Grid lcont;
    Grid lmol;

    //get the grid dimensions
    lcont.get_dim(p.box_x,p.box_y,traj.ibox,p.APS);
    lmol.get_dim(p.box_x,p.box_y,traj.ibox,p.APS);

    //generate name of output grid file for mols
    string mol_file_name = add_tag(p.lc_file_name,"_mol");

    //set the output file name for grid
    lcont.set_output(p.lc_file_name,p.out_data_format);
    lmol.set_output(mol_file_name,p.out_data_format);

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

        lip_cont(traj,s,p,param,lcont,lmol,custom_sel);

        traj.set_beta_lf();

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect contacts from mpi processes and compute the average
    perf.log_time(finalize_analysis(traj,s,p,lcont,lmol),"Fin Ana");

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

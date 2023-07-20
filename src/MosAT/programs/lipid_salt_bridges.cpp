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
#include "MosAT/program_variables/pv_lipid_salt_bridges.h"  //This has the variables specific to the analysis program
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the number of lipid contacts and adds it to the grid                               //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void lip_salt(Trajectory &traj,system_variables &s,program_variables &p,Index &param,Grid &saltb,Param &lip_charges,Param &prot_charges)
{
    int    i        = 0;                      //standard variable used in loops
    int    j        = 0;                      //standard variable used in loops
    int    k        = 0;                      //standard variable used in loops
    int    l        = 0;                      //standard variable used in loops
    int    m        = 0;                      //standard variable used in loops
    int    tot      = 0;                      //total number of salt bridges for frame
    double hx       = 0;                      //head atom x-component 
    double hy       = 0;                      //head atom y-component

    //clear the current frame grids
    saltb.clean_frame();

    for(i=0; i<traj.target_leaflet.size(); i++) //loop over target leaflet atoms
    {
        //get the first and last atom of the current lipid
        int min = traj.t_lip_start(i);
        int max = traj.t_lip_end(i);

        //jump to the next lipid
        i = traj.next_target_lipid(i);

        for(j=0; j<param.index_s.size(); j+=3) //loop over lipid types
        {
            if(strcmp(traj.res_name[min].c_str(), param.index_s[j].c_str()) == 0) //lipid type is correct
            {
                int contacts = 0;

                for(k=0; k<lip_charges.main_size_y(); k++) //loop over lipid charges
                {
                    if(strcmp(traj.res_name[min].c_str(), lip_charges.param_main_s[k][0].c_str() ) == 0) //lipid type is correct
                    {
                        sv1d lip_charge_grp = lip_charges.get_column_sec_s(k,0);
                        dv1d lip_center     = traj.center(lip_charge_grp,min,max); 

                        for(l=0; l<traj.prot.size(); l++) //loop over protein atoms
                        {
                            //get the first and last atoms of the current residue
                            int min_p = traj.p_res_start(l);
                            int max_p = traj.p_res_end(l);
 
                            //jump to the next residue
                            l = traj.next_prot_res(l);

                            for(m=0; m<prot_charges.main_size_y(); m++) //loop over protein charges
                            {
                                if(strcmp(traj.res_name[min_p].c_str(), prot_charges.param_main_s[m][0].c_str() ) == 0) //protein residue is correct
                                {
                                    sv1d prot_charge_grp = prot_charges.get_column_sec_s(m,0);
                                    dv1d prot_center     = traj.center(prot_charge_grp,min_p,max_p);

                                    int l_charge = 0;
                                    if(strcmp(lip_charges.param_main_s[k][1].c_str(), "+" ) == 0) //plus charge
                                    {
                                        l_charge = 1;
                                    }
                                    else if(strcmp(lip_charges.param_main_s[k][1].c_str(), "-" ) == 0) //negative charge
                                    {
                                        l_charge = -1;
                                    }
                                    else 
                                    {
                                        if(s.world_rank == 0)
                                        {
                                            printf("acceptable lipid charges include + or - (found %s). Please check items in your selection card (%s). \n",lip_charges.param_main_s[k][1].c_str(),p.lipid_charge_file_name.c_str());
                                        }
                                        MPI_Finalize();
                                        exit(EXIT_SUCCESS);
                                    }

                                    int p_charge = 0;
                                    if(strcmp(prot_charges.param_main_s[m][1].c_str(), "+" ) == 0) //plus charge
                                    {
                                        p_charge = 1;
                                    }
                                    else if(strcmp(prot_charges.param_main_s[m][1].c_str(), "-" ) == 0) //negative charge
                                    {
                                        p_charge = -1;
                                    }
                                    else
                                    {
                                        if(s.world_rank == 0)
                                        {
                                            printf("acceptable protein charges include + or - (found %s). Please check items in your selection card (%s). \n",prot_charges.param_main_s[m][1].c_str(),p.protein_charge_file_name.c_str());
                                        }
                                        MPI_Finalize();
                                        exit(EXIT_SUCCESS);
                                    }

                                    if(l_charge == -1*p_charge)
                                    {
                                        double dist = traj.get_dist(lip_center,prot_center);

                                        if(dist < p.contact_cutoff)
                                        {
                                            contacts++;
                                            tot++;

                                            //check h-bonds in pymol
                                            if(p.b_test == 1)
                                            {
                                                printf("lipid %d %s residue %d %s distance %f \n",traj.res_nr[min],traj.res_name[min].c_str(),traj.res_nr[min_p],traj.res_name[min_p].c_str(),dist);
                                                printf("pseudoatom lip%d, pos=[%f,%f,%f]\n",tot,10*lip_center[0],10*lip_center[1],10*lip_center[2]);
                                                printf("pseudoatom prot%d, pos=[%f,%f,%f]\n",tot,10*prot_center[0],10*prot_center[1],10*prot_center[2]);
                                                printf("show spheres, lip%d \n",tot);
                                                printf("show spheres, prot%d \n",tot);
                                                if(l_charge == 1)
                                                {
                                                    printf("color blue, lip%d \n",tot);
                                                }
                                                else if(l_charge == -1)
                                                {
                                                    printf("color red, lip%d  \n",tot);
                                                }
                                                printf("color green, prot%d \n",tot);
                                                printf("show licorice, resi %d \n",traj.res_nr[min]);
                                                printf("show licorice, resi %d \n",traj.res_nr[min_p]);
                                                printf("distance dist%d, lip%d, prot%d, 10, mode=0 \n",tot,tot,tot);
                                                printf("sel pair%d, prot%d + lip%d \n",tot,tot,tot);
                                                printf("orient pair%d \n",tot);
                                                printf("\n");
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
                    if(strcmp(traj.atom_name[k].c_str(), param.index_s[j+1].c_str()) == 0) //mapping atom 1
                    {
                        hx = traj.r[k][0];
                        hy = traj.r[k][1];

                        saltb.stamp(hx,hy,p.radius,contacts);
                    }
                    else if(strcmp(traj.atom_name[k].c_str(), param.index_s[j+2].c_str()) == 0) //mapping atom 2
                    {
                        hx = traj.r[k][0];
                        hy = traj.r[k][1];

                        saltb.stamp(hx,hy,p.radius,contacts);
                    }
                }
            }
        }
    }

    //get the average for the current frame
    saltb.norm_frame();

    //add the current frame grid to long term sum
    saltb.add_frame();

    //now we print the single frame lpsb data
    if(p.b_stdev == 1)
    {
        saltb.write_frame(traj.get_frame_global());
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Collect the contacts and compute the average                                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,Grid &saltb)
{
    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nFinalizing analysis. This requires communicating data over the grid and could take some time depending on the resolution. \n");
    }

    //collect saltb and rho from all ranks
    saltb.collect_grid();

    //normalize saltb and write the <saltb> and rho to file
    if(s.world_rank == 0)
    {
        saltb.normalize();

        saltb.exclude_data(p.cutoff,1);

        saltb.write_grid();
        saltb.write_rho();
    }

    //compute standard deviation
    saltb.get_stdev(p.b_stdev,p.b_clean,traj);

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
    s.program_name = "Lipid Salt Bridges";

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
    add_argument_mpi_s(argc,argv,"-crd_1",  p.param_file_name,            "Selection card with Lipid types + mapping atoms (crd)",       s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-crd_2",  p.lipid_charge_file_name,     "Selection card with lipid charge groups (crd)",               s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-crd_3",  p.protein_charge_file_name,   "Selection card with protein charge groups (crd)",             s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lpsb",   p.lpsb_file_name,             "Output file with spatially resolved salt bridge count (dat)", s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,           "PDB file with sorted leaflets (pdb)",                         s.world_rank, s.cl_tags, &p.b_lf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,  "File with additional leaflet finder parameters (prm)",        s.world_rank, s.cl_tags, &p.b_lf_param,0);
    add_argument_mpi_s(argc,argv,"-pf_pdb", p.pf_pdb_file_name,           "PDB file with selected protein (pdb)",                        s.world_rank, s.cl_tags, &p.b_pf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-pf_prm", p.protein_finder_param_name,  "File with additional protein finder parameters (prm)",        s.world_rank, s.cl_tags, &p.b_pf_param,0);
    add_argument_mpi_d(argc,argv,"-APS",    &p.APS,                       "Area per grid box (nm^2)",                                    s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-r",      &p.radius,                    "Radius of the target atom (nm)",                              s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-cutoff", &p.cutoff,                    "Cutoff for excluding data (chi)",                             s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-bx",     &p.box_x,                     "Grid x dimension (nm)",                                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-by",     &p.box_y,                     "Grid y dimension (nm)",                                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-leaf",   &p.leaflet,                   "Which leaflet? (0:both 1:upper 2:lower)",                     s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-odf",    &p.out_data_format,           "What is the output data format? (0:matrix 1:vector)",         s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-cdist",  &p.contact_cutoff,            "The contact cutoff distance (nm)",                            s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-stdev",  &p.b_stdev,                   "Compute the STDEV and STEM? (0:no 1:yes)",                    s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-clean",  &p.b_clean,                   "Remove single frame files? (0:no 1:yes)",                     s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-test",   &p.b_test,                    "Print info for checking salt bridges? (0:no 1:yes)",          s.world_rank, s.cl_tags, nullptr,      0);
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
    check_extension_mpi(s.world_rank,"-crd_1",p.param_file_name,".crd");
    check_extension_mpi(s.world_rank,"-crd_2",p.lipid_charge_file_name,".crd");
    check_extension_mpi(s.world_rank,"-crd_3",p.protein_charge_file_name,".crd");
    check_extension_mpi(s.world_rank,"-lpsb",p.lpsb_file_name,".dat");

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

    //create index objects
    Index param;

    //read the index files
    param.get_index(p.param_file_name);

    //create parameter files
    Param lip_charges;
    Param prot_charges;

    //read parameter files
    lip_charges.get_param(p.lipid_charge_file_name,3,2,1);
    prot_charges.get_param(p.protein_charge_file_name,3,2,1);

    //run leaflet/proten finder
    traj.get_leaflets(p.leaflet,p.leaflet_finder_param_name,p.b_lf_param);
    traj.get_protein(p.protein_finder_param_name,p.b_pf_param);

    //print a pdb with distinguished leaflets/protein
    traj.write_leaflets(p.lf_pdb_file_name,p.b_lf_pdb);
    traj.write_protein(p.pf_pdb_file_name,p.b_pf_pdb);

    //print info about the target leaflet
    traj.get_leaflet_stats();

    //print info about the lipid selection
    traj.get_lipid_selection_stats(param.get_column_s(3,0),"-crd_1");

    //print info about the protein
    traj.get_prot_stats();

    //create a grid to hold lipid contacts
    Grid saltb;

    //get the grid dimensions
    saltb.get_dim(p.box_x,p.box_y,traj.ibox,p.APS);

    //set the output file name for grid
    saltb.set_output(p.lpsb_file_name,p.out_data_format);

    //print info about the grid
    saltb.print_dim();

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

        lip_salt(traj,s,p,param,saltb,lip_charges,prot_charges);

        traj.set_beta_lf();

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect contacts from mpi processes and compute the average
    perf.log_time(finalize_analysis(traj,s,p,saltb),"Fin Ana");

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

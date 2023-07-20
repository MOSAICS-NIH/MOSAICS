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
#include "MosAT/program_variables/pv_interdigitation.h"     //This has the variables specific to the analysis program
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
// This function computes the interdigitation (rank) and adds it to the grid                                 //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void interdigitate(Trajectory &traj,system_variables &s,program_variables &p,Index &param_1,Param &param_2,Grid &digi)
{
    //Strategy used here is to loop over the target leaflet and select lipids of the correct type. Then loop over the lipid atoms 
    //and find the tail atom and mapping atom. Once these atoms are identified, we loop over the opposing leaflet atoms and look for lipids of 
    //of the correct type. We then find contacts and look up the ranks. contacts. If a contact is found the rank is added such that the average 
    //rank, over the contacts, is computed. This average is then stamped to the grid around the mapping atom. This process is repeated for
    //the second tail and mapping atom.  

    int i = 0;                                //standard variable used in loops
    int j = 0;                                //standard variable used in loops
    int k = 0;                                //standard variable used in loops
    int l = 0;                                //standard variable used in loops
    int m = 0;                                //standard variable used in loops
    int n = 0;                                //standard variable used in loops
    int o = 0;                                //standard variable used in loops
    int q = 0;                                //standard variable used in loops

    //clear the current frame grids
    digi.clean_frame();

    for(i=0; i<traj.target_leaflet.size(); i++) //loop over the target membrane atoms
    {
        //get the first and last atom of the current lipid
        int min_1 = traj.t_lip_start(i);
        int max_1 = traj.t_lip_end(i);

        //jump to the next lipid
        i = traj.next_target_lipid(i);

        for(j=0; j<param_1.index_s.size(); j+=5) //loop over the lipid types
        {
            if(strcmp(traj.res_name[min_1].c_str(), param_1.index_s[j].c_str()) == 0) //lipid type is correct
            {
                for(k=min_1; k<=max_1; k++) //loop over current residue atoms
                {
                    //find the average for the first lipid tail
                    if(strcmp(traj.atom_name[k].c_str(), param_1.index_s[j+1].c_str()) == 0) //tail atom 1
                    {
                        for(l=min_1; l<=max_1; l++) //loop over current residue atoms
                        {
                            if(strcmp(traj.atom_name[l].c_str(), param_1.index_s[j+3].c_str()) == 0) //mapping atom 1
                            {
                                double hx = traj.r[l][0];
                                double hy = traj.r[l][1];

                                double contacts = 0.0;
                                double sum_rank = 0.0;

                                for(m=0; m<traj.opposing_leaflet.size(); m++) //loop over opposing leaflet atoms
                                {
                                    //get the first and last atom of the current lipid
                                    int min_2 = traj.o_lip_start(m);
                                    int max_2 = traj.o_lip_end(m);

                                    //jump to the next lipid
                                    m = traj.next_opposing_lipid(m);

                                    for(n=0; n<param_2.main_size_y(); n++) //loop over lipid types 
                                    {
                                        if(strcmp(traj.res_name[min_2].c_str(), param_2.param_main_s[n][0].c_str()) == 0) //lipid type is correct
                                        {
                                            for(o=min_2; o<=max_2; o++) //loop over current lipid atoms 
                                            {
                                                for(q=0; q<param_2.sec_size_y(n); q++) //loop over ranked atoms
                                                {
                                                    if(strcmp(traj.atom_name[o].c_str(), param_2.param_sec_s[n][q][0].c_str()) == 0) //atoms is a ranked atom
                                                    {
                                                        //compute the distance between the atoms
                                                        double dx = traj.r[k][0] - traj.r[o][0];
                                                        double dy = traj.r[k][1] - traj.r[o][1];
                                                        double dz = traj.r[k][2] - traj.r[o][2];

                                                        double distance = sqrt(dx*dx + dy*dy + dz*dz);

                                                        if(distance < p.contact_cutoff) //atoms make a contact
                                                        {
                                                            contacts = contacts + 1.0;
                                                            sum_rank = sum_rank + param_2.param_sec_d[n][q][1]; 
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }

                                if(contacts > 0)
                                {
                                    digi.stamp(hx,hy,p.radius,sum_rank/contacts);
                                }
                                else
                                {
                                    digi.stamp(hx,hy,p.radius,0.0);
                                } 
                            }
                        }
                    }

                    //find the average for the second lipid tail
                    if(strcmp(traj.atom_name[k].c_str(), param_1.index_s[j+2].c_str()) == 0) //tail atom 2
                    {
                        for(l=min_1; l<=max_1; l++) //loop over current residue atoms
                        {
                            if(strcmp(traj.atom_name[l].c_str(), param_1.index_s[j+4].c_str()) == 0) //mapping atom 2
                            {
                                double hx = traj.r[l][0];
                                double hy = traj.r[l][1];

                                double contacts = 0.0;
                                double sum_rank = 0.0;

                                for(m=0; m<traj.opposing_leaflet.size(); m++) //loop over opposing leaflet atoms
                                {
                                    //get the first and last atom of the current lipid
                                    int min_2 = traj.o_lip_start(m);
                                    int max_2 = traj.o_lip_end(m);

                                    //jump to the next lipid
                                    m = traj.next_opposing_lipid(m);

                                    for(n=0; n<param_2.main_size_y(); n++) //loop over lipid types 
                                    {
                                        if(strcmp(traj.res_name[min_2].c_str(), param_2.param_main_s[n][0].c_str()) == 0) //lipid type is correct
                                        {
                                            for(o=min_2; o<=max_2; o++) //loop over current lipid atoms 
                                            {
                                                for(q=0; q<param_2.sec_size_y(n); q++) //loop over ranked atoms
                                                {
                                                    if(strcmp(traj.atom_name[o].c_str(), param_2.param_sec_s[n][q][0].c_str()) == 0) //atoms is a ranked atom
                                                    {
                                                        //compute the distance between the atoms
                                                        double dx = traj.r[k][0] - traj.r[o][0];
                                                        double dy = traj.r[k][1] - traj.r[o][1];
                                                        double dz = traj.r[k][2] - traj.r[o][2];

                                                        double distance = sqrt(dx*dx + dy*dy + dz*dz);

                                                        if(distance < p.contact_cutoff) //atoms make a contact
                                                        {
                                                            contacts = contacts + 1.0;
                                                            sum_rank = sum_rank + param_2.param_sec_d[n][q][1];
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }

                                double avg = sum_rank/contacts;
                                if(contacts > 0)
                                {
                                    digi.stamp(hx,hy,p.radius,sum_rank/contacts);
                                }
                                else 
                                {
                                    digi.stamp(hx,hy,p.radius,0.0);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    //get the average for the current frame
    digi.norm_frame();

    //add digi for each gridpoint to the long term sum
    digi.add_frame();

    //now we print the single frame digi data
    if(p.b_stdev == 1)
    {
        digi.write_frame(traj.get_frame_global());
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// collect digi and compute average.                                                                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,Grid &digi)
{
    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nfinalizing analysis. This requires communication of the grid and could take some time. \n");
    }

    //collect digi and rho from all ranks
    digi.collect_grid();

    //normalize digi and write the <digi> and rho to file
    if(s.world_rank == 0)
    {
        digi.normalize();

        digi.exclude_data(p.cutoff,1);

        digi.write_grid();
        digi.write_rho();
    }

    //compute standard deviation
    digi.get_stdev(p.b_stdev,p.b_clean,traj);

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
    s.program_name = "Interdigitation";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                      s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                      s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                     s.world_rank, s.cl_tags, &p.b_print,   0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                                  s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                   s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                     s.world_rank, s.cl_tags, &p.b_lsq,     0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-crd_1",  p.param_1_file_name,          "Selection card with lipid types, mapping, and tail atoms (crd)",  s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-crd_2",  p.param_2_file_name,          "Selection card with opposing lipid types and atom ranking (crd)", s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-digi",   p.digi_file_name,             "Output file with spatially resolved average rank (dat)",          s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,           "PDB file with sorted leaflets (pdb)",                             s.world_rank, s.cl_tags, &p.b_lf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,  "File with additional leaflet finder parameters (prm)",            s.world_rank, s.cl_tags, &p.b_lf_param,0);
    add_argument_mpi_d(argc,argv,"-APS",    &p.APS,                       "Area per grid box (nm^2)",                                        s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-r",      &p.radius,                    "Radius of the target atom (nm)",                                  s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-cutoff", &p.cutoff,                    "Cutoff for excluding data (chi)",                                 s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-bx",     &p.box_x,                     "Grid x dimension (nm)",                                           s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-by",     &p.box_y,                     "Grid y dimension (nm)",                                           s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-leaf",   &p.leaflet,                   "Which leaflet? (0:both 1:upper 2:lower)",                         s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-odf",    &p.out_data_format,           "What is the output data format? (0:matrix 1:vector)",             s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-cdist",  &p.contact_cutoff,            "The contact cutoff distance (nm)",                                s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-stdev",  &p.b_stdev,                   "Compute the STDEV and STEM? (0:no 1:yes)",                        s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-clean",  &p.b_clean,                   "Remove single frame files? (0:no 1:yes)",                         s.world_rank, s.cl_tags, nullptr,      0);
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
    check_extension_mpi(s.world_rank,"-digi",p.digi_file_name,".dat");

    if(p.b_lf_pdb == 1)
    {
        check_extension_mpi(s.world_rank,"-lf_pdb",p.lf_pdb_file_name,".pdb");
    }
    if(p.b_lf_param == 1)
    {
        check_extension_mpi(s.world_rank,"-lf_prm",p.leaflet_finder_param_name,".prm");
    }

    //create index objects
    Index param_1;

    //read the index files
    param_1.get_index(p.param_1_file_name);

    //create parameter files
    Param param_2;

    //read parameter files
    param_2.get_param(p.param_2_file_name,2,1,2);

    //check the integrity of the parameter files
    if(param_2.check_file() == 0) //bad files. kill program 
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
    traj.get_lipid_selection_stats(param_1.get_column_s(5,0),"-crd_1");
    traj.get_lipid_selection_stats(param_2.get_column_s(0),  "-crd_2");

    //create a grid to hold interdigitation
    Grid digi;

    //get the grid dimensions
    digi.get_dim(p.box_x,p.box_y,traj.ibox,p.APS);

    //set the output file name for grid
    digi.set_output(p.digi_file_name,p.out_data_format);

    //print info about the grid
    digi.print_dim();
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

        interdigitate(traj,s,p,param_1,param_2,digi);

        traj.set_beta_lf();

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect interdigitation data from mpi ranks and compute the average
    perf.log_time(finalize_analysis(traj,s,p,digi),"Fin Ana");

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

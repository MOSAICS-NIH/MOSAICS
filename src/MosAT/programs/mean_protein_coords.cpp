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
#include "MosAT/program_variables/pv_mean_protein_coords.h" //This has the variables specific to the analysis program
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
// This function adds the protein coordinates to the running sum.                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void add_r(Trajectory &traj,system_variables &s,program_variables &p)
{
    int i = 0;
    int j = 0;

    for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
    {
        for(j=0; j<3; j++) //loop over 3-dimensions
        {
            p.r_avg[i][j] = p.r_avg[i][j] + traj.r[traj.prot[i]-1][j]/(double)traj.get_ef_frames();
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects data from nodes and prints the final coords to pdb file.                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p)
{
    int i = 0;
    int j = 0;
    int k = 0;

    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nCollecting data and computing mean coords. \n");
    }

    //now we collect r_avg from all nodes
    double avg = 0;
    double avgs[s.world_size];
    for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
    {
        for(j=0; j<3; j++) //loop over x,y and z
        {
            avg = p.r_avg[i][j];

            MPI_Gather(&avg, 1, MPI_DOUBLE, avgs, 1, MPI_DOUBLE, 0,MPI_COMM_WORLD);
            if(s.world_rank == 0)
            {
                avg = 0;
                for(k=0; k<s.world_size; k++)
                {
                    avg = avg + avgs[k];
                }
                p.r_avg[i][j] = avg;
            }
        }
    }

    //print average coords to pdb file
    if(s.world_rank == 0 && p.b_mean_dist == 0)
    {
        //create arrays to hold pdb stuff
        vector <int>    set_atom_nr(traj.prot.size(),0);                  //atom_nr for pdb files
        vector <int>    set_res_nr(traj.prot.size(),0);                   //res_nr for pdb files
        vector <string> set_atom_name(traj.prot.size());                  //atom name for pdb files
        vector <string> set_res_name(traj.prot.size());                   //res name for pdb files
        vector <double> set_beta(traj.prot.size(),0.0);                   //Beta factor for pdb files
        vector <double> set_weight(traj.prot.size(),0.0);                 //Weight used for pdb files
        vector <string> set_element(traj.prot.size());                    //Element collumn in pdb file
        vector <char>   set_chain_id(traj.prot.size());                   //Chain id for pdb file

        //create info needed for pdb and get atomic coords
        for(i=0; i<traj.prot.size(); i++) //loop over lipid atoms
        {
            set_atom_nr[i]   = i+1;
            set_res_nr[i]    = traj.res_nr[traj.prot[i]-1];
            set_atom_name[i] = traj.atom_name[traj.prot[i]-1];
            set_res_name[i]  = traj.res_name[traj.prot[i]-1];
            set_beta[i]      = 0;
            set_weight[i]    = 0;
            set_element[i]   = "ca";
            set_chain_id[i]  = 'A';
        }

        //open file for printing output
        FILE *mpc_file = fopen(p.mpc_file_name.c_str(), "w");
        if(mpc_file == NULL)
        {
            printf("failure opening %s (avg coord). Make sure the file exists. \n",p.mpc_file_name.c_str());
        }

        //write the protein to a pdb file
        write_frame_pdb(traj.ibox,traj.prot.size(),set_atom_nr,set_res_nr,set_res_name,set_atom_name,p.r_avg,traj.title,s.world_rank,&mpc_file,set_beta,set_weight,set_element,set_chain_id,1);

        fclose(mpc_file);
    }

    //now we broadcast the average coords_x,y,z so all cores have the average
    //this is needed only if we are computing the average dist from the average coords
    if(p.b_mean_dist == 1)
    {
        if(s.world_size > 1)
        {
            for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
            {
                for(j=0; j<3; j++)
                {
                    MPI_Bcast(&p.r_avg[i][j], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                }
                //printf("world_rank %3d r_avg[%4d][0] %10f r_avg[%4d][1] %10f r_avg[%4d][2] %10f \n",world_rank,i,r_avg[i][0],i,r_avg[i][1],i,r_avg[i][2]);
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //compute and return time spent in function
    return (clock() - s.t)/CLOCKS_PER_SEC;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the distance the coords are from the average coords.                               //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void prot_dist(Trajectory &traj,system_variables &s,program_variables &p,double avg_dist[],double largest_dist[])
{
    int i = 0;
    int j = 0;
    for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
    {
        double dx = traj.r[traj.prot[i]-1][0] - p.r_avg[i][0];
        double dy = traj.r[traj.prot[i]-1][1] - p.r_avg[i][1];
        double dz = traj.r[traj.prot[i]-1][2] - p.r_avg[i][2];

        double dist = sqrt(dx*dx + dy*dy + dz*dz);

        avg_dist[i] = avg_dist[i] + dist;

        if(dist > largest_dist[i])
        {
            largest_dist[i] = dist;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects data from nodes and prints the final coords to pdb file.                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void finalize_analysis_dist(Trajectory &traj,system_variables &s,program_variables &p,double avg_dist[],double largest_dist[])
{
    if(p.b_mean_dist == 1)
    {
        int i = 0;
        int j = 0;
        int k = 0;

        MPI_Barrier(MPI_COMM_WORLD);

        if(s.world_rank == 0)
        {
            printf("Setting beta factor to the average distance. \n");
        }

        //now we collect avg_dist from all nodes
        double avg = 0;
        double world_avg[s.world_size];
        if(s.world_size > 0)
        {
            for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
            {
                avg = avg_dist[i];

                MPI_Gather(&avg, 1, MPI_DOUBLE, world_avg, 1, MPI_DOUBLE, 0,MPI_COMM_WORLD);
                if(s.world_rank == 0)
                {
                    avg = 0;
                    for(j=0; j<s.world_size; j++)
                    {
                        avg = avg + world_avg[j];
                    }
                    avg_dist[i] = avg/traj.get_ef_frames();
                }
            }
        }

        //now we collect largest_dist from all nodes
        double largest = 0;
        double world_largest[s.world_size];
        if(s.world_size > 0)
        {
            for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
            {
                largest = largest_dist[i];

                MPI_Gather(&largest, 1, MPI_DOUBLE, world_largest, 1, MPI_DOUBLE, 0,MPI_COMM_WORLD);
                if(s.world_rank == 0)
                {
                    largest = 0;
                    for(j=0; j<s.world_size; j++)
                    {
                        if(world_largest[j] > largest)
                        {
                            largest = world_largest[j];
                        }
                    }
                    largest_dist[i] = largest;
                }
            }
        }

        //print average coords to pdb file
        if(s.world_rank == 0)
        {
            //create arrays to hold pdb stuff
            vector <int>    set_atom_nr(traj.prot.size(),0);                  //atom_nr for pdb files
            vector <int>    set_res_nr(traj.prot.size(),0);                   //res_nr for pdb files
            vector <string> set_atom_name(traj.prot.size());                  //atom name for pdb files
            vector <string> set_res_name(traj.prot.size());                   //res name for pdb files
            vector <double> set_beta(traj.prot.size(),0.0);                   //Beta factor for pdb files
            vector <double> set_weight(traj.prot.size(),0.0);                 //Weight used for pdb files
            vector <string> set_element(traj.prot.size());                    //Element collumn in pdb file
            vector <char>   set_chain_id(traj.prot.size());                   //Chain id for pdb file

            //create info needed for pdb and get atomic coords
            for(i=0; i<traj.prot.size(); i++) //loop over lipid atoms
            {
                set_atom_nr[i]   = i+1;
                set_res_nr[i]    = traj.res_nr[traj.prot[i]-1];
                set_atom_name[i] = traj.atom_name[traj.prot[i]-1];
                set_res_name[i]  = traj.res_name[traj.prot[i]-1];
                set_beta[i]      = 0;
                set_weight[i]    = 0;
                set_element[i]   = "ca";
                set_chain_id[i]  = 'A';
            }

            //set beta value to the average dist
            for(i=0; i<traj.prot.size(); i++)
            {
                set_beta[i] = avg_dist[i];
            }

            string mean_dist_file_name = add_tag(p.mpc_file_name,"_mean_dist");

            //open file for printing output
            FILE *mpc_file = fopen(mean_dist_file_name.c_str(), "w");
            if(mpc_file == NULL)
            {
                printf("failure opening %s (avg coord). Make sure the file exists. \n",mean_dist_file_name.c_str());
            }

            //write the protein to a pdb file
            write_frame_pdb(traj.ibox,traj.prot.size(),set_atom_nr,set_res_nr,set_res_name,set_atom_name,p.r_avg,traj.title,s.world_rank,&mpc_file,set_beta,set_weight,set_element,set_chain_id,1);

            fclose(mpc_file);

            //set beta value to the largest dist
            for(i=0; i<traj.prot.size(); i++)
            {
                set_beta[i] = largest_dist[i];
            }

            string largest_dist_file_name = add_tag(p.mpc_file_name,"_largest_dist");

            //open file for printing output
            FILE *largest_dist_file = fopen(largest_dist_file_name.c_str(), "w");
            if(largest_dist_file == NULL)
            {
                printf("failure opening %s (avg coord). Make sure the file exists. \n",largest_dist_file_name.c_str());
            }

            //write the protein to a pdb file
            write_frame_pdb(traj.ibox,traj.prot.size(),set_atom_nr,set_res_nr,set_res_name,set_atom_name,p.r_avg,traj.title,s.world_rank,&largest_dist_file,set_beta,set_weight,set_element,set_chain_id,1);

            fclose(largest_dist_file);
        }

        MPI_Barrier(MPI_COMM_WORLD);
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
    s.program_name = "Mean Protein Coords";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                   s.world_rank, nullptr,          1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                   s.world_rank, nullptr,          1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                  s.world_rank, &p.b_print,       0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                    s.world_rank, nullptr,          0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                               s.world_rank, nullptr,          0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                s.world_rank, nullptr,          0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                  s.world_rank, &p.b_lsq,         0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                    s.world_rank, nullptr,          0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",    s.world_rank, nullptr,          0);
    add_argument_mpi_s(argc,argv,"-mpc",    p.mpc_file_name,              "Output file with time average protein coords (pdb)",           s.world_rank, nullptr,          1);
    add_argument_mpi_s(argc,argv,"-pf_pdb", p.pf_pdb_file_name,           "PDB file with protein indicated (pdb)",                        s.world_rank, &p.print_prot_pdb,0);
    add_argument_mpi_s(argc,argv,"-pf_prm", p.protein_finder_param_name,  "File with additional protein finder parameters (prm)",         s.world_rank, &p.b_pf_param,    0);
    add_argument_mpi_i(argc,argv,"-dist",   &p.b_mean_dist,               "Compute the average dist from the average coords (0:no,1:yes)",s.world_rank, nullptr,          0);
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
    check_extension_mpi(s.world_rank,"-mpc",p.mpc_file_name,".pdb");
    if(p.print_prot_pdb == 1)
    {
        check_extension_mpi(s.world_rank,"-pf_pdb",p.pf_pdb_file_name,".pdb");
    }   
    if(p.b_pf_param == 1)
    {
        check_extension_mpi(s.world_rank,"-pf_prm",p.protein_finder_param_name,".prm");
    }

    //run protein finder
    traj.get_protein(p.protein_finder_param_name,p.b_pf_param);

    //write pdb with protein indicated by beta factor
    traj.write_protein(p.pf_pdb_file_name,p.print_prot_pdb);

    //print info about the protein
    traj.get_prot_stats();

    //set protein size
    p.prot_size = traj.prot.size();

    //allocate memory for average coords
    //snew(p.r_avg,p.prot_size);
    p.r_avg = (rvec*)calloc(p.prot_size , sizeof(rvec));

    //allocate memory for average dist
    double avg_dist[p.prot_size];
    double largest_dist[p.prot_size];
    init_darray(avg_dist,p.prot_size);
    init_darray(largest_dist,p.prot_size);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //print info about the worlk load distribution
    traj.workload();

    s.t = clock();
    //read read frames of the trajector and perform analysis
    for(traj.current_frame=0; traj.current_frame<traj.get_num_frames(); traj.current_frame++)
    {
        traj.read_traj_frame();

        traj.do_fit();

        add_r(traj,s,p);

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //collect coords from mpi processes and compute average
    perf.log_time(finalize_analysis(traj,s,p),"Fin Ana");

    //print that analysis is beginning
    traj.report_progress();

    //loop over trajectory a second time to compute the average distance from average coords
    if(p.b_mean_dist == 1)
    {
        for(traj.current_frame=0; traj.current_frame<traj.get_num_frames(); traj.current_frame++)
        {
            traj.read_traj_frame();

            traj.do_fit();

            prot_dist(traj,s,p,avg_dist,largest_dist);
        }
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect distances from mpi processes and compute average
    finalize_analysis_dist(traj,s,p,avg_dist,largest_dist);

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

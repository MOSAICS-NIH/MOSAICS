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
#include "MosAT/program_variables/pv_group_contacts.h"       //This has the variables specific to the analysis program
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
#include "headers/atom_select.h"                             //This has routines used for making atom selections using a selection text
#include "headers/param.h"                                   //This has routines used for reading complex parameter data

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function checks if a contact is formed                                                               //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int check_contact(Trajectory &traj,program_variables &p,int a1,int a2)
{
    int result = 0;               //tells if the pair under investigation make a contact
    int i      = 0;               //standard variable used in loops

    rvec m;                       //difference vector between the atoms

    for(i=0; i<3; i++) //loop over 3 dimensions
    {
        m[i] = traj.r[a1][i] - traj.r[a2][i];
    }

    double dist = sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);   //distance between atoms

    if(dist < p.cdist)
    {
        result = 1;

        //check contacts in pymol
        if(p.b_test == 1)
        {
            printf("atom_1 %10d %10s %10d %10s \n",traj.atom_nr[a1],traj.atom_name[a1].c_str(),traj.res_nr[a1],traj.res_name[a1].c_str());
            printf("atom_2 %10d %10s %10d %10s \n",traj.atom_nr[a2],traj.atom_name[a2].c_str(),traj.res_nr[a2],traj.res_name[a2].c_str());
            printf("sel a1, resi %d & resn %s & name %s \n",traj.res_nr[a1]%10000,traj.res_name[a1].c_str(),traj.atom_name[a1].c_str());
            printf("sel a2, resi %d & resn %s & name %s \n",traj.res_nr[a2]%10000,traj.res_name[a2].c_str(),traj.atom_name[a2].c_str());
            printf("dist (a1), (a2) \n");
            printf("show licorice, resi %d \n",traj.res_nr[a1]);
            printf("show licorice, resi %d \n",traj.res_nr[a2]);
            printf("sel pair, a1 + a2 \n");
            printf("orient pair \n");
            printf("\n");
        }
    }

    return result; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function takes in the index data and tags the appropriate atoms for easy checks                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void tag_group_atoms(Trajectory &traj,system_variables &s,program_variables &p,iv1d &target_1,Index &n1)
{
    int i = 0;    //standard variable used in loops
    int j = 0;    //standard variable used in loops
    int k = 0;    //standard variable used in loops
    int l = 0;    //standard variable used in loops

    int pos = 0;
    string tag;

    for(i=0; i<n1.index_s.size(); i++)
    {
        target_1[n1.index_i[i]-1] = 1;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function checks if a set of atom pairs has been encountered befor to prevent double counting         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int check_pairs(int a1,int a2,iv1d &a1_found,iv1d &a2_found)
{
    int i = 0;
    int found = 0;

    for(i=0; i<a1_found.size(); i++)
    {
        if(a1 == a1_found[i] && a2 == a2_found[i])
        {
            found = 1;
        }
	else if(a2 == a1_found[i] && a1 == a2_found[i])
        {
            found = 1;
        }
    }

    return found;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the number of contacts between 2 groups                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_contacts(Trajectory &traj,system_variables &s,program_variables &p,iv1d &target_1,iv1d &target_2,iv1d &contact_count)
{
    int    i        = 0;                      //standard variable used in loops
    int    j        = 0;                      //standard variable used in loops
    int    k        = 0;                      //standard variable used in loops
    int contacts    = 0;                      //number of h-bonds for current frame

    iv1d a1_found(0,0); //used to keep track of contacts encountered so we dont double count
    iv1d a2_found(0,0); //used to keep track of contacts encountered so we dont double count

    for(i=0; i<traj.atoms(); i++) //loop over atoms
    {
        if(target_1[i] == 1) //residue atom is an acceptor
        {
            int atom_1 = i;

            for(j=0; j<traj.atoms(); j++) //loop over atoms
            {
                if(target_2[j] == 1)
                {
                    int atom_2 = j;

                    if(traj.res_nr[atom_1] != traj.res_nr[atom_2])
                    {    
                        if(check_contact(traj,p,atom_1,atom_2) == 1)
                        {
                            if(check_pairs(atom_1,atom_2,a1_found,a2_found) == 0)
                            {
                                a1_found.push_back(atom_1);
                                a2_found.push_back(atom_2);

                                contacts = contacts + 1;
                            }
                        }
                    }
                }
            }
        }
    }

    contact_count.push_back(contacts);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Collect the contacts and compute the average                                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,iv1d &contact_count)  
{
    int    i        = 0;                      //standard variable used in loops
    int    j        = 0;                      //standard variable used in loops
    int    k        = 0;                      //standard variable used in loops
    int    l        = 0;                      //standard variable used in loops

    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nFinalizing analysis. This requires communicating contacts data and could take some time. \n");
    }

    collect_iv1d(s.world_size,s.world_rank,contact_count);

    if(s.world_rank == 0)
    {
        FILE *contact_file = fopen(p.contact_file_name.c_str(),"w");
        fprintf(contact_file," %9s   %10s \n","#step","#contacts");
        for(i=0; i<contact_count.size(); i++)
        {
            fprintf(contact_file," %9d   %10d \n",i,contact_count[i]);
        }
        fclose(contact_file);

        //report average:
        double avg = 0.0;
        for(i=0; i<contact_count.size(); i++)
        {
            avg = avg + (double)contact_count[i];
        }
        avg = avg/(double)contact_count.size();
        printf("\naverage contacts per frame: %f \n",avg);
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
    s.program_name = "Group Contacts";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                   s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                   s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                  s.world_rank, s.cl_tags, &p.b_print,   0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                    s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                               s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                  s.world_rank, s.cl_tags, &p.b_lsq,     0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                    s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",    s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-cont",   p.contact_file_name,          "Output file with contact count (dat)",                         s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-test",   &p.b_test,                    "Print info for checking contacts? (0:no 1:yes)",               s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-n1",     p.n1_file_name,               "Index for group 1 atoms (ndx)",                                s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-n2",     p.n2_file_name,               "Index for group 2 atoms (ndx)",                                s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-cdist",  &p.cdist,                     "Distance cutoff for counting contacts (nm)",                   s.world_rank, s.cl_tags, nullptr,      1);
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
    check_extension_mpi(s.world_rank,"-n1",p.n1_file_name,".ndx");
    check_extension_mpi(s.world_rank,"-n2",p.n2_file_name,".ndx");

    //create index objects
    Index n1;                //holds atom group 1 
    Index n2;                //holds atom group 2

    //read the index files
    n1.get_index(p.n1_file_name);
    n2.get_index(p.n2_file_name);

    //create structures for tagging target atoms
    iv1d target_1(traj.atoms(),0);       //tags all atoms for group 1
    iv1d target_2(traj.atoms(),0);       //tags all atoms for group 2

    //tag the target atoms
    tag_group_atoms(traj,s,p,target_1,n1);
    tag_group_atoms(traj,s,p,target_2,n2);

    iv1d  contact_count(0,0);   //store number of contacts for each frame

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

        get_contacts(traj,s,p,target_1,target_2,contact_count);

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect contacts from mpi processes and compute the average
    perf.log_time(finalize_analysis(traj,s,p,contact_count),"Fin Ana");

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

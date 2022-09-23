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

#include "xdr/include/xdrfile_xtc.h"                            //used to read xtc files 
#include "xdr/include/xdr_seek.h"                               //used to get and set the file position in xtc and trr files
#include "xdr/include/xdrfile_trr.h"                            //used to read trr files
#include "xdr/include/xdrfile.h"                                //used to read xtc and trr files
#include "xdr/include/trr_header.h"                             //used to read the header info of trr files
#include "headers/multi_dim_vec.h"                              //This defines multidimensional vectors
#include "headers/switch.h"                                     //This defines a switch (on, off)
#include "headers/file_reader.h"                                //This has basic routines for reading text files
#include "headers/vector_mpi.h"                                 //This has routines for collecting vector data
#include "headers/mosat_routines.h"                            //This is where most of the functions called in main are located
#include "headers/file_naming.h"                                //This has routines for added tags to an existing file name    
#include "headers/file_naming_mpi.h"                            //This has routines for added tags to an existing file name (mpi)
#include "headers/command_line_args_mpi.h"                      //This has routines for adding command line arguments
#include "MosAT/program_variables/pv_protein_lipid_contacts.h" //This has the variables specific to the analysis program
#include "headers/array.h"                                      //This has routines used for working with arrays
#include "headers/performance.h"                                //This has a class for logging performance data
#include "headers/index.h"                                      //This has a class for working with index files
#include "headers/traj.h"                                       //This has a class for working with the trajectory
#include "headers/leaflet_finder.h"                             //This has routines used to find leaflets in membrane simulations
#include "headers/protein_finder.h"                             //This has routines used to find protein atoms
#include "headers/sol_finder.h"                                 //This has routines used to find the solvent
#include "headers/grid.h"                                       //This has routines used for working with a grid
#include "headers/protein.h"                                    //This has routines used for working with protein data
#include "headers/force_serial.h"                               //This has routines used for forcing the code to run on a single mpi process
#include "headers/param.h"                                      //This has routines used for reading complex parameter data

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function counts how many contacts are made between the lipids and each protein residue.              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void prot_lip_contacts(Trajectory &traj,system_variables &s,program_variables &p,Param &param,int contacts_sum[])
{
    int i        = 0;             //loop variable
    int j        = 0;             //loop variable
    int k        = 0;             //loop variable
    int l        = 0;             //loop variable
    int m        = 0;             //loop variable
    int n        = 0;             //loop variable

    for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
    {
        for(j=0; j<traj.target_leaflet.size(); j++) //loop over target leaflet atoms
        {
            for(k=0; k<param.main_size_y(); k++) //loop over lipid types 
            {
                if(strcmp(traj.res_name[traj.target_leaflet[j]-1].c_str(), param.param_main_s[k][0].c_str()) == 0) //residue is a target lipid
                {
                    for(l=0; l<param.sec_size_y(k); l++) //loop over target atoms
                    {
                        if(strcmp(traj.atom_name[traj.target_leaflet[j]-1].c_str(), param.param_sec_s[k][l][0].c_str()) == 0) //atom is a target atom
                        {
                            double dif_x = traj.r[traj.prot[i]-1][0] - traj.r[traj.target_leaflet[j]-1][0];
                            double dif_y = traj.r[traj.prot[i]-1][1] - traj.r[traj.target_leaflet[j]-1][1];
                            double dif_z = traj.r[traj.prot[i]-1][2] - traj.r[traj.target_leaflet[j]-1][2];

                            double dist = sqrt(dif_x*dif_x + dif_y*dif_y + dif_z*dif_z);

                            if(dist < p.contact_cutoff)
                            {
                                contacts_sum[traj.res_nr[traj.prot[i]-1]-1] = contacts_sum[traj.res_nr[traj.prot[i]-1]-1] + 1;
                            }
                        }
                    }
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects contact data from the ranks and prints the final results to a text and pdb file.   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,int contacts_sum[],int num_residues)
{
    int i = 0;                                       //standard variable used in loops
    int j = 0;                                       //standard variable used in loops
    int k = 0;                                       //standard variable used in loops
    FILE *contacts_file;                             //file for writing contacts to text 
    FILE *contacts_pdb_file;                         //file for writing contacts to pdb

    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nFinalizing analysis. This could take some time. \n");
    }

    //open files for writing
    if(s.world_rank == 0)
    {
        contacts_file = fopen(p.contacts_file_name.c_str(), "w");
        if(contacts_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",p.contacts_file_name.c_str());
        }
        else 
        {
            fprintf(contacts_file,"# %15s %15s \n","column 1:","res id");
            fprintf(contacts_file,"# %15s %15s \n","column 2:","Contacts formed per trajectory frame");
        }

        string contacts_pdb_file_name = chop_and_add_tag(p.contacts_file_name,".pdb");

        contacts_pdb_file = fopen(contacts_pdb_file_name.c_str(), "w");
        if(contacts_pdb_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",contacts_pdb_file_name.c_str());
        }
    }

    //now we collect the contacts from all nodes
    int contact = 0;
    int world_contacts[s.world_size];
    for(i=0; i<num_residues; i++)
    {
        contact = contacts_sum[i];

        MPI_Gather(&contact, 1, MPI_INT, world_contacts, 1, MPI_INT, 0,MPI_COMM_WORLD);
        if(s.world_rank == 0)
        {
            contact = 0;
            for(j=0; j<s.world_size; j++)
            {
                contact = contact + world_contacts[j];
            }
            contacts_sum[i] = contact;

            fprintf(contacts_file,"%15d %15f \n",i+1,(double)contacts_sum[i]/(double)traj.get_ef_frames());
        }
    }

    //now we print the contacts data to a pdb file as a B-factor
    if(s.world_rank == 0)
    {
        //here we set the beta factors
        for(i=0; i<traj.atoms(); i++)
        {
            traj.beta[i] = (double)contacts_sum[traj.res_nr[i] - 1]/(double)traj.get_ef_frames();
        }
        write_frame_pdb(traj.box_ref,traj.atoms(),traj.atom_nr,traj.res_nr,traj.res_name,traj.atom_name,traj.r_ref,traj.title,s.world_rank,&contacts_pdb_file,traj.beta,traj.weight,traj.element,traj.chain_id,0);

        fclose(contacts_file);
        fclose(contacts_pdb_file);
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
    s.program_name = "Protein Lipid Contacts";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                     s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                     s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                    s.world_rank, &p.b_print,   0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                      s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                                 s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                  s.world_rank, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                    s.world_rank, &p.b_lsq,     0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                      s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",      s.world_rank, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-crd",    p.param_file_name,            "Selection card with target lipids and atoms (crd)",              s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ct",     p.contacts_file_name,         "Output file with the number of contacts for each residue (dat)", s.world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-leaf",   &p.leaflet,                   "Which leaflet? (0:both 1:upper 2:lower)",                        s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,           "PDB file with sorted leaflets (B-factor) (pdb)",                 s.world_rank, &p.b_lf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,  "File with additional leaflet finder parameters (prm)",           s.world_rank, &p.b_lf_param,0);
    add_argument_mpi_s(argc,argv,"-pf_pdb", p.pf_pdb_file_name,           "PDB file with selected protein (B-factor) (pdb)",                s.world_rank, &p.b_pf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-pf_prm", p.protein_finder_param_name,  "File with additional protein finder parameters (prm)",           s.world_rank, &p.b_pf_param,0);
    add_argument_mpi_d(argc,argv,"-cdist",  &p.contact_cutoff,            "The contact cutoff distance (nm)",                               s.world_rank, nullptr,      1);
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
    check_extension_mpi(s.world_rank,"-crd",p.param_file_name,".crd");
    check_extension_mpi(s.world_rank,"-ct",p.contacts_file_name,".dat");

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
    Param param;

    //read parameter files
    param.get_param(p.param_file_name,2,1,1);

    //check the integrity of the parameter files
    if(param.check_file() == 0) //bad files. kill program 
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
    traj.get_lipid_selection_stats(param.get_column_s(0),"-crd");

    //run protein finder
    traj.get_protein(p.protein_finder_param_name,p.b_pf_param);

    //write pdb with protein indicated by beta factor
    traj.write_protein(p.pf_pdb_file_name,p.b_pf_pdb);

    //print info about the protein
    traj.get_prot_stats();

    int num_residues = traj.res_nr[traj.atoms()-1]; //How many residues does the system have
    int contacts_sum[num_residues];                 //Stores the number of contacts for each residue
    init_iarray(contacts_sum,num_residues);
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

        prot_lip_contacts(traj,s,p,param,contacts_sum);

        traj.set_beta_lf();

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect data from mpi processes and compute averages
    perf.log_time(finalize_analysis(traj,s,p,contacts_sum,num_residues),"Fin Ana");

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

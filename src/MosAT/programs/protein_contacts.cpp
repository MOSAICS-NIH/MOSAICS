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
#include "MosAT/program_variables/pv_protein_contacts.h"        //This has the variables specific to the analysis program
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
// This function adds the protein coordinates to the running sum.                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void prot_coords(Trajectory &traj,system_variables &s,program_variables &p,dv1d &prot_x,dv1d &prot_y,dv1d &prot_z)
{
    int i = 0;

    for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
    {
        prot_x[i] = prot_x[i] + traj.r[traj.prot[i]-1][0];
        prot_y[i] = prot_y[i] + traj.r[traj.prot[i]-1][1];
        prot_z[i] = prot_z[i] + traj.r[traj.prot[i]-1][2];
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function counts how many contacts are made between the lipids and each protein residue.              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void prot_lip_contacts(Trajectory &traj,system_variables &s,program_variables &p,Index &group_1,dv1d &contacts,iv1d &contacts_sum)
{
    int i        = 0;             //loop variable
    int j        = 0;             //loop variable

    for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
    {
        for(j=0; j<group_1.index_i.size(); j++)
        { 
            double dx = traj.r[traj.prot[i]-1][0] - traj.r[group_1.index_i[j]-1][0];
            double dy = traj.r[traj.prot[i]-1][1] - traj.r[group_1.index_i[j]-1][1];
            double dz = traj.r[traj.prot[i]-1][2] - traj.r[group_1.index_i[j]-1][2];

            double dist = sqrt(dx*dx + dy*dy + dz*dz);

            if(dist < p.contact_cutoff)
            {
                contacts[i] = contacts[i] + 1.0;
                contacts_sum[traj.res_nr[traj.prot[i]-1]-1] = contacts_sum[traj.res_nr[traj.prot[i]-1]-1] + 1;
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects contact data from the ranks and prints the final results to a text and pdb file.   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,dv1d &contacts,iv1d &contacts_sum,int num_residues,dv1d &prot_x,dv1d &prot_y,dv1d &prot_z)
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

    //collect data
    collect_and_sum_dv1d(s.world_size,s.world_rank,contacts);

    //collect atomic protein coordinates from mpi cores   
    collect_and_sum_dv1d(s.world_size,s.world_rank,prot_x);
    collect_and_sum_dv1d(s.world_size,s.world_rank,prot_y);
    collect_and_sum_dv1d(s.world_size,s.world_rank,prot_z);

    //now we print the contacts data to a pdb file as a B-factor
    if(s.world_rank == 0)
    {
        //normalize protein atom data 
        for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
        {
            prot_x[i] = prot_x[i]/(double)traj.get_ef_frames();
            prot_y[i] = prot_y[i]/(double)traj.get_ef_frames();
            prot_z[i] = prot_z[i]/(double)traj.get_ef_frames();
        }

        int this_size = traj.prot.size();

        //allocate memory to store info needed for writing out a pdb
        iv1d    this_atom_nr(this_size,0);                  //atom number used in pdb file
        iv1d    this_res_nr(this_size,1);                   //res number used in pdb file
        sv1d    this_atom_name(this_size);                  //atom name used in pdb file
        sv1d    this_res_name(this_size);                   //res name used in pdb file
        dv1d    this_beta(this_size,1.0);                   //beta factor used in pdb file
        dv1d    this_weight(this_size,1.0);                 //weight used ub  pdb file
        sv1d    this_element(this_size,"ele");              //element collumn used in pdb file
        cv1d    this_chain_id(this_size,'A');               //chain id used in pdb file

        //set the atom number
        for(i=0; i<this_size; i++) //loop over atoms
        {
            this_atom_nr[i] = i+1;
        }

        //set the atom names
        for(i=0; i<this_size; i++) //loop over protein atoms 
        {
            this_atom_name[i] = traj.atom_name[traj.prot[i]-1];
        }

        //set the residue names
        for(i=0; i<this_size; i++) //loop over protein atoms
        {
            this_res_name[i] = traj.res_name[traj.prot[i]-1];
        }

        //set residue numbers
        for(i=0; i<this_size; i++) //loop over protein atoms
        {
            this_res_nr[i] = traj.res_nr[traj.prot[i]-1];
        }

        //allocate memory for the coordinates
        rvec *this_r;
        this_r = (rvec *)calloc(this_size , sizeof(*this_r));

        for(i=0; i<this_size; i++) //loop over protein atoms
        {
            this_r[i][0] = prot_x[i];
            this_r[i][1] = prot_y[i];
            this_r[i][2] = prot_z[i];
        }

        //set the b factors
        for(i=0; i<this_size; i++) //loop over protein atoms
        {
            this_beta[i] = contacts[i]/(double)traj.get_ef_frames(); 
        }

        FILE *pdb_file;
        string pdb_file_name = chop_and_add_tag(p.contacts_file_name,".pdb");
        pdb_file             = fopen(pdb_file_name.c_str(), "w");
        if(pdb_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",pdb_file_name.c_str());
        }
        else
        {
            write_frame_pdb(traj.box_ref,this_size,this_atom_nr,this_res_nr,this_res_name,this_atom_name,this_r,traj.title,s.world_rank,&pdb_file,this_beta,this_weight,this_element,this_chain_id,0);
            fclose(pdb_file);
        }
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
    s.program_name = "Protein Contacts";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                     s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                     s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                    s.world_rank, s.cl_tags, &p.b_print,   0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                      s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                                 s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                  s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                    s.world_rank, s.cl_tags, &p.b_lsq,     0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                      s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",      s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-n1",     p.group_1_file_name,          "Index file with group 1 atoms (ndx)",                            s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-pc",     p.contacts_file_name,         "Output file with the number of contacts for each atom (dat)",    s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-pf_pdb", p.pf_pdb_file_name,           "PDB file with selected protein (B-factor) (pdb)",                s.world_rank, s.cl_tags, &p.b_pf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-pf_prm", p.protein_finder_param_name,  "File with additional protein finder parameters (prm)",           s.world_rank, s.cl_tags, &p.b_pf_param,0);
    add_argument_mpi_d(argc,argv,"-cdist",  &p.contact_cutoff,            "The contact cutoff distance (nm)",                               s.world_rank, s.cl_tags, nullptr,      1);
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
    check_extension_mpi(s.world_rank,"-pc",p.contacts_file_name,".dat");

    if(p.b_pf_pdb == 1)
    {
        check_extension_mpi(s.world_rank,"-pf_pdb",p.pf_pdb_file_name,".pdb");
    }
    if(p.b_pf_param == 1)
    {
        check_extension_mpi(s.world_rank,"-pf_prm",p.protein_finder_param_name,".prm");
    }

    //create index for each group of atoms
    Index group_1;

    //read the index files
    group_1.get_index(p.group_1_file_name);

    //run protein finder
    traj.get_protein(p.protein_finder_param_name,p.b_pf_param);

    //write pdb with protein indicated by beta factor
    traj.write_protein(p.pf_pdb_file_name,p.b_pf_pdb);

    //print info about the protein
    traj.get_prot_stats();

    int num_residues = traj.res_nr[traj.atoms()-1]; //How many residues does the system have
    iv1d contacts_sum(num_residues,0);              //Stores the number of contacts for each residue
    dv1d contacts(traj.prot.size(),0.0);            //Stores the number of contacts for each protein atom

    //allocate memory for average prot coords
    dv1d prot_x(traj.prot.size(),0.0);              //holds the protein x_coordinates
    dv1d prot_y(traj.prot.size(),0.0);              //holds the protein x_coordinates
    dv1d prot_z(traj.prot.size(),0.0);              //holds the protein x_coordinates
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

        prot_lip_contacts(traj,s,p,group_1,contacts,contacts_sum);

        prot_coords(traj,s,p,prot_x,prot_y,prot_z);

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect data from mpi processes and compute averages
    perf.log_time(finalize_analysis(traj,s,p,contacts,contacts_sum,num_residues,prot_x,prot_y,prot_z),"Fin Ana");

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

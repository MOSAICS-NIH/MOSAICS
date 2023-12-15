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
#include "MosAT/program_variables/pv_contact_rmsf.h"         //This has the variables specific to the analysis program
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
#include "headers/grid_io.h"                                 //This has routines used for reading in grid data
#include "headers/binding_events.h"                          //This has routines used for reading binding events files

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function stores the res_id for the bound lipid at each frame. Also checks if the be file and traj    //
// are compatible and sets the min, max, and resi_size arguments.                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_bound_lipids(Trajectory &traj,system_variables &s,program_variables &p,Binding_events &events,iv1d &bound_lipids)
{
    int i = 0;    //standard variable used in loops
    int j = 0;    //standard variable used in loops

    //check if the number of frames match
    if(traj.get_ef_frames() != events.ef_frames)
    {
        if(s.world_rank == 0)
        {
            printf("Frame count in the binding events file (%d) does not match the number of frames being analyzed in the trajectory (%d) \n",events.ef_frames,traj.get_ef_frames());
        }
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }

    //set the residue id for each frame that has a bound lipids
    //also set the min, max, and resi_size arguments
    for(i=0; i<events.lipid_nr.size(); i++) //loop over binding events
    {
        int min = traj.res_start[events.res_nr[i]-1];
        int max = traj.res_end[events.res_nr[i]-1];

        if(strcmp(traj.res_name[min].c_str(), p.target_res.c_str()) == 0)
        {
            p.min       = min;
            p.max       = max;
            p.resi_size = p.max - p.min + 1;

            for(j=events.bind_i[i]; j<=events.bind_f[i]; j++) //loop over bound frames
            {
                bound_lipids[j] = events.res_nr[i];
            }
        }
    }

    //check if any target lipids were found
    if(p.resi_size == 0)
    {
        if(s.world_rank == 0)
        {
            printf("No target residues (%s) were found in the binding events file (%s). \n",p.target_res.c_str(),p.be_file_name.c_str());
        }
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function searches for a  bound lipid in the current frame                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int find_lipid(Trajectory &traj,system_variables &s,program_variables &p,iv1d &bound_lipids,int global_frame)
{
    int i      = 0;   //standard variable used in loops
    int status = 0;   //tells if the current step contains a bound lipid or not

    if(bound_lipids[global_frame] != -1) //a lipid is bound
    {
        p.resi = bound_lipids[global_frame];
        p.min  = traj.res_start[p.resi-1];
        p.max  = traj.res_end[p.resi-1];
        status = 1;
    }

    return status;     
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the contact profiles and stores lipid/protein coordinates                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_contacts(Trajectory &traj,system_variables &s,program_variables &p,dv1d &x_coord,dv1d &y_coord,
                  dv1d &z_coord,dv1d &prot_x,dv1d &prot_y,dv1d &prot_z,FILE **tmp_file,Selection &this_sel,
		  iv1d &bound_lipids)
{
    int i = 0;                              //standard variable used in loops
    int j = 0;                              //standard variable used in loops

    if(p.b_be == 0 || find_lipid(traj,s,p,bound_lipids,traj.get_frame_global()) == 1)
    {
        //report frame number 
        fprintf(*tmp_file,"%s %d \n","frame",traj.get_frame_global());

        //Construct contact matrix
        for(i=p.min; i<=p.max; i++) //loop over target residue atoms
        {
            //report target atom number
            fprintf(*tmp_file," %d",i-p.min);

            for(j=0; j<this_sel.sel.size(); j++) //loop over the selected atoms
            {
                double dx = traj.r[i][0] - traj.r[this_sel.sel[j]-1][0];
                double dy = traj.r[i][1] - traj.r[this_sel.sel[j]-1][1];
                double dz = traj.r[i][2] - traj.r[this_sel.sel[j]-1][2];

                double dist = sqrt(dx*dx + dy*dy + dz*dz);

                if(dist < p.contact_cutoff)
                {
                    //report the protein atom number
                    fprintf(*tmp_file," %d",j);
                }
            }
            //add the end tag for the current target atom
            fprintf(*tmp_file," %s \n","end");
        }

        //store lipid coordinates 
        for(i=p.min; i<=p.max; i++) //loop over target residue atoms
        {
            x_coord[i-p.min] = x_coord[i-p.min] + traj.r[i][0];
            y_coord[i-p.min] = y_coord[i-p.min] + traj.r[i][1];
            z_coord[i-p.min] = z_coord[i-p.min] + traj.r[i][2];
        }

        //store protein coordinates
        for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
        {
            prot_x[i] = prot_x[i] + traj.r[traj.prot[i]-1][0];
            prot_y[i] = prot_y[i] + traj.r[traj.prot[i]-1][1];
            prot_z[i] = prot_z[i] + traj.r[traj.prot[i]-1][2];
        }

        p.norm_factor = p.norm_factor + 1;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Collect and report data                                                                                   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,dv1d &pf,dv1d &x_coord,
                         dv1d &y_coord,dv1d &z_coord,dv1d &prot_x,dv1d &prot_y,dv1d &prot_z,Selection &this_sel)
{
    int i           = 0;                              //standard variable used in loops
    int j           = 0;                              //standard variable used in loops
    int k           = 0;                              //standard variable used in loops
    int target_atom = -1;                             //store the target atom index for accessing data in tmp_contacts

    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nFinalizing analysis. This requires reading temporary contact data files and could take some time. \n");
    }

    //allocate memory for computing the rmsf
    dv2d contacts_mean(p.resi_size,dv1d(this_sel.sel.size(),0.0)); //stores the average contact profile for each target atom
    dv1d rmsf(p.resi_size,0.0);                                    //holds the sum of squares for each target atom
    iv1d contacts_global(this_sel.sel.size(),0);                   //holds the contacts profile for each target atom

    //collect the normalization factor 
    collect_and_sum_int(s.world_size,s.world_rank,&p.norm_factor);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Compute average contacts profile                                                                          //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //generate name of the temp contacts file
    string tag           = "_" + to_string(s.world_rank) + "_contacts.dat";
    string tmp_file_name = chop_and_add_tag(p.crmsf_file_name,tag);

    //make an index and read in contacts data
    Index tmp_contacts;
    tmp_contacts.get_index(tmp_file_name);

    //remove temp contacts file
    remove(tmp_file_name.c_str());

    for(j=0; j<tmp_contacts.index_s.size(); j++) //loop over items in contacts data file
    {
        if(strcmp("frame", tmp_contacts.index_s[j].c_str()) == 0) 
        {
            j = j + 2;
            target_atom = tmp_contacts.index_i[j];
        }
        else if(strcmp("end", tmp_contacts.index_s[j].c_str()) == 0)
        {
            if(j+1 < tmp_contacts.index_s.size()) //dont read past end of file
            {
                if(strcmp("frame", tmp_contacts.index_s[j+1].c_str()) == 0)
                {
                } 
                else 
                {
                    j = j + 1;
                    target_atom = tmp_contacts.index_i[j];
                } 
            }
        }
        else 
        {
            contacts_mean[target_atom][tmp_contacts.index_i[j]] = contacts_mean[target_atom][tmp_contacts.index_i[j]] + 1.0;
        }
    }

    //collect contacts mean
    for(i=0; i<contacts_mean.size(); i++) //loop over target atoms
    {
        collect_and_sum_dv1d(s.world_size,s.world_rank,contacts_mean[i]);
    }

    //normalize contacts mean
    if(s.world_rank == 0)
    { 
        //normalize the contacts mean
        for(i=0; i<p.resi_size; i++) //loop over the target atoms
        {
            for(j=0; j<this_sel.sel.size(); j++) //loop over the selected protein atoms
            {
                contacts_mean[i][j] = contacts_mean[i][j]/(double)p.norm_factor;
            }
        }
    }

    //distribute contacts mean
    for(i=0; i<contacts_mean.size(); i++) //loop over target atoms
    {
        broadcast_dv1d(s.world_size,s.world_rank,contacts_mean[i]);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Compute RMSF                                                                                              //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(j=0; j<tmp_contacts.index_s.size(); j++) //loop over items in contacts data file
    {
        if(strcmp("frame", tmp_contacts.index_s[j].c_str()) == 0)
        {
            j = j + 2;
            target_atom = tmp_contacts.index_i[j];
        }
        else if(strcmp("end", tmp_contacts.index_s[j].c_str()) == 0)
        {
            //compute delta
            double n_common = 0.0;
            double n_dif    = 0.0;
            double S        = 0.0;
            double delta    = 0.0;

            for(k=0; k<this_sel.sel.size(); k++) //loop over selected protein atoms
            {
                double largest  = 0.0; 
                double smallest = 0.0;

                if((double)contacts_global[k] > contacts_mean[target_atom][k])
                {
                    largest  = (double)contacts_global[k];
                    smallest = contacts_mean[target_atom][k];
                }
                else 
                {
                    largest  = contacts_mean[target_atom][k];
                    smallest = (double)contacts_global[k];
                }

                n_common = n_common + smallest; 
                n_dif    = n_dif    + largest - smallest; 
            }
         
            if(n_common + n_dif > 0.0)
            {
                S     = n_common/(n_common + n_dif);
                delta = 1.0 - S; 
                rmsf[target_atom] = rmsf[target_atom] + pow(delta,2.0);
            }
            else 
            {
                printf("Found zero for denominator \n");
            }
           
            //reset current profile 
            for(k=0; k<this_sel.sel.size(); k++) 
            {
                contacts_global[k] = 0;
            }

            //jump to next target atom or frame
            if(j+1 < tmp_contacts.index_s.size()) //dont read past end of file
            {
                if(strcmp("frame", tmp_contacts.index_s[j+1].c_str()) == 0)
                {
                }
                else
                {
                    j = j + 1;
                    target_atom = tmp_contacts.index_i[j];
                }
            }
        }
        else
        {
            contacts_global[tmp_contacts.index_i[j]] = 1;
        }            
    }

    //collect sum_o_squares data
    collect_and_sum_dv1d(s.world_size,s.world_rank,rmsf);

    if(s.world_rank == 0)
    { 
        for(i=0; i<p.resi_size; i++) //loop over target atoms
        {
            rmsf[i] = sqrt(rmsf[i]/(double)p.norm_factor);
        }
    }    

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Write data to PDB files                                                                                   //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(s.world_rank == 0)
    {
        printf("\nComputing time average coordinates over %d frames. \n",p.norm_factor);
    }

    //collect atomic coordinates from mpi cores   
    collect_and_sum_dv1d(s.world_size,s.world_rank,x_coord);
    collect_and_sum_dv1d(s.world_size,s.world_rank,y_coord);
    collect_and_sum_dv1d(s.world_size,s.world_rank,z_coord);
    collect_and_sum_dv1d(s.world_size,s.world_rank,prot_x);
    collect_and_sum_dv1d(s.world_size,s.world_rank,prot_y);
    collect_and_sum_dv1d(s.world_size,s.world_rank,prot_z);

    //project rmsf onto target residue atoms
    if(s.world_rank == 0)
    {
        //normalize target atom data
        for(i=0; i<x_coord.size(); i++) //loop over target atoms
        {
            x_coord[i] = x_coord[i]/(double)p.norm_factor;
            y_coord[i] = y_coord[i]/(double)p.norm_factor;
            z_coord[i] = z_coord[i]/(double)p.norm_factor;
        }

        //normalize protein atom data 
        for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
        {
            prot_x[i] = prot_x[i]/(double)p.norm_factor;
            prot_y[i] = prot_y[i]/(double)p.norm_factor;
            prot_z[i] = prot_z[i]/(double)p.norm_factor;
        }

        int this_size = prot_x.size() + x_coord.size();

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
        for(j=0; j<this_size; j++) //loop over atoms
        {
            this_atom_nr[j] = j+1;
        }

        //set the atom names
        for(j=0; j<traj.prot.size(); j++) //loop over protein atoms 
        {
            this_atom_name[j] = traj.atom_name[traj.prot[j]-1];
        }
        for(j=p.min; j<=p.max; j++) //loop over lipid atoms
        {
            this_atom_name[traj.prot.size()+j-p.min] = traj.atom_name[j];
        }

        //set the residue names
        for(j=0; j<traj.prot.size(); j++) //loop over protein atoms
        {
            this_res_name[j] = traj.res_name[traj.prot[j]-1];
        }
        for(j=p.min; j<=p.max; j++) //loop over lipid atoms
        {
            this_res_name[traj.prot.size()+j-p.min] = traj.res_name[j];
        }

        //set residue numbers
        for(j=0; j<traj.prot.size(); j++) //loop over protein atoms
        {
            this_res_nr[j] = traj.res_nr[traj.prot[j]-1];
        }
        for(j=p.min; j<=p.max; j++) //loop over lipid atoms
        {
            this_res_nr[traj.prot.size()+j-p.min] = this_res_nr[traj.prot.size()-1] + 1;
        }

        //allocate memory for the coordinates
        rvec *this_r;
        this_r = (rvec *)calloc(this_size , sizeof(*this_r));

        for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
        {
            this_r[i][0] = prot_x[i];
            this_r[i][1] = prot_y[i];
            this_r[i][2] = prot_z[i];
        }
        for(i=0; i<x_coord.size(); i++) //loop over target atoms
        {
            this_r[traj.prot.size()+i][0] = x_coord[i];
            this_r[traj.prot.size()+i][1] = y_coord[i];
            this_r[traj.prot.size()+i][2] = z_coord[i];
        }

        //set b-factor to the rmsf
        for(i=0; i<x_coord.size(); i++) //loop over target atoms
        {
            this_beta[traj.prot.size()+i] = rmsf[i]; 
        }
        for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
        {
            this_beta[i] = 0.0;
        }

        //open pdb file
        FILE *pdb_file;
        string pdb_file_name = chop_and_add_tag(p.crmsf_file_name,"_rmsf.pdb");
        pdb_file             = fopen(pdb_file_name.c_str(), "w");
        if(pdb_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",pdb_file_name.c_str());
        }
        else
        {
            write_frame_pdb(traj.box,this_size,this_atom_nr,this_res_nr,this_res_name,this_atom_name,this_r,traj.title,s.world_rank,&pdb_file,this_beta,this_weight,this_element,this_chain_id,i);
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
    s.program_name = "Contact RMSF";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                                       s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                                       s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                                      s.world_rank, s.cl_tags, &p.b_print,     0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                                        s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                                                   s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                                    s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                                      s.world_rank, s.cl_tags, &p.b_lsq,       0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                                        s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",                        s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_s(argc,argv,"-crmsf",  p.crmsf_file_name,            "Filename for output PDBs with contact RMSF projected onto the protein (pdb)",      s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_s(argc,argv,"-pf_pdb", p.pf_pdb_file_name,           "PDB file with selected protein (pdb) ",                                            s.world_rank, s.cl_tags, &p.b_pf_pdb,    0);
    add_argument_mpi_s(argc,argv,"-pf_prm", p.protein_finder_param_name,  "File with additional protein finder parameters (prm)",                             s.world_rank, s.cl_tags, &p.b_pf_param,  0);
    add_argument_mpi_i(argc,argv,"-resi",   &p.resi,                      "The target residue id",                                                            s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_d(argc,argv,"-cdist",  &p.contact_cutoff,            "The contact cutoff distance (nm)",                                                 s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_s(argc,argv,"-sel",    p.selection_text_file_name,   "Input file with the atom selection text (sel)",                                    s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_s(argc,argv,"-be",     p.be_file_name,               "Input binding events file (be)",                                                   s.world_rank, s.cl_tags, &p.b_be,        0);
    add_argument_mpi_i(argc,argv,"-x",      &p.target_x,                  "The target x lattice point used with the binding events file",                     s.world_rank, s.cl_tags, &p.b_x,         0);
    add_argument_mpi_i(argc,argv,"-y",      &p.target_y,                  "The target y lattice point used with the binding events file",                     s.world_rank, s.cl_tags, &p.b_y,         0);
    add_argument_mpi_s(argc,argv,"-type",   p.target_res,                 "Lipid type to select from the bound lipids timeline",                              s.world_rank, s.cl_tags, &p.b_target_res,0);
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
    //reset the clock
    s.t = clock();

    //check file extensions                                                                                     
    check_extension_mpi(s.world_rank,"-crmsf",p.crmsf_file_name,".pdb");
    check_extension_mpi(s.world_rank,"-sel",p.selection_text_file_name,".sel");

    if(p.b_be == 1)
    {
        check_extension_mpi(s.world_rank,"-be",p.be_file_name,".be");
    }
    if(p.b_pf_pdb == 1)
    {
        check_extension_mpi(s.world_rank,"-pf_pdb",p.pf_pdb_file_name,".pdb");
    }
    if(p.b_pf_param == 1)
    {
        check_extension_mpi(s.world_rank,"-pf_prm",p.protein_finder_param_name,".prm");
    }

    //check if everything nessesary for timeline was provided
    if(p.b_be == 1 || p.b_x == 1 || p.b_y == 1 || p.b_target_res == 1)
    {
        if(p.b_be == 0 || p.b_x == 0 || p.b_y == 0 || p.b_target_res == 0)
        {
            if(s.world_rank == 0)
            {
                printf("Arguments -be, -x, -y, and -type are used together to specify binding lipids dynamically. Please provide an argument for each of these or alternatively use the -resi argument alone. \n");
            }
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }
    }
 
    //run proten finder
    traj.get_protein(p.protein_finder_param_name,p.b_pf_param);

    //print a pdb with distinguished protein
    traj.write_protein(p.pf_pdb_file_name,p.b_pf_pdb);

    //print info about the target protein
    traj.get_prot_stats();
 
    //read in binding events
    Binding_events events;
    if(p.b_be == 1)
    { 
        int result = events.get_info(p.be_file_name); 
        if(result == 1)
        {
            result = events.get_binding_events_xy(p.be_file_name,p.target_x,p.target_y);
        }
        else 
        {
            if(s.world_rank == 0)
            {
                printf("Could not open binding events file %s \n",p.be_file_name.c_str());
            }
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }
    }

    //allocate memory to hold bound lipids 
    iv1d bound_lipids(traj.get_ef_frames(),-1);

    //get info about the target lipid or lipids
    if(p.b_be == 0)
    {
        p.min       = traj.res_start[p.resi-1];       //first atom of the target residue
        p.max       = traj.res_end[p.resi-1];         //last atom of the target residue
        p.resi_size = p.max - p.min + 1;              //how many atoms are in the target residue
    }
    else 
    {
        get_bound_lipids(traj,s,p,events,bound_lipids);
    }

    //allocate memory for RMSF data
    dv1d rmsf(traj.prot.size(),0.0);                //holds the contact rmsf data for each target atom 

    //allocate memory for target residue coords
    dv1d x_coord(p.resi_size,0.0);                    //holds the x coordinates 
    dv1d y_coord(p.resi_size,0.0);                    //holds the x coordinates 
    dv1d z_coord(p.resi_size,0.0);                    //holds the x coordinates 

    //allocate memory for average prot coords
    dv1d prot_x(traj.prot.size(),0.0);              //holds the protein x_coordinates
    dv1d prot_y(traj.prot.size(),0.0);              //holds the protein x_coordinates
    dv1d prot_z(traj.prot.size(),0.0);              //holds the protein x_coordinates

    //open files for writing temporary contact profiles data
    FILE *tmp_file;
    string tag           = "_" + to_string(s.world_rank) + "_contacts.dat";
    string tmp_file_name = chop_and_add_tag(p.crmsf_file_name,tag);
    tmp_file             = fopen(tmp_file_name.c_str(), "w");
    if(tmp_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",tmp_file_name.c_str());
    }

    //create a object to hold an atom selection
    Selection this_sel;

    //select the atoms
    this_sel.get_selection(traj,p.selection_text_file_name);

    //generate pdb file name for highlighting the selection 
    string pdb_filename = chop_and_add_tag(p.selection_text_file_name,".pdb");

    //highlight the selected atoms
    this_sel.highlight_sel(traj,pdb_filename);

    //log time spent preparing for analysis
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Other");

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

        get_contacts(traj,s,p,x_coord,y_coord,z_coord,prot_x,prot_y,prot_z,&tmp_file,this_sel,bound_lipids);

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //close tmp contacts files
    fclose(tmp_file);

    MPI_Barrier(MPI_COMM_WORLD);
 
    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //analyze and report data
    perf.log_time(finalize_analysis(traj,s,p,rmsf,x_coord,y_coord,z_coord,prot_x,prot_y,prot_z,this_sel),"Fin Ana");

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

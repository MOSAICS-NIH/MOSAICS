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
#include "MosAT/program_variables/pv_lipid_protein_min_dist.h"     //This has the variables specific to the analysis program
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
// This function computes the mean coords for each target atom and the protein                               //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_mean_coords(Trajectory &traj,system_variables &s,program_variables &p,dv1d &x_coord,dv1d &y_coord,
                     dv1d &z_coord,dv1d &prot_x,dv1d &prot_y,dv1d &prot_z,iv1d &bound_lipids)
{
    int i   = 0;                                                //standard variable used in loops

    if(p.b_be == 0 || find_lipid(traj,s,p,bound_lipids,traj.get_frame_global()) == 1)
    {  
        //store lipid coords
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
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the minimum distance between each target atom and the protein                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_min_dist(Trajectory &traj,system_variables &s,program_variables &p,dv1d &min_dist,dv1d &x_coord,dv1d &y_coord,
                  dv1d &z_coord,dv1d &prot_x,dv1d &prot_y,dv1d &prot_z,dv1d &contacts,dv1d &contacts_lip,dv1d &contacts_freq,
                  iv1d &bound_lipids)
{
    int i   = 0;                              //standard variable used in loops
    int j   = 0;                              //standard variable used in loops

    if(p.b_be == 0 || find_lipid(traj,s,p,bound_lipids,traj.get_frame_global()) == 1)
    {
        //Conpute the min dist 
        for(i=p.min; i<=p.max; i++) //loop over target residue atoms
        {
            int    num_contacts     = 0;
            int    num_contacts_lip = 0;
            double this_min_dist    = 999999.9;
 
            for(j=0; j<traj.prot.size(); j++) //loop over protein atoms
            {
                double dx = traj.r[i][0] - traj.r[traj.prot[j]-1][0];
                double dy = traj.r[i][1] - traj.r[traj.prot[j]-1][1];
                double dz = traj.r[i][2] - traj.r[traj.prot[j]-1][2];

                double dist = sqrt(dx*dx + dy*dy + dz*dz);

                if(dist < this_min_dist)
                {
                    this_min_dist = dist;
                }

                //count the number of contacts
                if(dist < p.contact_cutoff)
                {
                    num_contacts++;
                }
            }
            if(num_contacts > 0) //target atom is in contact with the protein
            {
                contacts_freq[i-p.min] = contacts_freq[i-p.min] + 1.0;
            }

            //count contacts with lipids
            for(j=0; j<traj.target_leaflet.size(); j++) //loop over target leaflet
            {
                if(traj.res_nr[i] != traj.res_nr[traj.target_leaflet[j]-1]) //dont count contacts with self
                {
                    double dx = traj.r[i][0] - traj.r[traj.target_leaflet[j]-1][0];
                    double dy = traj.r[i][1] - traj.r[traj.target_leaflet[j]-1][1];
                    double dz = traj.r[i][2] - traj.r[traj.target_leaflet[j]-1][2];

                    double dist = sqrt(dx*dx + dy*dy + dz*dz);

                    //count the number of contacts
                    if(dist < p.contact_cutoff)
                    {
                        num_contacts_lip++;
                    }
                }
            }

            min_dist[i-p.min]     = min_dist[i-p.min]     + this_min_dist;
            contacts[i-p.min]     = contacts[i-p.min]     + (double)num_contacts;
            contacts_lip[i-p.min] = contacts_lip[i-p.min] + (double)num_contacts_lip;
        }

        //increase the normalization factor
        p.norm_factor = p.norm_factor + 1;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function counts contacts between the protein atoms and other atoms                                   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_prot_contacts(Trajectory &traj,system_variables &s,program_variables &p,dv1d &prot_target,dv1d &prot_other,
                       iv1d &bound_lipids)
{
    int i            = 0;                              //standard variable used in loops
    int j            = 0;                              //standard variable used in loops
    int k            = 0;                              //standard variable used in loops
    int l            = 0;                              //standard variable used in loops
    int counter_prot = 0;                              //count protein residues as they are encountered
    int counter_lip  = 0;                              //count lipid residues as they are encountered

    if(p.b_be == 0 || find_lipid(traj,s,p,bound_lipids,traj.get_frame_global()) == 1)
    {
        //Count contacts with target resi 
        for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
        {
            int num_contacts       = 0;
            int num_contacts_other = 0;

            //count contacts with the target residue
            for(j=p.min; j<=p.max; j++) //loop over target residue atoms
            {
                double dx = traj.r[traj.prot[i]-1][0] - traj.r[j][0];
                double dy = traj.r[traj.prot[i]-1][1] - traj.r[j][1];
                double dz = traj.r[traj.prot[i]-1][2] - traj.r[j][2];

                double dist = sqrt(dx*dx + dy*dy + dz*dz);

                if(dist < p.contact_cutoff)
                {
                    num_contacts++;
                }
            }
            prot_target[i] = prot_target[i] + (double)num_contacts;
        }

        //compute center of mass of lipid and protein residues
        dv2d lipid_centers   = traj.get_centers_target_lf();
        dv2d protein_centers = traj.get_centers_prot();

        //compute number of contacts with non-target lipids
        for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
        {
            int min_p = traj.p_res_start(i);  
            int max_p = traj.p_res_end(i);

            counter_lip = 0;

            //count contacts with lipids excluding the target residue
            for(j=0; j<traj.target_leaflet.size(); j++) //loop over target leaflet
            {
                //get the first and last atom of the current lipid
                int min_l = traj.t_lip_start(j);
                int max_l = traj.t_lip_end(j);

                //jump to the next lipid
                j = traj.next_target_lipid(j);

                double dx_center = protein_centers[counter_prot][0] - lipid_centers[counter_lip][0];
                double dy_center = protein_centers[counter_prot][1] - lipid_centers[counter_lip][1];
                double dz_center = protein_centers[counter_prot][2] - lipid_centers[counter_lip][2];

                double dist_center = sqrt(dx_center*dx_center + dy_center*dy_center + dz_center*dz_center);

                if(dist_center < p.screen_dist)
                {
                    if(traj.res_nr[min_l] != p.resi) //dont count contacts with target resi
                    {
                        for(k=min_p; k<=max_p; k++) //loop over current residue atoms
                        {
                            for(l=min_l; l<=max_l; l++) //loop over current lipid
                            {     
                                double dx = traj.r[k][0] - traj.r[l][0];
                                double dy = traj.r[k][1] - traj.r[l][1];
                                double dz = traj.r[k][2] - traj.r[l][2];

                                double dist = sqrt(dx*dx + dy*dy + dz*dz);

                                if(dist < p.contact_cutoff)
                                {
                                    prot_other[i+k-min_p] = prot_other[i+k-min_p] + 1.0;
                                }
                            }
                        } 
                    }
                }
                counter_lip++;
            }

            //jump to the next residue
            i = traj.next_prot_res(i);

            counter_prot++;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Collect and report data                                                                                   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,dv1d &min_dist,dv1d &x_coord,
                         dv1d &y_coord,dv1d &z_coord,dv1d &prot_x,dv1d &prot_y,dv1d &prot_z,dv1d &contacts,dv1d &contacts_lip,
                         dv1d &prot_target,dv1d &prot_other,dv1d &contacts_freq)
{
    int i   = 0;                              //standard variable used in loops
    int j   = 0;                              //standard variable used in loops

    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nFinalizing analysis. This requires communicating data and could take some time. \n");
    }

    //collect the normalization factor 
    collect_and_sum_int(s.world_size,s.world_rank,&p.norm_factor);

    //collect min dist data from mpi cores   
    collect_and_sum_dv1d(s.world_size,s.world_rank,min_dist);
    collect_and_sum_dv1d(s.world_size,s.world_rank,contacts);
    collect_and_sum_dv1d(s.world_size,s.world_rank,contacts_lip);
    collect_and_sum_dv1d(s.world_size,s.world_rank,prot_target);
    collect_and_sum_dv1d(s.world_size,s.world_rank,prot_other);
    collect_and_sum_dv1d(s.world_size,s.world_rank,x_coord);
    collect_and_sum_dv1d(s.world_size,s.world_rank,y_coord);
    collect_and_sum_dv1d(s.world_size,s.world_rank,z_coord);
    collect_and_sum_dv1d(s.world_size,s.world_rank,prot_x);
    collect_and_sum_dv1d(s.world_size,s.world_rank,prot_y);
    collect_and_sum_dv1d(s.world_size,s.world_rank,prot_z);
    collect_and_sum_dv1d(s.world_size,s.world_rank,contacts_freq);

    //write summed contact matrix to file
    if(s.world_rank == 0)
    {
        printf("\nComputing time average coordinates over %d frames. \n",p.norm_factor);

    	//normalize target atom data
        for(i=0; i<min_dist.size(); i++) //loop over target atoms
        {
            min_dist[i]      = min_dist[i]/(double)p.norm_factor;
            contacts[i]      = contacts[i]/(double)p.norm_factor;
            contacts_lip[i]  = contacts_lip[i]/(double)p.norm_factor;
            x_coord[i]       = x_coord[i]/(double)p.norm_factor;
            y_coord[i]       = y_coord[i]/(double)p.norm_factor;
            z_coord[i]       = z_coord[i]/(double)p.norm_factor;
            contacts_freq[i] = contacts_freq[i]/(double)p.norm_factor;
        }

        //normalize protein atom data 
        for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
        {
            prot_target[i] = prot_target[i]/(double)p.norm_factor;
            prot_other[i]  = prot_other[i]/(double)p.norm_factor;
            prot_x[i]      = prot_x[i]/(double)p.norm_factor;
            prot_y[i]      = prot_y[i]/(double)p.norm_factor;
            prot_z[i]      = prot_z[i]/(double)p.norm_factor;
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

        //set b-factor to the min dist
        for(i=0; i<x_coord.size(); i++) //loop over target atoms
        {
            this_beta[traj.prot.size()+i] = min_dist[i]; 
        }

        //open pdb file
        FILE *pdb_file;
        string pdb_file_name = chop_and_add_tag(p.lpmd_file_name,"_min_dist.pdb");
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

        //set b-factor to the number of contacts
        for(i=0; i<x_coord.size(); i++) //loop over target atoms
        {
            this_beta[traj.prot.size()+i] = contacts[i];
        }

        pdb_file_name = chop_and_add_tag(p.lpmd_file_name,"_num_contacts.pdb");
        pdb_file      = fopen(pdb_file_name.c_str(), "w");
        if(pdb_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",pdb_file_name.c_str());
        }
        else
        {
            write_frame_pdb(traj.box,this_size,this_atom_nr,this_res_nr,this_res_name,this_atom_name,this_r,traj.title,s.world_rank,&pdb_file,this_beta,this_weight,this_element,this_chain_id,i);
            fclose(pdb_file);
        }

        //set b-factor to the number of lipid contacts
        for(i=0; i<x_coord.size(); i++) //loop over target atoms
        {
            this_beta[traj.prot.size()+i] = contacts_lip[i];
        }

        pdb_file_name = chop_and_add_tag(p.lpmd_file_name,"_num_contacts_lip.pdb");
        pdb_file      = fopen(pdb_file_name.c_str(), "w");
        if(pdb_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",pdb_file_name.c_str());
        }
        else
        {
            write_frame_pdb(traj.box,this_size,this_atom_nr,this_res_nr,this_res_name,this_atom_name,this_r,traj.title,s.world_rank,&pdb_file,this_beta,this_weight,this_element,this_chain_id,i);
            fclose(pdb_file);
        }

        //set b-factor to the number of protein contacts with target 
        for(i=0; i<x_coord.size(); i++) //loop over target atoms
        {
            this_beta[traj.prot.size()+i] = 0.0;
        }
        for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
        {
            this_beta[i] = prot_target[i];
        }

        pdb_file_name = chop_and_add_tag(p.lpmd_file_name,"_num_contacts_prot_target.pdb");
        pdb_file      = fopen(pdb_file_name.c_str(), "w");
        if(pdb_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",pdb_file_name.c_str());
        }
        else
        {
            write_frame_pdb(traj.box,this_size,this_atom_nr,this_res_nr,this_res_name,this_atom_name,this_r,traj.title,s.world_rank,&pdb_file,this_beta,this_weight,this_element,this_chain_id,i);
            fclose(pdb_file);
        }

        //set b-factor to the number of protein contacts with other lipids 
        for(i=0; i<x_coord.size(); i++) //loop over target atoms
        {
            this_beta[traj.prot.size()+i] = 0.0;
        }
        for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
        {
            this_beta[i] = prot_other[i];
        }

        pdb_file_name = chop_and_add_tag(p.lpmd_file_name,"_num_contacts_prot_other.pdb");
        pdb_file      = fopen(pdb_file_name.c_str(), "w");
        if(pdb_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",pdb_file_name.c_str());
        }
        else
        {
            write_frame_pdb(traj.box,this_size,this_atom_nr,this_res_nr,this_res_name,this_atom_name,this_r,traj.title,s.world_rank,&pdb_file,this_beta,this_weight,this_element,this_chain_id,i);
            fclose(pdb_file);
        }

        //set b-factor to the fraction of frames the target atom was in contact with the protein
        for(i=0; i<x_coord.size(); i++) //loop over target atoms
        {
            this_beta[traj.prot.size()+i] = contacts_freq[i];
        }
        for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
        {
            this_beta[i] = 0.0;
        }

        pdb_file_name = chop_and_add_tag(p.lpmd_file_name,"_contacts_freq.pdb");
        pdb_file      = fopen(pdb_file_name.c_str(), "w");
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

    //distribute average lipid coords to mpi world
    broadcast_dv1d(s.world_size,s.world_rank,x_coord);
    broadcast_dv1d(s.world_size,s.world_rank,y_coord);
    broadcast_dv1d(s.world_size,s.world_rank,z_coord);

    MPI_Barrier(MPI_COMM_WORLD);

    //compute and return time spent in function
    return (clock() - s.t)/CLOCKS_PER_SEC;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the RMSF for each target atom                                                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_rmsf(Trajectory &traj,system_variables &s,program_variables &p,dv1d &x_coord,dv1d &y_coord,dv1d &z_coord,
              dv1d &rmsf,iv1d &bound_lipids)
{
    int i   = 0;                              //standard variable used in loops

    if(p.b_be == 0 || find_lipid(traj,s,p,bound_lipids,traj.get_frame_global()) == 1)
    {
        //Conpute the min dist 
        for(i=p.min; i<=p.max; i++) //loop over target residue atoms
        {
            double dx = traj.r[i][0] - x_coord[i-p.min];
            double dy = traj.r[i][1] - y_coord[i-p.min];
            double dz = traj.r[i][2] - z_coord[i-p.min];

            double delta = dx*dx + dy*dy + dz*dz;

            rmsf[i-p.min] = rmsf[i-p.min] + delta;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Analyze and report RMSF data                                                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_rmsf(Trajectory &traj,system_variables &s,program_variables &p,dv1d &rmsf,dv1d &x_coord,
                     dv1d &y_coord,dv1d &z_coord,dv1d &prot_x,dv1d &prot_y,dv1d &prot_z)
{
    int i   = 0;                              //standard variable used in loops
    int j   = 0;                              //standard variable used in loops

    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("Finalizing RMSF calculation. This requires communicating RMSF data and could take some time. \n");
    }

    collect_and_sum_dv1d(s.world_size,s.world_rank,rmsf);

    //write summed contact matrix to file
    if(s.world_rank == 0)
    {
        //normalize RMSF data
        for(i=0; i<rmsf.size(); i++) //loop over target atoms
        {
            rmsf[i] = sqrt(rmsf[i]/(double)p.norm_factor);
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

        //set b-factor to the RMSF
        for(i=0; i<x_coord.size(); i++) //loop over target atoms
        {
            this_beta[traj.prot.size()+i] = rmsf[i];
        }

        FILE *pdb_file;
        string pdb_file_name = chop_and_add_tag(p.lpmd_file_name,"_rmsf.pdb");
        pdb_file      = fopen(pdb_file_name.c_str(), "w");
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
    s.program_name = "Lipid Protein Min Dist";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                                s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                                s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                               s.world_rank, s.cl_tags, &p.b_print,     0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                                 s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                                            s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                             s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                               s.world_rank, s.cl_tags, &p.b_lsq,       0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                                 s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",                 s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_s(argc,argv,"-lpmd",   p.lpmd_file_name,             "Base filename for output PDBs with data projected onto the lipids (pdb)",   s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_s(argc,argv,"-pf_pdb", p.pf_pdb_file_name,           "PDB file with selected protein (pdb) ",                                     s.world_rank, s.cl_tags, &p.b_pf_pdb,    0);
    add_argument_mpi_s(argc,argv,"-pf_prm", p.protein_finder_param_name,  "File with additional protein finder parameters (prm)",                      s.world_rank, s.cl_tags, &p.b_pf_param,  0);
    add_argument_mpi_d(argc,argv,"-cdist",  &p.contact_cutoff,            "The contact cutoff distance (nm)",                                          s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_i(argc,argv,"-resi",   &p.resi,                      "The target residue id",                                                     s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,           "PDB file with sorted leaflets (pdb)",                                       s.world_rank, s.cl_tags, &p.b_lf_pdb,    0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,  "File with additional leaflet finder parameters (prm) ",                     s.world_rank, s.cl_tags, &p.b_lf_param,  0);
    add_argument_mpi_i(argc,argv,"-leaf",   &p.leaflet,                   "Which leaflet? (0:both 1:upper 2:lower)",                                   s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_d(argc,argv,"-screen", &p.screen_dist,               "Screen lipids whose centers are within this disatnce (nm) of the protein",  s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_s(argc,argv,"-be",     p.be_file_name,               "Input binding events file (be)",                                            s.world_rank, s.cl_tags, &p.b_be,        0);
    add_argument_mpi_i(argc,argv,"-x",      &p.target_x,                  "The target x lattice point used with the binding events file",              s.world_rank, s.cl_tags, &p.b_x,         0);
    add_argument_mpi_i(argc,argv,"-y",      &p.target_y,                  "The target y lattice point used with the binding events file",              s.world_rank, s.cl_tags, &p.b_y,         0);
    add_argument_mpi_s(argc,argv,"-type",   p.target_res,                 "Lipid type to select from the bound lipids timeline",                       s.world_rank, s.cl_tags, &p.b_target_res,0);
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
    check_extension_mpi(s.world_rank,"-lpmd",p.lpmd_file_name,".pdb");

    if(p.b_be == 1)
    {
        check_extension_mpi(s.world_rank,"-be",p.be_file_name,".be");
    }
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

    //check if everything nessesary for timeline was provided
    if(p.b_be == 1 || p.b_target_res == 1)
    {
        if(p.b_be == 0 || p.b_target_res == 0)
        {
            if(s.world_rank == 0)
            {
                printf("Arguments -be, and -type are used together to specify binding lipids dynamically. Please provide an argument for each of these or alternatively use the -resi argument alone. \n");
            }
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }
        else
        {
            if(p.b_x == 1 && p.b_y == 0)
            {
                printf("A lattice point was specified for the x-direction but not y. Please include the y-direction if analyzing binding events for a lattice point. \n");
                MPI_Finalize();
                exit(EXIT_SUCCESS);
            }
            else if(p.b_x == 0 && p.b_y == 1)
            {
                printf("A lattice point was specified for the y-direction but not x. Please include the x-direction if analyzing binding events for a lattice point. \n");
                MPI_Finalize();
                exit(EXIT_SUCCESS);
            }
        }
    }

    //run leaflet/proten finder
    traj.get_leaflets(p.leaflet,p.leaflet_finder_param_name,p.b_lf_param);
    traj.get_protein(p.protein_finder_param_name,p.b_pf_param);

    //print a pdb with distinguished leaflets/protein
    traj.write_leaflets(p.lf_pdb_file_name,p.b_lf_pdb);
    traj.write_protein(p.pf_pdb_file_name,p.b_pf_pdb);

    //print info about the target leaflet/protein
    traj.get_leaflet_stats();
    traj.get_prot_stats();

    //read in binding events
    Binding_events target_events;
    if(p.b_be == 1)
    {
        int result = 0;
        if(p.b_x==1 && p.b_y==1)
        {
            result = target_events.get_info(p.be_file_name);
            if(result == 1)
            {
                result = target_events.get_binding_events_xy(p.be_file_name,p.target_x,p.target_y);
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
        else
        {
            result = target_events.get_binding_events_bin(p.be_file_name);
            if(result == 0)
            {
                if(s.world_rank == 0)
                {
                    printf("Could not open binding events file %s \n",p.be_file_name.c_str());
                }
                MPI_Finalize();
                exit(EXIT_SUCCESS);
            }
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
        get_bound_lipids(traj,s,p,target_events,bound_lipids);
    }

    //allocate memory for min dist and contacts data
    dv1d min_dist(p.resi_size,0.0);                 //holds the min dist data for each target atom 
    dv1d contacts(p.resi_size,0.0);                 //holds the number of contacts for each target atom
    dv1d contacts_lip(p.resi_size,0.0);             //holds the number of contacts for each target atom
    dv1d prot_target(traj.prot.size(),0.0);         //holds the number of contacts with target atoms for each protein atom 
    dv1d prot_other(traj.prot.size(),0.0);          //holds the number of contacts with other lipid lipids (not the target resi) and each protein atom
    dv1d contacts_freq(p.resi_size,0.0);            //holds the frequency that a contact is formed with the protein for each target atom

    //allocate memory for target residue coords
    dv1d x_coord(p.resi_size,0.0);                  //holds the x coordinates 
    dv1d y_coord(p.resi_size,0.0);                  //holds the x coordinates 
    dv1d z_coord(p.resi_size,0.0);                  //holds the x coordinates 

    //allocate memory for average prot coords
    dv1d prot_x(traj.prot.size(),0.0);              //holds the protein x_coordinates
    dv1d prot_y(traj.prot.size(),0.0);              //holds the protein x_coordinates
    dv1d prot_z(traj.prot.size(),0.0);              //holds the protein x_coordinates

    //allocate memory to hold rmsf for the lipid
    dv1d rmsf(p.resi_size,0.0);
 
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

        get_min_dist(traj,s,p,min_dist,x_coord,y_coord,z_coord,prot_x,prot_y,prot_z,contacts,contacts_lip,contacts_freq,bound_lipids);

        get_prot_contacts(traj,s,p,prot_target,prot_other,bound_lipids);

        get_mean_coords(traj,s,p,x_coord,y_coord,z_coord,prot_x,prot_y,prot_z,bound_lipids);

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //analyze and report data
    perf.log_time(finalize_analysis(traj,s,p,min_dist,x_coord,y_coord,z_coord,prot_x,prot_y,prot_z,contacts,contacts_lip,prot_target,prot_other,contacts_freq),"Fin Ana");

    MPI_Barrier(MPI_COMM_WORLD);

    //begin calculation of RMSF
    if(s.world_rank == 0)
    {	     
        printf("\nBeginning Calculation of RMSF for target residue. \n\n");
    }

    //print info about the worlk load distribution
    traj.workload();

    s.t = clock();
    //compute the rmsf for the lipid
    for(traj.current_frame=0; traj.current_frame<traj.get_num_frames(); traj.current_frame++)
    {
        traj.read_traj_frame();

        traj.do_fit();

        get_rmsf(traj,s,p,x_coord,y_coord,z_coord,rmsf,bound_lipids);

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent computing RMSF
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"RMSF");

    //analyze and report RMSF data
    perf.log_time(finalize_rmsf(traj,s,p,rmsf,x_coord,y_coord,z_coord,prot_x,prot_y,prot_z),"Fin RMSF");

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

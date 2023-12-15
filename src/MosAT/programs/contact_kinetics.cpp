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
#include "MosAT/program_variables/pv_contact_kinetics.h"     //This has the variables specific to the analysis program
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
#include "headers/parallel.h"                                //This has routines used for distributing the workload over various coordinates
#include "headers/histo.h"                                   //This has routines used for making a histogram
#include "headers/contacts.h"                                //This has routines used for reading and writing contacts matrices data files
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
// This function reports information about the contact                                                       //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void report_contact_info(Trajectory &traj,system_variables &s,program_variables &p,int contact_count,int atom_number_1,int atom_number_2,double distance)
{
    printf("select atom_a, (resi %d & resn %s & name %s) \n",traj.res_nr[atom_number_1-1]%10000,traj.res_name[atom_number_1-1].c_str(),traj.atom_name[atom_number_1-1].c_str());
    printf("select atom_b, (resi %d & resn %s & name %s) \n",traj.res_nr[atom_number_2-1]%10000,traj.res_name[atom_number_2-1].c_str(),traj.atom_name[atom_number_2-1].c_str());
    printf("select %d, (resi %d & resn %s & name %s) or ((resi %d & resn %s & name %s)) \n",contact_count,traj.res_nr[atom_number_1-1]%10000,traj.res_name[atom_number_1-1].c_str(),traj.atom_name[atom_number_1-1].c_str(),traj.res_nr[atom_number_2-1]%10000,traj.res_name[atom_number_2-1].c_str(),traj.atom_name[atom_number_2-1].c_str());
    printf("dist (atom_a), (atom_b) \n");
    printf("distance %f \n",distance);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function lets you check the mapping of contact data onto the time average coords                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void check_mapping(string lipid_atom_name,int protein_atom_nr,double *dist,Trajectory &traj,int i,int j)
{
    if(strcmp(traj.atom_name[i].c_str(), lipid_atom_name.c_str()) == 0 && traj.prot[j] == protein_atom_nr )
    {
        *dist = 0.0;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reports the frame number when a target contact is present                                   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void report_contact_frame(Trajectory &traj,system_variables &s,program_variables &p,Index &report_pairs,int atom_nr_1,int atom_nr_2)
{
    int i = 0;  //standard variable used in loops

    if(p.b_report == 1)
    {
        for(i=0; i<report_pairs.index_s.size(); i+=2)
        {
            if( (atom_nr_1 == report_pairs.index_i[i] && atom_nr_2 == report_pairs.index_i[i+1]) || (atom_nr_1 == report_pairs.index_i[i+1] && atom_nr_2 == report_pairs.index_i[i]) )
            {
                printf("Contact identified at trajectory frame %10d between %s %d (%d) atom %s %d (%d) and %s %d (%d) atom %s %d (%d). \n",traj.get_frame_full(),traj.res_name[atom_nr_1-1].c_str(),traj.res_nr[atom_nr_1-1],traj.res_nr[atom_nr_1-1]%10000,traj.atom_name[atom_nr_1-1].c_str(),traj.atom_nr[atom_nr_1-1],traj.atom_nr[atom_nr_1-1]%100000,traj.res_name[atom_nr_2-1].c_str(),traj.res_nr[atom_nr_2-1],traj.res_nr[atom_nr_2-1]%10000,traj.atom_name[atom_nr_2-1].c_str(),traj.atom_nr[atom_nr_2-1],traj.atom_nr[atom_nr_2-1]%100000); 
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the contact matrix for the current trajectory frame                                //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void binding_contacts(Trajectory &traj,system_variables &s,program_variables &p,Index &report_pairs,
		      Contacts &cont,iv1d &bound_lipids)
{
    //Note: We add a profile even if a lipid is not bound in the timeline. in this case it will be empty but is
    //Note: needed for reading back in contact profiles since a global frame number is used to lookup a contact profile.

    int i             = 0;                                                //standard variable used in loops
    int j             = 0;                                                //standard variable used in loops
    int k             = 0;                                                //standard variable used in loops
    int contact_count = 0;                                                //count the contacts as they are encountered

    iv2d this_profile(0,iv1d(0,0));        //hold contact info for the current frame

    if(p.b_be == 0 || find_lipid(traj,s,p,bound_lipids,traj.get_frame_global()) == 1)
    {
    	for(i=p.min; i<=p.max; i++) //loop over target residue atoms
        {
            for(j=0; j<traj.prot.size(); j++) //loop over protein atoms
            {
                double dx = traj.r[i][0] - traj.r[traj.prot[j]-1][0];
                double dy = traj.r[i][1] - traj.r[traj.prot[j]-1][1];
                double dz = traj.r[i][2] - traj.r[traj.prot[j]-1][2];

                double dist = sqrt(dx*dx + dy*dy + dz*dz);

                //check the mapping of contact data onto time average coords
                //check_mapping("H18R",1000,&dist,traj,i,j);

                if(dist < p.contact_cutoff)
                {
                    iv1d this_contact(2,0);      //holds atom id's for the current contact
                    this_contact[0] = i-p.min;
                    this_contact[1] = j;

                    this_profile.push_back(this_contact);

                    if(p.test == 1)
                    {
                        int atom_number_1 = traj.atom_nr[i];
                        int atom_number_2 = traj.prot[j];

                        report_contact_info(traj,s,p,contact_count,atom_number_1,atom_number_2,dist);
                    }

                    report_contact_frame(traj,s,p,report_pairs,traj.atom_nr[i],traj.atom_nr[traj.prot[j]-1]);

                    contact_count++;
                }
            }
        }

        //increase the normalization factor
        p.norm_factor = p.norm_factor + 1;
    }

    //add the current profile to tmp files
    cont.add_profile(this_profile,traj.get_frame_global());
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
    int j   = 0;                                                //standard variable used in loops

    if(p.b_be == 0 || find_lipid(traj,s,p,bound_lipids,traj.get_frame_global()) == 1)
    {
        //Conpute the min dist 
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
// Analyze contact matrices and compile data                                                                 //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,dv1d &x_coord,dv1d &y_coord,
		         dv1d &z_coord,dv1d &prot_x,dv1d &prot_y,dv1d &prot_z,iv3d &filter,iv2d &bound,dv2d &dwell_t,
			 dv2d &events,Contacts &cont,iv3d &filter_history,iv1d &bound_lipids)
{
    //Note: The strategy used here is to check for a new bound lipid using the timeline and then reinitialize the filter to 0 when a new lipid is found.
    //Note: this approach ensures that dwell times are recorded when a new lipid binds since the old contacts are kicked off. 
    //Note: Since a contact profile is made for every frame (some are empty if a target lipid was not bound) we do not need to screen for a bound lipid on 
    //Note: the main loop that reads and merges temporary contact profiles.  

    int i         = 0;                                          //standard variable used in loops
    int j         = 0;                                          //standard variable used in loops
    int k         = 0;                                          //standard variable used in loops
    int l         = 0;                                          //standard variable used in loops
    int prev_resi = -1;                                         //keep track of previous resi so you know when a new one is bound

    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    //collect the normalization factor 
    collect_and_sum_int(s.world_size,s.world_rank,&p.norm_factor);

    //merge tmp contact matrices to a single file
    cont.merge_profiles();

    //print info about the worlk load distribution
    traj.workload_prot();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Scan over trajectory, filter data, and record residence times                                             //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(s.world_rank == 0)
    { 
        printf("Computing contact lifetimes \n");
        printf("-----------------------------------------------------------------------------------------------------------------------------------\n");
    }

    //reset counter for printing progress + time estimates
    s.counter = 0;

    //get the time for reporting progress
    clock_t this_time = clock();

    for(i=0; i<traj.get_ef_frames(); i++) //loop over analyzed frames
    {
        int pos = i%(2*p.range + 1); //Tells current position in the noise filter

        //check for a new bound lipid
        if(p.b_be == 1 && find_lipid(traj,s,p,bound_lipids,i) == 1 && p.resi != prev_resi)
        {
            //reset the filter
            for(j=0; j<p.size_y; j++) //loop over y
            {
                for(k=traj.prot_start; k<=traj.prot_end; k++) //loop over my protein atoms
                {
                    for(l=0; l<(2*p.range + 1); l++) //loop over flanking frames
                    {
                        filter[l][j][k] = 0;
                    }
                }
            }

            //store the current resi
            prev_resi = p.resi;
        }

        //clear the current matrix
        cont.clear_profile(filter_history[pos],filter[pos]);

        //get the contact matrix for current frame
        filter_history[pos] = cont.get_profile_alt(i,filter[pos]);

        if(i >= (2*p.range)) //filter is full, compute frequencies
        {
            for(j=0; j<p.size_y; j++) //loop over y
            {
                for(k=traj.prot_start; k<=traj.prot_end; k++) //loop over my protein atoms
                {
                    double freq = 0.0;
                    
                    for(l=0; l<(2*p.range + 1); l++) //loop over flanking frames
                    {
                        if(filter[l][j][k] == 1)
                        {
                            freq = freq + 1.0;
                        }
                    }
                    freq = freq/((double)(2*p.range + 1));

                    if(freq >= 0.5 && (p.dump == 0 || i < traj.get_ef_frames()-1) )
                    {
                        bound[j][k-traj.prot_start] = bound[j][k-traj.prot_start] + 1;
                    }
                    else //contact is not formed
                    {
                        if(bound[j][k-traj.prot_start] > 0)
                        {
                            int time                      = bound[j][k-traj.prot_start];
                            dwell_t[k-traj.prot_start][j] = dwell_t[k-traj.prot_start][j] + (double)time*p.delta_t; 
                            events [k-traj.prot_start][j] = events [k-traj.prot_start][j] + 1.0; 
                        }
                        bound[j][k-traj.prot_start] = 0;
                    } 
                }
            }
        }

        //report progress
        time_stats(this_time,&s.counter,i,traj.get_ef_frames(),s.world_rank);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Collect residence time data                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(s.world_rank == 0)
    {
        printf("\nCollecting kinetics data and atomic coordinates from mpi cores\n");
    }

    //collect data
    collect_dv2d(s.world_size,s.world_rank,dwell_t);
    collect_dv2d(s.world_size,s.world_rank,events);

    //collect atomic coordinates from mpi cores   
    collect_and_sum_dv1d(s.world_size,s.world_rank,x_coord);
    collect_and_sum_dv1d(s.world_size,s.world_rank,y_coord);
    collect_and_sum_dv1d(s.world_size,s.world_rank,z_coord);
    collect_and_sum_dv1d(s.world_size,s.world_rank,prot_x);
    collect_and_sum_dv1d(s.world_size,s.world_rank,prot_y);
    collect_and_sum_dv1d(s.world_size,s.world_rank,prot_z);

    if(s.world_rank == 0)
    {
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Normalize data and extract relevant contacts                                                              //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        printf("\nNormalizing kinetics data and atomic coordinates. \n");

        //allocate memory for select contacts
        dv1d winners(0,0.0);          //holds the residence time for contacts with mean dwell time greater than zero
        dv1d winners_events(0,0.0);   //holds the number of binding events contributing to the mean dwell time
        iv1d atom_a(0,0);             //holds the atom number for atom a in the PDB with time average coordinates
        iv1d atom_b(0,0);             //holds the atom number for atom b in the PDB with time average coordinates
        iv1d atom_nr_a(0,0);          //holds the atom number for atom a in the reference file
        iv1d atom_nr_b(0,0);          //holds the atom number for atom b in the reference file

        //normalize data
        for(i=0; i<p.size_y; i++) //loop over y (target resi atoms)
        {
            for(j=0; j<p.size_x; j++) //loop over x (protein atoms)
            {
                if(events[j][i] > 0.0)
                {
                    dwell_t[j][i] = dwell_t[j][i]/events[j][i];
                }
                else 
                {
                    dwell_t[j][i] = 0.0;
                }

                if(dwell_t[j][i] > 0.0)
                {
                    winners.push_back(dwell_t[j][i]);
                    winners_events.push_back(events[j][i]);
                    atom_a.push_back(p.size_x + i + 1);
                    atom_b.push_back(j+1);
                    atom_nr_a.push_back(traj.atom_nr[p.min + i]);
                    atom_nr_b.push_back(traj.atom_nr[traj.prot[j]-1]);
                }
            }
        }

        printf("\nComputing time average coordinates over %d frames. \n",p.norm_factor);

        //normalize mean coords data
        for(i=0; i<x_coord.size(); i++) //loop over target atoms
        {
            x_coord[i] = x_coord[i]/(double)p.norm_factor;
            y_coord[i] = y_coord[i]/(double)p.norm_factor;
            z_coord[i] = z_coord[i]/(double)p.norm_factor;
        }

        //normalize protein coords
        for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
        {
            prot_x[i] = prot_x[i]/(double)p.norm_factor;
            prot_y[i] = prot_y[i]/(double)p.norm_factor;
            prot_z[i] = prot_z[i]/(double)p.norm_factor;
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Analyze dwell time distribution                                                                           //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        printf("\nReporting contact lifetime distribution \n");

        dv1d dwell_t_histo(0,0.0);

        for(i=0; i<p.size_y; i++) //loop over y
        {
            for(j=0; j<p.size_x; j++) //loop over x
            {
                if(dwell_t[j][i] > 0.0)
                {
                    dwell_t_histo.push_back(dwell_t[j][i]);                
                }
            }
        }

        string histo_file_name = chop_and_add_tag(p.ck_file_name,"_histo.dat");

        Histogram_d histo;
        histo.bin_data(dwell_t_histo,p.bin_width);
        histo.write_histo(histo_file_name,"Residence time (ps)");

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Write time average coordinates to file                                                                    //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

        //allocate memory for the curent cluster
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

        //open pdb file
        FILE *pdb_file;
        string pdb_file_name = chop_and_add_tag(p.ck_file_name,".pdb");
        pdb_file             = fopen(pdb_file_name.c_str(), "w");
        if(pdb_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",pdb_file_name.c_str());
        }
        else
        {
            printf("\nWriting time-average coordinates data to %s \n",pdb_file_name.c_str());

            write_frame_pdb(traj.box,this_size,this_atom_nr,this_res_nr,this_res_name,this_atom_name,this_r,traj.title,s.world_rank,&pdb_file,this_beta,this_weight,this_element,this_chain_id,0);
            fclose(pdb_file);
        } 

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Organize by largest dwell time                                                                            //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(winners.size() > 0)
        {
            int b_swap=0;
            for(b_swap=1; b_swap > 0; )
            {
                b_swap = 0;
                for(i=0; i<winners.size()-1; i++)
                {
                    if(winners[i] < winners[i+1])
                    {
                        double winners_tmp        = winners[i];
                        double winners_events_tmp = winners_events[i];
                        int    atom_a_tmp         = atom_a[i];
                        int    atom_b_tmp         = atom_b[i];
                        int    atom_nr_a_tmp      = atom_nr_a[i];
                        int    atom_nr_b_tmp      = atom_nr_b[i];
                        winners[i]                = winners[i+1];
                        winners_events[i]         = winners_events[i+1];
                        atom_a[i]                 = atom_a[i+1];
                        atom_b[i]                 = atom_b[i+1];
                        atom_nr_a[i]              = atom_nr_a[i+1];
                        atom_nr_b[i]              = atom_nr_b[i+1];
                        winners[i+1]              = winners_tmp;
                        winners_events[i+1]       = winners_events_tmp;
                        atom_a[i+1]               = atom_a_tmp;
                        atom_b[i+1]               = atom_b_tmp;
                        atom_nr_a[i+1]            = atom_nr_a_tmp;
                        atom_nr_b[i+1]            = atom_nr_b_tmp;

                        b_swap = 1;
                    }
                }
            }
        }
        else
        {
            printf("  No contacts found for the selected residue \n");
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Write PyMol select commands for identifying relevant contacts based on mean dwell time                    //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        FILE *contacts_file;
        string contacts_file_name = chop_and_add_tag(p.ck_file_name,"_contacts_dwell_t.pml");
        contacts_file             = fopen(contacts_file_name.c_str(), "w");
        if(contacts_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",contacts_file_name.c_str());
        }
        else
        {
            for(i=0; i<winners.size(); i++)
            {
                int    atom_number_1 = atom_a[i];
                int    atom_number_2 = atom_b[i];
                double freq          = (winners[i]*winners_events[i])/((double)p.norm_factor*p.delta_t); 
                double dash_rad      = (winners[i]/winners[0])*(p.max_dash_rad - p.min_dash_rad) + p.min_dash_rad;
                string dist_name     = "dist_" + to_string(i);

                if(dash_rad >= 0.0)
                {
                    fprintf(contacts_file,"\nContact %d: Residence time %f: Binding events %f: Frequency %f: Atom_nr b %d: Atom_nr a %d \n",i,winners[i],winners_events[i],freq,atom_nr_a[i],atom_nr_b[i]);
                    fprintf(contacts_file,"select atom_a, (resi %d & resn %s & name %s) \n",this_res_nr[atom_number_1-1]%10000,this_res_name[atom_number_1-1].c_str(),this_atom_name[atom_number_1-1].c_str());
                    fprintf(contacts_file,"select atom_b, (resi %d & resn %s & name %s) \n",this_res_nr[atom_number_2-1]%10000,this_res_name[atom_number_2-1].c_str(),this_atom_name[atom_number_2-1].c_str());
                    fprintf(contacts_file,"select %d, (resi %d & resn %s & name %s) or ((resi %d & resn %s & name %s)) \n",i,this_res_nr[atom_number_1-1]%10000,this_res_name[atom_number_1-1].c_str(),this_atom_name[atom_number_1-1].c_str(),this_res_nr[atom_number_2-1]%10000,this_res_name[atom_number_2-1].c_str(),this_atom_name[atom_number_2-1].c_str());
                    fprintf(contacts_file,"dist %s, (atom_a), (atom_b) \n",dist_name.c_str());
                    fprintf(contacts_file,"set dash_radius, %f, %s \n",dash_rad,dist_name.c_str());
                }
            }   
            fclose(contacts_file);
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Organize by contact frequency                                                                             //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(winners.size() > 0)
        {
            int b_swap=0;
            for(b_swap=1; b_swap > 0; )
            {
                b_swap = 0;
                for(i=0; i<winners.size()-1; i++)
                {
                    if((winners[i]*winners_events[i]) < (winners[i+1]*winners_events[i+1]))
                    {
                        double winners_tmp        = winners[i];
                        double winners_events_tmp = winners_events[i];
                        int    atom_a_tmp         = atom_a[i];
                        int    atom_b_tmp         = atom_b[i];
                        int    atom_nr_a_tmp      = atom_nr_a[i];
                        int    atom_nr_b_tmp      = atom_nr_b[i];
                        winners[i]                = winners[i+1];
                        winners_events[i]         = winners_events[i+1];
                        atom_a[i]                 = atom_a[i+1];
                        atom_b[i]                 = atom_b[i+1];
                        atom_nr_a[i]              = atom_nr_a[i+1];
                        atom_nr_b[i]              = atom_nr_b[i+1];
                        winners[i+1]              = winners_tmp;
                        winners_events[i+1]       = winners_events_tmp;
                        atom_a[i+1]               = atom_a_tmp;
                        atom_b[i+1]               = atom_b_tmp;
                        atom_nr_a[i+1]            = atom_nr_a_tmp;
                        atom_nr_b[i+1]            = atom_nr_b_tmp;

                        b_swap = 1;
                    }
                }
            }
        }
        else
        {
            printf("  No contacts found for the selected residue \n");
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Write PyMol select commands for identifying relevant contacts based on frequency                          //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        contacts_file_name = chop_and_add_tag(p.ck_file_name,"_contacts_freq.pml");
        contacts_file             = fopen(contacts_file_name.c_str(), "w");
        if(contacts_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",contacts_file_name.c_str());
        }
        else
        {
            for(i=0; i<winners.size(); i++)
            {
                int    atom_number_1 = atom_a[i];
                int    atom_number_2 = atom_b[i];
                double freq          = (winners[i]*winners_events[i])/((double)p.norm_factor*p.delta_t);
                double dash_rad      = freq*(p.max_dash_rad - p.min_dash_rad) + p.min_dash_rad;
                string dist_name     = "dist_" + to_string(i); 

                if(dash_rad >= 0.0)
                {
                    fprintf(contacts_file,"\nContact %d: Residence time %f: Binding events %f: Frequency %f: Atom_nr b %d: Atom_nr a %d \n",i,winners[i],winners_events[i],freq,atom_nr_a[i],atom_nr_b[i]);
                    fprintf(contacts_file,"select atom_a, (resi %d & resn %s & name %s) \n",this_res_nr[atom_number_1-1]%10000,this_res_name[atom_number_1-1].c_str(),this_atom_name[atom_number_1-1].c_str());
                    fprintf(contacts_file,"select atom_b, (resi %d & resn %s & name %s) \n",this_res_nr[atom_number_2-1]%10000,this_res_name[atom_number_2-1].c_str(),this_atom_name[atom_number_2-1].c_str());
                    fprintf(contacts_file,"select %d, (resi %d & resn %s & name %s) or ((resi %d & resn %s & name %s)) \n",i,this_res_nr[atom_number_1-1]%10000,this_res_name[atom_number_1-1].c_str(),this_atom_name[atom_number_1-1].c_str(),this_res_nr[atom_number_2-1]%10000,this_res_name[atom_number_2-1].c_str(),this_atom_name[atom_number_2-1].c_str());
                    fprintf(contacts_file,"dist %s, (atom_a), (atom_b) \n",dist_name.c_str());
                    fprintf(contacts_file,"set dash_radius, %f, %s \n",dash_rad,dist_name.c_str());
                }
            }
            fclose(contacts_file);
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
    s.program_name = "Contact Kinetics";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                                         s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                                         s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                                        s.world_rank, s.cl_tags, &p.b_print,     0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                                          s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                                                     s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                                      s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                                        s.world_rank, s.cl_tags, &p.b_lsq,       0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                                          s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",                          s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_s(argc,argv,"-ck",     p.ck_file_name,               "Output file with contacts kinetics data (dat)",                                      s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_s(argc,argv,"-pf_pdb", p.pf_pdb_file_name,           "PDB file with selected protein (pdb) ",                                              s.world_rank, s.cl_tags, &p.b_pf_pdb,    0);
    add_argument_mpi_s(argc,argv,"-pf_prm", p.protein_finder_param_name,  "File with additional protein finder parameters (prm)",                               s.world_rank, s.cl_tags, &p.b_pf_param,  0);
    add_argument_mpi_s(argc,argv,"-report", p.report_file_name,           "Selection card with contact pair. Report frames when this contact is present (crd)", s.world_rank, s.cl_tags, &p.b_report,    0);
    add_argument_mpi_d(argc,argv,"-cdist",  &p.contact_cutoff,            "The contact cutoff distance (nm)",                                                   s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_i(argc,argv,"-test",   &p.test,                      "Print PyMol select commands for identified contacts? (0:no, 1:yes)",                 s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_i(argc,argv,"-range",  &p.range,                     "Noise filter half-width ",                                                           s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_i(argc,argv,"-dump",   &p.dump,                      "Dump contacts on the last trajectory frame? (0:no, 1:yes) ",                         s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_d(argc,argv,"-dt",     &p.delta_t,                   "Effective time step between trajectory frames analyzed (ps) ",                       s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_d(argc,argv,"-bin",    &p.bin_width,                 "Residence time histogram bin width (ps) ",                                           s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_i(argc,argv,"-resi",   &p.resi,                      "Target residue for measuring contacts ",                                             s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_d(argc,argv,"-max",    &p.max_dash_rad,              "Maximum thickness of dash for PyMOL distance commands ",                             s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_d(argc,argv,"-min",    &p.min_dash_rad,              "Minimum thickness of dash for PyMOL distance commands ",                             s.world_rank, s.cl_tags, nullptr,        0);
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
    check_extension_mpi(s.world_rank,"-ck",p.ck_file_name,".dat");

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
    if(p.b_report == 1)
    {
        check_extension_mpi(s.world_rank,"-report",p.report_file_name,".crd");
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

    //read in atom selection parameter file
    Index report_pairs;
    if(p.b_report == 1)
    {
        report_pairs.get_index(p.report_file_name);
    }

    //run proten finder
    traj.get_protein(p.protein_finder_param_name,p.b_pf_param);

    //print a pdb with distinguished protein
    traj.write_protein(p.pf_pdb_file_name,p.b_pf_pdb);

    //print info about the protein
    traj.get_prot_stats();

    //set the parallelization scheme to the protein atoms for analyzing the filter
    traj.parallelize_by_prot(traj.prot.size());

    //read in binding events
    Binding_events target_events;
    if(p.b_be == 1)
    {
        int result = target_events.get_info(p.be_file_name);
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

    //store the contact matrix dimensions
    p.size_x = traj.prot.size();       //number of atoms from the protein
    p.size_y = p.resi_size;            //number of atoms from the resi selection

    //Estimate the memory requirements
    double mem_2 = (double)p.resi_size*3.0*8.0;                                            //x,y,z_coord
    double mem_3 = (double)traj.prot.size()*3.0*8.0;                                       //prot_x,y,z
    double mem_4 = (double)p.size_y*(double)traj.prot.size()*(double)(2*p.range + 1)*4.0;  //filter 
    double mem_5 = (double)p.size_y*(double)traj.my_prots*4.0;                             //bound
    double mem_6 = (double)p.size_y*(double)traj.my_prots*8.0*2.0;                         //dwell_t and events

    double mem_est = (mem_2 + mem_3 + mem_4 + mem_5 + mem_6)/1000000.0;

    if(s.world_rank == 0)
    {
        printf("Estimated memory required: %f MB \n\n",mem_est);
    }

    //allocate memory for target residue coords
    dv1d x_coord(p.resi_size,0.0);                                         //holds the x coordinates 
    dv1d y_coord(p.resi_size,0.0);                                         //holds the x coordinates 
    dv1d z_coord(p.resi_size,0.0);                                         //holds the x coordinates 

    //allocate memory for average prot coords
    dv1d prot_x(traj.prot.size(),0.0);                                     //holds the protein x_coordinates
    dv1d prot_y(traj.prot.size(),0.0);                                     //holds the protein x_coordinates
    dv1d prot_z(traj.prot.size(),0.0);                                     //holds the protein x_coordinates

    //allocate memory for analyzing the residence times
    iv3d filter (2*p.range + 1, iv2d(p.size_y, iv1d(traj.prot.size(),0))); //noise filter for the contact matrices
    iv2d bound  (p.size_y, iv1d(traj.my_prots, 0));                        //keeps track of how many frames a contacts has been around for 
    dv2d dwell_t(traj.my_prots, dv1d(p.size_y, 0.0));                      //keeps a sum of the residence time for each contact
    dv2d events (traj.my_prots, dv1d(p.size_y, 0.0));                      //keeps track of how many binding events there were for each contact
    iv3d filter_history (2*p.range + 1, iv2d(0, iv1d(2,0)));               //store a record of contacts in the noise filter

    //create object for working with contact matrices
    Contacts cont; 

    //initialize contacts
    cont.init(p.ck_file_name,p.size_y,traj.prot.size());

    //open files for writing temporary contact profiles data
    cont.prime_tmp_file();

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

        binding_contacts(traj,s,p,report_pairs,cont,bound_lipids);

        get_mean_coords(traj,s,p,x_coord,y_coord,z_coord,prot_x,prot_y,prot_z,bound_lipids);

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //close tmp contacts files
    cont.close_tmp_file();

    MPI_Barrier(MPI_COMM_WORLD);

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //analyze contact matrices
    perf.log_time(finalize_analysis(traj,s,p,x_coord,y_coord,z_coord,prot_x,prot_y,prot_z,filter,bound,dwell_t,events,cont,filter_history,bound_lipids),"Fin Ana");

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

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
#include "MosAT/program_variables//pv_dihedrals.h"           //This has the variables specific to the analysis program
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
#include "headers/histo.h"                                   //This has routines used for making a histogram

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function analyzes target residues and removes those specified by the user                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_new_target_residues(Trajectory &traj,system_variables &s,program_variables &p,string &target_res_file_name)
{
    string list_file_name = add_tag(target_res_file_name,"_ex");

    if(s.world_rank == 0) //rank 0 updates the target res file
    {
        int i = 0;    //standard variable used in loops
        int j = 0;    //standard variable used in loops

        sv1d new_list(0); //hold the new list of target residues

        Index initial_targets;
        Index ex;

        //read the index files
        initial_targets.get_index(target_res_file_name);
        ex.get_index(p.exclude_file_name);

        for(i=0; i<initial_targets.index_i.size(); i++) //loop over the target residues
        {
            int b_exclude       = 0;
            string residue_name = traj.res_name[traj.res_start[initial_targets.index_i[i]-1]];

            //check if residue is on the exclude list
            for(j=0; j<ex.index_s.size(); j++) //loop over exclude list
            {
                if(strcmp(residue_name.c_str(), ex.index_s[j].c_str()) == 0) //residue is on the exclude list
                {
                    b_exclude = 1;
                }
            }

            //generate new list of target residues
            if(b_exclude == 1)
            {
                string entry = "#";
                entry        = entry + initial_targets.index_s[i] + " " + "#" + residue_name;
                new_list.push_back(entry);
            }
            else
            {
                new_list.push_back(initial_targets.index_s[i]);
            }
        }

        //write new list to file
        FILE *list_file = fopen(list_file_name.c_str(), "w");
        if(list_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",list_file_name.c_str());
        }
        else
        {
            for(i=0; i<new_list.size(); i++) //loop over list
            {
                fprintf(list_file,"%s \n",new_list[i].c_str());
            }
            fclose(list_file);
        }
    }

    //replace old index with the new list
    target_res_file_name = list_file_name;
        
    MPI_Barrier(MPI_COMM_WORLD); //need a barrier to make sure all cores wait until the list is updated.
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function projects angle data onto the protein structure                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void project_to_b_factor(Trajectory &traj,system_variables &s,program_variables &p,Index &target_res,string filename,string tag,dv1d &data)
{
    int i = 0;    //standard variable used in loops
    int j = 0;    //standard variable used in loops

    //initialize b factor
    for(i=0; i<traj.atoms(); i++)
    {
        traj.beta[i] = -1.0;
    }

    //add data to the b factor
    for(i=0; i<target_res.index_s.size(); i++) //loop over target dihedrals
    {
        int min = traj.res_start[target_res.index_i[i]-1];
        int max = traj.res_end[target_res.index_i[i]-1];

        for(j=min; j<=max; j++) //loop over current residue atoms
        {
            traj.beta[j] = data[i];
        }
    }

    //write data to output file
    string pdb_file_name = chop_and_add_tag(filename,tag.c_str());
    FILE *pdb_file = fopen(pdb_file_name.c_str(), "w");
    if(pdb_file == NULL)
    {
        printf("failure opening %s (pdb single). Make sure the file exists. \n",pdb_file_name.c_str());
    }
    write_frame_pdb(traj.ibox,traj.atoms(),traj.atom_nr,traj.res_nr,traj.res_name,traj.atom_name,traj.r_ref,traj.title,s.world_rank,&pdb_file,traj.beta,traj.weight,traj.element,traj.chain_id,1);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reports data to a text file                                                                 //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void report_data(Index &target_res,string filename,string mode,string title,dv1d &data)
{
    int i = 0;     //standard variable used in loops

    FILE *stats_file = fopen(filename.c_str(),mode.c_str());

    fprintf(stats_file," %20s ",title.c_str());
    for(i=0; i<target_res.index_s.size(); i++) //loop over angles
    {
        fprintf(stats_file," %10f ",data[i]);
    }
    fprintf(stats_file,"\n");

    fclose(stats_file);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function takes 2 angles in degrees and returns the difference                                        //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double get_angle(double angle1,double angle2)
{
    double pi         = 3.14159;   //this times a circles radius gives it circumference
    double deg_to_rad = pi/180.0;  //used for convenient conversion between units
    double rad_to_deg = 180.0/pi;  //used for convenient conversion between units

    //convert angles to radians
    angle1 = deg_to_rad*angle1;
    angle2 = deg_to_rad*angle2;

    double x1      = cos(angle1);
    double y1      = sin(angle1);
    double x2      = cos(angle2);
    double y2      = sin(angle2);
    double u_dot_v = x1*x2 + y1*y2;
    double mag_v   = sqrt(x1*x1 + y1*y1);
    double mag_u   = sqrt(x2*x2 + y2*y2);
    double delta   = rad_to_deg*acos(u_dot_v/(mag_u*mag_v));

    if(angle1 == angle2)  
    {
        delta = 0.0;
    }

    return delta;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function takes 2 angles in degrees and returns the difference                                        //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double get_angle_alt(double angle1,double angle2)
{
    double pi         = 3.14159;   //this times a circles radius gives it circumference
    double deg_to_rad = pi/180.0;  //used for convenient conversion between units
    double rad_to_deg = 180.0/pi;  //used for convenient conversion between units

    //convert angles to radians
    angle1 = deg_to_rad*angle1;
    angle2 = deg_to_rad*angle2;

    double delta = rad_to_deg*atan2(sin(angle1-angle2),cos(angle1-angle2));

    if(angle1 == angle2)
    {
        delta = 0.0;
    }

    return delta;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the time average angle                                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double time_avg_angle(Trajectory &traj,dv1d &angle)
{
    int    i   = 0;
    double sum = 0.0;

    for(i=0; i<traj.get_ef_frames(); i++) //loop over trajectory frames
    {
        sum = sum + angle[i];
    }
    return sum/(double)traj.get_ef_frames();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the angle stdev                                                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double stdev_angle(Trajectory &traj,double mean,dv1d &angle)
{
    int    i             = 0;
    double sum_o_squares = 0.0;

    for(i=0; i<traj.get_ef_frames(); i++) //loop over trajectory frames
    {
        sum_o_squares = sum_o_squares + pow(angle[i]-mean,2);
    }
    return sqrt(sum_o_squares/(double)(traj.get_ef_frames()-1));
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the time average angle using projections on the unit circle                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double time_avg_angle_vec(Trajectory &traj,dv1d &dih)
{
    int    i          = 0;         //standard variable used in loops
    double pi         = 3.14159;   //this times a circles radius gives it circumference
    double deg_to_rad = pi/180.0;  //used for convenient conversion between units
    double rad_to_deg = 180.0/pi;  //used for convenient conversion between units
    double sum_x      = 0.0;
    double sum_y      = 0.0;
    double angle      = 0.0;

    for(i=0; i<traj.get_ef_frames(); i++) //loop over trajectory frames
    {
        angle = deg_to_rad*dih[i];
        sum_x = sum_x + cos(angle);
        sum_y = sum_y + sin(angle);
    }
    double avg_x = sum_x/(double)traj.get_ef_frames();
    double avg_y = sum_y/(double)traj.get_ef_frames();

    return rad_to_deg*atan2(avg_y,avg_x);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the angle stdev using projections on the unit circle                                //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double stdev_angle_vec(Trajectory &traj,double mean,dv1d &dih)
{
    int i = 0;   //standard variable used in loops

    double sum_dist = 0.0;
    for(i=0; i<traj.get_ef_frames(); i++) //loop over trajectory frames
    {
        double delta = get_angle(dih[i],mean);
        sum_dist     = sum_dist + delta*delta;
    }
    return sqrt(sum_dist/double(traj.get_ef_frames()-1));
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the dihedral angles                                                                //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void dihedrals(Trajectory &traj,system_variables &s,program_variables &p,Index &target_res,Index &def,dv2d &dih)
{
    //psi:  N:CA:CO:N      ----> N  : CA : ANGLE CO : N
    //phi:  CO:N:CA:CO     ----> CO : N  : ANGLE CA : CO
    //chi1: N:CA:CB:CG     ----> N  : CA : ANGLE CB : CG
    
    int     i = 0;       //standard variable used in loops
    int     j = 0;       //standard variable used in loops
    int     k = 0;       //standard variable used in loops
    int     l = 0;       //standard variable used in loops
    double pi = 3.14159; //this times a circles radius gives it circumference

    for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
    {
        //get the first and last atoms of the current res
        int min = traj.p_res_start(i);  
        int max = traj.p_res_end(i);    

        //store the current index in prot. used to find the previous res info
        int i_store = i;               

        //jump to the next residue
        i = traj.next_prot_res(i);

        for(j=0; j<target_res.index_s.size(); j++) //loop over the target residues
        {
            if(traj.res_nr[traj.prot[i]-1] == target_res.index_i[j]) //residue is a target res
            {
                int found_def = 0;
                for(k=0; k<def.index_s.size(); k+=5) //loop over dihedral definitions
                {
                    if(strcmp(traj.res_name[min].c_str(), def.index_s[k].c_str()) == 0) //residue type is correct
                    {
                        int atom_i = -1;
                        int atom_j = -1;
                        int atom_k = -1;
                        int atom_l = -1;
                        string ai  = def.index_s[k+1];
                        string aj  = def.index_s[k+2];
                        string ak  = def.index_s[k+3];
                        string al  = def.index_s[k+4];

                        found_def = 1;

                        if(p.type == 0) //phi
                        {
                            //atoms span res i and i-1
                            //ai = C,i-1
                            //aj = N,i
                            //ak = CA,i
                            //al = C,i

                            //get the first and last atom of res i-1
                            int min_2 = traj.p_res_start(i_store-1);
                            int max_2 = traj.p_res_end(i_store-1);

                            //find atom i
                            for(l=min_2; l<=max_2; l++) //loop over previous residue atoms
                            {
                                if(strcmp(traj.atom_name[l].c_str(), ai.c_str()) == 0) //atom i
                                {
                                    atom_i = l;
                                }
                            }

                            //find atoms j,k and l
                            for(l=min; l<=max; l++) //loop over current residue atoms
                            {
                                if(strcmp(traj.atom_name[l].c_str(), aj.c_str()) == 0) //atom j
                                {
                                    atom_j = l;
                                }
                                if(strcmp(traj.atom_name[l].c_str(), ak.c_str()) == 0) //atom k
                                {
                                    atom_k = l;
                                }
                                if(strcmp(traj.atom_name[l].c_str(), al.c_str()) == 0) //atom l
                                {
                                    atom_l = l;
                                }
                            }
                            //printf("dih %5d i %9d j %9d k %9d l %9d \n",target_res.index_i[j],atom_i,atom_j,atom_k,atom_l);
                            if(p.b_noisy == 1)
                            {
                                if(s.world_rank == 0)
                                { 
                                    printf("dih %5d ai %9s(%d) aj %9s(%d) ak %9s(%d) al %9s(%d) \n",target_res.index_i[j],traj.atom_name[atom_i].c_str(),traj.res_nr[atom_i],traj.atom_name[atom_j].c_str(),traj.res_nr[atom_j],traj.atom_name[atom_k].c_str(),traj.res_nr[atom_k],traj.atom_name[atom_l].c_str(),traj.res_nr[atom_l]);
                                }
                            }
                        }
                        else if(p.type == 1) //psi
                        {
                            //atoms span res i and i + 1
                            //ai = N,i
                            //aj = CA,i
                            //ak = C,i
                            //al = N,i+1 

                            //find the first and last atom of res i+1 
                            int min_2 = traj.p_res_start(i+1);
                            int max_2 = traj.p_res_end(i+1);

                            //find atoms i,j and k
                            for(l=min; l<=max; l++) //loop over current residue atoms
                            {
                                if(strcmp(traj.atom_name[l].c_str(), ai.c_str()) == 0) //atom i
                                {
                                    atom_i = l;
                                }
                                if(strcmp(traj.atom_name[l].c_str(), aj.c_str()) == 0) //atom j
                                {
                                    atom_j = l;
                                }
                                if(strcmp(traj.atom_name[l].c_str(), ak.c_str()) == 0) //atom k
                                {
                                    atom_k = l;
                                }
                            }

                            //find atom l
                            for(l=min_2; l<=max_2; l++) //loop over next residue atoms
                            {
                                if(strcmp(traj.atom_name[l].c_str(), al.c_str()) == 0) //atom l
                                {
                                    atom_l = l;
                                }
                            }
                            //printf("dih %5d i %9d j %9d k %9d l %9d max_1 %d prot_size %d \n",target_res.index_i[j],atom_i,atom_j,atom_k,atom_l,max,traj.prot.size());
                            if(p.b_noisy == 1)
                            {
                                if(s.world_rank == 0)
                                {
                                    printf("dih %5d ai %9s(%d) aj %9s(%d) ak %9s(%d) al %9s(%d) \n",target_res.index_i[j],traj.atom_name[atom_i].c_str(),traj.res_nr[atom_i],traj.atom_name[atom_j].c_str(),traj.res_nr[atom_j],traj.atom_name[atom_k].c_str(),traj.res_nr[atom_k],traj.atom_name[atom_l].c_str(),traj.res_nr[atom_l]);
                                }
                            }

                            if(max == traj.prot[traj.prot.size()-1]-1)
                            {
                                if(s.world_rank == 0)
                                {
                                    printf("Dihedral %6d runs off the protein end. i %9d j %9d k %9d l %9d \n",target_res.index_i[j],atom_i,atom_j,atom_k,atom_l);
                                }
                                MPI_Finalize();
                                exit(EXIT_SUCCESS);
                            }
                        }
                        else if(p.type == 2) //chi1
                        {
                            //atoms are all on res i
                            //ai = N,i
                            //aj = CA,i
                            //ak = CB,i
                            //al = CG,i

                            //find atoms i,j,k and l
                            for(l=min; l<=max; l++) //loop over current residue atoms
                            {
                                if(strcmp(traj.atom_name[l].c_str(), ai.c_str()) == 0) //atom i
                                {
                                    atom_i = l;
                                }
                                if(strcmp(traj.atom_name[l].c_str(), aj.c_str()) == 0) //atom j
                                {
                                    atom_j = l;
                                }
                                if(strcmp(traj.atom_name[l].c_str(), ak.c_str()) == 0) //atom k
                                {
                                    atom_k = l;
                                }
                                if(strcmp(traj.atom_name[l].c_str(), al.c_str()) == 0) //atom l
                                {
                                    atom_l = l;
                                }
                            }
                            if(p.b_noisy == 1)
                            {
                                if(s.world_rank == 0)
                                {
                                    printf("dih %5d ai %9s(%d) aj %9s(%d) ak %9s(%d) al %9s(%d) \n",target_res.index_i[j],traj.atom_name[atom_i].c_str(),traj.res_nr[atom_i],traj.atom_name[atom_j].c_str(),traj.res_nr[atom_j],traj.atom_name[atom_k].c_str(),traj.res_nr[atom_k],traj.atom_name[atom_l].c_str(),traj.res_nr[atom_l]);
                                }
                            } 
                        }

                        if(atom_i == -1 || atom_j == -1 || atom_k == -1 || atom_l == -1 )
                        { 
                            if(s.world_rank == 0)
                            {
                                printf("unable to find an atom for dihedral %5d %s. atom_i %9d atom_j %9d atom_k %9d atom_l %9d \n",target_res.index_i[j],traj.res_name[min].c_str(),atom_i,atom_j,atom_k,atom_l);
                            }
                             MPI_Finalize();
                             exit(EXIT_SUCCESS);
                        }

                        //compute the dihedral angle
                        double angle = traj.dihedral_angle(atom_i,atom_j,atom_k,atom_l);

                        //convert angle to degrees
                        angle = angle*180.0/pi;

                        //store the current angle
                        dih[j][traj.current_frame] = angle;
 
                        //prints Pymol compatible atom selection info for current angle
                        if(p.b_test == 1) 
                        {
                            printf("residue %d %s atom_i %s atom_j %s atom_k %s atom_l %s angle %f \n",traj.res_nr[min],traj.res_name[min].c_str(),traj.atom_name[atom_i].c_str(),traj.atom_name[atom_j].c_str(),traj.atom_name[atom_k].c_str(),traj.atom_name[atom_l].c_str(),angle);
                            printf("select atom_i, resi %d & resn %s & name %s \n",traj.res_nr[atom_i],traj.res_name[atom_i].c_str(),traj.atom_name[atom_i].c_str());
                            printf("select atom_j, resi %d & resn %s & name %s \n",traj.res_nr[atom_j],traj.res_name[atom_j].c_str(),traj.atom_name[atom_j].c_str());
                            printf("select atom_k, resi %d & resn %s & name %s \n",traj.res_nr[atom_k],traj.res_name[atom_k].c_str(),traj.atom_name[atom_k].c_str());
                            printf("select atom_l, resi %d & resn %s & name %s \n",traj.res_nr[atom_l],traj.res_name[atom_l].c_str(),traj.atom_name[atom_l].c_str());
                            printf("get_dihedral atom_i, atom_j, atom_k, atom_l \n");
                            printf("select quad, atom_i + atom_j + atom_k + atom_l \n");
                            printf("show licorice, quad \n");
                            printf("orient quad \n");
                        }
                    }
                }

                //dihedral angle was not defined 
                if(found_def == 0)
                { 
                    if(s.world_rank == 0)
                    {
                        printf("Unable to find a definition for residue %d %s. Please check your definitions list %s. \n",target_res.index_i[j],traj.res_name[min].c_str(),p.dihedrals_file_name.c_str());
                    }
                    MPI_Finalize();
                    exit(EXIT_SUCCESS);
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects the angles and writes them to file.                                                //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,Trajectory &traj_b,system_variables &s,program_variables &p,Index &target_res,
                         Index &target_res_cmp,dv2d &dih,dv2d &dih_b,dv2d &dih_cmp)
{
    int i             = 0;         //standard variable used in loops  
    int j             = 0;         //standard variable used in loops 
    double pi         = 3.14159;   //this times a circles radius gives it circumference
    double rad_to_deg = 180.0/pi;  //used for convenient conversion between units
    double deg_to_rad = pi/180.0;  //used for convenient conversion between units

    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank ==0)
    {
        printf("\nCollecting dihedrals data. This could take some time depending on the number of angles and frames analyzed. \n");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Collect dihedral angles from each mpi process                                                             //
    //                                                                                                           //
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    for(i=0; i<target_res.index_s.size(); i++) //loop over angles
    {
        collect_dv1d(s.world_size,s.world_rank,dih[i]);
      
        if(p.b_cmp == 1)
        {
            collect_dv1d(s.world_size,s.world_rank,dih_cmp[i]);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Copy dih_b over to a 1d vector to be used in fuctions for reporting data                                  //
    //                                                                                                           //
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    dv1d dih_b_prime(target_res.index_s.size(),0.0);
    for(i=0; i<target_res.index_s.size(); i++)
    {
        dih_b_prime[i] = dih_b[i][0];
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Write data to output file                                                                                 //
    //                                                                                                           //
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    if(s.world_rank == 0)
    {    
        //write angles vs time
        FILE *dih_file = fopen(p.dih_file_name.c_str(),"w");
        
        if(dih_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",p.dih_file_name.c_str());
        }
        else
        {
            printf("Writing dehedrals data to %s \n",p.dih_file_name.c_str());

            for(i=0; i<traj.get_ef_frames(); i++) //loop over trajectory frames
            {
                fprintf(dih_file," %10d ",i);
                for(j=0; j<target_res.index_s.size(); j++) //loop over angles
                {
                    fprintf(dih_file," %10.3f ",dih[j][i]);
                }
                fprintf(dih_file,"\n");
            }
            fclose(dih_file);
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Compute average and standard deviation of angles                                                          //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////    
        dv1d avg(target_res.index_s.size(),0.0);               //average angle
        dv1d stdev(target_res.index_s.size(),0.0);             //stdev angle
        dv1d avg_vec(target_res.index_s.size(),0.0);           //average angle from vec
        dv1d stdev_vec(target_res.index_s.size(),0.0);         //stdev angle from vec
        dv1d avg_dist_b(target_res.index_s.size(),0.0);        //average distance from target
        dv1d stdev_dist_b(target_res.index_s.size(),0.0);      //stdev of distance from target
        dv1d avg_vec_cmp(target_res.index_s.size(),0.0);       //average angle in -cmp
        dv1d stdev_vec_cmp(target_res.index_s.size(),0.0);     //stdev for angle in -cmp
        dv1d avg_dif_cmp(target_res.index_s.size(),0.0);       //average over instantaneous distances between vectors in -dih and -cmp
        dv1d dist_cmp(target_res.index_s.size(),0.0);          //distance between average vectors of -dih and -cmp
        dv1d avg_dist_b_alt(target_res.index_s.size(),0.0);    //average distance from target

        printf("Analyzing angles and computing averages/variance etc. \n");

        for(j=0; j<target_res.index_s.size(); j++) //loop over angles
        {
            avg[j]       = time_avg_angle(traj,dih[j]);             //compute average angle
            stdev[j]     = stdev_angle(traj,avg[j],dih[j]);         //compute stdev angle
            avg_vec[j]   = time_avg_angle_vec(traj,dih[j]);         //compute average angle from vec
            stdev_vec[j] = stdev_angle_vec(traj,avg_vec[j],dih[j]); //compute stdev angle from vec

            if(p.b_traj_b == 1)
            {
                //compute average deviation from target by measuring the angle between the vectors and the target vector
                double sum_x = 0.0;
                double sum_y = 0.0;
                for(i=0; i<traj.get_ef_frames(); i++) //loop over trajectory frames
                {
                    double delta = get_angle(dih[j][i],dih_b_prime[j]);

                    sum_x = sum_x + cos(deg_to_rad*delta);
                    sum_y = sum_y + sin(deg_to_rad*delta);
                }
                avg_dist_b[j] = rad_to_deg*atan2(sum_y,sum_x);

                //compute the standard deviation of the deviation from the target
                double sum_o_squares = 0.0;
                for(i=0; i<traj.get_ef_frames(); i++) //loop over trajectory frames
                {
                    double delta  = get_angle(dih[j][i],dih_b_prime[j]);
                    sum_o_squares = sum_o_squares + pow(delta-avg_dist_b[j],2);
                }
                stdev_dist_b[j] = sqrt(sum_o_squares/double(traj.get_ef_frames()-1));
            }

            if(p.b_cmp == 1)
            {
                avg_vec_cmp[j]   = time_avg_angle_vec(traj,dih_cmp[j]);             //compute average angle from vec
                stdev_vec_cmp[j] = stdev_angle_vec(traj,avg_vec_cmp[j],dih_cmp[j]); //compute stdev angle from vec 

                //compute the average difference between the angles in -res and -rescmp
                double sum_dist_cmp = 0.0;
                for(i=0; i<traj.get_ef_frames(); i++) //loop over trajectory frames
                {
                    sum_dist_cmp = sum_dist_cmp + get_angle(dih[j][i],dih_cmp[j][i]);
                }
                avg_dif_cmp[j] = sum_dist_cmp/double(traj.get_ef_frames());

                //compute the difference between the average angles in -res and -rescmp
                dist_cmp[j]  = get_angle(avg_vec[j],avg_vec_cmp[j]);
            }
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Write data to output files                                                                                //
        //                                                                                                           //
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
        printf("Writing angle statistics to files. \n");

        string stats_file_name = add_tag(p.dih_file_name,"_stats");

        report_data(target_res,stats_file_name,"w","res_id",target_res.index_d);
        report_data(target_res,stats_file_name,"a","avg",avg);
        report_data(target_res,stats_file_name,"a","stdev",stdev);
        report_data(target_res,stats_file_name,"a","avg_vec",avg_vec);
        report_data(target_res,stats_file_name,"a","stdev_vec",stdev_vec);
        if(p.b_traj_b == 1)
        {
            report_data(target_res,stats_file_name,"a","target",dih_b_prime);
            report_data(target_res,stats_file_name,"a","dist_target",avg_dist_b);
            report_data(target_res,stats_file_name,"a","stdev_target",stdev_dist_b);
        }
        if(p.b_cmp == 1)
        {
            report_data(target_res,stats_file_name,"a","avg_vec_cmp",avg_vec_cmp);
            report_data(target_res,stats_file_name,"a","stdev_vec_cmp",stdev_vec_cmp);
            report_data(target_res,stats_file_name,"a","avg_dif_cmp",avg_dif_cmp);
            report_data(target_res,stats_file_name,"a","dist_avg_vec_cmp",dist_cmp);
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Write data to the B factor                                                                                //
        //                                                                                                           //
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
        printf("Projecting data onto the protein coordinates. \n");

        project_to_b_factor(traj,s,p,target_res,p.dih_file_name,"_avg_vec.pdb",avg_vec);
        project_to_b_factor(traj,s,p,target_res,p.dih_file_name,"_stdev_vec.pdb",stdev_vec);
        if(p.b_traj_b == 1)
        {
            project_to_b_factor(traj,s,p,target_res,p.dih_file_name,"_target.pdb",dih_b_prime);
            project_to_b_factor(traj,s,p,target_res,p.dih_file_name,"_dist_target.pdb",avg_dist_b);
            project_to_b_factor(traj,s,p,target_res,p.dih_file_name,"_stdev_target.pdb",stdev_dist_b);
        }
        if(p.b_cmp == 1)
        {
            project_to_b_factor(traj,s,p,target_res,p.dih_file_name,"_avg_vec_cmp.pdb",avg_vec_cmp);
            project_to_b_factor(traj,s,p,target_res,p.dih_file_name,"_stdev_vec_cmp.pdb",stdev_vec_cmp);
            project_to_b_factor(traj,s,p,target_res,p.dih_file_name,"_avg_dif_cmp.pdb",avg_dif_cmp);
            project_to_b_factor(traj,s,p,target_res,p.dih_file_name,"_dist_avg_vec_cmp.pdb",dist_cmp);
            project_to_b_factor(traj,s,p,target_res_cmp,p.dih_file_name,"_dist_avg_vec_cmp_b.pdb",dist_cmp);
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Analyze deviation from target angles and get a probability distribution for each residue and each delta   //
        //                                                                                                           //
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
        if(p.b_traj_b == 1)
        {
            printf("Analyzing the frequency of deviation from the target angles. \n");

            dv2d freq(180,dv1d(target_res.index_s.size(),0.0)); //holds the frequency for each deviation and each resid

            for(i=0; i<target_res.index_s.size(); i++) //loop over target angles
            {
                dv1d delta(0,0.0); //holds deviation from target for each frame

                for(j=0; j<traj.get_ef_frames(); j++) //loop over frames analyzed
                {
                    delta.push_back(get_angle_alt(dih[i][j],dih_b[i][0]));
                }

                //bin the deviations data
                Histogram_d histo;
                histo.set_range(0.0,180.0);
                histo.bin_data(delta,1.0);

                for(j=0; j<histo.get_num_bins(); j++) //loop over the bins
                {
                    double probability = histo.bins[j]/(double)histo.get_size();
                    freq[j][i]         = probability; 
                }
            }

            //write frequencies to output file
            string freq_file_name = add_tag(p.dih_file_name,"_freq");
            FILE *freq_file       = fopen(freq_file_name.c_str(), "w");
            if(freq_file == NULL)
            {
                printf("failure opening %s. Make sure the file exists. \n",freq_file_name.c_str());
            }        
            else
            {
                for(i=0; i<180; i++) //loop over angles
                {
                    for(j=0; j<target_res.index_s.size(); j++) //loop over target residues
                    {
                        fprintf(freq_file," %6.4f ",freq[i][j]);
                    }
                    fprintf(freq_file,"\n");
                }
 
                fclose(freq_file);
            }
        }
    }

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
    s.program_name = "Dihedrals";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                                                     s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (tpr, gro)",                                                                     s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                                                    s.world_rank, s.cl_tags, &p.b_print,   0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                                                      s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                                                                 s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                                                  s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting",                                                                          s.world_rank, s.cl_tags, &p.b_lsq,     0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                                                      s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",                                      s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-pf_pdb", p.pf_pdb_file_name,           "PDB file with selected protein ",                                                                s.world_rank, s.cl_tags, &p.b_pf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-pf_prm", p.protein_finder_param_name,  "File with additional protein finder parameters ",                                                s.world_rank, s.cl_tags, &p.b_pf_param,0);
    add_argument_mpi_s(argc,argv,"-dih",    p.dih_file_name,              "Output data file file with dihedral angles (dat)",                                               s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-def",    p.dihedrals_file_name,        "Selection card with dihedral definitions (crd)",                                                 s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-res",    p.target_res_file_name,       "Selection card with target residue id's for measuring angles (crd) ",                            s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-type",   &p.type,                      "Which angle to measure? (0:phi, 1:psi, 2:chi1)",                                                 s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-noisy",  &p.b_noisy,                   "Print info about the selected residue/atoms? (0:no, 1:yes)",                                     s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-traj_b", p.traj_2_file_name,           "Input trajectory file with the desired angles (xtc, trr, pdb, gro)",                             s.world_rank, s.cl_tags, &p.b_traj_b,  0);
    add_argument_mpi_s(argc,argv,"-rescmp", p.target_res_cmp_file_name,   "Selection card with target residue id's that dihedral angles (-res) should be compared to (crd)",s.world_rank, s.cl_tags, &p.b_cmp,     0);
    add_argument_mpi_i(argc,argv,"-test",   &p.b_test,                    "Print info for checking angles? (0:no 1:yes)",                                                   s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-ex",     p.exclude_file_name,          "Selection card with residue names to be removed from target res list (crd)",                     s.world_rank, s.cl_tags, &p.b_ex,      0);
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

    //create a second trajectory
    Trajectory traj_b;

    //set second trajectory parameters
    if(p.b_traj_b == 1)
    {
        traj_b.set_block_parallel(on);
        traj_b.set_traj(p.traj_2_file_name);
        traj_b.set_ref(p.ref_file_name);
        traj_b.set_traj_w(p.out_file_name,p.b_print);
        traj_b.set_lsq(p.lsq_index_file_name,p.b_lsq,p.lsq_dim,p.lsq_ref);
        traj_b.set_res(1,0,1);

        perf.log_time(traj_b.build(),"Analyze Traj. B");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //check file extensions                                                                                     
    check_extension_mpi(s.world_rank,"-dih",p.dih_file_name,".dat");
    check_extension_mpi(s.world_rank,"-def",p.dihedrals_file_name,".crd");
    check_extension_mpi(s.world_rank,"-res",p.target_res_file_name,".crd");
    check_extension_mpi(s.world_rank,"-ex",p.exclude_file_name,".crd");
    if(p.b_pf_pdb == 1)
    {
        check_extension_mpi(s.world_rank,"-pf_pdb",p.pf_pdb_file_name,".pdb");
    }
    if(p.b_pf_param == 1)
    {
        check_extension_mpi(s.world_rank,"-pf_prm",p.protein_finder_param_name,".prm");
    }

    //generate new target res list that excludes specified residue types
    if(p.b_ex == 1)
    {
        get_new_target_residues(traj,s,p,p.target_res_file_name);
        if(p.b_cmp == 1)
        {
            get_new_target_residues(traj,s,p,p.target_res_cmp_file_name);
        }
    }

    //create index for residues making the dihedrals
    Index target_res;
    Index target_res_cmp;
    Index def;

    //read the index files
    target_res.get_index(p.target_res_file_name);
    if(p.b_cmp == 1)
    {
        target_res_cmp.get_index(p.target_res_cmp_file_name);
    }
    def.get_index(p.dihedrals_file_name);

    //run proten finder
    traj.get_protein(p.protein_finder_param_name,p.b_pf_param);
    traj_b.get_protein(p.protein_finder_param_name,p.b_pf_param);

    //print a pdb with distinguished protein
    traj.write_protein(p.pf_pdb_file_name,p.b_pf_pdb);

    //create structures to hold angles
    dv2d dih(target_res.index_s.size(),dv1d(traj.get_num_frames(),0.0));
    dv2d dih_b(target_res.index_s.size(),dv1d(1,0.0));
    dv2d dih_cmp(target_res_cmp.index_s.size(),dv1d(traj.get_num_frames(),0.0));

    if(p.b_cmp == 1)
    {
        if(target_res.index_s.size() != target_res_cmp.index_s.size()) //comparison is not possible
        {
            if(s.world_rank == 0)
            {
                printf("The user has opted to compare two sets of dihedral angles that are not compatible.  \n");
                printf("The size of each list must be equal. There are %d items in -res and %d in -rescmp \n",target_res.index_s.size(),target_res_cmp.index_s.size());
            }
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //print info about the worlk load distribution
    traj.workload();

    //print that analysis is beginning
    traj.report_progress();

    //measure target angles
    if(s.world_rank == 0 && p.b_traj_b == 1)
    {
        traj_b.current_frame = 0;
        traj_b.read_traj_frame();
        traj_b.do_fit();
        dihedrals(traj_b,s,p,target_res,def,dih_b);
    }

    s.t = clock();
    //read read frames of the trajectory and perform analysis
    for(traj.current_frame=0; traj.current_frame<traj.get_num_frames(); traj.current_frame++)
    {
        traj.read_traj_frame();

        traj.do_fit();

        dihedrals(traj,s,p,target_res,def,dih);

        if(p.b_cmp == 1)
        {
            dihedrals(traj,s,p,target_res_cmp,def,dih_cmp);
        }

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect dihedral angles from mpi processes and write output
    perf.log_time(finalize_analysis(traj,traj_b,s,p,target_res,target_res_cmp,dih,dih_b,dih_cmp),"Fin Ana");

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

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
#include "MosAT/program_variables/pv_h_bond_kinetics.h"      //This has the variables specific to the analysis program
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
#include "headers/histo.h"                                   //This has routines used for making a histogram
#include "headers/parallel.h"                                //This has routines used for distributing the workload over various coordinates
#include "headers/contacts.h"                                //This has routines used for reading and writing contacts matrices data files
#include "headers/grid_io.h"                                 //This has routines used for reading in grid data
#include "headers/binding_events.h"                          //This has routines used for reading binding events files

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads in the bonds list and builds a bonding list for each atom                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_bonds(Trajectory &traj,system_variables &s,program_variables &p,Index &bond,iv2d &bonds)
{
    int i = 0;  //standard variable used in loops
    int j = 0;  //standard variable used in loops

    /*                                                                                                                                                                  
     *   Notes:                                                                                                                                                           
     *   the strategy is to loop over the pairs of atoms found in the bond list (bond[]). For each pair, we loop over the bonding atoms (found in bonds[])               
     *   for the first atom in the pair. In this second loop, we check if the second atom in the pair (bond[]) is on already on the list. If not, then we                 
     *   add it. Before moving to the next pair, we loop over the bonding atoms (found in bonds[]) for the second atom in the pair. We then check if the first           
     *   atom in the pair (bond[]) is on already on the list. If not, then we add it. With this approach, each atom in the pair is given a list of bonded atoms.         
     *   We then check that each atom in the pair in on the other atoms list. This results in a 2d vector (bonds[]) where one dimension gives each atom in the system    
     *   and the second dimension lists every atom bonded to a given atom.                                                                                               
     */

    for(i=0; i<bond.index_i.size(); i+=2) //loop over bonds
    {
       int duplicate = 0;

       for(j=0; j<bonds[bond.index_i[i]-1].size(); j++) //loop over added bonds for the first atom in the pair
       {
           if(bonds[bond.index_i[i]-1][j] == bond.index_i[i+1]) //second atom in pair is already on the list
           {
               duplicate = 1;
           }
       }

       if(duplicate == 0) //atom not yet added
       {
           bonds[bond.index_i[i]-1].push_back(bond.index_i[i+1]);
       }
       else
       {
           printf("duplicate entry found in bonds list \n");
       }

       duplicate = 0;

       for(j=0; j<bonds[bond.index_i[i+1]-1].size(); j++) //loop over added bonds
       {
           if(bonds[bond.index_i[i+1]-1][j] == bond.index_i[i]) //atom already added
           {
               duplicate = 1;
           }
       }

       if(duplicate == 0) //atom not yet added
       {
           bonds[bond.index_i[i+1]-1].push_back(bond.index_i[i]);
       }
       else
       {
           printf("duplicate entry found in bonds list \n");
       }
    }
}

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
// This function checks if the acceptor/donor pair are on the exclude list                                   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int check_exclusions(Trajectory &traj,system_variables &s,program_variables &p,Index &exclude,int acceptor,int donor)
{
    int i         = 0; //standard variable used in loops
    int b_include = 1; //should the hydrogen bond be included

    for(i=0; i<exclude.index_s.size(); i+=2) //loop over excluded pairs
    {
        if(strcmp(traj.atom_name[acceptor].c_str(), exclude.index_s[i].c_str()) == 0 && strcmp(traj.atom_name[donor].c_str(), exclude.index_s[i+1].c_str()) == 0) //first acceptor second donor
        {
            b_include = 0;
        }
        else if(strcmp(traj.atom_name[acceptor].c_str(), exclude.index_s[i+1].c_str()) == 0 && strcmp(traj.atom_name[donor].c_str(), exclude.index_s[i].c_str()) == 0) //first donor second acceptor
        {
            b_include = 0;
        }
    }

    return b_include; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the number of lipid-prot h-bonds and adds it to the grid                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int check_h_bond(Trajectory &traj,program_variables &p,int acceptor,int donor,int h,Index &report)
{
    /*                                                                                                                                                                  
     *   Notes:                                                                                                                                                           
     *   Here we take in an atom index for the acceptor, donor, and hydrogen atoms. We then check if there is a hydrogen bond formed between them. The  
     *   function returns 1 if a h-bond is present. The function also prints PyMOL select commands used to visualize the bond in PyMOL when the -test 
     *   option is used. Similarly, the function prints information for selecting a hydrogen bond, like the frame in which it is detected, where the 
     *   bonds of interest are specified using the -report argument. This can be used to visualize hydrogen bonds that are unusual.                 
     */

    int result = 0;
    int pi     = 3.1415926535;
    int i      = 0;

    rvec m;
    rvec n;

    for(i=0; i<3; i++) //loop over 3 dimensions
    {
        m[i] = traj.r[acceptor][i] - traj.r[donor][i];
        n[i] = traj.r[h][i]        - traj.r[donor][i];
    }

    double angle = gmx_angle(m,n);
    double dist  = sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);

    //convert angle to degrees
    angle = angle*180/pi;

    if(dist < 0.35 && angle < 30)
    {
        result = 1;

        //check h-bonds in pymol
        if(p.b_test == 1)
        {
            printf("acceptor %10d %10s %10d %10s \n",traj.atom_nr[acceptor],traj.atom_name[acceptor].c_str(),traj.res_nr[acceptor],traj.res_name[acceptor].c_str());
            printf("donor    %10d %10s %10d %10s \n",traj.atom_nr[donor],   traj.atom_name[donor].c_str(),   traj.res_nr[donor],   traj.res_name[donor].c_str());
            printf("h        %10d %10s %10d %10s \n",traj.atom_nr[h],       traj.atom_name[h].c_str(),       traj.res_nr[h],       traj.res_name[h].c_str());
            printf("angle %f dist %f \n",angle,dist);
            printf("sel a, resi %d & resn %s & name %s \n",traj.res_nr[acceptor]%10000,traj.res_name[acceptor].c_str(),traj.atom_name[acceptor].c_str());
            printf("sel d, resi %d & resn %s & name %s \n",traj.res_nr[donor]%10000,traj.res_name[donor].c_str(),traj.atom_name[donor].c_str());
            printf("sel h, resi %d & resn %s & name %s \n",traj.res_nr[h]%10000,traj.res_name[h].c_str(),traj.atom_name[h].c_str());
            printf("get_angle (a), (d), (h), state=4  \n");
            printf("dist (a), (d) \n");
            printf("show licorice, resi %d \n",traj.res_nr[acceptor]);
            printf("show licorice, resi %d \n",traj.res_nr[donor]);
            printf("sel triad, a + d + h \n");
            printf("orient triad \n");
            printf("\n");
        }

        if(p.b_report == 1)
        {
            for(i=0; i<report.index_s.size(); i+=4)  //loop over target h-bonds
            {
                if(strcmp(traj.res_name[donor].c_str(), report.index_s[i].c_str()) == 0) //donor belongs to lipid
                {
                    if(strcmp(traj.atom_name[donor].c_str(), report.index_s[i+1].c_str()) == 0) //atom type for lipid
                    {
                        if(traj.res_nr[acceptor] == report.index_i[i+2]) //resid is correct
                        {
                            if(strcmp(traj.atom_name[acceptor].c_str(), report.index_s[i+3].c_str()) == 0) //atom type for the protein
                            {
                                printf("Target hydrogen bond identified on trajectory frame: %d \n",traj.get_frame_global());
                                printf("select lip,  resi %d and resname %s \n",traj.res_nr[donor]%10000,traj.res_name[donor].c_str());
                                printf("select prot, resi %d and name %s \n",traj.res_nr[acceptor]%10000,traj.atom_name[acceptor].c_str());
                            }
                        }
                    }
                }
                else if(strcmp(traj.res_name[acceptor].c_str(), report.index_s[i].c_str()) == 0) //acceptor belongs to lipid
                {
                    if(strcmp(traj.atom_name[acceptor].c_str(), report.index_s[i+1].c_str()) == 0) //atom type for lipid
                    {
                        if(traj.res_nr[donor] == report.index_i[i+2]) //resid is correct
                        {
                            if(strcmp(traj.atom_name[donor].c_str(), report.index_s[i+3].c_str()) == 0) //atom type for the protein
                            {
                                printf("Target hydrogen bond identified on trajectory frame: %d \n",traj.get_frame_global());
                                printf("select lip,  resi %d and resname %s \n",traj.res_nr[acceptor]%10000,traj.res_name[acceptor].c_str());
                                printf("select prot, resi %d and name %s \n",traj.res_nr[donor]%10000,traj.atom_name[donor].c_str());
                            }
                        }
                    }
                }
            }
        }
    }

    return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function counts the number of atoms for each lipid type                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void count_lipid_atoms(Trajectory &traj,system_variables &s,program_variables &p,sv1d &types,iv1d &types_count)
{
    int i   = 0;                              //standard variable used in loops
    int j   = 0;                              //standard variable used in loops
    int k   = 0;                              //standard variable used in loops

    types_count[0] = p.max - p.min + 1;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function get the atom names for each lipid type                                                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void name_lip_atoms(Trajectory &traj,system_variables &s,program_variables &p,sv1d &types,sv2d &atom_names_lip)
{
    int i   = 0;                              //standard variable used in loops
    int j   = 0;                              //standard variable used in loops
    int k   = 0;                              //standard variable used in loops

    for(i=p.min; i<=p.max; i++) //loop over target atoms
    {
        atom_names_lip[0][i-p.min] = traj.atom_name[i];
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function records the statistics for the h-bonds                                                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void store_h_bond(Trajectory &traj,system_variables &s,program_variables &p,int acceptor,int donor,sv1d &types,
                  iv1d &p_atom_id,sv1d &l_atom_type,dv2d &freq,vector <vector <rvec*>> &coords,string current_lip_t,
                  int min,int max,iv1d &types_count)
{
    /*                                                                                                                                                                  
     *   Notes:                       
     *   The goal here is to store coordinates for each h-bond type and the frequency so that the average coordinates can be shown for each bond
     *   type and the data ordered by the frequency of occurance. The algorith is complicated by having multiple lipid types posible but is explained here.                                                                                               
     *   The strategy here is to check if the current h-bond has been encountered before. This is done by comparing the protein atom id with those
     *   stored in p_atom_id[]. The lipid atom type is also compared to those in l_atom_type[]. This means the exact protein atom is important but 
     *   only the atom type is important (could come different lipid molecules). It is also important to consider tht the lipid atom type could be 
     *   common for multiple lipid types. For example, the carbonyl carbons (C21 and C31) are shared for POPE and POPG lipids. Thus a new bond is not detected
     *   when C21/POPG is detected if C21/POPE has already been encountered. This is ok though since the coords[] and frequency data are given data 
     *   for all lipid types automatically when a new bond is found, i.e. the new bond creates an entry for that bond type between all lipid types 
     *   even if the lipid type does not have the correct atom type (its frequency would never exceed 0 and no choords would be stored). Each time a 
     *   new bond is encountered, the protein_id is added to p_atom_id[] and the lipid atom name is added to l_atom_type[]. The frequency data is also 
     *   added to freq[]. Initially, freq is set to 1 only for the current lipid type (freq[bond][type]). Coordinates are also stored for the current 
     *   lipid and the protein, coordinates of 0 are stored for positions in coords[bond][type][x,y,z] of the incorrect lipid type. This approach lets
     *   us separate the data by the acceptor and donor types and the lipid type. In the case that a new bond is not found, its already on the list for 
     *   p_atom_id[] and l_atom_type[], then the position of the bond is found in p_atom_id and stored as duplicate_index. The lipid types are then looped 
     *   over and checked against the current lipid type. When the correct lipid type is found, then freq[bond][type] is updated for that type only. 
     *   The coordinates are also added to coords[bond][type][x,y,z] for the correct lipid type. The final computation of frequency and average coords 
     *   is done in finalize_analysis(). Note the protein atom is assumed to be the acceptor in the code and the lipid the donor. When this is not true, 
     *   the order of these arguments are reversed when calling the function.     
     */

    int i                   = 0;                        //standard variable used in loops
    int j                   = 0;                        //standard variable used in loops
    int new_bond            = 1;                        //tells whether the h-bond type is new or not
    int duplicate_index     = 0;                        //stores the index for the current h-bond assuming it is already on the list
    int    this_p_atom_id   = traj.atom_nr[acceptor];   //atom_id for the current protein atom
    string this_l_atom_type = traj.atom_name[donor];    //atom name for the current lipid atom

    //check for a new h-bond 
    for(i=0; i<p_atom_id.size(); i++) //loop over current h-bonds
    {
        if(p_atom_id[i] == this_p_atom_id && strcmp(l_atom_type[i].c_str(), this_l_atom_type.c_str()) == 0)
        {
            new_bond        = 0;
            duplicate_index = i;
            goto lp_end;
        }
    }
    lp_end:;

    if(new_bond == 1)
    {
        //update the p_atom_nr
        p_atom_id.push_back(this_p_atom_id);

        //update the l_atom_type
        l_atom_type.push_back(this_l_atom_type);

        //update the frequency
        dv1d this_freq(types.size(),0.0);
        for(i=0; i<types.size(); i++) //loop over the types
        {
            if(strcmp(types[i].c_str(), current_lip_t.c_str()) == 0)
            {
                this_freq[i] = 1.0;
            }
        }
        freq.push_back(this_freq);

        //update coords
        vector <rvec*> this_coords;

        for(i=0; i<types.size(); i++) //loop over lipid types 
        {
            int size_lip  = types_count[i];
            int size_prot = traj.prot.size();
            int this_size = size_prot + size_lip;
            rvec *this_r;
            this_r = (rvec *)calloc(this_size , sizeof(*this_r));

            if(strcmp(types[i].c_str(), current_lip_t.c_str()) == 0) //lipid type is correct
            {
                //store protein coords
                for(j=0; j<size_prot; j++) //loop over protein atoms
                {
                    this_r[j][0] = traj.r[traj.prot[j]-1][0];
                    this_r[j][1] = traj.r[traj.prot[j]-1][1];
                    this_r[j][2] = traj.r[traj.prot[j]-1][2];
                }

                //store lipid coords
                int count = 0;
                for(j=min; j<=max; j++) //loop over lipid atoms
                {
                    this_r[size_prot+count][0] = traj.r[j][0];
                    this_r[size_prot+count][1] = traj.r[j][1];
                    this_r[size_prot+count][2] = traj.r[j][2];
                    count++;
                }
            }
            else //lipid type is not correct. store zeros for coords
            {
                //store protein coords
                for(j=0; j<this_size; j++) //loop over protein atoms
                {
                    this_r[j][0] = 0.0;
                    this_r[j][1] = 0.0;
                    this_r[j][2] = 0.0;
                }
            }

            //store coords
            this_coords.push_back(this_r);
        }

        //store coords
        coords.push_back(this_coords);
    }
    else
    {
        //update frequency and add the coordinates
        for(i=0; i<types.size(); i++) //loop over lipid types 
        {
            if(strcmp(types[i].c_str(), current_lip_t.c_str()) == 0)
            {
                freq[duplicate_index][i] = freq[duplicate_index][i] + 1.0;

                //add protein coords
                for(j=0; j<traj.prot.size(); j++) //loop over protein atoms
                {
                    coords[duplicate_index][i][j][0] = coords[duplicate_index][i][j][0] + traj.r[traj.prot[j]-1][0];
                    coords[duplicate_index][i][j][1] = coords[duplicate_index][i][j][1] + traj.r[traj.prot[j]-1][1];
                    coords[duplicate_index][i][j][2] = coords[duplicate_index][i][j][2] + traj.r[traj.prot[j]-1][2];
                }

                int count = 0;
                for(j=min; j<=max; j++) //loop over lipid atoms
                {
                    coords[duplicate_index][i][traj.prot.size()+count][0] = coords[duplicate_index][i][traj.prot.size()+count][0] + traj.r[j][0];
                    coords[duplicate_index][i][traj.prot.size()+count][1] = coords[duplicate_index][i][traj.prot.size()+count][1] + traj.r[j][1];
                    coords[duplicate_index][i][traj.prot.size()+count][2] = coords[duplicate_index][i][traj.prot.size()+count][2] + traj.r[j][2];
                    count++;
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function takes in the accpetor and donor atom types and tags the appropriate atoms for easy checks   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void tag_h_bonds_atoms(Trajectory &traj,system_variables &s,program_variables &p,Index &lip_a,Index &lip_d,
                       Index &prot_a,Index &prot_d,iv1d &donors,iv1d &acceptors,iv1d &refined_sel)
{
    /*                                                                                                                                                                  
     *   Notes:                                                                                                                                                           
     *   Here we read the acceptor and donor files and tag these atoms in the acceptors[] and donors[] structures. This approach lets us check whether 
     *   an atom is an acceptor or donor very quickly. The approach uses a single acceptors[]/donors[] and does not have separate structures for lipids
     *   and the protein. This is okay since our later loops will span either the protein or lipids separately. We can thus be certain that a pair is 
     *   select with an atom from both the protein and lipid. We also consider the refined selection as created using a selection text but only for the 
     *   protein atoms. The user can easily fine tune the lipids by varying the -lip_a and -lip_d selection cards.  
     */

    int i = 0;    //standard variable used in loops
    int j = 0;    //standard variable used in loops
    int k = 0;    //standard variable used in loops
    int l = 0;    //standard variable used in loops

    int pos = 0;
    string tag;

    //do lipid atoms
    for(i=0; i<traj.atoms(); i++) //loop over system atoms
    {
        //get upper and lower range for current residue
        int min = traj.get_res_start(i);
        int max = traj.get_res_end(i);

        //prime loop for the next residue
        i = traj.next_residue(i);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // set lipid donors                                                                                          //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        while(lip_d.check_next_tag(&pos,tag)) //loop over all tags
        {
            pos++;                           //set position to first item after tag

            if(strcmp(traj.res_name[min].c_str(), tag.c_str()) == 0) //lipid type is correct
            {
                int    next_pos = pos;        //position next tag. end of array if not found 
                string next_tag;              //stores the next tag if one is found

                lip_d.check_next_tag(&next_pos,next_tag); //find end of current section

                for(j=min; j<=max; j++) //loop over current residue atoms
                {
                    for(k=pos; k<next_pos; k++) //loop over items in current section
                    {
                        if(strcmp(traj.atom_name[j].c_str(), lip_d.index_s[k].c_str()) == 0) //atom is an donor lipid atom
                        {
                            donors[j] = 1;
                        }
                    }
                }
            }
        }
        pos = 0;

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // set lipid acceptors                                                                                       //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        while(lip_a.check_next_tag(&pos,tag)) //loop over all tags
        {
            pos++;                           //set position to first item after tag

            if(strcmp(traj.res_name[min].c_str(), tag.c_str()) == 0) //lipid type is correct
            {
                int    next_pos = pos;        //position next tag. end of array if not found 
                string next_tag;              //stores the next tag if one is found

                lip_a.check_next_tag(&next_pos,next_tag); //find end of current section

                for(j=min; j<=max; j++) //loop over current residue atoms
                {
                    for(k=pos; k<next_pos; k++) //loop over items in current section
                    {
                        if(strcmp(traj.atom_name[j].c_str(), lip_a.index_s[k].c_str()) == 0) //atom is an acceptor lipid atom
                        {
                            acceptors[j] = 1;
                        }
                    }
                }
            }
        }
        pos = 0;
    }

    //do protein atoms
    for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
    {
        //get the first and last atom of the current residue
        int min = traj.p_res_start(i);
        int max = traj.p_res_end(i);

        //jump to the next residue
        i = traj.next_prot_res(i);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // set protein donors                                                                                        //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        while(prot_d.check_next_tag(&pos,tag)) //loop over all tags
        {
            pos++;                           //set position to first item after tag

            if(strcmp(traj.res_name[min].c_str(), tag.c_str()) == 0) //residue type is correct
            {
                int    next_pos = pos;        //position next tag. end of array if not found 
                string next_tag;              //stores the next tag if one is found

                prot_d.check_next_tag(&next_pos,next_tag); //find end of current section

                for(j=min; j<=max; j++) //loop over current residue atoms
                {
                    for(k=pos; k<next_pos; k++) //loop over items in current section
                    {
                        if(strcmp(traj.atom_name[j].c_str(), prot_d.index_s[k].c_str()) == 0) //atom is an donor protein atom
                        {
                            if(refined_sel[j] == 1)
                            {
                                donors[j] = 1;
                            }
                        }
                    }
                }
            }
        }
        pos = 0;

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // set protein acceptors                                                                                     //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        while(prot_a.check_next_tag(&pos,tag)) //loop over all tags
        {
            pos++;                           //set position to first item after tag

            if(strcmp(traj.res_name[min].c_str(), tag.c_str()) == 0) //residue type is correct
            {
                int    next_pos = pos;        //position next tag. end of array if not found 
                string next_tag;              //stores the next tag if one is found

                prot_a.check_next_tag(&next_pos,next_tag); //find end of current section

                for(j=min; j<=max; j++) //loop over current residue atoms
                {
                    for(k=pos; k<next_pos; k++) //loop over items in current section
                    {
                        if(strcmp(traj.atom_name[j].c_str(), prot_a.index_s[k].c_str()) == 0) //atom is an acceptor protein atom
                        {
                            if(refined_sel[j] == 1)
                            {
                                acceptors[j] = 1;
                            }
                        }
                    }
                }
            }
        }
        pos = 0;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // write pdb files with acceptors and donors highlighted by the B factor                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    FILE *acceptor_donor_pdb_file;

    string acceptor_pdb_file_name = chop_and_add_tag(p.hbk_file_name,"_acceptors.pdb");
    string donor_pdb_file_name    = chop_and_add_tag(p.hbk_file_name,"_donors.pdb");

    if(s.world_rank == 0)
    {
        acceptor_donor_pdb_file = fopen(acceptor_pdb_file_name.c_str(), "w");
        if(acceptor_donor_pdb_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",acceptor_pdb_file_name.c_str());
        }
        else
        {
            for(i=0; i<traj.atoms(); i++) //loop over atoms
            {
                traj.beta[i] = acceptors[i];
            }
            write_frame_pdb(traj.box_ref,traj.atoms(),traj.atom_nr,traj.res_nr,traj.res_name,traj.atom_name,traj.r_ref,traj.title,s.world_rank,&acceptor_donor_pdb_file,traj.beta,traj.weight,traj.element,traj.chain_id,0);
            fclose(acceptor_donor_pdb_file);
        }

        acceptor_donor_pdb_file = fopen(donor_pdb_file_name.c_str(), "w");
        if(acceptor_donor_pdb_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",donor_pdb_file_name.c_str());
        }
        else
        {
            for(i=0; i<traj.atoms(); i++) //loop over atoms
            {
                traj.beta[i] = donors[i];
            }
            write_frame_pdb(traj.box_ref,traj.atoms(),traj.atom_nr,traj.res_nr,traj.res_name,traj.atom_name,traj.r_ref,traj.title,s.world_rank,&acceptor_donor_pdb_file,traj.beta,traj.weight,traj.element,traj.chain_id,0);
            fclose(acceptor_donor_pdb_file);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the contact matrix for the current trajectory frame                                //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void binding_contacts(Trajectory &traj,system_variables &s,program_variables &p,
		      Index &lip_a,Index &lip_d,Index &prot_a,Index &prot_d,iv2d &bonds,sv1d &types,iv1d &p_atom_id,
		      sv1d &l_atom_type,dv2d &freq,vector <vector <rvec*>> &coords,iv1d &types_count,dv1d &b_factor_freq,
		      iv1d &refined_sel,Index &report,Contacts &cont,iv1d &bound_lipids,Index &exclude,iv1d &acceptors,iv1d &donors)
{
    /*                                                                                                                                                                  
     *   Notes:                                                                                                                                                           
     *   We add a profile even if a lipid is not bound in the timeline. in this case it will be empty but is needed for reading back in contact profiles 
     *   since a global frame number is used to lookup a contact profile.
     */

    int i = 0;                             //standard variable used in loops
    int j = 0;                             //standard variable used in loops
    int k = 0;                             //standard variable used in loops

    iv2d this_profile(0,iv1d(0,0));        //hold contact info for the current frame

    if(p.b_be == 0 || find_lipid(traj,s,p,bound_lipids,traj.get_frame_global()) == 1)
    {
        for(i=p.min; i<=p.max; i++) //loop over target residue atoms
        {	
            int contacts = 0;

            //get lip_d-prot_a h-bonds
            if(donors[i] == 1) //lipid atom is a donor
            {
                int donor    = i;
                int acceptor = 0;
                int h        = 0;

                for(j=0; j<bonds[donor].size(); j++) //loop over bonds 
                {
                    if(traj.atom_name[bonds[donor][j]-1].at(0) == 'H') //atom is a hydrogen
                    {
                        h = bonds[donor][j]-1;

                        for(k=0; k<traj.prot.size(); k++) //loop over protein atoms
                        {
                            if(acceptors[traj.prot[k]-1] == 1)
                            {
                                acceptor = traj.prot[k]-1;
    
                                if(check_exclusions(traj,s,p,exclude,acceptor,donor) == 1)
                                {
                                    //check for h-bond
                                    int result = check_h_bond(traj,p,acceptor,donor,h,report);

                                    contacts = contacts + result;

                                    if(result == 1)
                                    {
                                        iv1d this_contact(2,0);                //store the atom id's for the current contact
                                        this_contact[0] = i-p.min;
                                        this_contact[1] = k;

                                        this_profile.push_back(this_contact);

                                        store_h_bond(traj,s,p,acceptor,donor,types,p_atom_id,l_atom_type,freq,coords,traj.res_name[i],p.min,p.max,types_count);
                                        b_factor_freq[acceptor] = b_factor_freq[acceptor] + 1.0/traj.get_ef_frames();
                                    }
                                }
                            }
                        }
                    }
                }
            }


//printf("acceptor[%d] %s \n",traj.atom_name[i].c_str());	    //get lip_a-prot_d h-bonds
            if(acceptors[i] == 1) //lipid atom is an acceptor
            {
                int donor    = 0;
                int acceptor = i;
                int h        = 0;

//printf("acceptor found \n");

                for(j=0; j<traj.prot.size(); j++) //loop over protein atoms
                {
                    if(donors[traj.prot[j]-1] == 1)
                    {
                        donor = traj.prot[j]-1;

                        if(check_exclusions(traj,s,p,exclude,acceptor,donor) == 1)
                        {
                            for(k=0; k<bonds[donor].size(); k++) //loop over bonds 
                            {
                                if(refined_sel[bonds[donor][k]-1] == 1)
                                {
                                    if(traj.atom_name[bonds[donor][k]-1].at(0) == 'H') //atom is a hydrogen
                                    {
                                        h = bonds[donor][k]-1;

                                        //check for h-bond
                                        int result = check_h_bond(traj,p,acceptor,donor,h,report);

                                        contacts = contacts + result;

                                        if(result == 1)
                                        {
//printf("b_bond found \n");
                                            iv1d this_contact(2,0);                //store the atom id's for the current contact
                                            this_contact[0] = i-p.min;
                                            this_contact[1] = j;

                                            this_profile.push_back(this_contact);

                                            store_h_bond(traj,s,p,donor,acceptor,types,p_atom_id,l_atom_type,freq,coords,traj.res_name[i],p.min,p.max,types_count);
                                            b_factor_freq[donor] = b_factor_freq[donor] + 1.0/traj.get_ef_frames();
                                        }
                                    }
                                }
                            }
                        }
                    }
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
// Analyze contact matrices and compile data                                                                 //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,dv1d &x_coord,
                         dv1d &y_coord,dv1d &z_coord,dv1d &prot_x,dv1d &prot_y,dv1d &prot_z,
                         iv3d &filter,iv2d &bound,dv2d &dwell_t,dv2d &events,Contacts &cont,iv3d &filter_history,
			 iv1d &bound_lipids)
{
    /*                                                                                                                                                                  
     *   Notes:                                                                                                                                                           
     *   The strategy used here is to check for a new bound lipid using the timeline and then reinitialize the filter to 0 when a new lipid is found.
     *   this approach ensures that dwell times are recorded when a new lipid binds since the old contacts are kicked off. Since a contact profile is 
     *   made for every frame (some are empty if a target lipid was not bound) we do not need to screen for a bound lipid on the main loop that reads 
     *   and merges temporary contact profiles. 
     */

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

    for(i=0; i<traj.get_ef_frames(); i++)
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

        string histo_file_name = chop_and_add_tag(p.hbk_file_name,"_histo.dat");

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
        string pdb_file_name = chop_and_add_tag(p.hbk_file_name,".pdb");
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
        string contacts_file_name = chop_and_add_tag(p.hbk_file_name,"_hbonds_dwell_t.pml");
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

                if(dash_rad >= 0.0000005) //printf %f gives 6 decimal places. anything above 0.0000005 rounds up to 0.000001 giving a non zero radius. PyMOL shows a zero radius too big (maybe a default value)
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
        contacts_file_name = chop_and_add_tag(p.hbk_file_name,"_hbonds_freq.pml");
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

                if(dash_rad >= 0.0000005) //printf %f gives 6 decimal places. anything above 0.0000005 rounds up to 0.000001 giving a non zero radius. PyMOL shows a zero radius too big (maybe a default value)
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
// Collect the contacts and compute the average                                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis_cont(Trajectory &traj,system_variables &s,program_variables &p,sv1d &types,iv1d &p_atom_id,
                              sv1d &l_atom_type,dv2d &freq,vector <vector <rvec*>> &coords,iv1d &types_count,sv2d &atom_names_lip,
                              dv1d &b_factor_freq)
{
    /*                                                                                                                                                                  
     *   Notes:                                                                                                                                                           
     *   Here we collect data for the h-bonds statistics that was recorded in the store_h_bond() function. (like fin_ana() in lipid_h_bonds.cpp but no grid 
     *   stuff to finish). This data is collected using basic mpi_send/receive functions that are designed here to handle the complex form of 
     *   the data structures. Each data is stored by rank 0 in a new structure (this_p_atom_id[], this_l_atom_type[], this_freq[], and this_coords). 
     *   When collecting data, the new entries are added to these structures. Afterwards the data in this_p_atom_id[] and this_l_atom_type[] are compared
     *   to entries in p_atom_id[] and l_atom_type[] to see if a new bond was found for one of the non-rank 0 cores. If a new bond is found, then that entry 
     *   is added to the original structures (p_atom_id[], l_atom_type[], freq[], and coords[]). If the bond was already on the list, then we simply add
     *   this_freq[] ti freq[] and this_coords[] to coords[] at the correct index. This finally gives us a complete list for p_atom_id[], l_atom_type[], freq[], 
     *   and coords[]. Next, the data is normalized and written to output files. Here the coords are divided by data in freq[]. Before writing any data, the 
     *   data is organized by the freq[] so larger freqs come first. This is done for each lipid type independantly since output is written for each type 
     *   separately. For ordering the data we create a copy of the freq data (order_freq[]) and copy the initial index (order[]). We then shuffle both so 
     *   that order_freq[] is organized with larger numbers first. The data in order then points to a location in p_atom_id[] or coords[] that is used to 
     *   retreive coords when writing to the pdb file. Next, the data needed for a pdb is generated (res_name, element, B factor etc.) and the pdb is written. 
     *   Here, we loop over the bonds using order[] to pull data from freq[] and coords[]. For each bond, we set the B factor to the freq[]. And finally, the 
     *   data is written to the pdb file. Note, this whole process is in a loop over the types. So each lipid type receives its own pdb. In addition to the pdb
     *   file, a text file is generated that gives PyMOL select commands for each bond in the pdb file. And finally, in the last section, the frequency data is 
     *   collected for how often a protein atoms participates in h-bonds. This data is written to the B factor of a pdb with coords from the last analyzed frame
     *   by rank 0.      
     */

    int    i        = 0;                      //standard variable used in loops
    int    j        = 0;                      //standard variable used in loops
    int    k        = 0;                      //standard variable used in loops
    int    l        = 0;                      //standard variable used in loops

    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nFinalizing H-bond statistical analysis. This requires communicating data and could take some time. \n");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze h-bonds statistics                                                                                //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //collect p_atom_id
    if(s.world_size > 0)
    {
        if(s.world_rank == 0)
        {
            printf("\nCollecting H-bond statistics data. \n");
        }

        for(i=1; i<s.world_size; i++)
        {
            if(s.world_rank == 0)
            {
                printf("  Working on rank %d. \n",i);
            }

            iv1d this_p_atom_id(0,0);             //holds the p_atom_id data received
            sv1d this_l_atom_type(0,"");          //holds the l_atom_type data received
            dv2d this_freq(0,dv1d(0,0.0));        //holds the freq data received
            vector <vector <rvec*>> this_coords;  //holds the coords data received

            //collect p_atom_id
            if(s.world_rank == 0)
            { 
                int size = 0;
                MPI_Recv(&size, 1, MPI_INT, i, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                int rcv[size];            
                MPI_Recv(rcv, size, MPI_INT, i, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                for(j=0; j<size; j++) //loop over received atom_ids 
                {
                    this_p_atom_id.push_back(rcv[j]);
                }
            }
            else if(s.world_rank == i)
            {
                int size = p_atom_id.size(); 
                MPI_Send(&size, 1, MPI_INT, 0, 13, MPI_COMM_WORLD);
                int snd[size];
                for(j=0; j<size; j++)
                {
                    snd[j] = p_atom_id[j];
                }
                MPI_Send(snd, size, MPI_INT, 0, 13, MPI_COMM_WORLD);
            }

            //collect l_atom_type
            if(s.world_rank == 0)
            {
                int size = 0;
                MPI_Recv(&size, 1, MPI_INT, i, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                for(j=0; j<size; j++) //loop over each string in list
                {
                    int string_size = 0;
                    MPI_Recv(&string_size, 1, MPI_INT, i, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    char this_string[string_size];
                    MPI_Recv(this_string, string_size, MPI_CHAR, i, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    
                    string current_string = "";
                    for(k=0; k<string_size; k++) //loop over current string
                    {
                        current_string = current_string + this_string[k];
                    }

                    this_l_atom_type.push_back(current_string);
                }
            }
            else if(s.world_rank == i)
            {
                int size = l_atom_type.size();
                MPI_Send(&size, 1, MPI_INT, 0, 13, MPI_COMM_WORLD);

                for(j=0; j<size; j++) //loop over each string in list
                {
                    int string_size = l_atom_type[j].length();
                    MPI_Send(&string_size, 1, MPI_INT, 0, 13, MPI_COMM_WORLD);
                    char this_string[string_size];
                    
                    for(k=0; k<string_size; k++) //loop over current string
                    {
                        this_string[k] = l_atom_type[j][k];
                    }
                    MPI_Send(this_string, string_size, MPI_CHAR, 0, 13, MPI_COMM_WORLD);
                }
            }

            //collect freq
            if(s.world_rank == 0)
            {
                int size_x = 0;
                int size_y = 0;
                MPI_Recv(&size_x, 1, MPI_INT, i, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&size_y, 1, MPI_INT, i, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                double this_freq_ary[size_x][size_y];
                MPI_Recv(this_freq_ary, size_x*size_y, MPI_DOUBLE, i, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                for(j=0; j<size_x; j++) //loop over the h-bonds
                {
                    dv1d current_freq(size_y,0.0);
                    for(k=0; k<size_y; k++) //loop over the lipid types
                    {
                        current_freq[k] = this_freq_ary[j][k];
                    }
                    this_freq.push_back(current_freq);
                }
            }
            else if(s.world_rank == i)
            {
                int size_x = freq.size(); 
                int size_y = types.size();
                int size   = size_x*size_y;
                MPI_Send(&size_x, 1, MPI_INT, 0, 13, MPI_COMM_WORLD);
                MPI_Send(&size_y, 1, MPI_INT, 0, 13, MPI_COMM_WORLD);
                double this_freq[size_x][size_y];
                for(j=0; j<size_x; j++) //loop over the h-bonds
                {
                    for(k=0; k<size_y; k++) //loop over the lipid types
                    {
                        this_freq[j][k] = freq[j][k];
                    }
                }
                MPI_Send(this_freq, size, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD);
            }

            //collect coords
            int my_size = coords.size();
            int world_size[s.world_size];
            MPI_Allgather(&my_size, 1, MPI_INT, world_size, 1, MPI_INT,MPI_COMM_WORLD);

            for(j=0; j<world_size[i]; j++) //loop over current cores h-bonds
            {
                vector <rvec*> this_bond_coords(0);

                for(k=0; k<types.size(); k++) //loop over lipid types
                {
                    if(s.world_rank == 0)
                    {
                        int num_atoms = 0; 
                        MPI_Recv(&num_atoms, 1, MPI_INT, i, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                        double current_coords[num_atoms][3];
                        MPI_Recv(current_coords, num_atoms*3, MPI_DOUBLE, i, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
 
                        rvec *current_coords_rvec;
                        current_coords_rvec = (rvec *)calloc(num_atoms , sizeof(*current_coords_rvec));
                        for(l=0; l<num_atoms; l++) //loop over atoms
                        {
                            current_coords_rvec[l][0] = current_coords[l][0];
                            current_coords_rvec[l][1] = current_coords[l][1];
                            current_coords_rvec[l][2] = current_coords[l][2];
                        }
                        this_bond_coords.push_back(current_coords_rvec);
                    }
                    else if(s.world_rank == i)
                    {
                        int num_atoms = traj.prot.size() + types_count[k];
                        MPI_Send(&num_atoms, 1, MPI_INT, 0, 13, MPI_COMM_WORLD);
                        double current_coords[num_atoms][3];
                        for(l=0; l<num_atoms; l++) //loop over atoms
                        {
                            current_coords[l][0] = coords[j][k][l][0];
                            current_coords[l][1] = coords[j][k][l][1];
                            current_coords[l][2] = coords[j][k][l][2];
                        }
                        MPI_Send(current_coords, num_atoms*3, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD);
                    }
                }
                if(s.world_rank == 0) 
                {
                    this_coords.push_back(this_bond_coords);
                }
            }

            //add new h-bond data to the list
            if(s.world_rank == 0)
            {
                for(j=0; j<this_p_atom_id.size(); j++) //loop over new h-bonds
                {
                    int new_bond = 1;
                    int pos      = -1;

                    for(k=0; k<p_atom_id.size(); k++) //loop over old h-bonds
                    {
                        if(p_atom_id[k] == this_p_atom_id[j] && strcmp(l_atom_type[k].c_str(), this_l_atom_type[j].c_str()) == 0)
                        {
                            new_bond = 0; 
                            pos = k;                           
                        }
                    }

                    if(new_bond == 1)
                    {
                        p_atom_id.push_back(this_p_atom_id[j]);
                        l_atom_type.push_back(this_l_atom_type[j]);
                        freq.push_back(this_freq[j]);

                        vector <rvec*> current_bond_coords;
                        for(k=0; k<types.size(); k++) //loop over lipid types
                        {
                            int num_atoms = traj.prot.size()+types_count[k];
                            rvec *this_r;
                            this_r = (rvec *)calloc(num_atoms , sizeof(*this_r));
                            for(l=0; l<num_atoms; l++) //loop over atoms
                            {
                                this_r[l][0] = this_coords[j][k][l][0];
                                this_r[l][1] = this_coords[j][k][l][1];
                                this_r[l][2] = this_coords[j][k][l][2];
                            }
                            current_bond_coords.push_back(this_r);
                        }
                        coords.push_back(current_bond_coords);
                    }
                    else
                    {
                        for(k=0; k<types.size(); k++) //loop over lipid types
                        {
                            freq[pos][k] = freq[pos][k] + this_freq[j][k];  
                        }
                        for(k=0; k<types.size(); k++) //loop over lipid types
                        {
                            for(l=0; l<traj.prot.size()+types_count[k]; l++) //loop over atoms
                            {
                                coords[pos][k][l][0] = coords[pos][k][l][0] + this_coords[j][k][l][0];
                                coords[pos][k][l][1] = coords[pos][k][l][1] + this_coords[j][k][l][1];
                                coords[pos][k][l][2] = coords[pos][k][l][2] + this_coords[j][k][l][2];
                            }
                        } 
                    }
                }
            }
        }
    }

    if(s.world_rank == 0)
    {
        printf("\nNormalizing statistical data and writing output files. \n");

        //normalize coords
        int size_prot = traj.prot.size();
        for(i=0; i<p_atom_id.size(); i++) //loop over h-bonds
        {
            for(j=0; j<types.size(); j++) //loop over lipid types 
            {
                int size_lip  = types_count[j];
                int size_prot = traj.prot.size();
                int this_size = size_prot + size_lip;

                if(freq[i][j] > 0.0)
                {
                    for(k=0; k<this_size; k++) //loop over atoms
                    {
                        coords[i][j][k][0] = coords[i][j][k][0]/freq[i][j];
                        coords[i][j][k][1] = coords[i][j][k][1]/freq[i][j];
                        coords[i][j][k][2] = coords[i][j][k][2]/freq[i][j];
                    }
                }
            }
        }

        //write coords to pdb
        for(i=0; i<types.size(); i++) //loop over lipid types 
        {
            int size_lip  = types_count[i];
            int size_prot = traj.prot.size();
            int this_size = size_prot + size_lip;

            //allocate memory for shuffling data
            iv1d order(p_atom_id.size(),0);
            dv1d order_freq(p_atom_id.size(),0.0);
            for(j=0; j<p_atom_id.size(); j++) //loop over h-bonds
            {
                order[j] = j;
                order_freq[j] = freq[j][i];
            }

            //shuffle data by the freq value
            int b_swap = 0;
            for(b_swap=1; b_swap > 0; )
            {
                b_swap = 0;
                for(j=0; j<p_atom_id.size()-1; j++)
                {
                    if(order_freq[j] < order_freq[j+1])
                    {
                        int order_freq_tmp = order_freq[j];
                        int order_tmp      = order[j];
                        order_freq[j]      = order_freq[j+1];
                        order[j]           = order[j+1];
                        order_freq[j+1]    = order_freq_tmp; 
                        order[j+1]         = order_tmp;
                        b_swap             = 1; 
                    }
                }
            }

            //allocate memory to store info needed for writing out a pdb
            iv1d    this_atom_nr(this_size,0);                  //atom number used in pdb file
            iv1d    this_res_nr(this_size,1);                   //res number used in pdb file
            sv1d    this_atom_name(this_size);                  //atom name used in pdb file
            sv1d    this_res_name(this_size,types[i]);          //res name used in pdb file
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
            for(j=0; j<types_count[i]; j++) //loop over lipid atoms
            {
                this_atom_name[traj.prot.size()+j] = atom_names_lip[i][j];
            }

            //set the residue names
            for(j=0; j<traj.prot.size(); j++) //loop over protein atoms
            {
                this_res_name[j] = traj.res_name[traj.prot[j]-1];
            }
            for(j=0; j<types_count[i]; j++) //loop over lipid atoms
            {
                this_res_name[traj.prot.size()+j] = types[i];
            }

            //set residue numbers
            for(j=0; j<traj.prot.size(); j++) //loop over protein atoms
            {
                this_res_nr[j] = traj.res_nr[traj.prot[j]-1];
            }
            for(j=0; j<types_count[i]; j++) //loop over lipid atoms
            {
                this_res_nr[traj.prot.size()+j] = this_res_nr[traj.prot.size()-1] + 1;
            }        

            //open pdb file
            FILE *pdb_file;
            string tag = "_" + types[i] + ".pdb";
            string pdb_file_name = chop_and_add_tag(p.hbk_file_name,tag.c_str()); 
            pdb_file = fopen(pdb_file_name.c_str(), "w");
            if(pdb_file == NULL)
            {
                printf("failure opening %s. Make sure the file exists. \n",pdb_file_name.c_str());
            }
            else
            {
                printf("  Writing projections of h-bond statistics onto time average coordinates to %s. \n",pdb_file_name.c_str());

                for(j=0; j<p_atom_id.size(); j++) //loop over h-bonds
                {
                    if(freq[order[j]][i] > 0.0)
                    {
                        //set B-factor to the freq
                        for(k=0; k<traj.prot.size()+types_count[i]; k++) //loop over protein atoms
                        {
                            this_beta[k] = freq[order[j]][i]/traj.get_ef_frames();
                        }

                        //write the average coords to a pdb file
                        write_frame_pdb(traj.ibox,this_size,this_atom_nr,this_res_nr,this_res_name,this_atom_name,coords[order[j]][i],traj.title,s.world_rank,&pdb_file,this_beta,this_weight,this_element,this_chain_id,j);
                    }
                }
                
                //close file
                fclose(pdb_file);
            }

            //report frequency data
            FILE *freq_file;
            string freq_tag = "_" + types[i] + "_freq.dat";
            string freq_file_name = chop_and_add_tag(p.hbk_file_name,freq_tag.c_str());
            freq_file = fopen(freq_file_name.c_str(), "w");
            if(freq_file == NULL)
            {
                printf("failure opening %s. Make sure the file exists. \n",freq_file_name.c_str());
            }
            else
            {
                printf("  Writing h-bond statistics to %s. \n",freq_file_name.c_str());

                fprintf(freq_file," %12s %12s %-s \n","h-bond","frequency","pymol_selection");
                for(j=0; j<p_atom_id.size(); j++) //loop over h-bonds
                {
                    string select_tag = "select hb, name " + l_atom_type[order[j]] + " or (resid " + to_string(traj.res_nr[p_atom_id[order[j]]-1]) + " and name " + traj.atom_name[p_atom_id[order[j]]-1] + " )";
                    fprintf(freq_file," %12d %12f %-s \n",j+1,freq[order[j]][i]/traj.get_ef_frames(),select_tag.c_str());
                }
                fclose(freq_file);
            }
        }
    }

    //write h-bond frequencies for each protein atom to the b-factor
    //collect B-factor data, i.e., the h-bond frequency for each protein atom
    collect_and_sum_dv1d(s.world_size,s.world_rank,b_factor_freq);

    if(s.world_rank == 0)
    {
        //open pdb file for writing
        FILE *b_factor_file;
        string b_factor_tag = "_freq.pdb";
        string b_factor_file_name = chop_and_add_tag(p.hbk_file_name,b_factor_tag.c_str());
        b_factor_file = fopen(b_factor_file_name.c_str(), "w");
        if(b_factor_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",b_factor_file_name.c_str());
        }
        else
        {
            printf("  Writing h-bond frequency data for each protein atom to %s. \n",b_factor_file_name.c_str());

            write_frame_pdb(traj.ibox,traj.atoms(),traj.atom_nr,traj.res_nr,traj.res_name,traj.atom_name,traj.r,traj.title,s.world_rank,&b_factor_file,b_factor_freq,traj.weight,traj.element,traj.chain_id,1);
            fclose(b_factor_file);
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
    s.program_name = "H-Bond Kinetics";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                          s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                          s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                         s.world_rank, s.cl_tags, &p.b_print,     0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                           s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                                      s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                       s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                         s.world_rank, s.cl_tags, &p.b_lsq,       0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                           s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",           s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_s(argc,argv,"-lip_a",  p.lip_a_file_name,            "Selection card with lipid h-bond acceptor atom types (crd)",          s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_s(argc,argv,"-lip_d",  p.lip_d_file_name,            "Selection card with lipid h-bond donor atom types (crd)",             s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_s(argc,argv,"-prot_a", p.prot_a_file_name,           "Selection card with protein h-bond acceptor atom types (crd)",        s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_s(argc,argv,"-prot_d", p.prot_d_file_name,           "Selection card with protein h-bond donor atom types (crd)",           s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_s(argc,argv,"-bond",   p.bond_file_name,             "Selection card identifying protein bonding pairs (crd)",              s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_s(argc,argv,"-hbk",    p.hbk_file_name,              "Output file with h-bonds kinetics data (dat)",                        s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_s(argc,argv,"-pf_pdb", p.pf_pdb_file_name,           "PDB file with selected protein (pdb)",                                s.world_rank, s.cl_tags, &p.b_pf_pdb,    0);
    add_argument_mpi_s(argc,argv,"-pf_prm", p.protein_finder_param_name,  "File with additional protein finder parameters (prm)",                s.world_rank, s.cl_tags, &p.b_pf_param,  0);
    add_argument_mpi_i(argc,argv,"-test",   &p.b_test,                    "Print info for checking hydrogen bonds? (0:no 1:yes)",                s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_s(argc,argv,"-sel",    p.selection_text_file_name,   "Input file with the atom selection text (sel)",                       s.world_rank, s.cl_tags, &p.b_sel_text,  0);
    add_argument_mpi_s(argc,argv,"-report", p.report_file_name,           "Selection card with atom ids for h-bond to report when found",        s.world_rank, s.cl_tags, &p.b_report,    0);
    add_argument_mpi_i(argc,argv,"-range",  &p.range,                     "Noise filter half-width ",                                            s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_i(argc,argv,"-dump",   &p.dump,                      "Dump contacts on the last trajectory frame? (0:no, 1:yes) ",          s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_d(argc,argv,"-dt",     &p.delta_t,                   "Effective time step between trajectory frames analyzed (ps) ",        s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_d(argc,argv,"-bin",    &p.bin_width,                 "Residence time histogram bin width (ps) ",                            s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_i(argc,argv,"-resi",   &p.resi,                      "Target residue for measuring contacts ",                              s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_d(argc,argv,"-max",    &p.max_dash_rad,              "Maximum thickness of dash for PyMOL distance commands ",              s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_d(argc,argv,"-min",    &p.min_dash_rad,              "Minimum thickness of dash for PyMOL distance commands ",              s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_s(argc,argv,"-be",     p.be_file_name,               "Input binding events file (be)",                                      s.world_rank, s.cl_tags, &p.b_be,        0);
    add_argument_mpi_i(argc,argv,"-x",      &p.target_x,                  "The target x lattice point used with the binding events file",        s.world_rank, s.cl_tags, &p.b_x,         0);
    add_argument_mpi_i(argc,argv,"-y",      &p.target_y,                  "The target y lattice point used with the binding events file",        s.world_rank, s.cl_tags, &p.b_y,         0);
    add_argument_mpi_s(argc,argv,"-type",   p.target_res,                 "Lipid type to select from the bound lipids timeline",                 s.world_rank, s.cl_tags, &p.b_target_res,0);
    add_argument_mpi_s(argc,argv,"-ex",     p.exclude_file_name,          "Selection card with interaction pairs to exclude (crd)",              s.world_rank, s.cl_tags, &p.b_exclude,   0);
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
    check_extension_mpi(s.world_rank,"-hbk",p.hbk_file_name,".dat");
    check_extension_mpi(s.world_rank,"-lip_a",p.lip_a_file_name,".crd");
    check_extension_mpi(s.world_rank,"-lip_d",p.lip_d_file_name,".crd");
    check_extension_mpi(s.world_rank,"-prot_a",p.prot_a_file_name,".crd");
    check_extension_mpi(s.world_rank,"-prot_d",p.prot_d_file_name,".crd");
    check_extension_mpi(s.world_rank,"-bond",p.prot_d_file_name,".crd");

    if(p.b_be == 1)
    {
        check_extension_mpi(s.world_rank,"-be",p.be_file_name,".be");
    }
    if(p.b_exclude == 1)
    {
        check_extension_mpi(s.world_rank,"-ex",p.exclude_file_name,".crd");
    }
    if(p.b_pf_pdb == 1)
    {
        check_extension_mpi(s.world_rank,"-pf_pdb",p.pf_pdb_file_name,".pdb");
    }
    if(p.b_pf_param == 1)
    {
        check_extension_mpi(s.world_rank,"-pf_prm",p.protein_finder_param_name,".prm");
    }
    if(p.b_sel_text == 1)
    {
        check_extension_mpi(s.world_rank,"-sel",p.selection_text_file_name,".sel");
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

    //create index objects
    Index lip_a;
    Index lip_d;
    Index prot_a;
    Index prot_d;
    Index bond;
    Index report;
    Index exclude; 

    //read the index files
    lip_a.get_index(p.lip_a_file_name);
    lip_d.get_index(p.lip_d_file_name);
    prot_a.get_index(p.prot_a_file_name);
    prot_d.get_index(p.prot_d_file_name);
    bond.get_index(p.bond_file_name);
    if(p.b_report == 1)
    {
        report.get_index(p.report_file_name);
    }
    if(p.b_exclude == 1)
    {
        exclude.get_index(p.exclude_file_name);
    }

    //run proten finder
    traj.get_protein(p.protein_finder_param_name,p.b_pf_param);

    //print a pdb with distinguished protein
    traj.write_protein(p.pf_pdb_file_name,p.b_pf_pdb);

    //print info about the protein
    traj.get_prot_stats();

    //use bond list to create a list of bonds for each atom. Easier to look up this way.  
    iv2d bonds(traj.atoms(),iv1d(0,0));
    get_bonds(traj,s,p,bond,bonds);

    //set the parallelization scheme to the protein atoms for analyzing noise filter
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
    cont.init(p.hbk_file_name,p.size_y,traj.prot.size());

    //open files for writing temporary contact profiles data
    cont.prime_tmp_file();

    //allocate memory to hold h-bond stats
    sv1d types(1);                                  //stores the lipid types being analyzed
    iv1d types_count(types.size(),0);               //stores the number of atoms for each lipid type
    iv1d p_atom_id(0,0);                            //stores the protein atom id for each h-bond
    sv1d l_atom_type(0);                            //stores the lipid atom type for each h-bond 
    dv2d freq(0,dv1d(types.size(),0.0));            //stores the frequency that the h-bond occurs in the trajectory
    vector <vector <rvec*>> coords;                 //stores the atomic coordinates of the protein and a lipid

    //set the lipid type
    types[0] = traj.res_name[p.min];

    //count the atoms in each lipid type
    count_lipid_atoms(traj,s,p,types,types_count);

    //allocate memory to store atom names for lipids and protein
    sv2d atom_names_lip(0,sv1d(0));
    int i = 0;
    for(i=0; i<types.size(); i++) //loop over lipid types 
    {
        sv1d this_atom_names(types_count[i],"");
        atom_names_lip.push_back(this_atom_names);
    }

    //get lipid atom names
    name_lip_atoms(traj,s,p,types,atom_names_lip);

    dv1d b_factor_freq(traj.atoms(),0.0); //store the number of times a protein atom formed an h-bond as b-factor

    //create object to hold protein atom refinement
    iv1d refined_sel(traj.atoms(),1);

    if(p.b_sel_text == 1)
    {
        //create a object to hold an atom selection
        Selection this_sel;

        //select the atoms
        this_sel.get_selection(traj,p.selection_text_file_name);

        //generate pdb file name for highlighting the selection 
        string pdb_filename = chop_and_add_tag(p.selection_text_file_name,".pdb");

        //highlight the selected atoms
        this_sel.highlight_sel(traj,pdb_filename);

        //refine the selection
        refined_sel = this_sel.tag_atoms(traj);
    }

    iv1d donors(traj.atoms(),0);       //tags all donor atoms
    iv1d acceptors(traj.atoms(),0);    //tags all acceptor atoms 

    //tag the acceptor and donor atoms
    tag_h_bonds_atoms(traj,s,p,lip_a,lip_d,prot_a,prot_d,donors,acceptors,refined_sel);

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

	binding_contacts(traj,s,p,lip_a,lip_d,prot_a,prot_d,bonds,types,p_atom_id,l_atom_type,freq,coords,types_count,b_factor_freq,refined_sel,report,cont,bound_lipids,exclude,acceptors,donors);

        get_mean_coords(traj,s,p,x_coord,y_coord,z_coord,prot_x,prot_y,prot_z,bound_lipids);

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //close tmp contacts files
    cont.close_tmp_file();

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect contacts from mpi processes and compute the average
    perf.log_time(finalize_analysis(traj,s,p,x_coord,y_coord,z_coord,prot_x,prot_y,prot_z,filter,bound,dwell_t,events,cont,filter_history,bound_lipids),"Fin Ana");

    //print h-bond statistics to file
    perf.log_time(finalize_analysis_cont(traj,s,p,types,p_atom_id,l_atom_type,freq,coords,types_count,atom_names_lip,b_factor_freq),"Fin Ana cont.");

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

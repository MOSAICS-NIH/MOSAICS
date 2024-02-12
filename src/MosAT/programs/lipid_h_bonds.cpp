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
#include "MosAT/program_variables/pv_lipid_h_bonds.h"        //This has the variables specific to the analysis program
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

    int result = 0;               //tells if the pair under investigation make a hydrogen bond
    int pi     = 3.1415926535;    //its pi!
    int i      = 0;               //standard variable used in loops

    rvec m;                       //difference vector between the acceptor and donor atoms
    rvec n;                       //difference vector between the donor and hydrogen atoms

    for(i=0; i<3; i++) //loop over 3 dimensions
    {
        m[i] = traj.r[acceptor][i] - traj.r[donor][i];
        n[i] = traj.r[h][i]        - traj.r[donor][i];
    }

    double angle = gmx_angle(m,n);                            //angle between acceptor, donor, and hydrogen atoms
    double dist  = sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);   //distance between acceptor and donor atoms

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
    int i = 0;                      //standard variable used in loops
    int j = 0;                      //standard variable used in loops
    int k = 0;                      //standard variable used in loops

    for(i=0; i<types.size(); i++) //loop over lipid types
    {
        for(j=0; j<traj.target_leaflet.size(); j++) //loop over target leaflet atoms
        {
            //get the first and last atom of the current lipid
            int min = traj.t_lip_start(j);
            int max = traj.t_lip_end(j);

            //jump to the next lipid
            j = traj.next_target_lipid(j);

            if(strcmp(traj.res_name[min].c_str(), types[i].c_str()) == 0) //lipid type is correct
            {
                types_count[i] = max - min + 1;
                goto next_type;
            }
        }
        next_type:;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function get the atom names for each lipid type                                                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void name_lip_atoms(Trajectory &traj,system_variables &s,program_variables &p,sv1d &types,sv2d &atom_names_lip)
{
    int i = 0;                      //standard variable used in loops
    int j = 0;                      //standard variable used in loops
    int k = 0;                      //standard variable used in loops

    for(i=0; i<types.size(); i++) //loop over lipid types
    {
        for(j=0; j<traj.target_leaflet.size(); j++) //loop over target leaflet atoms
        {
            //get the first and last atom of the current lipid
            int min = traj.t_lip_start(j);
            int max = traj.t_lip_end(j);

            //jump to the next lipid
            j = traj.next_target_lipid(j);

            if(strcmp(traj.res_name[min].c_str(), types[i].c_str()) == 0) //lipid type is correct
            {
                for(k=min; k<=max; k++) //loop over current lipid atoms
                {
                    atom_names_lip[i][k-min] = traj.atom_name[k];
                }
                goto next_type;
            }
        }
        next_type:;
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

        //add frequency data for the new bond. set to 1 for the correct lipid type
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
        vector <rvec*> this_coords; //holds the coordinates of the current lipid/protein
 
        for(i=0; i<types.size(); i++) //loop over lipid types 
        {
            int size_lip  = types_count[i];        //size of the lipid type 
            int size_prot = traj.prot.size();      //size of the protein
            int this_size = size_prot + size_lip;  //overall size
            rvec *this_r;                          //hold coordinates for the current lipid type/protein combo
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
    for(i=0; i<traj.target_leaflet.size(); i++) //loop over target leaflet atoms
    {
        //get the first and last atom of the current lipid
        int min = traj.t_lip_start(i);
        int max = traj.t_lip_end(i);

        //jump to the next lipid
        i = traj.next_target_lipid(i);

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

    string acceptor_pdb_file_name = chop_and_add_tag(p.lphb_file_name,"_acceptors.pdb");
    string donor_pdb_file_name    = chop_and_add_tag(p.lphb_file_name,"_donors.pdb");

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
// This function computes the number of lipid-prot h-bonds and adds it to the grid                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void lip_h_bonds(Trajectory &traj,system_variables &s,program_variables &p,Index &param,Index &lip_a,Index &lip_d,
                 Index &prot_a,Index &prot_d,iv2d &bonds,Grid &hb,sv1d &types,iv1d &p_atom_id,sv1d &l_atom_type,
                 dv2d &freq,vector <vector <rvec*>> &coords,iv1d &types_count,dv1d &b_factor_freq,iv1d &refined_sel,
		 Index &report,Index &exclude,iv1d &acceptors,iv1d &donors)
{
    /*                                                                                                                                                                  
     *   Notes:                                                                                                                                                           
     *   Here we loop over the lipids and count hydrogen bonds with the protein. The routine works by scanning over the lipids and checking each 
     *   against the lipids listed in param[] (-crd). If a target lipid is found, then the number of hydrogen bonds is counted for it and eventually
     *   stamped around the mapping atoms specified in param[] (-crd). The routine has two branches of code. First, we scan over the current lipid atoms 
     *   and check if the atoms are an acceptor. If so, then we loop over the protein atoms and check if any of them are a donor. The second branch 
     *   checks if the lipid atom is a donor and then loops over the protein atoms and checks if they are an acceptor. Both branches loop over the bonded
     *   atoms of the donor and look for hydrogen atoms. Once a acceptor, donor, and hydrogen are identified, their indices are passed into the 
     *   check_h_bond() function. If the bond checks out, then the number of bonds for the current lipid is increased. Similarly, the B factor is increased for 
     *   the protien atom participating in the bonds. This will be used to determine the frequency that the protein atom form h-bonds. The statistics are 
     *   also recored for the h-bond using the store_h_bond() function. This records the time average coords for each h-bond type and arranges the data by 
     *   the frequency of occurance.    
     */

    int    i        = 0;                      //standard variable used in loops
    int    j        = 0;                      //standard variable used in loops
    int    k        = 0;                      //standard variable used in loops
    int    l        = 0;                      //standard variable used in loops
    int    m        = 0;                      //standard variable used in loops
    int    n        = 0;                      //standard variable used in loops
    double hx       = 0;                      //head atom x-component 
    double hy       = 0;                      //head atom y-component

    //clear the current frame grids
    hb.clean_frame();

    for(i=0; i<traj.target_leaflet.size(); i++) //loop over target leaflet atoms
    {
        //get the first and last atom of the current lipid
        int min_l = traj.t_lip_start(i);
        int max_l = traj.t_lip_end(i);

        //jump to the next lipid
        i = traj.next_target_lipid(i);

        for(j=0; j<param.index_s.size(); j+=3) //loop over lipid types
        {
            if(strcmp(traj.res_name[min_l].c_str(), param.index_s[j].c_str()) == 0) //lipid type is correct
            {
                int contacts = 0;  //number of h-bonds for current lipid

                for(k=min_l; k<=max_l; k++) //loop over current lipid atoms
                {
                    if(acceptors[k] == 1) //lipid atom is an acceptor
                    {
                        int acceptor = k; 

                        for(l=0; l<traj.prot.size(); l++) //loop over protein atoms
                        {
                            //get the first and last atom of the current residue
                            int min_p = traj.p_res_start(l);
                            int max_p = traj.p_res_end(l);

                            //jump to the next residue
                            l = traj.next_prot_res(l);        

                            for(m=min_p; m<=max_p; m++) //loop over current residue atoms
                            {
                                if(donors[m] == 1)
                                {
                                    int donor = m;
  
                                    for(n=0; n<bonds[donor].size(); n++) //loop over bonds 
                                    {
                                        if(traj.atom_name[bonds[donor][n]-1].at(0) == 'H') //atom is a hydrogen
                                        {
                                            int h = bonds[donor][n]-1;

                                            if(check_exclusions(traj,s,p,exclude,acceptor,donor) == 1)
                                            {
                                                //check for h-bond
                                                int result = check_h_bond(traj,p,acceptor,donor,h,report);

                                                contacts = contacts + result;

                                                if(result == 1)
                                                {
                                                    store_h_bond(traj,s,p,donor,acceptor,types,p_atom_id,l_atom_type,freq,coords,traj.res_name[k],min_l,max_l,types_count);
                                                    b_factor_freq[donor] = b_factor_freq[donor] + 1.0/traj.get_ef_frames();
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if(donors[k] == 1) //lipid atom is a donor
                    {
                        int donor = k;

                        for(l=0; l<traj.prot.size(); l++) //loop over protein atoms
                        {
                            //get the first and last atom of the current residue
                            int min_p = traj.p_res_start(l);
                            int max_p = traj.p_res_end(l);

                            //jump to the next residue
                            l = traj.next_prot_res(l);

                            for(m=min_p; m<=max_p; m++) //loop over current residue atoms
                            {
                                if(acceptors[m] == 1)
                                {
                                    int acceptor = m;

                                    for(n=0; n<bonds[donor].size(); n++) //loop over bonds 
                                    {
                                        if(traj.atom_name[bonds[donor][n]-1].at(0) == 'H') //atom is a hydrogen
                                        {
                                            int h = bonds[donor][n]-1;

                                            if(check_exclusions(traj,s,p,exclude,acceptor,donor) == 1)
                                            {
                                                //check for h-bond
                                                int result = check_h_bond(traj,p,acceptor,donor,h,report);

                                                contacts = contacts + result;

                                                if(result == 1)
                                                {
                                                    store_h_bond(traj,s,p,acceptor,donor,types,p_atom_id,l_atom_type,freq,coords,traj.res_name[k],min_l,max_l,types_count);
                                                    b_factor_freq[acceptor] = b_factor_freq[acceptor] + 1.0/traj.get_ef_frames();
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                //add the number of contacts to the grid
                for(k=min_l; k<=max_l; k++) //loop over current lipid atoms
                {
                    if(strcmp(traj.atom_name[k].c_str(), param.index_s[j+1].c_str()) == 0) //mapping atom 1
                    {
                        hx = traj.r[k][0]; 
                        hy = traj.r[k][1];

                        hb.stamp(hx,hy,p.radius,contacts);
                    }
                    else if(strcmp(traj.atom_name[k].c_str(), param.index_s[j+2].c_str()) == 0) //mapping atom 2
                    {
                        hx = traj.r[k][0];
                        hy = traj.r[k][1];

                        hb.stamp(hx,hy,p.radius,contacts);
                    }
                }
            }
        }
    }

    //get the average for the current frame
    hb.norm_frame();

    //add the current frame grid to long term sum
    hb.add_frame();

    //now we print the single frame lpsb data
    if(p.b_stdev == 1)
    {
        hb.write_frame(traj.get_frame_global());
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Collect the contacts and compute the average                                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,Grid &hb,sv1d &types,iv1d &p_atom_id,
                         sv1d &l_atom_type,dv2d &freq,vector <vector <rvec*>> &coords,iv1d &types_count,sv2d &atom_names_lip,
                         dv1d &b_factor_freq)
{
    /*                                                                                                                                                                  
     *   Notes:                                                                                                                                                           
     *   Here we wrap up grid based analysis for the h-bonds like usual. Then, we collect data for the h-bonds statistics that was recorded in the 
     *   store_h_bond() function. This data is collected using basic mpi_send/receive functions that are designed here to handle the complex form of 
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
        printf("\nFinalizing analysis. This requires communicating data over the grid and could take some time depending on the resolution. \n");
    }

    //collect hb and rho from all ranks
    hb.collect_grid();

    //normalize hb and write the <hb> and rho to file
    if(s.world_rank == 0)
    {
        hb.normalize();

        hb.exclude_data(p.cutoff,1);

        hb.write_grid();
        hb.write_rho();
    }

    //compute standard deviation
    hb.get_stdev(p.b_stdev,p.b_clean,traj);

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
                printf("Working on rank %d. \n",i);
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
        printf("\nNormalizing data and writing output files. \n\n");

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
            string pdb_file_name = chop_and_add_tag(p.lphb_file_name,tag.c_str()); 
            pdb_file = fopen(pdb_file_name.c_str(), "w");
            if(pdb_file == NULL)
            {
                printf("failure opening %s. Make sure the file exists. \n",pdb_file_name.c_str());
            }
            else
            {
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
            string freq_file_name = chop_and_add_tag(p.lphb_file_name,freq_tag.c_str());
            freq_file = fopen(freq_file_name.c_str(), "w");
            if(freq_file == NULL)
            {
                printf("failure opening %s. Make sure the file exists. \n",freq_file_name.c_str());
            }
            else
            {
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
    if(s.world_rank == 0)
    {
        printf("Writing h-bond frequency data for each protein atom. \n");
    }

    //collect B-factor data, i.e., the h-bond frequency for each protein atom
    collect_and_sum_dv1d(s.world_size,s.world_rank,b_factor_freq);

    if(s.world_rank == 0)
    {
        //open pdb file for writing
        FILE *b_factor_file;
        string b_factor_tag = "_freq.pdb";
        string b_factor_file_name = chop_and_add_tag(p.lphb_file_name,b_factor_tag.c_str());
        b_factor_file = fopen(b_factor_file_name.c_str(), "w");
        if(b_factor_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",b_factor_file_name.c_str());
        }
        else
        {
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
    s.program_name = "Lipid H-Bonds";

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
    add_argument_mpi_s(argc,argv,"-crd",    p.param_file_name,            "Lipid types selection card + mapping atoms (crd)",             s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lip_a",  p.lip_a_file_name,            "Selection card with lipid h-bond acceptor atom types (crd)",   s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lip_d",  p.lip_d_file_name,            "Selection card with lipid h-bond donor atom types (crd)",      s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-prot_a", p.prot_a_file_name,           "Selection card with protein h-bond acceptor atom types (crd)", s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-prot_d", p.prot_d_file_name,           "Selection card with protein h-bond donor atom types (crd)",    s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-bond",   p.bond_file_name,             "Selection card identifying protein bonding pairs (crd)",       s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lphb",   p.lphb_file_name,             "Output file with spatially resolved H-bond count (dat)",       s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,           "PDB file with sorted leaflets (pdb)",                          s.world_rank, s.cl_tags, &p.b_lf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,  "File with additional leaflet finder parameters (prm)",         s.world_rank, s.cl_tags, &p.b_lf_param,0);
    add_argument_mpi_s(argc,argv,"-pf_pdb", p.pf_pdb_file_name,           "PDB file with selected protein (pdb)",                         s.world_rank, s.cl_tags, &p.b_pf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-pf_prm", p.protein_finder_param_name,  "File with additional protein finder parameters (prm)",         s.world_rank, s.cl_tags, &p.b_pf_param,0);
    add_argument_mpi_s(argc,argv,"-sf_pdb", p.sf_pdb_file_name,           "PDB file with selected sol (pdb)",                             s.world_rank, s.cl_tags, &p.b_sf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-sf_prm", p.solvent_finder_param_name,  "File with additional solvent finder parameters (prm)",         s.world_rank, s.cl_tags, &p.b_sf_param,0);
    add_argument_mpi_d(argc,argv,"-APS",    &p.APS,                       "Area per grid box (nm^2)",                                     s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-r",      &p.radius,                    "Radius of the target atom (nm)",                               s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-cutoff", &p.cutoff,                    "Cutoff for excluding data (chi)",                              s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-bx",     &p.box_x,                     "Grid x dimension (nm)",                                        s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-by",     &p.box_y,                     "Grid y dimension (nm)",                                        s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-leaf",   &p.leaflet,                   "Which leaflet? (0:both 1:upper 2:lower)",                      s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-odf",    &p.out_data_format,           "What is the output data format? (0:matrix 1:vector)",          s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-stdev",  &p.b_stdev,                   "Compute the STDEV and STEM? (0:no 1:yes)",                     s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-clean",  &p.b_clean,                   "Remove single frame files? (0:no 1:yes)",                      s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-test",   &p.b_test,                    "Print info for checking hydrogen bonds? (0:no 1:yes)",         s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-sel",    p.selection_text_file_name,   "Input file with the atom selection text (sel)",                s.world_rank, s.cl_tags, &p.b_sel_text,0);
    add_argument_mpi_s(argc,argv,"-report", p.report_file_name,           "Selection card with atom ids for h-bond to report when found", s.world_rank, s.cl_tags, &p.b_report,  0);
    add_argument_mpi_s(argc,argv,"-ex",     p.exclude_file_name,          "Selection card with interaction pairs to exclude (crd)",       s.world_rank, s.cl_tags, &p.b_exclude, 0);
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
    check_extension_mpi(s.world_rank,"-crd",p.param_file_name,".crd");
    check_extension_mpi(s.world_rank,"-lip_a",p.lip_a_file_name,".crd");
    check_extension_mpi(s.world_rank,"-lip_d",p.lip_d_file_name,".crd");
    check_extension_mpi(s.world_rank,"-prot_a",p.prot_a_file_name,".crd");
    check_extension_mpi(s.world_rank,"-prot_d",p.prot_d_file_name,".crd");
    check_extension_mpi(s.world_rank,"-bond",p.prot_d_file_name,".crd");

    if(p.b_exclude == 1)
    {
        check_extension_mpi(s.world_rank,"-ex",p.exclude_file_name,".crd");
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
    if(p.b_sf_pdb == 1)
    {
        check_extension_mpi(s.world_rank,"-sf_pdb",p.sf_pdb_file_name,".pdb");
    }
    if(p.b_sf_param == 1)
    {
        check_extension_mpi(s.world_rank,"-sf_prm",p.solvent_finder_param_name,".prm");
    }
    if(p.b_sel_text == 1)
    {
        check_extension_mpi(s.world_rank,"-sel",p.selection_text_file_name,".sel");
    }
    if(p.b_report == 1)
    {
        check_extension_mpi(s.world_rank,"-report",p.report_file_name,".crd");
    }

    //create index objects
    Index lip_a;  
    Index lip_d; 
    Index prot_a;
    Index prot_d;  
    Index bond;
    Index param;
    Index report;
    Index exclude;

    //read the index files
    lip_a.get_index(p.lip_a_file_name);
    lip_d.get_index(p.lip_d_file_name);
    prot_a.get_index(p.prot_a_file_name);
    prot_d.get_index(p.prot_d_file_name);
    bond.get_index(p.bond_file_name);
    param.get_index(p.param_file_name);
    if(p.b_report == 1)
    {
        report.get_index(p.report_file_name);
    }
    if(p.b_exclude == 1)
    {
        exclude.get_index(p.exclude_file_name);
    }

    //run leaflet/proten/solvent finder
    traj.get_leaflets(p.leaflet,p.leaflet_finder_param_name,p.b_lf_param);
    traj.get_protein(p.protein_finder_param_name,p.b_pf_param);
    traj.get_solvent(p.solvent_finder_param_name,p.b_sf_param);

    //print a pdb with distinguished leaflets/protein/solvent
    traj.write_leaflets(p.lf_pdb_file_name,p.b_lf_pdb);
    traj.write_protein(p.pf_pdb_file_name,p.b_pf_pdb);
    traj.write_sol(p.sf_pdb_file_name,p.b_sf_pdb);

    //print info about the target leaflet
    traj.get_leaflet_stats();

    //print info about the lipid selection
    traj.get_lipid_selection_stats(param.get_column_s(3,0),"-crd");

    //print info about the protein
    traj.get_prot_stats();

    //print info about the water
    traj.get_sol_stats();

    //create a grid to hold lipid contacts
    Grid hb;

    //get the grid dimensions
    hb.get_dim(p.box_x,p.box_y,traj.ibox,p.APS);

    //set the output file name for grid
    hb.set_output(p.lphb_file_name,p.out_data_format);

    //print info about the grid
    hb.print_dim();

    //use bond list to create a list of bonds for each atom. Easier to look up this way.  
    iv2d bonds(traj.atoms(),iv1d(0,0));
    get_bonds(traj,s,p,bond,bonds);

    //allocate memory to hold h-bond stats
    sv1d types = param.get_column_s(3,0);           //stores the lipid types being analyzed
    iv1d types_count(types.size(),0);               //stores the number of atoms for each lipid type
    iv1d p_atom_id(0,0);                            //stores the protein atom id for each h-bond
    sv1d l_atom_type(0);                            //stores the lipid atom type for each h-bond 
    dv2d freq(0,dv1d(types.size(),0.0));            //stores the frequency that the h-bond occurs in the trajectory
    vector <vector <rvec*>> coords;                 //stores the atomic coordinates of the protein and a lipid
    int i = 0;                                      //standard variable used in loops

    //count the atoms in each lipid type
    count_lipid_atoms(traj,s,p,types,types_count);

    //allocate memory to store atom names for lipids and protein
    sv2d atom_names_lip(0,sv1d(0));
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

        lip_h_bonds(traj,s,p,param,lip_a,lip_d,prot_a,prot_d,bonds,hb,types,p_atom_id,l_atom_type,freq,coords,types_count,b_factor_freq,refined_sel,report,exclude,acceptors,donors);

        traj.set_beta_lf();

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect contacts from mpi processes and compute the average
    perf.log_time(finalize_analysis(traj,s,p,hb,types,p_atom_id,l_atom_type,freq,coords,types_count,atom_names_lip,b_factor_freq),"Fin Ana");

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

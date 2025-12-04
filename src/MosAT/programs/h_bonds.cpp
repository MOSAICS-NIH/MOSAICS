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
#include "MosAT/program_variables/pv_h_bonds.h"              //This has the variables specific to the analysis program
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
int check_h_bond(Trajectory &traj,program_variables &p,int acceptor,int donor,int h)
{
    /*                                                                                                                                                                  
     *   Notes:                                                                                                                                                           
     *   Here we take in an atom index for the acceptor, donor, and hydrogen atoms. We then check if there is a hydrogen bond formed between them. The  
     *   function returns 1 if a h-bond is present. The function also prints PyMOL select commands used to visualize the bond in PyMOL when the -test 
     *   option is used. 
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
    }

    return result; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function takes in the accpetor and donor atom types and tags the appropriate atoms for easy checks   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void tag_h_bonds_atoms(Trajectory &traj,system_variables &s,program_variables &p,Index &acc,Index &don,iv1d &donors,
		       iv1d &acceptors,Index &n1,string &file_name)
{
    int i = 0;    //standard variable used in loops
    int j = 0;    //standard variable used in loops
    int k = 0;    //standard variable used in loops
    int l = 0;    //standard variable used in loops

    int pos = 0;
    string tag;

    iv1d index_tag(traj.atoms(),0);    //tag each atom in index for fast checking
    for(i=0; i<n1.index_s.size(); i++)
    {
        index_tag[n1.index_i[i]-1] = 1;
    }

    //do protein atoms
    for(i=0; i<traj.atoms(); i++) //loop over atoms
    {
        //get the first and last atom of the current residue
        int min = traj.get_res_start(i);
        int max = traj.get_res_end(i);

        //jump to the next residue
        i = traj.next_residue(i); 

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // set donors                                                                                                //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        while(don.check_next_tag(&pos,tag)) //loop over all tags
        {
            pos++;                           //set position to first item after tag

            if(strcmp(traj.res_name[min].c_str(), tag.c_str()) == 0) //residue type is correct
            {
                int    next_pos = pos;        //position next tag. end of array if not found 
                string next_tag;              //stores the next tag if one is found

                don.check_next_tag(&next_pos,next_tag); //find end of current section

                for(j=min; j<=max; j++) //loop over current residue atoms
                {
                    for(k=pos; k<next_pos; k++) //loop over items in current section
                    {
                        if(strcmp(traj.atom_name[j].c_str(), don.index_s[k].c_str()) == 0) //atom is an donor atom
                        {
                            if(index_tag[j] == 1)
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
        // set acceptors                                                                                             //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        while(acc.check_next_tag(&pos,tag)) //loop over all tags
        {
            pos++;                           //set position to first item after tag

            if(strcmp(traj.res_name[min].c_str(), tag.c_str()) == 0) //residue type is correct
            {
                int    next_pos = pos;        //position next tag. end of array if not found 
                string next_tag;              //stores the next tag if one is found

                acc.check_next_tag(&next_pos,next_tag); //find end of current section

                for(j=min; j<=max; j++) //loop over current residue atoms
                {
                    for(k=pos; k<next_pos; k++) //loop over items in current section
                    {
                        if(strcmp(traj.atom_name[j].c_str(), acc.index_s[k].c_str()) == 0) //atom is an acceptor atom
                        {
                            if(index_tag[j] == 1)
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

    string acceptor_pdb_file_name = chop_and_add_tag(file_name,"_acceptors.pdb");
    string donor_pdb_file_name    = chop_and_add_tag(file_name,"_donors.pdb");

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
// This function checks if a set of a,d,h has been encountered befor to prevent double counting              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int check_pairs(int a,int d,int h,iv1d &a_found,iv1d &d_found,iv1d &h_found)
{
    int i = 0;
    int found = 0;    

    for(i=0; i<a_found.size(); i++)
    {
        if(a == a_found[i] && d == d_found[i] && h == h_found[i])
        {
            found = 1; 
        }
    }

    return found; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the number of lipid-prot h-bonds and adds it to the grid                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_h_bonds(Trajectory &traj,system_variables &s,program_variables &p,Index &acc,Index &don,iv2d &bonds,
		 Index &exclude,iv1d &acceptors_1,iv1d &donors_1,iv1d &acceptors_2,iv1d &donors_2,iv1d &h_bond_count)
{
    int    i        = 0;                      //standard variable used in loops
    int    j        = 0;                      //standard variable used in loops
    int    k        = 0;                      //standard variable used in loops
    int contacts    = 0;                      //number of h-bonds for current frame

    iv1d a_found(0,0); //used to keep track of h-bonds encountered so we dont double count
    iv1d d_found(0,0); //used to keep track of h-bonds encountered so we dont double count
    iv1d h_found(0,0); //used to keep track of h-bonds encountered so we dont double count

    for(i=0; i<traj.atoms(); i++) //loop over atoms
    {
        if(acceptors_1[i] == 1) //residue atom is an acceptor
        {
            int acceptor = i;

            for(j=0; j<traj.atoms(); j++) //loop over atoms
            {
                if(donors_2[j] == 1)
                {
                    int donor = j;

                    if(donor != acceptor)
                    {    
                        for(k=0; k<bonds[donor].size(); k++) //loop over bonds 
                        {
                            if(traj.atom_name[bonds[donor][k]-1].at(0) == 'H') //atom is a hydrogen
                            {
                                int h = bonds[donor][k]-1;

                                if(check_exclusions(traj,s,p,exclude,acceptor,donor) == 1)
                                {
                                    if(check_h_bond(traj,p,acceptor,donor,h) == 1)
                                    {
                                        if(check_pairs(acceptor,donor,h,a_found,d_found,h_found) == 0)
                                        {
                                            a_found.push_back(acceptor);
                                            d_found.push_back(donor);
                                            h_found.push_back(h);

                                            contacts = contacts + 1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if(donors_1[i] == 1) //atom is a donor
        {
            int donor = i;
 
            for(j=0; j<traj.atoms(); j++) //loop over atoms
            {
                if(acceptors_2[j] == 1)
                {
                    int acceptor = j;

                    if(donor != acceptor)
                    {
                        for(k=0; k<bonds[donor].size(); k++) //loop over bonds 
                        {
                            if(traj.atom_name[bonds[donor][k]-1].at(0) == 'H') //atom is a hydrogen
                            {
                                int h = bonds[donor][k]-1;

                                if(check_exclusions(traj,s,p,exclude,acceptor,donor) == 1)
                                {
                                    if(check_h_bond(traj,p,acceptor,donor,h) == 1)
                                    {
                                        if(check_pairs(acceptor,donor,h,a_found,d_found,h_found) == 0)
                                        {
                                            a_found.push_back(acceptor);
                                            d_found.push_back(donor);
                                            h_found.push_back(h);
                                        
                                            contacts = contacts + 1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    h_bond_count.push_back(contacts);

    //printf("contacts %d \n",contacts); 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Collect the contacts and compute the average                                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,iv1d &h_bond_count)  
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
        printf("\nFinalizing analysis. This requires communicating h-bond data and could take some time. \n");
    }

    collect_iv1d(s.world_size,s.world_rank,h_bond_count);

    if(s.world_rank == 0)
    {
        FILE *hb_file = fopen(p.hb_file_name.c_str(),"w");
        fprintf(hb_file," %9s   %10s \n","#step","#h_bonds");
        for(i=0; i<h_bond_count.size(); i++)
        {
            fprintf(hb_file," %9d   %10d \n",i,h_bond_count[i]);
        }
        fclose(hb_file);
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
    s.program_name = "H-Bonds";

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
    add_argument_mpi_s(argc,argv,"-acc",    p.acc_file_name,              "Selection card with h-bond acceptor atom types (crd)",         s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-don",    p.don_file_name,              "Selection card with h-bond donor atom types (crd)",            s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-bond",   p.bond_file_name,             "Selection card identifying protein bonding pairs (crd)",       s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-hb",     p.hb_file_name,               "Output file with H-bond count (dat)",                          s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-test",   &p.b_test,                    "Print info for checking hydrogen bonds? (0:no 1:yes)",         s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-ex",     p.exclude_file_name,          "Selection card with interaction pairs to exclude (crd)",       s.world_rank, s.cl_tags, &p.b_exclude, 0);
    add_argument_mpi_s(argc,argv,"-n1",     p.n1_file_name,               "Index for group 1 atoms (ndx)",                                s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-n2",     p.n2_file_name,               "Index for group 2 atoms (ndx)",                                s.world_rank, s.cl_tags, nullptr,      1);
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
    check_extension_mpi(s.world_rank,"-acc",p.acc_file_name,".crd");
    check_extension_mpi(s.world_rank,"-don",p.don_file_name,".crd");
    check_extension_mpi(s.world_rank,"-bond",p.bond_file_name,".crd");

    if(p.b_exclude == 1)
    {
        check_extension_mpi(s.world_rank,"-ex",p.exclude_file_name,".crd");
    }

    //create index objects
    Index n1;                //holds atom group 1 
    Index n2;                //holds atom group 2
    Index acc;               //holds acceptor atom types
    Index don;               //holds donor atom types
    Index bond;              //holds bonds list
    Index exclude;

    //read the index files
    n1.get_index(p.n1_file_name);
    n2.get_index(p.n2_file_name);
    acc.get_index(p.acc_file_name);
    don.get_index(p.don_file_name);
    bond.get_index(p.bond_file_name);
    if(p.b_exclude == 1)
    {
        exclude.get_index(p.exclude_file_name);
    }

    //use bond list to create a list of bonds for each atom. Easier to look up this way.  
    iv2d bonds(traj.atoms(),iv1d(0,0));
    get_bonds(traj,s,p,bond,bonds);

    //create structures for tagging acceptor and donor atoms
    iv1d donors_1(traj.atoms(),0);       //tags all donor atoms for group 1
    iv1d acceptors_1(traj.atoms(),0);    //tags all acceptor atoms for group 1
    iv1d donors_2(traj.atoms(),0);       //tags all donor atoms for group 2
    iv1d acceptors_2(traj.atoms(),0);    //tags all acceptor atoms for group 2

    //tag the acceptor and donor atoms
    tag_h_bonds_atoms(traj,s,p,acc,don,donors_1,acceptors_1,n1,p.n1_file_name);
    tag_h_bonds_atoms(traj,s,p,acc,don,donors_2,acceptors_2,n2,p.n2_file_name);

    iv1d  h_bond_count(0,0);   //store number of h-bonds for each frame

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

        get_h_bonds(traj,s,p,acc,don,bonds,exclude,acceptors_1,donors_1,acceptors_2,donors_2,h_bond_count);

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect contacts from mpi processes and compute the average
    perf.log_time(finalize_analysis(traj,s,p,h_bond_count),"Fin Ana");

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

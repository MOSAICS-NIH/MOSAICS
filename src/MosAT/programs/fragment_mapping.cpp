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
#include "MosAT/program_variables/pv_fragment_mapping.h"     //This has the variables specific to the analysis program
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function checks if a string is a filename                                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int check_for_filename(string this_string,string extension)
{
    int i                = 0;                     //standard variable used in loops
    int this_string_size = this_string.length();  //how long is the filename
    int extension_size   = extension.length();    //how long is the extension
    int is_filename      = 1;                     //was a filename provided

    if(extension_size >= this_string_size)
    {
        is_filename = 0;
    }
    else
    {
        for(i=0; i<extension_size; i++)
        {
            if(this_string[this_string_size-1-i] != extension[extension_size-1-i])
            {
		is_filename = 0;
            }
        }
    }
    return is_filename;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function finds the geometric center of the atom selection                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
dv1d get_center(Trajectory &traj,system_variables &s,program_variables &p,Param &param,dv1d &this_x,dv1d &this_y,dv1d &this_z)
{
    int i = 0;                //standard variable used in loops
    double cx = 0.0;          //center in x
    double cy = 0.0;          //center in y
    double cz = 0.0;          //center in z
    dv1d   r_center(3,0.0);   //the geometric center

    for(i=0; i<this_x.size(); i++) //loop over molecule atoms
    {
        cx = cx + this_x[i];
        cy = cy + this_y[i];
        cz = cz + this_z[i];
    }
    r_center[0] = cx/(double)this_x.size();
    r_center[1] = cy/(double)this_x.size();
    r_center[2] = cz/(double)this_x.size();

    return r_center;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function counts the number of molecules of the target type                                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int count_molecules(Trajectory &traj,system_variables &s,program_variables &p,Param &param)
{
    int i = 0;             //standard variable used in loops
    int count = 0;         //how many molecules of the target type

    for(i=0; i<traj.atoms(); i++) //loop over system atoms
    {
        //get the first and last atom of the current residue
        int min = traj.get_res_start(i);
        int max = traj.get_res_end(i);

        //jump to the next residue
        i = traj.next_residue(i);

        if(strcmp(traj.res_name[min].c_str(), p.target.c_str() ) == 0) //residue matches target type 
        {
            count++;
	}
    }
 
    return count; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function counts the number of atoms in the target molecule                                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int count_atoms(Trajectory &traj,system_variables &s,program_variables &p,Param &param)
{
    int i = 0;             //standard variable used in loops
    int count = 0;         //how many atoms in the target molecule

    for(i=0; i<param.main_size_y(); i++) //loop over fragments
    {
        count = count + param.sec_size_y(i);
    }

    return count; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function rotates the fragment to minimize distance between connecting atoms                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void rotate_fragment(Trajectory &traj,system_variables &s,program_variables &p,Param &param,dv1d &this_x,dv1d &this_y,
		     dv1d &this_z,int bond_aa_index,int bond_cg_index,sv1d &this_frag,dv1d &r_map)
{
    //Note: The strategy is as follows. We store the fragment coords in an rvec that can be used with GROMACS routines. 
    //We then store the reference molecule in an rvec. In this case, the reference only contains 1 atom, however, the 
    //revec must match the one for the fragment in size. We thus add dummy atoms to the ref rvec with bogus coords in 
    //addition to the real atom. The real atom is then selected by setting its mass value to 1 and all other atoms to 0. 
    //The bonding atoms mass in the fragment rvec is set to 1. This approach fits the fragment to the mapping atom while 
    //minimizing distance between the connecting atoms. 

    int i = 0;                //standard variable used in loops

    rvec *r_tmp;                          //holds current fragment coords
    rvec *r_ref_tmp;                      //holds ref coords
    int frag_size = this_x.size();        //how many atoms in current fragment

    //allocate memory for holding coords
    r_tmp     = (rvec*)calloc(frag_size , sizeof(rvec));
    r_ref_tmp = (rvec*)calloc(frag_size , sizeof(rvec));

    //get fragment and ref coords 
    for(i=0; i<this_x.size(); i++) //loop over fragment atoms
    {
        r_tmp[i][0] = this_x[i];
        r_tmp[i][1] = this_y[i];
        r_tmp[i][2] = this_z[i];

        r_ref_tmp[i][0] = 0.0;
        r_ref_tmp[i][1] = 0.0;
        r_ref_tmp[i][2] = 0.0;
    }

    //set coordinates for the cg-mapping atom
    r_ref_tmp[bond_aa_index][0] = traj.r[bond_cg_index][0];
    r_ref_tmp[bond_aa_index][1] = traj.r[bond_cg_index][1];
    r_ref_tmp[bond_aa_index][2] = traj.r[bond_cg_index][2];

    //compute shifting factor 
    double dx = -r_map[0];
    double dy = -r_map[1];
    double dz = -r_map[2];

    //center fragment at origin 
    for(i=0; i<this_x.size(); i++) //loop over fragment atoms
    {
        r_tmp[i][0] = r_tmp[i][0] + dx;
        r_tmp[i][1] = r_tmp[i][1] + dy;
        r_tmp[i][2] = r_tmp[i][2] + dz;

        r_ref_tmp[i][0] = r_ref_tmp[i][0] + dx;
        r_ref_tmp[i][1] = r_ref_tmp[i][1] + dy;
        r_ref_tmp[i][2] = r_ref_tmp[i][2] + dz;
    } 

    //tag atoms for fitting
    real   dummy_mass[frag_size];
    for(i=0; i<frag_size; i++)
    {
        dummy_mass[i] = 0.0;
    }
    dummy_mass[bond_aa_index] = 1.0;

    //do lsq fitting
    do_fit_ndim(3, frag_size, dummy_mass, r_ref_tmp, r_tmp);

    //move fragment center back to original position
    for(i=0; i<this_x.size(); i++) //loop over fragment atoms
    {
        r_tmp[i][0] = r_tmp[i][0] - dx;
        r_tmp[i][1] = r_tmp[i][1] - dy;
        r_tmp[i][2] = r_tmp[i][2] - dz;

        r_ref_tmp[i][0] = r_ref_tmp[i][0] - dx;
        r_ref_tmp[i][1] = r_ref_tmp[i][1] - dy;
        r_ref_tmp[i][2] = r_ref_tmp[i][2] - dz;
    }

    //copy rotated coords back to vectors
    for(i=0; i<this_x.size(); i++) //loop over fragment atoms
    {
        this_x[i] = r_tmp[i][0];
        this_y[i] = r_tmp[i][1];
        this_z[i] = r_tmp[i][2];
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function uses distance measurements to make a bonds list                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void frag_mapping(Trajectory &traj,system_variables &s,program_variables &p,Param &param,sv2d &mapping_atoms)
{
    int     i = 0;     //standard variable used in loops
    int     j = 0;     //standard variable used in loops
    int     k = 0;     //standard variable used in loops
    int     l = 0;     //standard variable used in loops
    int count_res = -1;    //Index the residues
    int count_frag = 0;    //Index the fragments

    //create rvec to hold all-atom coordinates
    rvec *r_aa;
    r_aa = (rvec*)calloc(p.num_atoms*p.num_mols, sizeof(rvec));

    for(i=0; i<traj.atoms(); i++) //loop over system atoms
    {
        //get the first and last atom of the current residue
        int min = traj.get_res_start(i);
        int max = traj.get_res_end(i);

        //jump to the next residue
        i = traj.next_residue(i);

        if(strcmp(traj.res_name[min].c_str(), p.target.c_str() ) == 0) //residue matches target type 
        {
            //count the residues as they are encountered
            count_res++;
            count_frag = 0;

            for(j=0; j<param.main_size_y(); j++) //loop over fragments
            {
                sv1d this_frag = param.get_column_sec_s(j,0); //get atom names for current fragment
                dv1d this_x    = param.get_column_sec_d(j,1); //get atom x-coord for current fragment
                dv1d this_y    = param.get_column_sec_d(j,2); //get atom y-coord for current fragment
                dv1d this_z    = param.get_column_sec_d(j,3); //get atom z-coord for current fragment

                //convert units to nm
                for(k=0; k<this_x.size(); k++) //loop over atoms
                {
                    this_x[k] = this_x[k]/10.0;
                    this_y[k] = this_y[k]/10.0;
                    this_z[k] = this_z[k]/10.0;
		}

                dv1d this_center = get_center(traj,s,p,param,this_x,this_y,this_z); //compute the fragments's geometric center

                //define variables for holing mapping atoms coords and index for connecting(bond) atoms
                int bond_cg_index = -1;     //index for cg-bonding atom
                int bond_aa_index = -1;     //index for frag-bonding atom
                dv1d r_map(3,0.0);          //holds coordinates of mapping atom or center
                dv1d map_center(3,0.0);     //geometric center of mapping atoms

                //find coords for mapping atom or the center when multiple map atoms are specified
                for(l=0; l<mapping_atoms[j].size(); l++) //loop over mapping atoms
                {
                    int found = 0;          //tells if the mapping atom was found

                    for(k=min; k<=max; k++) //loop over cg atoms
                    {
                        if(strcmp(traj.atom_name[k].c_str(), mapping_atoms[j][l].c_str() ) == 0) //cg atom matches mapping atom for fragment
                        {
                            map_center[0] = map_center[0] + traj.r[k][0];
                            map_center[1] = map_center[1] + traj.r[k][1];
                            map_center[2] = map_center[2] + traj.r[k][2];
                            
                            found = 1; 
                        }
                    }

                    if(found == 0)
                    {
                        printf("mapping atom not found %s \n",mapping_atoms[j][l].c_str());
                        MPI_Finalize();
                        exit(EXIT_SUCCESS);
		    }
                }
                r_map[0] = map_center[0]/(double)mapping_atoms[j].size();
                r_map[1] = map_center[1]/(double)mapping_atoms[j].size();
                r_map[2] = map_center[2]/(double)mapping_atoms[j].size();

                //compute shifting factor 
                double dx = r_map[0] - this_center[0];
                double dy = r_map[1] - this_center[1];
                double dz = r_map[2] - this_center[2];

                //shift the fragment
                for(k=0; k<this_x.size(); k++) //loop over atoms
                {
                    this_x[k] = this_x[k] + dx;
                    this_y[k] = this_y[k] + dy;
                    this_z[k] = this_z[k] + dz;
                }

                //find cg-bonding atom
                int found = 0;
                for(k=min; k<=max; k++) //loop over current residue
                {
                    if(strcmp(traj.atom_name[k].c_str(), param.param_main_s[j][1].c_str() ) == 0) //cg atom matches bonding cg-atom 
                    {
                        bond_cg_index = k; 
                        found = 1;
                        break;
                    }
                }
                if(found == 0)
                {
                    printf("CG connecting atom not found %s \n",param.param_main_s[j][1].c_str());
                    MPI_Finalize();
                    exit(EXIT_SUCCESS);
                }
                found = 0;

                //find aa-bonding atom 
                for(k=0; k<this_x.size(); k++) //loop over fragment
                {
                    if(strcmp(this_frag[k].c_str(), param.param_main_s[j][2].c_str() ) == 0) // frag-atom matches bonding atom   
	            {
                        bond_aa_index = k;
                        found = 1; 
			break;
	            }		    
                }
                if(found == 0)
                {
                    printf("Fragment connecting atom not found %s \n",param.param_main_s[j][2].c_str());
                    MPI_Finalize();
                    exit(EXIT_SUCCESS);
                }
                found = 0;

                //rotate fragment
                rotate_fragment(traj,s,p,param,this_x,this_y,this_z,bond_aa_index,bond_cg_index,this_frag,r_map);

		//add rotated fragment to rvec that will be written to pdb
                for(k=0; k<this_x.size(); k++) //loop over fragment atoms
                {
                    int current_index = count_res*p.num_atoms + count_frag + k;

                    r_aa[current_index][0] = this_x[k];
                    r_aa[current_index][1] = this_y[k];
                    r_aa[current_index][2] = this_z[k];
                } 

		count_frag = count_frag + this_x.size(); 
	    }
	}
    }

    int this_size = p.num_mols*p.num_atoms;

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
	this_atom_nr[j] = this_atom_nr[j]%99999;
    }

    //set the atom names
    int this_count = -1;
    for(i=0; i<p.num_mols; i++) //loop over molecules
    {
        for(j=0; j<param.main_size_y(); j++) //loop over fragments
        {
            for(k=0; k<param.sec_size_y(j); k++) //loop over fragment atoms
            {
                this_count++;
                this_atom_name[this_count] = param.param_sec_s[j][k][0].c_str();
            }
        }
    }

    //set the residue names
    for(j=0; j<this_size; j++) //loop over the atoms
    {
        this_res_name[j] = p.target;
    }

    //set residue numbers
    for(j=0; j<this_size; j++) //loop over protein atoms
    {
        this_res_nr[j] = floor(((double)j)/((double)p.num_atoms)) + 1;
        this_res_nr[j] = this_res_nr[j]%9999;
    }

    //open pdb file
    FILE *pdb_file;
    pdb_file             = fopen(p.frag_file_name.c_str(), "w");
    if(pdb_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",p.frag_file_name.c_str());
    }
    else
    {
        write_frame_pdb(traj.box,this_size,this_atom_nr,this_res_nr,this_res_name,this_atom_name,r_aa,traj.title,s.world_rank,&pdb_file,this_beta,this_weight,this_element,this_chain_id,i);
        fclose(pdb_file);
    }
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
    s.program_name = "Fragment Mapping";

    //force program to run in serial?
    enum Switch serial         = on;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                      s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                      s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                     s.world_rank, s.cl_tags, &p.b_print,   0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                                  s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                   s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                     s.world_rank, s.cl_tags, &p.b_lsq,     0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-frag",   p.frag_file_name,             "Output data file (gro, pdb) with all-atom molecules",             s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-target", p.target,                     "Residue type that is converted",                                  s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-crd",    p.param_file_name,            "Selection card with lipid atom info (crd)",                       s.world_rank, s.cl_tags, nullptr,      1);
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
    check_extension_mpi(s.world_rank,"-frag",p.frag_file_name,".pdb");

    //create parameter files
    Param param;

    //read parameter files
    param.get_param(p.param_file_name,4,3,4);

    //check the integrity of the parameter files
    if(param.check_file() == 0) //bad files. kill program 
    {
        MPI_Finalize();
        exit(EXIT_SUCCESS); 
    }

    p.num_mols  = count_molecules(traj,s,p,param);
    p.num_atoms = count_atoms(traj,s,p,param); 

    //extract mapping atoms from selection card
    int i = 0;                           //standard variable used in loops
    int j = 0;                           //standard variable used in loops
    sv2d mapping_atoms(0,sv1d(0));       //hold the mapping atoms for each fragment
    for(i=0; i<param.main_size_y(); i++) //loop over fragments
    {
       if(check_for_filename(param.param_main_s[i][0],".crd") == 1) //filename provided
       {
           Index mapping_atoms_index;
	   mapping_atoms_index.get_index(param.param_main_s[i][0].c_str());
        
           sv1d this_mapping_atoms(0);
           for(j=0; j<mapping_atoms_index.index_s.size(); j++) //loop over mapping atoms
           {
               this_mapping_atoms.push_back(mapping_atoms_index.index_s[j].c_str());
           }
           mapping_atoms.push_back(this_mapping_atoms);
       }
       else 
       {
           sv1d this_mapping_atoms(0);
           this_mapping_atoms.push_back(param.param_main_s[i][0].c_str());
           mapping_atoms.push_back(this_mapping_atoms);
       }  
    } 
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

        frag_mapping(traj,s,p,param,mapping_atoms);

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

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

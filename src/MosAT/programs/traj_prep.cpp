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
#include "MosAT/program_variables/pv_traj_prep.h"            //This has the variables specific to the analysis program
#include "headers/array.h"                                   //This has routines used for working with arrays
#include "headers/performance.h"                             //This has a class for logging performance data
#include "headers/index.h"                                   //This has a class for working with index files
#include "headers/traj.h"                                    //This has a class for working with the trajectory
#include "headers/leaflet_finder.h"                          //This has routines used to find leaflets in membrane simulations
#include "headers/protein_finder.h"                          //This has routines used to find protein atoms
#include "headers/sol_finder.h"                              //This has routines used to find the solvent
#include "headers/grid.h"                                    //This has routines used for working with a grid
#include "headers/protein.h"                                 //This has routines used for working with protein data
#include "headers/fit.h"                                     //This has routines used for fitting data
#include "headers/parallel.h"                                //This has routines for different parallelization schemes
#include "headers/force_serial.h"                            //This has routines used for forcing the code to run on a single mpi process
#include "headers/param.h"                                   //This has routines used for reading complex parameter data
#include "headers/grid_lt.h"                                 //This has routines used for reading in grid data
#include "headers/voronoi.h"                                 //This has routines used for computing voronoi diagrams
#include "headers/atom_select.h"                             //This has routines used for making atom selections using a selection text

/* 
Notes: Here, we have introduced the ability to remove some atoms from the trajectory. This feature works by 
       selecting the atoms one wishes to keep with atom_select. We thus read in an index file here with the 
       atoms to keep. The program then creates a second Trajectory object that works with the slimmed down 
       molecular system. The process is automated by extracting atomic coordinates from the reference file 
       for the target atoms. The new slimmed ref file is then provided to the second trajectory object as 
       both the -ref and -traj. By itself, this causes problems because there will be only a single frame. 
       The second traj will then assign no frames to all mpi cores other than rank 0. To fix this problem, 
       we create an ifo file for the new reference file. This info file is given the number of atoms for the 
       slimmed system, but the number of frames from the original trajectory. The position of each frame is 
       then set to zero. This way, we could loop over the number of frames specified, however, it would just 
       read the same ref file multiple times. Tricking the trajectory object into thinking there are multiple 
       frames creates another problem though. When there are multiple frames, the code will look for the time and 
       step in the title line of gro files. If this info is not found the program will hang. We thus add this 
       info to the title of the temporary ref file. The problem is further complicated by cross io. Here, we 
       have to copy many variables from the full sized trajectory like the coords, velocities, forces, box, 
       box dimensions, and stuff like chain id, etc.  
*/

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function generates an info file so the slimmed traj thinks there are multiple frames                 //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void gen_info_file(Trajectory &traj,system_variables &s,program_variables &p,int num_atoms,int num_frames,string info_file_name)
{
    if(s.world_rank == 0)
    {
        int i = 0;     //standard variable used in loop
        int64_t this_filesize = 0;
        int this_num_frames = 0;

        //read first 2 lines from existing info file
        FILE *info_file = fopen64(info_file_name.c_str(), "r");
        if(info_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",info_file_name.c_str());
        }
        char my_string[20];   //used to read in each item in the info file
        int line_count = 0;   //count the lines as they are read
        while(fscanf(info_file, "%s,", my_string) == 1)
        {
            if(line_count == 0) //file size
            {
                this_filesize = atol(my_string);
            }
            else if(line_count == 2) //frames
            {
                this_num_frames = atoi(my_string);
            }
            line_count++;
        }

        FILE *info_file_copy = fopen64("temp_ref.gro.info", "w");
        fprintf(info_file_copy," %ld \n", this_filesize);
        fprintf(info_file_copy," %d \n",num_atoms);
        fprintf(info_file_copy," %d \n",this_num_frames);
        for(i=0; i<this_num_frames; i++)
        {
            fprintf(info_file_copy," %ld \n",0);
        }
        fclose(info_file_copy);
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function creates a reference file with only the atoms specified in the index file                    //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_ref_slim(Trajectory &traj,system_variables &s,program_variables &p,Index &slim_index)
{
    if(s.world_rank == 0)
    {
        int i = 0;                                       //standard variable used in loops
        int num_atoms_slim = slim_index.index_i.size();  //how many atoms are there after removing other
        iv1d atom_nr_slim(num_atoms_slim,-1);            //hold the atom id's
        iv1d res_nr_slim(num_atoms_slim,-1);             //hold the res id's
        sv1d atom_name_slim(num_atoms_slim);             //hold the atom names
        sv1d res_name_slim(num_atoms_slim);              //hold the res names
        rvec *r_slim;                                    //store the coordinates
        rvec *v_slim;                                    //store the velocities (not used)

        //allocate memory for coordinates
        r_slim = (rvec *)calloc(num_atoms_slim , sizeof(*r_slim));
         
        for(i=0; i<slim_index.index_i.size(); i++) //loop index atoms
        {
            atom_nr_slim[i]   = traj.atom_nr[slim_index.index_i[i]-1];   
            res_nr_slim[i]    = traj.res_nr[slim_index.index_i[i]-1];   
            atom_name_slim[i] = traj.atom_name[slim_index.index_i[i]-1];   
            res_name_slim[i]  = traj.res_name[slim_index.index_i[i]-1];
            r_slim[i][0]      = traj.r_ref[slim_index.index_i[i]-1][0];
            r_slim[i][1]      = traj.r_ref[slim_index.index_i[i]-1][1];
            r_slim[i][2]      = traj.r_ref[slim_index.index_i[i]-1][2];
        }

        //make sure the title has time and step or the program will hang
        char this_title[200];
	string title_s = "step 0 time 0 \n";
        strcpy(this_title, title_s.c_str());

        //write out the new reference file
        string slim_ref_file_name = "temp_ref.gro";
        FILE *this_file = fopen64(slim_ref_file_name.c_str(), "w");
        write_frame_gro(traj.ibox,num_atoms_slim,atom_nr_slim,res_nr_slim,res_name_slim,atom_name_slim,r_slim,this_title,s.world_rank, &this_file,3,0,v_slim);
        fclose(this_file);

        //free memory
        free(r_slim);
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function creates a reference file with only the atoms specified in the index file                    //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void extract_atoms(Trajectory &traj,system_variables &s,program_variables &p,Index &slim_index,Trajectory &traj_slim)
{
    int i = 0;                                       //standard variable used in loops
    int j = 0;                                       //standard variable used in loops
    int counter = 0;                                 //track atoms as they are added

    //get the box
    for(i=0; i<3; i++)
    {
        for(j=0; j<3; j++)
        { 
            traj_slim.box[i][j] = traj.box[i][j];
        }
    }

    //get the title
    for(i=0; i<200; i++) 
    {
        traj_slim.title[i] = traj.title[i]; 
    }

    //get other variables
    traj_slim.step          = traj.step;
    traj_slim.time          = traj.time;
    traj_slim.bF            = traj.bF;
    traj_slim.bV            = traj.bV;
    traj_slim.bBox          = traj.bBox;
    traj_slim.global_frame  = traj.global_frame;             
    traj_slim.box_dimension = traj.box_dimension;  
 
    for(i=0; i<slim_index.index_i.size(); i++) //loop index atoms
    {
        //coords
        traj_slim.r[counter][0] = traj.r[slim_index.index_i[i]-1][0]; 
        traj_slim.r[counter][1] = traj.r[slim_index.index_i[i]-1][1];
        traj_slim.r[counter][2] = traj.r[slim_index.index_i[i]-1][2];

        //velocities
        if(traj.bV == 1)
        {
            traj_slim.v[counter][0] = traj.v[slim_index.index_i[i]-1][0];
            traj_slim.v[counter][1] = traj.v[slim_index.index_i[i]-1][1];
            traj_slim.v[counter][2] = traj.v[slim_index.index_i[i]-1][2];
        }
        //forces
        if(traj.bF == 1)
        {
            traj_slim.f[counter][0] = traj.f[slim_index.index_i[i]-1][0];
            traj_slim.f[counter][1] = traj.f[slim_index.index_i[i]-1][1];
            traj_slim.f[counter][2] = traj.f[slim_index.index_i[i]-1][2];
        }
        traj_slim.beta[counter]     = traj.beta[slim_index.index_i[i]-1];
        traj_slim.weight[counter]   = traj.weight[slim_index.index_i[i]-1];
        traj_slim.element[counter]  = traj.element[slim_index.index_i[i]-1];
        traj_slim.chain_id[counter] = traj.chain_id[slim_index.index_i[i]-1];
        counter++;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function checks the secondary selection cards for the correct extension                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void check_extension_recipe(string arg,string filename,string tag)
{
    int i             = 0;                  //standard variable used in loops
    int bad_extension = 0;                  //tells if the correct extension was provided
    int filename_size = filename.length();  //how long is the filename
    int tag_size      = tag.length();       //how long is the tag

    //make sure tag is not longer than filename
    if(tag_size > filename_size)
    {
        bad_extension = 1;
    }
    else //check for proper extension
    {
        for(i=filename_size-tag_size; i<filename_size; i++)
        {
            if(filename[i] != tag[i + tag_size - filename_size])
            {
                bad_extension = 1;
                break;
            }
        }
    }

    //report error and terminate program
    if(bad_extension == 1)
    {
        printf("The filename (%s) provided via the selection card (%s) requires the %s extension. \n",filename.c_str(),arg.c_str(),tag.c_str());
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function checks the arguments provided for a given operation                                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void check_args(Trajectory &traj,system_variables &s,program_variables &p,Param &param,sv1d &op_args,int i)
{
    int j = 0;    //standard variable used in loops
    int k = 0;    //standard variable used in loops

    for(j=0; j<param.sec_size_y(i); j++)
    {
        int found = 0;

        for(k=0; k<op_args.size(); k++) //loop over acceptable arguments
        {
            if(strcmp(op_args[k].c_str(), param.param_sec_s[i][j][0].c_str() ) == 0) //dimension
            {
                found = 1;
            }
        }

        //operation not found. terminate program
        if(found == 0)
        {
            if(s.world_rank == 0)
            {
                printf("Unrecognized argument (%s) identified in operation (%s, %s). Will terminate program. \n",param.param_sec_s[i][j][0].c_str(),param.param_main_s[i][0].c_str(),param.param_main_s[i][1].c_str());
            }
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function checks the recipe and reports any unrecognized arguments from the user                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void check_recipe(Trajectory &traj,system_variables &s,program_variables &p,Param &param)
{
    int i = 0;    //standard variable used in loops
    int j = 0;    //standard variable used in loops

    sv1d operations(0);  //holds supported operations

    operations.push_back("fit");     //least squares fitting
    operations.push_back("trans");   //translate system to a target point
    operations.push_back("center");  //center atom selection
    operations.push_back("wrap");    //put atoms/molecules in box
    operations.push_back("mend");    //fix broken molecules
    operations.push_back("time");    //set trajectory time
    operations.push_back("rot");     //rotate the system
    operations.push_back("shift");   //shift the system by a target amount

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check for unrecognized operation                                                                          //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(i=0; i<param.main_size_y(); i++) //loop over recipe
    {
        int found = 0;
        for(j=0; j<operations.size(); j++) //loop over operations
        {
            if(strcmp(operations[j].c_str(), param.param_main_s[i][0].c_str() ) == 0) //operation found
            {
                found = 1;
            } 
        }

        //operation not found. terminate program
        if(found == 0)
        {
            if(s.world_rank == 0)
            {
                printf("Unrecognized operation (%s) identified in recipe (%s). Will terminate program. \n",param.param_main_s[i][0].c_str(),p.param_file_name.c_str());
            }
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check secondary selection cards for the .crd extension                                                    //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(i=0; i<param.main_size_y(); i++) //loop over recipe
    {
        check_extension_recipe(p.param_file_name,param.param_main_s[i][1],".crd");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check each operation for unrecognized parameters                                                          //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(i=0; i<param.main_size_y(); i++) //loop over recipe
    {
        sv1d op_args(0);    //hold acceptable arguments for the operator

        if(strcmp("fit", param.param_main_s[i][0].c_str() ) == 0) //least squares fitting
        {
            op_args.push_back("lsq_d");    //argument for least squares fitting dimension
            op_args.push_back("lsq_r");    //argument for least squares fitting reference structure
            op_args.push_back("group");    //argument for least squares fitting atom selection

            check_args(traj,s,p,param,op_args,i);

        }
	else if(strcmp("trans", param.param_main_s[i][0].c_str() ) == 0) //translate system
        {
            op_args.push_back("x");        //target x-coordinate
            op_args.push_back("y");        //target y-coordinate
            op_args.push_back("z");        //target z-coordinate
            op_args.push_back("group");    //group of atoms whose center is placed at the target

            check_args(traj,s,p,param,op_args,i);
        }
        else if(strcmp("center", param.param_main_s[i][0].c_str() ) == 0) //center atoms selection
        {
            op_args.push_back("target");   //target (ref,frame_0,dynamic)
            op_args.push_back("group");    //group of atoms for centering            

            check_args(traj,s,p,param,op_args,i);
        }
        else if(strcmp("wrap", param.param_main_s[i][0].c_str() ) == 0) //put molecules in box
        {
            op_args.push_back("type");     //type (residues vs atom selection)
            op_args.push_back("group");    //group of atoms for wrapping             

            check_args(traj,s,p,param,op_args,i);
        }
        else if(strcmp("mend", param.param_main_s[i][0].c_str() ) == 0) //fix broken molecules
        {
            op_args.push_back("type");     //type (residues vs molecules from bonds)
            op_args.push_back("bonds");    //bonds list              
            op_args.push_back("cutoff");   //cutoff distance for the mend_mol_bonds option              

            check_args(traj,s,p,param,op_args,i);
        }
        else if(strcmp("time", param.param_main_s[i][0].c_str() ) == 0) //traj time
        {
            op_args.push_back("dt");       //time step between trajectory frames
            op_args.push_back("t0");       //the initial time
           
            check_args(traj,s,p,param,op_args,i); 
        }
        else if(strcmp("rot", param.param_main_s[i][0].c_str() ) == 0) //rotate system with matrix
        {
            op_args.push_back("matrix");   //filename with rotation matrix

            check_args(traj,s,p,param,op_args,i);
        }
        else if(strcmp("shift", param.param_main_s[i][0].c_str() ) == 0) //shift the system by a target amount
        {
            op_args.push_back("dx");   //distance to move in x-direction (nm)
            op_args.push_back("dy");   //distance to move in y-direction (nm)
            op_args.push_back("dz");   //distance to move in z-direction (nm)

            check_args(traj,s,p,param,op_args,i);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check the value of each operations parameters                                                             //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(i=0; i<param.main_size_y(); i++) //loop over recipe
    {
        if(strcmp("fit", param.param_main_s[i][0].c_str() ) == 0) //least squares fitting
        {
            for(j=0; j<param.sec_size_y(i); j++) //loop over current operations parameters
            {
                int kill = 0;

                if(strcmp("lsq_d", param.param_sec_s[i][j][0].c_str() ) == 0) //lsq_d
                {
                    if(param.param_sec_i[i][j][1] != 2 && param.param_sec_i[i][j][1] != 3 )
                    {
                        if(s.world_rank == 0)
                        {
                            printf("Acceptable options for (%s) provided in (%s) are 2 or 3. \n",param.param_sec_s[i][j][0].c_str(),param.param_main_s[i][1].c_str());
                        }
                        kill = 1;
                    }
                }
                if(strcmp("lsq_r", param.param_sec_s[i][j][0].c_str() ) == 0) //lsq_r
                {
                    if(param.param_sec_i[i][j][1] != 0 && param.param_sec_i[i][j][1] != 1 )
                    {
                        if(s.world_rank == 0)
                        {
                            printf("Acceptable options for lsq_r are 0 or 1. \n");
                        }
                        kill = 1;
                    }
                }
                if(strcmp("group", param.param_sec_s[i][j][0].c_str() ) == 0) //group
                {
                    check_extension_recipe(param.param_main_s[i][1],param.param_sec_s[i][j][1],".ndx");
                }
 
                if(kill == 1) //terminate program
                {
                    MPI_Finalize();
                    exit(EXIT_SUCCESS);
                }
            }
        }

        if(strcmp("trans", param.param_main_s[i][0].c_str() ) == 0) //translation
        {
            for(j=0; j<param.sec_size_y(i); j++) //loop over current operations parameters
            {
                int kill = 0;

                if(strcmp("x", param.param_sec_s[i][j][0].c_str() ) == 0) //target x
                {
                    if(check_float(param.param_sec_s[i][j][1].c_str()) == 0)
                    {
                        if(s.world_rank == 0)
                        {
                            printf("Argument (%s) provided in (%s) requires a floating point number while (%s) was provided. \n",param.param_sec_s[i][j][0].c_str(),param.param_main_s[i][1].c_str(),param.param_sec_s[i][j][1].c_str());
                        }
                        kill = 1;
                    }	
                }
                if(strcmp("y", param.param_sec_s[i][j][0].c_str() ) == 0) //target y
                {
                    if(check_float(param.param_sec_s[i][j][1].c_str()) == 0)
                    {
                        if(s.world_rank == 0)
                        {
                            printf("Argument (%s) provided in (%s) requires a floating point number while (%s) was provided. \n",param.param_sec_s[i][j][0].c_str(),param.param_main_s[i][1].c_str(),param.param_sec_s[i][j][1].c_str());
                        }
                        kill = 1;
                    }
                }
                if(strcmp("z", param.param_sec_s[i][j][0].c_str() ) == 0) //target z
                {
                    if(check_float(param.param_sec_s[i][j][1].c_str()) == 0)
                    {
                        if(s.world_rank == 0)
                        {
                            printf("Argument (%s) provided in (%s) requires a floating point number while (%s) was provided. \n",param.param_sec_s[i][j][0].c_str(),param.param_main_s[i][1].c_str(),param.param_sec_s[i][j][1].c_str());
                        }
                        kill = 1;
                    }
                }
                if(strcmp("group", param.param_sec_s[i][j][0].c_str() ) == 0) //group
                {
                    check_extension_recipe(param.param_main_s[i][1],param.param_sec_s[i][j][1],".ndx");
                }

		if(kill == 1) //terminate program
                {
                    MPI_Finalize();
                    exit(EXIT_SUCCESS);
                }
            }
        }

        if(strcmp("wrap", param.param_main_s[i][0].c_str() ) == 0) //wrap molecules
        {
            for(j=0; j<param.sec_size_y(i); j++) //loop over current operations parameters
            {
                int kill = 0;

                if(strcmp("type", param.param_sec_s[i][j][0].c_str() ) == 0) //type of wrapping
                {
                    if(param.param_sec_i[i][j][1] != 0 && param.param_sec_i[i][j][1] != 1 )
                    {
                        if(s.world_rank == 0)
                        {
                            printf("Acceptable options for (%s) provided in (%s) are 0 or 1. \n",param.param_sec_s[i][j][0].c_str(),param.param_main_s[i][1].c_str());
                        }
                        kill = 1;
                    }
                }
                if(strcmp("group", param.param_sec_s[i][j][0].c_str() ) == 0) //group
                {
                    check_extension_recipe(param.param_main_s[i][1],param.param_sec_s[i][j][1],".ndx");
                }

                if(kill == 1) //terminate program
                {
                    MPI_Finalize();
                    exit(EXIT_SUCCESS);
                }
            }
        }

        if(strcmp("center", param.param_main_s[i][0].c_str() ) == 0) //center selection
        {
            for(j=0; j<param.sec_size_y(i); j++) //loop over current operations parameters
            {
                int kill = 0;

                if(strcmp("target", param.param_sec_s[i][j][0].c_str() ) == 0) //target for center
                {
                    if(param.param_sec_i[i][j][1] != 0 && param.param_sec_i[i][j][1] != 1 && param.param_sec_i[i][j][1] != 2)
                    {
                        if(s.world_rank == 0)
                        {
                            printf("Acceptable options for (%s) provided in (%s) are 0, 1, or 2. \n",param.param_sec_s[i][j][0].c_str(),param.param_main_s[i][1].c_str());
                        }
                        kill = 1;
                    }
                }
                if(strcmp("group", param.param_sec_s[i][j][0].c_str() ) == 0) //group
                {
                    check_extension_recipe(param.param_main_s[i][1],param.param_sec_s[i][j][1],".ndx");
                }

                if(kill == 1) //terminate program
                {
                    MPI_Finalize();
                    exit(EXIT_SUCCESS);
                }
            }
        }

        if(strcmp("mend", param.param_main_s[i][0].c_str() ) == 0) //fix broken molecules
        {
            for(j=0; j<param.sec_size_y(i); j++) //loop over current operations parameters
            {
                int kill = 0;

                if(strcmp("type", param.param_sec_s[i][j][0].c_str() ) == 0) //type of mending
                {
                    if(param.param_sec_i[i][j][1] != 0 && param.param_sec_i[i][j][1] != 1 )
                    {
                        if(s.world_rank == 0)
                        {
                            printf("Acceptable options for (%s) provided in (%s) are 0 or 1. \n",param.param_sec_s[i][j][0].c_str(),param.param_main_s[i][1].c_str());
                        }
                        kill = 1;
                    }
                }
                if(strcmp("cutoff", param.param_sec_s[i][j][0].c_str() ) == 0) //cutoff distance for mend_mols_bonds option 
                {
                    if(check_float(param.param_sec_s[i][j][1].c_str()) == 0)
                    {
                        if(s.world_rank == 0)
                        {
                            printf("Argument (%s) provided in (%s) requires a floating point number while (%s) was provided. \n",param.param_sec_s[i][j][0].c_str(),param.param_main_s[i][1].c_str(),param.param_sec_s[i][j][1].c_str());
                        }
                        kill = 1;
                    }
                }
                if(strcmp("bonds", param.param_sec_s[i][j][0].c_str() ) == 0) //bonds
                {
                    check_extension_recipe(param.param_main_s[i][1],param.param_sec_s[i][j][1],".crd");
                }

                if(kill == 1) //terminate program
                {
                    MPI_Finalize();
                    exit(EXIT_SUCCESS);
                }
            }
        }

        if(strcmp("time", param.param_main_s[i][0].c_str() ) == 0) //time
        {
            for(j=0; j<param.sec_size_y(i); j++) //loop over current operations parameters
            {
                int kill = 0;

                if(strcmp("dt", param.param_sec_s[i][j][0].c_str() ) == 0) //dt
                {
                    if(check_float(param.param_sec_s[i][j][1].c_str()) == 0)
                    {
                        if(s.world_rank == 0)
                        {
                            printf("Argument (%s) provided in (%s) requires a floating point number while (%s) was provided. \n",param.param_sec_s[i][j][0].c_str(),param.param_main_s[i][1].c_str(),param.param_sec_s[i][j][1].c_str());
                        }
                        kill = 1;
                    }
                }
                if(strcmp("t0", param.param_sec_s[i][j][0].c_str() ) == 0) //initial time
                {
                    if(check_float(param.param_sec_s[i][j][1].c_str()) == 0)
                    {
                        if(s.world_rank == 0)
                        {
                            printf("Argument (%s) provided in (%s) requires a floating point number while (%s) was provided. \n",param.param_sec_s[i][j][0].c_str(),param.param_main_s[i][1].c_str(),param.param_sec_s[i][j][1].c_str());
                        }
                        kill = 1;
                    }
                }

                if(kill == 1) //terminate program
                {
                    MPI_Finalize();
                    exit(EXIT_SUCCESS);
                }
            }
        }
        if(strcmp("rot", param.param_main_s[i][0].c_str() ) == 0) //translation
        {
            for(j=0; j<param.sec_size_y(i); j++) //loop over current operations parameters
            {
                int kill = 0;

                if(strcmp("matrix", param.param_sec_s[i][j][0].c_str() ) == 0) //rotation matrix
                {
                    check_extension_recipe(param.param_main_s[i][1],param.param_sec_s[i][j][1],".crd");
                }
            }
        }
        if(strcmp("shift", param.param_main_s[i][0].c_str() ) == 0) //shift by a target amount
        {
            for(j=0; j<param.sec_size_y(i); j++) //loop over current operations parameters
            {
                int kill = 0;

                if(strcmp("dx", param.param_sec_s[i][j][0].c_str() ) == 0) //delta x
                {
                    if(check_float(param.param_sec_s[i][j][1].c_str()) == 0)
                    {
                        if(s.world_rank == 0)
                        {
                            printf("Argument (%s) provided in (%s) requires a floating point number while (%s) was provided. \n",param.param_sec_s[i][j][0].c_str(),param.param_main_s[i][1].c_str(),param.param_sec_s[i][j][1].c_str());
                        }
                        kill = 1;
                    }
                }
                if(strcmp("dy", param.param_sec_s[i][j][0].c_str() ) == 0) //delta y
                {
                    if(check_float(param.param_sec_s[i][j][1].c_str()) == 0)
                    {
                        if(s.world_rank == 0)
                        {
                            printf("Argument (%s) provided in (%s) requires a floating point number while (%s) was provided. \n",param.param_sec_s[i][j][0].c_str(),param.param_main_s[i][1].c_str(),param.param_sec_s[i][j][1].c_str());
                        }
                        kill = 1;
                    }
                }
                if(strcmp("dz", param.param_sec_s[i][j][0].c_str() ) == 0) //delta z
                {
                    if(check_float(param.param_sec_s[i][j][1].c_str()) == 0)
                    {
                        if(s.world_rank == 0)
                        {
                            printf("Argument (%s) provided in (%s) requires a floating point number while (%s) was provided. \n",param.param_sec_s[i][j][0].c_str(),param.param_main_s[i][1].c_str(),param.param_sec_s[i][j][1].c_str());
                        }
                        kill = 1;
                    }
                }

                if(kill == 1) //terminate program
                {
                    MPI_Finalize();
                    exit(EXIT_SUCCESS);
                }
            }
        }
    }
}

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
           printf("duplicate entry found in bonds list. Line %d %d %d \n",(int)((double)i/2.0),bond.index_i[i],bond.index_i[i+1]);
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
           printf("duplicate entry found in bonds list. Line %d %d %d \n",(int)((double)i/2.0),bond.index_i[i],bond.index_i[i+1]);
       }
    }

    //print info about the bonds. (who is bonded to each atom)
    if(p.b_test == 1)
    {	   
        string this_file_name = chop_and_add_tag(p.param_file_name,"_bonds.dat");
        FILE *this_file = fopen(this_file_name.c_str(), "w");
        if(this_file == NULL)
        {
            printf("failure opening %s for writing. \n",this_file_name.c_str());
        }
        else
        {
            for(i=0; i<bonds.size(); i++)
            {
                fprintf(this_file," %5s %5s %6d: ",traj.res_name[i].c_str(),traj.atom_name[i].c_str(),i+1);
                for(j=0; j<bonds[i].size(); j++)
                {
                    fprintf(this_file," %d ",bonds[i][j]);
                }
                fprintf(this_file,"\n");
            }
            fclose(this_file);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns a list of atoms bonded to a target atom                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
iv1d find_bonded_atoms(int target,iv2d &bonds,iv1d &status)
{
    int i = 0;  //standard variable used in loops
    int j = 0;  //standard variable used in loops
   
    /*                                                                                                                                                                  
     *   Notes:
     *   This function uses recursion and returns a list of atoms that are bonded to an atom. Since recursion is used, it is able to return a list 
     *   of atoms bond to the initial atom but also the atoms bonded to those atoms and so one. This is used to find all atoms that are conneced to 
     *   each other by one or more bonds thus defining a molecule. This function is initiated in get_molecules(). 
     */
 
    iv1d atoms(0,0); //stores the atom number

    atoms.push_back(target);
    status[target-1] = 1;

    for(i=0; i<bonds[target-1].size(); i++) //loop over bonded atoms
    {
        int this_atom_nr = bonds[target-1][i];

        if(status[this_atom_nr-1] == 0)
        {
            iv1d new_atoms = find_bonded_atoms(this_atom_nr,bonds,status);
            for(j=0; j<new_atoms.size(); j++)
            {
                atoms.push_back(new_atoms[j]);
            }
        }
    } 
    return atoms;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads in the bonding data and identifies molecules                                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_molecules(Trajectory &traj,system_variables &s,program_variables &p,iv2d &bonds,iv1d &status,iv2d &molecules)
{
    int i = 0;  //standard variable used in loops
    int j = 0;  //standard variable used in loops

    /*                                                                                                                                                                  
     *   Notes:
     *   Here, we hope to take the list of bonded atoms provided for each atom in the system (bonds[]) and sort the atoms into molecules. This is 
     *   accomplished by looping over the atoms in the system. We then check if an atom has already been assigned to a molecule (status[]). If not, then
     *   we pass the atom into the find_bonded_atoms() function. This function then adds the received atom to a list (atoms[]) and notes that the atom 
     *   has now been added to a molecule (status[]). The atoms bonded to the received atom (bonds[]) are then looped over. If these atoms have not been
     *   added to a molecule, then they are passed into find_bonded_atoms() as well. This routine therefore uses recursion. In this way, every new atom
     *   encountered (status[]) is examined and the atoms bonded to it are identified and added to the list (atoms[]). Finally, the find_bonded_atoms() 
     *   function, when no more new atoms are found will return the list of bonded atoms (atoms[]) which is added to the atoms from the previous itteration 
     *   of find_bonded_atoms() via the push_back() function. Finally, the first itteration of find_bonded_atoms() will return the full list of atoms making
     *   the molecule. This list (atoms[]) is then added to molecules[] using the push_back() function.    
     */

    for(i=0; i<traj.atoms(); i++)
    {
        if(status[i] == 0)
        {
            iv1d this_molecule = find_bonded_atoms(i+1,bonds,status);
            molecules.push_back(this_molecule);            
        }
    }
 
    //print molecule info
    if(p.b_test == 1)
    {
        string this_file_name = chop_and_add_tag(p.param_file_name,"_mol.pml");
        FILE *this_file = fopen(this_file_name.c_str(), "w");
        if(this_file == NULL)
        {
            printf("failure opening %s for writing. \n",this_file_name.c_str());
        }
        else
        {
            int stop = molecules.size(); //by default we print info for all molecules
            if(p.stop > 0)               //user specified an number of molecules to print info for
            {
                stop = p.stop; 
            }
            
            //print info to select molecules in pymol
            for(i=0; i<stop; i++) //loop over molecules
            {
                fprintf(this_file,"select mol_%d, ",i);
                for(j=0; j<molecules[i].size()-1; j++) //loop over current molecule
                {
                    fprintf(this_file,"(id %d and resi %d and resn %s) or ",molecules[i][j]%100000,traj.res_nr[molecules[i][j]-1]%10000,traj.res_name[molecules[i][j]-1].c_str());
                }
                fprintf(this_file,"(id %d and resi %d and resn %s) \n",molecules[i][molecules[i].size()-1]%100000,traj.res_nr[molecules[i][molecules[i].size()-1]-1]%10000,traj.res_name[molecules[i][molecules[i].size()-1]-1].c_str());
            }
            fclose(this_file);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads in molecules data and creates a list of bonds for each molecule                       //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_molecules_bonds(Trajectory &traj,system_variables &s,program_variables &p,iv2d &bonds,iv2d &molecules,iv3d &molecules_bonds)
{
    int i = 0;  //standard variable used in loops
    int j = 0;  //standard variable used in loops
    int k = 0;  //standard variable used in loops
    int l = 0;  //standard variable used in loops

    /*                                                                                                                                                                  
     *   Notes:
     *   This function is used to process a list of molecules data, where each molecule is defined by a list of atoms connected to each other by 
     *   one or more bonds, and extract a list of bonds for each molecule. The function works by looping over the molecules found in molecules[].
     *   Then, for each moleculei, the set of atoms composing it are looped over. For each atom in the molecule, the list of bonded atoms is looped
     *   over (bonds[]). This process gives a series of atom pairs (one atom from the bonds[] list and another from molecules[]). Then the bonds for 
     *   the molecule are looped over (molecule_bonds[]) and the atom pair is searched for. If the pair is not found it is added to molecule_bonds[].
     *   This process is continued until all the bonds are identified for the current molecule and for all molecules. 
     */

    //print info to select molecules in pymol
    for(i=0; i<molecules.size(); i++) //loop over molecules
    {
        for(j=0; j<molecules[i].size()-1; j++) //loop over current molecule
        {
            int this_atom_nr = molecules[i][j];

            for(k=0; k<bonds[this_atom_nr-1].size(); k++) //loop over bonds for current atom
            {
                int this_atom_nr_b = bonds[this_atom_nr-1][k];

                int found = 0;

                for(l=0; l<molecules_bonds[i].size(); l++) //loop over current molecule bonds
                {
                    if(molecules_bonds[i][l][0] == bonds[this_atom_nr-1][k] && molecules_bonds[i][l][1] == this_atom_nr)
                    {
                        found = 1;
                    }
                    else if(molecules_bonds[i][l][1] == bonds[this_atom_nr-1][k] && molecules_bonds[i][l][0] == this_atom_nr)
                    {
                        found = 1;
                    }
                }

                if(found == 0)
                {    
                    iv1d this_bond(2,0);
                    this_bond[0] = this_atom_nr;
                    this_bond[1] = this_atom_nr_b;
                    molecules_bonds[i].push_back(this_bond);
                }
            }
        }
    }
 
    //print PyMOL select commands for each bond in each molecule
    if(p.b_test == 1)
    {
        string this_file_name = chop_and_add_tag(p.param_file_name,"_mol_bonds.pml");
        FILE *this_file = fopen(this_file_name.c_str(), "w");
        if(this_file == NULL)
        {
            printf("failure opening %s for writing. \n",this_file_name.c_str());
        }
        else
        {
            int stop = molecules_bonds.size(); //by default we print info for all molecules
            if(p.stop > 0)                     //user specified an number of molecules to print info for
            {
                stop = p.stop;
            }

            for(i=0; i<stop; i++) //loop over molecules
            {
                fprintf(this_file,"mol %d \n",i);
                for(j=0; j<molecules_bonds[i].size(); j++) //loop over current molecule bonds 
                {
                    fprintf(this_file,"select atom_a, (id %d and resi %d and resn %s) \n",molecules_bonds[i][j][0]%100000,traj.res_nr[molecules_bonds[i][j][0]-1]%10000,traj.res_name[molecules_bonds[i][j][0]-1].c_str());
                    fprintf(this_file,"select atom_b, (id %d and resi %d and resn %s) \n",molecules_bonds[i][j][1]%100000,traj.res_nr[molecules_bonds[i][j][1]-1]%10000,traj.res_name[molecules_bonds[i][j][1]-1].c_str());
                    fprintf(this_file,"dist dist_%d, (atom_a), (atom_b) \n",j);
                }
            }
            fclose(this_file);
        }
    }	     
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function translates the atom selection center to a target position                                   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void trans(Trajectory &traj,system_variables &s,program_variables &p,Index &atom_group,double target_x,double target_y,double target_z)
{
    int i = 0;                                //standard variable used in loops

    dv1d center = traj.center_i(atom_group.index_i);

    double dx = target_x - center[0];
    double dy = target_y - center[1];
    double dz = target_z - center[2];

    for(i=0; i<traj.atoms(); i++) //loop over system atoms
    {
        traj.r[i][0] = traj.r[i][0] + dx;
        traj.r[i][1] = traj.r[i][1] + dy;
        traj.r[i][2] = traj.r[i][2] + dz;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function places the atom selection center in the box center                                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void center(Trajectory &traj,system_variables &s,program_variables &p,Index &atom_group,int target)
{
    int i = 0;                                //standard variable used in loops

    dv1d center = traj.center_i(atom_group.index_i);

    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;

    if(target == 0) //ref
    {
        dx = (0.5*traj.box_ref[XX][XX]) - center[0];
        dy = (0.5*traj.box_ref[YY][YY]) - center[1];
        dz = (0.5*traj.box_ref[ZZ][ZZ]) - center[2];
    }
    else if(target == 1) //frame 0
    {
        dx = (0.5*traj.ibox[XX][XX]) - center[0];
        dy = (0.5*traj.ibox[YY][YY]) - center[1];
        dz = (0.5*traj.ibox[ZZ][ZZ]) - center[2];
    }
    else if(target == 2) //dynamic
    {
        dx = (0.5*traj.box[XX][XX]) - center[0];
        dy = (0.5*traj.box[YY][YY]) - center[1];
        dz = (0.5*traj.box[ZZ][ZZ]) - center[2];
    }

    for(i=0; i<traj.atoms(); i++) //loop over system atoms
    {
        traj.r[i][0] = traj.r[i][0] + dx;
        traj.r[i][1] = traj.r[i][1] + dy;
        traj.r[i][2] = traj.r[i][2] + dz;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function places residues back inside the box                                                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void wrap(Trajectory &traj,system_variables &s,program_variables &p)
{
    int i = 0;                                //standard variable used in loops
    int j = 0;                                //standard variable used in loops
    int k = 0;                                //standard variable used in loops
    int l = 0;                                //standard variable used in loops
    int m = 0;                                //standard variable used in loops

    for(i=0; i<traj.atoms(); i++) //loop over system atoms
    {
        //prime loop for the next residue
        i = traj.next_residue(i);

        //get upper and lower range for current residue
        int min = traj.get_res_start(i);
        int max = traj.get_res_end(i);

        iv1d this_res(0,0);

        for(j=min; j<=max; j++) //loop over residue atoms
        {
            this_res.push_back(traj.atom_nr[j]);     
        }
        dv1d com = traj.com_i(this_res);

        //fix jumps in x direction
        if(com[0] < 0.0)
        {
            for(j=min; j<=max; j++) //loop over residue atoms
            {
                traj.r[j][0] = traj.r[j][0] + traj.box[XX][XX];
            }
        }
        else if(com[0] > traj.box[XX][XX])
        {
            for(j=min; j<=max; j++) //loop over residue atoms
            {
                traj.r[j][0] = traj.r[j][0] - traj.box[XX][XX];
            }
        }

        //fix jumps in y direction
        if(com[1] < 0.0)
        {
            for(j=min; j<=max; j++) //loop over residue atoms
            {
                traj.r[j][1] = traj.r[j][1] + traj.box[YY][YY];
            }
        }
        else if(com[1] > traj.box[YY][YY])
        {
            for(j=min; j<=max; j++) //loop over residue atoms
            {
                traj.r[j][1] = traj.r[j][1] - traj.box[YY][YY];
            }
        }

        //fix jumps in z direction
        if(com[2] < 0.0)
        {
            for(j=min; j<=max; j++) //loop over residue atoms
            {
                traj.r[j][2] = traj.r[j][2] + traj.box[ZZ][ZZ];
            }
        }
        else if(com[2] > traj.box[ZZ][ZZ])
        {
            for(j=min; j<=max; j++) //loop over residue atoms
            {
                traj.r[j][2] = traj.r[j][2] - traj.box[ZZ][ZZ];
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function places a specific selection back inside the box                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void wrap_this(Trajectory &traj,system_variables &s,program_variables &p,Index &atom_group,iv1d &atom_sel)
{
    int i = 0;                                //standard variable used in loops
    int j = 0;                                //standard variable used in loops
    int k = 0;                                //standard variable used in loops
    int l = 0;                                //standard variable used in loops
    int m = 0;                                //standard variable used in loops

    dv1d center = traj.center_i(atom_group.index_i);

    //fix jumps in x direction
    if(center[0] < 0.0)
    {
        for(j=0; j<traj.atoms(); j++) //loop over atoms
        {
            if(atom_sel[j] == 1)
            {
                traj.r[j][0] = traj.r[j][0] + traj.box[XX][XX];
            }
        }
    }
    else if(center[0] > traj.box[XX][XX])
    {
        for(j=0; j<traj.atoms(); j++) //loop over atoms
        {
            if(atom_sel[j] == 1)
            {
                traj.r[j][0] = traj.r[j][0] - traj.box[XX][XX];
            }
        }
    }

    //fix jumps in y direction
    if(center[1] < 0.0)
    {
        for(j=0; j<traj.atoms(); j++) //loop over atoms
        {
            if(atom_sel[j] == 1)
            {
                traj.r[j][1] = traj.r[j][1] + traj.box[YY][YY];
            }
        }
    }
    else if(center[1] > traj.box[YY][YY])
    {
        for(j=0; j<traj.atoms(); j++) //loop over atoms
        {
            if(atom_sel[j] == 1)
            {
                traj.r[j][1] = traj.r[j][1] - traj.box[YY][YY];
            }
        }
    }

    //fix jumps in z direction
    if(center[2] < 0.0)
    {
        for(j=0; j<traj.atoms(); j++) //loop over atoms
        {
            if(atom_sel[j] == 1)
            {
                traj.r[j][2] = traj.r[j][2] + traj.box[ZZ][ZZ];
            }
        }
    }
    else if(center[2] > traj.box[ZZ][ZZ])
    {
        for(j=0; j<traj.atoms(); j++) //loop over atoms
        {
            if(atom_sel[j] == 1)
            {
                traj.r[j][2] = traj.r[j][2] - traj.box[ZZ][ZZ];
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This mends broken residues.                                                                               //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void mend_res(Trajectory &traj,system_variables &s,program_variables &p)
{
    int i = 0;     //standard variable used in loops
    int j = 0;     //standard variable used in loops
    int k = 0;     //standard variable used in loops
    int l = 0;     //standard variable used in loops

    for(i=0; i<traj.atoms(); i++) //loop over system atoms
    {
        int min = traj.get_res_start(i);  //get first atom of residue
        int max = traj.get_res_end(i);    //get last atom of residue

        i = traj.next_residue(i);         //prime i for next residue

        int    target_index     = -1;
        double largest_sum_dist = 0.0;

        //find atom furthest from others in the residue
        for(j=min; j<=max; j++) //loop over residue atoms
        {
            double sum_dist = 0.0;

            for(k=min; k<=max; k++) //loop over residue atoms
            {
                double dx = traj.r[j][0] - traj.r[k][0];
                double dy = traj.r[j][1] - traj.r[k][1];
                double dz = traj.r[j][2] - traj.r[k][2];

                sum_dist = sum_dist + (dx*dx + dy*dy + dz*dz);
            }

	    if(sum_dist >= largest_sum_dist)
            {   
                target_index = j;
            }
        }

        //shift atoms if broken
        for(j=min; j<=max; j++) //loop over residue atoms
        {
            double dx = traj.r[j][0] - traj.r[target_index][0];
            double dy = traj.r[j][1] - traj.r[target_index][1];
            double dz = traj.r[j][2] - traj.r[target_index][2];

            //shift in x
            if(dx < -0.5*traj.box[XX][XX])
            {
                traj.r[j][0] = traj.r[j][0] + traj.box[XX][XX];
            }
            else if(dx > 0.5*traj.box[XX][XX])
            {
                traj.r[j][0] = traj.r[j][0] - traj.box[XX][XX];
            }

            //shift in y
            if(dy < -0.5*traj.box[YY][YY])
            {
                traj.r[j][1] = traj.r[j][1] + traj.box[YY][YY];
            }
            else if(dy > 0.5*traj.box[YY][YY])
            {
                traj.r[j][1] = traj.r[j][1] - traj.box[YY][YY];
            }

            //shift in z
            if(dz < -0.5*traj.box[ZZ][ZZ])
            {
                traj.r[j][2] = traj.r[j][2] + traj.box[ZZ][ZZ];
            }
            else if(dz > 0.5*traj.box[ZZ][ZZ])
            {
                traj.r[j][2] = traj.r[j][2] - traj.box[ZZ][ZZ];
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This rotates the molecular system                                                                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void rotate_system(Trajectory &traj,system_variables &s,program_variables &p,Index &rotation_matrix)
{
    int i = 0;    //standard variable used in loops
  
    int xx = rotation_matrix.index_d[0];

    float x_vec[3];            //holds the top row of rotation matrix
    float y_vec[3];            //holds the middle row of rotation matrix
    float z_vec[3];            //holds the bottom row of rotation matrix

    x_vec[0] = (float)rotation_matrix.index_d[0];
    x_vec[1] = (float)rotation_matrix.index_d[1];
    x_vec[2] = (float)rotation_matrix.index_d[2];

    y_vec[0] = (float)rotation_matrix.index_d[3];
    y_vec[1] = (float)rotation_matrix.index_d[4];
    y_vec[2] = (float)rotation_matrix.index_d[5];

    z_vec[0] = (float)rotation_matrix.index_d[6];
    z_vec[1] = (float)rotation_matrix.index_d[7];
    z_vec[2] = (float)rotation_matrix.index_d[8];

    for(i=0; i<traj.atoms(); i++) //loop over system atoms
    {
        traj.r[i][0] = iprod(x_vec, traj.r[i]);
        traj.r[i][1] = iprod(y_vec, traj.r[i]);
        traj.r[i][2] = iprod(z_vec, traj.r[i]);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This shifts the molecular system by a target amount                                                       //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void shift_system(Trajectory &traj,system_variables &s,program_variables,double dx,double dy,double dz)
{
    int i = 0;    //standard variable used in loops

    for(i=0; i<traj.atoms(); i++) //loop over system atoms
    {
        traj.r[i][0] = traj.r[i][0] + dx;
        traj.r[i][1] = traj.r[i][1] + dy;
        traj.r[i][2] = traj.r[i][2] + dz;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Forward declaration of update_distances() since this function calls shift_coords() and shift_coords()     //
// also calls update_distances()                                                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void update_distances(Trajectory &traj,system_variables &s,program_variables &p,iv3d &molecules_bonds,dv2d &distances,int num_bonds,int i,int shifted_atom,double cutoff); 

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This mends a broken molecule.                                                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void shift_coords(Trajectory &traj,system_variables &s,program_variables &p,iv3d &molecules_bonds,dv2d &distances,int num_bonds,int i,double cutoff,int atom_nr_a,int atom_nr_b,int j)
{
    /*                                                                                                                                                                  
     *   Notes:
     *   This function takes in a pair of atoms and check their distances in each direction for breaks. If a break is found, then the atom with the 
     *   larger coordinate is shifted. The shifted atom is then passed to update_distances() function, which updates the distance between this atom
     *   and any atoms bonded to it (molecules_bonds[]). Since update_distances() calls shift_coords() is any bonded atoms are found the process will 
     *   be continued until no atom pairs are found with a break.  
     */  

    //shift in x-direction
    if(distances[j][0] < -cutoff) //shift atom_b
    {
        int shifted_atom = atom_nr_b;
        traj.r[atom_nr_b-1][0] = traj.r[atom_nr_b-1][0] - traj.box[XX][XX];
        update_distances(traj,s,p,molecules_bonds,distances,num_bonds,i,shifted_atom,cutoff);
    }
    else if(distances[j][0] > cutoff) //shift atom_a
    {
        int shifted_atom = atom_nr_a;
        traj.r[atom_nr_a-1][0] = traj.r[atom_nr_a-1][0] - traj.box[XX][XX];
        update_distances(traj,s,p,molecules_bonds,distances,num_bonds,i,shifted_atom,cutoff);
    }

    //shift in y-direction
    if(distances[j][1] < -cutoff) //shift atom_b
    {
        int shifted_atom = atom_nr_b;
        traj.r[atom_nr_b-1][1] = traj.r[atom_nr_b-1][1] - traj.box[YY][YY];
        update_distances(traj,s,p,molecules_bonds,distances,num_bonds,i,shifted_atom,cutoff);
    }
    else if(distances[j][1] > cutoff) //shift atom_a
    {
        int shifted_atom = atom_nr_a;
        traj.r[atom_nr_a-1][1] = traj.r[atom_nr_a-1][1] - traj.box[YY][YY];
        update_distances(traj,s,p,molecules_bonds,distances,num_bonds,i,shifted_atom,cutoff);
    }

    //shift in z-direction
    if(distances[j][2] < -cutoff) //shift atom_b
    {
        int shifted_atom = atom_nr_b;
        traj.r[atom_nr_b-1][2] = traj.r[atom_nr_b-1][2] - traj.box[ZZ][ZZ];
        update_distances(traj,s,p,molecules_bonds,distances,num_bonds,i,shifted_atom,cutoff);
    }
    else if(distances[j][2] > cutoff) //shift atom_a
    {
        int shifted_atom = atom_nr_a;
        traj.r[atom_nr_a-1][2] = traj.r[atom_nr_a-1][2] - traj.box[ZZ][ZZ];
        update_distances(traj,s,p,molecules_bonds,distances,num_bonds,i,shifted_atom,cutoff);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This mends a broken molecule.                                                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void update_distances(Trajectory &traj,system_variables &s,program_variables &p,iv3d &molecules_bonds,dv2d &distances,int num_bonds,int i,int shifted_atom,double cutoff)
{
    int j = 0;  //standard variable used in loops
    int k = 0;  //standard variable used in loops

    /*                                                                                                                                                                  
     *   Notes:
     *   This function takes an atom id and updates the distance between this atom and any atom bonded to it (molecules_bonds[]). The function 
     *   Then passes the pair of bonded atoms to shift_coords() which checks the new distance for a break and shifts the atom with the larger 
     *   coordinate. Since shift_coords() calls update_distance() if a break is found, the process will continue to shift atoms until no breaks 
     *   are found. 
     */

    for(j=0; j<num_bonds; j++)
    {
        if(molecules_bonds[i][j][0] == shifted_atom)
        {
            int atom_nr_a = shifted_atom;
            int atom_nr_b = molecules_bonds[i][j][1];
 
            distances[j][0] = traj.r[atom_nr_a-1][0] - traj.r[atom_nr_b-1][0];
            distances[j][1] = traj.r[atom_nr_a-1][1] - traj.r[atom_nr_b-1][1];
            distances[j][2] = traj.r[atom_nr_a-1][2] - traj.r[atom_nr_b-1][2];

            shift_coords(traj,s,p,molecules_bonds,distances,num_bonds,i,cutoff,atom_nr_a,atom_nr_b,j);
        }
        else if(molecules_bonds[i][j][1] == shifted_atom)
        {
            int atom_nr_a = shifted_atom;
            int atom_nr_b = molecules_bonds[i][j][0];

            distances[j][0] = traj.r[atom_nr_a-1][0] - traj.r[atom_nr_b-1][0];
            distances[j][1] = traj.r[atom_nr_a-1][1] - traj.r[atom_nr_b-1][1];
            distances[j][2] = traj.r[atom_nr_a-1][2] - traj.r[atom_nr_b-1][2];

            shift_coords(traj,s,p,molecules_bonds,distances,num_bonds,i,cutoff,atom_nr_a,atom_nr_b,j);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This mends a broken molecule.                                                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void mend_mols_bonds(Trajectory &traj,system_variables &s,program_variables &p,iv3d &molecules_bonds,double cutoff)
{
    int i = 0;  //standard variable used in loops
    int j = 0;  //standard variable used in loops

    /*                                                                                                                                                                  
     *   Notes:
     *   This function uses a list of molecule bonds to check for broken molecules. This can be done since we know the bonds composing the molecule.
     *   We can thus look for instances where the bond length is too large (say greater than 1nm) and shift the atoms when such lengths are encountered.
     *   The algorithm starts by looping over the molecules in molecules_bonds[]. Then, for each molecule, the number of bonds is determined and these 
     *   bonds are looped over. The initial distances in each direction (distances[]) are measured for the bonds. These must be updated anytime an atom 
     *   is shifted. Then, the bonds are loop over again and each one sends a pair of atoms into the shift_coords() function. This function examines the  
     *   dx, dy, and dz terms between the pair and shifts the atom with the larger coordinate if a jump is detected. If so, the update_distances() function 
     *   is called. This function takes in the id of the shifted atom and updates the dx, dy, and dz terms for any atoms bonded to it (molecules_bonds[]). For 
     *   each bonded atom encountered, the shift_coords() is again called. This again takes in a pair of atoms. This time we have the atom that was shifted and an 
     *   atom bonded to it. Shift_coords() again checks the dx, dy, and dz terms for a break and if found the atom with the larger coordinate is shifted and passed
     *   into update_distances() where the process continues. The algorithm thus uses recursion. Eventually, no breaks will be encountered the program returns 
     *   to mend_mols_bonds() at which point the next molecule is examined. This outcome assumes that the broken molecules are a result of jumps across a 
     *   periodic boundary and that the molecules can be made whole again by shifting atoms by the box dimensions. If this is not true, the program will hang.       
     */

    for(i=0; i<molecules_bonds.size(); i++) //loop over molecules
    {
        int num_bonds = molecules_bonds[i].size();

        dv2d distances(num_bonds,dv1d(3,-1.0));

        //measure initial distances
        for(j=0; j<num_bonds; j++) //loop over current molecules bonds
        {
            int atom_nr_a = molecules_bonds[i][j][0];
            int atom_nr_b = molecules_bonds[i][j][1];

            distances[j][0] = traj.r[atom_nr_a-1][0] - traj.r[atom_nr_b-1][0];
            distances[j][1] = traj.r[atom_nr_a-1][1] - traj.r[atom_nr_b-1][1];
            distances[j][2] = traj.r[atom_nr_a-1][2] - traj.r[atom_nr_b-1][2];
        }

        //check for broken bond
        for(j=0; j<num_bonds; j++) //loop over current molecules bonds
        {
            int atom_nr_a = molecules_bonds[i][j][0];
            int atom_nr_b = molecules_bonds[i][j][1];

            shift_coords(traj,s,p,molecules_bonds,distances,num_bonds,i,cutoff,atom_nr_a,atom_nr_b,j);
        }
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
    s.program_name = "Traj Prep";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",               s.world_rank, s.cl_tags, &p.b_print,   1);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                 s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                            s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                             s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                               s.world_rank, s.cl_tags, &p.b_lsq,     0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                 s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)", s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-crd",    p.param_file_name,            "Selection card with list of operations (crd)",              s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-test",   &p.b_test,                    "Print info for testing molecule definitions (0:no, 1:yes)", s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-stop",   &p.stop,                      "Stop reporting -test info after this many molecules",       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-slim",   p.slim_index_file_name,       "Index for atoms to keep (ndx)",                             s.world_rank, s.cl_tags, &p.b_slim,    0);
    conclude_input_arguments_mpi(argc,argv,s.world_rank,s.program_name,s.cl_tags);

    //create a trajectory
    Trajectory traj; 

    //set trajectory parameters
    traj.set_block_parallel(on);
    traj.set_traj(p.in_file_name);
    traj.set_traj_w(p.out_file_name,p.b_print);
    traj.set_ref(p.ref_file_name);
    traj.set_lsq(p.lsq_index_file_name,p.b_lsq,p.lsq_dim,p.lsq_ref);
    traj.set_res(p.stride,p.start_frame,p.end_frame);

    //analyze the trajectory (log time spent) 
    perf.log_time(traj.build(),"Analyze Trajectory");

    //create a trajectory
    Trajectory traj_slim;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //check file extensions                                                                                     
    check_extension_mpi(s.world_rank,"-crd",p.param_file_name,".crd");

    //create parameter files
    Param param;

    //read parameter file
    param.get_param(p.param_file_name,2,1,2);

    //check the integrity of the parameter file
    if(param.check_file() == 0) //bad file. kill program 
    {
        MPI_Finalize();
        exit(EXIT_SUCCESS); 
    }

    //check recipe for unrecognizable parameters
    check_recipe(traj,s,p,param);

    int i = 0;  //standard variable used in loops
    int j = 0;  //standard variable used in loops
    int k = 0;  //standard variable used in loops

    //allocate memory for atom selections and bonds lists
    vector <Index> param_index(param.main_size_y());                //holds an index for selecting atoms each operation 
    iv2d param_sel(param.main_size_y(),iv1d(traj.atoms(),0));       //tags the selected atoms in param_index for fast identification of selected atoms 
    iv2d bonds(traj.atoms(),iv1d(0,0));                             //holds a list of bonded atoms for each atom in the system

    //read index files for least squres fitting operations
    for(i=0; i<param.main_size_y(); i++) //loop over operations
    {
        for(j=0; j<param.sec_size_y(i); j++) //loop over arguments for current operation
        {   
            if(strcmp("group", param.param_sec_s[i][j][0].c_str() ) == 0) //index file name
            {   
                param_index[i].get_index(param.param_sec_s[i][j][1].c_str()); //read the index file

                for(k=0; k<param_index[i].index_s.size(); k++) //loop over group of atoms
                {
                    param_sel[i][param_index[i].index_i[k]-1] = 1; //tag selected atom
                }
            }
            else if(strcmp("bonds", param.param_sec_s[i][j][0].c_str() ) == 0) //bonds list 
            {
                Index bond;
                bond.get_index(param.param_sec_s[i][j][1].c_str());
                get_bonds(traj,s,p,bond,bonds);
            }
            else if(strcmp("matrix", param.param_sec_s[i][j][0].c_str() ) == 0) //rotation matrix
            {
                 param_index[i].get_index(param.param_sec_s[i][j][1].c_str()); //read the index file

                 if(param_index[i].index_s.size() != 9)
                 {
                     if(s.world_rank == 0) 
                     {
                         printf("rotation matrix (%s) should contain 9 items \n",param.param_sec_s[i][j][1].c_str());
                     }
                     MPI_Finalize();
                     exit(EXIT_SUCCESS);
                 }
		 for(k=0; k<param_index[i].index_s.size(); k++) //check that matrix has floating point numbers only
                 {
                     if(check_float(param_index[i].index_s[k].c_str()) == 0)
                     {
                         if(s.world_rank == 0)
                         {
                             printf("rotation matrix (%s) should contain floating point numbers only \n",param.param_sec_s[i][j][1].c_str());
                         }
                         MPI_Finalize();
                         exit(EXIT_SUCCESS); 
                     }
                 } 
            }
        }
    }

    //estimate memory requirements to hold atom selection and bonds data
    double mem_sel = (double)param.main_size_y()*(double)traj.atoms()*4.0/1000000.0;
    if(s.world_rank == 0)
    {
        printf("Estimated memory usage: %f MB \n\n",mem_sel);
    }

    //get molecules from bonds list
    iv1d status(traj.atoms(),0);                                    //tells if an atom has been added to a molecule or not
    iv2d molecules(0,iv1d(0,0));                                    //stores a list of atoms for each molecule in the system
    get_molecules(traj,s,p,bonds,status,molecules);

    //get list of bonds for each molecule
    iv3d molecules_bonds(molecules.size(),iv2d(0,iv1d(2,0)));       //stores a list of bonds for each molecule in the system
    get_molecules_bonds(traj,s,p,bonds,molecules,molecules_bonds);

    //create index for removing atoms from trajectory
    Index slim_index;

    if(p.b_slim == 1)
    {
        //read index file
        slim_index.get_index(p.slim_index_file_name);

        //get reference file with atoms selected in slim_index
        get_ref_slim(traj,s,p,slim_index);

        //make an info file so it will think the ref_slim has multiple frames. otherwise, only rank 0 opens a file for writing when building. 
        string this_info_file_name = p.in_file_name + ".info";
        gen_info_file(traj,s,p,slim_index.index_i.size(),traj.get_ef_frames(),this_info_file_name);

        //set up new trajectory object
        traj_slim.set_block_parallel(on);
        traj_slim.set_traj("temp_ref.gro");
        traj_slim.set_traj_w(p.out_file_name,p.b_print);
	traj_slim.set_ref("temp_ref.gro");
        traj_slim.set_lsq(p.lsq_index_file_name,p.b_lsq,p.lsq_dim,p.lsq_ref);
        traj_slim.set_res(p.stride,p.start_frame,p.end_frame);

        //analyze the trajectory (log time spent) 
        perf.log_time(traj_slim.build(),"Analyze Traj Slim");
    }
    MPI_Barrier(MPI_COMM_WORLD);
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

	for(i=0; i<param.main_size_y(); i++) //loop over operations
        {
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                                                           //
            // Least Squares Fitting                                                                                     //
            //                                                                                                           //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            if(strcmp("fit", param.param_main_s[i][0].c_str() ) == 0) //perform least squares fitting
            {
                int this_lsq_d = 3;    //default value for lsq_d 
                int this_lsq_r = 0;    //default value for lsq_r

                for(j=0; j<param.sec_size_y(i); j++)
                {
                    if(strcmp("lsq_d", param.param_sec_s[i][j][0].c_str() ) == 0) //dimension
                    {
                        this_lsq_d = param.param_sec_i[i][j][1];
                    }
                    else if(strcmp("lsq_r", param.param_sec_s[i][j][0].c_str() ) == 0) //ref
                    {
                        this_lsq_r = param.param_sec_i[i][j][1];
                    }
                }
                
                traj.do_this_fit(param_index[i],this_lsq_d,this_lsq_r);
            }

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                                                           //
            // Translate selection to a target position                                                                  //
            //                                                                                                           //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            if(strcmp("trans", param.param_main_s[i][0].c_str() ) == 0) //translate system
            {
                double target_x = 0.0;    //default value for x    
                double target_y = 0.0;    //default value for y
                double target_z = 0.0;    //default value for z

                for(j=0; j<param.sec_size_y(i); j++)
                {
                    if(strcmp("x", param.param_sec_s[i][j][0].c_str() ) == 0) //target x
                    {
                        target_x = param.param_sec_d[i][j][1];
                    }
                    else if(strcmp("y", param.param_sec_s[i][j][0].c_str() ) == 0) //target y
                    {
                        target_y = param.param_sec_d[i][j][1];
                    }
                    else if(strcmp("z", param.param_sec_s[i][j][0].c_str() ) == 0) //target z
                    {
                        target_z = param.param_sec_d[i][j][1];
                    }
                }

                trans(traj,s,p,param_index[i],target_x,target_y,target_z);
            }

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                                                           //
            // Center a selection                                                                                        //
            //                                                                                                           //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            if(strcmp("center", param.param_main_s[i][0].c_str() ) == 0) //center selection
            {
                int target = 0;    //default value for target

                for(j=0; j<param.sec_size_y(i); j++)
                {
                    if(strcmp("target", param.param_sec_s[i][j][0].c_str() ) == 0) //target
                    {
                        target = param.param_sec_i[i][j][1];
                    }
                }

                center(traj,s,p,param_index[i],target);
            }

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                                                           //
            // Put residues or selection inside the box                                                                  //
            //                                                                                                           //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            if(strcmp("wrap", param.param_main_s[i][0].c_str() ) == 0) //put residues back inside the box
            {
                int type = 0;    //default value for type

                for(j=0; j<param.sec_size_y(i); j++)
                {
                    if(strcmp("type", param.param_sec_s[i][j][0].c_str() ) == 0) //type
                    {
                        type = param.param_sec_i[i][j][1];
                    }
                }

                if(type == 0) //residues
                {
                    wrap(traj,s,p);
                }
                else if(type == 1) //custom molecule
                {
                    wrap_this(traj,s,p,param_index[i],param_sel[i]);
                }
            }

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                                                           //
            // Fix broken molecules or residues                                                                          //
            //                                                                                                           //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            if(strcmp("mend", param.param_main_s[i][0].c_str() ) == 0) //fix broken molecules
            {
                int type      = 0;    //default value for type
                double cutoff = 1.0;  //cutoff distance for mend_mols_bonds option

                for(j=0; j<param.sec_size_y(i); j++)
                {
                    if(strcmp("type", param.param_sec_s[i][j][0].c_str() ) == 0) //type
                    {
                        type = param.param_sec_i[i][j][1];
                    }
                    else if(strcmp("cutoff", param.param_sec_s[i][j][0].c_str() ) == 0) //cutoff distance
                    {
                        cutoff = param.param_sec_d[i][j][1];
                    }
                }

                if(type == 0) //mend residues
                {
                    mend_res(traj,s,p);
                }
                else if(type == 1) //mend molecules
                {
                    mend_mols_bonds(traj,s,p,molecules_bonds,cutoff);
                }
            }

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                                                           //
            // Set the trajectory time                                                                                   //
            //                                                                                                           //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            if(strcmp("time", param.param_main_s[i][0].c_str() ) == 0) //set the trajectory time
            {
                double time_step = 1.0;    //default value for time_step     
                double t0        = 0.0;    //default value for t0 

                for(j=0; j<param.sec_size_y(i); j++)
                {
                    if(strcmp("dt", param.param_sec_s[i][j][0].c_str() ) == 0) //time step
                    {
                        time_step = param.param_sec_d[i][j][1];
                    }
		    if(strcmp("t0", param.param_sec_s[i][j][0].c_str() ) == 0) //initial time
                    {
                        t0 = param.param_sec_d[i][j][1];
                    }
                }
                traj.time = t0 + traj.get_frame_global()*time_step;
            }
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                                                           //
            // Rotate System                                                                                             //
            //                                                                                                           //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            if(strcmp("rot", param.param_main_s[i][0].c_str() ) == 0) //rotate the molecular system
            {
                for(j=0; j<param.sec_size_y(i); j++)
                {
                    if(strcmp("matrix", param.param_sec_s[i][j][0].c_str() ) == 0) //the rotation matrix
                    {
                        rotate_system(traj,s,p,param_index[i]);
                    }
                }
            }
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                                                           //
            // Shift System by a target amount                                                                           //
            //                                                                                                           //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            if(strcmp("shift", param.param_main_s[i][0].c_str() ) == 0) //shift the molecular system by a target amount
            {
                double dx = 0.0;        //shift this much in x direction
                double dy = 0.0;        //shift this much in y direction
                double dz = 0.0;        //shift this much in z direction

                for(j=0; j<param.sec_size_y(i); j++)
                {
                    if(strcmp("dx", param.param_sec_s[i][j][0].c_str() ) == 0) //dx
                    {   
                        dx = param.param_sec_d[i][j][1]; 
                    }
                    if(strcmp("dy", param.param_sec_s[i][j][0].c_str() ) == 0) //dy
                    {
                        dy = param.param_sec_d[i][j][1];
                    }
                    if(strcmp("dz", param.param_sec_s[i][j][0].c_str() ) == 0) //dz
                    {
                        dz = param.param_sec_d[i][j][1];
                    }
                }
                shift_system(traj,s,p,dx,dy,dz);
            }
        }

        traj.do_fit();

        if(p.b_slim == 0)
        {
            traj.write_traj_frame();
        }
        else if(p.b_slim == 1)
        {
            extract_atoms(traj,s,p,slim_index,traj_slim);
            traj_slim.write_traj_frame();
        }

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
        fflush(stdin);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //splice temporary traj file together (log time spent)
    if(p.b_slim == 0)
    {
        perf.log_time(traj.finalize_trajectory(),"Finalize Trajectory");
    }
    else if(p.b_slim == 1) 
    {
        perf.log_time(traj_slim.finalize_trajectory(),"Finalize Trajectory");
        remove("temp_ref.gro");
        remove("temp_ref.gro.info");
    }

    //print the performance stats
    perf.print_stats();

    //print closing statements
    print_closing(s.world_rank);

    //relinquish the mpi environment
    MPI_Finalize();

    return 0;
}

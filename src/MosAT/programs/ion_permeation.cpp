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
#include "MosAT/program_variables/pv_ion_permeation.h"       //This has the variables specific to the analysis program
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
// This function counts the number of target residues                                                        //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int count_targets(Trajectory &traj,system_variables &s,program_variables &p,iv1d &ion_res_nr)
{
    int       i = 0;     //standard variable used in loops
    int targets = 0;     //how many target residues are there

    for(i=0; i<traj.atoms(); i++) //loop over system atoms
    {
        //get the first and last atom of the current lipid
        int min = traj.get_res_start(i);
        int max = traj.get_res_end(i);

        //jump to the next lipid
        i = traj.next_residue(i);

        if(strcmp(traj.res_name[min].c_str(), p.target.c_str() ) == 0) //residue matches target type 
        {
            targets++;

            ion_res_nr.push_back(traj.res_nr[min]);
        }
    }

    //check if no targets were found
    if(targets == 0)
    {
        if(s.world_rank == 0)
        {
            printf("No residues found with type %s. Please provide a valid residue type with -target. \n",p.target.c_str());
        }
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }
    else
    {
        if(s.world_rank == 0)
        {
            printf("Target residues identified (%s): %d \n\n",p.target.c_str(),targets);
        }
    }

    return targets;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function uses distance measurements to make a bonds list                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double find_position_upper(Trajectory &traj,system_variables &s,program_variables &p,dv1d &this_center,Param &param,
                           int this_resi,string this_res_name)
{
    int     i = 0;     //standard variable used in loops
    int     j = 0;     //standard variable used in loops
    int     k = 0;     //standard variable used in loops

    dv1d top_dist(p.n_lipid,999999999.9);   //holds the top n closest lipids
    dv1d top_z(p.n_lipid,0.0);              //holds the z-coord of n closest lipids
    iv1d top_resi(p.n_lipid,0);             //holds the residue number of closest lipids
    sv1d top_resn(p.n_lipid);               //holds the residue name of closest lipids

    for(i=0; i<traj.target_leaflet.size(); i++) //loop over target leaflet
    {
        //jump to the next lipid
        i = traj.next_target_lipid(i);

        int min = traj.t_lip_start(i);         //first atom of the current lipid
        int max = traj.t_lip_end(i);           //last atom of the current lipid

        for(j=0; j<param.main_size_y(); j++) //loop over lipid types
        {
            if(strcmp(traj.res_name[min].c_str(), param.param_main_s[j][0].c_str() ) == 0) //lipid type is correct
            {   
                sv1d target_atoms = param.get_column_sec_s(j,0);
                dv1d lip_center   = traj.center(target_atoms,min,max);

		double dx = lip_center[0] - this_center[0];
                double dy = lip_center[1] - this_center[1];

                double dist = sqrt(dx*dx +dy*dy);
 
                double biggest_delta = 0.0;
                int    biggest_index = 0;
 
                dv1d delta(p.n_lipid,0.0);
                for(k=0; k<delta.size(); k++)
                {
                    delta[k] = top_dist[k] - dist;

                    if(delta[k] > biggest_delta)
                    {
                        biggest_delta = delta[k];
                        biggest_index = k;
                    }
                }

                if(dist < top_dist[biggest_index])
                {
                    top_dist[biggest_index] = dist;
                    top_z[biggest_index]    = lip_center[2]; 
                    top_resi[biggest_index] = traj.res_nr[min]; 
                    top_resn[biggest_index] = traj.res_name[min]; 
                } 
            }
        }
    }

    double avg = 0.0;
    for(i=0; i<top_z.size(); i++) //loop over z values
    {
        avg = avg + top_z[i];
    } 
    avg = avg/(double)top_z.size();

    return avg;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function uses distance measurements to make a bonds list                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double find_position_lower(Trajectory &traj,system_variables &s,program_variables &p,dv1d &this_center,Param &param,
                           int this_resi,string this_res_name)
{
    int     i = 0;     //standard variable used in loops
    int     j = 0;     //standard variable used in loops
    int     k = 0;     //standard variable used in loops

    dv1d top_dist(p.n_lipid,999999999.9);   //holds the top n closest lipids
    dv1d top_z(p.n_lipid,0.0);              //holds the z-coord of n closest lipids
    iv1d top_resi(p.n_lipid,0);             //holds the residue number of closest lipids
    sv1d top_resn(p.n_lipid);               //holds the residue name of closest lipids

    for(i=0; i<traj.opposing_leaflet.size(); i++) //loop over opposing leaflet
    {
        //jump to the next lipid
        i = traj.next_opposing_lipid(i);

        int min = traj.o_lip_start(i);         //first atom of the current lipid
        int max = traj.o_lip_end(i);           //last atom of the current lipid

        for(j=0; j<param.main_size_y(); j++) //loop over lipid types
        {
            if(strcmp(traj.res_name[min].c_str(), param.param_main_s[j][0].c_str() ) == 0) //lipid type is correct
            {
                sv1d target_atoms = param.get_column_sec_s(j,0);
                dv1d lip_center   = traj.center(target_atoms,min,max);

                double dx = lip_center[0] - this_center[0];
                double dy = lip_center[1] - this_center[1];

                double dist = sqrt(dx*dx +dy*dy);

                double biggest_delta = 0.0;
                int    biggest_index = 0;

                dv1d delta(p.n_lipid,0.0);
                for(k=0; k<delta.size(); k++)
                {
                    delta[k] = top_dist[k] - dist;

                    if(delta[k] > biggest_delta)
                    {
                        biggest_delta = delta[k];
                        biggest_index = k;
                    }
                }

                if(dist < top_dist[biggest_index])
                {
                    top_dist[biggest_index] = dist;
                    top_z[biggest_index]    = lip_center[2];
                    top_resi[biggest_index] = traj.res_nr[min];
                    top_resn[biggest_index] = traj.res_name[min];
                }
            }
        }
    }

    double avg = 0.0;
    for(i=0; i<top_z.size(); i++) //loop over z values
    {
        avg = avg + top_z[i];
    }
    avg = avg/(double)top_z.size();

    return avg;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function detemines the position of ions in the z-direction                                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void conductance(Trajectory &traj,system_variables &s,program_variables &p,dv2d &state,Param &param,dv2d &z_coord)
{
    int     i = 0;     //standard variable used in loops
    int     j = 0;     //standard variable used in loops
    int     k = 0;     //standard variable used in loops
    int count = -1;    //The position of residues in state

    for(i=0; i<traj.atoms(); i++) //loop over system atoms
    {
        //get the first and last atom of the current lipid
        int min = traj.get_res_start(i);
        int max = traj.get_res_end(i);

        //jump to the next lipid
        i = traj.next_residue(i);

        if(strcmp(traj.res_name[min].c_str(), p.target.c_str() ) == 0) //residue matches target type 
        {
            //get the position of the residue in state
            count++;

            //find the residue center
            sv1d target_atoms(0);
            for(j=min; j<=max; j++)
            {
                target_atoms.push_back(traj.atom_name[j]);
            }
            dv1d this_center = traj.center(target_atoms,min,max);

            //find local boundaries of membrane 
            double z_1 = find_position_upper(traj,s,p,this_center,param,traj.res_nr[min],traj.res_name[min]);
            double z_2 = find_position_lower(traj,s,p,this_center,param,traj.res_nr[min],traj.res_name[min]);

            //find position in z-direction
            if(this_center[2] >= z_1 + p.buf_up) //ion at the top buffer
            {
                state[traj.current_frame][count] = 0.0;
            }
	    else if(this_center[2] >= z_1 && this_center[2] < z_1 + p.buf_up) //ion at the top
            {
                state[traj.current_frame][count] = 1.0;
            }  
            else if(this_center[2] < z_1 && this_center[2] > z_2) //ion in the middle
            {
                state[traj.current_frame][count] = 2.0;
            }
            else if(this_center[2] <= z_2 && this_center[2] > z_2 - p.buf_low) //ion at the bottom
            {
                state[traj.current_frame][count] = 3.0;
            }
            else if(this_center[2] <= z_2 - p.buf_low) //ion at the bottom buffer
            {
                state[traj.current_frame][count] = 4.0;
            }

	    //store z-coords
	    z_coord[traj.current_frame][count] = this_center[2];
	}
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Collect states and write data to output file.                                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,dv2d &state,iv2d &history,
                         iv2d &history_frame,iv1d &ion_res_nr,dv2d &z_coord)
{
    int i = 0;              //standard variable used in loops
    int j = 0;              //standard variable used in loops 
    int k = 0;              //standard variable used in loops 
    int events_up = 0;      //How many upward events
    int events_down = 0;    //How many downward events

    iv2d perm_events(0,iv1d(6,0)); //stores info about each permeation event

    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nfinalizing analysis. This requires communication of the coordinates data and could take some time. \n");
    }

    //collect state and z-coord data
    collect_dv2d(s.world_size,s.world_rank,state);
    collect_dv2d(s.world_size,s.world_rank,z_coord);

    if(s.world_rank == 0)
    {
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                                 //
        // Substrate moving across membrane from top to bottom                                                             //
        //                                                                                                                 //
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        for(i=0; i<traj.get_ef_frames(); i++) //loop over trajectory
        {
            for(j=0; j<p.num_targets; j++) //loop over target ions
            {
                if(state[i][j] == 0.0) //at the top buffer
                {
                    history      [j][0] = 1;
                    history      [j][1] = 0;
                    history      [j][2] = 0;
                    history      [j][3] = 0;
                    history      [j][4] = 0;
		    history_frame[j][0] = i;
                    history_frame[j][1] = i;
                    history_frame[j][2] = i;
                    history_frame[j][3] = i;
                    history_frame[j][4] = i;
                }
                else if(state[i][j] == 1.0) //at the top
                {
                    if(history[j][0] == 1)
                    {
                        history      [j][1] = 1;
                        history_frame[j][1] = i;
                    }
                    else
                    {
                        history      [j][1] = 0;
                        history_frame[j][1] = i;
                    }
                }
                else if(state[i][j] == 2.0) //in the middle
                {
                    if(history[j][1] == 1)
                    {
                        history      [j][2] = 1;
			history_frame[j][2] = i;
                    }
		    else  
                    {
                        history      [j][2] = 0;
			history_frame[j][2] = i;
                    } 
                }
                else if(state[i][j] == 3.0) //at the bottom
                {
                    if(history[j][2] == 1)
                    {
                        history      [j][3] = 1;
                        history_frame[j][3] = i;
                    }
                    else
                    {
                        history      [j][3] = 0;
                        history_frame[j][3] = i;
                    }
                }
                else if(state[i][j] == 4.0) //at the bottom buffer
                {   
                    if(history[j][0] == 1 && history[j][1] == 1 && history[j][2] == 1 && history[j][3] == 1) //traversed all previous 4 segments 
                    {   
                        events_down++;

                        //store info about the event
			iv1d this_event(6,0);
			this_event[0] = j;
                        this_event[1] = history_frame[j][0]; 
                        this_event[2] = history_frame[j][1];
                        this_event[3] = history_frame[j][2];
                        this_event[4] = history_frame[j][3];
                        this_event[5] = history_frame[j][4];
                        perm_events.push_back(this_event);
                    }
                    history      [j][0] = 0;
                    history      [j][1] = 0;
                    history      [j][2] = 0;
                    history      [j][3] = 0;
                    history      [j][4] = 1;
		    history_frame[j][0] = i;
                    history_frame[j][1] = i;
                    history_frame[j][2] = i; 
                    history_frame[j][3] = i; 
                    history_frame[j][4] = i; 
                }
            } 
        }
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                                 //
        // Substrate moving across membrane from bottom to top                                                             //
        //                                                                                                                 //
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        for(i=0; i<traj.get_ef_frames(); i++) //loop over trajectory
        {
            for(j=0; j<p.num_targets; j++) //loop over target ions
            {
                if(state[i][j] == 0.0) //at the top buffer
                {
                    if(history[j][4] == 1 && history[j][3] == 1 && history[j][2] == 1 && history[j][1] == 1) //traversed all previous 2 segments 
                    {
                        events_up++;

                        //store info about the event
                        iv1d this_event(6,0);
                        this_event[0] = j;
                        this_event[1] = history_frame[j][0];
                        this_event[2] = history_frame[j][1];
                        this_event[3] = history_frame[j][2];
                        this_event[4] = history_frame[j][3];
                        this_event[5] = history_frame[j][4];
                        perm_events.push_back(this_event);
                    }
                    history      [j][0] = 1;
                    history      [j][1] = 0;
                    history      [j][2] = 0;
                    history      [j][3] = 0;
                    history      [j][4] = 0;

		    history_frame[j][0] = i;
                    history_frame[j][1] = i;
                    history_frame[j][2] = i;
                    history_frame[j][3] = i;
                    history_frame[j][4] = i;
                }

                else if(state[i][j] == 1.0) //at the top
                {
                    if(history[j][2] == 1)
                    {
                        history      [j][1] = 1;
                        history_frame[j][1] = i;
                    }
                    else
                    {
                        history      [j][1] = 0;
                        history_frame[j][1] = i;
                    }
                }
                else if(state[i][j] == 2.0) //in the middle
                {
                    if(history[j][3] == 1)
                    {
                        history      [j][2] = 1;
			history_frame[j][2] = i;
                    }
                    else 
                    {
                        history      [j][2] = 0;
			history_frame[j][2] = i;
                    }
                }
                else if(state[i][j] == 3.0) //at the bottom
                {
                    if(history[j][4] == 1)
                    {
                        history      [j][3] = 1;
                        history_frame[j][3] = i;
                    }
                    else
                    {
                        history      [j][3] = 0;
                        history_frame[j][3] = i;
                    }
                }
                else if(state[i][j] == 4.0) //at the bottom buffer
                {
                    history      [j][0] = 0;
                    history      [j][1] = 0;
                    history      [j][2] = 0;
                    history      [j][3] = 0;
                    history      [j][4] = 1;
		    history_frame[j][0] = i;
                    history_frame[j][1] = i;
                    history_frame[j][2] = i;
                    history_frame[j][3] = i;
                    history_frame[j][4] = i;
                }
            }
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                                 //
        // Print info about permeation events                                                                              //
        //                                                                                                                 //
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	printf("\n");
        printf("Upward events: %d \n",events_up);
        printf("Downward events: %d \n",events_down);
        printf("\n");

	//print info about permeation events
        printf("downward events:\n");
        printf("%10s %10s %10s %10s %12s %12s %13s %-10s \n","event","ion","duration","frame_top","frame_middle","frame_bottom","frame_bot_buf","pymol");
        printf("%10s-%10s-%10s-%10s-%12s-%12s-%13s-%10s---------- \n","----------","----------","----------","----------","------------","------------","-------------","----------");
        for(i=0; i<events_down; i++) //loop over permeation events
        {
            int duration = perm_events[i][4] - perm_events[i][1];
            printf("%10d %10d %10d %10d %12d %12d %13d select ion_%d, resi %d and resn %s \n",i,perm_events[i][0],duration,perm_events[i][1],perm_events[i][2],perm_events[i][3],perm_events[i][4],i,ion_res_nr[perm_events[i][0]]%10000,p.target.c_str());
        }
        printf("\nupward events:\n");
        printf("%10s %10s %10s %10s %12s %12s %13s %-10s \n","event","ion","duration","frame_top","frame_middle","frame_bottom","frame_bot_buf","pymol");
        printf("%10s-%10s-%10s-%10s-%12s-%12s-%13s-%10s---------- \n","----------","----------","----------","----------","------------","------------","-------------","----------");
	for(i=events_down; i<perm_events.size(); i++) //loop over permeation events
        {
            int duration = perm_events[i][4] - perm_events[i][1];
            printf("%10d %10d %10d %10d %12d %12d %13d select ion_%d, resi %d and resn %s \n",i,perm_events[i][0],duration,perm_events[i][1],perm_events[i][2],perm_events[i][3],perm_events[i][4],i,ion_res_nr[perm_events[i][0]]%10000,p.target.c_str());
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                                 //
        // Print z-coord data                                                                                              //
        //                                                                                                                 //
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        FILE *z_file;
        z_file = fopen(p.conductance_file_name.c_str(), "w");
        if(z_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",p.conductance_file_name.c_str());
        }
        else
        {
            fprintf(z_file," %10s ","#frame");
            for(i=0; i<p.num_targets; i++)
            {
                string tag = "#" + to_string(i);
                fprintf(z_file," %10s ",tag.c_str());
            }
            fprintf(z_file,"\n");	

	    for(i=0; i<traj.get_ef_frames(); i++) //loop over trajectory frames
            {
                fprintf(z_file," %10d ",i*p.stride);
                for(j=0; j<p.num_targets; j++) //loop over ions
                {
                    fprintf(z_file," %10f ",z_coord[i][j]);
                }
                fprintf(z_file,"\n");
            }
            fclose(z_file);
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                                 //
        // Print z-coord data for permeating ions only                                                                     //
        //                                                                                                                 //
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        string this_file_name = chop_and_add_tag(p.conductance_file_name,"_perm.dat");
 
        z_file = fopen(this_file_name.c_str(), "w");
        if(z_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",this_file_name.c_str());
        }
        else
        {
            fprintf(z_file," %10s ","#frame");
            for(i=0; i<p.num_targets; i++)
            {
                for(j=0; j<perm_events.size(); j++) //loop over permeation events
                {
                    if(perm_events[j][0] == i)
                    {
                        string tag = "#" + to_string(i);
                        fprintf(z_file," %10s ",tag.c_str());
                    }
                }
            }
            fprintf(z_file,"\n");

            for(i=0; i<traj.get_ef_frames(); i++) //loop over trajectory frames
            {
                fprintf(z_file," %10d ",i*p.stride);
                for(j=0; j<p.num_targets; j++) //loop over ions
                {
                    for(k=0; k<perm_events.size(); k++) //loop over permeation events
                    {
                        if(perm_events[k][0] == j)
                        {
                            fprintf(z_file," %10f ",z_coord[i][j]);
                        }
                    }
                }
                fprintf(z_file,"\n");
            }
            fclose(z_file);
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
    s.program_name = "Ion Permeation";

    //force program to run in serial?
    enum Switch serial         = off;

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
    add_argument_mpi_s(argc,argv,"-cond",   p.conductance_file_name,      "Output data file (dat) with conductance data",                    s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-target", p.target,                     "Residue type that will be monitored",                             s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,           "PDB file with sorted leaflets (pdb)",                             s.world_rank, s.cl_tags, &p.b_lf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,  "File with additional leaflet finder parameters (prm) ",           s.world_rank, s.cl_tags, &p.b_lf_param,0);
    add_argument_mpi_s(argc,argv,"-crd",    p.param_file_name,            "Selection card with lipid atom info (crd)",                       s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-buf_u",  &p.buf_up,                    "Distance for upper buffer (nm)",                                  s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-buf_l",  &p.buf_low,                   "Distance for lower buffer (nm)",                                  s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-n",      &p.n_lipid,                   "How many nearest lipids when computing average z-coord?",         s.world_rank, s.cl_tags, nullptr,      0);
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
    check_extension_mpi(s.world_rank,"-cond",p.conductance_file_name,".dat");
    if(p.b_lf_pdb == 1)
    {
        check_extension_mpi(s.world_rank,"-lf_pdb",p.lf_pdb_file_name,".pdb");
    }
    if(p.b_lf_param == 1)
    {
        check_extension_mpi(s.world_rank,"-lf_prm",p.leaflet_finder_param_name,".prm");
    }

    //create parameter files
    Param param;

    //read parameter files
    param.get_param(p.param_file_name,2,1,1);

    //check the integrity of the parameter files
    if(param.check_file() == 0) //bad files. kill program 
    {
        MPI_Finalize();
        exit(EXIT_SUCCESS); 
    }

    //set target leaflet to upper
    p.leaflet = 1; 

    //run leaflet finder
    traj.get_leaflets(p.leaflet,p.leaflet_finder_param_name,p.b_lf_param);

    //print a pdb with distinguished leaflets
    traj.write_leaflets(p.lf_pdb_file_name,p.b_lf_pdb);

    //print info about the target leaflet
    traj.get_leaflet_stats();

    iv1d ion_res_nr(0,0);                                            //store the residue number of each ion

    //count target residues
    p.num_targets = count_targets(traj,s,p,ion_res_nr);

    //report memory estimates
    int mem_ion_res_nr    = p.num_targets*4; 
    int mem_state         = traj.get_num_frames()*p.num_targets*8;
    int mem_history       = p.num_targets*5*4; 
    int mem_history_frame = p.num_targets*5*4;
    int mem_z_coord       = traj.get_ef_frames()*p.num_targets*8; 
    int mem_total         = mem_ion_res_nr + mem_state + mem_history + mem_history_frame + mem_z_coord; 
   
    if(s.world_rank == 0)
    {
        printf("Estimated memory: %f (MB) \n\n",(double)mem_total/1000000.0);
    }

    //create structures to store data
    dv2d state(traj.get_num_frames(),dv1d(p.num_targets,0.0));       //hold the state of each ion for each trajectory frame. 
    iv2d history(p.num_targets, iv1d(5,0));                          //the state history of the ion       
    iv2d history_frame(p.num_targets, iv1d(5,0));                    //store the frame when each event occurs      
    dv2d z_coord(traj.get_num_frames(),dv1d(p.num_targets,0.0));     //hold the z-coord of each ion for each trajectory frame. 

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

        conductance(traj,s,p,state,param,z_coord);

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //splice temporary traj file together (log time spent)
    perf.log_time(traj.finalize_trajectory(),"Finalize Trajectory");

    //collect data and cound permeation events
    perf.log_time(finalize_analysis(traj,s,p,state,history,history_frame,ion_res_nr,z_coord),"Finalize Analysis");

    //print the performance stats
    perf.print_stats();

    //print closing statements
    print_closing(s.world_rank);

    //relinquish the mpi environment
    MPI_Finalize();

    return 0;
}

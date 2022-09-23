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
#include "headers/mosat_routines.h"                         //This is where most of the functions called in main are located
#include "headers/file_naming.h"                             //This has routines for added tags to an existing file name    
#include "headers/file_naming_mpi.h"                         //This has routines for added tags to an existing file name (mpi)
#include "headers/command_line_args_mpi.h"                   //This has routines for adding command line arguments
#include "MosAT/program_variables/pv_lipid_mixing.h"        //This has the variables specific to the analysis program
#include "headers/array.h"                                   //This has routines used for working with arrays
#include "headers/performance.h"                             //This has a class for logging performance data
#include "headers/index.h"                                   //This has a class for working with index files
#include "headers/traj.h"                                    //This has a class for working with the trajectory
#include "headers/leaflet_finder.h"                          //This has routines used to find leaflets in membrane simulations
#include "headers/protein_finder.h"                          //This has routines used to find protein atoms
#include "headers/sol_finder.h"                              //This has routines used to find the solvent
#include "headers/grid.h"                                    //This has routines used for working with a grid
#include "headers/protein.h"                                 //This has routines used for working with protein data
#include "headers/parallel.h"                                //This has routines for different parallelization schemes
#include "headers/force_serial.h"                            //This has routines used for forcing the code to run on a single mpi process
#include "headers/param.h"                                   //This has routines used for reading complex parameter data

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Check if the arguments were provided to perform a distance projection                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void check_proj(Trajectory &traj,system_variables &s,program_variables &p)
{
    if(p.b_m1 == 1 || p.b_m2 == 1 || p.b_APS == 1 || p.b_r == 1 || p.b_cutoff == 1) //one of the arguments was set
    {
        if(p.b_m1 == 0 || p.b_m2 == 0 || p.b_APS == 0 || p.b_r == 0) //one of the arguments was not set 
        {
            if(s.world_rank == 0)
            {
                printf("To project the solvation number onto the XY plane you must provide the -m1, -m2, -APS, and -r tags. \n");
            }
            //terminate program
            MPI_Finalize();
            exit(EXIT_SUCCESS);           
        }
        else //do the distance projection
        {
            p.b_dist_proj = 1;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints to file binding events as they are encountered.                                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void report_binding_events(Trajectory &traj,system_variables &s,program_variables &p,int time,int lip_nr_1,int res_nr_1,int lip_count,string res_name_lip,int res_nr_lip)
{
    if(p.b_report_binding_events == 1)
    {
        int bind_i    = traj.current_frame - time;      //first frame binding was noted 
        int bind_f    = traj.current_frame - 1;         //last time binding was noted
        int num_lines = 0;                              //num lines in binding events file

        //create file name 
        string tag = "_" + to_string(lip_nr_1) + ".be";
        string binding_events_file_name = chop_and_add_tag(p.mix_file_name,tag);

        //check if header information has been printed
        FILE *binding_events_file = fopen(binding_events_file_name.c_str(), "r");
        if(binding_events_file == NULL)
        {
            //file does not yet exits. print header
        }
        else //file exists but could be empty. count lines
        {
            int c = 0;
            do
            {
              c = fgetc(binding_events_file);
              if (c == '\n')
              {
                  num_lines++;
                  break;       //header must exist. exit loop
              }
            }
            while (c != EOF);

            //close file
            fclose(binding_events_file);
        }

        //open file for writing (append)
        binding_events_file = fopen(binding_events_file_name.c_str(), "a");
        if(binding_events_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",binding_events_file_name.c_str());
        }
        else //file exists or was created succesfully
        {
            //print header information
            if(num_lines == 0)
            {
                fprintf(binding_events_file," lip_nr %10d res_nr %10d ef_dt(ps) %10f ef_frames %d num_lipids_2 %10d num_lip_1 %10d junk %10d APS(nm^2) %10f \n\n",lip_nr_1,res_nr_1,p.ef_dt,traj.get_ef_frames(),p.num_lipids_2,p.num_lipids_1,-1,p.APS);
                fprintf(binding_events_file," %10s %10s %10s %15s %15s %20s \n","lipid","res_nr","res_name","bind_i(frame)","bind_f(frame)","dwell time(frames)");
                fprintf(binding_events_file," %10s-%10s-%10s-%15s-%15s-%20s \n","----------","----------","----------","---------------","---------------","--------------------");
            }

            //report binding events
            fprintf(binding_events_file," %10d %10d %10s %15d %15d %20d \n",lip_count,res_nr_lip,res_name_lip.c_str(),bind_i,bind_f,time);

            //close file
            fclose(binding_events_file);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function finds first shell neighbors for the lipids and keeps track of ther precent of lipids to     //
// lipids to have visited this sehll                                                                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void update_neighbors(Trajectory &traj,system_variables &s,program_variables &p,Param &param_1,Param &param_2,vector <vector <double>> &neighbors,
                      vector <double> &percent_visited,vector <vector <int>> &bound,int dwell_times[],Grid &nbrs,iv3d &contacts)

{
    //Notes: Here the strategy is to loop over the target leaflet atoms and look for lipid types of the correct type. Once a target lipid is
    //found the mpi process checks if the lipid is one of its assigned lipid. If so, The geometric center is calculated. Then, the target leaflet
    //is again looped over and target lipids of the second type are searched for. Once a lipid is found of the correct type the same geometric center
    //is found for the lipid and the distance between the 2 centers is computed. If the 2 centers are close enough then the second lipid is counted as 
    //being in the first lipids first shell (neighbors[][] is set to 1). As the second loop progresses we also count how many lipids are in the first 
    //lipids first shell (num_neighbors++). And if the lipids are found to not be close enough but were previously then the dwell time is computedi and 
    //recorded (dwell_times[]++). Once the second loop has completed the number of neighbors is stamped to the grid. At this point the focus becomes computing
    //the percentage of lipid_2 to have been in lipid_1's first shell. Because we want rank 0 to have this inormation for all lipids we do collect this data. 
    //At the same time we collect the number of neighbor each lipid_1 has for the current frame. Then rank zero computes the average percent_visited over lipid_1
    //and the average current number of neighbors over lipid_1. The spread in both of these quantities over lipid_1 is computed next. And finally, the data is written
    //to the ouptu file.     

    int    i           = 0;                                  //standard variable used in loops
    int    j           = 0;                                  //standard variable used in loops
    int    k           = 0;                                  //standard variable used in loops
    int    l           = 0;                                  //standard variable used in loops
    int    m           = 0;                                  //standard variable used in loops
    int    n           = 0;                                  //standard variable used in loops
    int    o           = 0;                                  //standard variable used in loops
    int    prev_lip_1  = -1;                                 //last lipid (1) encountered
    int    prev_lip_2  = -1;                                 //last lipid (2) encountered
    int    lip_count_1 = -1;                                 //count lipids (1) as they are encountered 
    int    lip_count_2 = -1;                                 //count lipids (2) as they are encountered
    int    pos         = traj.current_frame%(2*p.range + 1); //where to store data in contacts[][][]
    int res_nr_lip[p.num_lipids_2];                          //Store the res_nr for reporting binding events
    string res_name_lip[p.num_lipids_2];                     //Store the res_name for reporting binding events
    double num_neighbors[p.num_lipids_1];                    //how many first shell neighbors in the current frame        
    init_darray(num_neighbors,p.num_lipids_1);

    if(p.b_dist_proj == 1)
    {
        //clear the current frame grids
        nbrs.clean_frame();
    }

    for(i=0; i<traj.target_leaflet.size(); i++) //loop over target leaflet 
    {
        //get the first and last atom of the current lipid
        int min = traj.t_lip_start(i);
        int max = traj.t_lip_end(i);

        //jump to the next lipid
        i = traj.next_target_lipid(i);
  
        for(j=0; j<param_1.main_size_y(); j++) //loop over lipid types 1
        {
            if(strcmp(traj.res_name[min].c_str(), param_1.param_main_s[j][0].c_str()) == 0) //lipid type is correct
            {
                //check for a new lipid
                if(traj.res_nr[traj.target_leaflet[i]-1] != prev_lip_1)
                {
                   lip_count_1++;
                   prev_lip_1 = traj.res_nr[traj.target_leaflet[i]-1];
                }

                if(lip_count_1 >= traj.lipid_start && lip_count_1 <= traj.lipid_end) //each core does their contributing lipids
                {
                    //extract atoms involved in geometric center and compute the center
                    sv1d target_atoms_1 = param_1.get_column_sec_s(j,0);

                    //compute the center for the lipid
                    dv1d r_center_1 = traj.center(target_atoms_1,min,max);

                    //reset the lipid 2 index
                    lip_count_2 = -1;
                    prev_lip_2  = -1;

                    for(m=0; m<traj.target_leaflet.size(); m++) //loop over target leaflet
                    {
                        //get the first and last atom of the current lipid
                        int min_2 = traj.t_lip_start(m);
                        int max_2 = traj.t_lip_end(m);

                        //jump to the next lipid
                        m = traj.next_target_lipid(m);

                        if(traj.res_nr[traj.target_leaflet[m]-1] != traj.res_nr[traj.target_leaflet[i]-1]) //dont count lipid i as being in its own first shell
                        {
                            for(n=0; n<param_2.main_size_y(); n++) //loop over lipid types 2 
                            {
                                if(strcmp(traj.res_name[min_2].c_str(), param_2.param_main_s[n][0].c_str()) == 0) //lipid type is correct
                                {
                                    //check for a new lipid
                                    if(traj.res_nr[traj.target_leaflet[m]-1] != prev_lip_2)
                                    {
                                       lip_count_2++;
                                       prev_lip_2                = traj.res_nr[traj.target_leaflet[m]-1];
                                       res_nr_lip[lip_count_2]   = traj.res_nr[min_2];
                                       res_name_lip[lip_count_2] = traj.res_name[min_2];
                                    }

                                    //extract atoms involved in geometric center and compute the center
                                    sv1d target_atoms_2 = param_2.get_column_sec_s(j,0);

                                    //compute the center for the lipid
                                    dv1d r_center_2 = traj.center(target_atoms_2,min_2,max_2);

                                    //compute distance between two lipid midpoints
                                    double dx = r_center_1[0] - r_center_2[0];
                                    double dy = r_center_1[1] - r_center_2[1];  

                                    double dist = sqrt(dx*dx + dy*dy);

                                    double cutoff_dist = param_1.param_main_d[n][2];

                                    //update neighbors data and compute dwell time
                                    if(dist < cutoff_dist) //lipid is in the first shell
                                    {
                                        contacts[lip_count_1][lip_count_2][pos] = 1;
                                    }
                                    else //lipid is not in the first shell
                                    {
                                        contacts[lip_count_1][lip_count_2][pos] = 0;          
                                    }
                                }
                            }
                        }
                    }

                    //check that window is full. if so count neighbors etc. 
                    if(traj.current_frame >= (2*p.range)) //window is full
                    {
                        for(m=0; m<p.num_lipids_2; m++) //loop over lipids_2
                        {
                            int count = 0;

                            for(n=0; n<(2*p.range + 1); n++) //loop over window
                            { 
                                if(contacts[lip_count_1][m][n] == 1)  
                                {
                                    count++;
                                }        
                            }
        
                            double freq = count/(double)(2*p.range + 1);

                            if(freq >= p.window_cutoff)
                            {
                                //update neighbors data and compute dwell time
                                neighbors[lip_count_1][m]  = 1.0;
                                num_neighbors[lip_count_1] = num_neighbors[lip_count_1] + 1.0;
                                bound[lip_count_1][m]      = bound[lip_count_1][m] + 1;
                            }
                            else //lipid is not in the first shell
                            {
                                if(bound[lip_count_1][m] > 0) //lipid was in the first shell. add dwell time
                                {
                                   int time = bound[lip_count_1][m];
                                   dwell_times[time-1] = dwell_times[time-1] + 1;

                                   //report binding event for checking
                                   report_binding_events(traj,s,p,time,lip_count_1,traj.res_nr[min],m,res_name_lip[m],res_nr_lip[m]);
                                }
                                bound[lip_count_1][m] = 0;
                            }
                        }
 
                        if(p.b_dist_proj == 1)
                        {
                            //loop over current residue atoms of lip_1 and add neighbors count to grid
                            for(m=min; m<=max; m++)
                            {
                                if(strcmp(traj.atom_name[m].c_str(), p.map_1.c_str()) == 0) //atom is mapping atom 1
                                {
                                    double rx = traj.r[m][0];    //mapping atom x-coord
                                    double ry = traj.r[m][1];    //mapping atom y-coord

                                    nbrs.stamp(rx,ry,p.radius,num_neighbors[lip_count_1]);
                                }
                                if(strcmp(traj.atom_name[m].c_str(), p.map_2.c_str()) == 0) //atom is mapping atom 2
                                {
                                    double rx = traj.r[m][0];    //mapping atom x-coord
                                    double ry = traj.r[m][1];    //mapping atom y-coord

                                    nbrs.stamp(rx,ry,p.radius,num_neighbors[lip_count_1]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    //need to create a new file here. if window size > 1 then the file is not 
    //opened on current_frame == 0. This way a new file is alway created when running program.
    FILE *mix_file;
    if(s.world_rank == 0 && traj.current_frame == 0)
    {
        mix_file = fopen(p.mix_file_name.c_str(), "w");
        if(mix_file == NULL)
        {
            printf("failure opening %s. \n",p.mix_file_name.c_str());
        }
        else
        {
            fprintf(mix_file," %15s %15s %15s %15s %15s \n","#frame","#mix_frac_avg","#mix_frac_stdev","#sol_num_avg","#sol_num_stdev");
            fclose(mix_file);
        }
    }

    if(traj.current_frame >= (2*p.range)) //window is full
    {
        if(p.b_dist_proj == 1)
        {
            //add the current frame grid to long term sum
            nbrs.add_frame_direct();
        }

        //compute percent visited and write data to output file
        if(traj.current_frame%p.mix_stride == 0) //write to file percent visited data on this frame
        {
            //now update percent visited
            for(i=0; i<p.num_lipids_1; i++) //loop over lipids 1
            {
                double visited       = 0;
                double my_neighbors  = num_neighbors[i];

                for(j=0; j<p.num_lipids_2; j++) //loop over lipids 2
                {
                    visited = visited + neighbors[i][j];
                }

                double world_visited[s.world_size];
                double world_neighbors[s.world_size];

                MPI_Gather(&visited,      1, MPI_DOUBLE, world_visited,   1, MPI_DOUBLE, 0,MPI_COMM_WORLD);
                MPI_Gather(&my_neighbors, 1, MPI_DOUBLE, world_neighbors, 1, MPI_DOUBLE, 0,MPI_COMM_WORLD);

                if(s.world_rank == 0)
                {
                    visited      = 0;
                    my_neighbors = 0;

                    for(j=0; j<s.world_size; j++) //loop over world 
                    {
                        visited      = visited      + world_visited[j];
                        my_neighbors = my_neighbors + world_neighbors[j];
                    }

                    percent_visited[i] = 100.0*visited/(double)p.num_lipids_2;
                    num_neighbors[i]   = my_neighbors;
                }
            }

            if(s.world_rank == 0)
            {
                //compute percent visited averaged over lipid 1
                double percent_visited_avg = 0;
                double num_neighbors_avg   = 0;

                for(i=0; i<p.num_lipids_1; i++) //loop over lipids 1
                {
                    percent_visited_avg = percent_visited_avg + percent_visited[i];
                    num_neighbors_avg   = num_neighbors_avg   + num_neighbors[i];
                }
                percent_visited_avg = percent_visited_avg/(double)p.num_lipids_1;
                num_neighbors_avg   = num_neighbors_avg/(double)p.num_lipids_1;

                //compute standard deviation of percent visited over lipid 1
                double percent_visited_stdev = 0;
                double num_neighbors_stdev   = 0;

                for(i=0; i<p.num_lipids_1; i++) //loop over lipids 1
                {
                    percent_visited_stdev = percent_visited_stdev + pow((percent_visited[i] - percent_visited_avg),2);
                    num_neighbors_stdev   = num_neighbors_stdev   + pow((num_neighbors[i] - num_neighbors_avg),2);
                }
                percent_visited_stdev = percent_visited_stdev/(double)(p.num_lipids_1-1);
                percent_visited_stdev = sqrt(percent_visited_stdev);
                num_neighbors_stdev   = num_neighbors_stdev/(double)(p.num_lipids_1-1);
                num_neighbors_stdev   = sqrt(num_neighbors_stdev);

                //write data to output file
                mix_file = fopen(p.mix_file_name.c_str(), "a");
                if(mix_file == NULL)
                {
                    printf("failure opening %s. \n",p.mix_file_name.c_str());
                }
                else
                {
                    fprintf(mix_file," %15d %15f %15f %15f %15f \n",traj.current_frame,percent_visited_avg,percent_visited_stdev,num_neighbors_avg,num_neighbors_stdev);
                    fclose(mix_file);
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects dwell_times from each mpi process and computes the average dwell time              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analyisis(Trajectory &traj,system_variables &s,program_variables &p,vector <vector <double>> &neighbors,
                      vector <double> &percent_visited,vector <vector <int>> &bound,int dwell_times[],Grid &nbrs)
{
    int i = 0;
    int j = 0;
    int world_count[s.world_size];
    init_iarray(world_count,s.world_size);

    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nfinalizing analysis. This requires communication of dwell times and could take some time. \n");
    }

    if(p.b_dist_proj == 1)
    {
        //collect nbrs and rho from all ranks
        nbrs.collect_grid();

        //normalize gnbrs and write the <gnbrs> and rho to file
        if(s.world_rank == 0)
        {
            nbrs.normalize();

            nbrs.exclude_data(p.cutoff,1);

            nbrs.write_grid();
            nbrs.write_rho();
        }
    }

    //collect the dwell_times from mpi processes
    for(i=0; i<traj.get_ef_frames(); i++) //loop over dwell times
    {
        int count = dwell_times[i];
        MPI_Gather(&count,      1, MPI_INT, world_count,   1, MPI_INT, 0,MPI_COMM_WORLD);

        if(s.world_rank == 0)
        {
            count = 0;
            for(j=0; j<s.world_size; j++) //loop over world
            {
                count = count + world_count[j];
            }
            dwell_times[i] = count;
        }
    }

    //compute the average dwell time and write to file the dwell time distribution
    if(s.world_rank == 0)
    {
        string freq_file_name = add_tag(p.mix_file_name,"_freq");
        FILE *freq_file = fopen(freq_file_name.c_str(), "w");
        if(freq_file == NULL)
        {

        }
        else
        {
            //count total number of binding events
            int binding_events = 0;
            for(i=0; i<traj.get_ef_frames(); i++) //loop over dwell times
            {
                binding_events = binding_events + dwell_times[i];
            }

            //compute average dwell time
            double freq = 0;
            double avg_dwell_time = 0;
            fprintf(freq_file,"%20s %15s \n","#dwell time (ps)","#probability");
            for(i=0; i<traj.get_ef_frames(); i++) //loop over dwell times
            {
                freq = (double)dwell_times[i]/(double)binding_events;
                avg_dwell_time = avg_dwell_time + freq*(i+1)*p.ef_dt;

                fprintf(freq_file,"%20f %15f \n",(i+1)*p.ef_dt,freq);
            }
            printf("Average dwell time %f ps \n",avg_dwell_time);
            fclose(freq_file);
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
    s.program_name = "Lipid Mixing";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                           s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                           s.world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                            s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                                       s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                        s.world_rank, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                          s.world_rank, &p.b_lsq,     0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                            s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",            s.world_rank, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-crd_1",  p.param_1_file_name,          "Selection card 1 (crd)",                                               s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-crd_2",  p.param_2_file_name,          "Selection card 2 (crd)",                                               s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-mix",    p.mix_file_name,              "Output data file with mixing fraction and solvation number (dat)",     s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,           "PDB file with sorted leaflets (pdb)",                                  s.world_rank, &p.b_lf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,  "File with additional leaflet finder parameters (prm)",                 s.world_rank, &p.b_lf_param,0);
    add_argument_mpi_i(argc,argv,"-leaf",   &p.leaflet,                   "Which leaflet? (0:both 1:upper 2:lower)",                              s.world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-mix_s",  &p.mix_stride,                "How often to print mixing data to output file (steps)",                s.world_rank, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-dt",     &p.ef_dt,                     "Effective time step between anayzed frames (counting stide, ps) ",     s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-m1",     p.map_1,                      "Name of mapping atom 1 ",                                              s.world_rank, &p.b_m1,      0);
    add_argument_mpi_s(argc,argv,"-m2",     p.map_2,                      "Name of mapping atom 2 ",                                              s.world_rank, &p.b_m2,      0);
    add_argument_mpi_d(argc,argv,"-APS",    &p.APS,                       "Area per grid box (nm^2)",                                             s.world_rank, &p.b_APS,     0);
    add_argument_mpi_d(argc,argv,"-r",      &p.radius,                    "Radius of the target atom (nm)",                                       s.world_rank, &p.b_r,       0);
    add_argument_mpi_d(argc,argv,"-cutoff", &p.cutoff,                    "Cutoff for excluding data (chi)",                                      s.world_rank, &p.b_cutoff,  0);
    add_argument_mpi_d(argc,argv,"-bx",     &p.box_x,                     "Grid x dimension (nm)",                                                s.world_rank, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-by",     &p.box_y,                     "Grid y dimension (nm)",                                                s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-odf",    &p.out_data_format,           "What is the output data format? (0:matrix 1:vector)",                  s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-range",  &p.range,                     "Noise filter half width (frames)",                                     s.world_rank, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-freq",   &p.window_cutoff,             "Noise filter significance threshold (omega)",                          s.world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-w_bind", &p.b_report_binding_events,   "Report binding events for each lipid's solvation shell? (0:no 1:yes)", s.world_rank, nullptr,      0);
    conclude_input_arguments_mpi(argc,argv,s.world_rank,s.program_name);

    //create a trajectory
    Trajectory traj; 

    //set trajectory parameters
    traj.set_block_parallel(off);
    traj.set_traj(p.in_file_name);
    traj.set_ref(p.ref_file_name);
    traj.set_lsq(p.lsq_index_file_name,p.b_lsq,p.lsq_dim,p.lsq_ref);
    traj.set_res(p.stride,p.start_frame,p.end_frame);

    //analyze the trajectory (log time spent) 
    perf.log_time(traj.build(),"Analyze Trajectory");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //check file extensions                                                                                     
    check_extension_mpi(s.world_rank,"-crd_1",p.param_1_file_name,".crd");
    check_extension_mpi(s.world_rank,"-crd_2",p.param_2_file_name,".crd");
    check_extension_mpi(s.world_rank,"-mix",p.mix_file_name,".dat");

    if(p.b_lf_pdb == 1)
    {
        check_extension_mpi(s.world_rank,"-lf_pdb",p.lf_pdb_file_name,".pdb");
    }
    if(p.b_lf_param == 1)
    {
        check_extension_mpi(s.world_rank,"-lf_prm",p.leaflet_finder_param_name,".prm");
    }

    //create parameter files
    Param param_1;
    Param param_2;

    //read parameter files
    param_1.get_param(p.param_1_file_name,3,1,1);
    param_2.get_param(p.param_2_file_name,2,1,1);

    //check the integrity of the parameter files
    if(param_1.check_file() == 0 || param_2.check_file() == 0) //bad files. kill program 
    {
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }

    //run leaflet finder
    traj.get_leaflets(p.leaflet,p.leaflet_finder_param_name,p.b_lf_param);

    //print a pdb with distinguished leaflets
    traj.write_leaflets(p.lf_pdb_file_name,p.b_lf_pdb);

    //print info about the target leaflet
    traj.get_leaflet_stats();

    //print info about the lipid selection
    traj.get_lipid_selection_stats(param_1.get_column_s(0),"-crd_1");
    traj.get_lipid_selection_stats(param_2.get_column_s(0),"-crd_2");

    //count the number of lipids in the target leaflet(s) of the target type
    p.num_lipids_1 = traj.count_target_lipids_type(param_1.get_column_s(0));
    p.num_lipids_2 = traj.count_target_lipids_type(param_2.get_column_s(0));

    //create arrays for holding data
    iv3d contacts(p.num_lipids_1, iv2d(p.num_lipids_2,iv1d(2*p.range + 1,0))); //keep record of which lipid were in contact with each other
    dv2d neighbors(p.num_lipids_1, dv1d(p.num_lipids_2,0.0));                  //a record of whol (lip_2) has been a in lip_1 first shell
    dv1d percent_visited(p.num_lipids_1,0.0);                                  //the percentace of lip_2 to have been in lip_1's first shell
    iv2d bound(p.num_lipids_1, iv1d(p.num_lipids_2,0.0));                      //A record of how many frames a lip_2 has been in lip_1's first shell
    int dwell_times[traj.get_ef_frames()];                                     //the number of lipids having a dwell time in the first shell
    init_iarray(dwell_times,traj.get_ef_frames());

    //parallelize by lipids
    traj.parallelize_by_lipid(p.num_lipids_1);

    //check if the solvation number is projected onto XY
    check_proj(traj,s,p);

    //create a grid to hold the number of neighbors
    Grid nbrs;

    if(p.b_dist_proj == 1)
    {
        //get the grid dimensions
        nbrs.get_dim(p.box_x,p.box_y,traj.ibox,p.APS);

        //set the output file name for grid
        string nbrs_file_name = add_tag(p.mix_file_name,"_neighbors");
        nbrs.set_output(nbrs_file_name,p.out_data_format);

        //print info about the grid
        nbrs.print_dim();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //print info about the worlk load distribution
    traj.workload_lipid();

    //print that analysis is beginning
    traj.report_progress();

    s.t = clock();
    //read read frames of the trajector and perform analysis
    for(traj.current_frame=0; traj.current_frame<traj.get_num_frames(); traj.current_frame++)
    {
        traj.read_traj_frame();

        traj.do_fit();

        update_neighbors(traj,s,p,param_1,param_2,neighbors,percent_visited,bound,dwell_times,nbrs,contacts);

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect num neighbors and write grid to file. Also collect dwell times and write distribution to file. 
    perf.log_time(finalize_analyisis(traj,s,p,neighbors,percent_visited,bound,dwell_times,nbrs),"Fin Ana");

    //print the performance stats
    perf.print_stats();

    //print closing statements
    print_closing(s.world_rank);

    //relinquish the mpi environment
    MPI_Finalize();

    return 0;
}

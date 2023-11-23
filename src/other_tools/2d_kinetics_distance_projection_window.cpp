
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <mpi.h>
#include <sstream>
#include <fstream>

using namespace std;

#include "../headers/multi_dim_vec.h"
#include "../headers/file_reader.h"
#include "../headers/vector_mpi.h"
#include "../headers/common_routines_mpi.h"
#include "../headers/common_routines.h"
#include "../headers/file_naming.h"
#include "../headers/binding_events.h"
#include "../headers/grid_lt.h"
#include "../headers/fit.h"
#include "../headers/command_line_args_mpi.h"
#include "../headers/array.h"
#include "../headers/file_naming_mpi.h"
#include "../headers/performance.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Main function which executes all other functions                                                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, const char * argv[]) 
{
    //Here we define some variables used throughout
    string in_file_name;                //Name of input binding events file 
    string mask_file_name;              //Name of the input file with protein mask
    string base_file_name_o;            //Base file name of output files
    string time_line_file_name;         //Name of the time line file
    string binding_events_file_name;    //Name of the binding_events output file
    string out_file_name;               //Name of the output files with the grid point selection and the big mask
    int i             = 0;              //General variable used in loops
    int j             = 0;              //General variable used in loops
    int k             = 0;              //General variable used in loops
    int q             = 0;              //General variable used in loops
    int target_x      = 0;              //The target grid point x when making a rectangular selection
    int target_y      = 0;              //The target grid point y when making a rectangular selection
    int range_x       = 0;              //The half width of x in the rectangular selection
    int range_y       = 0;              //The half width of y in the rectangular selection
    int iterations    = 1;              //How many times should the grid selection be moved
    int world_size    = 0;              //Size of the mpi world
    int world_rank    = 0;              //Rank in the mpi world
    int invert        = 0;              //Invert rectangular selection (select everything outside rectangle)
    int odf           = 0;              //Data file format
    int start         = 0;              //Start at this iteration
    int cutoff        = 1;              //Must have this many grid points for initial binding
    int width         = 0;              //How many frames of each side of current frame to analyze?
    int counter       = 0;              //How many times the "program run time" been displayed
    double cell_size  = 1.0;            //Distance between grid points
    double res        = 1.0;            //How much does the grid selection move with each iteration 
    double range      = 1.0;            //The half width of the grid point selection                   
    double dt         = 0.0;            //Time step used for converting frames to time. set equal to ef_dt
    double time_mask  = 0.0;            //Time spent generating selection masks
    double time_grid  = 0.0;            //Time spent reading in binding events data
    double time_main  = 0.0;            //Time spent in the main loop
    clock_t t;                          //Keeps the time for testing performance
    clock_t t_grid;                     //Keeps the time for testing performance (grid)
    clock_t t_mask;                     //Keeps the time for testing performance (mask)
    sv1d cl_tags;                       //Holds a list of command line tags for the program

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set up the mpi environment                                                                                //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    MPI_Init(NULL, NULL);;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);      //get the world size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);      //get the process rank

    //create object for logging performance data
    Performance perf;

    //take the initial time
    t = clock();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set program name and print info                                                                           //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    string program_name = "2D Kinetics Distance Projection Window";
    if(world_rank == 0)
    {
        print_credits(argc,argv,program_name);
    }

    string program_description = "2D Kinetics Distance Projection Window is an analysis tool that lets the user compute the lipid residence time as a function of distance from the protein surface, i.e., in thin shells around the protein.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments_mpi(argc,argv,world_rank,program_description);
    add_argument_mpi_s(argc,argv,"-d"     , in_file_name,               "Input binding events file (be)"                                                       , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-mask"  , mask_file_name,             "Input data file with protein mask (dat)"                                              , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o"     , base_file_name_o,           "Filename for deriving output data filenames (dat)"                                    , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-x"     , &target_x,                  "Rectangle center x (grid point)"                                                      , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-y"     , &target_y,                  "Rectangle center y (grid point)"                                                      , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-rx"    , &range_x,                   "Rectangle half width x (grid points)  "                                               , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-ry"    , &range_y,                   "Rectangle half width y (grid points)  "                                               , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-invert", &invert,                    "Invert rectangular selection? (0:no 1:yes)"                                           , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-iter"  , &iterations,                "How many iterations to perform"                                                       , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-res"   , &res,                       "Distance between each iteration (nm)"                                                 , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-range" , &range,                     "Half width of the grid selection shell (nm)"                                          , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-start" , &start,                     "Start at this iteration"                                                              , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-cutoff", &cutoff,                    "The lipid must occupy this many grid points in the mask befor being counted as bound" , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-odf"   , &odf,                       "Data file format (0:matrix 1:vector)"                                                 , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-n"     , &width,                     "Noise filter half width (frames)"                                                     , world_rank, cl_tags, nullptr,      1);
    conclude_input_arguments_mpi(argc,argv,world_rank,program_name,cl_tags);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension_mpi(world_rank,"-d",in_file_name,".be");
    check_extension_mpi(world_rank,"-o",base_file_name_o,".dat");
    check_extension_mpi(world_rank,"-mask",mask_file_name,".dat");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read in protein masking data                                                                              //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Grid_lt prot_mask;
    prot_mask.set_format(odf);
    prot_mask.get_grid(mask_file_name);
    if(world_rank == 0)
    {
        prot_mask.print_dim(1);
        printf("\n");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read binding events file 0_0 to get the header info                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Binding_events events_ref;
    int result    = events_ref.get_info(in_file_name);

    if(result == 0)
    {
        if(world_rank == 0)
        {
            printf("unable to open binding events file %s \n",in_file_name.c_str());
        }
        MPI_Finalize();
        return 0;
    }
    else 
    {
	result = events_ref.get_binding_events_xy(in_file_name,0,0);
        events_ref.get_binding_timeline();
    }

    cell_size = sqrt(events_ref.APS);
    dt        = events_ref.ef_dt; //use units of ps

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Distribute the workload across the cores                                                                  //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int my_iterations = count_workload(world_size,world_rank,iterations);

    //create array to hold each mpi processes iterations; Used for communication
    int world_iterations_ary[world_size];
    MPI_Allgather(&my_iterations, 1,MPI_INT,world_iterations_ary, 1, MPI_INT, MPI_COMM_WORLD );

    //allocate memory for world_iterations and copy data from the array
    iv1d world_iterations(world_size,0);
    for(i=0; i<world_size; i++)
    {
        world_iterations[i] = world_iterations_ary[i];
    }

    //print stats for distributing the iterations and distribute the iterations to each core
    int my_i = 0;
    int my_f = 0;
    int world_i[world_size];
    int world_f[world_size];
    get_workload(&my_i,&my_f,world_rank,world_iterations,iterations,world_i,world_f);

    my_i = my_i + start;
    my_f = my_f + start;
    for(i=0; i<world_size; i++)
    {
        world_i[i] = world_i[i] + start;
        world_f[i] = world_f[i] + start;
    }
    print_workload_stats(world_rank,world_i,world_f,world_iterations,world_size,"iterations","init","fin");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Make arrays for storing dwell times etc.                                                                  //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //local arrays
    double my_dwell_times[my_iterations];  //holds dwell times for the cores iterations
    double my_r2[my_iterations];           //holds r2 for the cores iterations
    double my_koff[my_iterations];         //holds koff for the cores iterations
    double my_num_events[my_iterations];   //holds number of binding events for the cores iterations

    //global arrays
    double world_dwell_time[iterations];   //holds the dwell times for all iterations
    double world_r2[iterations];           //holds the r2 for all iterations
    double world_koff[iterations];         //holds the koff for all iterations
    double world_num_events[iterations];   //holds the number of binding events for all iterations

    //log time spent performing main analysis
    perf.log_time((clock() - t)/CLOCKS_PER_SEC,"Other");

    //reset the clock
    t = clock();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Move the grid selection with each iteration (q)                                                           //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(my_iterations > 0)
    {
        for(q=my_i; q<=my_f; q++) //loop over the iterations
        {
            //print update so the user knows some progress has been made. 
            if(world_rank == 0)
            {
                printf("Building a selection map for iteration %d \n\n",q);
            }

            int effective_iteration = q - my_i;                //this is the index for storing local dwell times

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                                                           //
            // Create mask small                                                                                         //
            //                                                                                                           //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            Grid_lt small_mask;
            small_mask.set_format(odf);
            small_mask.get_grid(mask_file_name);
            small_mask.init_grid(0.0);              //initialize grid to zero
            small_mask.set_nan(0);                  //initialize nan tags

            t_mask = clock();
            small_mask.distance_projection(invert,target_x,target_y,range_x,range_y,cell_size,q,range,res,prot_mask.grid,0);
            time_mask = time_mask + clock() - t_mask;

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                                                           //
            // Write small grid to output file                                                                           //
            //                                                                                                           //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            out_file_name = add_tag(base_file_name_o, "_" + to_string(q) + "_small_mask");
            small_mask.write_grid(out_file_name);

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                                                           //
            // Create a binding events object for the small mask                                                         //
            //                                                                                                           //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            Binding_events events_small;
            result = events_small.get_info(in_file_name);
            int first_small = 1;

            //reset the clock
            t_grid = clock();

            if(world_rank == 0)
            {
                printf("Looping over the grid and building a binding timeline \n");
                printf("-----------------------------------------------------------------------------------------------------------------------------------\n");
            }
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                                                           //
            // Open select files for the current grid selecton (q) and make a time line                                  //
            //                                                                                                           //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            int mask_count = 0;   
	    for(i=0; i<prot_mask.size_y(); i++) //loop over y
            {
                for(j=0; j<prot_mask.size_x(); j++) //loop over x
                {
                    if(small_mask.grid[j][i][2][0] == 1) //check current iteration (q) mask
                    {
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        //                                                                                                           //
                        // Read in binding events                                                                                    //
                        //                                                                                                           //
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        result = events_small.get_binding_events_xy(in_file_name,j,i);

                        if(result == 1)
                        { 
		            if(events_small.lipid_nr.size() > 0)
                            {		
                                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                //                                                                                                           //
                                // Update bound_time_line_ij                                                                                 //
                                //                                                                                                           //
                                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                if(first_small == 1)
                                {
                                    events_small.get_binding_timeline();
                                    first_small = 0;
                                }
                                else
                                {
                                    events_small.stamp_to_binding_timeline();
                                }
                                //printf("Working on grid point x = %5d, y = %5d time_line_lipids %10d lipid_nr.size() %d \n",i,j,events_small.time_line_res_nr.size(),events_small.lipid_nr.size());
                            }
                        }
			else 
		        {
                            printf("failed to open binding events file %s \n",in_file_name);
                        }	
 
                        mask_count++;	 
                    }

		    //report progress and estimated time to completion
                    int current_step = i*prot_mask.size_x() + j + 1;
                    int my_steps     = prot_mask.size_x()*prot_mask.size_y();
                    ot_time_stats(t_grid,&counter,current_step,my_steps,world_rank,"lattice point");
                }
            }
            time_grid = time_grid + clock() - t_grid;

	    counter = 0;

            if(mask_count > 0) //only analyze kinetics if a selection mask was found
            {
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // analyze the time lines and make a binding events file for output                                          //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                Binding_events events;
                events.num_lipids = events_ref.num_lipids;
                events.ef_frames  = events_ref.ef_frames;
                events.x_i        = -1;
                events.y_i        = -1;
                events.num_g_x    = events_ref.num_g_x;
                events.num_g_y    = events_ref.num_g_y;
                events.ef_dt      = events_ref.ef_dt;
                events.APS        = events_ref.APS;
                events.size_timeline();

                for(i=0; i<events.num_lipids; i++) //loop over time line lipids
                {   
                    for(j=0; j<events.ef_frames; j++) //loop over time line frames
                    {
                        if(j >= width && j < events.ef_frames - width)  
                        {
                            int bound_count = 0;

                            for(k=j-width; k<=j+width; k++) //loop over window
                            {
                                if(events_small.bound_time_line[k][i] >= cutoff) //lipid is bound
                                {
                                    bound_count = bound_count + 1;
                                }
                            }
                            double percent = (double)bound_count/(double)(2*width + 1);  

                            if(percent > 0.5) //lipid is bound
                            {
                                events.bound_time_line[j][i] = 1;
                            }
                            else
                            {
                                events.bound_time_line[j][i] = 0;
                            }
                        }
                    }
                }
                events.time_line_lipid_nr = events_small.time_line_lipid_nr;
                events.time_line_res_nr   = events_small.time_line_res_nr;
                events.time_line_res_name = events_small.time_line_res_name;
                events.binding_events_from_timeline();
 
                if(events.lipid_nr.size() > 0)
                {
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //                                                                                                           //
                    // Sort binding events by bind_i (in perparation of writing binding_events file)                             //
                    //                                                                                                           //
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    events.organize_events(0);

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //                                                                                                           //
                    // Print current iteration (q) binding events and timeline to output file                                    //
                    //                                                                                                           //
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    binding_events_file_name = chop_and_add_tag(base_file_name_o, "_" + to_string(q) + ".be");
                    time_line_file_name      = chop_and_add_tag(base_file_name_o, "_" + to_string(q) + "_time_line.dat");
                    events.write_binding_events_bin(binding_events_file_name);
                    events.write_time_line(time_line_file_name);

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //                                                                                                           //
                    // Compute average dwell time                                                                                //
                    //                                                                                                           //
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    double avg_dwell_t = 0;     //The average dwell time in the selection for the current iteration (q)

                    for(i=0; i<events.lipid_nr.size(); i++) //loop over binding events
                    {
                        avg_dwell_t = avg_dwell_t + events.dwell_t[i];
                    }
                    avg_dwell_t = avg_dwell_t/(double)events.lipid_nr.size();

                    //store local average dwell time
                    my_dwell_times[effective_iteration] = avg_dwell_t;

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //                                                                                                           //
                    // Sort binding events by dwell_t (in perparation of computing koff)                                         //
                    //                                                                                                           //
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    events.organize_events(1);

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //                                                                                                           //
                    // Compute frequencies                                                                                       //
                    //                                                                                                           //
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    int largest_dwell_t = events.dwell_t[0];
                    double dwell_time_freq[largest_dwell_t-1];
                    init_darray(dwell_time_freq,largest_dwell_t-1);

                    for(i=0; i<largest_dwell_t; i++) //loop over dwell times
                    {
                        for(j=0; j<events.bind_i.size(); j++) //loop over binding events
                        {
                            if(events.dwell_t[j] == i+1) //dwell time matches
                            {
                                dwell_time_freq[i] = dwell_time_freq[i] + 1;
                            }
                        }
                        dwell_time_freq[i] = dwell_time_freq[i]/(double)events.bind_i.size();
                        //printf("dwell_time_freq %10.8f \n",dwell_time_freq[i]);
                    }

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //                                                                                                           //
                    // Compute probability of lasting t or longer. (add probabilities for t and longer)                          //
                    //                                                                                                           //
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    double nan = -999999;
                    int count = 0;
                    for(i=0; i<largest_dwell_t; i++) //loop over dwell times
                    {
                        if(dwell_time_freq[i] > 0) //inclue data in fit
                        {
                            for(j=i+1; j<largest_dwell_t; j++) //add probabilities
                            {
                                dwell_time_freq[i] = dwell_time_freq[i] + dwell_time_freq[j];
                            }
                            count++;
                        }
                        else //exclue data from fit
                        {
                            dwell_time_freq[i] = nan;
                        }
                    }

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //                                                                                                           //
                    // Get Ln(freq) and time                                                                                     //
                    //                                                                                                           //
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //create array to hold ln(frequencies) excluding Ln(0)
                    double log_freq[count];
                    double time[count];

                    //reset count
                    count = 0;

                    //fill array excluding Ln(0)
                    for(i=0; i<largest_dwell_t; i++) //loop over dwell times
                    {
                        if(dwell_time_freq[i] != nan)
                        {
                            log_freq[count] = log(dwell_time_freq[i]);
                            time[count]     = (i+1)*dt;
                            count++;
                        }
                    }

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //                                                                                                           //
                    // Compute koff                                                                                              //
                    //                                                                                                           //
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    double slope = 0;
                    double yint  = 0;
                    double r2    = 0;
                    double koff  = 0;

                    //do least squares regression
                    least_squares_regression(count,time,log_freq,&slope,&yint,&r2);
                    koff = -slope;

                    //store local r2
                    my_r2[effective_iteration] = r2;

                    //store local koff
                    my_koff[effective_iteration] = koff;

                    //store local num events
                    my_num_events[effective_iteration] = (double)events.lipid_nr.size();  

                    //printf("iteration %d r2 %f koff %f dt %f lipid_nr.size() %d count %d \n",q,r2,koff,dt,events.lipid_nr.size(),count,time);
                }
                else  //no binding events detected
                {
                    printf("Iteration %d was unable to detect any binding events. Perhaps the selection is too small posibly because it is close to the grid boundaries. Setting kintetics data to -1 for this iteration. \n",q);
                    my_dwell_times[effective_iteration] = -1/dt;
                    my_r2         [effective_iteration] = -1;
                    my_koff       [effective_iteration] = -1;
                    my_num_events [effective_iteration] = 0.0;
                }
            }
            else //too many itterations (probably)
            {
                printf("Iteration %d was unable to generate a selection mask. Perhaps the selection has moved beyond the grid boundaries. Setting kintetics data to -1 for this iteration. \n",q);
                my_dwell_times[effective_iteration] = -1.0/dt;
                my_r2         [effective_iteration] = -1.0;
                my_koff       [effective_iteration] = -1.0;
                my_num_events [effective_iteration] = -1.0;
            }

            if(world_rank == 0)
            {
                printf("\n");
            } 
	}
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //get time spent in main loop
    time_main = clock() - t;

    //log time spent performing main analysis
    perf.log_time(time_mask/CLOCKS_PER_SEC,"Selection Masks");
    perf.log_time(time_grid/CLOCKS_PER_SEC,"Binding Timelines");
    perf.log_time((time_main - time_mask - time_grid)/CLOCKS_PER_SEC,"Load Imbalance");

    //reset the clock
    t = clock();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Collect dwell times, r2 and koff values and print them to output files                                    //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(world_rank == 0)
    {
        printf("Collecting grid data and writing output files \n");
        printf("-----------------------------------------------------------------------------------------------------------------------------------\n");
    }

    collect_local_values(world_size,world_rank,world_dwell_time,my_dwell_times,world_iterations,"dwell times (ps)");
    collect_local_values(world_size,world_rank,world_r2,my_r2,world_iterations,"r2");
    collect_local_values(world_size,world_rank,world_koff,my_koff,world_iterations,"koff (ps^-1)");
    collect_local_values(world_size,world_rank,world_num_events,my_num_events,world_iterations,"number of binding events");

    if(world_rank == 0)
    {
        //create name for output files
        string dwell_time_file_name = add_tag(base_file_name_o, "_dwell_time");
        string r2_file_name         = add_tag(base_file_name_o, "_r2");
        string koff_file_name       = add_tag(base_file_name_o, "_koff");
        string num_events_file_name = add_tag(base_file_name_o, "_num_events");

        //open output files
        FILE *dwell_time_file = fopen(dwell_time_file_name.c_str(), "w");
        FILE *r2_file         = fopen(r2_file_name.c_str(), "w");
        FILE *koff_file       = fopen(koff_file_name.c_str(), "w");
        FILE *events_file     = fopen(num_events_file_name.c_str(), "w");

        //check that files were successfully opened. needed if directory does not exist. 
        if(dwell_time_file == NULL)
        {
            printf("Failure opening %s for writting \n",dwell_time_file_name.c_str());
        }
        else if(dwell_time_file == NULL)
        {
            printf("Failure opening %s for writting \n",r2_file_name.c_str());
        }
        else if(dwell_time_file == NULL)
        {
            printf("Failure opening %s for writting \n",koff_file_name.c_str());
        }
        else if(events_file == NULL)
        {
            printf("Failure opening %s for writting \n",num_events_file_name.c_str());
        }
        else
        {
            //print dwell times, r2 and koff to output file
            fprintf(dwell_time_file,"# %15s %20s \n","distance (nm)","dwell time (ps)");
            fprintf(r2_file,"# %15s %20s \n","distance (nm)","r2");
            fprintf(koff_file,"# %15s %20s \n","distance (nm)","koff (ps^-1)");
            fprintf(events_file,"# %15s %20s \n","distance (nm)","binding events");
            for(i=0; i<iterations; i++)
            {
                fprintf(dwell_time_file,"  %15f %20f \n",(double)(i+start)*res,world_dwell_time[i]*dt);
                fprintf(r2_file,"  %15f %20f \n",(double)(i+start)*res,world_r2[i]);
                fprintf(koff_file,"  %15f %20f \n",(double)(i+start)*res,world_koff[i]);
                fprintf(events_file,"  %15f %20f \n",(double)(i+start)*res,world_num_events[i]);
            }

            //close output files 
            fclose(dwell_time_file);
            fclose(r2_file);
            fclose(koff_file);
            fclose(events_file);
        }
    }

    //log time spent performing main analysis
    perf.log_time((clock() - t)/CLOCKS_PER_SEC,"Collect Data");

    //print the performance stats
    perf.print_stats();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Print final output and end protram                                                                        //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(world_rank == 0)
    {
        std::cout << "\nFormatting completed successfully" << "\n\n";
    }

    MPI_Finalize();

    return 0;
}



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
#include "../headers/binding_events_common_routines.h"
#include "../headers/common_routines_mpi.h"
#include "../headers/common_routines.h"
#include "../headers/file_naming.h"
#include "../headers/binding_events.h"
#include "../headers/command_line_args_mpi.h"
#include "../headers/array.h"
#include "../headers/file_naming_mpi.h"
#include "../headers/performance.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// The main function performing analysis                                                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, const char * argv[]) 
{
    //Here we define some variables used throughout
    string in_file_name;          //Name of the binding events file being read
    int i            = 0;         //Standard variable used in loops
    int j            = 0;         //Standard variable used in loops
    int k            = 0;         //Standard variable used in loops
    int l            = 0;         //Standard variable used in loops
    int world_size   = 0;         //How many mpi ranks
    int world_rank   = 0;         //Rank of the mpi process
    int counter      = 0;         //How many times the "program run time" been displayed
    int grid_counter = 0;         //Count lattice points as they are encountered
    double dt        = 0;         //Time step used for output 
    clock_t t;                    //Keeps the time for testing performance
    sv1d cl_tags;                 //Holds a list of command line tags for the program

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
    string program_name = "Binding Contributors";
    if(world_rank == 0)
    {
        print_credits(argc,argv,program_name);
    }

    string program_description = "Binding Contributors is an analysis tool that lets the user compute the average number of repeat visits and the time between visits from binding events data.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments_mpi(argc,argv,world_rank,program_description);
    add_argument_mpi_s(argc,argv,"-d"       , in_file_name,             "Input binding events file (be)"         , world_rank, cl_tags, nullptr,      1);
    conclude_input_arguments_mpi(argc,argv,world_rank,program_name,cl_tags);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension_mpi(world_rank,"-d",in_file_name,".be");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read binding events file 0_0 to get the header info                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Binding_events events;
    int result = events.get_info(in_file_name);

    if(result == 0)
    {
        if(world_rank == 0)
        {
            printf("unable to open binding events file %s \n",in_file_name.c_str());
        }
        MPI_Finalize();
        return 0;
    }

    dt = events.ef_dt;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Distribute the workload across the cores                                                                  //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int my_num_g = count_workload(world_size,world_rank,events.num_g_x*events.num_g_y);

    //create array to hold each mpi processes my_num_g; Used for communication
    int world_num_g_ary[world_size];
    MPI_Allgather(&my_num_g, 1,MPI_INT,world_num_g_ary, 1, MPI_INT, MPI_COMM_WORLD );

    //allocate memory for world_num_g and copy data from the array
    iv1d world_num_g(world_size,0);
    for(i=0; i<world_size; i++)
    {
        world_num_g[i] = world_num_g_ary[i];
    }

    //print stats for distributing the grid and distribute the grid to each core
    int my_gi = 0;
    int my_gf = 0;
    iv1d world_gi(world_size);
    iv1d world_gf(world_size);
    get_grid_points_alt(&my_gi,&my_gf,world_rank,world_size,world_num_g,events.num_g_x,events.num_g_y,world_gi,world_gf);
    print_workload_stats_alt(world_rank,world_gi,world_gf,world_num_g,world_size);
    MPI_Barrier(MPI_COMM_WORLD);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Allocate memory for the grid                                                                              //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    dv1d avg_visits_local(my_num_g,0.0);
    dv1d avg_off_time_local(my_num_g,0.0);

    if(world_rank == 0)
    {
        printf("Reading binding events data and analyzing repeat visits. \n");
        printf("-----------------------------------------------------------------------------------------------------------------------------------\n");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // loop over the grid                                                                                        //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(i=0; i<=events.num_g_x; i++) //loop over x
    {
        for(j=0; j<events.num_g_y; j++) //loop over y
        {
            if(grid_counter >= my_gi && grid_counter <= my_gf)
            {
		int ef_g = grid_counter - my_gi; 
 
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // Read in binding events                                                                                    //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                result       = events.get_binding_events_xy(in_file_name,i,j);

                if(result == 1) //binding events file exists
                { 
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //                                                                                                           //
                    // Characterize repeated binding                                                                             //
                    //                                                                                                           //
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    events.reduce_list();

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //                                                                                                           //
                    // Get the average off time and average number of repeat visits for the grid point                           //
                    //                                                                                                           //
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    double average_num_repeats = 0;    //the average number of repeat visits for the grid point
                    double average_off_time    = 0;    //the average time between repeat visits for the grid point

                    for(k=0; k<events.repeats_res_nr.size(); k++) //loop over repeats lipids
                    {
                        average_num_repeats = average_num_repeats + events.repeats[k];
                        average_off_time    = average_off_time + events.repeats_avg_off_time[k];       
                    }
                    average_num_repeats = average_num_repeats/events.repeats_res_nr.size();
                    average_off_time    = average_off_time/events.repeats_res_nr.size();

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //                                                                                                           //
                    // Add averages to the grid.                                                                                 //
                    //                                                                                                           //
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    avg_visits_local[ef_g]   = average_num_repeats;
                    avg_off_time_local[ef_g] = average_off_time*dt;          
                }

                //report progress and estimated time to completion
                int current_step = ef_g + 1;
                int my_steps     = my_gf - my_gi + 1;
                ot_time_stats(t,&counter,current_step,my_steps,world_rank,"lattice point");
            }
            grid_counter++;
        }
    }

    //log time spent performing main analysis
    perf.log_time((clock() - t)/CLOCKS_PER_SEC,"Main Loop");

    //take the initial time
    t = clock();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Collect data from cores and write the data to output files                                                //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    dv2d avg_visits_global(events.num_g_y,dv1d(events.num_g_x,0.0));
    dv2d avg_off_time_global(events.num_g_y,dv1d(events.num_g_x,0.0));
    iv2d nan(events.num_g_y,iv1d(events.num_g_x,0));

    gather_grid_d_gp_alt(world_size,world_rank,world_gi,world_gf,events.num_g_x,events.num_g_y,world_num_g,avg_visits_local,avg_visits_global);
    gather_grid_d_gp_alt(world_size,world_rank,world_gi,world_gf,events.num_g_x,events.num_g_y,world_num_g,avg_off_time_local,avg_off_time_global);

    if(world_rank == 0)
    {
        string nv_out_file_name = chop_and_add_tag(in_file_name, "_num_visits.dat");
        string ot_out_file_name = chop_and_add_tag(in_file_name, "_off_time.dat");

        write_grid_to_file(events.num_g_x,events.num_g_y,nan,nv_out_file_name,avg_visits_global);
        write_grid_to_file(events.num_g_x,events.num_g_y,nan,ot_out_file_name,avg_off_time_global);
    }

    //log time spent collecting data
    perf.log_time((clock() - t)/CLOCKS_PER_SEC,"Collect Data");

    //print the performance stats
    perf.print_stats();

    if(world_rank == 0)
    {
        std::cout << "\nFormatting completed successfully" << "\n\n";
    }

    MPI_Finalize();

    return 0;
}


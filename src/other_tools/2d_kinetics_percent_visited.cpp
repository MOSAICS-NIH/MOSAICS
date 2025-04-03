
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
#include "../headers/file_naming.h"
#include "../headers/binding_events_common_routines.h"
#include "../headers/common_routines_mpi.h"
#include "../headers/common_routines.h"
#include "../headers/binding_events.h"
#include "../headers/index.h"                             //This has a class for working with index files
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
    string base_file_name;        //Base file name for the binding events files
    string lip_t_file_name;       //Lipid types for analysis
    string out_file_name;         //Name of the output file
    int i            = 0;         //Standard variable used in loops
    int j            = 0;         //Standard variable used in loops
    int k            = 0;         //Standard variable used in loops
    int l            = 0;         //Standard variable used in loops
    int m            = 0;         //Standard variable used in loops
    int n            = 0;         //Standard variable used in loops
    int freq         = 1;         //How often to report percent visited
    int world_size   = 0;         //How many mpi ranks
    int world_rank   = 0;         //Rank of the mpi process
    int num_lipids   = 0;         //Number of target lipids
    int b_num_lip    = 0;         //Was the number of lipids specified?
    int counter      = 0;         //How many times the "program run time" been displayed
    int grid_counter = 0;         //Count lattice points as they are encountered
    int grid         = 0;         //Do analysis on a grid? 
    int b_x          = 0;         //Was a grid point in x-direction specified
    int b_y          = 0;         //Was a grid point in y-direction specified
    int x            = 0;         //Grid point in x-direction
    int y            = 0;         //Grid point in y-direction
    double cutoff    = 0.0;       //Exclude data with cutoff less than this
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
    string program_name = "2d Kinetics percent Visited";
    if(world_rank == 0)
    {
        print_credits(argc,argv,program_name);
    }

    string program_description = "2d Kinetics Percent Visited is an analysis tool that lets the user compute the percentage of lipids to visit each grid point using binding events data.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments_mpi(argc,argv,world_rank,program_description);
    add_argument_mpi_s(argc,argv,"-d"       , in_file_name,               "Input binding events file (be)              "                        , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o"       , out_file_name,              "Output filename used to derive names for percent visited data (dat)" , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-freq"    , &freq,                      "How often to report the percent visited (frames)"                    , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-crd"     , lip_t_file_name,            "Selection card with lipid types (crd)"                               , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-lipids"  , &num_lipids,                "How many lipids are there of the target type in the target leaflet?" , world_rank, cl_tags, &b_num_lip,   0);
    add_argument_mpi_i(argc,argv,"-grid"    , &grid,                      "Do analysis on grid data? (0:no, 1:yes)"                             , world_rank, cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-x"       , &x,                         "Grid point in x-direction"                                           , world_rank, cl_tags, &b_x,         0);
    add_argument_mpi_i(argc,argv,"-y"       , &y,                         "Grid point in y-direction"                                           , world_rank, cl_tags, &b_y,         0);
    add_argument_mpi_d(argc,argv,"-cutoff"  , &cutoff,                    "Exclude data with a dwell time smaller than this (ps)"               , world_rank, cl_tags, nullptr,      0);
    conclude_input_arguments_mpi(argc,argv,world_rank,program_name,cl_tags);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension_mpi(world_rank,"-d",in_file_name,".be");
    check_extension_mpi(world_rank,"-o",out_file_name,".dat");
    check_extension_mpi(world_rank,"-crd",lip_t_file_name,".crd");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check if a lattice point was specified                                                                    //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(world_rank == 0)
    {
        if(b_x == 1 && b_y == 0)
        {
            printf("A lattice point was specified for the x-direction but not y. Please include the y-direction if analyzing binding events for a lattice point. \n");
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }
        else if(b_x == 0 && b_y == 1)
        {
            printf("A lattice point was specified for the y-direction but not x. Please include the x-direction if analyzing binding events for a lattice point. \n");
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Create index to hold lipid types                                                                          //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //create index for lipid types
    Index lip_t;

    //read the index files
    lip_t.get_index(lip_t_file_name);

    if(grid == 1)
    {
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Check if a lattice point was specified                                                                    //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(world_rank == 0)
        {
            if(b_x == 1 || b_y == 1)
            {
                printf("A lattice point was specified (-x -y) in combinatioin with the -grid option. These options are incompatible with each other. \n");
                MPI_Finalize();
                exit(EXIT_SUCCESS);
            }
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Read binding events file 0_0 to get the header info                                                       //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        Binding_events events_ref;
        int result = events_ref.get_info(in_file_name);

        if(result == 0)
        {
            if(world_rank == 0)
            {
                printf("unable to open binding events file %s \n",in_file_name.c_str());
            }
            MPI_Finalize();
            return 0;
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Distribute the workload across the cores                                                                  //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        int my_num_g = count_workload(world_size,world_rank,events_ref.num_g_x*events_ref.num_g_y);

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
        get_grid_points_alt(&my_gi,&my_gf,world_rank,world_size,world_num_g,events_ref.num_g_x,events_ref.num_g_y,world_gi,world_gf);
        print_workload_stats_alt(world_rank,world_gi,world_gf,world_num_g,world_size);
        MPI_Barrier(MPI_COMM_WORLD);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Allocate memory for the grid                                                                              //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        int size = 0;
        for(i=0; i<events_ref.ef_frames; i++) //loop over time line frames
        {
            if(i%freq == 0)
            {
                size++;
            }
        }
        dv2d percent_local(size,dv1d(my_num_g,0.0));

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Make binding events object for reading in data                                                            //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        Binding_events events;
        result = events.get_info(in_file_name);

        if(result == 0)
        {
            if(world_rank == 0)
            {
                printf("unable to open binding events file %s \n",in_file_name.c_str());
            }
            MPI_Finalize();
            return 0;
        }

        if(world_rank == 0)
        {
            printf("Reading binding events data and computing percentage of lipids to visit each lattice point. \n");
            printf("-----------------------------------------------------------------------------------------------------------------------------------\n");
        }

        //log time spent performing main analysis
        perf.log_time((clock() - t)/CLOCKS_PER_SEC,"Other");

        //take the initial time
        t = clock();

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
                    result = events.get_binding_events_xy(in_file_name,i,j);

                    if(result == 1) //binding events file exists
                    {
                        if(events.lipid_nr.size() > 0)
                        {
                            events.organize_events(0); //organize by b_i

                            if(b_num_lip == 0)
                            {
                                num_lipids = events.num_lipids;
                            }

                            int ef_frame = -1; //used for adding data to percent_local

                            for(k=0; k<events.ef_frames; k+=freq) //loop over frames
                            {
                                ef_frame++;

                                iv1d lipid_numbers(0,0);    //store lipid numbers as they are encountered

                                double   dt = events.ef_dt; //effective time step between frames
                                double time = dt*k;         //time for current frame

                                for(l=0; l<events.dwell_t.size(); l++) //loop over binding events
                                {
                                    for(m=0; m<lip_t.index_s.size(); m++) //loop over lipid types
                                    {
                                        if(strcmp(events.res_name[l].c_str(), lip_t.index_s[m].c_str()) == 0) //lipid type is correct
                                        {
                                            double dwell_time = (double)events.dwell_t[l]*dt;   //dwell time for the event 
                                            int      lipid_nr = (double)events.lipid_nr[l];     //lipid number for the event
                                            double     bind_i = dt*(double)events.bind_i[l];    //initial binding time for event

                                            if(dwell_time >= cutoff && bind_i <= time)
                                            {
                                                int found = 0;

                                                for(n=0; n<lipid_numbers.size(); n++) //loop over the lipid numbers
                                                {
                                                    if(lipid_numbers[n] == lipid_nr)
                                                    {
                                                        found = 1;
                                                    }
                                                }

                                                if(found == 0)
                                                {
                                                    lipid_numbers.push_back(lipid_nr);
                                                }
                                            }
                                        }
                                    }
                                }

                                int lipids_found     = lipid_numbers.size();
                                double percent_found = (double)lipids_found/(double)num_lipids;

                                //store percent visited
                                percent_local[ef_frame][ef_g] = percent_found;
                            }
                        }
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
        // Collect grids and write data to output files                                                              //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(world_rank == 0)
        {
            printf("\nCollecting grid data and writing output data. \n");
            printf("-----------------------------------------------------------------------------------------------------------------------------------\n");
        }

        int ef_frame = -1;
        for(i=0; i<events_ref.ef_frames; i++) //loop over time line frames
        {
            if(i%freq == 0)
            {
                if(world_rank == 0)
                {
                    printf("Working on frame %d \n",i);
                }

                ef_frame++;

                dv2d percent_global(events_ref.num_g_y,dv1d(events_ref.num_g_x,0.0));
                iv2d nan(events_ref.num_g_y,iv1d(events_ref.num_g_x,0));

                gather_grid_d_gp_alt(world_size,world_rank,world_gi,world_gf,events_ref.num_g_x,events_ref.num_g_y,world_num_g,percent_local[ef_frame],percent_global);

                if(world_rank == 0)
                {
                    string tag            = "_" + to_string(i);
                    string this_file_name = add_tag(out_file_name,tag);

                    write_grid_to_file(events_ref.num_g_x,events_ref.num_g_y,nan,this_file_name,percent_global);
                }
            }
        }

        //log time spent collecting data
        perf.log_time((clock() - t)/CLOCKS_PER_SEC,"Collect Data");
    }
    else if(grid == 0)
    {
        //take the initial time
        t = clock();

        if(world_rank == 0)
        {
            //open file for writing data
            FILE *out_file = fopen(out_file_name.c_str(), "w");

            //print header info
            fprintf(out_file,"%10s %12s %10s %10s %10s \n","#frame","#time(ps)","#visited","#total","#percent");

	    int result = 0;

	    Binding_events events;

            //read binding events data
            if(b_x==1 && b_y==1)
            {
                result = events.get_info(in_file_name);
                if(result == 1)
                {
                    result = events.get_binding_events_xy(in_file_name,x,y);
                }
            }
            else
            {
                result = events.get_binding_events_bin(in_file_name);
            }

            if(result == 0) 
            {
                printf("unable to open binding events file %s \n",in_file_name.c_str());
            }
            else if(result == 1) //binding events file exists
            {
                events.organize_events(0); //organize by b_i

                if(b_num_lip == 0)
                {
                    num_lipids = events.num_lipids;
                }

                int ef_frame = -1; //used for adding data to percent_local

                for(k=0; k<events.ef_frames; k+=freq) //loop over frames
                {
                    ef_frame++;

                    iv1d lipid_numbers(0,0);    //store lipid numbers as they are encountered

                    double   dt = events.ef_dt; //effective time step between frames
                    double time = dt*k;         //time for current frame

                    for(l=0; l<events.dwell_t.size(); l++) //loop over binding events
                    {
                        for(m=0; m<lip_t.index_s.size(); m++) //loop over lipid types
                        {
                            if(strcmp(events.res_name[l].c_str(), lip_t.index_s[m].c_str()) == 0) //lipid type is correct
                            {
                                double dwell_time = (double)events.dwell_t[l]*dt;   //dwell time for the event 
                                int      lipid_nr = (double)events.lipid_nr[l];     //lipid number for the event
                                double     bind_i = dt*(double)events.bind_i[l];    //initial binding time for event

                                if(dwell_time >= cutoff && bind_i <= time)
                                {
                                    int found = 0;

                                    for(n=0; n<lipid_numbers.size(); n++) //loop over the lipid numbers
                                    {
                                        if(lipid_numbers[n] == lipid_nr)
                                        {
                                            found = 1;
                                        }
                                    }

                                    if(found == 0)
                                    {
                                        lipid_numbers.push_back(lipid_nr);
                                    }
                                }
                            }
                        }
                    }

                    int lipids_found     = lipid_numbers.size();
                    double percent_found = (double)lipids_found/(double)num_lipids;

                    //store percent visited
                    fprintf(out_file,"%10d %12.2f %10d %10d %10f \n",k,(double)k*events.ef_dt,lipids_found,num_lipids,percent_found);
                }
            }
            //close output file
            fclose(out_file);
	}

        //log time spent performing main analysis
        perf.log_time((clock() - t)/CLOCKS_PER_SEC,"Main Ana");
    }    

    if(world_rank == 0)
    {
        printf("Number of lipids used for normalizing factor: %d \n",num_lipids);
    }

    //print the performance stats
    perf.print_stats();

    if(world_rank == 0)
    {
        std::cout << "\nFormatting completed successfully" << "\n\n";
    }

    MPI_Finalize();

    return 0;
}


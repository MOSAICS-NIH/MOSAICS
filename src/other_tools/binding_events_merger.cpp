
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <sstream>
#include <mpi.h>
#include <fstream>

using namespace std;

#include "../headers/multi_dim_vec.h"
#include "../headers/file_reader.h"
#include "../headers/vector_mpi.h"
#include "../headers/switch.h"
#include "../headers/common_routines.h"
#include "../headers/common_routines_mpi.h"
#include "../headers/file_naming.h"
#include "../headers/binding_events.h"
#include "../headers/command_line_args_mpi.h"
#include "../headers/array.h"
#include "../headers/force_serial.h"
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
    string binding_events_file_name_1;        //Name of input file 1
    string binding_events_file_name_2;        //Name of input file 2
    string out_file_name;                     //Name of the output file
    int i                    = 0;             //Standard variable used in loops
    int j                    = 0;             //Standard variable used in loops
    int k                    = 0;             //Standard variable used in loops
    int sort                 = 0;             //Sort binding events by bind_i
    int compound_lipid_count = 0;             //Add lipid count for both files
    int lipid_nr_offset      = 0;             //Offset the second events file lipid nr by this amount
    int do_time_line         = 0;             //Get final events from timeline? 
    int world_size           = 0;             //Size of the mpi world
    int world_rank           = 0;             //Rank in the mpi world
    int grid                 = 0;             //Merge BE files for a pair of grids
    int result               = 0;             //Was the BE file read correctly?
    int do_grid              = 0;             //Merge BE files for a grid?
    int counter              = 0;             //How many times the "program run time" been displayed
    int grid_counter         = 0;             //Count lattice points as they are encountered
    clock_t t;                                //Keeps the time for testing performance
    clock_t t_grid;                           //Keeps the time when looping over the grid
    sv1d cl_tags;                             //Holds a list of command line tags for the program

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

    //force program to run in serial?
    enum Switch serial         = on;

    //here we check if the program supports parallelization or not
    check_serial(world_rank,world_size,serial);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set program name and print info                                                                           //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    string program_name = "Binding Events Merger";

    print_credits(argc,argv,program_name);

    string program_description = "Binding Events Merger is an analysis tool used for reading in 2 binding events files and merging them into a single file.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments_mpi(argc,argv,world_rank,program_description);
    add_argument_mpi_s(argc,argv,"-d1"       , binding_events_file_name_1,"Input binding events file 1 (be)"                               , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-d2"       , binding_events_file_name_2,"Input binding events file 2 (be)"                               , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-sort"     , &sort,                     "Sort the binding events by initial binding time (0:no, 1:yes)"  , world_rank, cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-o"        , out_file_name,             "Output data file with merged binding events (be)"               , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-offset"   , &lipid_nr_offset,          "Offset lipid nr in -d2 by this much"                            , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-compnd"   , &compound_lipid_count,     "Compound lipid count? (0:no, 1:yes)"                            , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-tl"       , &do_time_line,             "Get events from time line? (0:no, 1:yes)"                       , world_rank, cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-grid"     , &do_grid,                  "Merge binding events files for a pair of grids? (0:no, 1:yes)"  , world_rank, cl_tags, nullptr,      0);
    conclude_input_arguments_mpi(argc,argv,world_rank,program_name,cl_tags);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension_mpi(world_rank,"-d1",binding_events_file_name_1,".be");
    check_extension_mpi(world_rank,"-d2",binding_events_file_name_2,".be");
    check_extension_mpi(world_rank,"-o",out_file_name,".be");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Merge binding events for either a grid or single BE files                                                 //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(do_grid == 1)
    {
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Analyze first BE file (grid)                                                                              //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        Binding_events events_1;
        result = events_1.get_info(binding_events_file_name_1);
        if(result == 0)
        {
            printf("unable to open binding events file %s \n",binding_events_file_name_1.c_str());    
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Analyze second BE file (grid)                                                                             //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        Binding_events events_2; 
        result = events_2.get_info(binding_events_file_name_2);
        if(result == 0)
        {
            printf("unable to open binding events file %s \n",binding_events_file_name_2.c_str());
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        } 
       
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Check compatability of grids                                                                              //
        //                                                                                                           //
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
        if(events_1.num_g_x != events_2.num_g_x || events_1.num_g_y != events_2.num_g_y)
        {
            printf("Binding events files are incompatible with each other \n");
            printf("%15s: num_g_x %5d num_g_y %5d \n",binding_events_file_name_1.c_str(),events_1.num_g_x,events_1.num_g_y);
            printf("%15s: num_g_x %5d num_g_y %5d \n",binding_events_file_name_2.c_str(),events_2.num_g_x,events_2.num_g_y);
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Loop over grid points and merge the data                                                                  //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        printf("Looping over the grid and merging binding events data. \n");
        printf("-----------------------------------------------------------------------------------------------------------------------------------\n");

        ofstream be_file_o;                                                 //file for writing out binding events data
        i64v2d grid_be_pos(events_1.num_g_x, i64v1d(events_1.num_g_y,0));   //store the position of each grid point in the be file 

        //open file for writing the final binding events data
        be_file_o.open(out_file_name, ios::out | ios::binary);
        int64_t current_pos_o = be_file_o.tellp();

        t_grid = clock();

        for(i=0; i<events_1.num_g_x; i++) //loop over x
        {
            for(j=0; j<events_1.num_g_y; j++) //loop over y
            {
                //read in binding events for the lattice point
                result = events_1.get_binding_events_xy(binding_events_file_name_1,i,j);
                result = events_2.get_binding_events_xy(binding_events_file_name_2,i,j);

                //add binding events to current list
		for(k=0; k<events_2.lipid_nr.size(); k++) //loop over second set of binding events
                {
                    events_1.lipid_nr.push_back(events_2.lipid_nr[k] + lipid_nr_offset);
                    events_1.res_nr.  push_back(events_2.res_nr[k]);
                    events_1.res_name.push_back(events_2.res_name[k]);
                    events_1.bind_i.  push_back(events_2.bind_i[k]);
                    events_1.bind_f.  push_back(events_2.bind_f[k]);
                    events_1.dwell_t. push_back(events_2.dwell_t[k]);
                }

                //compound the number of lipids
                if(compound_lipid_count == 1)
                {
                    events_1.num_lipids = events_1.num_lipids + events_2.num_lipids;
                }

                //store the position of the grid point 
                grid_be_pos[i][j] = current_pos_o;

                //write binding event data to the final file
                current_pos_o = events_1.write_binding_events_tmp(be_file_o,current_pos_o);

                //count lattice points
                grid_counter++;

                //report progress
                int current_step = grid_counter;
                int my_steps     = events_1.num_g_x*events_1.num_g_y;
                ot_time_stats(t_grid,&counter,current_step,my_steps,world_rank,"grid points");
            }
        } 

        //write out .info file for the binding events file
        string info_file_name = out_file_name + ".info"; 
        FILE *this_file = fopen(info_file_name.c_str(), "w");
        if(this_file == NULL)
        {
            printf("Could not open file %s. \n",info_file_name.c_str());
        }
        else
        {
            for(i=0; i<events_1.num_g_x; i++) //loop over x
            {
                for(j=0; j<events_1.num_g_y; j++) //loop over y 
                {
                    fprintf(this_file," %ld ",grid_be_pos[i][j]);
                }
                fprintf(this_file,"\n");
            }
            fclose(this_file);
        }
    }
    else //single binding events file
    {
        Binding_events events;
        result = events.get_binding_events_bin(binding_events_file_name_1);
        if(result == 0)
        {
            printf("unable to open binding events file %s \n",binding_events_file_name_1.c_str());
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }
        result = events.add_binding_events(binding_events_file_name_2,compound_lipid_count,lipid_nr_offset);
        if(result == 0)
        {
            printf("unable to open binding events file %s \n",binding_events_file_name_2.c_str());
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }


        if(events.res_nr.size() > 0)
        {
            if(do_time_line == 1)
            {
                events.get_binding_timeline();
                events.binding_events_from_timeline();
            }

            if(sort == 1)
            {
                events.organize_events(0);
            }

            events.write_binding_events_bin(out_file_name);
        }
    }

    //log time spent performing misc tasks
    perf.log_time((clock() - t)/CLOCKS_PER_SEC,"Main");

    //print the performance stats
    perf.print_stats();

    //relinquish the mpi environment
    MPI_Finalize();

    std::cout << "\nFormatting completed successfully" << "\n\n";

    return 0;
}



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
#include "../headers/binding_events.h"
#include "../headers/command_line_args_mpi.h"
#include "../headers/array.h"
#include "../headers/file_naming_mpi.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// The main function performing analysis                                                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, const char * argv[]) 
{
    //Here we define some variables used throughout
    FILE *out_file;               //File for writing data to
    string in_file_name;          //Name of the input file
    string out_file_name;         //Name of the output file
    string base_file_name;        //Name of the input file (base)
    int i          = 0;           //General variable used in loops
    int j          = 0;           //General variable used in loops
    int k          = 0;           //General variable used in loops
    int l          = 0;           //General variable used in loops
    int world_size = 0;           //Size of the mpi world
    int world_rank = 0;           //Rank in the mpi world

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set up the mpi environment                                                                                //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    MPI_Init(NULL, NULL);;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);      //get the world size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);      //get the process rank

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set program name/description and print info                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    string program_name = "Binding Occupancy";
    if(world_rank == 0)
    {
        print_credits(argc,argv,program_name);
    }

    string program_description = "Binding Occupancy is an analysis tool that lets the user compute the percentage of time in which a lipid is bound for each grid point. This is done using binding events data produced by 2D Kinetics.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments_mpi(argc,argv,world_rank,program_description);
    add_argument_mpi_s(argc,argv,"-d"      , base_file_name,             "Base filename for input binding events files"            , world_rank, nullptr,      1);
    conclude_input_arguments_mpi(argc,argv,world_rank,program_name);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read binding events file 0_0 to get the header info                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Binding_events events;
    in_file_name = base_file_name + "_" + to_string(0) + "_" + to_string(0) + ".be";
    int result   = events.get_binding_events(in_file_name);

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
    int my_num_g_x = count_workload(world_size,world_rank,events.num_g_x);

    //create array to hold each mpi processes my_num_g_x; Used for communication
    int world_num_g_x_ary[world_size];
    MPI_Allgather(&my_num_g_x, 1,MPI_INT,world_num_g_x_ary, 1, MPI_INT, MPI_COMM_WORLD );

    //allocate memory for world_num_g_x and copy data from the array
    iv1d world_num_g_x(world_size,0);
    for(i=0; i<world_size; i++)
    {
        world_num_g_x[i] = world_num_g_x_ary[i];
    } 

    //print stats for distributing the grid and distribute the grid to each core
    int my_xi = 0;
    int my_xf = 0;
    int world_xi[world_size];
    int world_xf[world_size];
    get_workload(&my_xi,&my_xf,world_rank,world_num_g_x,events.num_g_x,world_xi,world_xf);
    print_workload_stats(world_rank,world_xi,world_xf,world_num_g_x,world_size,"num_g_x","xi","xf");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Allocate memory for the grid                                                                              //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    dv2d occupancy_local(events.num_g_y,dv1d(my_num_g_x,0.0));

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // loop over the grid                                                                                        //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(i=my_xi; i<=my_xf; i++) //loop over x
    {
        int ef_x = i - my_xi;

        if(world_rank == 0)
        {   
            printf("Working on grid collumn %d (%d) \n",ef_x,i);
        }

        for(j=0; j<events.num_g_y; j++) //loop over y
        {
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                                                           //
            // Read in binding events                                                                                    //
            //                                                                                                           //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            in_file_name = base_file_name + "_" + to_string(i) + "_" + to_string(j) + ".be";
            result       = events.get_binding_events(in_file_name);

            if(result == 1) //binding events file exists
            {
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // create binding time line                                                                                  //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                events.get_binding_timeline();

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // compute average number bound at any given time                                                            //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                double percent_occupied = 0;  //what percentage of frames have a lipid bound
                int occupied            = 0;  //is there a lipid bound for the frame
                int count               = 0;  //total number of frames with a lipid bound 

                for(k=0; k<events.ef_frames; k++) //loop over frames
                {
                    occupied = 0;
                    for(l=0; l<events.num_lipids; l++) //loop over lipids
                    {
                        if(events.bound_time_line[k][l] == 1)
                        {
                            occupied  = 1;
                        }
                    }
                    if(occupied == 1)
                    {
                        count++;
                    }
                }
                percent_occupied = (double)count/(double)events.ef_frames;

                //add averages to the grid
                occupancy_local[j][ef_x] = percent_occupied;
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Collect data from cores and write the data to output files                                                //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    dv2d occupancy_global(events.num_g_y,dv1d(events.num_g_x,0.0));
    iv2d nan(events.num_g_y,iv1d(events.num_g_x,0));

    gather_grid_d_gp(world_size,world_rank,my_num_g_x,events.num_g_x,events.num_g_y,world_num_g_x,occupancy_local,occupancy_global);

    if(world_rank == 0)
    {
        string oc_out_file_name    = base_file_name + "_occupancy.dat";

        write_grid_to_file(events.num_g_x,events.num_g_y,nan,oc_out_file_name,occupancy_global);
    }

    if(world_rank == 0)
    {
        std::cout << "\nFormatting completed successfully" << "\n\n";
    }

    MPI_Finalize();

    return 0;
}



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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// The main function performing analysis                                                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, const char * argv[]) 
{
    //Here we define some variables used throughout
    string binding_events_file_name;  //Name of the binding events input file   
    string out_file_name;             //Name of the binding timeline output file
    int i          = 0;               //General variable used in loops
    int j          = 0;               //General variable used in loops
    int world_size = 0;               //Size of the mpi world
    int world_rank = 0;               //Rank in the mpi world
    int b_x        = 0;               //Was a grid point in x-direction specified
    int b_y        = 0;               //Was a grid point in y-direction specified
    int x          = 0;               //Grid point in x-direction
    int y          = 0;               //Grid point in y-direction
    int result     = 0;               //Tells if the binding events file was read sucressfully
    sv1d cl_tags;                     //Holds a list of command line tags for the program

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set up the mpi environment                                                                                //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    MPI_Init(NULL, NULL);;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);      //get the world size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);      //get the process rank

    //force program to run in serial?
    enum Switch serial         = on;

    //here we check if the program supports parallelization or not
    check_serial(world_rank,world_size,serial);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set program name and print info                                                                           //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    string program_name = "Binding Time Line";

    print_credits(argc,argv,program_name);

    string program_description = "Binding Time Line is an analysis tool used for reading in a binding events file and making a binding time line.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments_mpi(argc,argv,world_rank,program_description);
    add_argument_mpi_s(argc,argv,"-d"        , binding_events_file_name,  "Binding events file (be)"                              , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o"        , out_file_name,             "Output data file with binding events timeline (dat)"   , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-x"        , &x,                        "Grid point in x-direction"                             , world_rank, cl_tags, &b_x,         0);
    add_argument_mpi_i(argc,argv,"-y"        , &y,                        "Grid point in y-direction"                             , world_rank, cl_tags, &b_y,         0);
    conclude_input_arguments_mpi(argc,argv,world_rank,program_name,cl_tags);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension_mpi(world_rank,"-d",binding_events_file_name,".be");
    check_extension_mpi(world_rank,"-o",out_file_name,".dat");

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
    // Read in binding events                                                                                    //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Binding_events events;

    if(b_x == 1 && b_y == 1)
    {
        result = events.get_info(binding_events_file_name);
        if(result == 1)
        {
            result = events.get_binding_events_xy(binding_events_file_name,x,y);
        }
    }
    else 
    { 
        result = events.get_binding_events_bin(binding_events_file_name);
    }

    if(result == 1) //binding events file exists
    {
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Make bound_time_line_ij                                                                                   //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    	    events.get_binding_timeline();

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // compute average number bound and the percentage of frames where a lipid is bound                          //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double percent_occupied = 0;      //what percentage of the frames is there a any lipid bound
        double avg_bound        = 0;      //how many lipids are bound on average 
        int occupied            = 0;      //tells is any lipid is bound for the frame
        int count               = 0;      //count how many frames a lipid is bound

        for(i=0; i<events.ef_frames; i++) //loop over frames
        {
            occupied = 0;
            for(j=0; j<events.num_lipids; j++) //loop over lipids
            {
                if(events.bound_time_line[i][j] == 1)
                {
                    occupied  = 1;
                    avg_bound = avg_bound + events.bound_time_line[i][j];
                }
            }
            if(occupied == 1)
            {
                count++;
            }
        }
        avg_bound        = avg_bound/((double)count);
        percent_occupied = (double)count/(double)events.ef_frames;

        printf("avg_bound %f percent_occupied %f \n",avg_bound,percent_occupied);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Write bound time line to file                                                                             //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        events.write_time_line(out_file_name); 
    }
    else //binding events file does not exist 
    {
        printf("Could not find binding events file \n");
    }

    //relinquish the mpi environment
    MPI_Finalize();

    std::cout << "\nFormatting completed successfully" << "\n\n";

    return 0;
}


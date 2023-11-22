
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
#include "../headers/histo.h"
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
    int bin_width  = 1.0;             //Bin width
    int world_size = 0;               //Size of the mpi world
    int world_rank = 0;               //Rank in the mpi world
    int b_min      = 0;               //Did the user specify a minimum range for the histogram?
    int b_max      = 0;               //Did the user specify a maximum range for the histogram?
    int min        = 0;               //minimum value of histogram
    int max        = 0;               //maximum value of histogram
    double avg     = 0.0;             //Average number of bound lipids 
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
    string program_name = "Binding Lipids";

    print_credits(argc,argv,program_name);

    string program_description = "Binding Lipids is an analysis tool used for reading in a binding events file and counting the number of bound lipids at each time point.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments_mpi(argc,argv,world_rank,program_description);
    add_argument_mpi_s(argc,argv,"-d"        , binding_events_file_name,  "Input binding events file (be)"              , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o"        , out_file_name,             "Output data file with histogram (dat)"       , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-bin"      , &bin_width,                "Bin wdith (number of lipids) "               , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-min"      , &min,                      "Minimum value of histogram (num lipids)"     , world_rank, cl_tags, &b_min,       0);
    add_argument_mpi_i(argc,argv,"-max"      , &max,                      "Maximum value of histogram (num lipids)"     , world_rank, cl_tags, &b_max,       0);
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
    // Check if a range was specified properly for the histogram                                                 //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(b_min != b_max) //only single boundary was specified
    {
        if(world_rank == 0)
        {
            printf("Only a single boundary was specified for the histogram. Please specify either both or no boundaries. \n");
        }
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Make vector to hold lipid counts                                                                          //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    iv1d data(0,0);        //store how many lipid are bound for each frame

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read in binding events                                                                                    //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Binding_events events;

    int result = events.get_binding_events_bin(binding_events_file_name);

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
        // count the number of bound lipids for each frame                                                           //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        for(i=0; i<events.ef_frames; i++) //loop over frames
        {
            int count = 0;

            for(j=0; j<events.num_lipids; j++) //loop over lipids
            {
                if(events.bound_time_line[i][j] == 1)
                {
                    count  = count + 1;
                }
            }
            data.push_back(count);
            avg = avg + count; 
        }
        avg = avg/(double)data.size(); 

        printf("Average number of bound lipids: %f \n",avg);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Make histogram                                                                                            //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        Histogram_i histo;
	if(b_min == 1 && b_max == 1)
        {
            histo.set_range(min,max);
        }
        histo.bin_data(data,bin_width);
        histo.write_histo(out_file_name,"number of lipids");
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


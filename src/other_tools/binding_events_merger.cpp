
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
    string binding_events_file_name_1;        //Name of input file 1
    string binding_events_file_name_2;        //Name of input file 2
    string out_file_name;                     //Name of the output file
    int i                    = 0;             //Standard variable used in loops
    int sort                 = 0;             //Sort binding events by bind_i
    int compound_lipid_count = 0;             //Add lipid count for both files
    int lipid_nr_offset      = 0;             //Offset the second events file lipid nr by this amount
    int do_time_line         = 0;             //Get final events from timeline? 
    int world_size           = 0;             //Size of the mpi world
    int world_rank           = 0;             //Rank in the mpi world

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
    string program_name = "Binding Events Merger";

    print_credits(argc,argv,program_name);

    string program_description = "Binding Events Merger is an analysis tool used for reading in 2 binding events files and merging them into a single file.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments_mpi(argc,argv,world_rank,program_description);
    add_argument_mpi_s(argc,argv,"-d1"       , binding_events_file_name_1,"Input binding events file 1 (be)"                               , world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-d2"       , binding_events_file_name_2,"Input binding events file 2 (be)"                               , world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-sort"     , &sort,                     "Sort the binding events by initial binding time (0:no, 1:yes)"  , world_rank, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-o"        , out_file_name,             "Output data file with merged binding events (be)"               , world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-offset"   , &lipid_nr_offset,          "Offset lipid nr in -d2 by this much"                            , world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-compnd"   , &compound_lipid_count,     "Compound lipid count? (0:no, 1:yes)"                            , world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-tl"       , &do_time_line,             "Get events from time line? (0:no, 1:yes)"                       , world_rank, nullptr,      0);
    conclude_input_arguments_mpi(argc,argv,world_rank,program_name);

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
    // Read in binding events                                                                                    //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Binding_events events;
    int result = events.get_binding_events(binding_events_file_name_1);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read second binding events file                                                                           //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    result = events.add_binding_events(binding_events_file_name_2,compound_lipid_count,lipid_nr_offset);

    if(events.res_nr.size() > 0)
    {
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Get events from timeline                                                                                  //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(do_time_line == 1)
        {
            events.get_binding_timeline();
            events.binding_events_from_timeline();
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Organize the data                                                                                         //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(sort == 1)
        {
            events.organize_events(0);
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Write data to output file                                                                                 //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        events.write_binding_events(out_file_name);
    }

    //relinquish the mpi environment
    MPI_Finalize();

    printf("\n");
    std::cout << "\nFormatting completed successfully" << "\n\n";

    return 0;
}


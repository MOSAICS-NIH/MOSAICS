
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
#include "headers/index.h"                             

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
    string card_file_name;            //Name of the selection card with resid's to be excluded
    int i          = 0;               //General variable used in loops
    int j          = 0;               //General variable used in loops
    int k          = 0;               //General variable used in loops
    int world_size = 0;               //Size of the mpi world
    int world_rank = 0;               //Rank in the mpi world
    int b_x        = 0;               //Was a grid point in x-direction specified
    int b_y        = 0;               //Was a grid point in y-direction specified
    int x          = 0;               //Grid point in x-direction
    int y          = 0;               //Grid point in y-direction
    int result     = 0;               //Tells if the binding events file was read sucressfully
    int threshold  = 0;               //Cutoff for mending fragmented binding events
    int b_card     = 0;               //Was a selection card provided or not?
    double cutoff  = 0;               //Exclude data with cutoff less than this
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
    string program_name = "Binding Extractor";

    print_credits(argc,argv,program_name);

    string program_description = "Binding Extractpr is an analysis tool used for reading in a binding events file and removing some of the events.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments_mpi(argc,argv,world_rank,program_description);
    add_argument_mpi_s(argc,argv,"-d"        , binding_events_file_name,  "Binding events file (be)"                                     , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o"        , out_file_name,             "Output data file with binding events (be)"                    , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-x"        , &x,                        "Grid point in x-direction"                                    , world_rank, cl_tags, &b_x,         0);
    add_argument_mpi_i(argc,argv,"-y"        , &y,                        "Grid point in y-direction"                                    , world_rank, cl_tags, &b_y,         0);
    add_argument_mpi_d(argc,argv,"-cutoff"   , &cutoff,                   "Exclude data with a dwell time smaller than this (ps)"        , world_rank, cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-repair"   , &threshold,                "Maximum allowed size (frames) for mending fragmented events"  , world_rank, cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-crd"      , card_file_name,            "Selection card (crd) with resid's to be excluded"             , world_rank, cl_tags, &b_card,      0);
    conclude_input_arguments_mpi(argc,argv,world_rank,program_name,cl_tags);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension_mpi(world_rank,"-d",binding_events_file_name,".be");
    check_extension_mpi(world_rank,"-o",out_file_name,".be");
    if(b_card == 1)
    {
        check_extension_mpi(world_rank,"-crd",card_file_name,".crd");
    } 

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read the selection card                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Index param;
    if(b_card == 1)
    {
        param.get_index(card_file_name);
    }

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
        // sort events by dwell time (largest first)                                                                 //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        events.organize_events(1);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // remove events with dwell time shorter than cutoff                                                         //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(cutoff > 0)
        {
            for(i=events.lipid_nr.size()-1; i>=0; i--) //loop over binding events
            {
                if((double)events.dwell_t[i]*events.ef_dt < cutoff)
                {
                    events.dwell_t.pop_back();
                    events.lipid_nr.pop_back();
                    events.bind_i.pop_back();
                    events.bind_f.pop_back();
                    events.res_nr.pop_back();
                    events.res_name.pop_back();
                }
            }
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Remove unwanted binding events                                                                            //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        events.get_binding_timeline();              //make a timeline        
        events.suppress_timeline_noise(threshold);  //mend fragmented events

        //exclude lipids from list
        if(b_card == 1)
        {
            for(i=0; i<events.num_lipids; i++) //loop over lipids
            {
                for(j=0; j<param.index_s.size(); j++) //loop over list of residue ids
                {
                    if(events.time_line_res_nr[i] == param.index_i[j]) //resid is correct
                    {
                        for(k=0; k<events.ef_frames; k++) //loop over frames
                        {
                            events.bound_time_line[k][i] = 0;
                        }
                    }
                }
            }
        }
	events.binding_events_from_timeline();      //generate binding events from mended timeline

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Write binding events to file                                                                              //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        events.write_binding_events_bin(out_file_name);

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


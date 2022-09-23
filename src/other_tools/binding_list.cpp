
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

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// The main function performing analysis                                                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, const char * argv[]) 
{
    //Here we define some variables used throughout
    string binding_events_file_name; //name of the binding events file
    int i           = 0;                 //Standard variable used in loops
    int mode        = 0;                 //Controls the mode for sorting output data
    int dt_mode     = 0;                 //controls whether time is in frames or ns
    int world_size  = 0;                 //Size of the mpi world
    int world_rank  = 0;                 //Rank in the mpi world
    double dt       = 0.0;               //dt used for final output

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
    // Set program name/description and print info                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    string program_name = "Binding List";

    print_credits(argc,argv,program_name);

    string program_description = "Binding List is an analysis tool used for organizing and displaying contents of a binding events file.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments_mpi(argc,argv,world_rank,program_description);
    add_argument_mpi_s(argc,argv,"-d"     , binding_events_file_name,   "Input binding events file (be)"                 , world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-mode"  , &mode,                      "Order events by (0:bind_i,1:dwell_t,2:repeats)" , world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-t"     , &dt_mode,                   "Report time in (0:ns,1:frames)"                 , world_rank, nullptr,      1);
    conclude_input_arguments_mpi(argc,argv,world_rank,program_name);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension_mpi(world_rank,"-d",binding_events_file_name,".be");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read in binding events                                                                                    //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Binding_events events;

    int result = events.get_binding_events(binding_events_file_name);

    if(result == 1) //binding events file exists
    {
        printf("x_i %d y_i %d ef_dt(ps) %f ef_frames %d num_lipids %d num_g_x %d num_g_y %d APS(nm^2) %f \n",events.x_i,events.y_i,events.ef_dt,events.ef_frames,events.num_lipids,events.num_g_x,events.num_g_y,events.APS);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Set dt                                                                                                    //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(dt_mode == 0) //time
        {
            dt = events.ef_dt; //convert to ps
        }
        else if(dt_mode == 1) //frames
        {
            dt = 1; 
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Characterize repeated binding                                                                             //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        events.reduce_list();

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Organize the binding events                                                                               //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        events.organize_events(mode);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Print info about the binding events                                                                       //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(mode == 0 || mode == 1) //print by initial binding time or dwell time
        {
            printf("\n");
            if(dt_mode == 1) //frames
            {
                printf(" %10s %10s %10s %10s %15s %15s %20s \n","rank","lipid_nr","res_nr","res_name","bind_i(frame)","bind_f(frame)","dwell time(frames)");
            }
            else if(dt_mode == 0) //ps
            {
                printf(" %10s %10s %10s %10s %15s %15s %20s \n","rank","lipid_nr","res_nr","res_name","bind_i(ps)","bind_f(ps)","dwell time(ps)");
            }
            printf(" %10s-%10s-%10s-%10s-%15s-%15s-%20s \n","----------","----------","----------","----------","---------------","---------------","--------------------");
            for(i=0; i<events.dwell_t.size(); i++)
            {
                printf(" %10d %10d %10d %10s %15.1f %15.1f %20.1f \n",i,events.lipid_nr[i],events.res_nr[i],events.res_name[i].c_str(),(double)events.bind_i[i]*dt,(double)events.bind_f[i]*dt,(double)events.dwell_t[i]*dt);
            }
        }
        else if(mode == 2) //print by repeated visits
        {
            printf("\n");
            if(dt_mode == 1) //frames
            {
                printf(" %10s %10s %10s %10s %20s %20s \n","rank","res_nr","res_name","visits","dwell_time(frames)","off_time(frames)");
            }
            else if(dt_mode == 0) //ps
            {
                printf(" %10s %10s %10s %10s %20s %20s \n","rank","res_nr","res_name","visits","dwell_time(ps)","off_time(ps)");
            }
            printf(" %10s-%10s-%10s-%10s-%10s-%10s \n","----------","----------","----------","----------","--------------------","--------------------");
            for(i=0; i<events.repeats_res_nr.size(); i++)
            {
                printf(" %10d %10d %10s %10d %20.1f %20.1f \n",i,events.repeats_res_nr[i],events.repeats_res_name[i].c_str(),events.repeats[i]+1,events.repeats_avg_dwell_time[i]*dt,(double)events.repeats_avg_off_time[i]*dt);
            }
        }
    }
    else //binding events file does not exist 
    {
        printf("Could not find binding events file \n");
    }

    std::cout << "\nFormatting completed successfully" << "\n\n";

    return 0;
}


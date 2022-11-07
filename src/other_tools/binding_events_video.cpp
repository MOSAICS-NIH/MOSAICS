
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
#include "../headers/binding_events.h"
#include "../headers/grid_lt.h"
#include "../headers/fit.h"
#include "../headers/command_line_args_mpi.h"
#include "../headers/array.h"
#include "../headers/file_naming_mpi.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Main function which executes all other functions                                                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, const char * argv[]) 
{
    //Here we define some variables used throughout
    string in_file_name;                //Name of input file 
    string base_file_name_i;            //Base file name of input binding events files
    string base_file_name_o;            //Base file name of output files
    string be_file_name;                //Name of the binding_events input file
    string out_file_name;               //Name of the output files with the grid point selection and the big mask
    string rho_file_name;               //Name of the rho file
    int i                 = 0;          //General variable used in loops
    int j                 = 0;          //General variable used in loops
    int k                 = 0;          //General variable used in loops
    int l                 = 0;          //General variable used in loops
    int world_size        = 0;          //Size of the mpi world
    int world_rank        = 0;          //Rank in the mpi world
    int odf               = 0;          //Rho data file format
    int stride            = 0;          //Skip stride frames
    int begin             = 0;          //Start on this frame
    int end               = 0;          //End on this frame
    int b_end             = 0;          //Did the user provide an end frame?
    int b_select          = 0;          //Was a binding events file provided for highlighting lipids?
    int b_rho             = 0;          //Did the user specify a rho input file name?
    double cell_size      = 1;          //Distance between grid points
    double dt             = 0;          //Time step used for converting frames to time. set equal to ef_dt
    double cutoff         = 0;          //Cutoff for excluding data
    double avg_rho        = 0;          //The average lipid density over the grid

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
    // Set program name and print info                                                                           //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    string program_name = "Binding Events Video";
    if(world_rank == 0)
    {
        print_credits(argc,argv,program_name);
    }

    string program_description = "Binding Events Video is an analysis tool that lets the user produce a video of the lipid dynamics from binding events data produced by 2D Kinetics. This tool has the option of highlighting select lipids as specified via a binding events file.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments_mpi(argc,argv,world_rank,program_description);
    add_argument_mpi_s(argc,argv,"-d"       , base_file_name_i,           "Base filename for input binding events files"                             , world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-be"      , be_file_name,               "Binding events file for highlighting select lipids (be)"                  , world_rank, &b_select,    0);
    add_argument_mpi_s(argc,argv,"-o"       , base_file_name_o,           "Base filename for output data files with noise filtered voronoi diagrams" , world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-stride"  , &stride,                    "Skip stride frames"                                                       , world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-b"       , &begin,                     "Start at this frame"                                                      , world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e"       , &end,                       "End on this frame"                                                        , world_rank, &b_end,       0);
    add_argument_mpi_s(argc,argv,"-rho"     , rho_file_name,              "Input data file with sample count (dat)"                                  , world_rank, &b_rho,       0);
    add_argument_mpi_d(argc,argv,"-cutoff"  , &cutoff,                    "Cutoff for excluding grid data (chi)"                                     , world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-odf"     , &odf,                       "Data file format for sample count (0:matrix 1:vector)"                    , world_rank, nullptr,      0);
    conclude_input_arguments_mpi(argc,argv,world_rank,program_name);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(b_select == 1)
    {
        check_extension_mpi(world_rank,"-be",be_file_name,".be");
    }
    if(b_rho == 1)
    {
        check_extension_mpi(world_rank,"-rho",rho_file_name,".dat");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Create Grid for excluding insignificant data                                                              //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Grid_lt rho;

    if(b_rho == 1)
    {
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Set the grid format                                                                                       //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        rho.set_format(odf);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Read in grid data                                                                                         //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        rho.get_grid(rho_file_name);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Compute average rho                                                                                       //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        for(i=0; i<rho.size_x(); i++) //loop over x
        {
            for(j=0; j<rho.size_y(); j++) //loop over y
            {
                avg_rho = avg_rho + rho.grid[i][j][2][0];
            }
        }
        avg_rho = avg_rho/(rho.size_x()*rho.size_y());

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Exclude insignificant data                                                                                //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        rho.exclude_grid_data(cutoff,avg_rho,rho.grid);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read binding events file 0_0 to get the header info                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Binding_events events_ref;
    in_file_name = base_file_name_i + "_" + to_string(0) + "_" + to_string(0) + ".be";
    int result    = events_ref.get_binding_events(in_file_name);

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
        events_ref.get_binding_timeline();
    }

    cell_size = sqrt(events_ref.APS);
    dt        = events_ref.ef_dt/1000;

    if(b_end == 0)
    {
        end = events_ref.ef_frames; 
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Distribute the workload across the cores                                                                  //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int my_num_g_x = count_workload(world_size,world_rank,events_ref.num_g_x);

    //create array to hold each mpi processes num_g_x; Used for communication
    int world_num_g_x_ary[world_size];
    MPI_Allgather(&my_num_g_x, 1,MPI_INT,world_num_g_x_ary, 1, MPI_INT, MPI_COMM_WORLD );

    //allocate memory for world_num_g_x and copy data from the array
    iv1d world_num_g_x(world_size,0);
    for(i=0; i<world_size; i++)
    {
        world_num_g_x[i] = world_num_g_x_ary[i];
    }

    //print stats for distributing the num_g_x and distribute the num_g_x to each core
    int my_xi = 0;
    int my_xf = 0;
    int world_xi[world_size];
    int world_xf[world_size];
    get_workload(&my_xi,&my_xf,world_rank,world_num_g_x,events_ref.num_g_x,world_xi,world_xf);
    print_workload_stats(world_rank,world_xi,world_xf,world_num_g_x,world_size,"num_g_x","init","fin");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set the effective number of frames taking into account the stride                                         //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int ef_frames = (int)events_ref.ef_frames/stride + 1;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read in binding events file of interest                                                                   //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Binding_events target_events;
    if(b_select == 1)
    {
        result = target_events.get_binding_events(be_file_name);

        if(result == 1)
        {
            //get binding timeline for the target be file
            target_events.get_binding_timeline();
        }
        else
        {
            printf("Could not find specified binding events file. Please check the name and try again. \n");
            MPI_Finalize();
            return 0;
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check compatability of grids                                                                              //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int fail = 0;

    if(b_rho == 1)
    {
        if(events_ref.num_g_x != rho.size_x() || events_ref.num_g_y != rho.size_y())
        {
            fail = 1;
        }
    }
    if(b_select == 1)
    {
        if(events_ref.num_g_x != target_events.num_g_x || events_ref.num_g_y != target_events.num_g_y)
        {
            fail = 1;
        }
    }

    if(fail == 1)
    {
        if(world_rank == 0)
        {
            printf("Grid data is incompatible. \n");
            printf("%20f: size_x %10d size_y %10d \n","binding_events -d",events_ref.num_g_x,events_ref.num_g_y);
            printf("%20f: size_x %10d size_y %10d \n","rho data -rho",rho.size_x(),rho.size_y());
            printf("%20f: size_x %10d size_y %10d \n","selection -be",target_events.num_g_x,target_events.num_g_y);
        }
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Get lipid blobs                                                                                           //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    events_ref.get_blobs(base_file_name_i,my_xi,my_xf,my_num_g_x,stride,ef_frames,world_rank);

    MPI_Barrier(MPI_COMM_WORLD);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Collect data and write grids to output                                                                    //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(i=0; i<events_ref.ef_frames; i+=stride)
    {
        if(i >= begin && i <= end) //check start and end frame condition
        {
            if(world_rank == 0)
            {
                printf("Writing grid data. Working on frame %d \n",i);
                fflush(stdin);
            }

            int this_frame = (int)(i/stride);

            events_ref.get_blobs_frame(this_frame,world_size,world_rank);

            if(world_rank == 0)
            {
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // Exclude insignificant data                                                                                //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                if(b_rho == 1)
                {
                    for(j=0; j<rho.size_x(); j++) //loop over x
                    {
                        for(k=0; k<rho.size_y(); k++) //loop over y
                        {
                            if(rho.grid[j][k][2][1] == 1)
                            {
                                events_ref.blob_nan_frame[k][j] = 1;
                            }
                        }
                    }
                }

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // Highlight target lipids                                                                                   //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                if(b_select == 1) //was a binding events file supplied for highlighting lipids?
                {
                    for(j=0; j<events_ref.num_g_x; j++) //loop over x
                    {
                        for(k=0; k<events_ref.num_g_y; k++) //loop over y
                        {
                            if(events_ref.blob_nan_frame[k][j] == 0) //check that grid point is a lipid
                            {
                                if(target_events.bound_time_line[this_frame*stride][events_ref.blob_frame[k][j]] == 1) //check if lipid is a target lipid
                                {
                                    events_ref.blob_frame[k][j] = -1*events_ref.blob_frame[k][j];
                                }
                            }
                        }
                    }
                }

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // Write grid to file                                                                                        //
                //                                                                                                           //
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
                out_file_name = base_file_name_o + "_" + to_string(i) + ".dat";
                events_ref.write_blobs_frame(out_file_name);
            }

            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Print final output and end program                                                                        //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(world_rank == 0)
    {
        std::cout << "\nFormatting completed successfully" << "\n\n";
    }

    MPI_Finalize();

    return 0;
}


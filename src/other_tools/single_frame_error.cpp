
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
#include "../headers/grid_lt.h"
#include "../headers/common_routines_mpi.h"
#include "../headers/common_routines.h"
#include "../headers/command_line_args_mpi.h"
#include "../headers/array.h"
#include "../headers/file_naming_mpi.h"

int main(int argc, const char * argv[]) 
{
    //Here we define some variables used throughout
    string in_file_name;               //Name of input file with mean 
    string base_file_name;             //Base name of input files with single frame 
    string frame_file_name;            //Name of current single frame file
    string out_file_name;              //Name of output file with error
    string mask_file_name;             //Name of file with mask data
    int i                 = 0;         //General variable used in loops
    int j                 = 0;         //General variable used in loops
    int k                 = 0;         //General variable used in loops
    int current_frame     = 0;         //Current single frame working on
    int odf               = 0;         //Data file format
    int num_frames        = 0;         //How many single frame grids?
    int winner            = -1;        //Frame with the smallest error                   
    int world_rank        = 0;         //Rank of the core
    int world_size        = 0;         //How many cores in world
    double nan            = 0.0;       //Value added to grid when NaN is encountered
    double smallest_error = 999999999; //Smallest error of all single frames
    sv1d cl_tags;                      //Holds a list of command line tags for the program

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
    string program_name = "Single Frame Error";

    if(world_rank == 0)
    {
        print_credits(argc,argv,program_name);
    }

    string program_description = "Single Frame Error is an analysis tool that compares single frame grid data to the spatially resolved time average and finds the frame that most closesly matches this average.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments_mpi(argc,argv,world_rank,program_description);
    add_argument_mpi_s(argc,argv,"-d"       , in_file_name,               "Spatially resolved time average observable of interest (dat)"      , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-mask"    , mask_file_name,             "Mask file with selection of grid region of interest (dat)"         , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-base"    , base_file_name,             "Base filename for the single frame grid data"                      , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o"       , out_file_name,              "Output data file with the error for each trajectory frame (dat)"   , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-frames"  , &num_frames,                "How many single frame grids are there? "                           , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-odf"     , &odf,                       "Data file format (0:matrix 1:vector) "                             , world_rank, cl_tags, nullptr,      1);
    conclude_input_arguments_mpi(argc,argv,world_rank,program_name,cl_tags);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension_mpi(world_rank,"-d",in_file_name,".dat");
    check_extension_mpi(world_rank,"-mask",mask_file_name,".dat");
    check_extension_mpi(world_rank,"-o",out_file_name,".dat");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Distribute the workload across the cores                                                                  //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int my_frames = count_workload(world_size,world_rank,num_frames);

    //create array to hold each mpi processes frames; Used for communication
    int world_frames_ary[world_size];
    MPI_Allgather(&my_frames, 1,MPI_INT,world_frames_ary, 1, MPI_INT, MPI_COMM_WORLD );

    //allocate memory for world_frames and copy data from the array
    iv1d world_frames(world_size,0);
    for(i=0; i<world_size; i++)
    {
        world_frames[i] = world_frames_ary[i];
    }

    //print stats for distributing the frames and distribute the frames to each core
    int my_i = -1;
    int my_f = -1;
    int world_i[world_size];
    int world_f[world_size];
    get_workload(&my_i,&my_f,world_rank,world_frames,num_frames,world_i,world_f);
    print_workload_stats(world_rank,world_i,world_f,world_frames,world_size,"frames","init","fin");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Create arrays for holding local/global error                                                              //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double local_error[my_frames];
    double global_error[num_frames];

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Create Grids                                                                                              //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Grid_lt data;
    Grid_lt mask;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set the grid format                                                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    data.set_format(odf);
    mask.set_format(odf);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read in grid data                                                                                         //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    data.get_grid(in_file_name);
    mask.get_grid(mask_file_name);

    //check that grid files are compatible
    int check = comp_grid(data,mask);

    if(check == 1)
    {
        if(world_rank == 0)
        {
            data.print_dim(1);
            mask.print_dim(0);
            printf("\n");
        }
    }
    else
    {
        if(world_rank == 0)
        {
            printf("Input files are not compatible. \n");

            data.print_dim(1);
            mask.print_dim(0);
        }
        MPI_Finalize();
        return 0;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Loop over single frame data                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(current_frame=my_i; current_frame<=my_f; current_frame++)
    {
        int effective_frame        = current_frame - my_i;
        int count                  = 0;
        double current_frame_error = 0;

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Set the file name for the current frame                                                                   //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        string frame_file_name = base_file_name + to_string(current_frame) + ".dat";

        if(world_rank == 0 && current_frame%10 == 0)
        {
            printf("Working on frame %s \n",frame_file_name.c_str());
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Read in single frame data                                                                                 //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        Grid_lt frame_data;
        frame_data.set_format(odf);
        frame_data.get_grid(frame_file_name);

        //check that grid files are compatible
        int check_1 = comp_grid(frame_data,data);

        if(check_1 == 0)
        {
            if(world_rank == 0)
            {
                printf("Input files are not compatible. \n");

                frame_data.print_dim(1);
            }
            MPI_Finalize();
            return 0;
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Compute the error                                                                                         //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        for(j=0; j<data.size_x(); j++) //loop over x
        {
            for(k=0; k<data.size_y(); k++) //loop over y
            {
                if(data.grid[j][k][2][1] == 0 && frame_data.grid[j][k][2][1] == 0 && mask.grid[j][k][2][0] == 1) //data present for both and mask passes
                {
                    double weight       = mask.grid[j][k][2][0];
                    current_frame_error = current_frame_error + weight*pow(data.grid[j][k][2][0]-frame_data.grid[j][k][2][0],2);
                    count++;
                }
            }
        }

        local_error[effective_frame] = sqrt(current_frame_error/(double)count);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Collect local errors                                                                                      //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    collect_local_values(world_size,world_rank,global_error,local_error,world_frames,"error");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the error                                                                                         //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(world_rank == 0)
    {
        FILE *out_file = fopen(out_file_name.c_str(), "w");
        if(out_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",out_file_name.c_str());
            exit(EXIT_SUCCESS);
        }
        else 
        {
            fprintf(out_file,"# %20s %10s \n","trajectory frame","error");
            for(i=0; i<num_frames; i++)
            {
                if(global_error[i] < smallest_error)
                {
                    smallest_error = global_error[i];
                    winner = i;
                }
                fprintf(out_file,"  %20d %10f \n",i,global_error[i]);
            }
            fclose(out_file);
        }

        printf("Frame %6d has the smallest error %f \n",winner,smallest_error);
    }

    if(world_rank == 0)
    {
        std::cout << "\nFormatting completed successfully" << "\n\n";
    }

    MPI_Finalize();

    return 0;
}


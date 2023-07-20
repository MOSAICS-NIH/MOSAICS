
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
#include "../headers/grid_lt.h"
#include "../headers/common_routines_mpi.h"
#include "../headers/common_routines.h"
#include "../headers/command_line_args_mpi.h"
#include "../headers/array.h"
#include "../headers/file_naming_mpi.h"

int main(int argc, const char * argv[]) 
{
    //Here we define some variables used throughout
    FILE *out_file;               //File for writing data
    string in_file_name_mask;     //Name of input data file
    string in_file_name_data;     //Name of input mask file
    string out_file_name;         //Name of output file
    string out_file_name_mask;    //Name of the output file with masking data
    int i                 = 0;    //General variable used in loops
    int j                 = 0;    //General variable used in loops
    int odf               = 0;    //Data file format
    int iterations        = 1;    //How many times to move selection
    int target_x          = 0;    //The target grid point x when making a rectangular selection
    int target_y          = 0;    //The target grid point y when making a rectangular selection
    int range_x           = 0;    //The half width of x in the rectangular selection
    int range_y           = 0;    //The half width of y in the rectangular selection
    int invert            = 0;    //Invert rectangular selection (select everything outside rectangle)
    int current_iteration = 0;    //The iteration currently being worked on
    int world_size        = 0;    //Size of the mpi world
    int world_rank        = 0;    //Rank in the mpi world
    int cumulative        = 0;    //Include all grid points up to dist?
    int b_highlight_prot  = 0;    //Highlight prot in mask files?
    double cell_size      = 1;    //Distance between grid points
    double res            = 1;    //How much does the grid selection move with each iteration
    double range          = 1;    //The half width of the grid point selection
    double nan            = 0.0;  //Value added to grid when NaN is encountered
    double APS            = 0;    //Area per square used for grid in analysis
    sv1d cl_tags;                 //Holds a list of command line tags for the program

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
    string program_name = "Grid Distance Projection";

    if(world_rank == 0)
    {
        print_credits(argc,argv,program_name);
    }

    string program_description = "Grid Distance Projection is an analysis tool that lets the user project grid data as a function of distance from the protein surface.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments_mpi(argc,argv,world_rank,program_description);
    add_argument_mpi_s(argc,argv,"-d"     , in_file_name_data,          "Input data file with the spatially resolved observable of interest (dat)"                , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-mask"  , in_file_name_mask,          "Input protein mask file (dat)"                                                           , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o"     , out_file_name,              "Output data file with observable as a function of distance to the protein surface (dat)" , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-x"     , &target_x,                  "Rectangle center x (grid point)"                                                         , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-y"     , &target_y,                  "Rectangle center y (grid point)"                                                         , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-rx"    , &range_x,                   "Rectangle half width x (grid points)  "                                                  , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-ry"    , &range_y,                   "Rectangle half width y (grid points)  "                                                  , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-invert", &invert,                    "Invert rectangular selection? (0:no 1:yes)"                                              , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-iter"  , &iterations,                "How many iterations to perform"                                                          , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-res"   , &res,                       "Distance moved between each iteration (nm)"                                              , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-range" , &range,                     "Half width of the grid selection shell (nm)"                                             , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-APS"   , &APS,                       "Area per grid square (nm^2)"                                                             , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-odf"   , &odf,                       "Data file format (0:matrix 1:vector)"                                                    , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-inc_p" , &b_highlight_prot,          "Highlight protein in mask files? (0:no 1:yes)"                                           , world_rank, cl_tags, nullptr,      0);
    conclude_input_arguments_mpi(argc,argv,world_rank,program_name,cl_tags);

    cell_size = sqrt(APS);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension_mpi(world_rank,"-d",in_file_name_data,".dat");
    check_extension_mpi(world_rank,"-mask",in_file_name_mask,".dat");
    check_extension_mpi(world_rank,"-o",out_file_name,".dat");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Distribute the workload across the cores                                                                  //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int my_iterations = count_workload(world_size,world_rank,iterations);

    //create array to hold each mpi processes iterations; Used for communication
    int world_iterations_ary[world_size];
    MPI_Allgather(&my_iterations, 1,MPI_INT,world_iterations_ary, 1, MPI_INT, MPI_COMM_WORLD );

    //allocate memory for world_iterations and copy data from the array
    iv1d world_iterations(world_size,0);
    for(i=0; i<world_size; i++)
    {
        world_iterations[i] = world_iterations_ary[i];
    }

    //print stats for distributing the iterations and distribute the iterations to each core
    int my_i = -1;
    int my_f = -1;
    int world_i[world_size];
    int world_f[world_size];
    get_workload(&my_i,&my_f,world_rank,world_iterations,iterations,world_i,world_f);
    print_workload_stats(world_rank,world_i,world_f,world_iterations,world_size,"iterations","init","fin");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Create arrays for holding averages                                                                        //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double local_averages[my_iterations];
    double global_averages[iterations];

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Create Grids                                                                                              //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Grid_lt data;
    Grid_lt init_mask;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set the grid format                                                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    data.set_format(odf);
    init_mask.set_format(odf);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read in grid data                                                                                         //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    data.get_grid(in_file_name_data);
    init_mask.get_grid(in_file_name_mask);

    //check that grid files are compatible
    int check = comp_grid(data,init_mask);

    if(check == 1)
    {
        if(world_rank == 0)
        {
            data.print_dim(1);
            init_mask.print_dim(0);
        }
    }
    else
    {
        if(world_rank == 0)
        {
            printf("Input files are not compatible. \n");

            data.print_dim(1);
            init_mask.print_dim(0);
        }
        MPI_Finalize();
        return 0;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Move the grid selection with each iteration                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(current_iteration=my_i; current_iteration<=my_f; current_iteration++) //loop over the iterations
    {
        //print update so the user knows some progress has been made.
        if(world_rank == 0)
        {
            printf("Working on iteration %5d \n",current_iteration);
        }

        int effective_iteration = current_iteration - my_i;                              //this is the index for storing the average

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Make a mask for the current iteration                                                                     //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        Grid_lt mask;
        mask.set_format(odf);
        mask.get_grid(in_file_name_mask);
        mask.init_grid(0.0);              //initialize grid to zero
        mask.set_nan(0);                  //initialize nan tags

        string tag = "_" + to_string(current_iteration) + "_mask";
        out_file_name_mask = add_tag(out_file_name,tag);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Compute the grid selection for the current iteration                                                      //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        mask.distance_projection(invert,target_x,target_y,range_x,range_y,cell_size,current_iteration,range,res,init_mask.grid,cumulative);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Highlight protein in the mask                                                                             //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(b_highlight_prot == 1)
        {
            for(i=0; i<data.size_x(); i++) //loop over x
            {
                for(j=0; j<data.size_y(); j++) //loop over y
                {
                    if(init_mask.grid[i][j][2][0] == 1) //grid points belongs to the protein
                    {
                        mask.grid[i][j][2][0] = -1;
                    }
                }
            }
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Write grid to output file                                                                                 //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        mask.write_grid(out_file_name_mask);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Compute average over mask                                                                                 //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double avg = 0;     //average over the mask
        int count = 0;      //How many grid points make the mask
        for(i=0; i<data.size_x(); i++) //loop over x
        {
            for(j=0; j<data.size_y(); j++) //loop over y
            {
                if(mask.grid[i][j][2][0] == 1) //grid points belongs to the protein
                {
                    avg = avg + data.grid[i][j][2][0];
                    count++;
                }
            }
        }
        avg = avg/(double)count;

        local_averages[effective_iteration] = avg;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Collect averages                                                                                          //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    collect_local_values(world_size,world_rank,global_averages,local_averages,world_iterations,"averages");

    if(world_rank == 0)
    {
        out_file      = fopen(out_file_name.c_str(), "w");
        fprintf(out_file,"# %15s %25s \n","distance (nm)","Observable of interest");
        for(i=0; i<iterations; i++)
        {
           fprintf(out_file,"  %15f %25f \n",(double)i*res,global_averages[i]);
        }
        fclose(out_file);
    }

    if(world_rank == 0)
    {
        std::cout << "\nFormatting completed successfully" << "\n\n";
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}


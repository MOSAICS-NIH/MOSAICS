
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
    FILE *out_file;                 //File for writing data
    string in_file_name_mask;       //Name of input data file
    string out_file_name;           //Name of output file
    string out_file_name_mask;      //Name of the output file with masking data
    string rho_file_name_1;         //Name of lipid density file 1
    string rho_file_name_2;         //Name of lipid density file 2
    int i                   = 0;    //General variable used in loops
    int j                   = 0;    //General variable used in loops
    int odf                 = 0;    //Data file format
    int iterations          = 1;    //How many times to move selection
    int target_x            = 0;    //The target grid point x when making a rectangular selection
    int target_y            = 0;    //The target grid point y when making a rectangular selection
    int range_x             = 0;    //The half width of x in the rectangular selection
    int range_y             = 0;    //The half width of y in the rectangular selection
    int invert              = 0;    //Invert rectangular selection (select everything outside rectangle)
    int current_iteration   = 0;    //The iteration currently being worked on
    int world_size          = 0;    //Size of the mpi world
    int world_rank          = 0;    //Rank in the mpi world
    int cumulative          = 0;    //Include all grid points up to dist?
    double cell_size        = 1;    //Distance between grid points
    double res              = 1;    //How much does the grid selection move with each iteration
    double range            = 1;    //The half width of the grid point selection
    double APS              = 0;    //Area per square used for grid in analysis
    double ratio_bulk       = 0;    //Ratio in the bulk

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
    string program_name = "2d Enrichment Distance Projection";

    if(world_rank == 0)
    {
        print_credits(argc,argv,program_name);
    }

    string program_description = "2d Enrichment Distance Projection is an analysis tools that lets the user project enrichment data as a function of distance from the protein surface.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments_mpi(argc,argv,world_rank,program_description);
    add_argument_mpi_s(argc,argv,"-rho_A" , rho_file_name_1,            "Input data file with sample count for lipids A (dat)"    , world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-rho_B" , rho_file_name_2,            "Input data file with sample count for lipids B (dat)"    , world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-mask"  , in_file_name_mask,          "Input protein mask file (dat)"                           , world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o"     , out_file_name,              "Output data file with projected enrichment factor (dat)" , world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-x"     , &target_x,                  "Rectangle center x (grid point)"                         , world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-y"     , &target_y,                  "Rectangle center y (grid point)"                         , world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-rx"    , &range_x,                   "Rectangle half width x (grid points)"                    , world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-ry"    , &range_y,                   "Rectangle half width y (grid points)"                    , world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-invert", &invert,                    "Invert rectangular selection? (0:no 1:yes)"              , world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-iter"  , &iterations,                "How many iterations to perform?"                         , world_rank, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-res"   , &res,                       "Distance moved between each iteration (nm)"              , world_rank, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-range" , &range,                     "Half width of the grid selection shell (nm)"             , world_rank, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-APS"   , &APS,                       "Area per grid square (nm^2)"                             , world_rank, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-ratio" , &ratio_bulk,                "Ratio of lip A to lip B in bulk?"                        , world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-odf"   , &odf,                       "Data file format (0:matrix 1:vector)"                    , world_rank, nullptr,      1);
    conclude_input_arguments_mpi(argc,argv,world_rank,program_name);

    cell_size = sqrt(APS);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension_mpi(world_rank,"-rho1",rho_file_name_1,".dat");
    check_extension_mpi(world_rank,"-rho2",rho_file_name_2,".dat");
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
    int my_i = 0;
    int my_f = 0;
    int world_i[world_size];
    int world_f[world_size];
    get_workload(&my_i,&my_f,world_rank,world_iterations,iterations,world_i,world_f);
    print_workload_stats(world_rank,world_i,world_f,world_iterations,world_size,"iterations","init","fin");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Create arrays for holding averages                                                                        //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double density_local_1_local[my_iterations];
    double density_local_2_local[my_iterations];
    double enrichment_local[my_iterations];

    double density_local_1_global[iterations];
    double density_local_2_global[iterations];
    double enrichment_global[iterations];

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Create Grids                                                                                              //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Grid_lt rho_1;
    Grid_lt rho_2;
    Grid_lt init_mask;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set the grid format                                                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    rho_1.set_format(odf);
    rho_2.set_format(odf);
    init_mask.set_format(odf);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read in grid data                                                                                         //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    rho_1.get_grid(rho_file_name_1);
    rho_2.get_grid(rho_file_name_2);
    init_mask.get_grid(in_file_name_mask);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check that grid files are compatible                                                                      //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int check_1 = comp_grid(init_mask,rho_1);
    int check_2 = comp_grid(init_mask,rho_2);

    if(check_1 == 1 && check_2 == 1)
    {
        if(world_rank == 0)
        {
            init_mask.print_dim(1);
            rho_1.print_dim(0);
            rho_2.print_dim(0);
        }
    }
    else
    {
        if(world_rank == 0)
        {
            printf("Input files are not compatible. \n");

            init_mask.print_dim(1);
            rho_1.print_dim(0);
            rho_2.print_dim(0);
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

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Compute the grid selection for the current iteration                                                      //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        mask.distance_projection(invert,target_x,target_y,range_x,range_y,cell_size,current_iteration,range,res,init_mask.grid,cumulative);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Write grid to output file                                                                                 //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        string tag = "_" + to_string(current_iteration) + "_mask";
        out_file_name_mask = add_tag(out_file_name,tag);

        mask.write_grid(out_file_name_mask);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Compute rho_1 sum over mask                                                                               //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double density_local_1 = 0;     //loval density of lipid 1
        for(i=0; i<rho_1.size_y(); i++) //loop over y
        {
            for(j=0; j<rho_1.size_x(); j++) //loop over x
            {
                int rectangle_pass = check_rectangular_pass_vector(invert,target_x,target_y,range_x,range_y,j,i);

                if(rectangle_pass == 1)
                {
                    if(mask.grid[j][i][2][0] == 1) //grid point is in mask
                    {
                        density_local_1 = density_local_1 + rho_1.grid[j][i][2][0];
                    }
                }
            }
        }

        density_local_1_local[effective_iteration] = density_local_1;

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Compute rho_2 sum over mask                                                                               //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double density_local_2 = 0;     //loval density of lipid 1
        for(i=0; i<rho_2.size_y(); i++) //loop over y
        {
            for(j=0; j<rho_2.size_x(); j++) //loop over x
            {
                int rectangle_pass = check_rectangular_pass_vector(invert,target_x,target_y,range_x,range_y,j,i);

                if(rectangle_pass == 1)
                {
                    if(mask.grid[j][i][2][0] == 1) //grid point is in the mask
                    {
                        density_local_2 = density_local_2 + rho_2.grid[j][i][2][0];
                    }
                }
            }
        }

        density_local_2_local[effective_iteration] = density_local_2;

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Compute the percent enrichment                                                                            //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double ratio                          = density_local_1/density_local_2;
        enrichment_local[effective_iteration] = 100*(ratio - ratio_bulk)/ratio_bulk;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Collect data from all cores and write to output file                                                      //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    collect_local_values(world_size,world_rank,density_local_1_global,density_local_1_local,world_iterations,"Sample count A");
    collect_local_values(world_size,world_rank,density_local_2_global,density_local_2_local,world_iterations,"Sample count B");
    collect_local_values(world_size,world_rank,enrichment_global,     enrichment_local,     world_iterations,"Enrichment");

    if(world_rank == 0)
    {
        out_file = fopen(out_file_name.c_str(), "w");
        if(out_file == NULL) //file does not exist
        {
            printf("failure opening %s. Make sure the file exists. \n",out_file_name.c_str());
        }
        else
        {
            fprintf(out_file,"#%15s %20s %20s %10s \n","distance (nm)","sample_count_A","sample_count_B","enrich (%)");

            for(i=0; i<iterations; i++)
            {
                double distance = i*res;
                fprintf(out_file," %15f %20.0f %20.0f %10f \n",distance,density_local_1_global[i],density_local_2_global[i],enrichment_global[i]);
            }
        }
    }

    if(world_rank == 0)
    {
        std::cout << "\nFormatting completed successfully" << "\n\n";
    }

    MPI_Finalize();
    return 0;
}


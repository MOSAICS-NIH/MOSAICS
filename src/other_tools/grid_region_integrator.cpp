
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <sstream>
#include <fstream>

using namespace std;

#include "../headers/multi_dim_vec.h"
#include "../headers/file_reader.h"
#include "../headers/grid_lt.h"
#include "../headers/common_routines.h"
#include "../headers/command_line_args.h"
#include "../headers/array.h"
#include "../headers/file_naming.h"

int main(int argc, const char * argv[]) 
{
    //Here we define some variables used throughout
    string in_file_name;          //Name of input file 
    string mask_file_name;        //Name of mask file
    string list_file_name;        //Name of output file containing values added
    int i                = 0;     //General variable used in loops
    int j                = 0;     //General variable used in loops
    int target_x         = 0;     //The target grid point x
    int target_y         = 0;     //The target grid point y
    int range_x          = 0;     //How many grid points on each side (x) to include
    int range_y          = 0;     //How many grid points on each side (y) to include
    int count            = 0;     //How many points fall in the range
    int odf              = 0;     //Data file format
    int b_mask           = 0;     //Was a masking file provided?
    int invert           = 0;     //Invert rectangular selection?
    int b_list           = 0;     //Write a list of integrated data?
    double avg           = 0;     //Average value reported
    double sum_o_squares = 0;     //Sum of squares for computing standard deviation
    double stdev         = 0;     //The standard deviation
    sv1d cl_tags;                 //Holds a list of command line tags for the program

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set program name and print info                                                                           //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    string program_name = "Grid Region Integrator";

    print_credits(argc,argv,program_name);

    string program_description = "Grid Region Integrator is a tool that lets the user integrate over a selection (rectangular/mask) of the grid.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments(argc,argv,program_description);
    add_argument_s(argc,argv,"-d"     , in_file_name,               "Input data file to be averaged (dat)"                                                , cl_tags, nullptr,      1);
    add_argument_s(argc,argv,"-mask"  , mask_file_name,             "Mask file used to select grid points in average (dat)"                               , cl_tags, &b_mask,      0);
    add_argument_s(argc,argv,"-list"  , list_file_name,             "Output data file with the value of each lattice point included in the average (dat)" , cl_tags, &b_list,      0);
    add_argument_i(argc,argv,"-x"     , &target_x,                  "Rectangle center x (grid point)"                                                     , cl_tags, nullptr,      1);
    add_argument_i(argc,argv,"-y"     , &target_y,                  "Rectangle center y (grid point)"                                                     , cl_tags, nullptr,      1);
    add_argument_i(argc,argv,"-rx"    , &range_x,                   "Rectangle half width x (grid points)  "                                              , cl_tags, nullptr,      1);
    add_argument_i(argc,argv,"-ry"    , &range_y,                   "Rectangle half width y (grid points)  "                                              , cl_tags, nullptr,      1);
    add_argument_i(argc,argv,"-invert", &invert,                    "Invert rectangular selection? (0:no 1:yes)"                                          , cl_tags, nullptr,      1);
    add_argument_i(argc,argv,"-odf"   , &odf,                       "Data file format (0:matrix 1:vector)"                                                , cl_tags, nullptr,      1);
    conclude_input_arguments(argc,argv,program_name,cl_tags);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension("-d",in_file_name,".dat");
    if(b_mask == 1)
    {
        check_extension("-mask",mask_file_name,".dat");
    }
    if(b_list == 1)
    {
        check_extension("-list",list_file_name,".dat");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Create Grids                                                                                              //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Grid_lt data;
    Grid_lt mask;
    Grid_lt selection;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set the grid format                                                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    data.set_format(odf);
    mask.set_format(odf);
    selection.set_format(odf);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Create Vector to store added values                                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    dv1d list(0,0.0);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read in grid data                                                                                         //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    data.get_grid(in_file_name);
    if(b_mask == 1)
    {
        mask.get_grid(mask_file_name);
    }
    else
    {
        mask.get_grid(in_file_name);
        mask.init_grid(1.0);
        mask.set_nan(0);
    }
    selection.get_grid(in_file_name);
    selection.init_grid(0.0);
    selection.set_nan(0);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Print info about grids                                                                                    //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    data.print_dim(1);
    if(b_mask == 1)
    {
        mask.print_dim(0);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Compute the average                                                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(i=0; i<data.size_x(); i++) //loop over x
    {
        for(j=0; j<data.size_y(); j++) //loop over y
        {
            if(mask.grid[i][j][2][0] == 1.0) //passes mask check
            { 
                int rectangle_pass = check_rectangular_pass_vector(invert,target_x,target_y,range_x,range_y,i,j);

                if(rectangle_pass == 1) //passes rectangle check
                {
                    if(data.grid[i][j][2][1] == 0) //data present
                    {
                        selection.grid[i][j][2][0] = 1;
                        selection.grid[i][j][2][1] = 0;
                        avg = avg + data.grid[i][j][2][0];
                        count++;
                        list.push_back(data.grid[i][j][2][0]);
                    }
                }
            }
        }
    }
    avg = avg/(double)count;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Compute the spread                                                                                        //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(i=0; i<data.size_x(); i++) //loop over x
    {
        for(j=0; j<data.size_y(); j++) //loop over y
        {
            if(mask.grid[i][j][2][0] == 1.0) //passes mask check
            {
                int rectangle_pass = check_rectangular_pass_vector(invert,target_x,target_y,range_x,range_y,i,j);

                if(rectangle_pass == 1) //passes rectangle check
                {
                    if(data.grid[i][j][2][1] == 0) //data present
                    {
                        sum_o_squares = sum_o_squares + pow(data.grid[i][j][2][0] - avg,2);
                    }
                }
            }
        }
    }
    stdev = sqrt(sum_o_squares/count);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Write grid selection to output file                                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    selection.write_grid("selection.dat");

    if(b_list == 1)
    {
        FILE *list_file = fopen(list_file_name.c_str(), "w");
        fprintf(list_file," %s \n","#values of <F^ij> included in average");
        for(i=0; i<list.size(); i++)
        {
            fprintf(list_file," %f \n",list[i]);
        }
        fclose(list_file);
    }

    printf("\n");
    printf("%8s %d \n","x:",target_x);
    printf("%8s %d \n","y:",target_y);
    printf("%8s %d \n","rx:",range_x);
    printf("%8s %d \n","ry:",range_y);
    printf("%8s %d \n","count:",count);
    printf("%8s %f \n","sum:",avg*(double)count);
    printf("%8s %f \n","average:",avg);
    printf("%8s %f \n","stdev:",stdev);
    std::cout << "\nFormatting completed successfully" << "\n\n";

    return 0;
}


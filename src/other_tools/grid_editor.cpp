
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
    string in_file_name;          //Name of input file data file
    string out_file_name;         //Name of output file
    string mask_file_name;        //Name of mask file
    int i           = 0;          //General variable used in loops
    int j           = 0;          //General variable used in loops
    int k           = 0;          //General variable used in loops
    int l           = 0;          //General variable used in loops
    int odf         = 0;          //Data file format
    int invert      = 0;          //Invert rectangular selection?
    int b_mask      = 0;          //Was a masking file provided?
    double val      = 0.0;        //Value to set select grid point equal to
    double nan      = 0.0;        //Value added to grid when NaN is encountered
    double target_x = 0;          //The target grid point x when making a rectangular selection
    double target_y = 0;          //The target grid point y when making a rectangular selection
    double range_x  = 0;          //The half width of x in the rectangular selection
    double range_y  = 0;          //The half width of y in the rectangular selection
    sv1d cl_tags;                 //Holds a list of command line tags for the program

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set program name and print info                                                                           //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    string program_name = "Grid Editor";

    print_credits(argc,argv,program_name);

    string program_description = "Grid Editor is an analysis tool that lets the user set the value for a selection (rectangular or mask) of the grid.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments(argc,argv,program_description);
    add_argument_s(argc,argv,"-d"     , in_file_name,               "Input file with grid data to be edited (dat)" , cl_tags, nullptr,      1);
    add_argument_s(argc,argv,"-o"     , out_file_name,              "Output data file with edited grid data (dat)" , cl_tags, nullptr,      1);
    add_argument_d(argc,argv,"-x"     , &target_x,                  "Rectangle center x (grid point)"              , cl_tags, nullptr,      1);
    add_argument_d(argc,argv,"-y"     , &target_y,                  "Rectangle center y (grid point)"              , cl_tags, nullptr,      1);
    add_argument_d(argc,argv,"-rx"    , &range_x,                   "Rectangle half width x (grid points)  "       , cl_tags, nullptr,      1);
    add_argument_d(argc,argv,"-ry"    , &range_y,                   "Rectangle half width y (grid points)  "       , cl_tags, nullptr,      1);
    add_argument_d(argc,argv,"-val"   , &val,                       "Value to set grid points to"                  , cl_tags, nullptr,      1);
    add_argument_i(argc,argv,"-odf"   , &odf,                       "Data file format (0:matrix 1:vector)"         , cl_tags, nullptr,      1);
    add_argument_s(argc,argv,"-mask"  , mask_file_name,             "Mask file for selecting lattice points (dat)" , cl_tags, &b_mask,      0);
    conclude_input_arguments(argc,argv,program_name,cl_tags);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension("-d",in_file_name,".dat");
    check_extension("-o",out_file_name,".dat");

    if(b_mask == 1)
    {
        check_extension("-mask",mask_file_name,".dat");
    }

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

    data.print_dim(1);
    if(b_mask == 1)
    {
        mask.print_dim(0);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set select points to val                                                                                  //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(i=0; i<data.size_x(); i++) //loop over x
    {
        for(j=0; j<data.size_y(); j++) //loop over y
        {
            int rectangle_pass = check_rectangular_pass_vector(invert,target_x,target_y,range_x,range_y,i,j);

            if(rectangle_pass == 1)
            {
                if(b_mask == 1)
                {
                    if(mask.grid[i][j][2][0] == 1)
                    {
                        data.grid[i][j][2][0] = val;
                        data.grid[i][j][2][1] = 0;
                    }
                }
                else
                {
                    data.grid[i][j][2][0] = val;
                    data.grid[i][j][2][1] = 0;
                }
            }
        }             
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Write grid to output file                                                                                 //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    data.write_grid(out_file_name);

    std::cout << "\nFormatting completed successfully" << "\n\n";

    return 0;
}


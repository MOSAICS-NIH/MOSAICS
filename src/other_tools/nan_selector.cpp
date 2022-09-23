
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
    string out_file_name;         //Name of output mask file
    int i         = 0;            //General variable used in loops
    int j         = 0;            //General variable used in loops
    int k         = 0;            //General variable used in loops
    int l         = 0;            //General variable used in loops
    int target_x  = 0;            //The target grid point x
    int target_y  = 0;            //The target grid point y
    int odf       = 0;            //Data file format
    int mask_grew = 0;            //Did the mask grow in previous iteration

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set program name and print info                                                                           //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    string program_name = "Nan Selector";

    print_credits(argc,argv,program_name);

    string program_description = "Nan Selector is a tool used to convert a region of grid data containing NaN into a selection mask.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments(argc,argv,program_description);
    add_argument_s(argc,argv,"-d"     , in_file_name,               "Input grid data file (dat)"                  , nullptr,      1);
    add_argument_s(argc,argv,"-o"     , out_file_name,              "Output mask file (dat)"                      , nullptr,      1);
    add_argument_i(argc,argv,"-x"     , &target_x,                  "Seed x (grid point)"                         , nullptr,      1);
    add_argument_i(argc,argv,"-y"     , &target_y,                  "Seed y (grid point)"                         , nullptr,      1);
    add_argument_i(argc,argv,"-odf"   , &odf,                       "Data file format (0:matrix 1:vector)"        , nullptr,      1);
    conclude_input_arguments(argc,argv,program_name);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension("-d",in_file_name,".dat");
    check_extension("-o",out_file_name,".dat");

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
    mask.get_grid(in_file_name);
    mask.init_grid(0.0);
    mask.set_nan(0);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Print info about grids                                                                                    //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    data.print_dim(1);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set nan to -99                                                                                            //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(i=0; i<data.size_x(); i++) //loop over x
    {
        for(j=0; j<data.size_y(); j++) //loop over y
        {
            if(data.grid[i][j][2][1] == 1) //nan
            {
                data.grid[i][j][2][1] = 0;
                data.grid[i][j][2][0] = -99;
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Grow the mask                                                                                             //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    mask.grow_mask(data.grid,target_x,target_y,-99.0);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Write grid to output file                                                                                 //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    mask.write_grid(out_file_name);

    std::cout << "\nFormatting completed successfully" << "\n\n";

    return 0;
}



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
    string in_file_name_1;             //Name of input file 1
    string in_file_name_2;             //Name of input file 2
    string out_file_name;              //Name of output file
    int i      = 0;                    //General variable used in loops
    int j      = 0;                    //General variable used in loops
    int odf    = 0;                    //Data file format 
    double nan = -999999999999;        //Number to add to grid when NaN is encountered
    sv1d cl_tags;                      //Holds a list of command line tags for the program

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set program name/description and print info                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    string program_name = "Delta Plot";

    print_credits(argc,argv,program_name);

    string program_description = "Delta Plot is an analysis tool used for computing the difference between 2 grids.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments(argc,argv,program_description);
    add_argument_s(argc,argv,"-d1"    , in_file_name_1,             "Input grid data file 1 (subtract me) (dat)"  , cl_tags, nullptr,      1);
    add_argument_s(argc,argv,"-d2"    , in_file_name_2,             "Input grid data file 2 (dat)"                , cl_tags, nullptr,      1);
    add_argument_s(argc,argv,"-o"     , out_file_name,              "Output grid data file with delta (dat) "     , cl_tags, nullptr,      1);
    add_argument_i(argc,argv,"-odf"   , &odf,                       "Data file format (0:matrix 1:vector)"        , cl_tags, nullptr,      1);
    conclude_input_arguments(argc,argv,program_name,cl_tags);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension("-d1",in_file_name_1,".dat");
    check_extension("-d2",in_file_name_2,".dat");
    check_extension("-o",out_file_name,".dat");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Create Grids                                                                                              //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Grid_lt data_1;
    Grid_lt data_2;
    Grid_lt delta;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set the grid format                                                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    data_1.set_format(odf);
    data_2.set_format(odf);
    delta.set_format(odf);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read in grid data                                                                                         //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    data_1.get_grid(in_file_name_1);
    data_2.get_grid(in_file_name_2);
    delta.get_grid(in_file_name_1);    //read in first grid to get dimensions
    delta.init_grid(0.0);              //initialize grid to zero
    delta.set_nan(0);                  //initialize nan tags

    int check = comp_grid(data_1,data_2);

    if(check == 1)
    {
        data_1.print_dim(1);
        data_2.print_dim(0);
    }
    else
    {
        printf("Input files are not compatible. \n");
        data_1.print_dim(1);
        data_2.print_dim(0);
        return 0;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Compute the difference                                                                                    //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(i=0; i<delta.size_x(); i++) //loop over x
    {
        for(j=0; j<delta.size_y(); j++) //loop over y
        {
            delta.grid[i][j][0][0] = data_1.grid[i][j][0][0];  //set x value
            delta.grid[i][j][1][0] = data_1.grid[i][j][1][0];  //sey y value

            if(data_1.grid[i][j][2][1] == 1 || data_2.grid[i][j][2][1] == 1) //NaN encountered
            {
                delta.grid[i][j][2][0] = nan;
                delta.grid[i][j][2][1] = 1;
            }
            else
            {
                delta.grid[i][j][2][0] = data_2.grid[i][j][2][0] - data_1.grid[i][j][2][0];
                delta.grid[i][j][2][1] = 0;
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Write delta to output file                                                                                //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    delta.write_grid(out_file_name);

    std::cout << "\nFormatting completed successfully" << "\n\n";

    return 0;
}


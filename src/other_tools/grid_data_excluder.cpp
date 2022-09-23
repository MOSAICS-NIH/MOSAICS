
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
    string in_file_name;          //Name of input file 1
    string rho_file_name;         //Name of the rho file 
    string out_file_name;         //Name of the output data file
    int i          = 0;           //General variable used in loops
    int j          = 0;           //General variable used in loops
    int odf        = 0;           //Data file format
    double avg_rho = 0;           //The average lipid density over the grid
    double cutoff  = 0;           //Cutoff for excluding data
    double nan     = 0.0;         //Value added to grid when NaN is encountered

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set program name/description and print info                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    string program_name = "Grid Data Excluder";

    print_credits(argc,argv,program_name);

    string program_description = "Grid Data Excluder is an analysis tool used for excluding insignificant grid data based on the sample count data.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments(argc,argv,program_description);
    add_argument_s(argc,argv,"-d"     , in_file_name,               "Input data file with grid data (dat)   "     , nullptr,      1);
    add_argument_s(argc,argv,"-rho"   , rho_file_name,              "Input data file with sample count (dat)  "   , nullptr,      1);
    add_argument_s(argc,argv,"-o"     , out_file_name,              "Output data file with grid data (dat)"       , nullptr,      1);
    add_argument_d(argc,argv,"-cutoff", &cutoff,                    "Cutoff for excluding data (chi)"             , nullptr,      1);
    add_argument_i(argc,argv,"-odf"   , &odf,                       "Data file format (0:matrix 1:vector)"        , nullptr,      1);
    conclude_input_arguments(argc,argv,program_name);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension("-d",in_file_name,".dat");
    check_extension("-rho",rho_file_name,".dat");
    check_extension("-o",out_file_name,".dat");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Create Grids                                                                                              //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Grid_lt rho;
    Grid_lt data;
    Grid_lt nan_orig; 

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set the grid format                                                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    rho.set_format(odf);
    data.set_format(odf);
    nan_orig.set_format(odf);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read in grid data                                                                                         //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    rho.get_grid(rho_file_name);
    data.get_grid(in_file_name);
    nan_orig.get_grid(in_file_name);

    //check that grids are compatible
    int check = comp_grid(data,rho);
    if(check == 0)
    {
        printf("Input files are not compatible. ");

        data.print_dim(1);
        rho.print_dim(0);

        return 0;
    }
    else
    {
        data.print_dim(1);
        rho.print_dim(0);
    }

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
    data.exclude_grid_data(cutoff,avg_rho,rho.grid);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Retain any initial nan grid points                                                                        //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(i=0; i<rho.size_x(); i++) //loop over x
    {
        for(j=0; j<rho.size_y(); j++) //loop over y
        {
            if(nan_orig.grid[i][j][2][1] == 1)
            {
                data.grid[i][j][2][1] = nan_orig.grid[i][j][2][1];
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


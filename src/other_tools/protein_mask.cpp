
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
    string rho_file_name;         //Name of input file data file
    string out_file_name;         //Name of output file
    int i          = 0;           //General variable used in loops
    int j          = 0;           //General variable used in loops
    int odf        = 0;           //Data file format
    double nan     = 0.0;         //Value added to grid when NaN is encountered
    double cutoff  = 0;           //Cutoff for excluding data
    double avg_rho = 0;           //The average lipid density over the grid
    sv1d cl_tags;                 //Holds a list of command line tags for the program

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set program name and print info                                                                           //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    string program_name = "Protein Mask";

    print_credits(argc,argv,program_name);

    string program_description = "Protein Mask is an analysis tool used for making a mask of the protein based on lipid sample count.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments(argc,argv,program_description);
    add_argument_s(argc,argv,"-d"     , rho_file_name,              "Input data file with sample count data (dat)"                              , cl_tags, nullptr,      1);
    add_argument_s(argc,argv,"-o"     , out_file_name,              "Output data file with mask indicating regions with low sample count (dat)" , cl_tags, nullptr,      1);
    add_argument_d(argc,argv,"-cutoff", &cutoff,                    "Cutoff for selecting regions with a low sample count (chi)"                , cl_tags, nullptr,      1);
    add_argument_i(argc,argv,"-odf"   , &odf,                       "Data file format (0:matrix 1:vector)"                                      , cl_tags, nullptr,      1);
    conclude_input_arguments(argc,argv,program_name,cl_tags);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension("-d",rho_file_name,".dat");
    check_extension("-o",out_file_name,".dat");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Create Grids                                                                                              //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Grid_lt rho;

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

    rho.print_dim(1);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Compute average rho                                                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    avg_rho = rho.get_average();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Get Protein mask                                                                                          //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(i=0; i<rho.size_x(); i++) //loop over x
    {
        for(j=0; j<rho.size_y(); j++) //loop over y
        {
            if(rho.grid[i][j][2][0] > cutoff*avg_rho)
            {
                rho.grid[i][j][2][0] = 0;
            }
            else
            {
                rho.grid[i][j][2][0] = 1;
            }
            rho.grid[i][j][2][1] = 0;
        }
    }            

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Write grid to output file                                                                                 //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    rho.write_grid(out_file_name);

    std::cout << "\nFormatting completed successfully" << "\n\n";

    return 0;
}


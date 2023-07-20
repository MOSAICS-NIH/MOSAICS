
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
    string in_file_name;          //Name of input mask file
    string out_file_name;         //Name of output file
    int i                 = 0;    //General variable used in loops
    int j                 = 0;    //General variable used in loops
    int odf               = 0;    //Data file format
    int current_iteration = 1;    //The current iteration
    int target_x          = 0;    //The target grid point x when making a rectangular selection
    int target_y          = 0;    //The target grid point y when making a rectangular selection
    int range_x           = 0;    //The half width of x in the rectangular selection
    int range_y           = 0;    //The half width of y in the rectangular selection
    int invert            = 0;    //Invert rectangular selection (select everything outside rectangle)
    int cumulative        = 1;    //Include all grid points up to dist?
    int b_prot            = 0;    //Copy protein mask onto new mask?
    double cell_size      = 1;    //Distance between grid points
    double res            = 1;    //How much does the grid selection move with each iteration
    double range          = 1;    //The half width of the grid point selection
    double nan            = 0.0;  //Value added to grid when NaN is encountered
    double APS            = 0;    //Area per square used for grid in analysis
    double dist           = 0;    //How far from protein to make mask
    sv1d cl_tags;                 //Holds a list of command line tags for the program

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set program name/description and print info                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    string program_name = "Mask Maker";

    print_credits(argc,argv,program_name);

    string program_description = "Mask Maker is an analysis tool used for making a mask around the protein.";
          
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments(argc,argv,program_description);
    add_argument_s(argc,argv,"-d"        , in_file_name,            "Input protein mask file (dat)"                    , cl_tags, nullptr,      1);
    add_argument_s(argc,argv,"-o"        , out_file_name,           "Output data file with the generated mask (dat)"   , cl_tags, nullptr,      1);
    add_argument_i(argc,argv,"-x"        , &target_x,               "Rectangle center x (grid point)"                  , cl_tags, nullptr,      1);
    add_argument_i(argc,argv,"-y"        , &target_y,               "Rectangle center y (grid point)"                  , cl_tags, nullptr,      1);
    add_argument_i(argc,argv,"-rx"       , &range_x,                "Rectangle half width x (grid points)  "           , cl_tags, nullptr,      1);
    add_argument_i(argc,argv,"-ry"       , &range_y,                "Rectangle half width y (grid points)  "           , cl_tags, nullptr,      1);
    add_argument_i(argc,argv,"-invert"   , &invert,                 "Invert rectangular selection (0:no,1:yes)"        , cl_tags, nullptr,      1);
    add_argument_d(argc,argv,"-dist"     , &dist,                   "How far from the protein does mask extend (nm)"   , cl_tags, nullptr,      1);
    add_argument_d(argc,argv,"-APS"      , &APS,                    "Area per grid point (nm^2)"                       , cl_tags, nullptr,      1);
    add_argument_i(argc,argv,"-copy"     , &b_prot,                 "Copy the protein mask onto new mask (0:no 1:yes)" , cl_tags, nullptr,      1);
    add_argument_i(argc,argv,"-odf"      , &odf,                    "Data file format (0:matrix 1:vector)"             , cl_tags, nullptr,      1);
    conclude_input_arguments(argc,argv,program_name,cl_tags);

    cell_size = sqrt(APS);

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
    Grid_lt init_mask;
    Grid_lt mask;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set the grid format                                                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    init_mask.set_format(odf);
    mask.set_format(odf);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read in grid data                                                                                         //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    init_mask.get_grid(in_file_name);
    mask.get_grid(in_file_name);
    mask.init_grid(0.0);              //initialize grid to zero
    mask.set_nan(0);                  //initialize nan tags

    init_mask.print_dim(1);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Compute the grid selection                                                                                //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    range             = 0;
    current_iteration = 1;
    res               = dist;
    cumulative        = 1;

    mask.distance_projection(invert,target_x,target_y,range_x,range_y,cell_size,current_iteration,range,res,init_mask.grid,cumulative);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Copy Protein mask                                                                                         //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(b_prot == 1)
    {
        for(i=0; i<mask.size_x(); i++) //loop over x
        {
            for(j=0; j<mask.size_y(); j++) //loop over y
            {
                if(init_mask.grid[i][j][2][0] == 1) //protein
                {
                    mask.grid[i][j][2][0] = 1;
                }
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Write the mask to output file                                                                             //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    mask.write_grid(out_file_name);

    std::cout << "\nFormatting completed successfully" << "\n\n";

    return 0;
}


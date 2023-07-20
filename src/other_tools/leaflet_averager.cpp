
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
    string in_file_name_1;         //Name of input file 1
    string in_file_name_2;         //Name of input file 2
    string rho_file_name_1;        //Name of the rho file 1
    string rho_file_name_2;        //Name of the rho file 2
    string out_file_name;          //Name of the output data file
    int i          = 0;            //General variable used in loops
    int j          = 0;            //General variable used in loops
    int odf        = 0;            //Data file format
    int nan        = 0;            //Number to be added to grid when NaN is encountered
    double avg_rho = 0;            //The average lipid density over the grid
    double cutoff  = 0;            //Cutoff for excluding data
    sv1d cl_tags;                  //Holds a list of command line tags for the program

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set program name/description and print info                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    string program_name = "Leaflet Averager";

    print_credits(argc,argv,program_name);

    string program_description = "Leaflet Averager is an analysis tool used for averaging grid data over 2 leaflets.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments(argc,argv,program_description);
    add_argument_s(argc,argv,"-d1"    , in_file_name_1,             "Input data file with grid data for leaflet 1 (dat)"    , cl_tags, nullptr,      1);
    add_argument_s(argc,argv,"-d2"    , in_file_name_2,             "Input data file with grid data for leaflet 2 (dat)"    , cl_tags, nullptr,      1);
    add_argument_s(argc,argv,"-rho1"  , rho_file_name_1,            "Input data file with sample count for leaflet 1 (dat)" , cl_tags, nullptr,      1);
    add_argument_s(argc,argv,"-rho2"  , rho_file_name_2,            "Input data file with sample count for leaflet 2 (dat)" , cl_tags, nullptr,      1);
    add_argument_s(argc,argv,"-o"     , out_file_name,              "Output data file with averaged grid data (dat)"        , cl_tags, nullptr,      1);
    add_argument_d(argc,argv,"-cutoff", &cutoff,                    "Cutoff for excluding data (chi)"                       , cl_tags, nullptr,      1);
    add_argument_i(argc,argv,"-odf"   , &odf,                       "Data file format (0:matrix 1:vector)"                  , cl_tags, nullptr,      1);
    conclude_input_arguments(argc,argv,program_name,cl_tags);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension("-d1",in_file_name_1,".dat");
    check_extension("-d2",in_file_name_2,".dat");
    check_extension("-rho1",rho_file_name_1,".dat");
    check_extension("-rho2",rho_file_name_2,".dat");
    check_extension("-o",out_file_name,".dat");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Create Grids                                                                                              //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Grid_lt data_1;
    Grid_lt data_2;
    Grid_lt rho_1;
    Grid_lt rho_2;
    Grid_lt rho;
    Grid_lt weight_1;
    Grid_lt weight_2;
    Grid_lt avg;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set the grid format                                                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    data_1.set_format(odf);
    data_2.set_format(odf);
    rho_1.set_format(odf);
    rho_2.set_format(odf);
    rho.set_format(odf);
    weight_1.set_format(odf);
    weight_2.set_format(odf);
    avg.set_format(odf);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read in grid data                                                                                         //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    data_1.get_grid(in_file_name_1);
    data_2.get_grid(in_file_name_2);
    rho_1.get_grid(rho_file_name_1);
    rho_2.get_grid(rho_file_name_2);
    rho.get_grid(rho_file_name_1);        //read in first grid to get dimensions
    weight_1.get_grid(rho_file_name_1);   //read in first grid to get dimensions
    weight_2.get_grid(rho_file_name_1);   //read in first grid to get dimensions
    avg.get_grid(in_file_name_1);         //read in first grid to get dimensions
    rho.init_grid(0.0);                   //initialize grid to zero
    weight_1.init_grid(0.0);              //initialize grid to zero
    weight_2.init_grid(0.0);              //initialize grid to zero
    avg.init_grid(0.0);                   //initialize grid to zero
    rho.set_nan(0);                       //initialize nan tags
    weight_1.set_nan(0);                  //initialize nan tags
    weight_2.set_nan(0);                  //initialize nan tags
    avg.set_nan(0);                       //initialize nan tags

    int check_1 = comp_grid(data_1,data_2);
    int check_2 = comp_grid(data_1,rho_1);
    int check_3 = comp_grid(data_1,rho_2);

    if(check_1 == 1 && check_2 == 1 && check_3 == 1)
    {
        data_1.print_dim(1);
        data_2.print_dim(0);
        rho_1.print_dim(0);
        rho_2.print_dim(0);
    }
    else
    {
        printf("Input files are not compatible. \n");

        data_1.print_dim(1);
        data_2.print_dim(0);
        rho_1.print_dim(0);
        rho_2.print_dim(0);

        return 0;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Get the combined lipid density                                                                            //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(i=0; i<rho.size_x(); i++) //loop over x
    {
        for(j=0; j<rho.size_y(); j++) //loop over y
        {
            rho.grid[i][j][0][0] = rho_1.grid[i][j][0][0];
            rho.grid[i][j][1][0] = rho_1.grid[i][j][1][0];
            rho.grid[i][j][2][0] = rho_1.grid[i][j][2][0] + rho_2.grid[i][j][2][0];

            //set the NaN number flags
            if(rho_1.grid[i][j][2][1] == 1 && rho_2.grid[i][j][2][1] == 1)
            {
                rho.grid[i][j][2][1] == 1;
            }
            else 
            {
                rho.grid[i][j][2][1] = 0;
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Compute weights for each leaflet                                                                          //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(i=0; i<weight_1.size_x(); i++) //loop over x
    {
        for(j=0; j<weight_1.size_y(); j++) //loop over y
        {
            double norm = rho_1.grid[i][j][2][0] + rho_2.grid[i][j][2][0];

            if(norm > 0)
            {
                weight_1.grid[i][j][2][0] = rho_1.grid[i][j][2][0]/norm;
                weight_2.grid[i][j][2][0] = rho_2.grid[i][j][2][0]/norm;
            }
            else
            {
                weight_1.grid[i][j][2][0] = 0;
                weight_2.grid[i][j][2][0] = 0;
            } 
        } 
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Compute average rho                                                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    avg_rho = rho.get_average();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Average Leaflets                                                                                          //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(i=0; i<avg.size_x(); i++) //loop over x
    {
        for(j=0; j<avg.size_y(); j++) //loop over y
        {
            avg.grid[i][j][0][0] = data_1.grid[i][j][0][0];
            avg.grid[i][j][1][0] = data_1.grid[i][j][1][0];
            avg.grid[i][j][2][0] = weight_1.grid[i][j][2][0]*data_1.grid[i][j][2][0] + weight_2.grid[i][j][2][0]*data_2.grid[i][j][2][0];
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Exclude insignificant data and write average to output file                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    avg.exclude_grid_data(cutoff,avg_rho,rho.grid);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Write grid to output file                                                                                 //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    avg.write_grid(out_file_name);

    std::cout << "\nFormatting completed successfully" << "\n\n";

    return 0;
}


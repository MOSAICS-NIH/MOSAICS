
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
#include "../headers/file_naming.h"
#include "../headers/grid_lt.h"
#include "../headers/common_routines.h"
#include "../headers/file_io.h"
#include "../headers/command_line_args.h"
#include "../headers/array.h"

int main(int argc, const char * argv[]) 
{
    //Here we define some variables used throughout
    FILE *fe_file;                         //File for reading in free energy data
    char my_string[200];                   //String to hold read in data entries
    string in_file_name;                   //Name of input file data file
    string rho_file_name;                  //Name of the density file
    string out_file_name;                  //Name of output file
    string fe_file_name;                   //Name of the data file with free energy values
    int i                  = 0;            //General variable used in loops
    int j                  = 0;            //General variable used in loops
    int k                  = 0;            //General variable used in loops
    int l                  = 0;            //General variable used in loops
    int number_of_lines_fe = 0;            //Number of lines in free energy file 
    int capacity_fe        = 0;            //Number of items in the fe file
    int items_per_line_fe  = 0;            //Items per line in the fe file
    int odf                = 0;            //Data file format
    double lipids_tot      = 0;            //How many lipids in the grid
   
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set program name/description and print info                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    string program_name = "Bilayer Free Energy";

    print_credits(argc,argv,program_name);

    string program_description = "Bilayer Free Energy is an analysis tool that takes a free energy profile and thickness data and generates a free energy projection in XY.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments(argc,argv,program_description);
    add_argument_s(argc,argv,"-d"     , in_file_name,               "Time average membrane thickness grid data (dat)"                     , nullptr,      1);
    add_argument_s(argc,argv,"-rho"   , rho_file_name,              "Normalized lipid density grid data (dat)"                            , nullptr,      1);
    add_argument_s(argc,argv,"-fe"    , fe_file_name,               "Free energy profile (used to look up cost based on thickness) (dat)" , nullptr,      1);
    add_argument_s(argc,argv,"-o"     , out_file_name,              "Output grid data file with free energy penaly (dat) "                , nullptr,      1);
    add_argument_i(argc,argv,"-odf"   , &odf,                       "Data file format (0:matrix 1:vector)"                                , nullptr,      1);
    conclude_input_arguments(argc,argv,program_name);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension("-d",in_file_name,".dat");
    check_extension("-rho",rho_file_name,".dat");
    check_extension("-fe",fe_file_name,".dat");
    check_extension("-o",out_file_name,".dat");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read in free energy lookup data                                                                           //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Data_file fe_lookup;
    fe_lookup.set_headers(7);
    fe_lookup.get_data(fe_file_name);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Create Grids                                                                                              //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Grid_lt data;
    Grid_lt rho;
    Grid_lt fe;
    Grid_lt fe_pair;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set the grid format                                                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    data.set_format(odf);
    rho.set_format(odf);
    fe.set_format(odf);
    fe_pair.set_format(odf);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read in grid data                                                                                         //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    data.get_grid(in_file_name);
    rho.get_grid(rho_file_name);
    fe.get_grid(in_file_name);          //for fe we copy the dimensions of data
    fe.init_grid(0.0);                  //initialize grid to 0
    fe.set_nan(0);                      //initialize nan tags
    fe_pair.get_grid(in_file_name);     //for fe we copy the dimensions of data
    fe_pair.init_grid(0.0);             //initialize grid to 0
    fe_pair.set_nan(0);                 //initialize nan tags

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Print info about grids                                                                                    //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int check = comp_grid(data,rho);
    if(check == 1)
    {
        data.print_dim(1);
        rho.print_dim(0);
    }
    else 
    {
        printf("Input files are not compatible. \n");
        data.print_dim(1);
        rho.print_dim(0);
        return 0;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Compute free energy for each grid point                                                                   //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(i=0; i<data.size_x(); i++) //loop over x
    {
        for(j=0; j<data.size_y(); j++) //loop over  y
        {
            //set the x and y values for the fe grid
            fe.grid[i][j][0][0]      = data.grid[i][j][0][0];
            fe.grid[i][j][1][0]      = data.grid[i][j][1][0];
            fe_pair.grid[i][j][0][0] = data.grid[i][j][0][0];
            fe_pair.grid[i][j][1][0] = data.grid[i][j][1][0];

            if(data.grid[i][j][2][1] == 0 && rho.grid[i][j][2][1] == 0) //thickness data present
            {
                int    found_thickness = 0;
                double thickness       = data.grid[i][j][2][0];
                double lipids          = rho.grid[i][j][2][0];
                lipids_tot             = lipids_tot + lipids;

                for(k=0; k<fe_lookup.size_y(); k++) //loop over fe data
                {
                    if(fe_lookup.data_d[k][0] <= thickness && fe_lookup.data_d[k+1][0] > thickness)
                    { 
                        if(strcmp(fe_lookup.data_s[k][4].c_str(), "NaN") == 0)
                        {
                            printf("grid point %d %d is out of range thickness %f \n",i,j,thickness);
                        }
                        fe.grid[i][j][2][0]       = (lipids)*fe_lookup.data_d[k][4];
                        fe.grid[i][j][2][1]       = 0;
                        fe_pair.grid[i][j][2][0]  = fe_lookup.data_d[k][4];
                        fe_pair.grid[i][j][2][1]  = 0;
                        found_thickness           = 1;
                    }
                }

                //check if the thickness was found in the fe_lookup
                if(found_thickness == 0)
                {
                    printf("Could not find thickness %f in fe_lookup for gridpoint x %d y %d \n",thickness,i,j);
                    fe.grid[i][j][2][1] = 1;
                    fe_pair.grid[i][j][2][1] = 1;
                }
            }
            else //no thickness data present
            {
                fe.grid[i][j][2][0]      = -999999999;
                fe.grid[i][j][2][1]      = 1;
                fe_pair.grid[i][j][2][0] = -999999999;
                fe_pair.grid[i][j][2][1] = 1;
            }
        }
    }            

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Write grid to output file                                                                                 //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    fe.write_grid(out_file_name);

    string fe_pair_file_name = add_tag(out_file_name,"_pair");
    fe_pair.write_grid(fe_pair_file_name);

    printf("lipids_tot %f \n",lipids_tot);

    std::cout << "\nFormatting completed successfully" << "\n\n";

    return 0;
}



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
#include "../headers/file_io.h"
#include "../headers/command_line_args.h"
#include "../headers/array.h"
#include "../headers/file_naming.h"

int main(int argc, const char * argv[]) 
{
    //Here we define some variables used throughout
    string in_file_name_1;                 //Name of input file data file 1
    string in_file_name_2;                 //Name of input file data file 2
    string out_file_name;                  //Name of output file
    int i  = 0;                            //General variable used in loops
    int j  = 0;                            //General variable used in loops
    int h1 = 0;                            //Number of header lines for data file 1
    int h2 = 0;                            //Number of header lines for data file 2
    sv1d cl_tags;                          //Holds a list of command line tags for the program
   
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set program name/description and print info                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    string program_name = "Data Averager";

    print_credits(argc,argv,program_name);

    string program_description = "Data Averager is an analysis tool used for averaging the contents of 2 data files.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments(argc,argv,program_description);
    add_argument_s(argc,argv,"-d1"    , in_file_name_1,             "Data file 1  "                               , cl_tags, nullptr,      1);
    add_argument_s(argc,argv,"-d2"    , in_file_name_2,             "Data file 2  "                               , cl_tags, nullptr,      1);
    add_argument_i(argc,argv,"-h1"    , &h1,                        "Number of header lines in data file 1  "     , cl_tags, nullptr,      0);
    add_argument_i(argc,argv,"-h2"    , &h2,                        "Number of header lines in data file 2  "     , cl_tags, nullptr,      0);
    add_argument_s(argc,argv,"-o"     , out_file_name,              "Output data file "                           , cl_tags, nullptr,      1);
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
    // Create data files                                                                                         //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Data_file data_1;
    Data_file data_2;
    Data_file data_avg;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set number of header lines                                                                                //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    data_1.set_headers(h1);
    data_2.set_headers(h2);
    data_avg.set_headers(h1);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read data files                                                                                           //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    data_1.get_data(in_file_name_1);
    data_2.get_data(in_file_name_2);
    data_avg.get_data(in_file_name_1);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Compute average for the data                                                                              //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for(i=0; i<data_1.size_y(); i++) //loop over rows
    {
        for(j=0; j<data_1.size_x(); j++) //loop over columns
        {
            data_avg.data_d[i][j] = 0.5*(data_1.data_d[i][j] + data_2.data_d[i][j]);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Write grid to output file                                                                                 //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    data_avg.write_data(out_file_name,1);

    std::cout << "\nFormatting completed successfully" << "\n\n";

    return 0;
}


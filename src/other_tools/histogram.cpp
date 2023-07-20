
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
#include "../headers/common_routines.h"
#include "../headers/file_io.h"
#include "../headers/histo.h"
#include "../headers/command_line_args.h"
#include "../headers/array.h"
#include "../headers/file_naming.h"

int main(int argc, const char * argv[]) 
{
    //Here we define some variables used throughout
    FILE *out_file;               //File for writing data
    string in_file_name;          //Name of input file 1
    string out_file_name;         //Name of output file
    int i                = 0;     //General variable used in loops
    int j                = 0;     //General variable used in loops
    int col              = 1;     //What collumn do we use data for
    double bin_width     = 1;     //bin width
    double largest_value = 0;     //largest value in data set
    double percent       = 0;     //Percent reported for the bin
    sv1d cl_tags;                 //Holds a list of command line tags for the program

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set program name/description and print info                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    string program_name = "Histogram";

    print_credits(argc,argv,program_name);

    string program_description = "Histogram is an anaylsis tool that bins a data set to produce a histogram.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments(argc,argv,program_description);
    add_argument_s(argc,argv,"-d",   in_file_name,               "Input data file"                         , cl_tags, nullptr,      1);
    add_argument_s(argc,argv,"-o",   out_file_name,              "Output file"                             , cl_tags, nullptr,      1);
    add_argument_i(argc,argv,"-col", &col,                       "What collumn to work with"               , cl_tags, nullptr,      1);
    add_argument_d(argc,argv,"-bin", &bin_width,                 "Bin size"                                , cl_tags, nullptr,      1);
    conclude_input_arguments(argc,argv,program_name,cl_tags);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension("-d",in_file_name,".dat");
    check_extension("-o",out_file_name,".dat");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read in data                                                                                              //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Data_file data;
    data.set_headers(0);
    data.get_data(in_file_name);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Copy data to vector                                                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    dv1d my_data(0,0.0);
    for(i=0; i<data.data_d.size(); i++) //loop over data
    {
        if(strcmp(data.data_s[i][col].c_str(), "NaN") == 0)
        {
        }
        else
        {
            my_data.push_back(data.data_d[i][col]);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Make histogram                                                                                            //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Histogram_d histo;
    histo.bin_data(my_data,bin_width);
    histo.write_histo(out_file_name,"binned data (same units as the provided data)");

    std::cout << "\nFormatting completed successfully" << "\n\n";

    return 0;
}



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
#include "../headers/index.h"  
#include "../headers/param.h"

int main(int argc, const char * argv[]) 
{
    //Here we define some variables used throughout
    string in_file_name;                   //Name of input file with data to be averaged
    string out_file_name;                  //Name of output file
    int i      = 0;                         //General variable used in loops
    int j      = 0;                         //General variable used in loops
    int k      = 0;                         //General variable used in loops
    int range  = 0;                         //How many frames on each side
    int col    = 0;                         //Which column holds the target data
    int header = 0;                         //How many header lines does the file contain
    sv1d cl_tags;                           //Holds a list of command line tags for the program
   
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set program name/description and print info                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    string program_name = "Running Average";

    print_credits(argc,argv,program_name);

    string program_description = "Running Average is an analysis tool used for averaging the contents of a timeline data over a window around each data point.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments(argc,argv,program_description);
    add_argument_s(argc,argv,"-d"     , in_file_name,               "Input file with data to average (dat)"       , cl_tags, nullptr,      1);
    add_argument_s(argc,argv,"-o"     , out_file_name,              "Output data file "                           , cl_tags, nullptr,      1);
    add_argument_i(argc,argv,"-range" , &range,                     "How many frames on each to average over?"    , cl_tags, nullptr,      1);
    add_argument_i(argc,argv,"-col"   , &col,                       "Which column holds target data to average?"  , cl_tags, nullptr,      1);
    add_argument_i(argc,argv,"-header", &header,                    "How many header lines does the file contain?", cl_tags, nullptr,      1);
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
    // Create data files                                                                                         //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Data_file data;
    data.set_headers(header);
    data.get_data(in_file_name.c_str());

    Data_file data_avg;
    data_avg.set_headers(header);
    data_avg.get_data(in_file_name.c_str());

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Compute running average                                                                                   //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(i=range; i<data.size_y()-range; i++) //loop over rows
    {
        double avg = 0.0;
        for(j=i-range; j<=i+range; j++) //loop over window
        {
            avg = avg + data.data_d[j][col];
        }
        avg = avg/(double)(2*range + 1);

       data_avg.data_d[i][col] = avg; 
    }

    data_avg.double_to_string();

    //take care of frames outside the window
    for(i=0; i<range; i++)
    {
        data_avg.data_s[i][col] = "NaN";
    }
    for(i=data.size_y()-range; i<data.size_y(); i++)
    {
        data_avg.data_s[i][col] = "NaN";
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


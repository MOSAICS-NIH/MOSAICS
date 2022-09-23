
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
    int i      = 0;               //General variable used in loops
    int j      = 0;               //General variable used in loops
    int height = 0;               //How many rows in pdb
    int width  = 0;               //How many columns in pdb
    string row;                   //Name the atom selection 
    string out_file_name;         //Name of file to write commands to
    string end_text;              //Any ending text to print

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set program name/description and print info                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    string program_name = "Mean Coords Row Selector";

    print_credits(argc,argv,program_name);

    string program_description = "Mean Coords Row Selector prints PyMOL commands letting the user select the rows and collumns of the lipids produced by Mean Lipid Coords.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments(argc,argv,program_description);
    add_argument_i(argc,argv,"-height",&height,                    "How many rows in the pdb?"                         , nullptr,      1);
    add_argument_i(argc,argv,"-width", &width,                     "How many columnds in the pdb?"                     , nullptr,      1);
    add_argument_s(argc,argv,"-o",     out_file_name,              "Name of output file with selection commands (pml)" , nullptr,      1);
    add_argument_s(argc,argv,"-e",     end_text,                   "Any selection text to be added (string)"           , nullptr,      0);
    conclude_input_arguments(argc,argv,program_name);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension("-o",out_file_name,".pml");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Print the pymol command                                                                                   //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    FILE *out_file = fopen(out_file_name.c_str(),"a");

    fprintf(out_file,"#columns \n");

    //print columns
    for(i=0; i<width; i++) 
    {   
        string col = "col_" + to_string(i+1);

        fprintf(out_file,"select %s, resi ",col.c_str());

        for(j=0; j<height; j++)
        {
            int lipid = (i) + j*width;

            fprintf(out_file,"%d",lipid);
            if(j < height-1)
            {
                fprintf(out_file,"+");
            }
            else
            {
                fprintf(out_file," %s \n",end_text.c_str());
            }
        }
    }            

    //print rows
    fprintf(out_file,"#rows \n");
    for(i=0; i<height; i++) 
    {
        string row = "row_" + to_string(i+1);

        fprintf(out_file,"select %s, resi ",row.c_str());

        for(j=0; j<width; j++)
        {
            int lipid = j + i*width;

            fprintf(out_file,"%d",lipid);
            if(j < width-1)
            {
                fprintf(out_file,"+");
            }
            else
            {
                fprintf(out_file," %s \n",end_text.c_str());
            }
        }
    }

    fprintf(out_file,"\n");
    fclose(out_file);

    std::cout << "\nFormatting completed successfully" << "\n\n";

    return 0;
}


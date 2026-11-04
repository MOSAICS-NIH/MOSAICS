
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
#include "../headers/file_io_variable.h"

int main(int argc, const char * argv[]) 
{
    //Here we define some variables used throughout
    string in_file_name_1;                 //Name of input file data file 1
    string in_file_name_2;                 //Name of input file data file 2
    string out_file_name;                  //Name of output file
    int i        = 0;                      //General variable used in loops
    int j        = 0;                      //General variable used in loops
    int print    = 1;                      //Tells if current line should be printed
    int remove   = 0;                      //Tells is the contact should be removed due to a cutoff
    int b_screen = 0;                      //Tells if we are removing duplicate bonds
    int mode     = 0;                      //Mode when -d2 is provided
    double scale = 1.0;                    //How much to scale dash radius by
    double cutoff = 0.0;                   //Cutoff value for removing contacts
    sv1d cl_tags;                          //Holds a list of command line tags for the program
   
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set program name/description and print info                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    string program_name = "Contact Editor";

    print_credits(argc,argv,program_name);

    string program_description = "Contact Editor is an analysis tool used for scanning PyMOL files generated with Contact Kinetics ";
    program_description = program_description + "or H-bond Kinetics and removing duplicate contacts, scaling dash radius, removing infrequent ";
    program_description = program_description + "contacts, or transfering contact frequency data. ";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments(argc,argv,program_description);
    add_argument_s(argc,argv,"-d1"    , in_file_name_1,             "Contacts file to remove duplicates from (pml) "                   , cl_tags, nullptr,      1);
    add_argument_s(argc,argv,"-d2"    , in_file_name_2,             "Contacts file to compare with (pml) "                             , cl_tags, &b_screen,    0);
    add_argument_s(argc,argv,"-o"     , out_file_name,              "Output contacts file (pml)"                                       , cl_tags, nullptr,      1);
    add_argument_d(argc,argv,"-scale" , &scale,                     "Scaling factor for dash radious"                                  , cl_tags, nullptr,      0);
    add_argument_d(argc,argv,"-cutoff", &cutoff,                    "Cutoff value for removing contacts"                               , cl_tags, nullptr,      0);
    add_argument_i(argc,argv,"-mode" ,  &mode,                      "What mode if -d2 is provided? (0:rm duplicates, 1:transfer freq)" , cl_tags, nullptr,      0);
    conclude_input_arguments(argc,argv,program_name,cl_tags);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension("-d1",in_file_name_1,".pml");
    check_extension("-o",out_file_name,".pml");
    if(b_screen == 1)
    {
        check_extension("-d2",in_file_name_2,".pml");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Create data files                                                                                         //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Data_file_vari data_1;
    Data_file_vari data_2;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read data files                                                                                           //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    data_1.get_data(in_file_name_1);
    if(b_screen == 1)
    { 
        data_2.get_data(in_file_name_2);
    }

    FILE *out_file = fopen(out_file_name.c_str(),"w");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // scale dash radius and Check for and remove duplicate bonds in the data                                    //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int cont_index = 0;
    for(i=0; i<data_1.data_s.size(); i++) //loop over file 1 rows
    {
        if(data_1.data_s[i].size() > 0) //skip empty lines
        {
            //scale dash radious
            if(data_1.data_s[i].size() > 1) //must have at least 2 lines
            {
                if(strcmp(data_1.data_s[i][1].c_str(), "dash_radius,") == 0)
                {
    			data_1.data_s[i][2] = to_string(scale*atof(data_1.data_s[i][2].c_str())) + ",";
                }
            }

            if(b_screen == 1)
            {
                if(strcmp(data_1.data_s[i][0].c_str(), "Contact") == 0)
                {
                    int duplicate = 0; 

                    int atom_a = stoi(data_1.data_s[i][12]);
                    int atom_b = stoi(data_1.data_s[i][15]);

                    for(j=0; j<data_2.data_s.size(); j++) //loop over file 2 rows
                    {
                        if(data_2.data_s[j].size() > 0) //skip empty lines
                        {
                            if(strcmp(data_2.data_s[j][0].c_str(), "Contact") == 0)
                            {
                                int cmp_atom_a = stoi(data_2.data_s[j][12]);
                                int cmp_atom_b = stoi(data_2.data_s[j][15]);

                                if( (atom_a == cmp_atom_a && atom_b == cmp_atom_b) || (atom_a == cmp_atom_b && atom_b == cmp_atom_a) )
                                {
                                    duplicate = 1; 
                                    printf("duplicate found \n");
                      
                                    if(mode == 1) //
                                    {
                                        printf("transfering frequency data for atoms %d %d: data_1 %s data_2 %s \n",cmp_atom_a,cmp_atom_b,data_1.data_s[i][9].c_str(),data_2.data_s[j][9].c_str());
                                        data_1.data_s[i][9] = data_2.data_s[j][9];
                                        print = 1;
                                    }
                                }                 
                            }
                        }
                    }

	            if(duplicate == 0)
                    {
                        if(mode == 0) //remove duplicate contacts
                        {
                            print = 1;

	                    //renumber contacts 
                            data_1.data_s[i+0][1] = to_string(cont_index) + ":";
                            data_1.data_s[i+3][1] = to_string(cont_index) + ",";
                            data_1.data_s[i+4][1] = "dist_" + to_string(cont_index) + ",";
                            data_1.data_s[i+5][3] = "dist_" + to_string(cont_index);
	                    cont_index++;
                        }
                    }
                    else 
                    {
                        if(mode == 0)
                        {
                            print = 0;
                        }
                    }
                }
            }

            //remove contacts based on cutoff
            if(strcmp(data_1.data_s[i][0].c_str(), "Contact") == 0)
            {
                if(atof(data_1.data_s[i][9].c_str()) >= cutoff)
                {
                    remove = 0; 
                }
                else
                {
                    remove = 1; 
                }
            } 
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Report data                                                                                               //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(print == 1 && remove == 0)
        {
            for(j=0; j<data_1.data_s[i].size(); j++)
            {
                fprintf(out_file,"%s ",data_1.data_s[i][j].c_str());
            }
            fprintf(out_file,"\n");
        }
        else //print to screen data removed
        {
            for(j=0; j<data_1.data_s[i].size(); j++)
            {
                printf("%s ",data_1.data_s[i][j].c_str());
            }
            printf("\n");
        }
    }
    fclose(out_file);

    std::cout << "\nFormatting completed successfully" << "\n\n";

    return 0;
}


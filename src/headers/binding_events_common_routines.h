
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <vector>

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes the grid to an output file.                                                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_grid_to_file(int num_g_x,int num_g_y,iv2d &nan,string my_file_name,dv2d &data)
{
    int j=0;
    int k=0;

    FILE *my_file = fopen(my_file_name.c_str(), "w");
    if(my_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",my_file_name.c_str());
    }

    for(k=0; k<num_g_y; k++) //loop over y-dimension
    {
        for(j=0; j<num_g_x; j++) //loop over x-dimension
        {
                if(nan[k][j] == 0)
                {
                    fprintf(my_file," %10.6f",data[k][j]);
                }
                else //data excluded
                {
                    fprintf(my_file," %10s ","NaN");
                }
        }
        fprintf(my_file,"\n");
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This functionc reads a binding events file                                                                //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int get_binding_events(string in_file_name,iv1d &lipid_nr,iv1d &res_nr,sv1d &res_name,iv1d &bind_i,iv1d &bind_f,
                        iv1d &dwell_t,int *x_i,int *y_i,double *ef_dt,int *ef_frames,int *num_lipids,int *num_g_x,
                        int *num_g_y,double *APS)
{
    int number_of_lines = 0;      //Number of lines in input files
    int items_per_line  = 0;      //How many item is a single line
    int k               = 0;      //Standard variable used in loops
    int l               = 0;      //Standard variable used in loops
    int m               = 0;      //Standard variable used in loops
    int outcome         = 0;      //Return whether the file was read succesfully or not
    char my_string[200];          //String to hold read in data entries

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Here we open files for reading/writing                                                                    //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    FILE *in_file = fopen(in_file_name.c_str(), "r");
    if(in_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",in_file_name.c_str());
    }
    else
    {
        outcome = 1;

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Get information about the data files                                                                      //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        characterize_file(in_file_name,&number_of_lines,&items_per_line,4);

        //printf("in_file_name %20s number_of_lines %10d items_per_line %10d \n",in_file_name.c_str(),number_of_lines,items_per_line);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Read in binding events                                                                                    //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        for(k=0; k<number_of_lines+4; k++) //loop over rows
        {
            if(k==0) //first line containts x_i,y_i,ef_dt,ef_frames,num_lipids,num_g_x and num_g_y
            {
                for(l=0; l<16; l++) //loop over header items
                {
                    for(m=0; m<200; m++)
                    {
                        my_string[m] = ' ';
                    }
                    fscanf(in_file, "%s,", my_string);

                    if(l == 1) //x_i
                    {
                        *x_i = atoi(my_string);
                    }
                    if(l == 3) //y_i
                    {
                        *y_i = atoi(my_string);
                    }
                    if(l == 5) //ef_dt
                    {
                        *ef_dt = atof(my_string);
                    }
                    if(l == 7) //ef_frames
                    {
                        *ef_frames = atoi(my_string);
                    }
                    if(l == 9) //num_lipids
                    {
                        *num_lipids = atoi(my_string);
                    }
                    if(l == 11) //num_g_x
                    {
                        *num_g_x = atoi(my_string);
                    }
                    if(l == 13) //num_g_y
                    {
                        *num_g_y = atoi(my_string);
                    }
                    if(l == 15) //APS
                    {
                        *APS = atof(my_string);
                        next_line(in_file);
                    }
                }
            }
            else if(k>3) //first 4 lines are header information
            {
                for(l=0; l<items_per_line; l++) //loop over collumns
                {
                    for(m=0; m<200; m++)
                    {
                        my_string[m] = ' ';
                    }
                    fscanf(in_file, "%s,", my_string);

                    if(l == 0) //lipid number
                    {
                        lipid_nr.push_back(atoi(my_string));
                    }
                    if(l == 1) //res number
                    {
                        res_nr.push_back(atoi(my_string));
                    }
                    if(l == 2) //res name
                    {
                        res_name.push_back(strdup(my_string));
                    }
                    if(l == 3) //bind_i
                    {
                        bind_i.push_back(atoi(my_string));
                    }
                    if(l == 4) //bind_f
                    {
                        bind_f.push_back(atoi(my_string));
                    }
                    if(l == 5) //dwell_t
                    {
                        dwell_t.push_back(atoi(my_string));
                    }
                }
            }
            else //read the line
            {
                next_line(in_file);
            }
        }
        fclose(in_file);
    }
    return outcome;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function takes in binding evets vectors and fills binding timeline vectors                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_binding_timeline(int num_lipids,int ef_frames,iv1d &lipid_nr,iv1d &res_nr,sv1d &res_name,iv1d &bind_i,iv1d &bind_f,
                          iv1d &time_line_lipid_nr,iv1d &time_line_res_nr,sv1d &time_line_res_name,iv2d &bound_time_line)
{
    int i = 0;
    int j = 0;

    //initialize the time line
    for(i=0; i<num_lipids; i++) //loop over lipids
    {
        for(j=0; j<ef_frames; j++) //loop over frames
        {
            bound_time_line[j][i] = 0;
        }
    }

    //now add bound lipids to the time line
    for(i=0; i<res_nr.size(); i++) //loop over binding events
    {
        for(j=bind_i[i]; j<=bind_f[i]; j++) //loop over bound frames
        {
            time_line_lipid_nr[lipid_nr[i]] = lipid_nr[i];
            time_line_res_nr[lipid_nr[i]]   = res_nr[i];
            time_line_res_name[lipid_nr[i]] = res_name[i];

            bound_time_line[j][lipid_nr[i]] = 1;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function takes in binding evets vectors and adds to a binding timeline vector                        //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void add_to_binding_timeline(int num_lipids,int ef_frames,iv1d &lipid_nr,iv1d &res_nr,sv1d &res_name,iv1d &bind_i,iv1d &bind_f,
                             iv1d &time_line_lipid_nr,iv1d &time_line_res_nr,sv1d &time_line_res_name,iv2d &bound_time_line)
{
    int i = 0;
    int j = 0;

    for(i=0; i<res_nr.size(); i++) //loop over binding events
    {
        for(j=bind_i[i]; j<=bind_f[i]; j++) //loop over bound frames
        {
            time_line_lipid_nr[lipid_nr[i]] = lipid_nr[i];
            time_line_res_nr[lipid_nr[i]]   = res_nr[i];
            time_line_res_name[lipid_nr[i]] = res_name[i];

            bound_time_line[j][lipid_nr[i]] = 1;
        }
    }
}



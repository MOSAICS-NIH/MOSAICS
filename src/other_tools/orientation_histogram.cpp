
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
    FILE *out_file;               //File for writing data
    string in_file_name;          //Name of input file 
    string out_file_name;         //Name of output file
    int i                = 0;     //General variable used in loops
    int j                = 0;     //General variable used in loops
    int k                = 0;     //General variable used in loops
    int l                = 0;     //General variable used in loops
    int next             = 0;     //Used for controlling output printing
    double res_theta     = 1;     //bin size for theta
    double res_phi       = 1;     //bin size for phi
    double p_res_theta   = 0.005; //printing resultion
    double biggest_phi   = 0;     //biggest value of phi
    double biggest_theta = 0;     //biggest value of theta
    double percent       = 0;     //percent of data in bin for theta/phi
    double line_value    = 1;     //what value to print for theta contour lines

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set program name/description and print info                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    string program_name = "Orientation Histogram";

    print_credits(argc,argv,program_name);

    string program_description = "Orientation Histogram is an analysis tool used for representing protein tilt angle data as a polar histogram.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments(argc,argv,program_description);
    add_argument_s(argc,argv,"-d",         in_file_name,               "Input data file"                         , nullptr,      1);
    add_argument_s(argc,argv,"-o",         out_file_name,              "Output file"                             , nullptr,      1);
    add_argument_d(argc,argv,"-res_t",     &res_theta,                 "Bin size for theta"                      , nullptr,      1);
    add_argument_d(argc,argv,"-res_p",     &res_phi,                   "Bin size for phi"                        , nullptr,      1);
    add_argument_d(argc,argv,"-res",       &p_res_theta,               "Printing resolution"                     , nullptr,      1);
    add_argument_d(argc,argv,"-line_v",    &line_value,                "What value to use for contour lines"     , nullptr,      1);
    conclude_input_arguments(argc,argv,program_name);

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
    Data_file angles;
    angles.set_headers(1);
    angles.get_data(in_file_name);
    //angles.show_data();
 
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Determine the biggest phi and theta                                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(i=0; i<angles.data_d.size(); i++) //loop over angles
    {
        if(angles.data_d[i][2] > biggest_phi)
        {
            biggest_phi = angles.data_d[i][2];
        }
        if(angles.data_d[i][1] > biggest_theta)
        {
            biggest_theta = angles.data_d[i][1];
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Determine how frequent each theta and phi values are                                                      //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int num_bins_theta = ceil((biggest_theta)/res_theta);    //How many bins are theta values grouped into
    int num_bins_phi = ceil((360)/res_phi);                  //How many bins are phi values grouped into
    //printf("num_bins_theta %10d num_bins_phi %10d res_theta %10.3f res_phi %10.3f biggest_theta %10.3f biggest_phi %10.3f \n",num_bins_theta,num_bins_phi,res_theta,res_phi,biggest_theta,biggest_phi);

    iv2d count(num_bins_theta,iv1d(num_bins_phi,0));

    //now we read through the data and count how frequent each theta/phi is
    for(i=0; i<angles.data_d.size(); i++) //loop over angles
    {
        for(k=0; k<num_bins_theta; k++) //loop over bins theta
        {
            for(l=0; l<num_bins_phi; l++) //loop over bins phi
            {
                if(angles.data_d[i][1] >= k*res_theta && angles.data_d[i][1] < (k+1)*res_theta && angles.data_d[i][2] >= l*res_phi && angles.data_d[i][2] < (l+1)*res_phi)
                {
                    count[k][l] = count[k][l] + 1;
                }
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Convert the counts into percentages                                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int total = 0;      //this is the toal number frames in traj
    for(i=0; i<num_bins_theta; i++) //loop over bins theta
    {
        for(j=0; j<num_bins_phi; j++) //loop over bins phi
        {
            total = total + count[i][j];
        }
    }

    for(i=0; i<num_bins_theta; i++) //loop over bins theta
    {
        for(j=0; j<num_bins_phi; j++) //loop over bins phi
        {
            percent = 100*(double)count[i][j]/(double)total;
            //printf(" %10d %10d %10.3f \n",i,count[i][j],percent);
        }
    }   

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Print data to the output file                                                                             //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int max_theta      = ceil(biggest_theta);  //the is the maximum value of theta from the data
    int num_cont_lines = (int)max_theta;       //the number of contour lines indicating the value of theta
    //printf("num_cont_lines %d max_theta %d \n",num_cont_lines,max_theta);

    int num_bins_theta_p = ceil(((double)max_theta)/p_res_theta);  //how many printing bins are needed to go from zero to theta max. 

    out_file = fopen(out_file_name.c_str(), "w");

    for(j=-num_bins_theta_p; j<=num_bins_theta_p; j++) //loop over y direction
    {
        for(i=-num_bins_theta_p; i<=num_bins_theta_p; i++) //loop over x direction
        {
            if(i == 0 && j == 0) //print a dot at the origin
            {
                fprintf(out_file,"%10.6f ",line_value);
            }
            else
            {
                double x    = i*p_res_theta;   //distance in x from origin
                double y    = j*p_res_theta;   //distnace in y from origin
                double dist = sqrt(x*x + y*y); //distance from origin. represents the theta value on polar plot
               
                //measure angle from [0,1]. use cos_angle = v_dot_y/(mag_v*mag_y)
                double y_vec[2]  = {0, 1};
                double mag_y     = 1;
                double mag_v     = dist;
                double v_dot_y   = x*0 + y*1;
                double cos_angle = v_dot_y/(mag_v*mag_y);
                double pi        = 3.141592654;
                double angle     = acos(cos_angle)*(180/pi);    //angle from y-axis in degrees. represents phi on polar plot

                if(i < 0) //convert to a 0-360 scale
                {
                    angle = 2*(180 - angle) + angle;
                }

                //given we have a distance (theta value) and an angle we must find the bins index that have these values             
                int phi_i   = ceil(angle/res_phi) - 1;        //index for bin holding phi = angle
                int theta_i = ceil(dist/res_theta) - 1;       //index for bin holding theta = distance
                if(angle == 0) //needed because ceil(0) = 0
                {
                    phi_i = 0;
                }

                //put lines indicating the theta value
                for(k=0; k<max_theta; k++) //loop over theta values
                {
                    if(dist > (k - 0.05) && dist < (k + 0.05) ) //gives the lines a thickness of 0.1
                    {
                        fprintf(out_file,"%10.6f ",line_value);
                        next = 1;
                    }
                }

                //write the percentage of frames having theta and phi values
                if(next == 0)
                {
                    if(phi_i < num_bins_phi && theta_i < num_bins_theta)
                    {
                        if(count[theta_i][phi_i] > 0)
                        { 
                            fprintf(out_file,"%10.6f ",100*(double)count[theta_i][phi_i]/(double)total);
                        }
                        else 
                        {
                            fprintf(out_file,"%10s ","NaN");
                        }
                    }
                    else
                    {
                        fprintf(out_file,"%10s ","NaN");
                    }
                }
                next = 0;
            }
        }
        fprintf(out_file,"\n");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Close files and end program                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    fclose(out_file);

    std::cout << "\nFormatting completed successfully" << "\n\n";

    return 0;
}


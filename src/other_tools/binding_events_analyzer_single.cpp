
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <sstream>
#include <mpi.h>
#include <fstream>

using namespace std;

#include "../headers/multi_dim_vec.h"
#include "../headers/file_reader.h"
#include "../headers/vector_mpi.h"
#include "../headers/switch.h"
#include "../headers/common_routines.h"
#include "../headers/common_routines_mpi.h"
#include "../headers/file_naming.h"
#include "../headers/binding_events.h"
#include "../headers/fit.h"
#include "../headers/index.h"
#include "../headers/command_line_args_mpi.h"
#include "../headers/array.h"
#include "../headers/force_serial.h"
#include "../headers/file_naming_mpi.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This is the main function of the program.                                                                 //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, const char * argv[]) 
{
    //Here we define some variables used throughout
    FILE *out_file;                    //File for writing frequencies to
    string binding_events_file_name;   //Name of binding events file 
    string out_file_name;              //Name of the output file
    string lip_t_file_name;            //Name of the lipid types file
    int i            = 0;              //General variable used in loops
    int j            = 0;              //General variable used in loops
    int k            = 0;              //General variable used in loops
    int x            = 0;              //Grid point in x-direction
    int y            = 0;              //Grid point in y-direction
    int b_x          = 0;              //Was a grid point in x-direction specified
    int b_y          = 0;              //Was a grid point in y-direction specified
    int result       = 0;              //result of reading in binding events file
    int b_write_hist = 0;              //Write frequencies?
    int world_size   = 0;              //Size of the mpi world
    int world_rank   = 0;              //Rank in the mpi world
    int threshold    = 0;              //Cutoff for mending fragmented binding events
    int report       = 0;              //Report a list of the binding events?
    double slope     = 0.0;            //slope of LnP vs time
    double yint      = 0.0;            //ying of LnP vs time
    double r2        = 0.0;            //Correlation coeficient in linear regression
    double koff      = 0.0;            //The koff value
    double cutoff    = 0.0;            //Exclude data with cutoff less than this
    double bin_width = 1.0;            //bin width (default is 1ps)
    double percent   = 0.0;            //Percent reported for the bin
    sv1d cl_tags;                      //Holds a list of command line tags for the program

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set up the mpi environment                                                                                //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    MPI_Init(NULL, NULL);;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);      //get the world size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);      //get the process rank

    //force program to run in serial?
    enum Switch serial         = on;

    //here we check if the program supports parallelization or not
    check_serial(world_rank,world_size,serial);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set program name/description and print info                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    string program_name = "Binding Events Analyzer Single";

    print_credits(argc,argv,program_name);

    string program_description = "Binding Events Analyzer Single is an analysis tool used for reading a binding events file and computing the average dwell time, and koff, etc.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments_mpi(argc,argv,world_rank,program_description);
    add_argument_mpi_s(argc,argv,"-d"        , binding_events_file_name,"Input binding events file (be)"                               , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-histo"    , out_file_name,           "Output data file with dwell time histogram (dat)"             , world_rank, cl_tags, &b_write_hist,0);
    add_argument_mpi_s(argc,argv,"-crd"      , lip_t_file_name,         "Selection card with lipid types (crd)"                        , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-x"        , &x,                      "Grid point in x-direction"                                    , world_rank, cl_tags, &b_x,         0);
    add_argument_mpi_i(argc,argv,"-y"        , &y,                      "Grid point in y-direction"                                    , world_rank, cl_tags, &b_y,         0);
    add_argument_mpi_d(argc,argv,"-cutoff"   , &cutoff,                 "Exclude data with a dwell time smaller than this (ps)"        , world_rank, cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-bin"      , &bin_width,              "Bin width (ps)"                                               , world_rank, cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-repair"   , &threshold,              "Maximum allowed size (frames) for mending fragmented events"  , world_rank, cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-report"   , &report,                 "Print a list of the dwell times (0:no, 1:yes)"                , world_rank, cl_tags, nullptr,      0);
    conclude_input_arguments_mpi(argc,argv,world_rank,program_name,cl_tags);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension_mpi(world_rank,"-d",binding_events_file_name,".be");
    check_extension_mpi(world_rank,"-crd",lip_t_file_name,".crd");

    if(b_write_hist == 1)
    {
        check_extension_mpi(world_rank,"-histo",out_file_name,".dat");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read in lipid types                                                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Index lip_t;
    lip_t.get_index(lip_t_file_name);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check if a lattice point was specified                                                                    //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(world_rank == 0)
    {
        if(b_x == 1 && b_y == 0)
        {
            printf("A lattice point was specified for the x-direction but not y. Please include the y-direction if analyzing binding events for a lattice point. \n");
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }
        else if(b_x == 0 && b_y == 1)
        {
            printf("A lattice point was specified for the y-direction but not x. Please include the x-direction if analyzing binding events for a lattice point. \n");
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read in binding events                                                                                    //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Binding_events events;
  
    if(b_x==1 && b_y==1)
    {
        result = events.get_info(binding_events_file_name);
        if(result == 1)
        {
            result = events.get_binding_events_xy(binding_events_file_name,x,y);
        }
    }
    else
    {
        result = events.get_binding_events_bin(binding_events_file_name);
    }

    if(result == 1) //bind events file exists
    { 
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Mend any fragmented events                                                                                //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        events.suppress_timeline_noise_be_alt(threshold);  //mend fragmented events

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // remove events with dwell time shorter than cutoff_dt                                                      //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        events.screen_events(cutoff);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Compute average dwell time                                                                                //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double dwell_time = 0.0;
        int ef_num_events = 0;

        for(i=0; i<events.dwell_t.size(); i++) //loop over binding events
        {
            for(j=0; j<lip_t.index_s.size(); j++) //loop over lipid types
            {
                if(strcmp(events.res_name[i].c_str(), lip_t.index_s[j].c_str()) == 0) //lipid type is correct 
                {
                    dwell_time    = dwell_time + (double)events.dwell_t[i];  
                    ef_num_events = ef_num_events + 1;

                    if(report == 1)
                    {
                        printf("Dwell time %5d: %f (ps)\n",ef_num_events,(double)events.dwell_t[i]*events.ef_dt);
                    }
                }
            }
        }
        dwell_time = dwell_time/(double)ef_num_events;

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Compute variance in dwell time                                                                            //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double stdev         = 0.0;
        double sum_o_squares = 0.0;

        for(i=0; i<events.dwell_t.size(); i++) //loop over binding events
        {
            for(j=0; j<lip_t.index_s.size(); j++) //loop over lipid types
            {
                if(strcmp(events.res_name[i].c_str(), lip_t.index_s[j].c_str()) == 0) //lipid type is correct
                {
                    sum_o_squares = sum_o_squares + pow((double)events.dwell_t[i]-dwell_time,2);
                }
            }   
        }
        stdev = sqrt(sum_o_squares/(double)(ef_num_events-1));

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Compute effective largest dwell_t                                                                         //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        int largest_dwell_t = 0;
        for(i=0; i<events.dwell_t.size(); i++) //loop over binding events
        {
            for(j=0; j<lip_t.index_s.size(); j++) //loop over lipid types
            {
                if(strcmp(events.res_name[i].c_str(), lip_t.index_s[j].c_str()) == 0) //lipid type is correct
                {
                    if(events.dwell_t[i] > largest_dwell_t)
                    {
                        largest_dwell_t = events.dwell_t[i];
                    }
                }
            }
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Compute frequencies                                                                                       //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double dwell_time_freq[largest_dwell_t];       //holds how many binding events had dwell time i
        init_darray(dwell_time_freq,largest_dwell_t);  

        for(i=0; i<largest_dwell_t; i++) //loop over dwell times 
        {
            for(j=0; j<events.dwell_t.size(); j++) //loop over binding events
            {
                if(events.dwell_t[j] == i+1) //dwell time matches
                {
                    for(k=0; k<lip_t.index_s.size(); k++) //loop over lipid types
                    {
                        if(strcmp(events.res_name[j].c_str(), lip_t.index_s[k].c_str()) == 0) //lipid type is correct
                        {
                            dwell_time_freq[i] = dwell_time_freq[i] + 1;
                        }
                    }
                }
            }
            dwell_time_freq[i] = dwell_time_freq[i]/(double)ef_num_events;   
            //printf("dwell_time_freq %10.3f \n",dwell_time_freq[i]);
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Compute probability of lasting t or longer. (add probabilities for t and longer)                          //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double nan = -999999;       //assign a negative value to freq for frequenies with no data  
        int count = 0;              //how many dwell times have a positive frequency
        for(i=0; i<largest_dwell_t; i++) //loop over dwell times
        {
            if(dwell_time_freq[i] > 0) //inclue data in fit
            {
                for(j=i+1; j<largest_dwell_t; j++) //add probabilities
                {
                    dwell_time_freq[i] = dwell_time_freq[i] + dwell_time_freq[j];
                }
                count++;
            }
            else //exclue data from fit
            {
                dwell_time_freq[i] = nan;
            }
        } 

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Get Ln(freq) and time                                                                                     //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //create array to hold ln(frequencies) excluding Ln(0)
        double log_freq[count];         //holds the log freq for dwell times having positive freq
        double time[count];             //holds time for dwell times having positive freq

        //reset count
        count = 0;

        //fill array excluding Ln(0)
        for(i=0; i<largest_dwell_t; i++) //loop over dwell times
        {
            if(dwell_time_freq[i] != nan) //dont include dwell times with no samples
            {
                log_freq[count] = log(dwell_time_freq[i]);
                time[count]     = (i+1)*events.ef_dt;
                count++;
            }
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Compute koff                                                                                              //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //do least squares regression
        least_squares_regression(count,time,log_freq,&slope,&yint,&r2);
        koff = -slope; 

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Print output                                                                                              //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        printf(" %25s %-10f \n","effective dt(ps)",events.ef_dt);
        printf(" %25s %-10f \n","largest dwell time(ps)",(double)largest_dwell_t*events.ef_dt);
        printf(" %25s %-10d \n","binding events",ef_num_events);
        printf(" %25s %-10f +/- %10f \n","average dwell time(ps)",dwell_time*events.ef_dt,stdev*events.ef_dt);
        printf(" %25s %-10f \n","r^2",r2);
        printf(" %25s %-10f \n","koff(ps^-1)",koff);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Here we have code for making a probability histogram                                                      //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(b_write_hist == 1) //write log freq to file
        {
            printf("\nBinning the data and writting a histogram to %s. \n",out_file_name.c_str());

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                                                           //
            // Bin the data                                                                                              //
            //                                                                                                           //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            int num_bins = ceil((events.ef_dt*(double)largest_dwell_t)/bin_width) + 1;
            iv1d bin_count(num_bins,0);
            for(i=0; i<events.dwell_t.size(); i++) //loop over data
            {
                for(j=0; j<lip_t.index_s.size(); j++) //loop over lipid types
                {
                    if(strcmp(events.res_name[i].c_str(), lip_t.index_s[j].c_str()) == 0) //lipid type is correct
                    {
                        for(k=0; k<num_bins; k++)
                        {
                            if(events.ef_dt*(double)events.dwell_t[i] >= (double)k*bin_width && events.ef_dt*(double)events.dwell_t[i] < (double)(k+1)*bin_width)
                            {
                                bin_count[k] = bin_count[k] + 1;
                            }
                        }
                    }
                }
            }

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                                                           //
            // Normalize histogram and write data to output file                                                         //
            //                                                                                                           //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            out_file = fopen(out_file_name.c_str(), "w");
            if(out_file == NULL)
            {
                printf("failure opening %s. Make sure the file exists. \n",out_file_name.c_str());
            }
            else
            {
                double sum_p = 0;
                fprintf(out_file," %10s %10s %10s %10s \n","#bin","time(ps)","count","probability");
                for(i=0; i<num_bins; i++)
                {
                    percent = (double)bin_count[i]/(double)ef_num_events;
                    sum_p = sum_p + percent;
                    fprintf(out_file," %10d %10.3f %10d %10.4f \n",i,(double)i*bin_width,bin_count[i],percent);
                }
                fclose(out_file);
                printf("probability sums to %f. \n",sum_p);
            }
        }
    }
    else //binding events file does not exits
    {
        printf("Could not find binding events file \n"); 
    }

    //relinquish the mpi environment
    MPI_Finalize();

    std::cout << "\nFormatting completed successfully" << "\n\n";

    return 0;
}


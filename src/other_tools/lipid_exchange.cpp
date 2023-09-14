
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
#include "../headers/binding_events.h"
#include "../headers/fit.h"
#include "../headers/index.h"
#include "../headers/histo.h"
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
    string binding_events_file_name;   //Name of binding events file 
    string leaving_lipid_file_name;    //Lipid type leaving
    string replacing_lipid_file_name;  //Lipid type replacing
    string histo_file_name;            //Name of the histogram file
    int i                 = 0;         //General variable used in loops
    int j                 = 0;         //General variable used in loops
    int k                 = 0;         //General variable used in loops
    int bin_width         = 1;         //The bin width 
    int b_histo           = 0;         //Make histogram for the time between the leaving lipid and the entering lipid?
    int b_min             = 0;         //Did the user specify a minimum range for the histogram?
    int b_max             = 0;         //Did the user specify a maximum range for the histogram?
    int min               = 0;         //minimum value of histogram
    int max               = 0;         //maximum value of histogram
    int world_size        = 0;         //Size of the mpi world
    int world_rank        = 0;         //Rank in the mpi world
    int b_self            = 1;         //Count events when leaving and replacing lipids are one         
    double lipid_fraction = 0;         //What percentage of the lipids are the target type
    double slope          = 0;         //slope of LnP vs time
    double yint           = 0;         //ying of LnP vs time
    double r2             = 0;         //Correlation coeficient in linear regression
    double koff           = 0;         //The koff value
    double cutoff         = 0;         //Exclude data with cutoff less than this
    sv1d cl_tags;                      //Holds a list of command line tags for the program

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set up the mpi environment                                                                                //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    MPI_Init(NULL, NULL);
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
    string program_name = "Lipid Exchange";

    print_credits(argc,argv,program_name);

    string program_description = "Lipid Exchange is an analysis tool used for reading a binding events file containing multiple lipid types and computing the lipid exchange frequency.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments_mpi(argc,argv,world_rank,program_description);
    add_argument_mpi_s(argc,argv,"-d"        , binding_events_file_name, "Input binding events file (be)"                                                           , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-l_in"     , replacing_lipid_file_name,"Selection card with replacing lipid types (crd)"                                          , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-l_out"    , leaving_lipid_file_name,  "Selection card with outgoing lipid types (crd)"                                           , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-l_frac"   , &lipid_fraction,          "Fraction of lipids in the target leaflet that are of the type given in -l_in (0-1)"       , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-cutoff"   , &cutoff,                  "Exclude binding events with a dwell time smaller than this (ps)"                          , world_rank, cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-histo"    , histo_file_name,          "Output data file with the exchange duration histogram (dat)"                              , world_rank, cl_tags, &b_histo,     0);
    add_argument_mpi_i(argc,argv,"-bin"      ,&bin_width,                "Bin width for the exchange duration histogram (frames)"                                   , world_rank, cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-min"      , &min,                     "Minimum value of histogram (num lipids)"                                                  , world_rank, cl_tags, &b_min,       0);
    add_argument_mpi_i(argc,argv,"-max"      , &max,                     "Maximum value of histogram (num lipids)"                                                  , world_rank, cl_tags, &b_max,       0);
    add_argument_mpi_i(argc,argv,"-self"     , &b_self,                  "Count exchanges when the outgoing and replacing lipids are the same lipid? (0:no, 1:yes)" , world_rank, cl_tags, nullptr,      1);
    conclude_input_arguments_mpi(argc,argv,world_rank,program_name,cl_tags);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension_mpi(world_rank,"-d",binding_events_file_name,".be");
    check_extension_mpi(world_rank,"-l_in",replacing_lipid_file_name,".crd");
    check_extension_mpi(world_rank,"-l_out",leaving_lipid_file_name,".crd");

    if(b_histo == 1)
    {
        check_extension_mpi(world_rank,"-histo",histo_file_name,".dat");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check if a range was specified properly for the histogram                                                 //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(b_min != b_max) //only single boundary was specified
    {
        if(world_rank == 0)
        {
            printf("Only a single boundary was specified for the histogram. Please specify either both or no boundaries. \n");
        }
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read in lipid types                                                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Index replacing_lipid;
    Index leaving_lipid;
    replacing_lipid.get_index(replacing_lipid_file_name);
    leaving_lipid.get_index(leaving_lipid_file_name);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read in binding events                                                                                    //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Binding_events events;
    int result = events.get_binding_events(binding_events_file_name);

    if(result == 1) //bind events file exists
    { 
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // sort events by dwell time (largest first)                                                                 //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        events.organize_events(1);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // remove events with dwell time shorter than cutoff                                                         //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(cutoff > 0)
        {
            for(i=events.lipid_nr.size()-1; i>=0; i--) //loop over binding events
            {   
                if((double)events.dwell_t[i]*events.ef_dt < cutoff)
                {   
                    events.dwell_t.pop_back(); 
                    events.lipid_nr.pop_back();
                    events.bind_i.pop_back();  
                    events.bind_f.pop_back();  
                    events.res_nr.pop_back();  
                    events.res_name.pop_back(); 
                }
            }
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // characterize replacement binding                                                                          //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        int target_replacements = 0;       //total number of replacements by the correct lipid type
        int total               = 0;       //total number of lipid replacements
        iv1d delta_time(0,0);              //store the time between entering lipid and leaving lipid

        for(i=0; i<events.bind_i.size()-1; i++) //loop over binding events
        {
            //printf("here 1 res_name %s leaving_lipid %s \n",events.repeats_res_name[l].c_str(),leaving_lipid.c_str());
            for(j=0; j<leaving_lipid.index_s.size(); j++) //loop over leaving lipids
            {
                if(strcmp(events.res_name[i].c_str(), leaving_lipid.index_s[j].c_str()) == 0) //leaving lipid is of correct type
                {
                    int closest   = 999999;
                    int closest_i = -1;
                    int delta_t   = 0;

                    for(k=0; k<events.bind_i.size()-1; k++) //loop over binding events
                    {
                        int delta = abs(events.bind_i[k] - events.bind_f[i]);

                        if(delta < closest)
                        {
                            closest   = delta;
                            closest_i = k;
                            delta_t   = events.bind_i[k] - events.bind_f[i];
                        }
                    }

                    if(closest_i != i || b_self == 1) //dont count if the leaving and incoming lipids are one
                    {
                        //check if exchange was with a target lipid
                        for(k=0; k<replacing_lipid.index_s.size(); k++) //loop over replacing lipids
                        {
                            if(strcmp(events.res_name[closest_i].c_str(), replacing_lipid.index_s[k].c_str()) == 0) //replacing lipid is correct type
                            {
                                target_replacements++;
                                break;
                            }
                        }
                        total++;

                        if(b_histo == 1)
                        {
                            delta_time.push_back(delta_t);
                        }
                    }
                }
            }
        }
        double percent = (double)target_replacements/(double)(total) - lipid_fraction;

        printf(" %30s %f \n","observed exchange frequency:",(double)target_replacements/(double)(total));
        printf(" %30s %f \n","relative exchange frequency:",percent);
        printf(" %30s %d \n", "number of exchanges:",total);

        if(b_histo == 1)
        {
            printf("Binning data \n");
            Histogram_i histo;
            if(b_min == 1 && b_max == 1)
            {
                histo.set_range(min,max);
            }
            histo.bin_data(delta_time,bin_width);
            histo.write_histo(histo_file_name,"exchange duration (frames)"); 
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


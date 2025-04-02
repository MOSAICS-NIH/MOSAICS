
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <mpi.h>
#include <sstream>
#include <fstream>

using namespace std;

#include "../headers/multi_dim_vec.h"
#include "../headers/file_reader.h"
#include "../headers/vector_mpi.h"
#include "../headers/binding_events_common_routines.h"
#include "../headers/common_routines_mpi.h"
#include "../headers/common_routines.h"
#include "../headers/file_naming.h"
#include "../headers/binding_events.h"
#include "../headers/command_line_args_mpi.h"
#include "../headers/array.h"
#include "../headers/file_naming_mpi.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// The main function performing analysis                                                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, const char * argv[]) 
{
    //Here we define some variables used throughout
    string be_file_name;          //Input file name for the binding events files
    int i            = 0;         //Standard variable used in loops
    int j            = 0;         //Standard variable used in loops
    int k            = 0;         //Standard variable used in loops
    int l            = 0;         //Standard variable used in loops
    int m            = 0;         //Standard variable used in loops
    int world_size   = 0;         //How many mpi ranks
    int world_rank   = 0;         //Rank of the mpi process
    int stride       = 1;         //Skip frames 
    double cutoff    = 0.0;       //Lipids must be in the shell for this long to be counted
    sv1d cl_tags;                 //Holds a list of command line tags for the program

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set up the mpi environment                                                                                //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    MPI_Init(NULL, NULL);;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);      //get the world size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);      //get the process rank

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set program name and print info                                                                           //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    string program_name = "Binding Events Number Lipids Visited";
    if(world_rank == 0)
    {
        print_credits(argc,argv,program_name);
    }

    string program_description = "Binding Events Number Lipids Visited is a program that analyzes a binding events file and quantifies how many lipids have visited the region of interest as a function of time.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments_mpi(argc,argv,world_rank,program_description);
    add_argument_mpi_s(argc,argv,"-be"      , be_file_name,               "Input file name for binding events file (be)"          , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-stride"  , &stride,                    "Skip frames when counting visiting lipids vs time"     , world_rank, cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-cutoff"  , &cutoff,                    "Must be in the shell for this long to be counted (ps)" , world_rank, cl_tags, nullptr,      1);
    conclude_input_arguments_mpi(argc,argv,world_rank,program_name,cl_tags);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension_mpi(world_rank,"-be",be_file_name,".be");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Perform the analysis                                                                                      //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(world_rank == 0)
    {
        Binding_events events;   //used to store binding events data

        int result = events.get_binding_events_bin(be_file_name);

        if(result == 1) //be file successfully read
        {
            events.organize_events(0); //organize by b_i

            printf(" %20s %10s %10s \n","#Time:","#Lipids","#Percent"); 

    	    for(i=0; i<events.ef_frames; i+=stride) //loop over frames
            {
                iv1d lipid_numbers(0,0);    //store lipid numbers as they are encountered

    	        double   dt = events.ef_dt; //effective time step between frames
    	        double time = dt*i;         //time for current frame

    	        for(k=0; k<events.dwell_t.size(); k++) //loop over binding events
                {
                    double dwell_time = (double)events.dwell_t[k]*dt;   //dwell time for the event 
                    int      lipid_nr = (double)events.lipid_nr[k];     //lipid number for the event
                    double     bind_i = dt*(double)events.bind_i[k];    //initial binding time for event

                    if(dwell_time >= cutoff && bind_i <= time)
                    {
                        int found = 0;

                        for(l=0; l<lipid_numbers.size(); l++) //loop over the lipid numbers
                        {
                            if(lipid_numbers[l] == lipid_nr)
                            {
                                found = 1;
                            }
                        }

                        if(found == 0)
                        {
                            lipid_numbers.push_back(lipid_nr);
                        }
                    } 
                }
    
    	        int lipids_found = lipid_numbers.size();
    	        int total_lipids = events.num_lipids;
    	        double percent_found = (double)lipids_found/(double)total_lipids;

                printf(" %20f %10d %10f \n",time,lipids_found,total_lipids,percent_found);
            }
        }
        else
        {
            printf("unable to open binding events file %s \n",be_file_name.c_str());
        } 
    }
    
    if(world_rank == 0)
    {
        std::cout << "\nFormatting completed successfully" << "\n\n";
    }

    MPI_Finalize();

    return 0;
}


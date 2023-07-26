
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
#include "../headers/file_naming.h"
#include "../headers/common_routines_mpi.h"
#include "../headers/common_routines.h"
#include "../headers/binding_events.h"
#include "../headers/grid_lt.h"
#include "../headers/fit.h"
#include "../headers/histo.h"
#include "../headers/command_line_args_mpi.h"
#include "../headers/array.h"
#include "../headers/index.h"
#include "../headers/file_naming_mpi.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Main function which executes all other functions                                                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, const char * argv[]) 
{
    //Here we define some variables used throughout
    string in_file_name;                //Name of input file 
    string base_file_name_i;            //Base file name of input binding events files
    string be_file_name;                //Name of the binding_events input file
    string out_file_name;               //Name of the output files with the grid point selection and the big mask
    string rho_file_name;               //Name of the rho file
    string leaving_lipid_file_name;     //Name of the index file with leaving lipid types
    int i                 = 0;          //General variable used in loops
    int j                 = 0;          //General variable used in loops
    int k                 = 0;          //General variable used in loops
    int l                 = 0;          //General variable used in loops
    int world_size        = 0;          //Size of the mpi world
    int world_rank        = 0;          //Rank in the mpi world
    int odf               = 0;          //Rho data file format
    int stride            = 0;          //Skip stride frames
    int begin             = 0;          //Start on this frame
    int end               = 0;          //End on this frame
    int b_end             = 0;          //Did the user provide an end frame?
    int b_rho             = 0;          //Did the user specify a rho input file name?
    int b_min             = 0;          //Did the user specify a minimum range for the histogram?
    int b_max             = 0;          //Did the user specify a maximum range for the histogram?
    double cell_size      = 1;          //Distance between grid points
    double dt             = 0;          //Time step used for converting frames to time. set equal to ef_dt
    double cutoff         = 0;          //Cutoff for excluding data
    double avg_rho        = 0;          //The average lipid density over the grid
    double ex             = 0;          //Exclude binding events shorter than this (ns)
    double bin_width      = 1.0;        //Bin width for histogram
    double min            = 0.0;        //minimum value of histogram
    double max            = 0.0;        //maximum value of histogram
    sv1d cl_tags;                       //Holds a list of command line tags for the program

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
    string program_name = "Lipid Exchange Distances";
    if(world_rank == 0)
    {
        print_credits(argc,argv,program_name);
    }

    string program_description = "Lipid Exchange Distances is an analysis tool that lets the user compute the distance between exchanging lipids. This is plotted as a probability histogram.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments_mpi(argc,argv,world_rank,program_description);
    add_argument_mpi_s(argc,argv,"-d"       , base_file_name_i,           "Base filename for input binding events files"                            , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-be"      , be_file_name,               "Binding events file for the target region (be)"                          , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o"       , out_file_name,              "Output data files with exchange distance histogram (dat)"                , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-stride"  , &stride,                    "Skip stride frames"                                                      , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-b"       , &begin,                     "Start at this frame"                                                     , world_rank, cl_tags, &b_end,       0);
    add_argument_mpi_i(argc,argv,"-e"       , &end,                       "End on this frame"                                                       , world_rank, cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-rho"     , rho_file_name,              "Input data file with sample count (dat)"                                 , world_rank, cl_tags, &b_rho,       0);
    add_argument_mpi_d(argc,argv,"-cutoff"  , &cutoff,                    "Cutoff for excluding grid data in noise filtered voronoi diagrams (chi)" , world_rank, cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-odf"     , &odf,                       "Data file format for rho (0:matrix 1:vector)"                            , world_rank, cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-crd"     , leaving_lipid_file_name,    "Selection card with leaving lipid types (crd)"                           , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-ex"      , &ex,                        "Exclude binding events with a dwell time smaller than this (ps)"         , world_rank, cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-bin"     , &bin_width,                 "Bin width (nm)"                                                          , world_rank, cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-min"     , &min,                       "Minimum value of histogram (nm)"                                         , world_rank, cl_tags, &b_min,       0);
    add_argument_mpi_d(argc,argv,"-max"     , &max,                       "Maximum value of histogram (nm)"                                         , world_rank, cl_tags, &b_max,       0);
    conclude_input_arguments_mpi(argc,argv,world_rank,program_name,cl_tags);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension_mpi(world_rank,"-be",be_file_name,".be");
    check_extension_mpi(world_rank,"-o",out_file_name,".dat");
    check_extension_mpi(world_rank,"-crd",leaving_lipid_file_name,".crd");
    if(b_rho == 1)
    {
        check_extension_mpi(world_rank,"-rho",rho_file_name,".dat");
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
    // Read in lipid types from index file                                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Index leaving_lipid;
    leaving_lipid.get_index(leaving_lipid_file_name);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Create Grid for excluding insignificant data                                                              //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Grid_lt rho;

    if(b_rho == 1)
    {
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Set the grid format                                                                                       //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        rho.set_format(odf);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Read in grid data                                                                                         //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        rho.get_grid(rho_file_name);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Compute average rho                                                                                       //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        for(i=0; i<rho.size_x(); i++) //loop over x
        {
            for(j=0; j<rho.size_y(); j++) //loop over y
            {
                avg_rho = avg_rho + rho.grid[i][j][2][0];
            }
        }
        avg_rho = avg_rho/(rho.size_x()*rho.size_y());

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Exclude insignificant data                                                                                //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        rho.exclude_grid_data(cutoff,avg_rho,rho.grid);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read binding events file 0_0 to get the header info                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Binding_events events_ref;
    in_file_name = base_file_name_i + "_" + to_string(0) + "_" + to_string(0) + ".be";
    int result    = events_ref.get_binding_events(in_file_name);

    if(result == 0)
    {
        if(world_rank == 0)
        {
            printf("unable to open binding events file %s \n",in_file_name.c_str());
        }
        MPI_Finalize();
        return 0;
    }
    else
    {
        events_ref.get_binding_timeline();
    }

    cell_size     = sqrt(events_ref.APS);
    dt            = events_ref.ef_dt; //keep units of ps
    int ef_frames = (int)events_ref.ef_frames/stride + 1;

    if(b_end == 0)
    {
        end = events_ref.ef_frames;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Distribute the workload across the cores                                                                  //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int my_num_g_x = count_workload(world_size,world_rank,events_ref.num_g_x);

    //create array to hold each mpi processes num_g_x; Used for communication
    int world_num_g_x_ary[world_size];
    MPI_Allgather(&my_num_g_x, 1,MPI_INT,world_num_g_x_ary, 1, MPI_INT, MPI_COMM_WORLD );

    //allocate memory for world_num_g_x and copy data from the array
    iv1d world_num_g_x(world_size,0);
    for(i=0; i<world_size; i++)
    {
        world_num_g_x[i] = world_num_g_x_ary[i];
    }

    //print stats for distributing the num_g_x and distribute the num_g_x to each core
    int my_xi = 0;
    int my_xf = 0;
    int world_xi[world_size];
    int world_xf[world_size];
    get_workload(&my_xi,&my_xf,world_rank,world_num_g_x,events_ref.num_g_x,world_xi,world_xf);
    print_workload_stats(world_rank,world_xi,world_xf,world_num_g_x,world_size,"num_g_x","init","fin");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Get lipid blobs                                                                                           //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    events_ref.get_blobs(base_file_name_i,my_xi,my_xf,my_num_g_x,stride,ef_frames,world_rank);

    MPI_Barrier(MPI_COMM_WORLD);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Create vector to hold distances                                                                           //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    dv1d distances(0,0.0);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read in binding events at interface                                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Binding_events events;
    result = events.get_binding_events(be_file_name);

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
        if(ex > 0)
        {
            for(i=events.lipid_nr.size()-1; i>=0; i--) //loop over binding events
            {
                if((double)events.dwell_t[i]*events.ef_dt < ex)
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
        //for(i=0; i<events.bind_i.size()-1; i++) //loop over binding events
        for(i=0; i<events.bind_i.size(); i++) //loop over binding events
	{
            for(j=0; j<leaving_lipid.index_s.size(); j++) //loop over leaving lipids
            {
                if(strcmp(events.res_name[i].c_str(), leaving_lipid.index_s[j].c_str()) == 0) //leaving lipid is of correct type
                {
                    int closest   = 999999;
                    int closest_i = -1;
                    int in_frame  = 0;
                    int out_frame = events.bind_f[i]/stride;
                    int in_nr     =  0;
                    int out_nr    = events.lipid_nr[i];

                    //for(k=0; k<events.bind_i.size()-1; k++) //loop over binding events
                    for(k=0; k<events.bind_i.size(); k++) //loop over binding events
		    {
                        int delta = abs(events.bind_i[k] - events.bind_f[i]);

                        if(delta < closest)
                        {
                            closest   = delta;
                            closest_i = k;
                            in_frame  = events.bind_i[k]/stride;
                            in_nr     = events.lipid_nr[k];
                        }
                    }

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //                                                                                                           //
                    // find center of leaving lipid                                                                              //
                    //                                                                                                           //
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    events_ref.get_blobs_frame(out_frame,world_size,world_rank);

                    dv1d center_1(2,0.0);
                    if(world_rank == 0)
                    {
                        if(b_rho == 1)
                        {
                            for(k=0; k<rho.size_x(); k++) //loop over x
                            {
                                for(l=0; l<rho.size_y(); l++) //loop over y
                                {
                                    if(rho.grid[k][l][2][1] == 1)
                                    {
                                        events_ref.blob_nan_frame[l][k] = 1;
                                    }
                                }
                            }
                        }

                        center_1 = events_ref.find_blobs_center(out_nr);

                        string tag = "_" + to_string(out_frame);
                        string this_file_name = add_tag(out_file_name,tag);
                        events_ref.write_blobs_frame(this_file_name);
                    }

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //                                                                                                           //
                    // find center of incoming lipid                                                                             //
                    //                                                                                                           //
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    events_ref.get_blobs_frame(in_frame,world_size,world_rank);

                    dv1d center_2(2,0.0);
                    if(world_rank == 0)
                    {
                        if(b_rho == 1)
                        {
                            for(k=0; k<rho.size_x(); k++) //loop over x
                            {
                                for(l=0; l<rho.size_y(); l++) //loop over y
                                {
                                    if(rho.grid[k][l][2][1] == 1)
                                    {
                                        events_ref.blob_nan_frame[l][k] = 1;
                                    }
                                }
                            }
                        }

                        center_2 = events_ref.find_blobs_center(in_nr);

                        string tag = "_" + to_string(in_frame);
                        string this_file_name = add_tag(out_file_name,tag);
                        events_ref.write_blobs_frame(this_file_name);
                    }

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //                                                                                                           //
                    // find distance between centers                                                                             //
                    //                                                                                                           //
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    if(world_rank == 0)
                    {
                        if(closest_i != i) //dont count exchanges with self
                        {
                            dv1d dif(2,0.0);
                            dif[0] = center_1[0] - center_2[0];
                            dif[1] = center_1[1] - center_2[1];

                            double dist = cell_size*sqrt(dif[0]*dif[0] + dif[1]*dif[1]);

                            distances.push_back(dist);

                            //printf("in_frame %10d out_frame %10d in_nr %10d out_nr %10d dist %10f \n",in_frame,out_frame,in_nr,out_nr,dist);
                            printf("in_frame %10d out_frame %10d in_nr %10d out_nr %10d in_x %10f in_y %10f out_x %10f out_y %10f dist %10f\n",in_frame,out_frame,in_nr,out_nr,center_2[0],center_2[1],center_1[0],center_1[1],dist);
                        }
                    }
                }
            }
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Bin data and write probability distribution to file                                                       //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(world_rank == 0)
        {
            Histogram_d histo;
            if(b_min == 1 && b_max == 1)
            {
	        histo.set_range(min,max);
            }
            histo.bin_data(distances,bin_width);
            histo.write_histo(out_file_name,"distance (nm)");
        }
    }
    else //binding events file does not exits
    {
        printf("Could not find binding events file \n");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Print final output and end program                                                                        //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(world_rank == 0)
    {
        std::cout << "\nFormatting completed successfully" << "\n\n";
    }

    MPI_Finalize();

    return 0;
}



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
#include "../headers/index.h"
#include "../headers/grid_lt.h"
#include "../headers/command_line_args_mpi.h"
#include "../headers/array.h"
#include "../headers/fit.h"
#include "../headers/file_naming_mpi.h"
#include "../headers/performance.h"        

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// The main function performing analysis                                                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, const char * argv[]) 
{
    //Here we define some variables used throughout
    string in_file_name;          //Name of the input file
    string out_file_name;         //Name of the output file
    string lip_t_file_name;       //Name of the lipid types file
    string base_file_name_o;      //Base file name of output files
    string rho_file_name;         //Name of the rho file
    int i                 = 0;    //Standard variable used in loops
    int j                 = 0;    //Standard variable used in loops
    int k                 = 0;    //Standard variable used in loops
    int l                 = 0;    //Standard variable used in loops
    int m                 = 0;    //Standard variable used in loops
    int world_size        = 0;    //Size of the mpi world
    int world_rank        = 0;    //Rank in the mpi world
    int b_rho             = 0;    //Did the user specify a rho input file name?
    int odf               = 0;    //Rho data file format
    int threshold         = 0;    //Cutoff for mending fragmented binding events
    int counter           = 0;    //How many times the "program run time" been displayed
    int grid_counter      = 0;    //Count lattice points as they are encountered
    int b_ref             = 0;    //Did the user provide a ref be file? 
    double slope          = 0.0;  //slope of LnP vs time
    double yint           = 0.0;  //ying of LnP vs time
    double r2             = 0.0;  //Correlation coeficient in linear regression
    double koff           = 0.0;  //The koff value
    double cutoff_dt      = 0.0;  //Exclude data with cutoff less than this
    double cutoff         = 0.0;  //Cutoff for excluding data
    double avg_rho        = 0.0;  //The average lipid density over the grid
    clock_t t;                    //Keeps the time for testing performance
    clock_t t1;                   //Keeps the time for testing performance
    sv1d cl_tags;                 //Holds a list of command line tags for the program

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set up the mpi environment                                                                                //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    MPI_Init(NULL, NULL);;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);      //get the world size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);      //get the process rank

    //create object for logging performance data
    Performance perf;

    //take the initial time
    t = clock();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set program name/description and print info                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    string program_name = "Binding Events Analyzer";

    if(world_rank == 0)
    {
        print_credits(argc,argv,program_name);
    }

    string program_description = "Binding Events Analyzer is an analysis tool used for reading binding events files and computing the mean dwell time, and koff, etc. This data is projected onto the XY plane.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments_mpi(argc,argv,world_rank,program_description);
    add_argument_mpi_s(argc,argv,"-d"       , in_file_name,               "Filename for input binding events file (be)"                       , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o"       , base_file_name_o,           "Filename for driving output data filenames (dat)"                  , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-crd"     , lip_t_file_name,            "Selection card with lipid types (crd)"                             , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-ex"      , &cutoff_dt,                 "Exclude data with dwell time less than this (ps)"                  , world_rank, cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-rho"     , rho_file_name,              "Input data file with sample count (dat)"                           , world_rank, cl_tags, &b_rho,       0);
    add_argument_mpi_d(argc,argv,"-cutoff"  , &cutoff,                    "Cutoff for excluding grid data (chi)"                              , world_rank, cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-odf"     , &odf,                       "Data file format for sample count data (0:matrix 1:vector)"        , world_rank, cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-repair"  , &threshold,                 "Maximum allowed size (frames) for mending fragmented events"       , world_rank, cl_tags, nullptr,      0);
    conclude_input_arguments_mpi(argc,argv,world_rank,program_name,cl_tags);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension_mpi(world_rank,"-crd",lip_t_file_name,".crd");
    check_extension_mpi(world_rank,"-d",in_file_name,".be");
    check_extension_mpi(world_rank,"-o",base_file_name_o,".dat");
    if(b_rho == 1)
    {
        check_extension_mpi(world_rank,"-rho",rho_file_name,".dat");
    }

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
    // Read binding events file to get the header info                                                           //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Binding_events events_ref;
    int result = events_ref.get_info(in_file_name); 

    if(result == 0)
    {
        if(world_rank == 0)
        {
            printf("unable to open binding events file %s \n",in_file_name.c_str());
        }
        MPI_Finalize();
        return 0;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Distribute the workload across the cores                                                                  //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int my_num_g = count_workload(world_size,world_rank,events_ref.num_g_x*events_ref.num_g_y);

    //create array to hold each mpi processes my_num_g; Used for communication
    int world_num_g_ary[world_size];
    MPI_Allgather(&my_num_g, 1,MPI_INT,world_num_g_ary, 1, MPI_INT, MPI_COMM_WORLD );

    //allocate memory for world_num_g and copy data from the array
    iv1d world_num_g(world_size,0);
    for(i=0; i<world_size; i++)
    {
        world_num_g[i] = world_num_g_ary[i];
    }

    //print stats for distributing the grid and distribute the grid to each core
    int my_gi = 0;
    int my_gf = 0;
    iv1d world_gi(world_size);
    iv1d world_gf(world_size);
    get_grid_points_alt(&my_gi,&my_gf,world_rank,world_size,world_num_g,events_ref.num_g_x,events_ref.num_g_y,world_gi,world_gf);
    print_workload_stats_alt(world_rank,world_gi,world_gf,world_num_g,world_size);
    MPI_Barrier(MPI_COMM_WORLD);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read in lipid types                                                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Index lip_t;
    lip_t.get_index(lip_t_file_name);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Allocate memory for the grid                                                                              //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    dv1d num_events_local(my_num_g,0.0);
    dv1d dwell_time_local(my_num_g,0.0);
    dv1d stdev_local     (my_num_g,0.0);
    dv1d r2_local        (my_num_g,0.0);
    dv1d koff_local      (my_num_g,0.0);
    dv1d largest_local   (my_num_g,0.0);

    //create some items to track time in each region
    clock_t t_read_be  = 0.0;   //time spent reading binding events
    clock_t t_suppress = 0.0;   //time spent mending fragmented events
    clock_t t_sweep    = 0.0;   //time spent removing events with too small a dwell time
    clock_t t_dwell_t  = 0.0;   //time spent computing dwell times
    clock_t t_variance = 0.0;   //time spent computing variance
    clock_t t_largest  = 0.0;   //time spent getting largest dwell time
    clock_t t_koff     = 0.0;   //time spent getting koff
    clock_t t_waiting  = 0.0;   //time spent waiting on other cores to finish

    //log time spent leading up to main analysis
    perf.log_time((clock() - t)/CLOCKS_PER_SEC,"Other");

    //get time for reporting progress
    t1 = clock();

    Binding_events events;
    result       = events.get_info(in_file_name);

    if(result == 1) 
    {
        if(world_rank == 0)
        {
            printf("Reading binding events data and computing dwell times. \n");
            printf("-----------------------------------------------------------------------------------------------------------------------------------\n");
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // loop over the grid                                                                                        //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        for(i=0; i<=events_ref.num_g_x; i++) //loop over x
        {
            for(j=0; j<events_ref.num_g_y; j++) //loop over y
            {
                if(grid_counter >= my_gi && grid_counter <= my_gf)
                {
                    int ef_count = grid_counter - my_gi;
 
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //                                                                                                           //
                    // Read in binding events                                                                                    //
                    //                                                                                                           //
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    t         = clock();
                    result    = events.get_binding_events_xy(in_file_name,i,j);
                    t_read_be = t_read_be + clock() - t;

                    if(events.lipid_nr.size() > 0)
                    {
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        //                                                                                                           //
                        // Mend any fragmented events                                                                                //
                        //                                                                                                           //
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        t = clock();
                        events.suppress_timeline_noise_be_alt(threshold);  //mend fragmented events
                        t_suppress = t_suppress + clock() - t;

                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        //                                                                                                           //
                        // remove events with dwell time shorter than cutoff_dt                                                      //
                        //                                                                                                           //
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        t = clock();
                        events.screen_events(cutoff_dt);
                        t_sweep = t_sweep + clock() - t;

                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        //                                                                                                           //
                        // Compute average dwell time                                                                                //
                        //                                                                                                           //
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        double dwell_time = 0.0;
                        int ef_num_events = 0;

                        t = clock();
                        for(k=0; k<events.dwell_t.size(); k++) //loop over binding events
                        {
                            for(l=0; l<lip_t.index_s.size(); l++) //loop over lipid types
                            {
                                if(strcmp(events.res_name[k].c_str(), lip_t.index_s[l].c_str()) == 0) //lipid type is correct
                                {
                                    dwell_time = dwell_time + (double)events.dwell_t[k];
                                    ef_num_events = ef_num_events + 1;
                                }
                            }
                        }
                        dwell_time                 = dwell_time/(double)ef_num_events;
            	        dwell_time_local[ef_count] = dwell_time*events.ef_dt;
                        num_events_local[ef_count] = (double)ef_num_events;
                        t_dwell_t = t_dwell_t + clock() - t;

                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        //                                                                                                           //
                        // Compute variance in dwell time                                                                            //
                        //                                                                                                           //
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        double stdev         = 0.0;
                        double sum_o_squares = 0.0;
 
                        t = clock();
                        for(k=0; k<events.dwell_t.size(); k++) //loop over binding events
                        {
                            for(l=0; l<lip_t.index_s.size(); l++) //loop over lipid types
                            {
                                if(strcmp(events.res_name[k].c_str(), lip_t.index_s[l].c_str()) == 0) //lipid type is correct
                                {
                                    sum_o_squares = sum_o_squares + pow((double)events.dwell_t[k]-dwell_time,2);
                                }
                            }
                        }
                        stdev = sqrt(sum_o_squares/(double)(ef_num_events-1));
                        stdev_local[ef_count] = stdev*events.ef_dt;
                        t_variance = t_variance + clock() - t;

                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        //                                                                                                           //
                        // Compute effective largest dwell_t                                                                         //
                        //                                                                                                           //
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        int largest_dwell_t = 0;
      
                        t = clock();
                        for(k=0; k<events.dwell_t.size(); k++) //loop over binding events
                        {
                            for(l=0; l<lip_t.index_s.size(); l++) //loop over lipid types
                            {
                                if(strcmp(events.res_name[k].c_str(), lip_t.index_s[l].c_str()) == 0) //lipid type is correct
                                {
                                    if(events.dwell_t[k] > largest_dwell_t)
                                    {
                                        largest_dwell_t = events.dwell_t[k];
                                    }
                                }
                            }
                        }
                        largest_local[ef_count] = (double)largest_dwell_t*events.ef_dt;
                        t_largest = t_largest + clock() - t;

                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        //                                                                                                           //
                        // Compute frequencies                                                                                       //
                        //                                                                                                           //
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        double dwell_time_freq[largest_dwell_t-1];       //holds how many binding events had dwell time i
                        init_darray(dwell_time_freq,largest_dwell_t-1);

                        t = clock();
                        for(k=0; k<largest_dwell_t; k++) //loop over dwell times
                        {
                            for(l=0; l<events.bind_i.size(); l++) //loop over binding events
                            {
                                if(events.dwell_t[l] == k+1) //dwell time matches
                                {
                                    for(m=0; m<lip_t.index_s.size(); m++) //loop over lipid types
                                    {
                                        if(strcmp(events.res_name[l].c_str(), lip_t.index_s[m].c_str()) == 0) //lipid type is correct
                                        {
                                            dwell_time_freq[k] = dwell_time_freq[k] + 1;
                                        }
                                    }
                                }
                            }
                            dwell_time_freq[k] = dwell_time_freq[k]/(double)ef_num_events;
                            //printf("dwell_time_freq %10.3f \n",dwell_time_freq[k]);
                        }

                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        //                                                                                                           //
                        // Compute probability of lasting t or longer. (add probabilities for t and longer)                          //
                        //                                                                                                           //
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        double nan = -999999;       //assign a negative value to freq for frequenies with no data
                        int count = 0;              //how many dwell times have a positive frequency
                        for(k=0; k<largest_dwell_t; k++) //loop over dwell times
                        {
                            if(dwell_time_freq[k] > 0) //inclue data in fit
                            {
                                for(l=k+1; l<largest_dwell_t; l++) //add probabilities
                                {
                                    dwell_time_freq[k] = dwell_time_freq[k] + dwell_time_freq[l];
                                }
                                count++;
                            }
                            else //exclue data from fit
                            {
                                dwell_time_freq[k] = nan;
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
                        for(k=0; k<largest_dwell_t; k++) //loop over dwell times
                        {
                            if(dwell_time_freq[k] != nan) //dont include dwell times with no samples
                            {
                                log_freq[count] = log(dwell_time_freq[k]);
                                time[count]     = (k+1)*events.ef_dt;
                                count++;
                            }
                        }

                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        //                                                                                                           //
                        // Compute koff                                                                                              //
                        //                                                                                                           //
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        least_squares_regression(count,time,log_freq,&slope,&yint,&r2);
                        koff                = -slope;
                        r2_local[ef_count]   = r2;
                        koff_local[ef_count] = koff;
                        t_koff = t_koff + clock() - t;
                    }
                    else //no binding events present 
                    {
                    }

                    //report progress and estimated time to completion
                    int current_step = ef_count + 1; 
                    int my_steps     = my_gf - my_gi + 1;
                    ot_time_stats(t1,&counter,current_step,my_steps,world_rank,"lattice point");
                }
                grid_counter++;
            }
        }
    }
    else 
    {
        printf("failure opening %s. Make sure the file exists. \n",in_file_name.c_str());
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }

    t = clock();
    MPI_Barrier(MPI_COMM_WORLD);
    t_waiting = t_waiting + clock() - t;

    //log time spent performing main analysis
    perf.log_time(t_read_be/CLOCKS_PER_SEC,"read be");
    perf.log_time(t_suppress/CLOCKS_PER_SEC,"mend be");
    perf.log_time(t_sweep/CLOCKS_PER_SEC,"exclude be");
    perf.log_time(t_dwell_t/CLOCKS_PER_SEC,"dwell time");
    perf.log_time(t_variance/CLOCKS_PER_SEC,"dwell time var.");
    perf.log_time(t_largest/CLOCKS_PER_SEC,"largest time");
    perf.log_time(t_koff/CLOCKS_PER_SEC,"koff");
    perf.log_time(t_waiting/CLOCKS_PER_SEC,"waiting");

    //reset the clock
    t = clock();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Collect grids                                                                                             //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    dv2d num_events_global(events_ref.num_g_y,dv1d(events_ref.num_g_x,0.0));
    dv2d dwell_time_global(events_ref.num_g_y,dv1d(events_ref.num_g_x,0.0));
    dv2d stdev_global     (events_ref.num_g_y,dv1d(events_ref.num_g_x,0.0));
    dv2d r2_global        (events_ref.num_g_y,dv1d(events_ref.num_g_x,0.0));
    dv2d koff_global      (events_ref.num_g_y,dv1d(events_ref.num_g_x,0.0));
    dv2d largest_global   (events_ref.num_g_y,dv1d(events_ref.num_g_x,0.0));
    iv2d nan              (events_ref.num_g_y,iv1d(events_ref.num_g_x,0));

    //exclude insignificant data
    if(b_rho == 1)
    {
        for(i=0; i<rho.size_x(); i++) //loop over x
        {
            for(j=0; j<rho.size_y(); j++) //loop over y
            {
                nan[j][i] = rho.grid[i][j][2][1];
            }
        }
    }

    gather_grid_d_gp_alt(world_size,world_rank,world_gi,world_gf,events_ref.num_g_x,events_ref.num_g_y,world_num_g,num_events_local,num_events_global);
    gather_grid_d_gp_alt(world_size,world_rank,world_gi,world_gf,events_ref.num_g_x,events_ref.num_g_y,world_num_g,dwell_time_local,dwell_time_global);
    gather_grid_d_gp_alt(world_size,world_rank,world_gi,world_gf,events_ref.num_g_x,events_ref.num_g_y,world_num_g,stdev_local,stdev_global);
    gather_grid_d_gp_alt(world_size,world_rank,world_gi,world_gf,events_ref.num_g_x,events_ref.num_g_y,world_num_g,r2_local,r2_global);
    gather_grid_d_gp_alt(world_size,world_rank,world_gi,world_gf,events_ref.num_g_x,events_ref.num_g_y,world_num_g,koff_local,koff_global);
    gather_grid_d_gp_alt(world_size,world_rank,world_gi,world_gf,events_ref.num_g_x,events_ref.num_g_y,world_num_g,largest_local,largest_global);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Write grid data to output files                                                                           //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(world_rank == 0)
    {
	out_file_name = add_tag(base_file_name_o,"_num_events");
        write_grid_to_file(events_ref.num_g_x,events_ref.num_g_y,nan,out_file_name,num_events_global);

        out_file_name = add_tag(base_file_name_o,"_dwell_time");
	write_grid_to_file(events_ref.num_g_x,events_ref.num_g_y,nan,out_file_name,dwell_time_global);

        out_file_name = add_tag(base_file_name_o,"_stdev");
        write_grid_to_file(events_ref.num_g_x,events_ref.num_g_y,nan,out_file_name,stdev_global);

        out_file_name = add_tag(base_file_name_o,"_r2");
        write_grid_to_file(events_ref.num_g_x,events_ref.num_g_y,nan,out_file_name,r2_global);

        out_file_name = add_tag(base_file_name_o,"_koff");
        write_grid_to_file(events_ref.num_g_x,events_ref.num_g_y,nan,out_file_name,koff_global);

        out_file_name = add_tag(base_file_name_o,"_largest_dwell_time");
        write_grid_to_file(events_ref.num_g_x,events_ref.num_g_y,nan,out_file_name,largest_global);
    }

    //log time spent performing main analysis
    perf.log_time((clock() - t)/CLOCKS_PER_SEC,"Collect Data");

    //print the performance stats
    perf.print_stats();

    if(world_rank == 0)
    {
        std::cout << "\nFormatting completed successfully" << "\n\n";
    }

    MPI_Finalize();

    return 0;
}


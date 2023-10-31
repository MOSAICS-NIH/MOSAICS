
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
// This function determines how much time has passed and gives an estimate of the time remaining.            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ot_time_stats(clock_t t,int *counter,int current_step,int my_steps,int world_rank,string my_string)
{
    if(world_rank == 0)
    {
        double seconds = 0;

        seconds = (clock() - t)/CLOCKS_PER_SEC;

        if((int)seconds/10 > *counter)
        {
            double percent_done         = ((double)(current_step)/(double)my_steps)*100.0;
            double estimated_total_time = 100.0*seconds/percent_done;
            double time_remaining       = estimated_total_time - seconds;

            if(current_step == my_steps)
            {
                estimated_total_time = 0.0;
            }

            int phr  = 0;
            int pmin = 0;
            int psec = 0;
            int lhr  = 0;
            int lmin = 0;
            int lsec = 0;

            phr  = (seconds)/(60*60);
            pmin = (seconds - (phr*60*60))/60;
            psec = seconds - (phr*60*60) - (pmin*60);

            lhr  = ((int)time_remaining)/(60*60);
            lmin = (time_remaining - (lhr*60*60))/60;
            lsec = time_remaining - (lhr*60*60) - (lmin*60);

            if(percent_done != 0)
            {
                printf("Finished %s %7d with %5.1f percent done overall in %2d hr %2d min %2d sec. Estimated time to completion is %2d hr %2d min %2d sec. \n",my_string.c_str(),current_step+1,percent_done,phr,pmin,psec,lhr,lmin,lsec);
            }
            *counter = *counter + 1;
        }
    }
}

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
    string base_file_name;        //Name of the input file (base)
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
    add_argument_mpi_s(argc,argv,"-d"       , base_file_name,             "Base filename for input binding events files "                     , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o"       , base_file_name_o,           "Base filename for output data files"                               , world_rank, cl_tags, nullptr,      1);
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
    // Read binding events file 0_0 to get the header info                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Binding_events events_ref;
    in_file_name = base_file_name + "_" + to_string(0) + "_" + to_string(0) + ".be";
    int result   = events_ref.get_binding_events(in_file_name);

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
    int my_num_g_x = count_workload(world_size,world_rank,events_ref.num_g_x);

    //create array to hold each mpi processes my_num_g_x; Used for communication
    int world_num_g_x_ary[world_size];
    MPI_Allgather(&my_num_g_x, 1,MPI_INT,world_num_g_x_ary, 1, MPI_INT, MPI_COMM_WORLD );

    //allocate memory for world_num_g_x and copy data from the array
    iv1d world_num_g_x(world_size,0);
    for(i=0; i<world_size; i++)
    {
        world_num_g_x[i] = world_num_g_x_ary[i];
    }

    //print stats for distributing the grid and distribute the grid to each core
    int my_xi = 0;
    int my_xf = 0;
    int world_xi[world_size];
    int world_xf[world_size];
    get_workload(&my_xi,&my_xf,world_rank,world_num_g_x,events_ref.num_g_x,world_xi,world_xf);
    print_workload_stats(world_rank,world_xi,world_xf,world_num_g_x,world_size,"num_g_x","xi","xf");

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
    dv2d num_events_local(events_ref.num_g_y,dv1d(my_num_g_x,0.0));
    dv2d dwell_time_local(events_ref.num_g_y,dv1d(my_num_g_x,0.0));
    dv2d stdev_local     (events_ref.num_g_y,dv1d(my_num_g_x,0.0));
    dv2d r2_local        (events_ref.num_g_y,dv1d(my_num_g_x,0.0));
    dv2d koff_local      (events_ref.num_g_y,dv1d(my_num_g_x,0.0));
    dv2d largest_local   (events_ref.num_g_y,dv1d(my_num_g_x,0.0));

    //log time spent leading up to main analysis
    perf.log_time((clock() - t)/CLOCKS_PER_SEC,"Other");

    //reset the clock
    t = clock();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // loop over the grid                                                                                        //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(i=my_xi; i<=my_xf; i++) //loop over x
    {
        int ef_x = i - my_xi;

        if(world_rank == 0)
        {   
            printf("Working on grid collumn %d (%d) \n",ef_x,i);
        }

        for(j=0; j<events_ref.num_g_y; j++) //loop over y
        {
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                                                           //
            // Read in binding events                                                                                    //
            //                                                                                                           //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            Binding_events events;
            in_file_name = base_file_name + "_" + to_string(i) + "_" + to_string(j) + ".be";
            result       = events.get_binding_events(in_file_name);

            if(result == 1)
            {
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // sort events by dwell time (largest first)                                                                 //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                events.organize_events(1);

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // remove events with dwell time shorter than cutoff_dt                                                      //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                if(cutoff_dt > 0.0)
                {
                    for(k=events.lipid_nr.size()-1; k>=0; k--) //loop over binding events
                    {
                        if((double)events.dwell_t[k]*events.ef_dt < cutoff_dt)
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
                // Make a timeline and mend any fragmented events                                                            //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                events.get_binding_timeline();              //make a timeline        
                events.suppress_timeline_noise(threshold);  //mend fragmented events
                events.binding_events_from_timeline();      //generate binding events from mended timeline

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // Compute average dwell time                                                                                //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                double dwell_time = 0.0;
                int ef_num_events = 0;

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
                dwell_time                = dwell_time/(double)ef_num_events;
                dwell_time_local[j][ef_x] = dwell_time*events.ef_dt;
                num_events_local[j][ef_x] = (double)ef_num_events;

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // Compute variance in dwell time                                                                            //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                double stdev         = 0.0;
                double sum_o_squares = 0.0;

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
                stdev_local[j][ef_x] = stdev*events.ef_dt;

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // Compute effective largest dwell_t                                                                         //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                int largest_dwell_t = 0;
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
                largest_local[j][ef_x] = (double)largest_dwell_t*events.ef_dt;

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // Compute frequencies                                                                                       //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                double dwell_time_freq[largest_dwell_t-1];       //holds how many binding events had dwell time i
                init_darray(dwell_time_freq,largest_dwell_t-1);

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
                r2_local[j][ef_x]   = r2;
                koff_local[j][ef_x] = koff;
            }
            else //no such binding events file 
            {
            }

            //report progress and estimated time to completion
            int current_step = events_ref.num_g_y*ef_x + j + 1; 
            int my_steps     = (my_xf - my_xi + 1)*events_ref.num_g_y;
	    ot_time_stats(t,&counter,current_step,my_steps,world_rank,"lattice point");
        }
    }

    //log time spent performing main analysis
    perf.log_time((clock() - t)/CLOCKS_PER_SEC,"Main Loop");

    //reset the clock
    t = clock();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Collect grids                                                                                             //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    dv2d num_events_global(events_ref.num_g_y,dv1d(events_ref.num_g_x,0.0));
    dv2d dwell_time_global(events_ref.num_g_y,dv1d(events_ref.num_g_x,0.0));
    dv2d stdev_global(events_ref.num_g_y,dv1d(events_ref.num_g_x,0.0));
    dv2d r2_global(events_ref.num_g_y,dv1d(events_ref.num_g_x,0.0));
    dv2d koff_global(events_ref.num_g_y,dv1d(events_ref.num_g_x,0.0));
    dv2d largest_global(events_ref.num_g_y,dv1d(events_ref.num_g_x,0.0));
    iv2d nan(events_ref.num_g_y,iv1d(events_ref.num_g_x,0));

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

    gather_grid_d_gp(world_size,world_rank,my_num_g_x,events_ref.num_g_x,events_ref.num_g_y,world_num_g_x,num_events_local,num_events_global);
    gather_grid_d_gp(world_size,world_rank,my_num_g_x,events_ref.num_g_x,events_ref.num_g_y,world_num_g_x,dwell_time_local,dwell_time_global);
    gather_grid_d_gp(world_size,world_rank,my_num_g_x,events_ref.num_g_x,events_ref.num_g_y,world_num_g_x,stdev_local,stdev_global);
    gather_grid_d_gp(world_size,world_rank,my_num_g_x,events_ref.num_g_x,events_ref.num_g_y,world_num_g_x,r2_local,r2_global);
    gather_grid_d_gp(world_size,world_rank,my_num_g_x,events_ref.num_g_x,events_ref.num_g_y,world_num_g_x,koff_local,koff_global);
    gather_grid_d_gp(world_size,world_rank,my_num_g_x,events_ref.num_g_x,events_ref.num_g_y,world_num_g_x,largest_local,largest_global);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Write grid data to output files                                                                           //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(world_rank == 0)
    {
        out_file_name = base_file_name_o + "_num_events.dat";
        write_grid_to_file(events_ref.num_g_x,events_ref.num_g_y,nan,out_file_name,num_events_global);

        out_file_name = base_file_name_o + "_dwell_time.dat";
        write_grid_to_file(events_ref.num_g_x,events_ref.num_g_y,nan,out_file_name,dwell_time_global);

        out_file_name = base_file_name_o + "_stdev.dat";
        write_grid_to_file(events_ref.num_g_x,events_ref.num_g_y,nan,out_file_name,stdev_global);

        out_file_name = base_file_name_o + "_r2.dat";
        write_grid_to_file(events_ref.num_g_x,events_ref.num_g_y,nan,out_file_name,r2_global);

        out_file_name = base_file_name_o + "_koff.dat";
        write_grid_to_file(events_ref.num_g_x,events_ref.num_g_y,nan,out_file_name,koff_global);

        out_file_name = base_file_name_o + "_largest_dwell_time.dat";
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


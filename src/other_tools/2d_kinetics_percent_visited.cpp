
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
#include "../headers/binding_events_common_routines.h"
#include "../headers/common_routines_mpi.h"
#include "../headers/common_routines.h"
#include "../headers/binding_events.h"
#include "../headers/index.h"                             //This has a class for working with index files
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
    string in_file_name;          //Name of the binding events file being read
    string base_file_name;        //Base file name for the binding events files
    string lip_t_file_name;       //Lipid types for analysis
    string out_file_name;         //Name of the output file
    int i          = 0;           //Standard variable used in loops
    int j          = 0;           //Standard variable used in loops
    int k          = 0;           //Standard variable used in loops
    int l          = 0;           //Standard variable used in loops
    int m          = 0;           //Standard variable used in loops
    int freq       = 1;           //How often to report percent visited
    int world_size = 0;           //How many mpi ranks
    int world_rank = 0;           //Rank of the mpi process
    int num_lipids = 0;           //Number of target lipids
    int b_num_lip  = 0;           //Was the number of lipids specified?
    double dt      = 0;           //Time step used for output     

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
    string program_name = "2d Kinetics percent Visited";
    if(world_rank == 0)
    {
        print_credits(argc,argv,program_name);
    }

    string program_description = "2d Kinetics Percent Visited is an analysis tool that lets the user compute the percentage of lipids to visit each grid point using binding events data.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments_mpi(argc,argv,world_rank,program_description);
    add_argument_mpi_s(argc,argv,"-d"       , base_file_name,             "Base filename for input binding events files"                        , world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o"       , out_file_name,              "Output filename used to derive names for percent visited data (dat)" , world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-freq"    , &freq,                      "How often to report the percent visited"                             , world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-crd"     , lip_t_file_name,            "Selection card with lipid types (crd)"                               , world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-lipids"  , &num_lipids,                "How many lipids are there of the target type in the target leaflet?" , world_rank, &b_num_lip,   0);
    conclude_input_arguments_mpi(argc,argv,world_rank,program_name);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension_mpi(world_rank,"-o",out_file_name,".dat");
    check_extension_mpi(world_rank,"-crd",lip_t_file_name,".crd");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Create index to hold lipid types                                                                          //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //create index for lipid types
    Index lip_t;

    //read the index files
    lip_t.get_index(lip_t_file_name);

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

    dt = events_ref.ef_dt/1000;

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
    // Allocate memory for the grid                                                                              //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int size = 0;
    for(i=0; i<events_ref.ef_frames; i++) //loop over time line frames
    {
        if(i%freq == 0)
        {
            size++;
        }
    }
    dv3d percent_local(size,dv2d(events_ref.num_g_y,dv1d(my_num_g_x,0.0)));

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

            if(result == 1) //binding events file exists
            {
                if(events.lipid_nr.size() > 0)
                {
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //                                                                                                           //
                    // Make bound_time_line_ij                                                                                   //
                    //                                                                                                           //
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    events.get_binding_timeline();

                    //check if the number of lipids was specified
                    if(b_num_lip == 0)
                    {
                        num_lipids = events.num_lipids;
                    }

                    iv1d visited(events.num_lipids, 0);  //tells if each lipid has visited previously
                    int ef_frame = -1;                   //used for adding data to percent_local

                    for(k=0; k<events.ef_frames; k++) //loop over time line frames
                    {
                        for(l=0; l<events.num_lipids; l++) //loop over time line lipids
                        {
                            for(m=0; m<lip_t.index_s.size(); m++) //loop over lipid types
                            {
                                if(strcmp(events.time_line_res_name[l].c_str(), lip_t.index_s[m].c_str()) == 0) //lipid type is correct  
                                {
                                    if(events.bound_time_line[k][l] == 1)
                                    {
                                        visited[l] = 1;
                                    }
                                }
                            }
                        }

                        if(k%freq == 0) //compute percent visited
                        {
                            ef_frame++;

                            int count = 0;

                            for(l=0; l<events.num_lipids; l++) //loop over lipids
                            {
                                if(visited[l] == 1)
                                {
                                    count++;
                                }                            
                            }
                            double percent = (double)count/(double)num_lipids;

                            //store percent visited
                            percent_local[ef_frame][j][ef_x] = percent;                    
                        }
                    }
                }
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Collect grids and write data to output files                                                              //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int ef_frame = -1;
    for(i=0; i<events_ref.ef_frames; i++) //loop over time line frames
    {
        if(i%freq == 0)
        {
            if(world_rank == 0)
            {
                printf("Collect grid data. Working on frame %d \n",i);
            }

            ef_frame++;

            dv2d percent_global(events_ref.num_g_y,dv1d(events_ref.num_g_x,0.0));
            iv2d nan(events_ref.num_g_y,iv1d(events_ref.num_g_x,0));

            gather_grid_d_gp(world_size,world_rank,my_num_g_x,events_ref.num_g_x,events_ref.num_g_y,world_num_g_x,percent_local[ef_frame],percent_global);

            if(world_rank == 0)
            {
                string tag            = "_" + to_string(i);
                string this_file_name = add_tag(out_file_name,tag);

                write_grid_to_file(events_ref.num_g_x,events_ref.num_g_y,nan,this_file_name,percent_global);
            }
        }
    }

    if(world_rank == 0)
    {
        std::cout << "\nFormatting completed successfully" << "\n\n";
    }

    MPI_Finalize();

    return 0;
}


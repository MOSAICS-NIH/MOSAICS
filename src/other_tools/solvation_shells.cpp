
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
#include "../headers/common_routines_mpi.h"
#include "../headers/common_routines.h"
#include "../headers/binding_events.h"
#include "../headers/grid_lt.h"
#include "../headers/fit.h"
#include "../headers/command_line_args_mpi.h"
#include "../headers/array.h"
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
    string base_file_name_o;            //Base file name of output files
    string out_file_name;               //Name of the output files with the grid point selection and the big mask
    string prot_file_name;              //Name of the protein mask file
    string mask_file_name;              //Name of the mask file
    FILE *first_shell_file;             //File for writing first shell binding events
    FILE *second_shell_file;            //File for writing second shell binding events
    FILE *third_shell_file;             //File for writing third shell binding events
    FILE *fourth_shell_file;            //File for writing fourth shell binding events
    FILE *fifth_shell_file;             //File for writing fifth shell binding events
    int i                 = 0;          //General variable used in loops
    int j                 = 0;          //General variable used in loops
    int k                 = 0;          //General variable used in loops
    int l                 = 0;          //General variable used in loops
    int m                 = 0;          //General variable used in loops
    int n                 = 0;          //General variable used in loops
    int world_size        = 0;          //Size of the mpi world
    int world_rank        = 0;          //Rank in the mpi world
    int odf               = 0;          //Data file format
    int stride            = 0;          //Skip stride frames
    int begin             = 0;          //Start on this frame
    int end               = 0;          //End on this frame
    int b_end             = 0;          //Did the user provide an end frame?
    int range_shell       = 0;          //How many frames on each side to analyze for assinging shells?
    int range_box         = 0;          //How many frames on each side to analyze for placing lipids in the box?
    int range_big         = 0;          //Which range value is biggest?;
    int target_x          = 0;          //The target grid point x when making a rectangular selection
    int target_y          = 0;          //The target grid point y when making a rectangular selection
    int range_x           = 0;          //The half width of x in the rectangular selection
    int range_y           = 0;          //The half width of y in the rectangular selection
    int invert            = 0;          //Invert rectangular selection (select everything outside rectangle)
    int b_mask            = 0;          //Did the user provide a mask?
    double cell_size      = 1;          //Distance between grid points
    double dt             = 0;          //Time step used for converting frames to time. set equal to ef_dt
    double cutoff_1       = 1;          //Cutoff distance for first shell lipids
    double cutoff_n       = 1;          //Cutoff distance for other shells

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
    string program_name = "Solvation Shells";
    if(world_rank == 0)
    {
        print_credits(argc,argv,program_name);
    }

    string program_description = "Solvation Shells is an analysis tool designed for characterizing the lipid dynamics within each solvation shell around the protein. This tool uses the binding events files produced by 2D Kinetics.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments_mpi(argc,argv,world_rank,program_description);
    add_argument_mpi_s(argc,argv,"-d"     , base_file_name_i,           "Base filename for input binding events files"                          , world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-prot"  , prot_file_name,             "Input file with protein mask (dat)"                                    , world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-mask"  , mask_file_name,             "Input file with selection mask (overrides rectangular selection, dat)" , world_rank, &b_mask,      0);
    add_argument_mpi_s(argc,argv,"-o"     , base_file_name_o,           "Base filename for output data files"                                   , world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-odf"   , &odf,                       "Data file format (0:matrix 1:vector)"                                  , world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-stride", &stride,                    "Skip stride frames"                                                    , world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-b"     , &begin,                     "Start at this frame"                                                   , world_rank, &b_end,       0);
    add_argument_mpi_i(argc,argv,"-e"     , &end,                       "End on this frame"                                                     , world_rank, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-dist_1", &cutoff_1,                  "Cutoff distance for first shell lipids (nm)"                           , world_rank, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-dist_n", &cutoff_n,                  "Cutoff distance for first outer shells (nm)"                           , world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-n_shel", &range_shell,               "Noise filter half width for assigning lipids to the shells? (frames)"  , world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-n_box" , &range_box,                 "Noise filter half width for assigning lipids to the box? (frames)"     , world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-x"     , &target_x,                  "Rectangle center x (grid point)"                                       , world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-y"     , &target_y,                  "Rectangle center y (grid point)"                                       , world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-rx"    , &range_x,                   "Rectangle half width x (grid points)  "                                , world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-ry"    , &range_y,                   "Rectangle half width y (grid points)  "                                , world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-invert", &invert,                    "Invert rectangular selection? (0:no 1:yes)"                            , world_rank, nullptr,      1);
    conclude_input_arguments_mpi(argc,argv,world_rank,program_name);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension_mpi(world_rank,"-prot",prot_file_name,".dat");
    if(b_mask == 1)
    {
        check_extension_mpi(world_rank,"-mask",mask_file_name,".dat");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Get the bigger range value                                                                                //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(range_shell > range_box)
    {
        range_big = range_shell;
    }
    else 
    {
        range_big = range_box; 
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read in protein masking data                                                                              //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Grid_lt prot_mask;
    prot_mask.set_format(odf);
    prot_mask.get_grid(prot_file_name);
    if(world_rank == 0)
    {
        prot_mask.print_dim(1);
        printf("\n");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read in masking data overriding rectangular selection                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Grid_lt mask;
    if(b_mask == 1)
    {
        mask.set_format(odf);
        mask.get_grid(mask_file_name);
        if(world_rank == 0)
        {
            mask.print_dim(0);
            printf("\n");
        }
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

    cell_size = sqrt(events_ref.APS);
    dt        = events_ref.ef_dt; //use units of ps

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
    // Set the effective number of frames taking into account the stride                                         //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int ef_frames = (int)events_ref.ef_frames/stride + 1;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Allocate memory to hold residue names and numbers                                                         //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    iv1d lipid_nr(0,0);
    iv1d res_nr(0,0);
    sv1d res_name(0);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read binding events files until a complete list of lipids is found                                        //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(world_rank == 0)
    {
        events_ref.get_complete_set(base_file_name_i,lipid_nr,res_nr,res_name);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Get lipid blobs                                                                                           //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    events_ref.get_blobs(base_file_name_i,my_xi,my_xf,my_num_g_x,stride,ef_frames,world_rank);

    MPI_Barrier(MPI_COMM_WORLD);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Open binding events files and print header information                                                    //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //create file names for binding events
    string first_shell_file_name  = base_file_name_o + "_first_shell.be";
    string second_shell_file_name = base_file_name_o + "_second_shell.be";
    string third_shell_file_name  = base_file_name_o + "_third_shell.be";
    string fourth_shell_file_name = base_file_name_o + "_fourth_shell.be";
    string fifth_shell_file_name  = base_file_name_o + "_fifth_shell.be";

    if(world_rank == 0)
    {
        //open file for writing 
        first_shell_file  = fopen(first_shell_file_name.c_str(), "w");
        second_shell_file = fopen(second_shell_file_name.c_str(), "w");
        third_shell_file  = fopen(third_shell_file_name.c_str(), "w");
        fourth_shell_file = fopen(fourth_shell_file_name.c_str(), "w");
        fifth_shell_file  = fopen(fifth_shell_file_name.c_str(), "w");

        //print header first shell
        if(first_shell_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",first_shell_file_name.c_str());
        }
        else
        {
            fprintf(first_shell_file," x_i %10d y_i %10d ef_dt(ps) %10f ef_frames %d num_lipids %10d num_g_x %10d num_g_y %10d APS(nm^2) %10f \n\n",-1,-1,events_ref.ef_dt,events_ref.ef_frames,events_ref.num_lipids,events_ref.num_g_x,events_ref.num_g_y,events_ref.APS);
            fprintf(first_shell_file," %10s %10s %10s %15s %15s %20s \n","lipid","res_nr","res_name","bind_i(frame)","bind_f(frame)","dwell time(frames)");
            fprintf(first_shell_file," %10s-%10s-%10s-%15s-%15s-%20s \n","----------","----------","----------","---------------","---------------","--------------------");
        }

        //print header second shell
        if(second_shell_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",second_shell_file_name.c_str());
        }
        else
        {
            fprintf(second_shell_file," x_i %10d y_i %10d ef_dt(ps) %10f ef_frames %d num_lipids %10d num_g_x %10d num_g_y %10d APS(nm^2) %10f \n\n",-1,-1,events_ref.ef_dt,events_ref.ef_frames,events_ref.num_lipids,events_ref.num_g_x,events_ref.num_g_y,events_ref.APS);
            fprintf(second_shell_file," %10s %10s %10s %15s %15s %20s \n","lipid","res_nr","res_name","bind_i(frame)","bind_f(frame)","dwell time(frames)");
            fprintf(second_shell_file," %10s-%10s-%10s-%15s-%15s-%20s \n","----------","----------","----------","---------------","---------------","--------------------");
        }

        //print header third shell
        if(third_shell_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",third_shell_file_name.c_str());
        }
        else
        {
            fprintf(third_shell_file," x_i %10d y_i %10d ef_dt(ps) %10f ef_frames %d num_lipids %10d num_g_x %10d num_g_y %10d APS(nm^2) %10f \n\n",-1,-1,events_ref.ef_dt,events_ref.ef_frames,events_ref.num_lipids,events_ref.num_g_x,events_ref.num_g_y,events_ref.APS);
            fprintf(third_shell_file," %10s %10s %10s %15s %15s %20s \n","lipid","res_nr","res_name","bind_i(frame)","bind_f(frame)","dwell time(frames)");
            fprintf(third_shell_file," %10s-%10s-%10s-%15s-%15s-%20s \n","----------","----------","----------","---------------","---------------","---------------");
        }

        //print header fourth shell
        if(fourth_shell_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",fourth_shell_file_name.c_str());
        }
        else
        {
            fprintf(fourth_shell_file," x_i %10d y_i %10d ef_dt(ps) %10f ef_frames %d num_lipids %10d num_g_x %10d num_g_y %10d APS(nm^2) %10f \n\n",-1,-1,events_ref.ef_dt,events_ref.ef_frames,events_ref.num_lipids,events_ref.num_g_x,events_ref.num_g_y,events_ref.APS);
            fprintf(fourth_shell_file," %10s %10s %10s %15s %15s %20s \n","lipid","res_nr","res_name","bind_i(frame)","bind_f(frame)","dwell time(frames)");
            fprintf(fourth_shell_file," %10s-%10s-%10s-%15s-%15s-%20s \n","----------","----------","----------","---------------","---------------","--------------------");
        }

        //print header fifth shell
        if(fifth_shell_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",fifth_shell_file_name.c_str());
        }
        else
        {
            fprintf(fifth_shell_file," x_i %10d y_i %10d ef_dt(ps) %10f ef_frames %d num_lipids %10d num_g_x %10d num_g_y %10d APS(nm^2) %10f \n\n",-1,-1,events_ref.ef_dt,events_ref.ef_frames,events_ref.num_lipids,events_ref.num_g_x,events_ref.num_g_y,events_ref.APS);
            fprintf(fifth_shell_file," %10s %10s %10s %15s %15s %20s \n","lipid","res_nr","res_name","bind_i(frame)","bind_f(frame)","dwell time(frames)");
            fprintf(fifth_shell_file," %10s-%10s-%10s-%15s-%15s-%20s \n","----------","----------","----------","---------------","---------------","--------------------");
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Collect data and write grids to output                                                                    //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    iv2d shells(2*range_big + 1,iv1d(events_ref.num_lipids,0));                                //holds the shell assignment for each lipid over the window
    iv2d box_window(2*range_big + 1,iv1d(events_ref.num_lipids));                              //holds the box assignment for each lipid over the window
    iv3d grid_window(2*range_big + 1,iv2d(events_ref.num_g_y,iv1d(events_ref.num_g_x,0)));     //holds the grids over the window
    iv3d nan_window(2*range_big + 1,iv2d(events_ref.num_g_y,iv1d(events_ref.num_g_x,0)));      //holds the nan data for the grids over the window
    iv2d bound(events_ref.num_lipids,iv1d(5,0));                                               //holds the bound state for each lipid and each shell
    iv2d age(2*range_big + 1,iv1d(events_ref.num_lipids,0));                                   //holds the age of each entry in the window

    if(world_rank == 0)
    {
        printf("Collecting grid data and assigning lipids to solvation shells. \n");
    }

    for(i=0; i<events_ref.ef_frames; i+=stride)
    {
        if(i >= begin && i <= end) //check start and end frame condition
        {
            if(world_rank == 0)
            {
                printf("Working on frame %d \n",i);
            }

            int this_frame = (int)(i/stride);

            events_ref.get_blobs_frame(this_frame,world_size,world_rank);

            if(world_rank == 0)
            {
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // Find first shell lipids                                                                                   //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                iv1d first_shell(0,0);

                for(j=0; j<events_ref.num_g_y; j++) //loop over y
                {
                    for(k=0; k<events_ref.num_g_x; k++) //loop over x
                    {
                        if(events_ref.blob_nan_frame[j][k] == 0) //check nan
                        {
                            int present = 0;

                            for(l=0; l<first_shell.size(); l++) //loop over first shell lipids
                            {
                                if(events_ref.blob_frame[j][k] == first_shell[l])
                                {
                                    present = 1;
                                }
                            }

                            if(present == 0)
                            {
                                int upper_x = k + (int)ceil(cutoff_1/cell_size);
                                int lower_x = k - (int)ceil(cutoff_1/cell_size);
                                int upper_y = j + (int)ceil(cutoff_1/cell_size);
                                int lower_y = j - (int)ceil(cutoff_1/cell_size);

                                if(upper_x > events_ref.num_g_x)
                                {
                                    upper_x = events_ref.num_g_x;
                                }
                                if(lower_x < 0)
                                {
                                    lower_x = 0;
                                }
                                if(upper_y > events_ref.num_g_y)
                                {
                                    upper_y = events_ref.num_g_y;
                                }
                                if(lower_y < 0)
                                {
                                    lower_y = 0;
                                }
 
                                for(l=lower_x; l<upper_x; l++) //loop over x
                                {
                                    for(m=lower_y; m<upper_y; m++) //loop over y 
                                    {
                                        if(prot_mask.grid[l][m][2][0] == 1) //check if grid point belongs to protein
                                        {
                                            double dx = (k - l)*cell_size;
                                            double dy = (j - m)*cell_size;

                                            double dist = sqrt(dx*dx + dy*dy);

                                            if(dist < cutoff_1)
                                            {
                                                first_shell.push_back(events_ref.blob_frame[j][k]);
                                                goto end_loop;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        end_loop:;
                    }
                }

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // Find second shell lipids                                                                                  //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                iv1d second_shell(0,0);

                for(j=0; j<events_ref.num_g_y; j++) //loop over y
                {
                    for(k=0; k<events_ref.num_g_x; k++) //loop over x
                    {
                        if(events_ref.blob_nan_frame[j][k] == 0) //check nan
                        {
                            int present = 0;

                            for(l=0; l<first_shell.size(); l++) //loop over first shell lipids
                            {
                                if(events_ref.blob_frame[j][k] == first_shell[l])
                                {
                                    present = 1;
                                }
                            }

                            for(l=0; l<second_shell.size(); l++) //loop over second shell lipids
                            {
                                if(events_ref.blob_frame[j][k] == second_shell[l])
                                {
                                    present = 1;
                                }
                            }

                            if(present == 0)
                            {
                                int upper_x = k + (int)ceil(cutoff_1/cell_size);
                                int lower_x = k - (int)ceil(cutoff_1/cell_size);
                                int upper_y = j + (int)ceil(cutoff_1/cell_size);
                                int lower_y = j - (int)ceil(cutoff_1/cell_size);

                                if(upper_x > events_ref.num_g_x)
                                {
                                    upper_x = events_ref.num_g_x;
                                }
                                if(lower_x < 0)
                                {
                                    lower_x = 0;
                                }
                                if(upper_y > events_ref.num_g_y)
                                {
                                    upper_y = events_ref.num_g_y;
                                }
                                if(lower_y < 0)
                                {
                                    lower_y = 0;
                                }

                                for(l=lower_x; l<upper_x; l++) //loop over x
                                {
                                    for(m=lower_y; m<upper_y; m++) //loop over y
                                    {
                                        if(events_ref.blob_nan_frame[m][l] == 0) //grid point contains a lipid
                                        {
                                            for(n=0; n<first_shell.size(); n++) //loop over first shell lipids
                                            {
                                                if(events_ref.blob_frame[m][l] == first_shell[n]) //check if grid point belongs to a first shell lipid
                                                {
                                                    double dx = (k - l)*cell_size;
                                                    double dy = (j - m)*cell_size;

                                                    double dist = sqrt(dx*dx + dy*dy);

                                                    if(dist < cutoff_n)
                                                    {
                                                        second_shell.push_back(events_ref.blob_frame[j][k]);
                                                        goto end_loop_1;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        end_loop_1:;
                    }
                }

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // Find third shell lipids                                                                                  //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                iv1d third_shell(0,0);

                for(j=0; j<events_ref.num_g_y; j++) //loop over y
                {
                    for(k=0; k<events_ref.num_g_x; k++) //loop over x
                    {
                        if(events_ref.blob_nan_frame[j][k] == 0) //check nan
                        {
                            int present = 0;

                            for(l=0; l<first_shell.size(); l++) //loop over first shell lipids
                            {
                                if(events_ref.blob_frame[j][k] == first_shell[l])
                                {
                                    present = 1;
                                }
                            }

                            for(l=0; l<second_shell.size(); l++) //loop over second shell lipids
                            {
                                if(events_ref.blob_frame[j][k] == second_shell[l])
                                {
                                    present = 1;
                                }
                            }

                            for(l=0; l<third_shell.size(); l++) //loop over third shell lipids
                            {
                                if(events_ref.blob_frame[j][k] == third_shell[l])
                                {
                                    present = 1;
                                }
                            }

                            if(present == 0)
                            {
                                int upper_x = k + (int)ceil(cutoff_1/cell_size);
                                int lower_x = k - (int)ceil(cutoff_1/cell_size);
                                int upper_y = j + (int)ceil(cutoff_1/cell_size);
                                int lower_y = j - (int)ceil(cutoff_1/cell_size);

                                if(upper_x > events_ref.num_g_x)
                                {
                                    upper_x = events_ref.num_g_x;
                                }
                                if(lower_x < 0)
                                {
                                    lower_x = 0;
                                }
                                if(upper_y > events_ref.num_g_y)
                                {
                                    upper_y = events_ref.num_g_y;
                                }
                                if(lower_y < 0)
                                {
                                    lower_y = 0;
                                }

                                for(l=lower_x; l<upper_x; l++) //loop over x
                                {
                                    for(m=lower_y; m<upper_y; m++) //loop over y
                                    {
                                        if(events_ref.blob_nan_frame[m][l] == 0) //grid point contains a lipid
                                        {
                                            for(n=0; n<second_shell.size(); n++) //loop over second shell lipids
                                            {
                                                if(events_ref.blob_frame[m][l] == second_shell[n]) //check if grid point belongs to a second shell lipid
                                                {
                                                    double dx = (k - l)*cell_size;
                                                    double dy = (j - m)*cell_size;

                                                    double dist = sqrt(dx*dx + dy*dy);

                                                    if(dist < cutoff_n)
                                                    {
                                                        third_shell.push_back(events_ref.blob_frame[j][k]);
                                                        goto end_loop_2;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        end_loop_2:;
                    }
                }

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // Find fourth shell lipids                                                                                  //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                iv1d fourth_shell(0,0);

                for(j=0; j<events_ref.num_g_y; j++) //loop over y
                {
                    for(k=0; k<events_ref.num_g_x; k++) //loop over x
                    {
                        if(events_ref.blob_nan_frame[j][k] == 0) //check nan
                        {
                            int present = 0;

                            for(l=0; l<first_shell.size(); l++) //loop over first shell lipids
                            {
                                if(events_ref.blob_frame[j][k] == first_shell[l])
                                {
                                    present = 1;
                                }
                            }

                            for(l=0; l<second_shell.size(); l++) //loop over second shell lipids
                            {
                                if(events_ref.blob_frame[j][k] == second_shell[l])
                                {
                                    present = 1;
                                }
                            }

                            for(l=0; l<third_shell.size(); l++) //loop over third shell lipids
                            {
                                if(events_ref.blob_frame[j][k] == third_shell[l])
                                {
                                    present = 1;
                                }
                            }

                            for(l=0; l<fourth_shell.size(); l++) //loop over fourth shell lipids
                            {
                                if(events_ref.blob_frame[j][k] == fourth_shell[l])
                                {
                                    present = 1;
                                }
                            }

                            if(present == 0)
                            {
                                int upper_x = k + (int)ceil(cutoff_1/cell_size);
                                int lower_x = k - (int)ceil(cutoff_1/cell_size);
                                int upper_y = j + (int)ceil(cutoff_1/cell_size);
                                int lower_y = j - (int)ceil(cutoff_1/cell_size);

                                if(upper_x > events_ref.num_g_x)
                                {
                                    upper_x = events_ref.num_g_x;
                                }
                                if(lower_x < 0)
                                {
                                    lower_x = 0;
                                }
                                if(upper_y > events_ref.num_g_y)
                                {
                                    upper_y = events_ref.num_g_y;
                                }
                                if(lower_y < 0)
                                {
                                    lower_y = 0;
                                }

                                for(l=lower_x; l<upper_x; l++) //loop over x
                                {
                                    for(m=lower_y; m<upper_y; m++) //loop over y
                                    {
                                        if(events_ref.blob_nan_frame[m][l] == 0) //grid point contains a lipid
                                        {
                                            for(n=0; n<third_shell.size(); n++) //loop over third shell lipids
                                            {
                                                if(events_ref.blob_frame[m][l] == third_shell[n]) //check if grid point belongs to a third shell lipid
                                                {
                                                    double dx = (k - l)*cell_size;
                                                    double dy = (j - m)*cell_size;

                                                    double dist = sqrt(dx*dx + dy*dy);

                                                    if(dist < cutoff_n)
                                                    {
                                                        fourth_shell.push_back(events_ref.blob_frame[j][k]);
                                                        goto end_loop_3;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        end_loop_3:;
                    }
                }

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // Find fifth shell lipids                                                                                   //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                iv1d fifth_shell(0,0);

                for(j=0; j<events_ref.num_g_y; j++) //loop over y
                {
                    for(k=0; k<events_ref.num_g_x; k++) //loop over x
                    {
                        if(events_ref.blob_nan_frame[j][k] == 0) //check nan
                        {
                            int present = 0;

                            for(l=0; l<first_shell.size(); l++) //loop over first shell lipids
                            {
                                if(events_ref.blob_frame[j][k] == first_shell[l])
                                {
                                    present = 1;
                                }
                            }

                            for(l=0; l<second_shell.size(); l++) //loop over second shell lipids
                            {
                                if(events_ref.blob_frame[j][k] == second_shell[l])
                                {
                                    present = 1;
                                }
                            }

                            for(l=0; l<third_shell.size(); l++) //loop over third shell lipids
                            {
                                if(events_ref.blob_frame[j][k] == third_shell[l])
                                {
                                    present = 1;
                                }
                            }

                            for(l=0; l<fourth_shell.size(); l++) //loop over fourth shell lipids
                            {
                                if(events_ref.blob_frame[j][k] == fourth_shell[l])
                                {
                                    present = 1;
                                }
                            }

                            for(l=0; l<fifth_shell.size(); l++) //loop over fifth shell lipids
                            {
                                if(events_ref.blob_frame[j][k] == fifth_shell[l])
                                {
                                    present = 1;
                                }
                            }

                            if(present == 0)
                            {
                                int upper_x = k + (int)ceil(cutoff_1/cell_size);
                                int lower_x = k - (int)ceil(cutoff_1/cell_size);
                                int upper_y = j + (int)ceil(cutoff_1/cell_size);
                                int lower_y = j - (int)ceil(cutoff_1/cell_size);

                                if(upper_x > events_ref.num_g_x)
                                {
                                    upper_x = events_ref.num_g_x;
                                }
                                if(lower_x < 0)
                                {
                                    lower_x = 0;
                                }
                                if(upper_y > events_ref.num_g_y)
                                {
                                    upper_y = events_ref.num_g_y;
                                }
                                if(lower_y < 0)
                                {
                                    lower_y = 0;
                                }

                                for(l=lower_x; l<upper_x; l++) //loop over x
                                {
                                    for(m=lower_y; m<upper_y; m++) //loop over y
                                    {
                                        if(events_ref.blob_nan_frame[m][l] == 0) //grid point contains a lipid
                                        {
                                            for(n=0; n<fourth_shell.size(); n++) //loop over fourth shell lipids
                                            {
                                                if(events_ref.blob_frame[m][l] == fourth_shell[n]) //check if grid point belongs to a fourth shell lipid
                                                {
                                                    double dx = (k - l)*cell_size;
                                                    double dy = (j - m)*cell_size;

                                                    double dist = sqrt(dx*dx + dy*dy);

                                                    if(dist < cutoff_n)
                                                    {
                                                        fifth_shell.push_back(events_ref.blob_frame[j][k]);
                                                        goto end_loop_4;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        end_loop_4:;
                    }
                }

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // Assign other lipids                                                                                       //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                iv1d other_shells(0,0);

                for(j=0; j<events_ref.num_lipids; j++) //loop over lipids
                {
                    int present = 0;

                    for(k=0; k<first_shell.size(); k++)
                    {
                        if(first_shell[k] == j)
                        {
                            present = 1;
                        }
                    }

                    for(k=0; k<second_shell.size(); k++)
                    {
                        if(second_shell[k] == j)
                        {
                            present = 1;
                        }
                    }

                    for(k=0; k<third_shell.size(); k++)
                    {
                        if(third_shell[k] == j)
                        {
                            present = 1;
                        }
                    }

                    for(k=0; k<fourth_shell.size(); k++)
                    {
                        if(fourth_shell[k] == j)
                        {
                            present = 1;
                        }
                    }

                    for(k=0; k<fifth_shell.size(); k++)
                    {
                        if(fifth_shell[k] == j)
                        {
                            present = 1;
                        }
                    }

                    if(present == 0)
                    {
                        other_shells.push_back(j);
                    }
                }

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // Store assignments                                                                                         //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                int pos = this_frame%(2*range_big + 1);
                
                for(j=0; j<events_ref.num_lipids; j++)
                {
                    age[pos][j] = 0;
                    for(k=0; k<2*range_big + 1; k++) //loop over window
                    {
                        age[k][j] = age[k][j] + 1;
                    }
 
                    for(k=0; k<first_shell.size(); k++)
                    {
                        if(first_shell[k] == j)
                        {
                            shells[pos][j] = 1;
                            goto end_loop_5;
                        }
                    }

                    for(k=0; k<second_shell.size(); k++)
                    {
                        if(second_shell[k] == j)
                        {
                            shells[pos][j] = 2;
                            goto end_loop_5;
                        }
                    }

                    for(k=0; k<third_shell.size(); k++)
                    {
                        if(third_shell[k] == j)
                        {
                            shells[pos][j] = 3;
                            goto end_loop_5;
                        }
                    }

                    for(k=0; k<fourth_shell.size(); k++)
                    {
                        if(fourth_shell[k] == j)
                        {
                            shells[pos][j] = 4;
                            goto end_loop_5;
                        }
                    }

                    for(k=0; k<fifth_shell.size(); k++)
                    {
                        if(fifth_shell[k] == j)
                        {
                            shells[pos][j] = 5;
                            goto end_loop_5;
                        }
                    }

                    for(k=0; k<other_shells.size(); k++)
                    {
                        if(other_shells[k] == j)
                        {
                            shells[pos][j] = 6;
                            goto end_loop_5;
                        }
                    }
                    end_loop_5:;
                }
                grid_window[pos] = events_ref.blob_frame;
                nan_window[pos]  = events_ref.blob_nan_frame; 

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // Find the center of each lipid                                                                             //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                dv1d centers_x(events_ref.num_lipids,0.0);
                dv1d centers_y(events_ref.num_lipids,0.0);
                iv1d centers_count(events_ref.num_lipids,0);

                for(j=0; j<events_ref.num_g_y; j++) //loop over y-dimension
                {
                    for(k=0; k<events_ref.num_g_x; k++) //loop over x-dimension
                    {
                        if(prot_mask.grid[k][j][2][0] == 0) //check that grid point is not the protein
                        {
                            if(events_ref.blob_nan_frame[j][k] == 0)
                            {
                                int lip_index = events_ref.blob_frame[j][k];

                                centers_x[lip_index]     = centers_x[lip_index] + (double)k*cell_size;
                                centers_y[lip_index]     = centers_y[lip_index] + (double)j*cell_size;
                                centers_count[lip_index] = centers_count[lip_index] + 1;
                            }
                        }
                    }
                }

                for(j=0; j<events_ref.num_lipids; j++) //loop over lipids
                {
                    if(centers_count[j] > 0)
                    {
                        centers_x[j] = centers_x[j]/(double)centers_count[j];
                        centers_y[j] = centers_y[j]/(double)centers_count[j];
                    }
                    else
                    {
                        centers_x[j] = -999999;
                        centers_y[j] = -999999;
                    }
                }

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // Determine which lipids are inside or outside the rectangular selection/mask                               //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                for(j=0; j<events_ref.num_lipids; j++)
                {
                    int closest_x = (int)ceil(centers_x[j]/cell_size);
                    int closest_y = (int)ceil(centers_y[j]/cell_size);
                    int rectangular_selection_pass = 0;

                    if(centers_x[j] == -999999.0 && centers_y[j] == -999999.0) //first frame may have no data (empty grid)
                    {
                        closest_x = -999999;
                        closest_y = -999999;
                    }
                    if(b_mask == 0) //no mask data was provided. use rectangular selection
                    {
                        rectangular_selection_pass = check_rectangular_pass_vector(invert,target_x,target_y,range_x,range_y,closest_x,closest_y);
                    }
                    else if(closest_x >= 0 && closest_x < mask.size_x() && closest_y >= 0 && closest_y < mask.size_y()) //center is inside the grid
                    {
                        if(mask.grid[closest_x][closest_y][2][0] == 1) //grid point is inside the mask
                        {
                            rectangular_selection_pass = 1;
                        }
                    }

                    box_window[pos][j] = rectangular_selection_pass;
                }

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // Analyze assignments over window                                                                           //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                iv1d final_assignment(events_ref.num_lipids,0);
                iv1d final_assignment_box(events_ref.num_lipids,0);

                if(this_frame >= 2*range_big + 1)
                {
                    for(j=0; j<events_ref.num_lipids; j++)
                    { 
                        int delta_shell = range_big - range_shell; 
                        int delta_box   = range_big - range_box; 

                        int count_1   = 0;
                        int count_2   = 0;
                        int count_3   = 0;
                        int count_4   = 0;
                        int count_5   = 0;
                        int count_6   = 0;
                        int count_box = 0;

                        //count frequency of shells and being in box
                        for(k=0; k<2*range_big + 1; k++)
                        {
                            if(age[k][j] > delta_shell && age[k][j] <= (2*range_big + 1) - delta_shell)
                            {
                                if(shells[k][j] == 1)
                                {
                                    count_1 = count_1 + 1;
                                }
                                else if(shells[k][j] == 2)
                                {
                                    count_2 = count_2 + 1;
                                }
                                else if(shells[k][j] == 3)
                                {
                                    count_3 = count_3 + 1;
                                }
                                else if(shells[k][j] == 4)
                                {
                                    count_4 = count_4 + 1;
                                }
                                else if(shells[k][j] == 5)
                                {
                                    count_5 = count_5 + 1;
                                }
                                else if(shells[k][j] == 6)
                                {
                                    count_6 = count_6 + 1;
                                }
                            }
                            
                            if(age[k][j] > delta_box && age[k][j] <= (2*range_big + 1) - delta_box)
                            {
                                if(box_window[k][j] == 1)
                                {
                                    count_box = count_box + 1;
                                }
                            }
                        }

                        final_assignment[j] = j;

                        if((double)count_box/(double)(2*range_box + 1) > 0.5)
                        {
                            final_assignment_box[j] = 1;
                        }

                        if(count_1 >= count_2 && count_1 >= count_3 && count_1 >= count_4 && count_1 >= count_5 && count_1 >= count_6 && final_assignment_box[j] == 1)
                        {
                            final_assignment[j] = -1;
                        }
                        else if(count_2 >= count_1 && count_2 >= count_3 && count_2 >= count_4 && count_2 >= count_5 && count_2 >= count_6 && final_assignment_box[j] == 1)
                        {
                            final_assignment[j] = -2;
                        }
                        else if(count_3 >= count_1 && count_3 >= count_2 && count_3 >= count_4 && count_3 >= count_5 && count_3 >= count_6 && final_assignment_box[j] == 1)
                        {
                            final_assignment[j] = -3;
                        }
                        else if(count_4 >= count_1 && count_4 >= count_2 && count_4 >= count_3 && count_4 >= count_5 && count_4 >= count_6 && final_assignment_box[j] == 1)
                        {
                            final_assignment[j] = -4;
                        }
                        else if(count_5 >= count_1 && count_5 >= count_2 && count_5 >= count_3 && count_5 >= count_4 && count_5 >= count_6 && final_assignment_box[j] == 1)
                        {
                            final_assignment[j] = -5;
                        }
                        else if(count_6 >= count_1 && count_6 >= count_2 && count_6 >= count_3 && count_6 >= count_4 && count_6 >= count_5 && final_assignment_box[j] == 1)
                        {
                            final_assignment[j] = -6;
                        }
                    }

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //                                                                                                           //
                    // Check for tansitions between shells and record dwell times                                                //
                    //                                                                                                           //
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    for(j=0; j<events_ref.num_lipids; j++) //loop over lipids
                    {
                        if(final_assignment[j] == -1) //first shell
                        {
                            bound[j][0] = bound[j][0] + 1;
                        }
                        else
                        {
                            if(bound[j][0] > 0) //record dwell time
                            {
                                int dwell_t = bound[j][0]*stride; 
                                int bind_f  = stride*(this_frame - range_big);
                                int bind_i  = bind_f - dwell_t;
 
                                //report binding event
                                fprintf(first_shell_file," %10d %10d %10s %15d %15d %20d \n",j,res_nr[j],res_name[j].c_str(),bind_i,bind_f,dwell_t);

                                bound[j][0] = 0;   
                            }
                        }

                        if(final_assignment[j] == -2) //second shell
                        {
                            bound[j][1] = bound[j][1] + 1;
                        }
                        else
                        {
                            if(bound[j][1] > 0) //record dwell time
                            {
                                int dwell_t = bound[j][1]*stride;
                                int bind_f  = stride*(this_frame - range_big);
                                int bind_i  = bind_f - dwell_t;

                                //report binding event
                                fprintf(second_shell_file," %10d %10d %10s %15d %15d %20d \n",j,res_nr[j],res_name[j].c_str(),bind_i,bind_f,dwell_t);

                                bound[j][1] = 0;
                            }
                        }

                        if(final_assignment[j] == -3) //third shell
                        {
                            bound[j][2] = bound[j][2] + 1;
                        }
                        else
                        {
                            if(bound[j][2] > 0) //record dwell time
                            {
                                int dwell_t = bound[j][2]*stride;
                                int bind_f  = stride*(this_frame - range_big);
                                int bind_i  = bind_f - dwell_t;

                                //report binding event
                                fprintf(third_shell_file," %10d %10d %10s %15d %15d %20d \n",j,res_nr[j],res_name[j].c_str(),bind_i,bind_f,dwell_t);

                                bound[j][2] = 0;
                            }
                        }

                        if(final_assignment[j] == -4) //fourth shell
                        {
                            bound[j][3] = bound[j][3] + 1;
                        }
                        else
                        {
                            if(bound[j][3] > 0) //record dwell time
                            {
                                int dwell_t = bound[j][3]*stride;
                                int bind_f  = stride*(this_frame - range_big);
                                int bind_i  = bind_f - dwell_t;

                                //report binding event
                                fprintf(fourth_shell_file," %10d %10d %10s %15d %15d %20d \n",j,res_nr[j],res_name[j].c_str(),bind_i,bind_f,dwell_t);

                                bound[j][3] = 0;
                            }
                        }

                        if(final_assignment[j] == -5) //fifth shell
                        {
                            bound[j][4] = bound[j][4] + 1;
                        }
                        else
                        {
                            if(bound[j][4] > 0) //record dwell time
                            {
                                int dwell_t = bound[j][4]*stride;
                                int bind_f  = stride*(this_frame - range_big);
                                int bind_i  = bind_f - dwell_t;

                                //report binding event
                                fprintf(fifth_shell_file," %10d %10d %10s %15d %15d %20d \n",j,res_nr[j],res_name[j].c_str(),bind_i,bind_f,dwell_t);

                                bound[j][4] = 0;
                            }
                        }
                    }

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //                                                                                                           //
                    // Get the old grid and highlight shells selection                                                           //
                    //                                                                                                           //
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    int prev_pos = (this_frame-range_big)%(2*range_big + 1);

                    events_ref.blob_frame = grid_window[prev_pos]; 
                    events_ref.blob_nan_frame  = nan_window[prev_pos];
 
                    for(j=0; j<events_ref.num_g_x; j++) //loop over x
                    {
                        for(k=0; k<events_ref.num_g_y; k++) //loop over y
                        {
                            if(events_ref.blob_nan_frame[k][j] == 0)
                            {
                                events_ref.blob_frame[k][j] = final_assignment[events_ref.blob_frame[k][j]];
                            }
                        }
                    }

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //                                                                                                           //
                    // Write grid to file                                                                                        //
                    //                                                                                                           //
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //set nan flags for the protein lattice points
                    for(j=0; j<events_ref.num_g_y; j++) //loop over y-dimension
                    {
                        for(k=0; k<events_ref.num_g_x; k++) //loop over x-dimension
                        {
                            if(prot_mask.grid[k][j][2][0] == 1)
                            {
                                events_ref.blob_nan_frame[j][k] = 1;
                            }
                        }
                    }

                    out_file_name = base_file_name_o + "_" + to_string(stride*(this_frame - range_big)) + ".dat";
                    events_ref.write_blobs_frame(out_file_name);
                }
            }

            MPI_Barrier(MPI_COMM_WORLD);            
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Close binding events files                                                                                //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(world_rank == 0)
    {
        fclose(first_shell_file);
        fclose(second_shell_file);
        fclose(third_shell_file);
        fclose(fourth_shell_file);
        fclose(fifth_shell_file);
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


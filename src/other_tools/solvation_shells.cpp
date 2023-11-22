
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
#include "../headers/file_naming.h"
#include "../headers/binding_events.h"
#include "../headers/grid_lt.h"
#include "../headers/fit.h"
#include "../headers/command_line_args_mpi.h"
#include "../headers/array.h"
#include "../headers/file_naming_mpi.h"     
#include "../headers/performance.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function checks if the grid point is an exterior point for the lipid                                 //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int check_exterior(int x,int y,int num_g_x,int num_g_y,iv2d &voro_frame,iv2d &voro_nan_frame)
{
    int i        = 0;      //standard variable used in loops
    int j        = 0;      //standard variable used in loops    
    int exterior = 0;      //is the grid point on the lipid surface or not
    int upper_x  = x + 1;  //upper bound for x
    int lower_x  = x - 1;  //lower bound for x
    int upper_y  = y + 1;  //upper bound for y
    int lower_y  = y - 1;  //lower bound for y

    //check that boundaries do not exceed those of the grid
    if(upper_x > num_g_x - 1)
    {
        upper_x = num_g_x - 1;
    }
    if(lower_x < 0)
    {
        lower_x = 0;
    }
    if(upper_y > num_g_y - 1)
    {
        upper_y = num_g_y - 1;
    }
    if(lower_y < 0)
    {
        lower_y = 0;
    }

    //check if the lattice point is on the surface
    for(i=lower_x; i<=upper_x; i++) //loop over x
    {
        for(j=lower_y; j<=upper_y; j++) //loop over y
        {
            if(voro_frame[j][i] != voro_frame[y][x] || voro_nan_frame[j][i] == 1)
            {
                exterior = 1; 
                goto end_ext;
            }
        }
    }
    end_ext:;

    return exterior; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function gets the upper and lower range for the lattice point                                        //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_range(double cutoff,double cell_size,int *upper_x,int *lower_x,int *upper_y,int *lower_y,int x,int y,int num_g_x,int num_g_y)
{
    *upper_x = x + (int)ceil(cutoff/cell_size);
    *lower_x = x - (int)ceil(cutoff/cell_size);
    *upper_y = y + (int)ceil(cutoff/cell_size);
    *lower_y = y - (int)ceil(cutoff/cell_size);

    if(*upper_x > num_g_x)
    {
        *upper_x = num_g_x;
    }
    if(*lower_x < 0)
    {
        *lower_x = 0;
    }
    if(*upper_y > num_g_y)
    {
        *upper_y = num_g_y;
    }
    if(*lower_y < 0)
    {
        *lower_y = 0;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Main function which executes all other functions                                                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, const char * argv[]) 
{
    //Here we define some variables used throughout
    string in_file_name;                //Name of input file 
    string base_file_name_o;            //Base file name of output files
    string out_file_name;               //Name of the output files with the grid point selection and the big mask
    string prot_file_name;              //Name of the protein mask file
    string mask_file_name;              //Name of the mask file
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
    int counter           = 0;          //How many times the "program run time" been displayed
    int shell_index       = 0;          //used to loop over the shells
    int test              = 0;          //print tessellation data with highlighted shells?
    int dump              = 0;          //dump lipids on the last frame?
    double cell_size      = 1;          //Distance between grid points
    double dt             = 0;          //Time step used for converting frames to time. set equal to ef_dt
    double cutoff_1       = 1;          //Cutoff distance for first shell lipids
    double cutoff_n       = 1;          //Cutoff distance for other shells
    clock_t t;                          //Keeps the time for testing performance
    sv1d cl_tags;                       //Holds a list of command line tags for the program
  
    //NOTES:
    //This program uses two noise filters. One tells if the lipid is inside the shell and another tells if it is 
    //inside the box. The halfwidth of the filters are chosen independent of each other. To enable this, we determine
    //which half-width (range_shell vs range_box) is larger. Then a noise filter is created with this larger range. 
    //The age of each entry in the filter is maintained as shown in the example below. We can then analye a subset of the filter alone for 
    //the smaller filter using the check if(age[k][j] > delta_shell && age[k][j] <= (2*range_big + 1) - delta_shell).
    //where delta_shell is range_big - range_shell and tells the difference in size of the shell filter and the largest filter. 
    //Assignments are only made once the larger filter is full. Still the smaller filter is considered for making the assignments. 

    // age[pos][lipid]
    //          |4 5 1 2 3| 
    //          |4 5 1 2 3|
    // lipids   |4 5 1 2 3|
    //          |4 5 1 2 3|
    //          |4 5 1 2 3|
    //             pos

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
    add_argument_mpi_s(argc,argv,"-d"     , in_file_name,               "Input binding events file (be)"                                        , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-prot"  , prot_file_name,             "Input file with protein mask (dat)"                                    , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-mask"  , mask_file_name,             "Input file with selection mask (overrides rectangular selection, dat)" , world_rank, cl_tags, &b_mask,      0);
    add_argument_mpi_s(argc,argv,"-o"     , base_file_name_o,           "Filename for deriving output data filenames (be)"                      , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-odf"   , &odf,                       "Data file format (0:matrix 1:vector)"                                  , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-stride", &stride,                    "Skip stride frames"                                                    , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-b"     , &begin,                     "Start at this frame"                                                   , world_rank, cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e"     , &end,                       "End on this frame"                                                     , world_rank, cl_tags, &b_end,       0);
    add_argument_mpi_d(argc,argv,"-dist_1", &cutoff_1,                  "Cutoff distance for first shell lipids (nm)"                           , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-dist_n", &cutoff_n,                  "Cutoff distance for first outer shells (nm)"                           , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-n_shel", &range_shell,               "Noise filter half width for assigning lipids to the shells? (frames)"  , world_rank, cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-n_box" , &range_box,                 "Noise filter half width for assigning lipids to the box? (frames)"     , world_rank, cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-x"     , &target_x,                  "Rectangle center x (grid point)"                                       , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-y"     , &target_y,                  "Rectangle center y (grid point)"                                       , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-rx"    , &range_x,                   "Rectangle half width x (grid points)  "                                , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-ry"    , &range_y,                   "Rectangle half width y (grid points)  "                                , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-invert", &invert,                    "Invert rectangular selection? (0:no 1:yes)"                            , world_rank, cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-test",   &test,                      "Write voronoi tessellations with shells highlighted? (0:no 1:yes)"     , world_rank, cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-dump",   &dump,                      "Dump lipids from each shell on the last frame? (0:no 1:yes)"           , world_rank, cl_tags, nullptr,      0);
    conclude_input_arguments_mpi(argc,argv,world_rank,program_name,cl_tags);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension_mpi(world_rank,"-d",in_file_name,".be");
    check_extension_mpi(world_rank,"-o",base_file_name_o,".be");
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
    else
    {
        result = events_ref.get_binding_events_xy(in_file_name,0,0);
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
    // Distribute the workload across the cores for getting tessellations                                        //
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
    // Set the effective number of frames taking into account the stride                                         //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int ef_frames = (int)events_ref.ef_frames/stride + 1;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Allocate memory to hold residue names and numbers                                                         //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    iv1d lipid_nr(0,0);    //holds the lipid number (complete set) 
    iv1d res_nr(0,0);      //holds the residue id (complete set)
    sv1d res_name(0);      //holds the residue names (complete set)

    //log time spent performing misc tasks
    perf.log_time((clock() - t)/CLOCKS_PER_SEC,"Other");

    //reset the clock
    t = clock();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read lipid info file or read binding events files until a complete list of lipids is found                //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(world_rank == 0)
    {
        events_ref.get_complete_set(in_file_name,lipid_nr,res_nr,res_name);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //log time spent getting a complete set
    perf.log_time((clock() - t)/CLOCKS_PER_SEC,"Get Set");

    //reset the clock
    t = clock();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Get lipid tessellations                                                                                   //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    perf.log_time(events_ref.get_tessellations(in_file_name,my_gi,my_gf,my_num_g,stride,ef_frames,world_rank),"Get Tessellations");

    MPI_Barrier(MPI_COMM_WORLD);

    //reset the clock
    t = clock();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Create binding events objects and set the header information                                              //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    vector <Binding_events> be_shells(5);     //store the binding events data for each shell

    //set the header info
    for(i=0; i<5; i++)
    {
        be_shells[i].lipid_nr.resize(0,0);
        be_shells[i].res_nr.resize(0,0);
        be_shells[i].res_name.resize(0);
        be_shells[i].bind_i.resize(0,0);
        be_shells[i].bind_f.resize(0,0);
        be_shells[i].dwell_t.resize(0,0);

    	be_shells[i].x_i          = -1;                                                                                                        
        be_shells[i].y_i          = -1;                                                                                                        
        be_shells[i].ef_frames    = events_ref.ef_frames;                                                                                               
        be_shells[i].num_lipids   = events_ref.num_lipids;                                                                                                
        be_shells[i].num_g_x      = events_ref.num_g_x;                                                                                                   
        be_shells[i].num_g_y      = events_ref.num_g_y;                                                                                                   
        be_shells[i].ef_dt        = events_ref.ef_dt;                                                                                                     
        be_shells[i].APS          = events_ref.APS;  
    }

    //create file names for binding events
    sv1d shell_file_names(5);                //store the filename for each shell's binding events file
    shell_file_names[0] = add_tag(base_file_name_o, "_first_shell");
    shell_file_names[1] = add_tag(base_file_name_o, "_second_shell");
    shell_file_names[2] = add_tag(base_file_name_o, "_third_shell");
    shell_file_names[3] = add_tag(base_file_name_o, "_fourth_shell");
    shell_file_names[4] = add_tag(base_file_name_o, "_fifth_shell");

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
        printf("\nCollecting grid data and assigning lipids to solvation shells. \n");
        printf("-----------------------------------------------------------------------------------------------------------------------------------\n");
    }

    for(i=0; i<events_ref.ef_frames; i+=stride) //loop over frames
    {
        if(i >= begin && i <= end) //check start and end frame condition
        {
            int this_frame = (int)(i/stride);

            events_ref.get_voro_frame(this_frame,world_size,world_rank);

            if(world_rank == 0)
            {
                iv1d current_shell(lipid_nr.size(),6);  //store the shell for each lipid

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // Find first shell lipids                                                                                   //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                for(j=0; j<events_ref.num_g_y; j++) //loop over y
                {
                    for(k=0; k<events_ref.num_g_x; k++) //loop over x
                    {
                        if(events_ref.voro_nan_frame[j][k] == 0) //check nan
                        {
                            if(current_shell[events_ref.voro_frame[j][k]] == 6) //not assigned yet
                            {
                                if(check_exterior(k,j,events_ref.num_g_x,events_ref.num_g_y,events_ref.voro_frame,events_ref.voro_nan_frame) == 1)
                                {
                                    int upper_x, lower_x, upper_y, lower_y; 

                                    get_range(cutoff_1,cell_size,&upper_x,&lower_x,&upper_y,&lower_y,k,j,events_ref.num_g_x,events_ref.num_g_y); 

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
                                                    current_shell[events_ref.voro_frame[j][k]] = 1;
                                                    goto end_loop;
                                                }
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
                // Find remaining shells                                                                                     //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                for(shell_index=2; shell_index<=5; shell_index++) //loop over the shells
                {                
                    int previous_shell = shell_index - 1; 
 
                    for(j=0; j<events_ref.num_g_y; j++) //loop over y
                    {
                        for(k=0; k<events_ref.num_g_x; k++) //loop over x
                        {
                            if(events_ref.voro_nan_frame[j][k] == 0) //check nan
                            {
                                if(current_shell[events_ref.voro_frame[j][k]] == 6) //not assigned yet
                                {
                                    if(check_exterior(k,j,events_ref.num_g_x,events_ref.num_g_y,events_ref.voro_frame,events_ref.voro_nan_frame) == 1)
                                    {
                                        int upper_x, lower_x, upper_y, lower_y;

                                        get_range(cutoff_n,cell_size,&upper_x,&lower_x,&upper_y,&lower_y,k,j,events_ref.num_g_x,events_ref.num_g_y);

                                        for(l=lower_x; l<upper_x; l++) //loop over x
                                        {
                                            for(m=lower_y; m<upper_y; m++) //loop over y
                                            {
                                                if(events_ref.voro_nan_frame[m][l] == 0) //grid point contains a lipid
                                                {
                                                    if(current_shell[events_ref.voro_frame[m][l]] == previous_shell) //a previous shell lipid
                                                    {
                                                        double dx = (k - l)*cell_size;
                                                        double dy = (j - m)*cell_size;

                                                        double dist = sqrt(dx*dx + dy*dy);

                                                        if(dist < cutoff_n)
                                                        {
                                                            current_shell[events_ref.voro_frame[j][k]] = shell_index;
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
                }

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // Store assignments and tessellations and monitor the age of each entry in the window                       //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                int pos = this_frame%(2*range_big + 1);
                
                for(j=0; j<events_ref.num_lipids; j++) //loop over lipids
                {
                    age[pos][j] = 0;
                    for(k=0; k<2*range_big + 1; k++) //loop over window
                    {
                        age[k][j] = age[k][j] + 1;
                    }

                    shells[pos][j] = current_shell[j];
                }
                grid_window[pos] = events_ref.voro_frame;
                nan_window[pos]  = events_ref.voro_nan_frame; 

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // Find the center of each lipid                                                                             //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                dv1d centers_x(events_ref.num_lipids,0.0);        //holds the center x-component for each lipid
                dv1d centers_y(events_ref.num_lipids,0.0);        //holds the center y-component for each lipid
                iv1d centers_count(events_ref.num_lipids,0);      //holds the number of lattice points occupied by each lipid

                for(j=0; j<events_ref.num_g_y; j++) //loop over y-dimension
                {
                    for(k=0; k<events_ref.num_g_x; k++) //loop over x-dimension
                    {
                        if(prot_mask.grid[k][j][2][0] == 0) //check that grid point is not the protein
                        {
                            if(events_ref.voro_nan_frame[j][k] == 0) //not a NaN lattice point
                            {
                                int lip_index = events_ref.voro_frame[j][k];

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
                        centers_x[j] = -999999.9;
                        centers_y[j] = -999999.9;
                    }
                }

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                                                                                                           //
                // Determine which lipids are inside or outside the rectangular selection/mask                               //
                //                                                                                                           //
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                for(j=0; j<events_ref.num_lipids; j++)
                {
                    int closest_x = (int)ceil(centers_x[j]/cell_size);   //find the lattice point in x that is closest to the lipid center
                    int closest_y = (int)ceil(centers_y[j]/cell_size);   //find the lattice point in y that is closest to the lipid center
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
                iv1d final_assignment(events_ref.num_lipids,0);        //holds the final shell assignment for each lipid after noise filtering
                iv1d final_assignment_box(events_ref.num_lipids,0);    //holds the final box assignment for each lipid after noise filtering

                if(this_frame >= 2*range_big + 1) //filter is full
                {
                    for(j=0; j<events_ref.num_lipids; j++) //loop over lipids
                    { 
                        int delta_shell = range_big - range_shell;   //tells difference in size between shell filter and largest filter
                        int delta_box   = range_big - range_box;     //tells difference in size between box filter and largest filter

                        int count_1   = 0;                           //how many timis is the lipid in shell 1 over the window
                        int count_2   = 0;                           //how many timis is the lipid in shell 2 over the window
                        int count_3   = 0;                           //how many timis is the lipid in shell 3 over the window
                        int count_4   = 0;                           //how many timis is the lipid in shell 4 over the window
                        int count_5   = 0;                           //how many timis is the lipid in shell 5 over the window
                        int count_6   = 0;                           //how many timis is the lipid in shell other over the window
                        int count_box = 0;                           //how many timis is the lipid in the box over the window

                        //count frequency of shells and being in box
                        for(k=0; k<2*range_big + 1; k++) //loop over filter
                        {
                            if(age[k][j] > delta_shell && age[k][j] <= (2*range_big + 1) - delta_shell)  //chooses frames to use in the shells filter (whether it is the bigger or smaller filter) 
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
                            
                            if(age[k][j] > delta_box && age[k][j] <= (2*range_big + 1) - delta_box) //chooses frames to use in the box filter (whether it is the bigger or smaller filter) 
                            {
                                if(box_window[k][j] == 1)
                                {
                                    count_box = count_box + 1;
                                }
                            }
                        }

                        final_assignment[j] = j;  //set equal to the lipid number for lipids that are outside the box

                        if((double)count_box/(double)(2*range_box + 1) > 0.5) //lipid is inside the box
                        {
                            final_assignment_box[j] = 1;
                        }

                        if(count_1 >= count_2 && count_1 >= count_3 && count_1 >= count_4 && count_1 >= count_5 && count_1 >= count_6 && final_assignment_box[j] == 1) //inside the box and shell 1
                        {
                            final_assignment[j] = -1;
                        }
                        else if(count_2 >= count_1 && count_2 >= count_3 && count_2 >= count_4 && count_2 >= count_5 && count_2 >= count_6 && final_assignment_box[j] == 1) //inside the box and shell 2
                        {
                            final_assignment[j] = -2;
                        }
                        else if(count_3 >= count_1 && count_3 >= count_2 && count_3 >= count_4 && count_3 >= count_5 && count_3 >= count_6 && final_assignment_box[j] == 1) //inside the box and shell 3
                        {
                            final_assignment[j] = -3;
                        }
                        else if(count_4 >= count_1 && count_4 >= count_2 && count_4 >= count_3 && count_4 >= count_5 && count_4 >= count_6 && final_assignment_box[j] == 1) //inside the box and shell 4
                        {
                            final_assignment[j] = -4;
                        }
                        else if(count_5 >= count_1 && count_5 >= count_2 && count_5 >= count_3 && count_5 >= count_4 && count_5 >= count_6 && final_assignment_box[j] == 1) //inside the box and shell 5
                        {
                            final_assignment[j] = -5;
                        }
                        else if(count_6 >= count_1 && count_6 >= count_2 && count_6 >= count_3 && count_6 >= count_4 && count_6 >= count_5 && final_assignment_box[j] == 1) //inside the box and shell 6 (other)
                        {
                            final_assignment[j] = -6;
                        }
                    }

                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //                                                                                                           //
                    // Check for tansitions between shells and record dwell times                                                //
                    //                                                                                                           //
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    for(j=0; j<5; j++) //loop over shells
	            { 
                        for(k=0; k<events_ref.num_lipids; k++) //loop over lipids
                        {
                            if(final_assignment[k] == -1*(j+1) && (dump == 0 || i < end) ) //belong to shell j
                            {
                                bound[k][j] = bound[k][j] + 1;
                            }
                            else
                            {
                                if(bound[k][j] > 0) //record dwell time
                                {
                                    int dwell_t = bound[k][j]*stride;
                                    int bind_f  = stride*(this_frame - range_big);
                                    int bind_i  = bind_f - dwell_t;

                                    //store binding events info
                                    be_shells[j].lipid_nr.push_back(k);
                                    be_shells[j].res_nr.push_back(res_nr[k]);
                                    be_shells[j].res_name.push_back(res_name[k]);
                                    be_shells[j].bind_i.push_back(bind_i);
                                    be_shells[j].bind_f.push_back(bind_f);
                                    be_shells[j].dwell_t.push_back(dwell_t);

                                    bound[k][j] = 0;
                                }
                            }
                        }
		    }

                    if(test == 1)
                    {
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        //                                                                                                           //
                        // Get the old grid and highlight shells selection                                                           //
                        //                                                                                                           //
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        int prev_pos = (this_frame-range_big)%(2*range_big + 1);

                        events_ref.voro_frame = grid_window[prev_pos]; 
                        events_ref.voro_nan_frame  = nan_window[prev_pos];
 
                        for(j=0; j<events_ref.num_g_x; j++) //loop over x
                        {
                            for(k=0; k<events_ref.num_g_y; k++) //loop over y
                            {
                                if(events_ref.voro_nan_frame[k][j] == 0) //not a NaN 
                                {
                                    events_ref.voro_frame[k][j] = final_assignment[events_ref.voro_frame[k][j]];
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
                                if(prot_mask.grid[k][j][2][0] == 1) //belongs to the protein
                                {
                                    events_ref.voro_nan_frame[j][k] = 1;
                                }
                            }
                        }

                        out_file_name = chop_and_add_tag(base_file_name_o, "_" + to_string(stride*(this_frame - range_big)) + ".dat");
                        events_ref.write_voro_frame(out_file_name);
                    }
                }
            }

            MPI_Barrier(MPI_COMM_WORLD);  

            //report time stats
            int current_step = i - begin + 1;
            int my_steps     = end - begin + 1;
            ot_time_stats(t,&counter,current_step,my_steps,world_rank,"frames");
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Write binding events to file                                                                              //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(world_rank == 0)
    {
        for(i=0; i<5; i++) //loop over the shells
        {
            be_shells[i].write_binding_events_bin(shell_file_names[i]);
        }
    }

    //log time spent in main loop
    perf.log_time((clock() - t)/CLOCKS_PER_SEC,"Assign Shells");

    //reset the clock
    t = clock();

    //print the performance stats
    perf.print_stats();

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


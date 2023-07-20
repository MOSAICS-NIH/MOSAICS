#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h> 
#include <mpi.h>
#include <fstream>
#include <vector>
using namespace std;

#include "xdr/include/xdrfile_xtc.h"                         //used to read xtc files 
#include "xdr/include/xdr_seek.h"                            //used to get and set the file position in xtc and trr files
#include "xdr/include/xdrfile_trr.h"                         //used to read trr files
#include "xdr/include/xdrfile.h"                             //used to read xtc and trr files
#include "xdr/include/trr_header.h"                          //used to read the header info of trr files
#include "headers/multi_dim_vec.h"                           //This defines multidimensional vectors
#include "headers/switch.h"                                  //This defines a switch (on, off)
#include "headers/file_reader.h"                             //This has basic routines for reading text files
#include "headers/vector_mpi.h"                              //This has routines for collecting vector data
#include "headers/mosat_routines.h"                         //This is where most of the functions called in main are located
#include "headers/file_naming.h"                             //This has routines for added tags to an existing file name    
#include "headers/file_naming_mpi.h"                         //This has routines for added tags to an existing file name (mpi)
#include "headers/command_line_args_mpi.h"                   //This has routines for adding command line arguments
#include "MosAT/program_variables/pv_2d_kinetics.h"         //This has the variables specific to the analysis program
#include "headers/array.h"                                   //This has routines used for working with arrays
#include "headers/performance.h"                             //This has a class for logging performance data
#include "headers/index.h"                                   //This has a class for working with index files
#include "headers/traj.h"                                    //This has a class for working with the trajectory
#include "headers/leaflet_finder.h"                          //This has routines used to find leaflets in membrane simulations
#include "headers/protein_finder.h"                          //This has routines used to find protein atoms
#include "headers/sol_finder.h"                              //This has routines used to find the solvent
#include "headers/grid.h"                                    //This has routines used for working with a grid
#include "headers/protein.h"                                 //This has routines used for working with protein data
#include "headers/fit.h"                                     //This has routines used for fitting data
#include "headers/parallel.h"                                //This has routines for different parallelization schemes
#include "headers/force_serial.h"                            //This has routines used for forcing the code to run on a single mpi process
#include "headers/param.h"                                   //This has routines used for reading complex parameter data
#include "headers/voronoi.h"                                 //This has routines used for computing voronoi diagrams

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function removes any pre-existing binging events files                                               //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double clear_be_files(Trajectory &traj,system_variables &s,program_variables &p)
{
    int i = 0;
    int j = 0;

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("Removing any previous binding events files. \n\n");
    }

    if(traj.my_num_g_x > 0)
    {
        for(i=traj.my_xi; i<=traj.my_xf; i++) //loop over grid x
        {
            for(j=0; j<p.num_g_y; j++) //loop over grid y
            {
                //create file name 
                string tag = "_" + to_string(i) + "_" + to_string(j) + ".be";
                string binding_events_file_name = chop_and_add_tag(p.k_file_name,tag);

                remove(binding_events_file_name.c_str());
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    int term = 0;

    if(s.world_rank == 0)
    {
        for(i=0; i<=p.num_g_x; i++) //loop over grid x
        {
            for(j=0; j<p.num_g_y; j++) //loop over grid y
            {
                //create file name 
                string tag = "_" + to_string(i) + "_" + to_string(j) + ".be";
                string binding_events_file_name = chop_and_add_tag(p.k_file_name,tag);

                FILE *test_file = fopen(binding_events_file_name.c_str(), "r");
                if(test_file == NULL)
                {
                }
                else 
                {
                    printf("Binding events file %s was not deleted. Please delete binding events file manually to prevent duplicate entries. \n",binding_events_file_name.c_str());
                    fclose(test_file);
                    term = 1;
                }
            }
        }
    }

    MPI_Bcast(&term, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(term == 1)
    {
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //compute and return time spent in function
    return (clock() - s.t)/CLOCKS_PER_SEC;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function deposits the lipids to rho_t for excluding data in later steps                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void update_density(Trajectory &traj,system_variables &s,program_variables &p,dv2d &rho_t,Param &param)
{
    double hx1 = 0;                           //mapping atom 1 x-coord
    double hx2 = 0;                           //mapping atom 2 x-coord
    double hy1 = 0;                           //mapping atom 1 y-coord
    double hy2 = 0;                           //mapping atom 2 y-coord
    int i      = 0;                           //standard variable used in loops
    int j      = 0;                           //standard variable used in loops
    int k      = 0;                           //standard variable used in loops

    for(i=0; i<traj.target_leaflet.size(); i++) //loop over target leaflet atoms
    {
        //get the first and last atom of the current lipid
        int min = traj.t_lip_start(i);
        int max = traj.t_lip_end(i);

        //jump to the next lipid
        i = traj.next_target_lipid(i);

        for(j=0; j<param.main_size_y(); j++) //loop over lipid types
        { 
            if(strcmp(traj.res_name[min].c_str(), param.param_main_s[j][0].c_str()) == 0) //lipid type is correct
            {
                for(k=min; k<=max; k++) //loop over current residue atoms
                {
                    if(strcmp(traj.atom_name[k].c_str(), param.param_main_s[j][1].c_str()) == 0) //map atom 1
                    {
                        hx1 = traj.r[k][0];
                        hy1 = traj.r[k][1];

                        //count the lipid for rho_t
                        add_to_grid_d(hx1,hy1,p.radius,p.cell_size,rho_t,1.0);
                    }
                    if(strcmp(traj.atom_name[k].c_str(), param.param_main_s[j][2].c_str()) == 0) //map atom 2
                    {
                        hx2 = traj.r[k][0];
                        hy2 = traj.r[k][1];

                        //count the lipid for rho_t
                        add_to_grid_d(hx2,hy2,p.radius,p.cell_size,rho_t,1.0);
                    }
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints the lipid type for each lipid given its index number                                 //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void print_lipid_info(Trajectory &traj,system_variables &s,program_variables &p,sv1d &res_name_lip,iv1d &res_nr_lip,int num_lipids)
{
    int i = 0;
   
    if(s.world_rank == 0)
    {
        //create file name 
        string tag = "_lipid_info";
        string info_file_name = add_tag(p.k_file_name,tag);

        printf("Writing lipid information to %s. \n\n",info_file_name.c_str());
        fflush(stdin);

        //check if header information has been printed
        FILE *info_file = fopen(info_file_name.c_str(), "w");
        if(info_file == NULL)
        {
            printf("Could not open file for writing (%s). \n",info_file_name.c_str());
        }
        else
        {
            fprintf(info_file," %10s %10s %10s \n","lipid","res_name","res_id");
            for(i=0; i<num_lipids; i++)
            {
                fprintf(info_file," %10d %10s %10d \n",i,res_name_lip[i].c_str(),res_nr_lip[i]);
            }
            fclose(info_file);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function stores information about the lipids like the residue name and number corresponding to a     //
// lipid number                                                                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void store_lipid_info(Trajectory &traj,system_variables &s,program_variables &p,Param &param,iv1d &res_nr_lip,sv1d &res_name_lip)
{
    int i         = 0;                                  //Standard variable used in loops
    int j         = 0;                                  //Standard variable used in loops
    int lip_count = -1;

    for(i=0; i<traj.target_leaflet.size(); i++) //loop over target leaflet atoms
    {
        //get the first and last atom of the current lipid
        int min = traj.t_lip_start(i);
        int max = traj.t_lip_end(i);

        //jump to the next lipid
        i = traj.next_target_lipid(i);

        //index target lipids and store res_id and res_name
        for(j=0; j<param.main_size_y(); j++) //loop over lipid types
        {
            if(strcmp(traj.res_name[min].c_str(), param.param_main_s[j][0].c_str()) == 0) //lipid type is correct
            {
                lip_count++;
                res_nr_lip[lip_count]   = traj.res_nr[min];
                res_name_lip[lip_count] = traj.res_name[min];
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes a voronoi diagram for the current trajectory frame                                 //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_voronoi(Trajectory &traj,system_variables &s,program_variables &p,Param &param,Grid_i &nan)
{
    int dynamic      = 0;                      //use a dynamic box for voronoi diagrams?
    int pbc          = 0;                      //account for periodic boundary conditions in voronoi diagrams?

    //get the voronoi diagram
    Grid_i voronoi;

    if(p.com == 1)
    {
        voronoi = voronoi_diagram_com(traj,p.APS,p.num_g_x,p.num_g_y,param,p.c_dist,p.voro_stamp_rad,p.v_prot,dynamic,pbc);
    }
    else 
    {
        voronoi = voronoi_diagram(traj,p.APS,p.num_g_x,p.num_g_y,param,p.c_dist,p.voro_stamp_rad,p.v_prot,dynamic,pbc);
    }

    //write the voronoi diagram to file
    string voronoi_file_name = add_tag(p.k_file_name,to_string(traj.get_frame_global()));
    voronoi.set_output(voronoi_file_name,p.out_data_format);
    voronoi.write_grid(nan);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function analyzes a voronoi diagram and stores binding events information                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void analyze_diagram(Trajectory &traj,system_variables &s,program_variables &p,
                   iv3d &bound,iv3d &filter,iv2d &voro,iv3d &bi,iv3d &bf,iv3d &be_lipid_nr,iv1d &res_nr_lip)
{
    int i           = 0;                                  //Standard variable used in loops
    int j           = 0;                                  //Standard variable used in loops
    int k           = 0;                                  //Standard variable used in loops
    int l           = 0;                                  //Standard variable used in loops
    int m           = 0;                                  //Standard variable used in loops
    int n           = 0;                                  //Standard variable used in loops
    int ef_x        = 0;                                  //Effective value of of x used for adding to arrays
    int pos         = traj.current_frame%(2*p.range + 1); //Tells current position in filter

    if(traj.my_num_g_x > 0) //If the user specifies too many cores some may have no grid points assigned.
    {
        for(i=traj.my_xi; i<=traj.my_xf; i++) //loop over grid x
        {
            //set effective x used for filling arrays
            ef_x = i - traj.my_xi;

            for(j=0; j<p.num_g_y; j++) //loop over grid y
            {
                //add the lipid number to the window (noise fileter)
                if(voro[j][i] == -1) //protein 
                {
                    filter[ef_x][j][pos] = -1;
                }
                else
                {
                    for(k=0; k<p.num_lipids; k++) //loop over lipids
                    {
                        if(res_nr_lip[k]  == voro[j][i])
                        {
                            filter[ef_x][j][pos] = k;
                        }
                    }
                }

                if(traj.current_frame >= (2*p.range)) //filter is full, compute frequencies
                {
                    for(k=0; k<p.num_lipids; k++) //loop over lipids
                    {
                        double freq = 0;
                        for(m=0; m<(2*p.range + 1); m++) //loop over flanking frames
                        {
                            if(filter[ef_x][j][m] == k) //lipid is a match
                            {
                                freq = freq + 1.0;
                            }
                        }
                        freq = freq/((double)(2*p.range + 1));

                        if(freq >= 0.5 && (p.dump == 0 || traj.current_frame < traj.get_ef_frames()-1) ) //lipid takes the grid point or last frame dump
                        {
                            bound[ef_x][j][k] = bound[ef_x][j][k] + 1;
                        }
                        else //lipid is not bound
                        {
                            if(bound[ef_x][j][k] > 0) //lipid was bound. add dwell time
                            {
                               int time   = bound[ef_x][j][k];
                               int bind_i = traj.current_frame - time;      //first frame binding was noted 
                               int bind_f = traj.current_frame - 1;         //last time binding was noted

                               //store binding event in memory
                               bi[j][ef_x].push_back(bind_i);
                               bf[j][ef_x].push_back(bind_f);
                               be_lipid_nr[j][ef_x].push_back(k);
                            }
                            bound[ef_x][j][k] = 0;
                        }
                    }
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads in voronoi digrams, filters them and records binding events data                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double get_dwell_times(Trajectory &traj,system_variables &s,program_variables &p,iv3d &bound,iv3d &filter,iv2d &voro,
                     iv3d &bi,iv3d &bf,iv3d &be_lipid_nr,iv1d &res_nr_lip,Param &param,sv1d &res_name_lip)
{
    char my_string[200];                     //Used for reading data from single frame files
    int j = 0;
    int k = 0;

    //wait until all cores have computed their voronoi diagrams
    MPI_Barrier(MPI_COMM_WORLD);

    //get the residue number and names
    store_lipid_info(traj,s,p,param,res_nr_lip,res_name_lip);

    //write lipid infor to file
    print_lipid_info(traj,s,p,res_name_lip,res_nr_lip,p.num_lipids);

    if(s.world_rank == 0)
    {
        printf("Computing residence times.\n");
        fflush(stdin);
    }

    //reset counter for printing progress + time estimates
    s.counter = 0;

    //record time spent analyzing diagrams
    s.t = clock();

    //read in voronoi diagrams and compute residence times
    for(traj.current_frame=0; traj.current_frame<traj.get_ef_frames(); traj.current_frame++) //loop over traj frames
    {
        string voronoi_file_name = add_tag(p.k_file_name,to_string(traj.current_frame));

        FILE *voro_file = fopen(voronoi_file_name.c_str(), "r");
        if(voro_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",voronoi_file_name.c_str());
        }
        else
        {
            for(k=0; k<p.num_g_y; k++) //loop over y-dimension
            {
                for(j=0; j<p.num_g_x; j++) //loop over x-dimension
                {
                    int result = fscanf(voro_file, "%s,", my_string);
                    voro[k][j] = atof(my_string);
                }
            }
            fclose(voro_file);
        }

        //compute residence times
        analyze_diagram(traj,s,p,bound,filter,voro,bi,bf,be_lipid_nr,res_nr_lip);

        //report progress
        time_stats(s.t,&s.counter,traj.current_frame,traj.get_ef_frames(),s.world_rank);
        fflush(stdin);
    }
    //compute and return time spent in function
    return (clock() - s.t)/CLOCKS_PER_SEC; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects and writes rho_t and writes binding events to output files                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,iv3d &bi,iv3d &bf,iv3d &be_lipid_nr,
                         iv1d &res_nr_lip,sv1d &res_name_lip,dv2d &rho_t,Grid_i &nan)
{
    int i        = 0;                              //general variable used in loops
    int j        = 0;                              //general variable used in loops
    int k        = 0;                              //general variable used in loops

    //record time when beginning analysis
    s.t = clock();

    //wait until all cores have computed their voronoi diagrams
    MPI_Barrier(MPI_COMM_WORLD);

    //collect rho_t
    collect_grid_d(s.world_size,s.world_rank,rho_t);

    if(s.world_rank == 0)
    {
       //create name of output files additional to koff
       string rho_file_name = add_tag(p.k_file_name,"_rho");

       printf("\nWriting rho to %s. \n\n",rho_file_name.c_str());
       fflush(stdin);

       //exclude insignificant data and write grid data to file
       double avg_rho_t = get_average_rho(rho_t);
       exclude_insignificant_data(0.0,avg_rho_t,rho_t,nan.grid);
       write_grid_to_file_d(p.out_data_format,p.cell_size,rho_file_name,nan.grid,rho_t);
    }

    if(s.world_rank == 0)
    {
        printf("Writing binding events files.\n");
        fflush(stdin);
    }

    //write binding events files
    for(i=traj.my_xi; i<=traj.my_xf; i++) //loop over grid x
    {
        //set effective x used for filling arrays
        int ef_x = i - traj.my_xi;

        for(j=0; j<p.num_g_y; j++) //loop over grid y
        {
            if(be_lipid_nr[j][ef_x].size() > 0) //only create files if binding events were sampled
            {
                //create file name 
                string tag = "_" + to_string(i) + "_" + to_string(j) + ".be";
                string binding_events_file_name = chop_and_add_tag(p.k_file_name,tag);

                //open file for writing
                FILE *binding_events_file = fopen(binding_events_file_name.c_str(), "w");
                if(binding_events_file == NULL)
                {
                    printf("failure opening %s for writing. \n",binding_events_file_name.c_str());
                    fflush(stdin);
                }
                else
                {
                    //print header 
                    fprintf(binding_events_file," x_i %10d y_i %10d ef_dt(ps) %10f ef_frames %d num_lipids %10d num_g_x %10d num_g_y %10d APS(nm^2) %10f \n\n",i,j,p.delta_t,traj.get_ef_frames(),p.num_lipids,p.num_g_x,p.num_g_y,p.APS);
                    fprintf(binding_events_file," %10s %10s %10s %15s %15s %20s \n","lipid","res_nr","res_name","bind_i(frame)","bind_f(frame)","dwell time(frames)");
                    fprintf(binding_events_file," %10s-%10s-%10s-%15s-%15s-%20s \n","----------","----------","----------","---------------","---------------","--------------------");
    
                    //print binding events
                    for(k=0; k<bi[j][ef_x].size(); k++) //loop over binding events
                    {
                        int res_nr      = res_nr_lip[be_lipid_nr[j][ef_x][k]];
                        string res_name = res_name_lip[be_lipid_nr[j][ef_x][k]];
                        int time        = bf[j][ef_x][k] - bi[j][ef_x][k] + 1; 

                        fprintf(binding_events_file," %10d %10d %10s %15d %15d %20d \n",be_lipid_nr[j][ef_x][k],res_nr,res_name.c_str(),bi[j][ef_x][k],bf[j][ef_x][k],time);
                    }
                }
                //close file
                fclose(binding_events_file);
            }
        }
    }

    //remove voronoi diagram files
    if(p.b_clean == 1)
    {
        if(s.world_rank == 0)
        {
            printf("\nRemoving voronoi diagram files. \n");
            fflush(stdin);

            for(traj.current_frame=0; traj.current_frame<traj.get_ef_frames(); traj.current_frame++) //loop over traj frames
            {
                string voronoi_file_name = add_tag(p.k_file_name,to_string(traj.get_frame_global()));

                remove(voronoi_file_name.c_str());

                FILE *test_file = fopen(voronoi_file_name.c_str(), "r");
                if(test_file == NULL)
                {
                }
                else
                {
                    printf("\nUnable to delete file %s.\n",voronoi_file_name.c_str());
                    fflush(stdin);
                    fclose(test_file);
                }
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //compute and return time spent in function
    return (clock() - s.t)/CLOCKS_PER_SEC;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This is the main function which executes the other functions.                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, const char * argv[]) 
{
    //Here we set up the MPI environment
    MPI_Init(NULL, NULL);
    
    //nowe we define the system and program variables and initialize them
    system_variables s;
    program_variables p;
    initialize_system_variables(&s);
    initialize_program_variables(&p);

    //create object for logging performance data
    Performance perf; 

    //set the name of your analysis program here
    s.program_name = "2D Kinetics";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                                   s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                                   s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                                    s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                                               s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                                s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                                  s.world_rank, s.cl_tags, &p.b_lsq,     0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                                    s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",                    s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-crd",    p.param_file_name,            "Selection card with target lipids for voronoi diagrame + mapping atoms (crd)", s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-k",      p.k_file_name,                "Output filename used to derive filenames for kinetics data (.dat)",            s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,           "PDB file with sorted leaflets (pdb)",                                          s.world_rank, s.cl_tags, &p.b_lf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,  "File with additional leaflet finder parameters (prm)",                         s.world_rank, s.cl_tags, &p.b_lf_param,0);
    add_argument_mpi_s(argc,argv,"-pf_pdb", p.pf_pdb_file_name,           "PDB file with protein marked (pdb)",                                           s.world_rank, s.cl_tags, &p.b_pf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-pf_prm", p.protein_finder_param_name,  "File with additional protein finder parameters (prm)",                         s.world_rank, s.cl_tags, &p.b_pf_param,0);
    add_argument_mpi_d(argc,argv,"-APS",    &p.APS,                       "Area per grid box (nm^2)",                                                     s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-r",      &p.radius,                    "Stamping radius when computing the sample count (nm)",                         s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-bx",     &p.box_x,                     "Grid x dimension (nm)",                                                        s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-by",     &p.box_y,                     "Grid y dimension (nm)",                                                        s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-leaf",   &p.leaflet,                   "Which leaflet? (0:both 1:upper 2:lower)",                                      s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-odf",    &p.out_data_format,           "What is the output data format? (0:matrix 1:vector)",                          s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-range",  &p.range,                     "Noise filter half width (frames)",                                             s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-dt",     &p.delta_t,                   "Effective time step between trajectory frames analyzed (ps) ",                 s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-voro_r", &p.voro_stamp_rad,            "Stamping radius used when computing Voronoi diagrams (nm)",                    s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-v_prot", &p.v_prot,                    "Include protein atoms in voronoi diagram? (0:no 1:yes)",                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-c_dist", &p.c_dist,                    "Distance cutoff for counting protein atoms in voronoi diagram (nm)",           s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-clean",  &p.b_clean,                   "Remove voronoi diagrams after writing binding events files? (0:no 1:yes)",     s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-dump",   &p.dump,                      "Dump bound lipids on last frame? (0:no 1:yes)",                                s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-com",    &p.com,                       "Use the lipid center of mass for voronoi tessellations? (0:no 1:yes)",         s.world_rank, s.cl_tags, nullptr,      0);
    conclude_input_arguments_mpi(argc,argv,s.world_rank,s.program_name,s.cl_tags);

    //create a trajectory
    Trajectory traj; 

    //set trajectory parameters
    traj.set_block_parallel(on);
    traj.set_traj(p.in_file_name);
    traj.set_ref(p.ref_file_name);
    traj.set_lsq(p.lsq_index_file_name,p.b_lsq,p.lsq_dim,p.lsq_ref);
    traj.set_res(p.stride,p.start_frame,p.end_frame);

    //analyze the trajectory (log time spent) 
    perf.log_time(traj.build(),"Analyze Trajectory");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //check file extensions                                                                                     
    check_extension_mpi(s.world_rank,"-crd",p.param_file_name,".crd");
    check_extension_mpi(s.world_rank,"-k",p.k_file_name,".dat");

    if(p.b_lf_pdb == 1)
    {
        check_extension_mpi(s.world_rank,"-lf_pdb",p.lf_pdb_file_name,".pdb");
    }
    if(p.b_lf_param == 1)
    {
        check_extension_mpi(s.world_rank,"-lf_prm",p.leaflet_finder_param_name,".prm");
    }
    if(p.b_pf_pdb == 1)
    {
        check_extension_mpi(s.world_rank,"-pf_pdb",p.pf_pdb_file_name,".pdb");
    }
    if(p.b_pf_param == 1)
    {
        check_extension_mpi(s.world_rank,"-pf_prm",p.protein_finder_param_name,".prm");
    }

    //create parameter file
    Param param;

    //read parameter file
    param.get_param(p.param_file_name,4,3,1);

    //check the integrity of the parameter files
    if(param.check_file() == 0) //bad files. kill program 
    {
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }

    //run leaflet/protein finder
    traj.get_leaflets(p.leaflet,p.leaflet_finder_param_name,p.b_lf_param);
    traj.get_protein(p.protein_finder_param_name,p.b_pf_param);

    //print a pdb with distinguished leaflets/protein
    traj.write_leaflets(p.lf_pdb_file_name,p.b_lf_pdb);
    traj.write_protein(p.pf_pdb_file_name,p.b_pf_pdb);

    //print info about the target leaflet
    traj.get_leaflet_stats();

    //print info about the lipid selection
    traj.get_lipid_selection_stats(param.get_column_s(0),"-crd");

    //print info about the protein
    traj.get_prot_stats();

    //count the number of lipids in the target leaflet(s) of the target type
    p.num_lipids = traj.count_target_lipids_type(param.get_column_s(0));

    //set the lenght of a cell for the grid
    p.cell_size = sqrt(p.APS);

    //counts grid points in x and y dimensions
    get_grid_size(p.box_x,p.box_y,traj.ibox,&p.num_g_x,&p.num_g_y,p.cell_size);

    //print information about the grid
    print_grid_stats(p.box_x,p.box_y,traj.ibox,p.num_g_x,p.num_g_y,p.cell_size,s.world_rank);

    //set the parallelization scheme to grid points
    traj.parallelize_by_grid(p.num_g_x);

    //compute memory requirements  
    double memory_filter    = ((double)traj.my_num_g_x*(double)p.num_g_y*(double)((2.0*(double)p.range + 1.0)*4.0))/1000000.0;
    double memory_bound     = ((double)traj.my_num_g_x*(double)p.num_g_y*(double)p.num_lipids*4.0)/1000000.0;

    if(s.world_rank == 0)
    {
        printf("Memory Estimates per core: bound[] %f MB     filter %f MB \n\n",memory_bound,memory_filter);
    }

    //create multi-dimensional vectors to store data 
    iv3d filter     (traj.my_num_g_x, iv2d(p.num_g_y, iv1d(2*p.range + 1,0)));        //contains a history of which lipid was closest to each grid point
    iv3d bound      (traj.my_num_g_x, iv2d(p.num_g_y, iv1d(p.num_lipids,0)));         //keeps track of how many frames a lipid has been bound for 

    //allocate other memory
    dv2d rho_t(p.num_g_y,dv1d(p.num_g_x,0.0));                                        //stores the total lipid density

    //store the binding events for each lattice point
    iv3d bi(p.num_g_y,iv2d(p.num_g_x,iv1d(0,0)));
    iv3d bf(p.num_g_y,iv2d(p.num_g_x,iv1d(0,0)));
    iv3d be_lipid_nr(p.num_g_y,iv2d(p.num_g_x,iv1d(0,0)));

    //create vectors to hold the residue number and name for a given lipid number 
    iv1d res_nr_lip(p.num_lipids,0);
    sv1d res_name_lip(p.num_lipids);

    //create grid for reading in voronoi diagrams
    iv2d voro(p.num_g_y,iv1d(p.num_g_x,0));

    //create a grid to use when writing voronoi diagrams to file
    Grid_i nan;

    //get the grid dimensions
    nan.set_dim(p.APS,p.num_g_x,p.num_g_y);

    //erase existing binding events files
    perf.log_time(clear_be_files(traj,s,p),"Cleaning");
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //print info about the worlk load distribution
    traj.workload();

    //print that analysis is beginning
    traj.report_progress();

    s.t = clock();
    //read read frames of the trajector and perform analysis
    for(traj.current_frame=0; traj.current_frame<traj.get_num_frames(); traj.current_frame++)
    {
        traj.read_traj_frame();

        traj.do_fit();

        get_voronoi(traj,s,p,param,nan);

        update_density(traj,s,p,rho_t,param);

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
        fflush(stdin);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //print info about the worlk load distribution
    if(s.world_rank == 0)
    {
        printf("\n");
    }
    traj.workload_grid();

    //read voronoi diagrams and compute residence times
    perf.log_time(get_dwell_times(traj,s,p,bound,filter,voro,bi,bf,be_lipid_nr,res_nr_lip,param,res_name_lip),"Dwell time");

    //write binding events files
    perf.log_time(finalize_analysis(traj,s,p,bi,bf,be_lipid_nr,res_nr_lip,res_name_lip,rho_t,nan),"Write BE");

    //print the performance stats
    perf.print_stats();

    //print closing statements
    print_closing(s.world_rank);

    //relinquish the mpi environment
    MPI_Finalize();

    return 0;
}

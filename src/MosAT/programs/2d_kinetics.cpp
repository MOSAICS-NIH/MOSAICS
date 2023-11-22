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
#include "headers/grid_lt.h"                                 //This has routines used for reading in grid data
#include "headers/voronoi.h"                                 //This has routines used for computing voronoi diagrams
#include "headers/binding_events.h"                          //This has routines used for reading in binding events files

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function tags frames whose tessellation data is to be stored in memory                               //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
iv2d tag_frames_mem(Trajectory &traj,system_variables &s,program_variables &p,int num_frames_ram)
{
    int i       = 0;                           //standard variable used in loops
    int counter = 0;                           //count the target frames as they are encountered

    iv2d mem_tags(traj.get_num_frames(),iv1d(2,0));     //tags which frames should be stored in memory

    for(traj.current_frame=0; traj.current_frame<traj.get_num_frames(); traj.current_frame++) //loop over core's frames
    {
        if(traj.current_frame > p.range && traj.current_frame <= p.range + num_frames_ram && traj.current_frame < traj.get_num_frames() - p.range - 1) //store in memory
        {
            //store tags	
            mem_tags[traj.current_frame][0] = 1;
	    mem_tags[traj.current_frame][1] = counter;

	    //count the frame
	    counter++;
        }
    }

    return mem_tags;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function tags frames whose tessellation data is to be stored in memory (global)                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
iv2d tag_frames_mem_global(Trajectory &traj,system_variables &s,program_variables &p,int num_frames_ram)
{
    int i       = 0;                           //standard variable used in loops
    int counter = 0;                           //count the target frames as they are encountered

    iv2d mem_tags_global(traj.get_ef_frames(),iv1d(2,0));     //tags which frames the current core should stored in memory (global)

    for(traj.current_frame=0; traj.current_frame<traj.get_num_frames(); traj.current_frame++) //loop over core's frames
    {
        if(traj.current_frame > p.range && traj.current_frame <= p.range + num_frames_ram && traj.current_frame < traj.get_num_frames() - p.range - 1) //store in memory
        {
            //store tags        
            mem_tags_global[traj.get_frame_global()][0] = 1; 
            mem_tags_global[traj.get_frame_global()][1] = counter;

            //count the frame
            counter++;
        }
    }

    return mem_tags_global;
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

        printf("\nWriting lipid information to %s. \n\n",info_file_name.c_str());
        fflush(stdin);

        //check if header information has been printed
        FILE *info_file = fopen(info_file_name.c_str(), "w");
        if(info_file == NULL)
        {
            printf("Could not open file for writing (%s). \n",info_file_name.c_str());
        }
        else
        {
            fprintf(info_file," %10s %10s %10s \n","#lipid","#res_name","#res_id");
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
void store_lipid_info(Trajectory &traj,system_variables &s,program_variables &p,Param &param,iv1d &res_nr_lip,
		      sv1d &res_name_lip,iv1d &lip_nr_res)
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
                res_nr_lip[lip_count]          = traj.res_nr[min];
                res_name_lip[lip_count]        = traj.res_name[min];
                lip_nr_res[traj.res_nr[min]-1] = lip_count;
	    }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes a voronoi diagram for the current trajectory frame                                 //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_voronoi(Trajectory &traj,system_variables &s,program_variables &p,Param &param,Grid_i &nan,Grid_lt &prot_mask,
		 iv2d &mem_tags,iv3d &voro_ram)
{
    int dynamic      = 0;                      //use a dynamic box for voronoi diagrams?
    int pbc          = 0;                      //account for periodic boundary conditions in voronoi diagrams?
    int i = 0;
    int j = 0;

    //get the voronoi diagram
    Grid_i voronoi;

    if(p.com == 1)
    {
        voronoi = voronoi_diagram_com(traj,p.APS,p.num_g_x,p.num_g_y,param,p.c_dist,p.voro_stamp_rad,p.v_prot,dynamic,pbc,prot_mask,p.b_prot_mask);
    }
    else 
    {
        voronoi = voronoi_diagram(traj,p.APS,p.num_g_x,p.num_g_y,param,p.c_dist,p.voro_stamp_rad,p.v_prot,dynamic,pbc,prot_mask,p.b_prot_mask);
    }

    if(mem_tags[traj.current_frame][0] == 0) //write tessellation data to hdd
    {
        string voronoi_file_name = add_tag(p.k_file_name,to_string(traj.get_frame_global()));
        ofstream voro_file(voronoi_file_name, ios::out | ios::binary);
        voro_file.write(reinterpret_cast<const char *>(&p.num_g_x), sizeof(int));
        voro_file.write(reinterpret_cast<const char *>(&p.num_g_y), sizeof(int));
        for(i=0; i<p.num_g_y; i++)
        {
            for(j=0; j<p.num_g_x; j++)
            {
                voro_file.write(reinterpret_cast<const char *>(&voronoi.grid[i][j]), sizeof(int));
            }
        }
        voro_file.close();

        //check file afterwards
        if(!voro_file)
        {
            printf("Failure writing voronoi file %s. \n",voronoi_file_name.c_str());
        }
    }
    else //store tessellation data in memory 
    {
        for(i=0; i<p.num_g_y; i++)
        {
            for(j=0; j<p.num_g_x; j++)
            {
                voro_ram[mem_tags[traj.current_frame][1]][j][i] = voronoi.grid[i][j];
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function analyzes a voronoi diagram and stores binding events information                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void analyze_diagram(Trajectory &traj,system_variables &s,program_variables &p,iv3d &bound,iv3d &filter,iv2d &voro,
		     iv3d &bi,iv3d &bf,iv3d &be_lipid_nr,iv1d &res_nr_lip,iv1d &lip_nr_res,int ef_frame,int ef_num_frames,
		     int real_frame)
{
    int i           = 0;                                  //Standard variable used in loops
    int j           = 0;                                  //Standard variable used in loops
    int k           = 0;                                  //Standard variable used in loops
    int l           = 0;                                  //Standard variable used in loops
    int m           = 0;                                  //Standard variable used in loops
    int n           = 0;                                  //Standard variable used in loops
    int pos         = ef_frame%(2*p.range + 1); //Tells current position in filter

    for(i=0; i<p.num_g_x; i++) //loop over grid x
    {
        for(j=0; j<p.num_g_y; j++) //loop over grid y
        {
            //add the lipid number to the window (noise fileter)
            if(voro[j][i] == -1) //protein 
            {
                filter[j][i][pos] = -1;
            }
            else
            {
                filter[j][i][pos] = lip_nr_res[voro[j][i]-1];
            }

            if(ef_frame >= (2*p.range)) //filter is full, compute frequencies
            {
                dv1d lipid_freq(p.num_lipids, 0.0);  //hold the frequency of each lipid

                for(k=0; k<(2*p.range + 1); k++) //loop over flanking frames
                {
                    int this_lipid = filter[j][i][k]; 

                    if(this_lipid > -1) //not occupied by the protein
                    {
                        lipid_freq[this_lipid] = lipid_freq[this_lipid] + 1.0;
                    }
                }

                for(k=0; k<p.num_lipids; k++) //loop over lipids
                {   
                    double freq = lipid_freq[k]/((double)(2*p.range + 1));

                    if(freq >= 0.5 && (p.dump == 0 || ef_frame < ef_num_frames-1) ) //lipid takes the grid point or last frame dump
                    {
                        bound[j][i][k] = bound[j][i][k] + 1;
                    }
                    else //lipid is not bound
                    {
                        if(bound[j][i][k] > 0) //lipid was bound. add dwell time
                        {
                           int time   = bound[j][i][k];
                           int bind_i = real_frame - time;      //first frame binding was noted 
                           int bind_f = real_frame - 1;         //last time binding was noted

                           //note: b_f will leave the last frame unbound if dump is used. Careful splicing! 
                           //note: the bi and bf are shifted here to the right by p.range. shouldnt matter though.

                           //store binding event in memory
                           bi[j][i].push_back(bind_i);
                           bf[j][i].push_back(bind_f);
                           be_lipid_nr[j][i].push_back(k);
                        }
                        bound[j][i][k] = 0;
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
                       iv3d &bi,iv3d &bf,iv3d &be_lipid_nr,iv1d &res_nr_lip,Param &param,sv1d &res_name_lip,iv1d &lip_nr_res,
		       iv2d &mem_tags_global,iv3d &voro_ram)
{
    char my_string[200];                     //Used for reading data from single frame files
    int i = 0;
    int j = 0;
    int k = 0;

    //wait until all cores have computed their voronoi diagrams
    MPI_Barrier(MPI_COMM_WORLD);

    //record time spent analyzing diagrams
    s.t = clock();

    //get the residue number and names
    store_lipid_info(traj,s,p,param,res_nr_lip,res_name_lip,lip_nr_res);

    //write lipid infor to file
    print_lipid_info(traj,s,p,res_name_lip,res_nr_lip,p.num_lipids);

    if(s.world_rank == 0)
    {
        printf("Computing residence times.\n");
	printf("-----------------------------------------------------------------------------------------------------------------------------------\n");
        fflush(stdin);
    }

    //reset counter for printing progress + time estimates
    s.counter = 0;

    int start = traj.get_global_frame_i() - p.range - 1; 
    int end   = traj.get_global_frame_f() + p.range + 1;

    if(start < 0)
    {
        start = 0;
    } 
    if(end >= traj.get_ef_frames()) 
    {
        end = traj.get_ef_frames()-1;
    } 

    //read in voronoi diagrams and compute residence times
    for(i=start; i<=end; i++) //loop over traj frames 
    {
        int ef_frame      = i - start;         //set current frame to the effective frame
        int ef_num_frames = end - start + 1;   //how many frames the core is responsible for

        if(mem_tags_global[i][0] == 0) //read from hdd
        { 
            string voronoi_file_name = add_tag(p.k_file_name,to_string(i));

            ifstream voro_file(voronoi_file_name, ios::out | ios::binary);
            if(!voro_file) 
            {
                printf("Could not find voronoi file %s. \n",voronoi_file_name.c_str());
            }
            else 
            {
                int size_x = 0;
                int size_y = 0;
                voro_file.read(reinterpret_cast<char *>(&size_x), sizeof(int));
                voro_file.read(reinterpret_cast<char *>(&size_y), sizeof(int));
                for(j=0; j<size_y; j++) //loop over y
                {
                    for(k=0; k<size_x; k++) //loop over x
                    {
                        voro_file.read(reinterpret_cast<char *>(&voro[j][k]), sizeof(int));
                    }
                }
                voro_file.close();
            }
        }
        else //pull tessellation from memory
        {
            for(j=0; j<p.num_g_y; j++) //loop over y
            {
                for(k=0; k<p.num_g_x; k++) //loop over x
                {
                    voro[j][k] = voro_ram[mem_tags_global[i][1]][k][j]; 
                }
            }    
        }

        //compute residence times
        analyze_diagram(traj,s,p,bound,filter,voro,bi,bf,be_lipid_nr,res_nr_lip,lip_nr_res,ef_frame,ef_num_frames,i);

        //report progress
	time_stats(s.t,&s.counter,ef_frame,ef_num_frames,s.world_rank);
        fflush(stdin);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //compute and return time spent in function
    return (clock() - s.t)/CLOCKS_PER_SEC; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects and writes rho_t and writes binding events to output files                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,iv3d &bi,iv3d &bf,iv3d &be_lipid_nr,
                         iv1d &res_nr_lip,sv1d &res_name_lip,dv2d &rho_t,Grid_i &nan,Performance &perf,iv2d &mem_tags)
{
    int i              = 0;                              //general variable used in loops
    int j              = 0;                              //general variable used in loops
    int k              = 0;                              //general variable used in loops
    int grid_count     = 0;                              //count grid points as they are encountered 
    int result         = 0;                              //record if the be file was read successfully
    int current_rank   = 0;                              //used to loop over mpi cores when mending BE files

    //record time when beginning analysis
    s.t = clock();

    //wait until all cores have computed their voronoi diagrams
    MPI_Barrier(MPI_COMM_WORLD);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Write out lipid density data                                                                              //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Remove voronoi tessellation files                                                                         //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //remove voronoi diagram files
    if(p.b_clean == 1)
    {
        if(s.world_rank == 0)
        {
            printf("Removing voronoi diagram files. \n");
        }

        for(traj.current_frame=0; traj.current_frame<traj.get_num_frames(); traj.current_frame++) //loop over traj frames
        {
            if(mem_tags[traj.current_frame][0] == 0) //tessellation was not stored in memory 
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

    //wait until all cores are here
    MPI_Barrier(MPI_COMM_WORLD);

    //log time spent writing tmp be files
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Clean Voro Files");

    //reset the clock
    s.t = clock();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Write temporary binding events files (partial traj, full grid)                                            //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(s.world_rank == 0)
    {
        printf("\nWriting temporary binding events files (partial trajectory, full grid).\n");
        printf("-----------------------------------------------------------------------------------------------------------------------------------\n");
    }

    //reset counter for printing progress + time estimates
    s.counter = 0;

    string tag = "_" + to_string(s.world_rank) + ".be";
    string binding_events_file_name = chop_and_add_tag(p.k_file_name,tag);
    ofstream be_file(binding_events_file_name, ios::out | ios::binary);

    i64v3d be_pos(s.world_size, i64v2d(p.num_g_y, i64v1d(p.num_g_x,0))); //store the position of each grid point in the be file 

    //write header info
    int num_ef_frames = traj.get_ef_frames();
    be_file.write(reinterpret_cast<const char *>(&p.delta_t),     sizeof(double)); //ef_dt
    be_file.write(reinterpret_cast<const char *>(&num_ef_frames), sizeof(int));    //ef_frames
    be_file.write(reinterpret_cast<const char *>(&p.num_lipids),  sizeof(int));    //num_lipids
    be_file.write(reinterpret_cast<const char *>(&p.num_g_x),     sizeof(int));    //number of grid points in x
    be_file.write(reinterpret_cast<const char *>(&p.num_g_y),     sizeof(int));    //number of grid points in y
    be_file.write(reinterpret_cast<const char *>(&p.APS),         sizeof(double)); //APS
    for(i=0; i<p.num_g_x; i++)
    {
        for(j=0; j<p.num_g_y; j++)
        {
            int num_events = bi[j][i].size();

            //store the current position in the file
            be_pos[s.world_rank][j][i] = be_file.tellp();

            be_file.write(reinterpret_cast<const char *>(&i),              sizeof(int));    //grid points index in x
            be_file.write(reinterpret_cast<const char *>(&j),              sizeof(int));    //grid points index in y
            be_file.write(reinterpret_cast<const char *>(&num_events),     sizeof(int));    //number of binding events
            for(k=0; k<bi[j][i].size(); k++) //loop over binding events
            {
                int res_nr      = res_nr_lip[be_lipid_nr[j][i][k]];
                string res_name = res_name_lip[be_lipid_nr[j][i][k]];
                int time        = bf[j][i][k] - bi[j][i][k] + 1;

		size_t this_res_name_size = res_name.size(); 

                be_file.write(reinterpret_cast<const char *>(&be_lipid_nr[j][i][k]), sizeof(int));                //lipid number
                be_file.write(reinterpret_cast<const char *>(&res_nr),               sizeof(int));                //residue id
                be_file.write(reinterpret_cast<const char *>(&this_res_name_size),   sizeof(this_res_name_size)); //size of residue name
		be_file.write(reinterpret_cast<const char *>(res_name.c_str()),      this_res_name_size);         //lipid number
                be_file.write(reinterpret_cast<const char *>(&bi[j][i][k]),          sizeof(int));                //bi
                be_file.write(reinterpret_cast<const char *>(&bf[j][i][k]),          sizeof(int));                //bf
                be_file.write(reinterpret_cast<const char *>(&time),                 sizeof(int));                //time
            }

            //report progress
            int ef_g     = i*p.num_g_y + j;
            int ef_num_g = p.num_g_x*p.num_g_y;
            time_stats(s.t,&s.counter,ef_g,ef_num_g,s.world_rank);
        }
    }
    be_file.close();

    MPI_Barrier(MPI_COMM_WORLD);

    //collect the position of each lattice point from other cores
    collect_lv3d(s.world_size,s.world_rank,be_pos);

    //broadcast the complete set of positions
    broadcast_lv3d(s.world_size,s.world_rank,be_pos);

    //log time spent writing tmp be files
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Write Tmp BE Files");

    //reset the clock
    s.t = clock();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Write temporary binding events files (full traj, partial grid)                                            //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //print info about the worlk load distribution
    if(s.world_rank == 0)
    {
        printf("\n");
    }
    traj.workload_grid_alt();

    if(s.world_rank == 0)
    {
        printf("Mending binding events data and writing temporary BE files (full trajectory, partial grid). \n");
        printf("-----------------------------------------------------------------------------------------------------------------------------------\n");
    }

    //record time spent analyzing diagrams
    s.t = clock();

    //reset counter for printing progress + time estimates
    s.counter = 0;

    //open file for writing out temporary BE file
    string this_tag          = "_" + to_string(s.world_rank) + "_tmp.be";
    string this_be_file_name = chop_and_add_tag(p.k_file_name,this_tag);
    ofstream tmp_be_file_o(this_be_file_name, ios::out | ios::binary);
    int64_t current_pos = tmp_be_file_o.tellp(); 

    //read in individual binding events files
    for(i=0; i<p.num_g_x; i++) //loop over x
    {
        for(j=0; j<p.num_g_y; j++) //loop over y 
        {
            if(grid_count >= traj.my_gi && grid_count <= traj.my_gf)
            {
                Binding_events be;

                for(k=0; k<s.world_size; k++) //loop over mpi cores
                {
                    //set the name of the current BE file to be read
                    this_tag          = "_" + to_string(k) + ".be";
                    this_be_file_name = chop_and_add_tag(p.k_file_name,this_tag);

                    //read in binding events
                    if(k==0)
                    {
                        result = be.get_binding_events_grid(this_be_file_name,be_pos[k],i,j);
                    }
                    else 
                    {
                        result = be.add_binding_events_grid(this_be_file_name,be_pos[k],i,j);
                    }
                }

                be.get_binding_timeline();
                be.binding_events_from_timeline();

                current_pos = be.write_binding_events_tmp(tmp_be_file_o,current_pos);

                //report progress
                int ef_g     = grid_count - traj.my_gi;
                int ef_num_g = traj.my_gf - traj.my_gi + 1; 
                time_stats(s.t,&s.counter,ef_g,ef_num_g,s.world_rank);
	    }
            grid_count++;
        }
    }
    //close file
    tmp_be_file_o.close();

    //reset grid_count
    grid_count = 0;

    MPI_Barrier(MPI_COMM_WORLD);

    //remove tmp files
    for(i=0; i<s.world_size; i++)
    {
        this_tag          = "_" + to_string(i) + ".be";
        this_be_file_name = chop_and_add_tag(p.k_file_name,this_tag);
        remove(this_be_file_name.c_str());
    }

    //log time spent writing tmp be files
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Mend BE Files");

    //reset the clock
    s.t = clock();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Merge binding events files                                                                                //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(s.world_rank == 0)
    {
        ifstream tmp_be_file_i;                          //file for reading in binding events data
        ofstream be_file_o;                              //file for writing out binding events data

        i64v2d grid_be_pos(p.num_g_x, i64v1d(p.num_g_y,0));   //store the position of each grid point in the be file 

        printf("\nMerging binding events data and writing the final BE file (full trajectory, full grid). \n");
        printf("-----------------------------------------------------------------------------------------------------------------------------------\n");

	//reset counter for printing progress + time estimates
        s.counter = 0;

        //keep track of count restarting for each rank
        int my_grid_count = 0;

        //open file for writing the final binding events data
        this_tag          = ".be";
        this_be_file_name = chop_and_add_tag(p.k_file_name,this_tag);
        be_file_o.open(this_be_file_name, ios::out | ios::binary);
        int64_t current_pos_o = be_file_o.tellp();

        //open initial temporary BE file for reading
        this_tag          = "_" + to_string(0) + "_tmp.be";
        this_be_file_name = chop_and_add_tag(p.k_file_name,this_tag);
        tmp_be_file_i.open(this_be_file_name, ios::out | ios::binary);
        int64_t current_pos_i = tmp_be_file_i.tellg();

        for(i=0; i<p.num_g_x; i++) //loop over x
        {
            for(j=0; j<p.num_g_y; j++) //loop over y 
            {
                Binding_events be;

                if(my_grid_count == traj.world_num_g[current_rank])
                {
                    //close current temporary BE file
                    tmp_be_file_i.close();

                    //remove tmp file
                    remove(this_be_file_name.c_str()); 

                    //change to the next rank
                    current_rank++;
                    my_grid_count = 0;
           
                    //open next temporary BE file 
                    this_tag          = "_" + to_string(current_rank) + "_tmp.be";
                    this_be_file_name = chop_and_add_tag(p.k_file_name,this_tag);
                    tmp_be_file_i.open(this_be_file_name, ios::out | ios::binary);
                    current_pos_i = tmp_be_file_i.tellg();
                }

                if(grid_count >= traj.world_gi[current_rank] && grid_count <= traj.world_gf[current_rank])
                {
                    current_pos_i = be.get_binding_events_tmp(tmp_be_file_i,current_pos_i);                     
                }
                else
                {
                    printf("This should not have failed! current_rank %4d grid_count %7d world_num_g %7d world_gi %7d world_gf %7d \n",current_rank,grid_count,traj.world_num_g[current_rank],traj.world_gi[current_rank],traj.world_gf[current_rank]);
                }

                //store the position of the grid point 
                grid_be_pos[i][j] = current_pos_o;
          
                //write binding event data to the final file
                current_pos_o = be.write_binding_events_tmp(be_file_o,current_pos_o);

                grid_count++;
                my_grid_count++;

                //report progress
                int ef_g     = grid_count;
                int ef_num_g = p.num_g_x*p.num_g_y;
                time_stats(s.t,&s.counter,ef_g,ef_num_g,s.world_rank);
            }
        }

        //remove last tmp file
        remove(this_be_file_name.c_str());

        //close output file
        be_file_o.close(); 

        //write out .info file for the binding events file
        string info_tag = ".be.info";
        string info_file_name = chop_and_add_tag(p.k_file_name,info_tag);
        FILE *this_file = fopen(info_file_name.c_str(), "w");
        if(this_file == NULL)
        {
            printf("Could not open file %s. \n",info_file_name.c_str());
        }
        else
        {
            for(i=0; i<p.num_g_x; i++) //loop over x
            {
                for(j=0; j<p.num_g_y; j++) //loop over y 
                {
                    fprintf(this_file," %ld ",grid_be_pos[i][j]);
                }
                fprintf(this_file,"\n");
            }
            fclose(this_file);
        }
    }

    //wait until all cores are here
    MPI_Barrier(MPI_COMM_WORLD);

    //log time spent writing BE files
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Write BE Files");

    //reset the clock
    s.t = clock();

    MPI_Barrier(MPI_COMM_WORLD);
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
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                                        s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                                        s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                                         s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                                                    s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                                     s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                                       s.world_rank, s.cl_tags, &p.b_lsq,       0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                                         s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",                         s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_s(argc,argv,"-crd",    p.param_file_name,            "Selection card with target lipids for voronoi diagrame + mapping atoms (crd)",      s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_s(argc,argv,"-k",      p.k_file_name,                "Output filename used to derive filenames for kinetics data (.dat)",                 s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,           "PDB file with sorted leaflets (pdb)",                                               s.world_rank, s.cl_tags, &p.b_lf_pdb,    0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,  "File with additional leaflet finder parameters (prm)",                              s.world_rank, s.cl_tags, &p.b_lf_param,  0);
    add_argument_mpi_s(argc,argv,"-pf_pdb", p.pf_pdb_file_name,           "PDB file with protein marked (pdb)",                                                s.world_rank, s.cl_tags, &p.b_pf_pdb,    0);
    add_argument_mpi_s(argc,argv,"-pf_prm", p.protein_finder_param_name,  "File with additional protein finder parameters (prm)",                              s.world_rank, s.cl_tags, &p.b_pf_param,  0);
    add_argument_mpi_d(argc,argv,"-APS",    &p.APS,                       "Area per grid box (nm^2)",                                                          s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_d(argc,argv,"-r",      &p.radius,                    "Stamping radius when computing the sample count (nm)",                              s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_d(argc,argv,"-bx",     &p.box_x,                     "Grid x dimension (nm)",                                                             s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_d(argc,argv,"-by",     &p.box_y,                     "Grid y dimension (nm)",                                                             s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_i(argc,argv,"-leaf",   &p.leaflet,                   "Which leaflet? (0:both 1:upper 2:lower)",                                           s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_i(argc,argv,"-odf",    &p.out_data_format,           "What is the output data format? (0:matrix 1:vector)",                               s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_i(argc,argv,"-range",  &p.range,                     "Noise filter half width (frames)",                                                  s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_d(argc,argv,"-dt",     &p.delta_t,                   "Effective time step between trajectory frames analyzed (ps) ",                      s.world_rank, s.cl_tags, nullptr,        1);
    add_argument_mpi_d(argc,argv,"-voro_r", &p.voro_stamp_rad,            "Stamping radius used when computing Voronoi diagrams (nm)",                         s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_i(argc,argv,"-v_prot", &p.v_prot,                    "Include protein atoms in voronoi diagram? (0:no 1:yes)",                            s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_d(argc,argv,"-c_dist", &p.c_dist,                    "Distance cutoff for counting protein atoms in voronoi diagram (nm)",                s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_i(argc,argv,"-clean",  &p.b_clean,                   "Remove voronoi diagrams after writing binding events files? (0:no 1:yes)",          s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_i(argc,argv,"-dump",   &p.dump,                      "Dump bound lipids on last frame? (0:no 1:yes)",                                     s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_i(argc,argv,"-com",    &p.com,                       "Use the lipid center of mass for voronoi tessellations? (0:no 1:yes)",              s.world_rank, s.cl_tags, nullptr,        0);
    add_argument_mpi_s(argc,argv,"-mask"  , p.prot_mask_file_name,        "Protein mask (matrix format) to be used in tessellation computations (dat) ",       s.world_rank, s.cl_tags, &p.b_prot_mask, 0);
    add_argument_mpi_d(argc,argv,"-mem",    &p.ram,                       "How much memory to use for storing tessellation data (MB). Rest goes on the disc.", s.world_rank, s.cl_tags, nullptr,        0);
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
    traj.parallelize_by_grid_alt(p.num_g_x,p.num_g_y);

    //count frames to store tessellation data
    int    grid_size      = p.num_g_x*p.num_g_y;
    double mem_per_grid   = (double)grid_size*4.0/1000000.0;
    int    num_frames_ram = floor(p.ram/mem_per_grid);
    if(num_frames_ram > traj.get_num_frames())
    {
        num_frames_ram = traj.get_num_frames();
    }

    //estimate memory requirements  
    double memory_filter = ((double)p.num_g_y*(double)p.num_g_x*(double)(2.0*(double)p.range + 1.0)*4.0)/1000000.0;  //filter
    double memory_bound  = ((double)p.num_g_y*(double)p.num_g_x*(double)p.num_lipids*4.0)/1000000.0;                 //bound
    double memory_voro   = ((double)num_frames_ram*(double)grid_size*4.0)/1000000.0;                                 //voro
    double memory_rho_t  = ((double)p.num_g_y*(double)p.num_g_x*8.0)/1000000.0;                                      //rho_t
    if(s.world_rank == 0)
    {
        printf("Memory Estimates per core (not including memory required to store binding events): \n");
        printf("  bound   %f MB \n",memory_bound);
        printf("  filter  %f MB \n",memory_filter);
        printf("  voro    %f MB \n",memory_voro);
        printf("  rho     %f MB \n",memory_rho_t);
        printf("  total   %f MB \n\n",memory_bound + memory_filter + memory_voro + memory_rho_t);
    }

    //allocate memory to hold tessellation data
    iv3d voro_ram(num_frames_ram, iv2d(p.num_g_x, iv1d(p.num_g_y,0)));                 //holds the voronoi tessellations for as many frames as possible

    //determine which frames to store in memory
    iv2d mem_tags        = tag_frames_mem       (traj,s,p,num_frames_ram);
    iv2d mem_tags_global = tag_frames_mem_global(traj,s,p,num_frames_ram);

    //create multi-dimensional vectors to store data 
    iv3d filter     (p.num_g_y,iv2d(p.num_g_x,iv1d(2*p.range + 1,0)));                //contains a history of which lipid was closest to each grid point
    iv3d bound      (p.num_g_y,iv2d(p.num_g_x,iv1d(p.num_lipids,0)));                 //keeps track of how many frames a lipid has been bound for 

    //allocate other memory
    dv2d rho_t(p.num_g_y,dv1d(p.num_g_x,0.0));                                        //stores the total lipid density

    //store the binding events for each lattice point
    iv3d bi(p.num_g_y,iv2d(p.num_g_x,iv1d(0,0)));
    iv3d bf(p.num_g_y,iv2d(p.num_g_x,iv1d(0,0)));
    iv3d be_lipid_nr(p.num_g_y,iv2d(p.num_g_x,iv1d(0,0)));

    //create vectors to hold the residue number and name for a given lipid number 
    iv1d res_nr_lip(p.num_lipids,0);
    sv1d res_name_lip(p.num_lipids);

    //store the lipid number for each residue for fast lookup
    iv1d lip_nr_res(traj.res_nr[traj.atoms()-1],0);

    //create grid for reading in voronoi diagrams
    iv2d voro(p.num_g_y,iv1d(p.num_g_x,0));

    //create a grid to use when writing voronoi diagrams to file
    Grid_i nan;

    //get the grid dimensions
    nan.set_dim(p.APS,p.num_g_x,p.num_g_y);

    //read in the protein mask
    Grid_lt prot_mask;
    if(p.b_prot_mask == 1)
    { 
        prot_mask.set_format(0);
        prot_mask.get_grid(p.prot_mask_file_name);
    }

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

        get_voronoi(traj,s,p,param,nan,prot_mask,mem_tags,voro_ram);

        update_density(traj,s,p,rho_t,param);

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
        fflush(stdin);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //read voronoi diagrams and compute residence times
    perf.log_time(get_dwell_times(traj,s,p,bound,filter,voro,bi,bf,be_lipid_nr,res_nr_lip,param,res_name_lip,lip_nr_res,mem_tags_global,voro_ram),"Dwell time");

    //write binding events files
    finalize_analysis(traj,s,p,bi,bf,be_lipid_nr,res_nr_lip,res_name_lip,rho_t,nan,perf,mem_tags);

    //print the performance stats
    perf.print_stats();

    //print closing statements
    print_closing(s.world_rank);

    //relinquish the mpi environment
    MPI_Finalize();

    return 0;
}

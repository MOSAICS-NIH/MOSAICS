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
#include <iomanip>
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
#include "MosAT/program_variables/pv_membrane_thickness.h"  //This has the variables specific to the analysis program
#include "headers/array.h"                                   //This has routines used for working with arrays
#include "headers/performance.h"                             //This has a class for logging performance data
#include "headers/index.h"                                   //This has a class for working with index files
#include "headers/traj.h"                                    //This has a class for working with the trajectory
#include "headers/leaflet_finder.h"                          //This has routines used to find leaflets in membrane simulations
#include "headers/protein_finder.h"                          //This has routines used to find protein atoms
#include "headers/sol_finder.h"                              //This has routines used to find the solvent
#include "headers/grid.h"                                    //This has routines used for working with a grid
#include "headers/protein.h"                                 //This has routines used for working with protein data
#include "headers/force_serial.h"                            //This has routines used for forcing the code to run on a single mpi process
#include "headers/param.h"                                   //This has routines used for reading complex parameter data

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function checks if dz makes sense given the leaflet selection                                        //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void validate_delta(program_variables &p,double dif_z,int lip_1,int lip_2,int global_frame)
{
    if(p.leaflet == 1 && dif_z > 0)
    {
        printf("found inconsistent delta z. Perhaps there is something wrong with the trajectory. lipid_1 %10d lipid_2 %10d global_frame %10d \n",lip_1,lip_2,global_frame);
    }
    else if(p.leaflet == 2 && dif_z < 0)
    {
        printf("found inconsistent delta z. Perhaps there is something wrong with the trajectory. lipid_1 %10d lipid_2 %10d global_frame %10d  \n",lip_1,lip_2,global_frame);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function checks if the distance is in a range of interest                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void check_dist(Trajectory &traj,program_variables &p,double min_dist_z,int min_1,int max_1,int min_dist_resnr,int *dump_frame)
{
    int i = 0;  //standard variable used in loops

    if(p.lower_cutoff != p.upper_cutoff) //both set to zero by default
    {
        if(min_dist_z > p.lower_cutoff && min_dist_z < p.upper_cutoff) //in range of interest
        {
            *dump_frame = 1;

            //set beta to 1 
            for(i=min_1; i<=max_1; i++) //loop over current residue
            {
                traj.beta[i] = 1;
            }

            //set beta for the second lipid
            int min_2 = traj.res_start[min_dist_resnr-1];         //index for the first atom of residue i
            int max_2 = traj.res_end[min_dist_resnr-1];           //index for the last atom + 1 of residue i
            for(i=min_2; i<=max_2; i++) //loop over current residue
            {
                traj.beta[i] = 1;
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes membrane thickness and adds to the grid                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void measure_thickness(Trajectory &traj,system_variables &s,program_variables &p,Param &param,Grid &rho,Grid &rho_t,Grid &dist,dv1d &thickness)
{
    int i          = 0;                              //standard variable used in loops
    int j          = 0;                              //standard variable used in loops
    int k          = 0;                              //standard variable used in loops
    int l          = 0;                              //standard variable used in loops
    int dump_frame = 0;                              //Used to decide if the current frame should be written to pdb

    //clear the current frame grids
    rho.clean_frame();
    rho_t.clean_frame();
    dist.clean_frame();

    //set beta values to zero. set if distance falls in special range
    for(i=0; i<traj.atoms(); i++)
    {
        traj.beta[i] = 0;
    }

    for(i=0; i<traj.target_leaflet.size(); i++) //loop over the target membrane atoms
    {
        //jump to the next lipid
        i = traj.next_target_lipid(i);

        for(j=0; j<param.main_size_y(); j++) //loop over lipid types        
        {
            if(strcmp(traj.res_name[traj.target_leaflet[i]-1].c_str(), param.param_main_s[j][0].c_str() ) == 0) //lipid type is correct
            {
                //get the first and last atom of the current lipid
                int    min_1          = traj.t_lip_start(i);
                int    max_1          = traj.t_lip_end(i);
                int    found_neighbor = 0;                                  //keeps track of whether a pair is found or not
                int    min_dist_resnr = -1;                                 //residue number of lipid with min dist_xy
                double min_dist_xy    = 9999;                               //min dist in xy between centers
                double min_dist_z     = 9999;                               //min dist in z between centers

                //compute the center for the lipid
                sv1d target_atoms_1 = param.get_column_sec_s(j,0);
                dv1d r_center_1     = traj.center(target_atoms_1,min_1,max_1);

                for(k=0; k<traj.opposing_leaflet.size(); k++) //loop over the opposing leaflet atoms
                {
                    //jump to the next lipid
                    k = traj.next_opposing_lipid(k);

                    for(l=0; l<param.main_size_y(); l++) //loop over lipid types 
                    {
                        if(strcmp(traj.res_name[traj.opposing_leaflet[k]-1].c_str(), param.param_main_s[l][0].c_str() ) == 0 || p.count_all == 1 ) //lipid type is correct
                        {
                             int min_2 = traj.o_lip_start(k);
                             int max_2 = traj.o_lip_end(k);

                            //compute the center for the lipid
                            sv1d target_atoms_2 = param.get_column_sec_s(l,0);
                            dv1d r_center_2     = traj.center(target_atoms_2,min_2,max_2);

                            //compute distance between centers
                            double dif_x = r_center_2[0] - r_center_1[0];
                            double dif_y = r_center_2[1] - r_center_1[1];
                            double dif_z = r_center_2[2] - r_center_1[2];

                            double distance_z  = sqrt(dif_z*dif_z);
                            double distance_xy = sqrt(dif_x*dif_x + dif_y*dif_y);

                            if(distance_xy < p.xy_cutoff)
                            {
                                found_neighbor = 1;

                                if(distance_xy < min_dist_xy)
                                {
                                    min_dist_xy    = distance_xy;
                                    min_dist_z     = distance_z;
                                    min_dist_resnr = traj.res_nr[traj.opposing_leaflet[k]-1];
                                }

                                //check that upper lipid is above lower lipid
                                validate_delta(p,dif_z,traj.res_nr[min_1],traj.res_nr[min_2],traj.get_frame_global());
                            }
                            goto end_loop;    
                        }
                    }
                    end_loop:;
                }

                //add min dist to grid if a neighbor was found
                if(found_neighbor == 1)
                {
                    //find the mapping atoms and add min dist to the grid
                    for(k=min_1; k<=max_1; k++) //loop over current residue
                    {
                        if(strcmp(traj.atom_name[k].c_str(), param.param_main_s[j][1].c_str()) == 0) //atom is mapping atom 1
                        {
                            double rx = traj.r[k][0];
                            double ry = traj.r[k][1];

                            //add the smallest distance to the grid
                            dist.stamp(rx,ry,p.radius,min_dist_z);
                        }
                        if(strcmp(traj.atom_name[k].c_str(), param.param_main_s[j][2].c_str()) == 0) //atom is mapping atom 2
                        {
                            double rx = traj.r[k][0];
                            double ry = traj.r[k][1];

                            //add the smallest distance to the grid
                            dist.stamp(rx,ry,p.radius,min_dist_z);
                        }
                    }

                    //store thickness measurement
                    thickness.push_back(min_dist_z);

                    //check if distance is in a range of interest
                    check_dist(traj,p,min_dist_z,min_1,max_1,min_dist_resnr,&dump_frame);
                }

                //add lipid to grid for computing rho_t
                for(l=min_1; l<=max_1; l++) //loop over current residue
                {
                    if(strcmp(traj.atom_name[l].c_str(), param.param_main_s[j][1].c_str()) == 0) //atom is mapping atom 1
                    {
                        double rx = traj.r[l][0];
                        double ry = traj.r[l][1];

                        rho_t.stamp(rx,ry,p.radius,1.0);
                    }
                    if(strcmp(traj.atom_name[l].c_str(), param.param_main_s[j][2].c_str()) == 0) //atom is mapping atom 2
                    {
                        double rx = traj.r[l][0];
                        double ry = traj.r[l][1];

                        rho_t.stamp(rx,ry,p.radius,1.0);
                    }
                }
            }
        }
    }

    //get the average for the current frame
    dist.norm_frame();
    rho_t.norm_frame();

    //add the current frame grid to long term sum
    dist.add_frame();
    rho_t.add_frame();

    //now we print the single frame thickness data
    if(p.b_stdev == 1)
    {
        dist.write_frame(traj.get_frame_global());
    }

    //dump pdb with highlighted lipid pairs of interest
    if(dump_frame == 1)
    {
        string dump_file_name = "frame_dump_" + to_string(traj.get_frame_global()) + ".pdb";
        write_data_pdb(traj.box,traj.atoms(),traj.atom_nr,traj.res_nr,traj.res_name,traj.atom_name,traj.r,traj.time,traj.step,traj.beta,traj.weight,traj.chain_id,traj.element,
                       traj.title,s.world_rank,dump_file_name);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Collect distances and compute average                                                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,Grid &rho,Grid &rho_t,Grid &dist,Grid &pairs_percent,dv1d &thickness)
{
    int i = 0;           //standard variable used in loops
    int j = 0;           //standard variable used in loops

    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nfinalizing analysis. This requires communication of the grid and could take some time. \n");
    }

    //collect dist and rho from all ranks
    dist.collect_grid();
    rho_t.collect_grid();

    //compute percentage of target lipids to pair up
    for(i=0; i<rho_t.num_x(); i++) //loop over x
    {
        for(j=0; j<rho_t.num_y(); j++) //loop over y
        {
            pairs_percent.grid[j][i] = dist.rho[j][i]/rho_t.rho[j][i];
            pairs_percent.rho[j][i]  = rho_t.rho[j][i];
        }
    }

    //normalize dist and write the <dist> and rho to file
    if(s.world_rank == 0)
    {
        dist.normalize();

        dist.exclude_data(p.cutoff,0);
        rho_t.exclude_data(p.cutoff,1);
        pairs_percent.exclude_data(p.cutoff,0);

        dist.write_grid();
        pairs_percent.write_grid();
        dist.write_rho();
        rho_t.write_grid();
    }

    //compute standard deviation
    dist.get_stdev(p.b_stdev,p.b_clean,traj);

    //here we collect the thickness measurements from each core
    collect_dv1d(s.world_size,s.world_rank,thickness);

    //now we find the free energy 
    if(s.world_rank == 0)
    {
        //now we find the largest thickness
        double largest_thickness = 0;
        for(i=0; i<thickness.size(); i++)
        {
            if(thickness[i] > largest_thickness)
            {
                largest_thickness = thickness[i];
            }
        }

        //now we bin the data
        int  num_bins = (int)ceil(largest_thickness/(p.bin_width));
        iv1d bins(num_bins,0);
        for(i=0; i<thickness.size(); i++) //loop over thickness measurements
        {
            for(j=0; j<num_bins; j++) //loop over num bins
            {
                if(thickness[i] >= (double)j*(p.bin_width) && thickness[i] < (double)(j+1)*(p.bin_width))
                {
                    //check the values being binned
                    //if(j == 300)
                    //{
                    //    printf("thickness[%d] %12.10f lower %12.10f upper %12.10f \n",i,thickness[i],(double)j*(p.bin_width),(double)(j+1)*(p.bin_width));
                    //}

                    bins[j] = bins[j] + 1;
                }
            }
        }

        //get the most populated bin (used to shift F minimum to zero)
        int biggest_bin = 0;
        for(i=0; i<num_bins; i++)
        {
            if(bins[i] > biggest_bin)
            {
                biggest_bin = bins[i];
            }
        }

        //print data to output file
        string f_file_name = add_tag(p.thk_file_name,"_free_energy");
        FILE *f_file = fopen(f_file_name.c_str(), "w");
        if(f_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",f_file_name.c_str());
        }
        else
        {
            double sum_prob    = 0;
            double avagadro    = 6.0221409 * pow(10,23);
            double k_boltzmann = 1.38064852 * pow(10,-23);
            double beta        = 1/(k_boltzmann*p.temp);

            fprintf(f_file," #%10s %15f \n","T",p.temp);
            fprintf(f_file," #%10s %15f \n","Beta",beta);
            fprintf(f_file," #%10s %15s \n","Column_1","Thickness (nm)");
            fprintf(f_file," #%10s %15s \n","Column_2","Count");
            fprintf(f_file," #%10s %15s \n","Column_3","Probability");
            fprintf(f_file," #%10s %15s \n","Column_4","F (kJ/mol)");
            fprintf(f_file," #%10s %15s \n","Column_5","F_shifted (kJ/mol)");

            for(i=0; i<num_bins; i++)
            {
                double probability = (double)bins[i]/(double)thickness.size();
                double F           = -(1/beta)*log(probability)*avagadro/1000.0;
                double F_shifted   = -(1/beta)*log((double)bins[i]/(double)biggest_bin)*avagadro/1000.0;
                if(bins[i] == 0)
                {
                    fprintf(f_file,"%10.4f %10d %10.8f %10s %10s \n",(double)i*p.bin_width,bins[i],probability,"NaN","NaN");
                }
                else
                {
                    fprintf(f_file,"%10.4f %10d %10.8f %10.8e %10.8e \n",(double)i*p.bin_width,bins[i],probability,F,F_shifted);
                }

                sum_prob = sum_prob + probability; 
            }
            fclose(f_file);
            printf("Summing probability over bins. sum_prob = %f \n",sum_prob);
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
    s.program_name = "Membrane Thickness";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                  s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                  s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                 s.world_rank, s.cl_tags, &p.b_print,   0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                   s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                              s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                               s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                 s.world_rank, s.cl_tags, &p.b_lsq,     0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                   s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",   s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-crd",    p.param_file_name,            "Selection card with lipid types and center atoms (crd)",      s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-thk",    p.thk_file_name,              "Output grid data with <thickness> (dat)",                     s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,           "PDB file with sorted leaflets (pdb)",                         s.world_rank, s.cl_tags, &p.b_lf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,  "File with additional leaflet finder parameters (prm)",        s.world_rank, s.cl_tags, &p.b_lf_param,0);
    add_argument_mpi_d(argc,argv,"-APS",    &p.APS,                       "Area per grid box (nm^2)",                                    s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-r",      &p.radius,                    "Radius of the target atom (nm)",                              s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-cutoff", &p.cutoff,                    "Cutoff for excluding data (chi)",                             s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-bx",     &p.box_x,                     "Grid x dimension (nm)",                                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-by",     &p.box_y,                     "Grid y dimension (nm)",                                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-leaf",   &p.leaflet,                   "Which leaflet? (0:both 1:upper 2:lower)",                     s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-odf",    &p.out_data_format,           "What is the output data format? (0:matrix 1:vector)",         s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-stdev",  &p.b_stdev,                   "Compute the STDEV and STEM? (0:no 1:yes)",                    s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-clean",  &p.b_clean,                   "Remove single frame files? (0:no 1:yes)",                     s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-xy",     &p.xy_cutoff,                 "How close in xy before counting lipid in min dist (nm)",      s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-width",  &p.bin_width,                 "Width of bins for free energy profile (nm)",                  s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-dump_u", &p.upper_cutoff,              "Dump frame with thickness below this value (nm)",             s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-dump_l", &p.lower_cutoff,              "Dump frame with thickness above this value (nm)",             s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-temp",   &p.temp,                      "Temperature of the simulation (k)",                           s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-all",    &p.count_all,                 "Use all lipids from opposing leaflet? (0:no 1:yes)",          s.world_rank, s.cl_tags, nullptr,      0);
    conclude_input_arguments_mpi(argc,argv,s.world_rank,s.program_name,s.cl_tags);

    //create a trajectory
    Trajectory traj; 

    //set trajectory parameters
    traj.set_block_parallel(on);
    traj.set_traj(p.in_file_name);
    traj.set_ref(p.ref_file_name);
    traj.set_traj_w(p.out_file_name,p.b_print);
    traj.set_lsq(p.lsq_index_file_name,p.b_lsq,p.lsq_dim,p.lsq_ref);
    traj.set_res(p.stride,p.start_frame,p.end_frame);

    //analyze the trajectory (log time spent) 
    perf.log_time(traj.build(),"Analyze Trajectory");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //check file extensions                                                                                     
    check_extension_mpi(s.world_rank,"-crd",p.param_file_name,".crd");
    check_extension_mpi(s.world_rank,"-thk",p.thk_file_name,".dat");

    if(p.b_lf_pdb == 1)
    {
        check_extension_mpi(s.world_rank,"-lf_pdb",p.lf_pdb_file_name,".pdb");
    }
    if(p.b_lf_param == 1)
    {
        check_extension_mpi(s.world_rank,"-lf_prm",p.leaflet_finder_param_name,".prm");
    }

    //create parameter files
    Param param;

    //read parameter files
    param.get_param(p.param_file_name,4,3,1);

    //check the integrity of the parameter files
    if(param.check_file() == 0) //bad files. kill program 
    {
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }

    //run leaflet finder
    traj.get_leaflets(p.leaflet,p.leaflet_finder_param_name,p.b_lf_param);

    //print a pdb with distinguished leaflets
    traj.write_leaflets(p.lf_pdb_file_name,p.b_lf_pdb);

    //print info about the target leaflet
    traj.get_leaflet_stats();

    //print info about the lipid selection
    traj.get_lipid_selection_stats(param.get_column_s(0),"-crd");

    //create a grid to hold distances etc.
    Grid rho;
    Grid rho_t;
    Grid dist;
    Grid pairs_percent;

    //get the grid dimensions
    rho.get_dim(p.box_x,p.box_y,traj.ibox,p.APS);
    rho_t.get_dim(p.box_x,p.box_y,traj.ibox,p.APS);
    dist.get_dim(p.box_x,p.box_y,traj.ibox,p.APS);
    pairs_percent.get_dim(p.box_x,p.box_y,traj.ibox,p.APS);

    //set the output file name for grid
    string rho_file_name     = add_tag(p.thk_file_name,"_rho");
    string rho_t_file_name   = add_tag(p.thk_file_name,"_rho_t");
    string percent_file_name = add_tag(p.thk_file_name,"_pairs_percent");
    rho.set_output(rho_file_name,p.out_data_format);
    rho_t.set_output(rho_t_file_name,p.out_data_format);
    dist.set_output(p.thk_file_name,p.out_data_format);
    pairs_percent.set_output(percent_file_name,p.out_data_format);

    //print info about the grid
    dist.print_dim();

    //make vector to hold distances
    dv1d thickness(0,0.0);                               //stores each thickness measurement
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

        measure_thickness(traj,s,p,param,rho,rho_t,dist,thickness);

        traj.set_beta_lf();

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect thickness data from mpi processes and compute average
    perf.log_time(finalize_analysis(traj,s,p,rho,rho_t,dist,pairs_percent,thickness),"Fin Ana");

    //splice temporary traj file together (log time spent)
    perf.log_time(traj.finalize_trajectory(),"Finalize Trajectory");

    //print the performance stats
    perf.print_stats();

    //print closing statements
    print_closing(s.world_rank);

    //relinquish the mpi environment
    MPI_Finalize();

    return 0;
}

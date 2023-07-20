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
#include "headers/linear_algebra.h"                          //This has linear algebra stuff
#include "headers/mosat_routines.h"                         //This is where most of the functions called in main are located
#include "headers/file_naming.h"                             //This has routines for added tags to an existing file name    
#include "headers/file_naming_mpi.h"                         //This has routines for added tags to an existing file name (mpi)
#include "headers/command_line_args_mpi.h"                   //This has routines for adding command line arguments
#include "MosAT/program_variables/pv_lipid_orientation.h"   //This has the variables specific to the analysis program
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the orientation vector for the lipids and adds it to the grid                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void lip_orientation(Trajectory &traj,system_variables &s,program_variables &p,Index &pairs,Grid &sum_x,Grid &sum_y,Grid &sum_z)
{
    //use cos_theta = u_dot_v/(mag_u*mag_v) to compute the angle between 2 vecors u and v

    int i     = 0;                              //standard variable used in loops
    int j     = 0;                              //standard variable used in loops
    int k     = 0;                              //standard variable used in loops
    int l     = 0;                              //standard variable used in loops
    int m     = 0;                              //standard variable used in loops
    double hx = 0;                              //mapping atom x coord
    double hy = 0;                              //mapping atom y coord
    double dx = 0;                              //vector connecting head and tail atoms (x component)
    double dy = 0;                              //vector connecting head and tail atoms (y component)
    double dz = 0;                              //vector connecting head and tail atoms (z component)

    //clear the current frame grids
    sum_x.clean_frame();
    sum_y.clean_frame();
    sum_z.clean_frame();

    for(i=0; i<traj.target_leaflet.size(); i++) //loop over target leaflet atoms
    {
        //get the first and last atom of the current lipid
        int min = traj.t_lip_start(i);
        int max = traj.t_lip_end(i);

        //jump to the next lipid
        i = traj.next_target_lipid(i);

        for(j=0; j<pairs.index_s.size(); j+=4) //loop over lipid types
        {
            if(strcmp(traj.res_name[traj.target_leaflet[i]-1].c_str(), pairs.index_s[j].c_str()) == 0) //residue is a target lipid
            {
                for(k=min; k<=max; k++) //loop over current residue atoms
                {
                    if(strcmp(traj.atom_name[k].c_str(), pairs.index_s[j+1].c_str()) == 0) //atom is head 
                    {
                        for(l=min; l<=max; l++) //loop over current residue atoms
                        {
                            if(strcmp(traj.atom_name[l].c_str(), pairs.index_s[j+2].c_str()) == 0) //atom is tail 
                            {
                                for(m=min; m<=max; m++) //loop over current residue atoms
                                {
                                    if(strcmp(traj.atom_name[m].c_str(), pairs.index_s[j+3].c_str()) == 0) //atom is mapping atom
                                    {
                                        hx = traj.r[m][0];
                                        hy = traj.r[m][1];

                                        dx = traj.r[l][0] - traj.r[k][0];
                                        dy = traj.r[l][1] - traj.r[k][1];
                                        dz = traj.r[l][2] - traj.r[k][2];

                                        sum_x.stamp(hx,hy,p.radius,dx);
                                        sum_y.stamp(hx,hy,p.radius,dy);
                                        sum_z.stamp(hx,hy,p.radius,dz);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    //get the average for the current frame
    sum_x.norm_frame();
    sum_y.norm_frame();
    sum_z.norm_frame();

    //add the current frame grid to long term sum
    sum_x.add_frame();
    sum_y.add_frame();
    sum_z.add_frame();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function does communication and computes the average orientation vector. It finally then             //
// characterizes this vector and writes the spherical angels etc. to file                                    //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,Grid &sum_x,Grid &sum_y,Grid &sum_z,Grid &p2,Grid &phi,Grid &theta,Grid &dist)
{
    int    i  = 0;                //standard variable used in loops
    int    j  = 0;                //standard variable used in loops
    int    k  = 0;                //standard variable used in loops
    int    l  = 0;                //standard variable used in loops
    double pi = 3.141592654;      //pi

    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nFinalizing analysis. This requires communicating data over the grid and could take some time depending on the resolution. \n");
    }

    //collect coords 
    sum_x.collect_grid();
    sum_y.collect_grid();
    sum_z.collect_grid();

    //rank 0 prints data to output file
    if(s.world_rank == 0)
    {
        sum_x.normalize();
        sum_y.normalize();
        sum_z.normalize();

        //characterize the average vector across the grid
        for(j=0; j<sum_x.num_x(); j++) //loop over x-dimension
        {
            for(k=0; k<sum_x.num_y(); k++) //loop over y-dimension
            {
                if(sum_x.rho[k][j] > 0) //data present
                {
                    //transfer vector to array for use with dot function
                    double vec[3];
                    vec[0] = sum_x.grid[k][j];
                    vec[1] = sum_y.grid[k][j];
                    vec[2] = sum_z.grid[k][j];

                    //compute theta
                    double z[3] = {0, 0, 1};
                    double mag_vec   = sqrt( dot(vec,vec) );
                    double vec_dot_z = dot(vec,z);
                    double cos_theta = vec_dot_z/mag_vec;
                    theta.grid[k][j] = acos(cos_theta)*(180/pi); //angle in degrees

                    //compute p2
                    p2.grid[k][j] = 0.5*( (3*(cos_theta)*(cos_theta)) - 1);

                    //compute phi
                    double proj_v[3] = {vec[0],vec[1],0};
                    double mag_proj_v = sqrt( dot(proj_v,proj_v) );
                    double y[3] = {0,1,0};
                    double mag_y = sqrt( dot(y,y) );
                    double proj_v_dot_y = dot(proj_v,y);
                    double cos_phi = proj_v_dot_y/(mag_proj_v*mag_y);
                    phi.grid[k][j] = acos(cos_phi)*(180/pi);

                    //if orientation vector is in quadrant 2 or 3 the anlgel (greater than 180 from y) should be corrected 
                    if(vec[0] < 0)
                    {
                        phi.grid[k][j] = 2*(180 - phi.grid[k][j]) + phi.grid[k][j];
                    }

                    //compute r
                    dist.grid[k][j] = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
                }
            }
        }

        //copy lipid density 
        p2.copy_rho(sum_x);
        phi.copy_rho(sum_x);         
        theta.copy_rho(sum_x);         
        dist.copy_rho(sum_x);            
        
        //exclude insignificant data
        sum_x.exclude_data(p.cutoff,0);
        sum_y.exclude_data(p.cutoff,0);
        sum_z.exclude_data(p.cutoff,0);
        p2.exclude_data(p.cutoff,0);
        phi.exclude_data(p.cutoff,0);
        theta.exclude_data(p.cutoff,1);
        dist.exclude_data(p.cutoff,0);

        //create file names for x,y and z
        string x_file_name     = add_tag(p.p2_file_name,"_x");
        string y_file_name     = add_tag(p.p2_file_name,"_y");
        string z_file_name     = add_tag(p.p2_file_name,"_z");
        string p2_file_name    = add_tag(p.p2_file_name,"_p2");
        string phi_file_name   = add_tag(p.p2_file_name,"_phi");
        string theta_file_name = add_tag(p.p2_file_name,"_theta");
        string dist_file_name  = add_tag(p.p2_file_name,"_dist");

        //set the output file name and format for grids
        sum_x.set_output(x_file_name,p.out_data_format);
        sum_y.set_output(y_file_name,p.out_data_format);
        sum_z.set_output(z_file_name,p.out_data_format);
        p2.set_output(p2_file_name,p.out_data_format);
        phi.set_output(phi_file_name,p.out_data_format);
        theta.set_output(theta_file_name,p.out_data_format);
        dist.set_output(dist_file_name,p.out_data_format);

        //write data to the output files 
        sum_x.write_grid();
        sum_y.write_grid();
        sum_z.write_grid();
        p2.write_grid();
        phi.write_grid();
        theta.write_grid();
        dist.write_grid();
        p2.write_rho();
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
    s.program_name = "Lipid Orientation";

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
    add_argument_mpi_s(argc,argv,"-pairs",  p.pairs_file_name,            "Selection card with orientational vector definitions (crd)",  s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-p2",     p.p2_file_name,               "Output file with spatially resolved tilt angle (dat)",        s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,           "PDB file with sorted leaflets (pdb)",                         s.world_rank, s.cl_tags, &p.b_lf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,  "File with additional leaflet finder parameters (prm)",        s.world_rank, s.cl_tags, &p.b_lf_param,0);
    add_argument_mpi_d(argc,argv,"-APS",    &p.APS,                       "Area per grid box (nm^2)",                                    s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-r",      &p.radius,                    "Radius of the target atom (nm)",                              s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-cutoff", &p.cutoff,                    "Cutoff for excluding data (chi)",                             s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-bx",     &p.box_x,                     "Grid x dimension (nm)",                                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-by",     &p.box_y,                     "Grid y dimension (nm)",                                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-leaf",   &p.leaflet,                   "Which leaflet? (0:both 1:upper 2:lower)",                     s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-odf",    &p.out_data_format,           "What is the output data format? (0:matrix 1:vector)",         s.world_rank, s.cl_tags, nullptr,      0);
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
    check_extension_mpi(s.world_rank,"-pairs",p.pairs_file_name,".crd");
    check_extension_mpi(s.world_rank,"-p2",p.p2_file_name,".dat");

    if(p.b_lf_pdb == 1)
    {
        check_extension_mpi(s.world_rank,"-lf_pdb",p.lf_pdb_file_name,".pdb");
    }
    if(p.b_lf_param == 1)
    {
        check_extension_mpi(s.world_rank,"-lf_prm",p.leaflet_finder_param_name,".prm");
    }

    //read in atom selection parameter file
    Index pairs;
    pairs.get_index(p.pairs_file_name);

    //run leaflet finder
    traj.get_leaflets(p.leaflet,p.leaflet_finder_param_name,p.b_lf_param);

    //print a pdb with distinguished leaflets
    traj.write_leaflets(p.lf_pdb_file_name,p.b_lf_pdb);

    //print info about the target leaflet
    traj.get_leaflet_stats();

    //print info about the lipid selection
    traj.get_lipid_selection_stats(pairs.get_column_s(4,0),"-pairs");

    //create grids to hold orientation vector etc.
    Grid sum_x;
    Grid sum_y;
    Grid sum_z;
    Grid p2;
    Grid phi;
    Grid theta;
    Grid dist;

    //get the grid dimensions
    sum_x.get_dim(p.box_x,p.box_y,traj.ibox,p.APS);
    sum_y.get_dim(p.box_x,p.box_y,traj.ibox,p.APS);
    sum_z.get_dim(p.box_x,p.box_y,traj.ibox,p.APS);
    p2.get_dim(p.box_x,p.box_y,traj.ibox,p.APS);
    phi.get_dim(p.box_x,p.box_y,traj.ibox,p.APS);
    theta.get_dim(p.box_x,p.box_y,traj.ibox,p.APS);
    dist.get_dim(p.box_x,p.box_y,traj.ibox,p.APS);

    //print info about the grid
    sum_x.print_dim();
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

        lip_orientation(traj,s,p,pairs,sum_x,sum_y,sum_z);

        traj.set_beta_lf();

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect data from mpi processes and characterize the mean orientation vector
    perf.log_time(finalize_analysis(traj,s,p,sum_x,sum_y,sum_z,p2,phi,theta,dist),"Fin Ana");

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

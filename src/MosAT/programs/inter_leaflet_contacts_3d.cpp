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

#include "xdr/include/xdrfile_xtc.h"                               //used to read xtc files 
#include "xdr/include/xdr_seek.h"                                  //used to get and set the file position in xtc and trr files
#include "xdr/include/xdrfile_trr.h"                               //used to read trr files
#include "xdr/include/xdrfile.h"                                   //used to read xtc and trr files
#include "xdr/include/trr_header.h"                                //used to read the header info of trr files
#include "headers/multi_dim_vec.h"                                 //This defines multidimensional vectors
#include "headers/switch.h"                                        //This defines a switch (on, off)
#include "headers/file_reader.h"                                   //This has basic routines for reading text files
#include "headers/vector_mpi.h"                                    //This has routines for collecting vector data
#include "headers/mosat_routines.h"                               //This is where most of the functions called in main are located
#include "headers/file_naming.h"                                   //This has routines for added tags to an existing file name    
#include "headers/file_naming_mpi.h"                               //This has routines for added tags to an existing file name (mpi)
#include "headers/command_line_args_mpi.h"                         //This has routines for adding command line arguments
#include "MosAT/program_variables/pv_inter_leaflet_contacts_3d.h" //This has the variables specific to the analysis program
#include "headers/array.h"                                         //This has routines used for working with arrays
#include "headers/performance.h"                                   //This has a class for logging performance data
#include "headers/index.h"                                         //This has a class for working with index files
#include "headers/traj.h"                                          //This has a class for working with the trajectory
#include "headers/leaflet_finder.h"                                //This has routines used to find leaflets in membrane simulations
#include "headers/protein_finder.h"                                //This has routines used to find protein atoms
#include "headers/sol_finder.h"                                    //This has routines used to find the solvent
#include "headers/grid.h"                                          //This has routines used for working with a grid
#include "headers/grid_3d.h"                                       //This has routines used for working with a 3d grid
#include "headers/protein.h"                                       //This has routines used for working with protein data
#include "headers/force_serial.h"                                  //This has routines used for forcing the code to run on a single mpi process
#include "headers/param.h"                                         //This has routines used for reading complex parameter data

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the inter-leaflet contacts and stamps it onto the grid                             //
//                                                                                                           //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_contacts(Trajectory &traj,system_variables &s,program_variables &p,Param &param_1,Param &param_2,Grid_3d &lfc,iv1d &contacts,dv1d &contacts_frame)
{
    int i               = 0;          //standard variable used in loops
    int j               = 0;          //standard variable used in loops
    int k               = 0;          //standard variable used in loops
    int l               = 0;          //standard variable used in loops
    int m               = 0;          //standard variable used in loops
    int n               = 0;          //standard variable used in loops
    int o               = 0;          //standard variable used in loops
    int q               = 0;          //standard variable used in loops
    int measurements    = 0;          //how many measurements for the single frame
    double contacts_avg = 0.0;        //the average number of contacts for the frame

    //clear the current frame grids
    lfc.clean_frame();

    for(i=0; i<traj.target_leaflet.size(); i++) //loop over the target leaflet atoms
    {
        //get the first and last atom of the current lipid
        int min = traj.t_lip_start(i);
        int max = traj.t_lip_end(i);

        //jump to the next lipid
        i = traj.next_target_lipid(i);

        for(j=0; j<param_1.main_size_y(); j++) //loop over lipid types 1
        {
            if(strcmp(traj.res_name[min].c_str(), param_1.param_main_s[j][0].c_str() ) == 0) //lipid 1 type is correct
            {
                int contact_count = 0;

                for(k=min; k<=max; k++) //loop over current residue atoms
                {
                    for(l=0; l<param_1.sec_size_y(j); l++) //loop over target atoms
                    {
                        if(strcmp(traj.atom_name[k].c_str(), param_1.param_sec_s[j][l][0].c_str() ) == 0) //atom is a target atom
                        {
                            for(m=0; m<traj.opposing_leaflet.size(); m++) //loop over the opposing membrane atoms
                            {
                                //get the first and last atom of the current lipid
                                int min_2 = traj.o_lip_start(m);
                                int max_2 = traj.o_lip_end(m);

                                //jump to the next lipid
                                m = traj.next_opposing_lipid(m);

                                for(n=0; n<param_2.main_size_y(); n++) //loop over lipid types 2
                                {
                                    if(strcmp(traj.res_name[min_2].c_str(), param_2.param_main_s[n][0].c_str() ) == 0) //lipid 2 type is correct
                                    {
                                        for(o=min_2; o<=max_2; o++) //loop over current residue atoms
                                        {
                                            for(q=0; q<param_2.sec_size_y(n); q++) //loop over target atoms
                                            {
                                                if(strcmp(traj.atom_name[o].c_str(), param_2.param_sec_s[n][q][0].c_str() ) == 0) //atom is a target atom
                                                {
                                                    //compute the distance between the atoms
                                                    double dx = traj.r[k][0] - traj.r[o][0];
                                                    double dy = traj.r[k][1] - traj.r[o][1];
                                                    double dz = traj.r[k][2] - traj.r[o][2];

                                                    double distance = sqrt(dx*dx + dy*dy + dz*dz);

                                                    if(distance < p.contact_cutoff) //contact
                                                    {
                                                        contact_count++;
                                                    }
                                                    goto end_loop_op;
                                                }
                                            }
                                            end_loop_op:;
                                        }
                                    }
                                }
                            }
                            goto end_loop_target;
                        }
                    }
                    end_loop_target:;
                }

                //find mapping atoms
                for(k=min; k<=max; k++) //loop over current residue atoms
                {
                    int b_stamp = 0;
                    double hx   = 0.0;
                    double hy   = 0.0;
                    double hz   = 0.0;

                    if(strcmp(traj.atom_name[k].c_str(), param_1.param_main_s[j][1].c_str() ) == 0) //mapping atom 1
                    {
                        hx      = traj.r[k][0];
                        hy      = traj.r[k][1];
                        hz      = traj.r[k][2];
                        b_stamp = 1;
                    }
                    else if(strcmp(traj.atom_name[k].c_str(), param_1.param_main_s[j][2].c_str() ) == 0) //mapping atom 2
                    {
                        hx      = traj.r[k][0];
                        hy      = traj.r[k][1];
                        hz      = traj.r[k][2];
                        b_stamp = 1;
                    }

                    if(b_stamp == 1)
                    {
                        lfc.stamp(hx,hy,hz,p.radius,contact_count);

                        //add inter leaflet contacts measurement to contacts 
                        if(hx > 0.15*traj.ibox[XX][XX] && hx < 0.85*traj.ibox[XX][XX] && hy > 0.15*traj.ibox[YY][YY] && hy < 0.85*traj.ibox[YY][YY]) //ignore measurements at boundaries
                        {
                            contacts.push_back(contact_count);
                            contacts_avg = contacts_avg + contact_count;
                            measurements++;
                        }
                    }
                }
            }
        }
    }

    //compute average number of contacts for the current frame
    contacts_avg = contacts_avg/(double)measurements;
    contacts_frame.push_back(contacts_avg);

    //get the average for the current frame
    lfc.norm_frame();

    //add the interleaflet contacts for each gridpoint to the long term sum
    lfc.add_frame();

    //now we print the single frame interleaflet contacts
    if(p.b_stdev == 1)
    {
        lfc.write_frame(traj.get_frame_global(),p.ex_val);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Collect inter-leaflet contacts and compute the average.                                                   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,Grid_3d &lfc,iv1d &contacts,dv1d &contacts_frame)
{
    int i = 0; //standard variable used in loops
    int j = 0; //standard variable used in loops 

    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nfinalizing analysis. This requires communication of the grid and could take some time. \n");
    }

    //collect lfc and rho from all ranks
    lfc.collect_grid();

    //normalize lfc and write the <lfc> and rho to file
    if(s.world_rank == 0)
    {
        lfc.normalize();

        lfc.exclude_data(p.cutoff,1);

        lfc.write_grid(p.ex_val);
        lfc.write_rho(p.ex_val);
    }

    //compute standard deviation
    lfc.get_stdev(p.b_stdev,p.b_clean,traj,p.ex_val);

    if(p.b_histo == 1) //analyze the free energy 
    {
        //here we collect the thickness measurements from each core
        collect_iv1d(s.world_size,s.world_rank,contacts);

        //now we find the free energy 
        if(s.world_rank == 0)
        {
            //now we find the largest number of contacts
            double largest_contacts = 0;
            for(i=0; i<contacts.size(); i++)
            {
                if(contacts[i] > largest_contacts)
                {
                    largest_contacts = contacts[i];
                }
            }

            //now we bin the data
            int  num_bins = (int)ceil(largest_contacts/(p.bin_width));
            iv1d bins(num_bins,0);
            for(i=0; i<contacts.size(); i++) //loop over contacts measurements
            {
                for(j=0; j<num_bins; j++) //loop over num bins
                {
                    if(contacts[i] >= (double)j*(p.bin_width) && contacts[i] < (double)(j+1)*(p.bin_width))
                    {
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
            string f_file_name = add_tag(p.lfc_file_name,"_free_energy");
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

                fprintf(f_file,"# %10s %15f \n","T",p.temp);
                fprintf(f_file,"# %10s %15f \n","Beta",beta);
                fprintf(f_file,"# %10s %15s \n","Column_1","Contacts");
                fprintf(f_file,"# %10s %15s \n","Column_2","Count");
                fprintf(f_file,"# %10s %15s \n","Column_3","Probability");
                fprintf(f_file,"# %10s %15s \n","Column_4","F (kJ/mol)");
                fprintf(f_file,"# %10s %15s \n","Column_5","F_shifted (kJ/mol)");

                for(i=0; i<num_bins; i++)
                {
                    double probability = (double)bins[i]/(double)contacts.size();
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
    }

    //here we collect the thickness measurements from each core
    collect_dv1d(s.world_size,s.world_rank,contacts_frame);

    if(s.world_rank == 0)
    {
        //print to file the average number of contacts for each frame
        string frame_file_name = chop_and_add_tag(p.lfc_file_name,"_contacts_frame") + ".dat";
        FILE *frame_file = fopen(frame_file_name.c_str(), "w");
        if(frame_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",frame_file_name.c_str());
        }
        else
        {
            fprintf(frame_file,"# column 1: frame number \n");
            fprintf(frame_file,"# column 2: contacts formed per chemical group (averaged over the target lipids) \n");
            for(i=0; i<contacts_frame.size(); i++) //loop over frames
            {
                fprintf(frame_file," %10d %10.4f \n",i,contacts_frame[i]);
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
    s.program_name = "Inter-Leaflet Contacts 3D";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                  s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                  s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                 s.world_rank, &p.b_print,   0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                   s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                              s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                               s.world_rank, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                 s.world_rank, &p.b_lsq,     0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                   s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",   s.world_rank, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-crd_1",  p.param_1_file_name,          "Selection card target leaflet (crd)",                         s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-crd_2",  p.param_2_file_name,          "Selection card opposing leaflet (crd)",                       s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lfc",    p.lfc_file_name,              "Output file with spatially resolved ILC (dx)",                s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,           "PDB file with sorted leaflets (pdb)",                         s.world_rank, &p.b_lf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,  "File with additional leaflet finder parameters (prm)",        s.world_rank, &p.b_lf_param,0);
    add_argument_mpi_d(argc,argv,"-APS",    &p.APS,                       "Area per grid box (nm^2)",                                    s.world_rank, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-r",      &p.radius,                    "Radius of the target atom (nm)",                              s.world_rank, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-cutoff", &p.cutoff,                    "Cutoff for excluding data (chi)",                             s.world_rank, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-bx",     &p.box_x,                     "Grid x dimension (nm)",                                       s.world_rank, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-by",     &p.box_y,                     "Grid y dimension (nm)",                                       s.world_rank, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-bz",     &p.box_z,                     "Grid z dimension (nm)",                                       s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-leaf",   &p.leaflet,                   "Which leaflet? (0:both 1:upper 2:lower)",                     s.world_rank, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-cdist",  &p.contact_cutoff,            "The contact cutoff distance (nm)",                            s.world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-stdev",  &p.b_stdev,                   "Compute the STDEV and STEM? (0:no 1:yes)",                    s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-clean",  &p.b_clean,                   "Remove single frame files? (0:no 1:yes)",                     s.world_rank, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-width",  &p.bin_width,                 "Width of bins for free energy profile (contacts)",            s.world_rank, &p.b_histo,   0);
    add_argument_mpi_d(argc,argv,"-temp",   &p.temp,                      "Temperature of the simulation (k)",                           s.world_rank, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-ex_val", &p.ex_val,                    "Set excluded lattice points to this value",                   s.world_rank, nullptr,      0);
    conclude_input_arguments_mpi(argc,argv,s.world_rank,s.program_name);

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
    check_extension_mpi(s.world_rank,"-crd_1",p.param_1_file_name,".crd");
    check_extension_mpi(s.world_rank,"-crd_2",p.param_2_file_name,".crd");
    check_extension_mpi(s.world_rank,"-lfc",p.lfc_file_name,".dx");

    if(p.b_lf_pdb == 1)
    {
        check_extension_mpi(s.world_rank,"-lf_pdb",p.lf_pdb_file_name,".pdb");
    }
    if(p.b_lf_param == 1)
    {
        check_extension_mpi(s.world_rank,"-lf_prm",p.leaflet_finder_param_name,".prm");
    }

    //create parameter files
    Param param_1;
    Param param_2; 

    //read parameter files
    param_1.get_param(p.param_1_file_name,4,3,1);
    param_2.get_param(p.param_2_file_name,2,1,1);

    //check the integrity of the parameter files
    if(param_1.check_file() == 0 || param_2.check_file() == 0) //bad files. kill program 
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
    traj.get_lipid_selection_stats(param_1.get_column_s(0),"-crd_1");
    traj.get_lipid_selection_stats(param_2.get_column_s(0),"-crd_2");

    //create a grid to hold enrichment etc.
    Grid_3d lfc;

    //get the grid dimensions
    lfc.get_dim(p.box_x,p.box_y,p.box_z,traj.ibox,p.APS);
    if(s.world_rank == 0) //attempt to estimate requirements befor a crash
    {
        float mem_d = 4.0*(float)lfc.num_x()*(float)lfc.num_y()*(float)lfc.num_z()*8.0;
        float mem_i = 2.0*(float)lfc.num_x()*(float)lfc.num_y()*(float)lfc.num_z()*4.0;
        float mem_t = 1.0*(mem_d + mem_i);
        printf("Estimated memory to hold the grid: %f (MB) \n\n",mem_t/1000000.0);
    }

    //set the output file name for grid
    lfc.set_output(p.lfc_file_name);

    //print info about the grid
    lfc.print_dim();

    //make vector to hold distances
    iv1d contacts(0,0);                               //stores each inter leaflet contacts measurements
    dv1d contacts_frame(0,0.0);                       //stores average number of contacts for each frame
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

        get_contacts(traj,s,p,param_1,param_2,lfc,contacts,contacts_frame);

        traj.set_beta_lf();

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect the contacts from each mpi process and compute the average
    perf.log_time(finalize_analysis(traj,s,p,lfc,contacts,contacts_frame),"Fin Ana");

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

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
#include "headers/mosat_routines.h"                          //This is where most of the functions called in main are located
#include "headers/file_naming.h"                             //This has routines for added tags to an existing file name    
#include "headers/file_naming_mpi.h"                         //This has routines for added tags to an existing file name (mpi)
#include "headers/command_line_args_mpi.h"                   //This has routines for adding command line arguments
#include "MosAT/program_variables/pv_contact_analysis.h"     //This has the variables specific to the analysis program
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
#include "headers/contacts_full.h"                           //This has routines used for reading and writing contacts matrices data files

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function counts the contacts between the 2 groups of atoms.                                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void contacts(Trajectory &traj,system_variables &s,program_variables &p,Index &group_1,Index &group_2,iv1d &cont,Contacts &cont_mat)
{
    int i        = 0;                                //standard variable used in loops
    int j        = 0;                                //standard variable used in loops
    int contacts = 0;                                //keep track of the number of contacts

    iv2d this_profile(0,iv1d(0,0));                  //hold contact info for the current frame

    for(i=0; i<group_1.index_i.size(); i++) //loop over the first group of atoms
    {
        for(j=0; j<group_2.index_i.size(); j++) //loop over the second group of atoms
        {
            int index_1 = group_1.index_i[i]-1;
            int index_2 = group_2.index_i[j]-1;

            double dx = traj.r[index_1][0] - traj.r[index_2][0];
            double dy = traj.r[index_1][1] - traj.r[index_2][1];
            double dz = traj.r[index_1][2] - traj.r[index_2][2];

	    double distance = sqrt(dx*dx + dy*dy + dz*dz); 

	    if(distance <= p.cdist)
            {
                iv1d this_contact(2,0);      //holds the positions int the indices for the current contact 
                this_contact[0] = i;
                this_contact[1] = j;
		this_profile.push_back(this_contact);

                contacts++;
            } 
	}	
    }
    cont[traj.current_frame] = contacts;

    //add the current profile to tmp files
    cont_mat.add_profile(this_profile,traj.get_frame_global());
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Collect contacts and write data to output file.                                                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,iv1d &cont,Contacts &cont_mat,Index &group_1,Index &group_2)
{
    int i = 0; //standard variable used in loops
    int j = 0; //standard variable used in loops 

    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nfinalizing analysis. This requires communication of data and could take some time. \n");
    }

    if(s.world_rank == 0)
    {
        double mem = (double)p.size_x*(double)p.size_y*4.0/1000000;
        printf("Memory estimate: %f MB \n",mem);
    }

    //merge tmp contact matrices to a single file
    cont_mat.merge_profiles();

    //here we collect the contact counts from each core
    collect_iv1d(s.world_size,s.world_rank,cont);

    if(s.world_rank == 0)
    {
        FILE *cont_file = fopen(p.cont_file_name.c_str(),"w");

        for(i=0; i<cont.size(); i++) //loop over frames
        {
            fprintf(cont_file," %10d %10d \n",i,cont[i]);
        }
        fclose(cont_file);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Get the frequency of each contact                                                                         //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        iv2d this_cont_mat(p.size_x,iv1d(p.size_y,0));
        for(i=0; i<traj.get_ef_frames(); i++) //loop over the frames
        {
	    cont_mat.get_profile_stamp(i,this_cont_mat);
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Create a list of contacts so we can arrange them with largest freq first                                  //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        FILE *contacts_file;
        string contacts_file_name = chop_and_add_tag(p.cont_file_name,"_contacts_freq.pml");
        contacts_file      = fopen(contacts_file_name.c_str(), "w");
        if(contacts_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",contacts_file_name.c_str());
        }
        else
        {
            iv2d contacts_list(0,iv1d(2,0));      //holds the pair of atom indexes for each contact
            dv1d contacts_list_freq(0,0.0);       //holds the frequency of each contact

            for(i=0; i<p.size_y; i++) //loop over y
            {
                for(j=0; j<p.size_x; j++) //loop over x
                {
                    double this_freq = (double)this_cont_mat[j][i]/(double)traj.get_ef_frames(); 
         
                    if(this_freq > 0.0)
                    {
                        iv1d this_cont(2,0);

                        this_cont[0] = j;
                        this_cont[1] = i;
                        contacts_list.push_back(this_cont);
                        contacts_list_freq.push_back(this_freq);
                    }
                }
            }

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                                                           //
            // Organize by largest frequency                                                                             //
            //                                                                                                           //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            if(contacts_list_freq.size() > 0)
            {
                int b_swap=0;
                for(b_swap=1; b_swap > 0; )
                {
                    b_swap = 0;
                    for(i=0; i<contacts_list_freq.size()-1; i++)
                    {
                        if(contacts_list_freq[i] < contacts_list_freq[i+1])
                        {
                            double contacts_list_freq_tmp  = contacts_list_freq[i];
                            int contacts_list_x_tmp        = contacts_list[i][0];
                            int contacts_list_y_tmp        = contacts_list[i][1];

                            contacts_list_freq[i] = contacts_list_freq[i+1];
                            contacts_list[i][0]   = contacts_list[i+1][0];
                            contacts_list[i][1]   = contacts_list[i+1][1];

                            contacts_list_freq[i+1] = contacts_list_freq_tmp;
                            contacts_list[i+1][0]   = contacts_list_x_tmp;
                            contacts_list[i+1][1]   = contacts_list_y_tmp;

                            b_swap = 1;
                        }
                    }
                }
            }

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                                                           //
            // Write PyMol select commands for identifying relevant contacts based on frequency                          //
            //                                                                                                           //
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
            for(i=0; i<contacts_list_freq.size(); i++)
            {
                int    atom_number_1 = group_1.index_i[contacts_list[i][0]];
                int    atom_number_2 = group_2.index_i[contacts_list[i][1]];
                double freq          = contacts_list_freq[i];
                double dash_rad      = freq*(p.max_dash_rad - p.min_dash_rad) + p.min_dash_rad;
                string dist_name     = "dist_" + to_string(i); 

                if(dash_rad >= 0.0000005) //printf %f gives 6 decimal places. anything above 0.0000005 rounds up to 0.000001 giving a non zero radius. PyMOL shows a zero radius too big (maybe a default value)
                {
                    fprintf(contacts_file,"\nContact %d: Frequency %f: Atom_nr b %d: Atom_nr a %d \n",i,freq,atom_number_1,atom_number_2);
                    fprintf(contacts_file,"select atom_a, (resi %d & resn %s & name %s) \n",traj.res_nr[atom_number_1-1]%10000,traj.res_name[atom_number_1-1].c_str(),traj.atom_name[atom_number_1-1].c_str());
                    fprintf(contacts_file,"select atom_b, (resi %d & resn %s & name %s) \n",traj.res_nr[atom_number_2-1]%10000,traj.res_name[atom_number_2-1].c_str(),traj.atom_name[atom_number_2-1].c_str());
                    fprintf(contacts_file,"select %d, (resi %d & resn %s & name %s) or ((resi %d & resn %s & name %s)) \n",i,traj.res_nr[atom_number_1-1]%10000,traj.res_name[atom_number_1-1].c_str(),traj.atom_name[atom_number_1-1].c_str(),traj.res_nr[atom_number_2-1]%10000,traj.res_name[atom_number_2-1].c_str(),traj.atom_name[atom_number_2-1].c_str());
                    fprintf(contacts_file,"dist %s, (atom_a), (atom_b) \n",dist_name.c_str());
                    fprintf(contacts_file,"set dash_radius, %f, %s \n",dash_rad,dist_name.c_str());
                }
            }
            fclose(contacts_file);
        }
    }

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
    s.program_name = "Contact Analysis";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                  s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (tpr, gro)",                                  s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                 s.world_rank, s.cl_tags, &p.b_print,   0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                   s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                              s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                               s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting",                                       s.world_rank, s.cl_tags, &p.b_lsq,     0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                   s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",   s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-n1",     p.group_1_file_name,          "Index file with group 1 atoms (ndx)",                         s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-n2",     p.group_2_file_name,          "Index file with group 2 atoms (ndx)",                         s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-cont",   p.cont_file_name,             "Output file name (dat) with contacts",                        s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-cdist",  &p.cdist,                     "Distance threshold (nm) for counting contacts",               s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-max",    &p.max_dash_rad,              "Maximum thickness of dash for PyMOL distance commands ",      s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-min",    &p.min_dash_rad,              "Minimum thickness of dash for PyMOL distance commands ",      s.world_rank, s.cl_tags, nullptr,      0);
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
    check_extension_mpi(s.world_rank,"-n1",p.group_1_file_name,".ndx");
    check_extension_mpi(s.world_rank,"-n2",p.group_1_file_name,".ndx");
    check_extension_mpi(s.world_rank,"-cont",p.cont_file_name,".dat");

    //create index for each group of atoms
    Index group_1;
    Index group_2;

    //read the index files
    group_1.get_index(p.group_1_file_name);
    group_2.get_index(p.group_2_file_name);

    //create vector to hold contacts
    iv1d cont(traj.get_num_frames(),0);

    //store the contact matrix dimensions
    p.size_x = group_1.index_i.size();
    p.size_y = group_2.index_i.size();

    //create object for working with contact matrices
    Contacts cont_mat; 

    //initialize contacts
    cont_mat.init(p.cont_file_name,p.size_x,p.size_y);

    //open files for writing temporary contact profiles data
    cont_mat.prime_tmp_file();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //print info about the worlk load distribution
    traj.workload();

    //print that analysis is beginning
    traj.report_progress();

    s.t = clock();
    //read read frames of the trajectory and perform analysis
    for(traj.current_frame=0; traj.current_frame<traj.get_num_frames(); traj.current_frame++)
    {
        traj.read_traj_frame();

        traj.do_fit();

        contacts(traj,s,p,group_1,group_2,cont,cont_mat);

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect the distances from each mpi process and write data to file
    perf.log_time(finalize_analysis(traj,s,p,cont,cont_mat,group_1,group_2),"Fin Ana");

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

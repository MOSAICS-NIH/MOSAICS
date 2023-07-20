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

#include "xdr/include/xdrfile_xtc.h"                            //used to read xtc files 
#include "xdr/include/xdr_seek.h"                               //used to get and set the file position in xtc and trr files
#include "xdr/include/xdrfile_trr.h"                            //used to read trr files
#include "xdr/include/xdrfile.h"                                //used to read xtc and trr files
#include "xdr/include/trr_header.h"                             //used to read the header info of trr files
#include "headers/multi_dim_vec.h"                              //This defines multidimensional vectors
#include "headers/switch.h"                                     //This defines a switch (on, off)
#include "headers/file_reader.h"                                //This has basic routines for reading text files
#include "headers/vector_mpi.h"                                 //This has routines for collecting vector data
#include "headers/mosat_routines.h"                             //This is where most of the functions called in main are located
#include "headers/file_naming.h"                                //This has routines for added tags to an existing file name    
#include "headers/file_naming_mpi.h"                            //This has routines for added tags to an existing file name (mpi)
#include "headers/command_line_args_mpi.h"                      //This has routines for adding command line arguments
#include "MosAT/program_variables/pv_surface_residue_finder.h"  //This has the variables specific to the analysis program
#include "headers/performance.h"                                //This has a class for logging performance data
#include "headers/array.h"                                      //This has routines used for working with arrays
#include "headers/index.h"                                      //This has a class for working with index files
#include "headers/traj.h"                                       //This has a class for working with the trajectory
#include "headers/leaflet_finder.h"                             //This has routines used to find leaflets in membrane simulations
#include "headers/protein_finder.h"                             //This has routines used to find protein atoms
#include "headers/sol_finder.h"                                 //This has routines used to find the solvent
#include "headers/grid.h"                                       //This has routines used for working with a grid
#include "headers/protein.h"                                    //This has routines used for working with protein data
#include "headers/force_serial.h"                               //This has routines used for forcing the code to run on a single mpi process
#include "headers/param.h"                                      //This has routines used for reading complex parameter data
#include "headers/histo.h"                                      //This has routines used for binning data and making a histogram
#include "headers/atom_select.h"                                //This has routines used for making atom selections using a selection text

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function counts how many lipids make contact with each protein residue.                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void protein_surface(Trajectory &traj,system_variables &s,program_variables &p,Param &param,int contacts_sum[],
                     iv1d &surface_count,iv1d &refined_sel)
{
    int i          = 0;                               //loop variable
    int j          = 0;                               //loop variable
    int k          = 0;                               //loop variable
    int l          = 0;                               //loop variable
    int m          = 0;                               //loop variable
    int n          = 0;                               //loop variable
    int prot_count = 0;

    //allocate memory to hold contact information for current frame
    iv1d contact(traj.prot.size(),0);  //number of contacts between the atom and lipids

    //count contacts with each atom
    for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
    {
        int min_p = traj.p_res_start(i);
        int max_p = traj.p_res_end(i);

        i = traj.next_prot_res(i);

        for(j=min_p; j<=max_p; j++) //loop over current residue atoms
        {
            if(refined_sel[j] == 1) //check refined selection
            {
                for(k=0; k<traj.target_leaflet.size(); k++) //loop over target leaflet atoms
                {
                    //get the first and last atom of the current lipid
                    int min = traj.t_lip_start(k);
                    int max = traj.t_lip_end(k);

                    //jump to the next lipid
                    k = traj.next_target_lipid(k);

                    //count contacts for target lipids 
                    for(l=0; l<param.main_size_y(); l++) //loop over lipid types 
                    {
                        if(strcmp(traj.res_name[min].c_str(), param.param_main_s[l][0].c_str()) == 0) //residue is a target lipid
                        {
                            for(m=min; m<=max; m++) //loop over current lipid atoms
                            {
                                for(n=0; n<param.sec_size_y(l); n++) //loop over target atoms
                                {
                                    if(strcmp(traj.atom_name[m].c_str(), param.param_sec_s[l][n][0].c_str()) == 0) //atom is a target atom
                                    {
                                        double dif_x = traj.r[j][0] - traj.r[m][0];
                                        double dif_y = traj.r[j][1] - traj.r[m][1];
                                        double dif_z = traj.r[j][2] - traj.r[m][2];

                                        double dist = sqrt(dif_x*dif_x + dif_y*dif_y + dif_z*dif_z);

                                        if(dist < p.contact_cutoff)
                                        {
                                            contact[prot_count] = 1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            prot_count++;
        }
    }

    //add contacts to long term sums and count contacts for the curren frame
    int count = 0;
    for(i=0; i<traj.prot.size(); i++)
    {
        contacts_sum[i] = contacts_sum[i] + contact[i];

        if(contact[i] == 1)
        {
            count++;
        }
    }
    surface_count.push_back(count);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects contact data from the ranks and prints the final surface results to a text and     //
// pdb file.                                                                                                 //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,int contacts_sum[],iv1d &surface_count)  
{
    int i = 0;                                       //standard variable used in loops
    int j = 0;                                       //standard variable used in loops
    int k = 0;                                       //standard variable used in loops

    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nFinalizing analysis. This could take some time. \n");
    }

    //collect contacts_sum from all cores
    collect_protein_res_data_i(s.world_rank,s.world_size,traj.prot.size(),contacts_sum);

    //collect surface count
    collect_iv1d(s.world_size,s.world_rank,surface_count);

    if(s.world_rank == 0)
    {
        dv1d percent(traj.prot.size(),0.0); 

        for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
        {
            percent[i] = 100.0*(double)contacts_sum[i]/(double)traj.get_ef_frames();
        }

        //set beta factor to surface
        for(i=0; i<traj.atoms(); i++) //loop over system atoms
        {
            traj.beta[i] = -1;

            for(j=0; j<traj.prot.size(); j++) //loop over protein atoms
            {
                if(i == traj.prot[j]-1) //atom is a protein atom
                {
                    traj.beta[traj.prot[j]-1] = percent[j];
                }
            }
        }

        //write scored atoms to pdb
        string ranks_file_name = add_tag(p.surface_pdb_file_name,"_ranks");
        write_data_pdb(traj.box_ref,traj.atoms(),traj.atom_nr,traj.res_nr,traj.res_name,traj.atom_name,traj.r_ref,traj.time,traj.step,traj.beta,traj.weight,traj.chain_id,traj.element,
                       traj.title,s.world_rank,ranks_file_name);

        //compute average surface atoms
        double avg = 0.0;
        for(i=0; i<surface_count.size(); i++) //loop over frames
        {
            avg = avg + (double)surface_count[i];
        }
        avg = avg/(double)surface_count.size();

        //compute stdev
        double stdev = 0.0;
        for(i=0; i<surface_count.size(); i++) //loop over frames
        {
            stdev = stdev + pow(((double)surface_count[i] - avg),2);
        }
        stdev = stdev/((double)surface_count.size()-1.0);
        stdev = sqrt(stdev);
 
        printf("There are %f +/- %f surface atoms \n",avg,stdev);

        //set beta factor to 0
        for(i=0; i<traj.atoms(); i++) //loop over system atoms
        {
            traj.beta[i] = 0;
        }

        //select surface atoms
        int target = (int)floor(avg + p.cutoff*stdev); 
        int largest_j  = 0;
        for(i=0; i<target; i++) //loop over num surface atoms
        {
            double largest = 0.0;
            for(j=0; j<traj.prot.size(); j++) //loop over protein atoms
            {
                if(percent[j] > largest)
                {
                    largest   = percent[j];
                    largest_j = j;
                }
            }

            traj.beta[traj.prot[largest_j]-1] = 1;
            percent[largest_j]                = -1 ;
        }

        //write surface to pdb
        write_data_pdb(traj.box_ref,traj.atoms(),traj.atom_nr,traj.res_nr,traj.res_name,traj.atom_name,traj.r_ref,traj.time,traj.step,traj.beta,traj.weight,traj.chain_id,traj.element,
                       traj.title,s.world_rank,p.surface_pdb_file_name);

        //make histogram of the number of surface atoms
        Histogram_i histo;
        histo.bin_data(surface_count,1);
        histo.write_histo(p.histo_file_name,"surface atoms count"); 

        //print select commands for surface atoms
        printf("\n");
        printf("select surface, id ");
        int count = 0;
        for(i=0; i<traj.atoms(); i++) //loop over system atoms
        { 
            if(traj.beta[i] == 1) //atom is a surface atom
            {
                if(count > 0)
                {
                    printf("+");      
                }
                printf("%d",i+1);
                count++;
            }
        }
        printf("\n\n");

        //print select commands for core atoms
        printf("select core, id ");
        count = 0;
        for(i=0; i<traj.atoms(); i++) //loop over system atoms
        {
            for(j=0; j<traj.prot.size(); j++) //loop over protein atoms
            {
                if(traj.prot[j]-1 == i && traj.beta[i] == 0) //atom is a core atom
                {
                    if(count > 0)
                    {
                        printf("+");
                    }
                    printf("%d",i+1);
                    count++;
                }
            }
        }
        printf("\n\n");
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
    s.program_name = "Surface Residue Finder";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                    s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                    s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                   s.world_rank, s.cl_tags, &p.b_print,   0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                     s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                                s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                 s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                   s.world_rank, s.cl_tags, &p.b_lsq,     0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                     s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",     s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-crd",    p.param_file_name,            "Selection card with lipid types and target atoms (crd)",        s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-histo",  p.histo_file_name,            "Output file with surface atoms histogram (dat)",                s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-bin",    &p.lsq_ref,                   "Bin width for the surface atoms histogram (atoms, int)",        s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-s_pdb",  p.surface_pdb_file_name,      "Output file with surface residues indicated by B-factor (pdb)", s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,           "PDB file with sorted leaflets (B-factor) (pdb)",                s.world_rank, s.cl_tags, &p.b_lf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,  "File with additional leaflet finder parameters (prm)",          s.world_rank, s.cl_tags, &p.b_lf_param,0);
    add_argument_mpi_s(argc,argv,"-pf_pdb", p.pf_pdb_file_name,           "PDB file with selected protein (B-factor) (pdb)",               s.world_rank, s.cl_tags, &p.b_pf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-pf_prm", p.protein_finder_param_name,  "File with additional protein finder parameters (prm)",          s.world_rank, s.cl_tags, &p.b_pf_param,0);
    add_argument_mpi_d(argc,argv,"-cutoff", &p.cutoff,                    "Cutoff for selecting the N surface atoms (chi)",                s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-cdist",  &p.contact_cutoff,            "The contact cutoff distance (nm)",                              s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-sel",    p.selection_text_file_name,   "Input file with the atom selection text (sel)",                 s.world_rank, s.cl_tags, &p.b_sel_text,0);
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
    check_extension_mpi(s.world_rank,"-histo",p.histo_file_name,".dat");
    check_extension_mpi(s.world_rank,"-s_pdb",p.surface_pdb_file_name,".pdb");

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
    if(p.b_sel_text == 1)
    {
        check_extension_mpi(s.world_rank,"-sel",p.selection_text_file_name,".sel");
    }

    //create parameter files
    Param param;

    //read parameter files
    param.get_param(p.param_file_name,2,1,1);

    //check the integrity of the parameter files
    if(param.check_file() == 0) //bad files. kill program 
    {
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }

    //run leaflet/proten/solvent finder
    traj.get_leaflets(p.leaflet,p.leaflet_finder_param_name,p.b_lf_param);
    traj.get_protein(p.protein_finder_param_name,p.b_pf_param);

    //print a pdb with distinguished leaflets/protein/solvent
    traj.write_leaflets(p.lf_pdb_file_name,p.b_lf_pdb);
    traj.write_protein(p.pf_pdb_file_name,p.b_pf_pdb);

    //print info about the target leaflet
    traj.get_leaflet_stats();

    //print info about the lipid selection
    traj.get_lipid_selection_stats(param.get_column_s(0),"-crd");

    //print info about the protein
    traj.get_prot_stats();

    //allocate memory for contacts 
    int contacts_sum[traj.prot.size()];           //Stores the contact counts
    init_iarray(contacts_sum,traj.prot.size());

    iv1d surface_count(0,0);                      //Holds the number of surface atoms for each frame

    //create object to hold protein atom refinement
    iv1d refined_sel(traj.atoms(),1);

    if(p.b_sel_text == 1)
    {
        //create a object to hold an atom selection
        Selection this_sel;

        //select the atoms
        this_sel.get_selection(traj,p.selection_text_file_name);

        //generate pdb file name for highlighting the selection 
        string pdb_filename = chop_and_add_tag(p.selection_text_file_name,".pdb");

        //highlight the selected atoms
        this_sel.highlight_sel(traj,pdb_filename);

        //refine the selection
        refined_sel = this_sel.tag_atoms(traj);
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

        protein_surface(traj,s,p,param,contacts_sum,surface_count,refined_sel);

        traj.set_beta_lf();

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect data from mpi processes and compute averages
    perf.log_time(finalize_analysis(traj,s,p,contacts_sum,surface_count),"Fin Ana");

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

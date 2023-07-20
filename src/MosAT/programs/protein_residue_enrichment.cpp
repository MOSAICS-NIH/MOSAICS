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

#include "xdr/include/xdrfile_xtc.h"                                //used to read xtc files 
#include "xdr/include/xdr_seek.h"                                   //used to get and set the file position in xtc and trr files
#include "xdr/include/xdrfile_trr.h"                                //used to read trr files
#include "xdr/include/xdrfile.h"                                    //used to read xtc and trr files
#include "xdr/include/trr_header.h"                                 //used to read the header info of trr files
#include "headers/multi_dim_vec.h"                                  //This defines multidimensional vectors
#include "headers/switch.h"                                         //This defines a switch (on, off)
#include "headers/file_reader.h"                                    //This has basic routines for reading text files
#include "headers/vector_mpi.h"                                     //This has routines for collecting vector data
#include "headers/mosat_routines.h"                                //This is where most of the functions called in main are located
#include "headers/file_naming.h"                                    //This has routines for added tags to an existing file name    
#include "headers/file_naming_mpi.h"                                //This has routines for added tags to an existing file name (mpi)
#include "headers/command_line_args_mpi.h"                          //This has routines for adding command line arguments
#include "MosAT/program_variables/pv_protein_residue_enrichment.h" //This has the variables specific to the analysis program
#include "headers/array.h"                                          //This has routines used for working with arrays
#include "headers/performance.h"                                    //This has a class for logging performance data
#include "headers/index.h"                                          //This has a class for working with index files
#include "headers/traj.h"                                           //This has a class for working with the trajectory
#include "headers/leaflet_finder.h"                                 //This has routines used to find leaflets in membrane simulations
#include "headers/protein_finder.h"                                 //This has routines used to find protein atoms
#include "headers/sol_finder.h"                                     //This has routines used to find the solvent
#include "headers/grid.h"                                           //This has routines used for working with a grid
#include "headers/protein.h"                                        //This has routines used for working with protein data
#include "headers/force_serial.h"                                   //This has routines used for forcing the code to run on a single mpi process
#include "headers/param.h"                                          //This has routines used for reading complex parameter data

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function counts how many lipids make contact with each protein residue.                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void protein_enrichment(Trajectory &traj,system_variables &s,program_variables &p,Param &param_1,Param &param_2,int contacts_sum_1[],int contacts_sum_2[])
{
    int i         = 0;                                //loop variable
    int j         = 0;                                //loop variable
    int k         = 0;                                //loop variable
    int l         = 0;                                //loop variable
    int m         = 0;                                //loop variable
    int n         = 0;                                //loop variable
    int count_r   = 0;                                //Count the residues as we find them
    int count_l_1 = 0;                                //Count the lipids as we find them
    int count_l_2 = 0;                                //Count the lipids as we find them

    //allocate memory to hold contact information for current frame
    iv2d contact_1(p.prot_res_count,iv1d(p.lip_count_1,0));  //number of contacts between the residue and lipids 1
    iv2d contact_2(p.prot_res_count,iv1d(p.lip_count_2,0));  //number of contacts between the residue and lipids 2

    //count contacts with each lipid type
    for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
    {
        int min_p = traj.p_res_start(i);  
        int max_p = traj.p_res_end(i);

        i = traj.next_prot_res(i);

        for(j=min_p; j<=max_p; j++) //loop over current residue atoms
        {
            //reset lipid count
            count_l_1 = 0;
            count_l_2 = 0;

            for(k=0; k<traj.target_leaflet.size(); k++) //loop over target leaflet atoms
            {
                //get the first and last atom of the current lipid
                int min = traj.t_lip_start(k);
                int max = traj.t_lip_end(k);

                //jump to the next lipid
                k = traj.next_target_lipid(k);

                //count contacts for lipid types 1
                for(l=0; l<param_1.main_size_y(); l++) //loop over lipid types 1
                {
                    if(strcmp(traj.res_name[min].c_str(), param_1.param_main_s[l][0].c_str()) == 0) //residue is a target lipid
                    {
                        for(m=min; m<=max; m++) //loop over current lipid atoms
                        {
                            for(n=0; n<param_1.sec_size_y(l); n++) //loop over target atoms
                            {
                                if(strcmp(traj.atom_name[m].c_str(), param_1.param_sec_s[l][n][0].c_str()) == 0) //atom is a target atom
                                {
                                    double dif_x = traj.r[j][0] - traj.r[m][0];
                                    double dif_y = traj.r[j][1] - traj.r[m][1];
                                    double dif_z = traj.r[j][2] - traj.r[m][2];

                                    double dist = sqrt(dif_x*dif_x + dif_y*dif_y + dif_z*dif_z);

                                    if(dist < p.contact_cutoff)
                                    {
                                        contact_1[count_r][count_l_1] = 1;
                                    }
                                }
                            }
                        }
                        count_l_1++;
                    }
                }

                //count contacts for lipid types 2
                for(l=0; l<param_2.main_size_y(); l++) //loop over lipid types 2
                {
                    if(strcmp(traj.res_name[min].c_str(), param_2.param_main_s[l][0].c_str()) == 0) //residue is a target lipid
                    {
                        for(m=min; m<=max; m++) //loop over current lipid atoms
                        {
                            for(n=0; n<param_2.sec_size_y(l); n++) //loop over target atoms
                            {
                                if(strcmp(traj.atom_name[m].c_str(), param_2.param_sec_s[l][n][0].c_str()) == 0) //atom is a target atom
                                {
                                    double dif_x = traj.r[j][0] - traj.r[m][0];
                                    double dif_y = traj.r[j][1] - traj.r[m][1];
                                    double dif_z = traj.r[j][2] - traj.r[m][2];

                                    double dist = sqrt(dif_x*dif_x + dif_y*dif_y + dif_z*dif_z);

                                    if(dist < p.contact_cutoff)
                                    {
                                        contact_2[count_r][count_l_2] = 1;
                                    }
                                }
                            }
                        }
                        count_l_2++;
                    }
                }
            }
        }
        count_r++;
    }

    //add contacts to long term sums
    for(i=0; i<p.prot_res_count; i++)
    {
        for(j=0; j<p.lip_count_1; j++)
        {
            contacts_sum_1[i] = contacts_sum_1[i] + contact_1[i][j];
        }
        for(j=0; j<p.lip_count_2; j++)
        {
            contacts_sum_2[i] = contacts_sum_2[i] + contact_2[i][j];
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects contact data from the ranks and prints the final enrichment results to a text and  //
// pdb file.                                                                                                 //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,int contacts_sum_1[],int contacts_sum_2[])
{
    int i       = 0;                                 //standard variable used in loops
    int j       = 0;                                 //standard variable used in loops
    int k       = 0;                                 //standard variable used in loops
    int count_r = 0;                                 //keeps track of residues as they are encountered
    double enrichment[p.prot_res_count];             //Enrichment for each residue
    init_darray(enrichment,p.prot_res_count);

    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nFinalizing analysis. This could take some time. \n");
    }

    //collect contact counts from mpi processes
    collect_protein_res_data_i(s.world_rank,s.world_size,p.prot_res_count,contacts_sum_1);
    collect_protein_res_data_i(s.world_rank,s.world_size,p.prot_res_count,contacts_sum_2);

    //compute percent enrichment and write data to ouput file
    if(s.world_rank == 0)
    {
        //open output file for writing data
        FILE *enrich_file = fopen(p.enrich_file_name.c_str(), "w");
        if(enrich_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",p.enrich_file_name.c_str());
        }
        else
        {
            fprintf(enrich_file,"%15s %15s %15s %15s %15s \n","Residue","Contacts_A","Contacts_B","Ratio","Enrichment (%)");
            for(i=0; i<p.prot_res_count; i++) //loop over residues
            {
                if(contacts_sum_1[i] == 0 || contacts_sum_2[i] == 0)
                {
                    fprintf(enrich_file,"%15d %15d %15d %15s %15s \n",i,contacts_sum_1[i],contacts_sum_2[i],"NaN","NaN");
                }
                else
                {
                    double ratio  = (double)contacts_sum_1[i]/(double)contacts_sum_2[i];
                    double enrich = 100*(ratio - p.bulk)/p.bulk;
                    enrichment[i] = enrich;
                    fprintf(enrich_file,"%15d %15d %15d %15.5f %15.5f \n",i,contacts_sum_1[i],contacts_sum_2[i],ratio,enrich);
                }
            }
            fclose(enrich_file);
        }

        //initiallize non protein atoms
        for(i=0; i<traj.atoms(); i++) //loop over system atoms
        {
            traj.beta[i] = 0;
        }

        //set beta factor to enrichment
        for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
        {
            int min_p = traj.p_res_start(i);
            int max_p = traj.p_res_end(i);

            i = traj.next_prot_res(i);

            if(contacts_sum_1[count_r] + contacts_sum_2[count_r] > p.cutoff*traj.get_ef_frames()) //residue is sufficiently solvated
            {
                for(j=min_p; j<=max_p; j++) //loop over curent res atoms
                {
                    traj.beta[j] = enrichment[count_r];
                }
            }
            else //exlude data
            {
                traj.beta[j] = 0;
            }
            count_r++;
        }
        count_r = 0;

        //write enrichment to pdb
        write_data_pdb(traj.box_ref,traj.atoms(),traj.atom_nr,traj.res_nr,traj.res_name,traj.atom_name,traj.r_ref,traj.time,traj.step,traj.beta,traj.weight,traj.chain_id,traj.element,
                       traj.title,s.world_rank,p.enrich_pdb_file_name);

        //create name of density files 
        string density_1_file_pdb_name = add_tag(p.enrich_pdb_file_name,"_A");
        string density_2_file_pdb_name = add_tag(p.enrich_pdb_file_name,"_B");

        
        //set beta factor to density 1
        for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
        {
            int min_p = traj.p_res_start(i);
            int max_p = traj.p_res_end(i);

            i = traj.next_prot_res(i);

            for(j=min_p; j<=max_p; j++) //loop over curent res atoms
            {
                traj.beta[j] = contacts_sum_1[count_r]/(double)traj.get_ef_frames();
            }
            count_r++;
        }
        count_r = 0;

        //write density_1 to pdb
        write_data_pdb(traj.box_ref,traj.atoms(),traj.atom_nr,traj.res_nr,traj.res_name,traj.atom_name,traj.r_ref,traj.time,traj.step,traj.beta,traj.weight,traj.chain_id,traj.element,
                       traj.title,s.world_rank,density_1_file_pdb_name);

        //set beta factor to density 2
        for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
        {
            int min_p = traj.p_res_start(i);
            int max_p = traj.p_res_end(i);

            i = traj.next_prot_res(i);

            for(j=min_p; j<=max_p; j++) //loop over curent res atoms
            {
                traj.beta[j] = contacts_sum_2[count_r]/(double)traj.get_ef_frames();
            }
            count_r++;
        }
        count_r = 0;

        //write density_2 to pdb
        write_data_pdb(traj.box_ref,traj.atoms(),traj.atom_nr,traj.res_nr,traj.res_name,traj.atom_name,traj.r_ref,traj.time,traj.step,traj.beta,traj.weight,traj.chain_id,traj.element,
                       traj.title,s.world_rank,density_2_file_pdb_name);
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
    s.program_name = "Protein Residue Enrichment";

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
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                                  s.world_rank, s.cl_tags, &p.b_print,   0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                                    s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                                               s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                                s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                                  s.world_rank, s.cl_tags, &p.b_lsq,     0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                                    s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",                    s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-crd_1",  p.param_1_file_name,          "Selection card for lipids A (crd)",                                            s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-crd_2",  p.param_2_file_name,          "Selection card for lipids B (crd)",                                            s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-enrich", p.enrich_file_name,           "Output file with percent enrichment for each protein residue (dat)",           s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-e_pdb",  p.enrich_pdb_file_name,       "Output file with percent enrichment projected on the protein structure (pdb)", s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-leaf",   &p.leaflet,                   "Which leaflet? (0:both 1:upper 2:lower)",                                      s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,           "PDB file with sorted leaflets (B-factor) (pdb)",                               s.world_rank, s.cl_tags, &p.b_lf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,  "File with additional leaflet finder parameters (prm)",                         s.world_rank, s.cl_tags, &p.b_lf_param,0);
    add_argument_mpi_s(argc,argv,"-pf_pdb", p.pf_pdb_file_name,           "PDB file with selected protein (B-factor) (pdb)",                              s.world_rank, s.cl_tags, &p.b_pf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-pf_prm", p.protein_finder_param_name,  "File with additional protein finder parameters (prm)",                         s.world_rank, s.cl_tags, &p.b_pf_param,0);
    add_argument_mpi_d(argc,argv,"-cutoff", &p.cutoff,                    "Cutoff for excluding data (chi)",                                              s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-cdist",  &p.contact_cutoff,            "The contact cutoff distance (nm)",                                             s.world_rank, s.cl_tags, nullptr,      1);
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
    check_extension_mpi(s.world_rank,"-crd_1",p.param_1_file_name,".crd");
    check_extension_mpi(s.world_rank,"-crd_2",p.param_2_file_name,".crd");
    check_extension_mpi(s.world_rank,"-enrich",p.enrich_file_name,".dat");
    check_extension_mpi(s.world_rank,"-e_pdb",p.enrich_pdb_file_name,".pdb");

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

    //create parameter files
    Param param_1;
    Param param_2;

    //read parameter files
    param_1.get_param(p.param_1_file_name,2,1,1);
    param_2.get_param(p.param_2_file_name,2,1,1);

    //check the integrity of the parameter files
    if(param_1.check_file() == 0 || param_2.check_file() == 0) //bad files. kill program 
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
    traj.get_lipid_selection_stats(param_1.get_column_s(0),"-crd_1");
    traj.get_lipid_selection_stats(param_2.get_column_s(0),"-crd_2");

    //print info about the protein
    traj.get_prot_stats();

    //count the number of protein residues and the number of lipids for each type
    p.prot_res_count = traj.get_num_res_prot();
    p.lip_count_1    = traj.count_target_lipids_type(param_1.get_column_s(0));
    p.lip_count_2    = traj.count_target_lipids_type(param_2.get_column_s(0));
    p.bulk = (double)p.lip_count_1/(double)p.lip_count_2;

    //allocate memory for the protein
    int contacts_sum_1[p.prot_res_count];           //Stores the density for lipid type 1
    int contacts_sum_2[p.prot_res_count];           //Stores the density for lipid type 2
    init_iarray(contacts_sum_1,p.prot_res_count);
    init_iarray(contacts_sum_2,p.prot_res_count);
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

        protein_enrichment(traj,s,p,param_1,param_2,contacts_sum_1,contacts_sum_2);

        traj.set_beta_lf();

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect data from mpi processes and compute average
    perf.log_time(finalize_analysis(traj,s,p,contacts_sum_1,contacts_sum_2),"Fin Ana");

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

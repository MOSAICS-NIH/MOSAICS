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
#include "MosAT/program_variables/pv_lipid_h_bonds.h"       //This has the variables specific to the analysis program
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
// This function computes the number of lipid-prot h-bonds and adds it to the grid                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int check_h_bond(Trajectory &traj,program_variables &p,int acceptor,int donor,int h)
{
    int result = 0;
    int pi     = 3.1415926535;
    int i      = 0;

    rvec m;
    rvec n;

    for(i=0; i<3; i++) //loop over 3 dimensions
    {
        m[i] = traj.r[acceptor][i] - traj.r[donor][i];
        n[i] = traj.r[h][i]        - traj.r[donor][i];
    }

    double angle = gmx_angle(m,n);
    double dist  = sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);

    //convert angle to degrees
    angle = angle*180/pi;
          
    if(dist < 0.35 && angle < 30)
    {
        result = 1;

        //check h-bonds in pymol
        if(p.b_test == 1)
        {
            printf("acceptor %10d %10s %10d %10s \n",traj.atom_nr[acceptor],traj.atom_name[acceptor].c_str(),traj.res_nr[acceptor],traj.res_name[acceptor].c_str());
            printf("donor    %10d %10s %10d %10s \n",traj.atom_nr[donor],   traj.atom_name[donor].c_str(),   traj.res_nr[donor],   traj.res_name[donor].c_str());
            printf("h        %10d %10s %10d %10s \n",traj.atom_nr[h],       traj.atom_name[h].c_str(),       traj.res_nr[h],       traj.res_name[h].c_str());
            printf("angle %f dist %f \n",angle,dist);
            printf("sel a, resi %d & resn %s & name %s \n",traj.res_nr[acceptor],traj.res_name[acceptor].c_str(),traj.atom_name[acceptor].c_str());
            printf("sel d, resi %d & resn %s & name %s \n",traj.res_nr[donor],traj.res_name[donor].c_str(),traj.atom_name[donor].c_str());
            printf("sel h, resi %d & resn %s & name %s \n",traj.res_nr[h],traj.res_name[h].c_str(),traj.atom_name[h].c_str());
            printf("get_angle (a), (d), (h), state=4  \n");
            printf("dist (a), (d) \n");
            printf("show licorice, resi %d \n",traj.res_nr[acceptor]);
            printf("show licorice, resi %d \n",traj.res_nr[donor]);
            printf("sel triad, a + d + h \n");
            printf("orient triad \n");
            printf("\n");
        }
    }

    return result; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the number of lipid-prot h-bonds and adds it to the grid                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void lip_h_bonds(Trajectory &traj,system_variables &s,program_variables &p,Index &param,Index &lip_a,Index &lip_d,Index &prot_a,Index &prot_d,iv2d &bonds,Grid &hb)
{
    int    i        = 0;                      //standard variable used in loops
    int    j        = 0;                      //standard variable used in loops
    int    k        = 0;                      //standard variable used in loops
    int    l        = 0;                      //standard variable used in loops
    int    m        = 0;                      //standard variable used in loops
    int    n        = 0;                      //standard variable used in loops
    int    o        = 0;                      //standard variable used in loops
    int    q        = 0;                      //standard variable used in loops
    double distance = 0;                      //how far between atoms
    double hx       = 0;                      //head atom x-component 
    double hy       = 0;                      //head atom y-component

    //clear the current frame grids
    hb.clean_frame();

    for(i=0; i<traj.target_leaflet.size(); i++) //loop over target leaflet atoms
    {
        //get the first and last atom of the current lipid
        int min = traj.t_lip_start(i);
        int max = traj.t_lip_end(i);

        //jump to the next lipid
        i = traj.next_target_lipid(i);

        for(j=0; j<param.index_s.size(); j+=3) //loop over lipid types
        {
            if(strcmp(traj.res_name[traj.target_leaflet[i]-1].c_str(), param.index_s[j].c_str()) == 0) //lipid type is correct
            {
                int contacts = 0;

                for(k=min; k<=max; k++) //loop over current residue atoms
                {
                    //get lip_d-prot_a h-bonds
                    for(l=0; l<lip_d.index_s.size(); l++) //loop over donor lipid atoms
                    {
                        if(strcmp(traj.atom_name[k].c_str(), lip_d.index_s[l].c_str()) == 0) //atom is a donor lipid atom
                        {
                            int donor    = k;
                            int acceptor = 0;
                            int h        = 0;

                            for(m=0; m<bonds[donor].size(); m++) //loop over bonds 
                            {
                                if(traj.atom_name[bonds[donor][m]-1].at(0) == 'H') //atom is a hydrogen
                                {
                                    h = bonds[donor][m]-1;

                                    for(o=0; o<traj.prot.size(); o++) //loop over protein atoms
                                    {
                                        for(q=0; q<prot_a.index_s.size(); q++) //loop over acceptor protein atoms
                                        {
                                            if(strcmp(traj.atom_name[traj.prot[o]-1].c_str(), prot_a.index_s[q].c_str()) == 0) //atom is an acceptor protein atom
                                            {
                                                acceptor = traj.prot[o]-1;

                                                //check for h-bond
                                                contacts = contacts + check_h_bond(traj,p,acceptor,donor,h);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    //get lip_a-prot_d h-bonds
                    for(l=0; l<lip_a.index_s.size(); l++) //loop over acceptor lipid atoms
                    {
                        if(strcmp(traj.atom_name[k].c_str(), lip_a.index_s[l].c_str()) == 0) //atom is an acceptor lipid atom
                        {
                            int donor    = 0;
                            int acceptor = k;
                            int h        = 0;
 
                            for(m=0; m<traj.prot.size(); m++) //loop over protein atoms
                            {
                                for(n=0; n<prot_d.index_s.size(); n++) //loop over donor protein atoms
                                {
                                    if(strcmp(traj.atom_name[traj.prot[m]-1].c_str(), prot_d.index_s[n].c_str()) == 0) //atom is a donor protein atom
                                    {
                                        donor = traj.prot[m]-1;

                                        for(o=0; o<bonds[donor].size(); o++) //loop over bonds 
                                        {
                                            if(traj.atom_name[bonds[donor][o]-1].at(0) == 'H') //atom is a hydrogen
                                            {
                                                //printf("HBOND acceptor %5s %5s donor %5s %5s bond %5s \n",traj.atom_name[acceptor].c_str(),traj.res_name[acceptor].c_str(),traj.atom_name[donor].c_str(),traj.res_name[donor].c_str(),traj.atom_name[bond.index_i[o+1]-1].c_str());

                                                h = bonds[donor][o]-1;

                                                //check if h-bond
                                                contacts = contacts + check_h_bond(traj,p,acceptor,donor,h);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                //add the number of contacts to the grid
                for(k=min; k<=max; k++) //loop over current residue atoms
                {
                    if(strcmp(traj.atom_name[k].c_str(), param.index_s[j+1].c_str()) == 0) //mapping atom 1
                    {
                        hx = traj.r[k][0]; 
                        hy = traj.r[k][1];

                        hb.stamp(hx,hy,p.radius,contacts);
                    }
                    else if(strcmp(traj.atom_name[k].c_str(), param.index_s[j+2].c_str()) == 0) //mapping atom 2
                    {
                        hx = traj.r[k][0];
                        hy = traj.r[k][1];

                        hb.stamp(hx,hy,p.radius,contacts);
                    }
                }
            }
        }
    }

    //get the average for the current frame
    hb.norm_frame();

    //add the current frame grid to long term sum
    hb.add_frame();

    //now we print the single frame lpsb data
    if(p.b_stdev == 1)
    {
        hb.write_frame(traj.get_frame_global());
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Collect the contacts and compute the average                                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,Grid &hb)
{
    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nFinalizing analysis. This requires communicating data over the grid and could take some time depending on the resolution. \n");
    }

    //collect hb and rho from all ranks
    hb.collect_grid();

    //normalize hb and write the <hb> and rho to file
    if(s.world_rank == 0)
    {
        hb.normalize();

        hb.exclude_data(p.cutoff,1);

        hb.write_grid();
        hb.write_rho();
    }

    //compute standard deviation
    hb.get_stdev(p.b_stdev,p.b_clean,traj);

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
    s.program_name = "Lipid H-Bonds";

    //force program to run in serial?
    enum Switch serial         = off;

    //here we check if the program supports parallelization or not
    check_serial(s.world_rank,s.world_size,serial);

    //print program name and input arguments
    name_and_record(s.world_rank,argc,argv,s.program_name); 

    //analyze the command line arguments 
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,               "Input trajectory file (xtc, trr, pdb, gro)",                   s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,              "Refference file (pdb, gro)",                                   s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-o",      p.out_file_name,              "Output trajectory file (xtc, trr, pdb, gro)",                  s.world_rank, &p.b_print,   0);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                    "Read every 'stride' frame",                                    s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,               "Skip frames before this number",                               s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                 "Skip frames after this number",                                s.world_rank, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-lsq",    p.lsq_index_file_name,        "Index for lsq fitting (ndx)",                                  s.world_rank, &p.b_lsq,     0);
    add_argument_mpi_i(argc,argv,"-lsq_d",  &p.lsq_dim,                   "Dimension for lsq fitting (3:x,y,z 2:x,y)",                    s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-lsq_r",  &p.lsq_ref,                   "Reference structure for lsq fitting (0:ref 1:first_frame)",    s.world_rank, nullptr,      0);
    add_argument_mpi_s(argc,argv,"-crd",    p.param_file_name,            "Lipid types selection card + mapping atoms (crd)",             s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lip_a",  p.lip_a_file_name,            "Selection card with lipid h-bond acceptor atom types (crd)",   s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lip_d",  p.lip_d_file_name,            "Selection card with lipid h-bond donor atom types (crd)",      s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-prot_a", p.prot_a_file_name,           "Selection card with protein h-bond acceptor atom types (crd)", s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-prot_d", p.prot_d_file_name,           "Selection card with protein h-bond donor atom types (crd)",    s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-bond",   p.bond_file_name,             "Selection card identifying protein bonding pairs (crd)",       s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lphb",   p.lphb_file_name,             "Output file with spatially resolved H-bond count (dat)",       s.world_rank, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,           "PDB file with sorted leaflets (pdb)",                          s.world_rank, &p.b_lf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,  "File with additional leaflet finder parameters (prm)",         s.world_rank, &p.b_lf_param,0);
    add_argument_mpi_s(argc,argv,"-pf_pdb", p.pf_pdb_file_name,           "PDB file with selected protein (pdb)",                         s.world_rank, &p.b_pf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-pf_prm", p.protein_finder_param_name,  "File with additional protein finder parameters (prm)",         s.world_rank, &p.b_pf_param,0);
    add_argument_mpi_s(argc,argv,"-sf_pdb", p.sf_pdb_file_name,           "PDB file with selected sol (pdb)",                             s.world_rank, &p.b_sf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-sf_prm", p.solvent_finder_param_name,  "File with additional solvent finder parameters (prm)",         s.world_rank, &p.b_sf_param,0);
    add_argument_mpi_d(argc,argv,"-APS",    &p.APS,                       "Area per grid box (nm^2)",                                     s.world_rank, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-r",      &p.radius,                    "Radius of the target atom (nm)",                               s.world_rank, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-cutoff", &p.cutoff,                    "Cutoff for excluding data (chi)",                              s.world_rank, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-bx",     &p.box_x,                     "Grid x dimension (nm)",                                        s.world_rank, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-by",     &p.box_y,                     "Grid y dimension (nm)",                                        s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-leaf",   &p.leaflet,                   "Which leaflet? (0:both 1:upper 2:lower)",                      s.world_rank, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-odf",    &p.out_data_format,           "What is the output data format? (0:matrix 1:vector)",          s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-stdev",  &p.b_stdev,                   "Compute the STDEV and STEM? (0:no 1:yes)",                     s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-clean",  &p.b_clean,                   "Remove single frame files? (0:no 1:yes)",                      s.world_rank, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-test",   &p.b_test,                    "Print info for checking hydrogen bonds? (0:no 1:yes)",         s.world_rank, nullptr,      0);
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
    check_extension_mpi(s.world_rank,"-crd",p.param_file_name,".crd");
    check_extension_mpi(s.world_rank,"-lip_a",p.lip_a_file_name,".crd");
    check_extension_mpi(s.world_rank,"-lip_d",p.lip_d_file_name,".crd");
    check_extension_mpi(s.world_rank,"-prot_a",p.prot_a_file_name,".crd");
    check_extension_mpi(s.world_rank,"-prot_d",p.prot_d_file_name,".crd");
    check_extension_mpi(s.world_rank,"-bond",p.prot_d_file_name,".crd");

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
    if(p.b_sf_pdb == 1)
    {
        check_extension_mpi(s.world_rank,"-sf_pdb",p.sf_pdb_file_name,".pdb");
    }
    if(p.b_sf_param == 1)
    {
        check_extension_mpi(s.world_rank,"-sf_prm",p.solvent_finder_param_name,".prm");
    }

    //create index objects
    Index lip_a;
    Index lip_d;
    Index prot_a;
    Index prot_d;
    Index bond;
    Index param;

    //read the index files
    lip_a.get_index(p.lip_a_file_name);
    lip_d.get_index(p.lip_d_file_name);
    prot_a.get_index(p.prot_a_file_name);
    prot_d.get_index(p.prot_d_file_name);
    bond.get_index(p.bond_file_name);
    param.get_index(p.param_file_name);

    //run leaflet/proten/solvent finder
    traj.get_leaflets(p.leaflet,p.leaflet_finder_param_name,p.b_lf_param);
    traj.get_protein(p.protein_finder_param_name,p.b_pf_param);
    traj.get_solvent(p.solvent_finder_param_name,p.b_sf_param);

    //print a pdb with distinguished leaflets/protein/solvent
    traj.write_leaflets(p.lf_pdb_file_name,p.b_lf_pdb);
    traj.write_protein(p.pf_pdb_file_name,p.b_pf_pdb);
    traj.write_sol(p.sf_pdb_file_name,p.b_sf_pdb);

    //print info about the target leaflet
    traj.get_leaflet_stats();

    //print info about the lipid selection
    traj.get_lipid_selection_stats(param.get_column_s(3,0),"-crd");

    //print info about the protein
    traj.get_prot_stats();

    //print info about the water
    traj.get_sol_stats();

    //create a grid to hold lipid contacts
    Grid hb;

    //get the grid dimensions
    hb.get_dim(p.box_x,p.box_y,traj.ibox,p.APS);

    //set the output file name for grid
    hb.set_output(p.lphb_file_name,p.out_data_format);

    //print info about the grid
    hb.print_dim();

    //use bond list to create a list of bonds for each atom. Easier to look up this way.  
    iv2d bonds(traj.atoms(),iv1d(0,0));
    int i = 0;
    int j = 0;
    for(i=0; i<bond.index_i.size(); i+=2) //loop over bonds
    {
       int duplicate = 0;

       for(j=0; j<bonds[bond.index_i[i]-1].size(); j++) //loop over added bonds
       {
           if(bonds[bond.index_i[i]-1][j] == bond.index_i[i+1]) //atom already added
           {
               duplicate = 1;
           }
       }
       
       if(duplicate == 0) //atom not yet added
       {
           bonds[bond.index_i[i]-1].push_back(bond.index_i[i+1]);
       }
       else 
       {  
           printf("duplicate entry found in bonds list \n");
       }

       duplicate = 0;

       for(j=0; j<bonds[bond.index_i[i+1]-1].size(); j++) //loop over added bonds
       {
           if(bonds[bond.index_i[i+1]-1][j] == bond.index_i[i]) //atom already added
           {
               duplicate = 1;
           }
       }

       if(duplicate == 0) //atom not yet added
       {
           bonds[bond.index_i[i+1]-1].push_back(bond.index_i[i]);
       }
       else 
       {
           printf("duplicate entry found in bonds list \n");
       }
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

        lip_h_bonds(traj,s,p,param,lip_a,lip_d,prot_a,prot_d,bonds,hb);

        traj.set_beta_lf();

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect contacts from mpi processes and compute the average
    perf.log_time(finalize_analysis(traj,s,p,hb),"Fin Ana");

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

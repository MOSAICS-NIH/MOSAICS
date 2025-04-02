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
#include "MosAT/program_variables/pv_mean_lipid_coords_3d.h"//This has the variables specific to the analysis program
#include "headers/array.h"                                   //This has routines used for working with arrays
#include "headers/performance.h"                             //This has a class for logging performance data
#include "headers/index.h"                                   //This has a class for working with index files
#include "headers/traj.h"                                    //This has a class for working with the trajectory
#include "headers/leaflet_finder.h"                          //This has routines used to find leaflets in membrane simulations
#include "headers/protein_finder.h"                          //This has routines used to find protein atoms
#include "headers/sol_finder.h"                              //This has routines used to find the solvent
#include "headers/grid.h"                                    //This has routines used for working with a grid
#include "headers/grid_3d.h"                                 //This has routines used for working with a grid
#include "headers/protein.h"                                 //This has routines used for working with protein data
#include "headers/force_serial.h"                            //This has routines used for forcing the code to run on a single mpi process
#include "headers/ref.h"                                     //This has routines used for reading in ref files (single frame)

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function dumps the grid of lipid coords to a pdb file.                                               //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void dump_grid_pdb(Trajectory &traj,system_variables &s,program_variables &p,vector<Grid_3d_d> &coords_x,vector <Grid_3d_d> &coords_y,
                   vector <Grid_3d_d> &coords_z,Grid_3d_i &nan,Grid_3d_d &rho,Ref &top)
{
    int i         = 0;                                              //standard variable used in loops
    int j         = 0;                                              //standard variable used in loops
    int k         = 0;                                              //standard variable used in loops
    int l         = 0;                                              //standard variable used in loops
    int count     = 0;                                              //used to determine the current atom when filling data structures                                                                
    int res_count = 0;                                              //used to determine the current res when filling data structures                                                               
    int size_grid = top.get_num_atoms()*p.ef_size_grid;                 //How many lines in the pdb required to hold all atoms

    //allocate memory for coords
    rvec *lip_coords;
    lip_coords = (rvec*)calloc(size_grid , sizeof(rvec));

    //create arrays to hold pdb stuff
    vector <int>    set_atom_nr(size_grid,0);                  //atom number used in pdb file
    vector <int>    set_res_nr(size_grid,0);                   //res number used in pdb file
    vector <string> set_atom_name(size_grid);                  //atom name used in pdb file
    vector <string> set_res_name(size_grid);                   //res name used in pdb file
    vector <double> lip_beta(size_grid);                       //beta factor used in pdb file
    vector <double> lip_weight(size_grid);                     //weight used ub  pdb file
    vector <string> lip_element(size_grid);                    //element collumn used in pdb file
    vector <char>   lip_chain_id(size_grid);                   //chain id used in pdb file

    for(i=0; i<p.num_g_z; i++) //loop over z
    {
        count     = 0;
        res_count = 0;

        //create info needed for pdb and get atomic coords
        for(j=0; j<p.num_g_y; j++) //loop over y
        {
            if(j%p.grid_stride == 0) //only print every grid_stride point
            {
                for(k=0; k<p.num_g_x; k++) //loop over x
                {
                    if(k%p.grid_stride == 0) //only print every grid_stride point
                    {
                        for(l=0; l<top.get_num_atoms(); l++) //loop over lipid atoms
                        {
                            if(nan.grid[i][j][k] == 0) //data is significant
                            {
                                //use average coords  
                                lip_coords[count][0] = coords_x[l].grid[i][j][k];
                                lip_coords[count][1] = coords_y[l].grid[i][j][k];
                                lip_coords[count][2] = coords_z[l].grid[i][j][k];
                            }
                            else //data is too rare
                            {
                                //set coords to the origin   
                                lip_coords[count][0] = 0;
                                lip_coords[count][1] = 0;
                                lip_coords[count][2] = 0;
                            }

                            //set the beta factor    
                            lip_beta[count] = rho.grid[i][j][k];

                            //set remaining pdb data                   
                            set_atom_nr[count]   = count;
                            set_res_nr[count]    = res_count;
                            set_atom_name[count] = top.atom_name[l];
                            set_res_name[count]  = top.res_name[l];
                            lip_weight[count]    = 0;
                            lip_element[count]   = "ca";
                            lip_chain_id[count]  = 'A';
                            count++;
                        }
                        res_count++;
                    }
                }
            }
        }

        //open file for printing output 
        FILE *mlc_file;
        string tag = "_" + to_string(i);
        string pdb_file_name = add_tag(p.mlc_file_name,tag);
        mlc_file = fopen(pdb_file_name.c_str(), "w");
        if(mlc_file == NULL)
        {
            printf("failure opening %s (pdb single). Make sure the file exists. \n",p.mlc_file_name.c_str());
        }
        else 
        {
            //write the average coords to a pdb file
            write_frame_pdb(traj.ibox,size_grid,set_atom_nr,set_res_nr,set_res_name,set_atom_name,lip_coords,top.title,s.world_rank,&mlc_file,lip_beta,lip_weight,lip_element,lip_chain_id,1);

            //close the pdb file
            fclose(mlc_file);
        }
    }
    //free memory not needed any more
    free(lip_coords);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function adds the lipid coords to the grid and keeps track of how many lipids were at each grid      //
// point so the average coords may be computed later and so that density data may be displayed and so that   //
// insignificant data may be identified and excluded from the final pdb.                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void lip_coords(Trajectory &traj,system_variables &s,program_variables &p,vector <Grid_3d_d> &sf_coords_x,vector <Grid_3d_d> &sf_coords_y,
                vector <Grid_3d_d> &sf_coords_z,vector <Grid_3d_d> &coords_x,vector <Grid_3d_d> &coords_y,vector <Grid_3d_d> &coords_z,Grid_3d_d &rho,
                Grid_3d_d &rho_t,Grid_3d_d &cell_count,Grid_3d_d &cell_count_t,Ref &top)
{
    int    i   = 0;                           //standard variable used in loops
    int    j   = 0;                           //standard variable used in loops
    int    k   = 0;                           //standard variable used in loops
    int    l   = 0;                           //standard variable used in loops

    //clear the current frame grids
    for(i=0; i<top.get_num_atoms(); i++) //loop over lipid atoms
    {
        sf_coords_x[i].clean_grid();
        sf_coords_y[i].clean_grid();
        sf_coords_z[i].clean_grid();
    }
    cell_count.clean_grid();
    cell_count_t.clean_grid();

    for(i=0; i<traj.target_leaflet.size(); i++) //loop over target leaflet atoms
    {
        //get the first and last atom of the current lipid
        int min = traj.t_lip_start(i);
        int max = traj.t_lip_end(i);

        //jump to the next lipid
        i = traj.next_target_lipid(i);

        for(j=min; j<=max; j++) //loop over current residue atoms
        {
            if(strcmp(traj.atom_name[j].c_str(), p.map_1.c_str()) == 0 || strcmp(traj.atom_name[j].c_str(), p.map_2.c_str()) == 0) //map atom 1 or 2 
            {
                double hx = traj.r[j][0];
                double hy = traj.r[j][1];
                double hz = traj.r[j][2];

                cell_count_t.stamp(hx,hy,hz,p.radius,1.0);

                //add the current lipids to the grid    
                if(strcmp(traj.res_name[traj.target_leaflet[i]-1].c_str(), p.target_lip.c_str()) == 0) //lipid type is correct
                {
                    for(k=min; k<=max; k++) //loop over current lipid atoms
                    {
                        for(l=0; l<top.get_num_atoms(); l++) //loop over reference file atoms
                        {
                            if(strcmp(traj.atom_name[k].c_str(), top.atom_name[l].c_str()) == 0) //atom was include in ref
                            {
                                //get the coords to be added to the grid
                                double x = traj.r[k][0];
                                double y = traj.r[k][1];
                                double z = traj.r[k][2];

                                //add coords to long term sum
                                sf_coords_x[l].stamp(hx,hy,hz,p.radius,x);
                                sf_coords_y[l].stamp(hx,hy,hz,p.radius,y);
                                sf_coords_z[l].stamp(hx,hy,hz,p.radius,z);
                            }
                        }
                    }
                    cell_count.stamp(hx,hy,hz,p.radius,1.0);
                }
            }
        }
    }

    //get the single frame average and add to the long term sums
    for(i=0; i<top.get_num_atoms(); i++) //loop over lipid atoms
    {
        //normalize current frame
        sf_coords_x[i].normalize(cell_count);
        sf_coords_y[i].normalize(cell_count);
        sf_coords_z[i].normalize(cell_count);

        //add to long term sums
        for(j=0; j<coords_x[i].num_z(); j++) //loop over z direction
        {
            for(k=0; k<coords_x[i].num_y(); k++) //loop over y direction
            {
                for(l=0; l<coords_x[i].num_x(); l++) //loop over x direction
                {
                    coords_x[i].grid[j][k][l] = coords_x[i].grid[j][k][l] + sf_coords_x[i].grid[j][k][l];
                    coords_y[i].grid[j][k][l] = coords_y[i].grid[j][k][l] + sf_coords_y[i].grid[j][k][l];
                    coords_z[i].grid[j][k][l] = coords_z[i].grid[j][k][l] + sf_coords_z[i].grid[j][k][l];
                }
            }
        }
    }

    //update rho and rho_t
    for(j=0; j<coords_x[0].num_z(); j++) //loop over z direction
    {
        for(k=0; k<coords_x[0].num_y(); k++) //loop over y direction
        {
            for(l=0; l<coords_x[0].num_x(); l++) //loop over x direction
            {
                //update rho
                if(cell_count.grid[j][k][l] > 0)
                {
                    rho.grid[j][k][l] = rho.grid[j][k][l] + 1.0;
                }

                //update rho_t
                if(cell_count_t.grid[j][k][l] > 0)
                {
                    rho_t.grid[j][k][l] = rho_t.grid[j][k][l] + 1.0;
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function normalizes the weights and does communication between nodes. It finally prints the average  //
// coords with density data. The average coords with spread data is printed later. Rho data is also printed  //
// here.                                                                                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,vector <Grid_3d_d> &coords_x,vector <Grid_3d_d> &coords_y,
                       vector <Grid_3d_d> &coords_z,Grid_3d_d &rho,Grid_3d_d &rho_t,Grid_3d_i &nan,Ref &top)
{
    int i = 0;                               //Standard variable used in loops
    int j = 0;                               //Standard variable used in loops
    int k = 0;                               //Standard variable used in loops

    MPI_Barrier(MPI_COMM_WORLD);

    //record time when beginning analysis
    s.t = clock();

    if(s.world_rank == 0)
    {
        printf("\nComputing Average coords.\n");
        printf("This requires communicating data over the grid and could take some time depending on the resolution. \n");
    }

    //collect coords and rho from all ranks
    if(s.world_size > 1)
    {
        rho.collect_grid();
        rho_t.collect_grid();

        for(i=0; i<top.get_num_atoms(); i++) //loop over lipid atoms
        {
            coords_x[i].collect_grid();
            coords_y[i].collect_grid();
            coords_z[i].collect_grid();
        }
    }

    //normalize coords and write the <coords> and rho to file
    if(s.world_rank == 0)
    {
        //normalize coords
        for(i=0; i<top.get_num_atoms(); i++) //loop over lipid atoms
        {
            coords_x[i].normalize(rho);
            coords_y[i].normalize(rho);
            coords_z[i].normalize(rho);
        }

        //set flags for excluding data in pdb file
        if(p.b_rho_t == 0)
        {
            exclude_data(rho_t,nan,p.cutoff,1);
        }
        else if(p.b_rho_t == 1)
        {
            exclude_data(rho,nan,p.cutoff,1);
        }
        else 
        {
            printf("Invalid value received for -rho_t. Please provide a value of 0 or 1. \n");
        }

        //write lipid density to output file
        rho_t.write_grid(nan,p.ex_val);
        rho.write_grid(nan,p.ex_val);

        //write aerage coords with density data to pdb
        dump_grid_pdb(traj,s,p,coords_x,coords_y,coords_z,nan,rho,top);
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
    s.program_name = "Mean Lipid Coords 3d";

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
    add_argument_mpi_s(argc,argv,"-param",  p.param_file_name,            "Ref file with target lipid (pdb, gro)",                       s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-mlc",    p.mlc_file_name,              "Output data file with average coords (pdb)",                  s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,           "PDB file with sorted leaflets (pdb)",                         s.world_rank, s.cl_tags, &p.b_lf_pdb,  0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,  "File with additional leaflet finder parameters (prm)",        s.world_rank, s.cl_tags, &p.b_lf_param,0);
    add_argument_mpi_s(argc,argv,"-m1",     p.map_1,                      "Name of mapping atom 1 ",                                     s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_s(argc,argv,"-m2",     p.map_2,                      "Name of mapping atom 2 ",                                     s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-APS",    &p.APS,                       "Area per grid box (nm^2)",                                    s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-r",      &p.radius,                    "Radius of the target atom (nm)",                              s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_d(argc,argv,"-cutoff", &p.cutoff,                    "Cutoff for excluding data (chi)",                             s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-bx",     &p.box_x,                     "Grid x dimension (nm)",                                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-by",     &p.box_y,                     "Grid y dimension (nm)",                                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-bz",     &p.box_z,                     "Grid Z dimension (nm)",                                       s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-leaf",   &p.leaflet,                   "Which leaflet? (0:both 1:upper 2:lower)",                     s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-g_strd", &p.grid_stride,               "How many grid points to skip in output pdb?",                 s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_d(argc,argv,"-ex_val", &p.ex_val,                    "Set excluded lattice points to this value",                   s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-rho_t",  &p.b_rho_t,                   "Use rho_t for excluding data (0:no, 1:yes)",                  s.world_rank, s.cl_tags, nullptr,      0);
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
    check_extension_mpi(s.world_rank,"-mlc",p.mlc_file_name,".pdb");

    if(p.b_lf_pdb == 1)
    {
        check_extension_mpi(s.world_rank,"-lf_pdb",p.lf_pdb_file_name,".pdb");
    }
    if(p.b_lf_param == 1)
    {
        check_extension_mpi(s.world_rank,"-lf_prm",p.leaflet_finder_param_name,".prm");
    }

    //run leaflet finder
    traj.get_leaflets(p.leaflet,p.leaflet_finder_param_name,p.b_lf_param);

    //print a pdb with distinguished leaflets
    traj.write_leaflets(p.lf_pdb_file_name,p.b_lf_pdb);

    //print info about the target leaflet
    traj.get_leaflet_stats();

    //create a ref file with lipid topology
    Ref top;
    top.get_ref(p.param_file_name,"-param");

    //get the lipid type
    p.target_lip = top.res_name[0];

    //print info about the lipid selection
    vector <string> lip_t;
    lip_t.push_back(p.target_lip);
    traj.get_lipid_selection_stats(lip_t,"-lip_t");

    //set the lenght of a cell for the grid
    p.cell_size = sqrt(p.APS);

    //counts grid points in x and y dimensions
    get_grid_size(p.box_x,p.box_y,p.box_z,traj.ibox,&p.num_g_x,&p.num_g_y,&p.num_g_z,p.cell_size);

    //print information about the grid
    print_grid_stats(p.box_x,p.box_y,p.box_z,traj.ibox,p.num_g_x,p.num_g_y,p.num_g_z,p.cell_size,s.world_rank);

    //estimate memoroy usage
    int lattice_size = p.num_g_z*p.num_g_y*p.num_g_x;
    double mem = 2.0*3.0*(double)top.get_num_atoms()*(double)lattice_size*8.0;   //vectors holding coords
    mem = mem + 4.0*(double)lattice_size*8.0;                                    //cell counts and rhos
    mem = mem + (double)lattice_size*4.0;                                        //nan
    mem = mem/1000000.0;                                                         //convert to mega bytes

    if(s.world_rank == 0)
    {
        printf("memory required for grid: %f (Mb) \n\n",mem);
        //printf("lip_num_atoms %d lattice_size %d num_g_x %d num_g_y %d num_g_z %d \n",top.get_num_atoms(),lattice_size,p.num_g_x,p.num_g_y,p.num_g_z);
    }

    //allocate memory for the grid to hold the coords, rho etc.
    vector <Grid_3d_d> sf_coords_x(top.get_num_atoms()); //holds the x coord of each atom
    vector <Grid_3d_d> sf_coords_y(top.get_num_atoms()); //holds the y coord of each atom
    vector <Grid_3d_d> sf_coords_z(top.get_num_atoms()); //holds the z coord of each atom 
    vector <Grid_3d_d> coords_x(top.get_num_atoms());    //holds the x coord of each atom
    vector <Grid_3d_d> coords_y(top.get_num_atoms());    //holds the y coord of each atom
    vector <Grid_3d_d> coords_z(top.get_num_atoms());    //holds the z coord of each atom
    Grid_3d_d cell_count;                            //hold the single frame normalization factor
    Grid_3d_d cell_count_t;                          //hold the single frame normalization factor for all lipids
    Grid_3d_d rho;                                   //holds the lipid density for target lipids
    Grid_3d_d rho_t;                                 //holds the lipid density for all lipids
    Grid_3d_i nan;                                   //flags which lattice points should be hidden 

    int i = 0;
    int j = 0;
    //set the dimensions of the grid
    for(i=0; i<top.get_num_atoms(); i++)
    {
        sf_coords_x[i].set_dim(p.APS,p.num_g_x,p.num_g_y,p.num_g_z);
        sf_coords_y[i].set_dim(p.APS,p.num_g_x,p.num_g_y,p.num_g_z);
        sf_coords_z[i].set_dim(p.APS,p.num_g_x,p.num_g_y,p.num_g_z);
        coords_x[i].set_dim(p.APS,p.num_g_x,p.num_g_y,p.num_g_z);
        coords_y[i].set_dim(p.APS,p.num_g_x,p.num_g_y,p.num_g_z);
        coords_z[i].set_dim(p.APS,p.num_g_x,p.num_g_y,p.num_g_z);
    }
    cell_count.set_dim(p.APS,p.num_g_x,p.num_g_y,p.num_g_z);
    cell_count_t.set_dim(p.APS,p.num_g_x,p.num_g_y,p.num_g_z);
    rho.set_dim(p.APS,p.num_g_x,p.num_g_y,p.num_g_z);
    rho_t.set_dim(p.APS,p.num_g_x,p.num_g_y,p.num_g_z);
    nan.set_dim(p.APS,p.num_g_x,p.num_g_y,p.num_g_z);

    //set the output file name for grid
    p.rho_file_name = chop_and_add_tag(p.mlc_file_name,"_rho.dx");
    p.rho_t_file_name = chop_and_add_tag(p.rho_file_name,"_rho_t.dx");
    rho_t.set_output(p.rho_t_file_name);
    rho.set_output(p.rho_file_name);

    //get grid size after accounting for grid stride
    p.size_x       = ceil((double)p.num_g_x/(double)p.grid_stride);     //How many grid points in x-direction to include in the pdb (after g_stride)
    p.size_y       = ceil((double)p.num_g_y/(double)p.grid_stride);     //How many grid points in y-direction to include in the pdb (after g_stride)
    p.ef_size_grid = p.size_x*p.size_y;                                 //How many grid points after accounting for grid stride
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

        lip_coords(traj,s,p,sf_coords_x,sf_coords_y,sf_coords_z,coords_x,coords_y,coords_z,rho,rho_t,cell_count,cell_count_t,top);

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect coords from mpi pocesses and compute average coords
    perf.log_time(finalize_analysis(traj,s,p,coords_x,coords_y,coords_z,rho,rho_t,nan,top),"Fin Ana");

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

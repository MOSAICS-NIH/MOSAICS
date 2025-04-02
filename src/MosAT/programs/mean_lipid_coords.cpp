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
#include "MosAT/program_variables/pv_mean_lipid_coords.h"   //This has the variables specific to the analysis program
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
#include "headers/ref.h"                                     //This has routines used for reading in ref files (single frame)

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function dumps the grid of lipid coords to a pdb file.                                               //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void dump_grid_pdb(Trajectory &traj,system_variables &s,program_variables &p,vector<Grid> &coords_x,vector <Grid> &coords_y,
                   vector <Grid> &coords_z,vector <Grid_d> &beta,Grid_i &nan,Ref &top,string &out_name)
{
    int i         = 0;                                              //standard variable used in loops
    int j         = 0;                                              //standard variable used in loops
    int k         = 0;                                              //standard variable used in loops
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

    //create info needed for pdb and get atomic coords
    for(j=0; j<p.num_g_y; j++) //loop over y
    {
        if(j%p.grid_stride == 0) //only print every grid_stride point
        {
            for(i=0; i<p.num_g_x; i++) //loop over x
            {
                if(i%p.grid_stride == 0) //only print every grid_stride point
                {
                    for(k=0; k<top.get_num_atoms(); k++) //loop over lipid atoms
                    {
                        if(nan.grid[j][i] == 0) //data is significant
                        {
                            //use average coords  
                            lip_coords[count][0] = coords_x[k].grid[j][i];
                            lip_coords[count][1] = coords_y[k].grid[j][i];
                            lip_coords[count][2] = coords_z[k].grid[j][i];
                        }
                        else //data is too rare
                        {
                            //set coords to the origin   
                            lip_coords[count][0] = 0;
                            lip_coords[count][1] = 0;
                            lip_coords[count][2] = 0;
                        }

                        //set the beta factor    
                        lip_beta[count] = beta[k].grid[j][i];

                        //set remaining pdb data                   
                        set_atom_nr[count]   = count;
                        set_res_nr[count]    = res_count;
                        set_atom_name[count] = top.atom_name[k];
                        set_res_name[count]  = top.res_name[k];
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
    mlc_file = fopen(out_name.c_str(), "w");
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
void lip_coords(Trajectory &traj,system_variables &s,program_variables &p,vector <Grid> &coords_x,vector <Grid> &coords_y,
                vector <Grid> &coords_z,Grid &rho_t,Ref &top)
{
    int    i   = 0;                           //standard variable used in loops
    int    j   = 0;                           //standard variable used in loops
    int    k   = 0;                           //standard variable used in loops
    int    l   = 0;                           //standard variable used in loops

    //clear the current frame grids
    for(i=0; i<top.get_num_atoms(); i++) //loop over lipid atoms
    {
        coords_x[i].clean_frame();
        coords_y[i].clean_frame();
        coords_z[i].clean_frame();
    }
    rho_t.clean_frame();

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

                rho_t.stamp(hx,hy,p.radius,1.0);

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
                                coords_x[l].stamp(hx,hy,p.radius,x);
                                coords_y[l].stamp(hx,hy,p.radius,y);
                                coords_z[l].stamp(hx,hy,p.radius,z);
                            }
                        }
                    }
                }
            }
        }
    }

    //get the single frame average and add to the long term sums
    for(i=0; i<top.get_num_atoms(); i++) //loop over lipid atoms
    {
        //normalize current frame
        coords_x[i].norm_frame();
        coords_y[i].norm_frame();
        coords_z[i].norm_frame();

        //add to long term sums
        coords_x[i].add_frame();
        coords_y[i].add_frame();
        coords_z[i].add_frame();
    }
    rho_t.norm_frame();
    rho_t.add_frame();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function normalizes the weights and does communication between nodes. It finally prints the average  //
// coords with density data. The average coords with spread data is printed later. Rho data is also printed  //
// here.                                                                                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,vector <Grid> &coords_x,vector <Grid> &coords_y,
                       vector <Grid> &coords_z,Grid &rho_t,vector <Grid_d> &beta,Grid_i &nan,Ref &top)
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
        //normalize and exclude insignificant data for coords
        for(i=0; i<top.get_num_atoms(); i++) //loop over lipid atoms
        {
            //copy the taget lipid density for the beta factor
            beta[i].grid = coords_x[i].rho;

            //normalize coords
            coords_x[i].normalize();
            coords_y[i].normalize();
            coords_z[i].normalize();
       
            //copy rho_t for excluding insignificant data
            if(p.b_rho_t == 1)
            {
                coords_x[i].copy_rho(rho_t);  
                coords_y[i].copy_rho(rho_t);
                coords_z[i].copy_rho(rho_t);
            }

            //exclude insignificant data
            coords_x[i].exclude_data(p.cutoff,0);
            coords_y[i].exclude_data(p.cutoff,0);
            coords_z[i].exclude_data(p.cutoff,0);
        }

        //write lipid density to output file
        coords_x[0].set_output(p.rho_file_name,p.out_data_format);
        coords_x[0].write_rho();

        //set flags for excluding data in pdb file
        nan.grid = coords_x[0].nan;

        //write aerage coords with density data to pdb
        dump_grid_pdb(traj,s,p,coords_x,coords_y,coords_z,beta,nan,top,p.mlc_file_name);
    }

    //now we broadcast coords_x,y,z so all cores have the average
    //this is needed only if we are computing the average dist from the average coords
    if(p.b_mean_dist == 1)
    {
        if(s.world_size > 1)
        {
            for(i=0; i<top.get_num_atoms(); i++) //loop over lipid atoms
            {
                coords_x[i].bcast_grid();
                coords_y[i].bcast_grid();
                coords_z[i].bcast_grid();
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //compute and return time spent in function
    return (clock() - s.t)/CLOCKS_PER_SEC;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the RMSD relative to the average coords of a given lattice point                   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double get_rmsd(vector <Grid> &coords_x,vector <Grid> &coords_y,vector <Grid> &coords_z,Ref &top,int target_x,int target_y,
		int min, int max,Trajectory &traj)
{
    //1. The routine for shifting the coords to the origin uses the atom numbers in the index[] to reference items in the massi[] array. This can be problematic
    //   since mosat starts numbering at 1. If an index is used (which it is not here) then the index would be shifted down by 1. However, when nullptr is given for 
    //   ind_cm then the center of mass is computed for atoms 0 to ncm (this is what is done here). If ind_reset != nullptr then atoms in ind_reset are centered. Otherwise,
    //   atoms from 0 up to nreset are centered. In the example used here, both optioins are nullptr so the center of mass is computed for the complete lipid and the complete
    //   lipid is centered at the origin. 
    //2. The reference lipid may be only part of the lipid. That is the user may choose to ommit some atoms in the reference structure (.gro file). When this is the case, the 
    //   program will build an rvec containing only atoms included in the ref molecule for fitting/rmsd calculations. 
    //3. For do_fit_ndim() the masses are used to select atoms included in the fit (all atoms with a non zero mass). In the present case all atoms are given a mass of 1. 
    //4. The rmsdev() function computes the rmsd over the entire molecule (no index provided). There is a second function (using an index) if it is desireable to consider only
    //   part of the molecule (not done here). 
    //5. see build/installation/include/gromacs/math/do_fit.h for notes about these functions
    //   void reset_x_ndim(int ndim, int ncm, const int *ind_cm,int nreset, const int *ind_reset,rvec x[], const real mass[]);
    //   void do_fit_ndim(int ndim, int natoms, real *w_rls, const rvec *xp, rvec *x);
    //   real rmsdev(int natoms, real mass[], rvec x[], rvec xp[]);

    int i = 0;                           //standard variable used in loops
    int j = 0;                           //standard variable used in loops
    double my_rmsd = 0.0;                //rmsd of lipid compared to mean coords

    rvec *r_tmp;                         //holds current lipid coords
    rvec *r_ref_tmp;                     //holds ref lipid coords
    int lip_size = top.get_num_atoms();  //how many atoms in current lipid

    //allocate memory for holding coords
    r_tmp     = (rvec*)calloc(lip_size , sizeof(rvec));
    r_ref_tmp = (rvec*)calloc(lip_size , sizeof(rvec));

    //get coords for ref lipid
    for(i=0; i<lip_size; i++) //loop over lipid atoms
    {
        r_ref_tmp[i][0] = coords_x[i].grid[target_y][target_x];
	r_ref_tmp[i][1] = coords_y[i].grid[target_y][target_x];
        r_ref_tmp[i][2] = coords_z[i].grid[target_y][target_x];
    }

    //get coords for current lipid frame
    for(i=min; i<=max; i++) //loop over current lipid atoms
    {
        for(j=0; j<top.get_num_atoms(); j++) //loop over reference file atoms
        {   
            if(strcmp(traj.atom_name[i].c_str(), top.atom_name[j].c_str()) == 0) //atom was include in ref
            {   
                r_tmp[j][0] = traj.r[i][0];
                r_tmp[j][1] = traj.r[i][1];
                r_tmp[j][2] = traj.r[i][2];
            }
        }
    }

    //center lipid/ref at origin
    real   dummy_mass[lip_size];
    for(i=0; i<lip_size; i++)
    {
        dummy_mass[i] = 1.0;
    }
    reset_x_ndim(3, lip_size, nullptr, lip_size, nullptr,r_tmp, dummy_mass);
    reset_x_ndim(3, lip_size, nullptr, lip_size, nullptr,r_ref_tmp, dummy_mass);

    //perform least squares fitting
    do_fit_ndim(3, lip_size, dummy_mass, r_ref_tmp, r_tmp);

    //measure rmsd
    my_rmsd = rmsdev(lip_size, dummy_mass, r_tmp, r_ref_tmp);

    //free memory
    free(r_tmp);
    free(r_ref_tmp);

    return my_rmsd; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function deposits the rmsd to the grid, which requires special stamping as the value can vary        // 
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void deposit_rmsd_to_grid_v(program_variables &p,double hx,double hy,vector <Grid> &coords_x,vector <Grid> &coords_y,
                            vector <Grid> &coords_z,Grid &rmsd,int min,int max,Trajectory &traj,Ref &top)
{
    int i       = 0;         //standard variable used in loops
    int j       = 0;         //standard variable used in loops
    int lower_x = 0;         //starting cell in the x direction needed to enclose the atom        
    int lower_y = 0;         //starting cell in the y direction needed to enclose the atom
    int upper_x = 0;         //end cell in the x direction needed to enclose the atom
    int upper_y = 0;         //end cell in the y direction needed to enclose the atom

    //determine the upper and lower cells needed to enclose the atom
    get_bounds(hx,hy,p.radius,p.cell_size,p.num_g_x,p.num_g_y,&lower_x,&lower_y,&upper_x,&upper_y);

    for(i=lower_x; i<upper_x; i++) //loop over the grid x-axis
    {
        for(j=lower_y; j<upper_y; j++) //loop over the grid y-axis
        {
            //compute the distance the cell point is from the mapping atom
            double dist_x = i*p.cell_size - hx;
            double dist_y = j*p.cell_size - hy;

            double dist = sqrt(dist_x*dist_x + dist_y*dist_y);

            if(dist <= p.radius)
            {
                //compute rmsd for each grid point
                double this_rmsd = get_rmsd(coords_x,coords_y,coords_z,top,i,j,min,max,traj);

                //add the rmsd to current lattice point
                rmsd.frame_grid[j][i] = rmsd.frame_grid[j][i] + this_rmsd;
                rmsd.cell_count[j][i] = rmsd.cell_count[j][i] + 1;
            }
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function deposits the distance from the mean to the grid which requires special stamping as the      // 
// value can vary.                                                                                           // 
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void deposit_dist_to_grid_v(program_variables &p,double hx,double hy,vector <Grid> &coords_x,vector <Grid> &coords_y,
                            vector <Grid> &coords_z,int current_atom_index,vector <Grid> &sum_dist,double rx,double ry,double rz)
{
    int i       = 0;         //standard variable used in loops
    int j       = 0;         //standard variable used in loops
    int lower_x = 0;         //starting cell in the x direction needed to enclose the atom        
    int lower_y = 0;         //starting cell in the y direction needed to enclose the atom
    int upper_x = 0;         //end cell in the x direction needed to enclose the atom
    int upper_y = 0;         //end cell in the y direction needed to enclose the atom

    //determine the upper and lower cells needed to enclose the atom
    get_bounds(hx,hy,p.radius,p.cell_size,p.num_g_x,p.num_g_y,&lower_x,&lower_y,&upper_x,&upper_y);

    for(i=lower_x; i<upper_x; i++) //loop over the grid x-axis
    {
        for(j=lower_y; j<upper_y; j++) //loop over the grid y-axis
        {
            //compute the distance the cell point is from the mapping atom
            double dist_x = i*p.cell_size - hx;
            double dist_y = j*p.cell_size - hy;

            double dist = sqrt(dist_x*dist_x + dist_y*dist_y);

            if(dist <= p.radius)
            {
                //compute the distance the atom is from the mean coord
                double dx = rx - coords_x[current_atom_index].grid[j][i];
                double dy = ry - coords_y[current_atom_index].grid[j][i];
                double dz = rz - coords_z[current_atom_index].grid[j][i];

                double r_dist = sqrt(dx*dx + dy*dy + dz*dz);

                //add the distance to the grid point
                sum_dist[current_atom_index].frame_grid[j][i] = sum_dist[current_atom_index].frame_grid[j][i] + r_dist;
                sum_dist[current_atom_index].cell_count[j][i] = sum_dist[current_atom_index].cell_count[j][i] + 1;
	    }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the distance each lipid atom is from the average coords and adds it to the long    //
// term sums.                                                                                                //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void lip_coords_dist(Trajectory &traj,system_variables &s,program_variables &p,vector <Grid> &coords_x,vector <Grid> &coords_y,
                vector <Grid> &coords_z,Grid &rho_t,Ref &top,vector <Grid> &sum_dist,Grid &rmsd)

{
    int    i   = 0;                           //standard variable used in loops
    int    j   = 0;                           //standard variable used in loops
    int    k   = 0;                           //standard variable used in loops
    int    l   = 0;                           //standard variable used in loops

    //clear the current frame grids
    for(i=0; i<top.get_num_atoms(); i++) //loop over lipid atoms
    {
        sum_dist[i].clean_frame();
    }
    rmsd.clean_frame();

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

                //add the current lipids to the grid    
                if(strcmp(traj.res_name[traj.target_leaflet[i]-1].c_str(), p.target_lip.c_str()) == 0) //lipid type is correct
                {
                    for(k=min; k<=max; k++) //loop over current lipid atoms
                    {
                        for(l=0; l<top.get_num_atoms(); l++) //loop over reference file atoms
                        {
                            if(strcmp(traj.atom_name[k].c_str(), top.atom_name[l].c_str()) == 0) //atom was include in ref
                            {
                                //get the coords. We will compute distances for each grid point when stamping. 
                                double x = traj.r[k][0];
                                double y = traj.r[k][1];
                                double z = traj.r[k][2];

                                //add distances to grid
				deposit_dist_to_grid_v(p,hx,hy,coords_x,coords_y,coords_z,l,sum_dist,x,y,z);
                            }
                        }
                    }

                    //put this outside so we dont do for multiple atoms, only once per mapping atom 
                    deposit_rmsd_to_grid_v(p,hx,hy,coords_x,coords_y,coords_z,rmsd,min,max,traj,top);
                }
            }
        }
    }

    for(i=0; i<top.get_num_atoms(); i++) //loop over lipid atoms
    {
        sum_dist[i].norm_frame();
        sum_dist[i].add_frame();
    }
    rmsd.norm_frame();
    rmsd.add_frame();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects sum_dist from all cores and computes the final average dist.                       //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void finalize_avg_dist(Trajectory &traj,system_variables &s,program_variables &p,vector <Grid> &coords_x,vector <Grid> &coords_y,
                       vector <Grid> &coords_z,Grid &rho_t,vector <Grid> &sum_dist,Grid &rmsd,Grid_i &nan,Ref &top,vector <Grid_d> &beta)
{
    int i = 0;                               //Standard variable used in loops

    MPI_Barrier(MPI_COMM_WORLD);

    if(s.world_rank == 0)
    {
        printf("\nComputing average distance from the mean coords.\n");
        printf("This requires communicating data over the grid and could take some time depending on the resolution. \n");
    }

    //collect sum_dist and rho from all ranks
    if(s.world_size > 1)
    {
        for(i=0; i<top.get_num_atoms(); i++) //loop over lipid atoms 
        {
            sum_dist[i].collect_grid();
        }
	rmsd.collect_grid();
    }

    //normalize sum_dist and write the <dist> to file
    if(s.world_rank == 0)
    {
        for(i=0; i<top.get_num_atoms(); i++) //loop over lipid atoms
        {
            sum_dist[i].normalize();

            if(p.b_rho_t == 1)
            {
                sum_dist[i].copy_rho(rho_t);
            }
            sum_dist[i].exclude_data(p.cutoff,0);

            beta[i].grid = sum_dist[i].grid;

            //create name for pdb file
            string out_file_name = chop_and_add_tag(p.mlc_file_name,"_dist.pdb");

            //write the average coords with average dist as a beta value to the pdb file
	    dump_grid_pdb(traj,s,p,coords_x,coords_y,coords_z,beta,nan,top,out_file_name);  
        }

	//now do the rmsd
	rmsd.normalize();
	if(p.b_rho_t == 1)
        {
            rmsd.copy_rho(rho_t);
        }
        rmsd.exclude_data(p.cutoff,0);
	for(i=0; i<top.get_num_atoms(); i++) //loop over lipid atoms
        {
            beta[i].grid = rmsd.grid;
        }	
        //create name for pdb file
        string out_file_name = chop_and_add_tag(p.mlc_file_name,"_rmsd.pdb");

        //write the average coords with average dist as a beta value to the pdb file
        dump_grid_pdb(traj,s,p,coords_x,coords_y,coords_z,beta,nan,top,out_file_name);
    }
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
    s.program_name = "Mean Lipid Coords";

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
    add_argument_mpi_s(argc,argv,"-param",  p.param_file_name,            "Ref file with target lipid (gro, pdb)",                       s.world_rank, s.cl_tags, nullptr,      1);
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
    add_argument_mpi_i(argc,argv,"-leaf",   &p.leaflet,                   "Which leaflet? (0:both 1:upper 2:lower)",                     s.world_rank, s.cl_tags, nullptr,      1);
    add_argument_mpi_i(argc,argv,"-dist",   &p.b_mean_dist,               "Compute the <dist> from <coords>? (0:no 1:yes)",              s.world_rank, s.cl_tags, nullptr,      0);
    add_argument_mpi_i(argc,argv,"-g_strd", &p.grid_stride,               "How many grid points to skip in output pdb?",                 s.world_rank, s.cl_tags, nullptr,      0);
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
    traj.get_lipid_selection_stats(lip_t,"-param");

    //set the lenght of a cell for the grid
    p.cell_size = sqrt(p.APS);

    //counts grid points in x and y dimensions
    get_grid_size(p.box_x,p.box_y,traj.ibox,&p.num_g_x,&p.num_g_y,p.cell_size);

    //allocate memory for the grid to hold the coords, rho etc. 
    vector <Grid> coords_x(top.get_num_atoms());    //holds the x coord of each atom
    vector <Grid> coords_y(top.get_num_atoms());    //holds the y coord of each atom
    vector <Grid> coords_z(top.get_num_atoms());    //holds the z coord of each atom
    vector <Grid> sum_dist(top.get_num_atoms());    //holds the average distance from meaan coord of each atom
    Grid rmsd;                                      //holds the average rmsd relative to the meaan coords at each grid point

    Grid rho_t;                                     //holds the lipid density for all lipids
    vector <Grid_d> beta(top.get_num_atoms());      //holds the beta value for each atom of the lipids
    Grid_i nan;                                     //flags which lattice points should be hidden 

    int i = 0;
    int j = 0;
    //set the dimensions of the grid
    for(i=0; i<top.get_num_atoms(); i++)
    {
        coords_x[i].get_dim(p.box_x,p.box_y,traj.ibox,p.APS);
        coords_y[i].get_dim(p.box_x,p.box_y,traj.ibox,p.APS);
        coords_z[i].get_dim(p.box_x,p.box_y,traj.ibox,p.APS);
	if(p.b_mean_dist == 1)
	{	
            sum_dist[i].get_dim(p.box_x,p.box_y,traj.ibox,p.APS);
        }
        beta[i].set_dim(p.APS,p.num_g_x,p.num_g_y);
    }
    if(p.b_mean_dist == 1)
    {
        rmsd.get_dim(p.box_x,p.box_y,traj.ibox,p.APS);
    } 
    rho_t.get_dim(p.box_x,p.box_y,traj.ibox,p.APS);
    nan.set_dim(p.APS,p.num_g_x,p.num_g_y);

    //print info about the grid
    coords_x[0].print_dim();

    //set the output file name for grid
    p.rho_file_name = chop_and_add_tag(p.mlc_file_name,".dat");

    //get grid size after accounting for grid stride
    p.size_x       = ceil((double)p.num_g_x/(double)p.grid_stride);     //How many grid points in x-direction to include in the pdb (after g_stride)
    p.size_y       = ceil((double)p.num_g_y/(double)p.grid_stride);     //How many grid points in y-direction to include in the pdb (after g_stride)
    p.ef_size_grid = p.size_x*p.size_y;                                 //How many grid points after accounting for grid stride

    //print memory estimates
    int mem_per_grid = (2*p.num_g_x*p.num_g_y*4) + (4*p.num_g_x*p.num_g_y*8);   //2 int grids (4 bytes) and 4 double grids (8 bytes)
    int mem_grid     = (4*mem_per_grid*top.get_num_atoms()) + (4*mem_per_grid); //4 grids for each atom and 1 for rho 
    int mem_grid_d   = 2*p.num_g_x*p.num_g_y*top.get_num_atoms()*8;             //beta and rmsd
    int mem_grid_i   = p.num_g_x*p.num_g_y*top.get_num_atoms()*4;               //nan
    double mem_tot   = (double)(mem_grid + mem_grid_d + mem_grid_i)/1000000.0;  //mega bytes
    if(s.world_rank == 0)
    {
        printf("Estimated memory for grid data: %f MB \n\n",mem_tot);
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

        lip_coords(traj,s,p,coords_x,coords_y,coords_z,rho_t,top);

        traj.write_traj_frame();

        time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
    }

    //log time spent in main loop
    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");

    //collect coords from mpi pocesses and compute average coords
    perf.log_time(finalize_analysis(traj,s,p,coords_x,coords_y,coords_z,rho_t,beta,nan,top),"Fin Ana");

    //loop over the trajectory a second time to get the average dist from mean
    if(p.b_mean_dist == 1)
    {
        //reset counter for printing progress + time estimates
        s.counter = 0;
        s.t = clock();

        if(s.world_rank == 0)
        {
            printf("\nComputing fluctuations.\n");
            printf("-----------------------------------------------------------------------------------------------------------------------------------\n");
            fflush(stdin);
        }

	for(traj.current_frame=0; traj.current_frame<traj.get_num_frames(); traj.current_frame++)
        {
            traj.read_traj_frame();

            traj.do_fit();

            lip_coords_dist(traj,s,p,coords_x,coords_y,coords_z,rho_t,top,sum_dist,rmsd);

            //report progress
            time_stats(s.t,&s.counter,traj.current_frame,traj.get_num_frames(),s.world_rank);
	}
        perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Compute Mean Dist.");

        s.t = clock();
        finalize_avg_dist(traj,s,p,coords_x,coords_y,coords_z,rho_t,sum_dist,rmsd,nan,top,beta);
        perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Fin. Analysis dist.");
    }

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

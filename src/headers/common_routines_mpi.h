
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <mpi.h>

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects the projection (iterations) data                                                   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void collect_local_values(int world_size,int world_rank,double world_values[],double my_values[],iv1d &world_iterations,
                            string values_string)
{
    int i = 0;                                  //general variable used in loops
    int j = 0;                                  //general variable used in loops
    int k = 0;                                  //general variable used in loops
    int current_i = 0;                          //keep track of global iteration value
    MPI_Status status;

    //rank 0 copies its local values to the global_values array
    if(world_rank == 0)
    {
        for(j=0; j<world_iterations[0]; j++) //loop over rank 0 iterations
        {
            world_values[j] = my_values[j];
            current_i++;
        }
    }

    if(world_size > 1) //only gather if world size is greater than 1
    {
        MPI_Barrier(MPI_COMM_WORLD);

        if(world_rank == 0)
        {
            printf("Collecting %s \n",values_string.c_str());
        }

        //now rank 0 collects local values from other cores
        for(i=1; i<world_size; i++) //loop over world
        {
            if(world_rank == 0) //receive data
            {
                for(j=0; j<world_iterations[i]; j++)
                {
                    //allocate memory to store incoming data
                    double value_rec;

                    //receive data
                    MPI_Recv(&value_rec, 1, MPI_DOUBLE, i, 13, MPI_COMM_WORLD, &status);

                    world_values[current_i] = value_rec;
                    current_i++;
                }
            }
            else if(world_rank == i) //send data
            {
                for(j=0; j<world_iterations[world_rank]; j++)
                {
                    MPI_Send(&my_values[j], 1, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD);
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function determines which grid points each core is responsible for.                                  //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int count_workload(int world_size,int world_rank,int num_g_x)
{
    int i=0;
    int my_points = 0;

    //count how many grid points each core is responsible for
    for(i=0; i<num_g_x; i++)
    {
        if(i%world_size == world_rank)
        {
            my_points++;
        }
    }
    return my_points;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function determines the upper and lower bounds for num_g_x for each core.                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_workload(int *my_xi,int *my_xf,int world_rank,iv1d &world_num_g_x,int num_g_x,int world_xi[],int world_xf[])
{
    int i = 0;
    int rank = 0;
    int count = 0;
    int first = 1;

    for(i=0; i<num_g_x; i++)
    {
        if(rank == world_rank)
        {
            if(first == 1)
            {
                *my_xi = i;
                first = 0;
            }
            *my_xf = i;
        }

        count++;
        if(count == world_num_g_x[rank])
        {
            rank++;
            count = 0;
        }
    }

    //collect xi and xf
    MPI_Allgather(my_xi, 1,MPI_INT,world_xi, 1, MPI_INT, MPI_COMM_WORLD );
    MPI_Allgather(my_xf, 1,MPI_INT,world_xf, 1, MPI_INT, MPI_COMM_WORLD );
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints statistics for each core describing how the grid was split between them.             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void print_workload_stats(int world_rank,int world_xi[],int world_xf[],iv1d &world_num_g_x,int world_size,string variable,string init,string fin)
{
    int i =0;

    if(world_rank == 0)
    {
        printf("Distributing the workload accross %d mpi processes. \n",world_size);
        printf(" %10s %10s %10s %10s \n","rank",variable.c_str(),init.c_str(),fin.c_str());
        printf("-%10s-%10s-%10s-%10s \n","----------","----------","----------","----------");
        for(i=0; i<world_size; i++)
        {
            printf(" %10d %10d %10d %10d \n",i,world_num_g_x[i],world_xi[i],world_xf[i]);
        }
        printf("\n");
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function gathers data (double) for the grid points. Used when parallelization scheme is by grid point//
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void gather_grid_d_gp(int world_size,int world_rank,int my_num_g_x,int num_g_x,int num_g_y,iv1d &world_num_g_x,dv2d &local_data,dv2d &global_data)
{
    if(world_size > 1) //only gather if world size is greater than 1
    {
        int i = 0;                                  //general variable used in loops
        int j = 0;                                  //general variable used in loops
        int k = 0;                                  //general variable used in loops
        int current_x = 0;                          //keep track of global x value
        MPI_Status status;                          //Used for mpi_recv

        MPI_Barrier(MPI_COMM_WORLD);

        //rank 0 copies its local data to the global_data array
        if(world_rank == 0)
        {
            for(j=0; j<my_num_g_x; j++) //loop over x
            {
                for(k=0; k<num_g_y; k++) //loop over y
                {
                    global_data[k][current_x] = local_data[k][j];
                }
                current_x++;
            }
        }

        //now rank 0 collects data from other cores
        for(i=1; i<world_size; i++) //loop over world
        {
            if(world_rank == 0) //receive data
            {
                //allocate memory to store incoming data
                double data_rec[num_g_y][world_num_g_x[i]];

                //receive data
                MPI_Recv(data_rec, world_num_g_x[i]*num_g_y, MPI_DOUBLE, i, 13, MPI_COMM_WORLD, &status);

                //add received data to the global data
                for(j=0; j<world_num_g_x[i]; j++) //loop over x
                {
                    for(k=0; k<num_g_y; k++) //loop over y
                    {
                        global_data[k][current_x] = data_rec[k][j];
                    }
                    current_x++;
                }
            }
            else if(world_rank == i) //send data
            {
                //copy data to array for sending
                double data_snd[num_g_y][my_num_g_x];
                for(j=0; j<my_num_g_x; j++) //loop over x
                {
                    for(k=0; k<num_g_y; k++) //loop over y
                    {
                        data_snd[k][j] = local_data[k][j];
                    }
                }

                MPI_Send(data_snd, my_num_g_x*num_g_y, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function determines the upper and lower bounds for num_g for each core.                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_grid_points_alt(int *my_gi,int *my_gf,int world_rank,int world_size,vector <int> &world_num_g,int num_g_x,int num_g_y,vector <int> &world_gi,vector <int> &world_gf)
{
    int i     = 0;
    int j     = 0;
    int rank  = 0;
    int count = 0;
    int first = 1;
    int world_gi_ary[world_size];
    int world_gf_ary[world_size];
    int grid_count = 0;

    for(i=0; i<num_g_x; i++)
    {
        for(j=0; j<num_g_y; j++)
        {
            if(rank == world_rank)
            {
                if(first == 1)
                {
                    *my_gi = grid_count;
                    first = 0;
                }
                *my_gf = grid_count;
            }

            count++;
            if(count == world_num_g[rank])
            {
                rank++;
                count = 0;
            }

            grid_count++;
        }
    }

    //collect xi and xf
    MPI_Allgather(my_gi, 1,MPI_INT,world_gi_ary, 1, MPI_INT, MPI_COMM_WORLD );
    MPI_Allgather(my_gf, 1,MPI_INT,world_gf_ary, 1, MPI_INT, MPI_COMM_WORLD );

    //copy world_xi/xf arrays to the vectors
    for(i=0; i<world_size; i++)
    {
        world_gi[i] = world_gi_ary[i];
        world_gf[i] = world_gf_ary[i];
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints statistics for each core describing how the grid was split between them.             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void print_workload_stats_alt(int world_rank,vector <int> &world_gi,vector <int> &world_gf,vector <int> &world_num_g,int world_size)
{
    int i =0;

    if(world_rank == 0)
    {
        printf("Distributing the workload (the grid) accross %d mpi processes. \n",world_size);
        printf(" %4s %8s %8s %8s \n","rank","num_g","gi","gf");
        printf("-%4s-%8s-%8s-%8s \n","--------","--------","--------","--------");
        for(i=0; i<world_size; i++)
        {
            printf(" %4d %8d %8d %8d \n",i,world_num_g[i],world_gi[i],world_gf[i]);
        }
        printf("\n");
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function gathers data (double) for the grid points. Used when parallelization scheme is by grid point//
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void gather_grid_d_gp_alt(int world_size,int world_rank,iv1d &world_gi,iv1d &world_gf,int num_g_x,int num_g_y,iv1d &world_num_g,dv1d &local_data,dv2d &global_data)
{
    int i            = 0;                       //general variable used in loops
    int j            = 0;                       //general variable used in loops
    int k            = 0;                       //general variable used in loops
    int grid_counter = 0;                       //count grid points as they are encountered
    MPI_Status status;                          //Used for mpi_recv

    MPI_Barrier(MPI_COMM_WORLD);

    //rank 0 copies its local data to the global_data array
    if(world_rank == 0)
    {
        for(j=0; j<num_g_x; j++) //loop over x
        {
            for(k=0; k<num_g_y; k++) //loop over y
            {
                if(grid_counter >= world_gi[0] && grid_counter <= world_gf[0])	
                {
                    int ef_count = grid_counter - world_gi[0];

                    global_data[k][j] = local_data[ef_count];
                }
                grid_counter++;
            }
        }
    }

    //now rank 0 collects data from other cores
    for(i=1; i<world_size; i++) //loop over world
    {
        if(world_rank == 0) //receive data
        {
            int size = world_num_g[i];

            //allocate memory to store incoming data
            double data_rec[size];

            //receive data
            MPI_Recv(data_rec, size, MPI_DOUBLE, i, 13, MPI_COMM_WORLD, &status);

	    grid_counter = 0;
 
            //add received data to the global data
            for(j=0; j<num_g_x; j++) //loop over x
            {
                for(k=0; k<num_g_y; k++) //loop over y
                {
                    if(grid_counter >= world_gi[i] && grid_counter <= world_gf[i])
                    {
                        int ef_count = grid_counter - world_gi[i];

                        global_data[k][j] = data_rec[ef_count];
                    }
                    grid_counter++;
                }
            }
        }
        else if(world_rank == i) //send data
        {
            //copy data to array for sending
            int size = world_num_g[i];

            double data_snd[size];
            for(j=0; j<size; j++) //loop over x
            {
                data_snd[j] = local_data[j];
            }

            MPI_Send(data_snd, size, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function determines how much time has passed and gives an estimate of the time remaining.            //
// Takes in the current step starting at 1 (not zero).                                                       //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ot_time_stats(clock_t t,int *counter,int current_step,int my_steps,int world_rank,string my_string)
{
    if(world_rank == 0)
    {
        double seconds = 0;

        seconds = (clock() - t)/CLOCKS_PER_SEC;

        if((int)seconds/10 > *counter)
        {
            double percent_done         = ((double)(current_step)/(double)my_steps)*100.0;
            double estimated_total_time = 100.0*seconds/percent_done;
            double time_remaining       = estimated_total_time - seconds;

            if(current_step == my_steps)
            {
                estimated_total_time = 0.0;
            }

            int phr  = 0;
            int pmin = 0;
            int psec = 0;
            int lhr  = 0;
            int lmin = 0;
            int lsec = 0;

            phr  = (seconds)/(60*60);
            pmin = (seconds - (phr*60*60))/60;
            psec = seconds - (phr*60*60) - (pmin*60);

            lhr  = ((int)time_remaining)/(60*60);
            lmin = (time_remaining - (lhr*60*60))/60;
            lsec = time_remaining - (lhr*60*60) - (lmin*60);

            if(percent_done != 0)
            {
                printf("Finished %s %7d with %5.1f percent done overall in %2d hr %2d min %2d sec. Estimated time to completion is %2d hr %2d min %2d sec. \n",my_string.c_str(),current_step+1,percent_done,phr,pmin,psec,lhr,lmin,lsec);
            }
            *counter = *counter + 1;
        }
    }
}


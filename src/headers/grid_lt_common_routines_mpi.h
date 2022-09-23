#include <mpi.h>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes how many iterations the core is responsible for                                    //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int count_my_iterations(int iterations,int world_rank,int world_size)
{
    int i = 0;
    int j = 0;
    int my_count = 0;

    for(i=0; i<iterations; i++) //loop over iterations
    {
        if(i%world_size == world_rank)
        {
            my_count++;
        }
    }
    return my_count;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Main function computes the initial and final iteration number for each core                               //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void distribute_work_load(int iterations,int world_rank,int world_size,int *my_i,int *my_f,int world_i[],int world_f[],int world_iterations[])
{
    int i = 0;
    int rank = 0;
    int first = 1;
    int count = 0;

    for(i=0; i<iterations; i++) //loop over iterations
    {
        if(rank == world_rank)
        {
            if(first == 1)
            {
                *my_i = i;
                first = 0;
            }
            *my_f = i;
        }

        count++;
        if(count == world_iterations[rank])
        {
            rank++;
            count = 0;
        }
    }

    //collect xi and xf
    MPI_Allgather(my_i, 1,MPI_INT,world_i, 1, MPI_INT, MPI_COMM_WORLD );
    MPI_Allgather(my_f, 1,MPI_INT,world_f, 1, MPI_INT, MPI_COMM_WORLD );
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints information about how the workload was distributed                                   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void print_workload_stats(int world_i[],int world_f[],int world_iterations[],int world_rank,int world_size)
{
    int i = 0;

    if(world_rank == 0)
    {
        printf("\n");
        printf(" %10s %10s %10s %10s \n","Rank","Iterations","Initial","Final");
        printf("-%10s-%10s-%10s-%10s \n","----------","----------","----------","----------");
        for(i=0; i<world_size; i++) //loop over cores
        {
            printf(" %10d %10d %10d %10d \n",i,world_iterations[i],world_i[i],world_f[i]);
        }
        printf("\n");
   }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects the dwell times from each core                                                     //
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


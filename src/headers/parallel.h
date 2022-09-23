
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function determines which grid points each core is responsible for.                                  //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int count_grid_points(int world_size,int world_rank,int num_g_x)
{
    int i         = 0;
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
void get_grid_points(int *my_xi,int *my_xf,int world_rank,int world_size,vector <int> &world_num_g_x,int num_g_x,vector <int> &world_xi,vector <int> &world_xf)
{
    int i     = 0;
    int rank  = 0;
    int count = 0;
    int first = 1;
    int world_xi_ary[world_size];
    int world_xf_ary[world_size];

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
    MPI_Allgather(my_xi, 1,MPI_INT,world_xi_ary, 1, MPI_INT, MPI_COMM_WORLD );
    MPI_Allgather(my_xf, 1,MPI_INT,world_xf_ary, 1, MPI_INT, MPI_COMM_WORLD );

    //copy world_xi/xf arrays to the vectors
    for(i=0; i<world_size; i++)
    {
        world_xi[i] = world_xi_ary[i];
        world_xf[i] = world_xf_ary[i];
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints statistics for each core describing how the grid was split between them.             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void print_workload_stats(int world_rank,vector <int> &world_xi,vector <int> &world_xf,vector <int> &world_num_g_x,int world_size)
{
    int i =0;

    if(world_rank == 0)
    {
        printf("Distributing the workload (the grid) accross %d mpi processes. \n",world_size);
        printf(" %4s %8s %8s %8s \n","rank","num_g_x","xi","xf");
        printf("-%4s-%8s-%8s-%8s \n","--------","--------","--------","--------");
        for(i=0; i<world_size; i++)
        {
            printf(" %4d %8d %8d %8d \n",i,world_num_g_x[i],world_xi[i],world_xf[i]);
        }
        printf("\n");
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function breaks up the workload by grid points (num_g_x)                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::parallelize_by_grid(int num_g_x)
{
    int i = 0;

    //compute how many grid points each frame is responsible for
    my_num_g_x = count_grid_points(world_size,world_rank,num_g_x);

    //create array to hold each mpi processes my_num_g_x; Used for communication
    int world_num_g_x_ary[world_size];
    MPI_Allgather(&my_num_g_x, 1,MPI_INT,world_num_g_x_ary, 1, MPI_INT, MPI_COMM_WORLD );

    //allocate memory for world_num_g_x and copy data from the array
    world_num_g_x.resize(world_size,0);
    for(i=0; i<world_size; i++)
    {
        world_num_g_x[i] = world_num_g_x_ary[i];
    }

    //allocate memory for world_xi and world_xf
    world_xi.resize(world_size,0);
    world_xf.resize(world_size,0);

    //compute each mpi processes start and end grid point
    get_grid_points(&my_xi,&my_xf,world_rank,world_size,world_num_g_x,num_g_x,world_xi,world_xf);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints the workload distribution for grid points                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::workload_grid()
{
    print_workload_stats(world_rank,world_xi,world_xf,world_num_g_x,world_size);
}













///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function determines how many lipids each core is responsible for                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int compute_workload_lip(int num_lipids_1,int world_rank,int world_size)
{
    int i         = 0;
    int my_lipids = 0;

    for(i=0; i<num_lipids_1; i++)
    {
        if(i%world_size == world_rank)
        {
            my_lipids = my_lipids + 1;
        }
    }

    return my_lipids; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function determines which lipids each core is responsible for                                        //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void distribute_work_load_lip(int num_lipids_1,int world_rank,int world_size,int my_lipids,vector <int> &world_lipids,int *lipid_start,
                              int *lipid_end,vector <int> &world_lipid_start,vector <int> &world_lipid_end)
{
    int i             = 0;
    int running_count = 0;
    int world_lipid_start_ary[world_size];
    int world_lipid_end_ary[world_size];

    for(i=0; i<world_size; i++)
    {
        if(world_rank == i)
        {
            *lipid_start = running_count;
            *lipid_end   = running_count + world_lipids[i] - 1;
        }
        running_count = running_count + world_lipids[i];
    }

    MPI_Allgather(lipid_start, 1, MPI_INT, world_lipid_start_ary, 1, MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(lipid_end,   1, MPI_INT, world_lipid_end_ary,   1, MPI_INT,MPI_COMM_WORLD);

    for(i=0; i<world_size; i++)
    {
        world_lipid_start[i] = world_lipid_start_ary[i];
        world_lipid_end[i]   = world_lipid_end_ary[i];
    }

    //printf("world_rank %3d num_lipids_1 %5d lipid_start %5d lipid_end %5d \n",world_rank,num_lipids_1,*lipid_start,*lipid_end);    
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints information about the work load distribution                                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void print_work_load_stats_lip(int world_rank,int world_size,vector <int> &world_lipids,vector <int> &world_lipid_start,vector <int> &world_lipid_end)
{
    int i = 0;

    if(world_rank == 0)
    {
        printf(" %10s %10s %20s %20s \n","Rank","Lipids","Lipid_Start","Lipid_End");
        printf(" %10s-%10s-%20s-%20s \n","----------","----------","--------------------","--------------------");
        for(i=0; i<world_size; i++)
        {
            printf(" %10d %10d %20d %20d \n",i,world_lipids[i],world_lipid_start[i],world_lipid_end[i]);
        }
        printf("\n");
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function breaks up the workload by lipid                                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::parallelize_by_lipid(int num_lipids_1)
{
    int i = 0;

    //count the number of lipids the mpi process is responsible for
    my_lipids = compute_workload_lip(num_lipids_1,world_rank,world_size);

    //create an array for collecting my_lipids. then copy to the vector
    int world_lipids_ary[world_size];
    MPI_Allgather(&my_lipids, 1, MPI_INT, world_lipids_ary, 1, MPI_INT,MPI_COMM_WORLD);

    //allocate memory for world_lipids and copy world_lipids_ary 
    world_lipids.resize(world_size,0);
    for(i=0; i<world_size; i++)
    {
        world_lipids[i] = world_lipids_ary[i];
    }

    //allocate memory for world_lipid_start and world_lipid_end
    world_lipid_start.resize(world_size,0);
    world_lipid_end.resize(world_size,0);

    //compute each mpi processes start and end lipid
    distribute_work_load_lip(num_lipids_1,world_rank,world_size,my_lipids,world_lipids,&lipid_start,&lipid_end,world_lipid_start,world_lipid_end);

    //print_work_load_stats(s.world_rank,s.world_size,world_lipids,world_lipid_start,world_lipid_end);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints the workload distribution for lipid                                                  //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::workload_lipid()
{
    print_work_load_stats_lip(world_rank,world_size,world_lipids,world_lipid_start,world_lipid_end);
}














///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function determines how many waters each core is responsible for                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int compute_workload_wat(int num_waters,int world_rank,int world_size)
{
    int i         = 0;
    int my_waters = 0;

    for(i=0; i<num_waters; i++)
    {
        if(i%world_size == world_rank)
        {
            my_waters = my_waters + 1;
        }
    }

    return my_waters;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function determines which waters each core is responsible for                                        //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void distribute_work_load_wat(int num_waters,int world_rank,int world_size,int my_waters,vector <int> &world_waters,int *water_start,
                              int *water_end,vector <int> &world_water_start,vector <int> &world_water_end)
{
    int i             = 0;
    int running_count = 0;
    int world_water_start_ary[world_size];
    int world_water_end_ary[world_size];

    for(i=0; i<world_size; i++)
    {
        if(world_rank == i)
        {
            *water_start = running_count;
            *water_end   = running_count + world_waters[i] - 1;
        }
        running_count = running_count + world_waters[i];
    }

    MPI_Allgather(water_start, 1, MPI_INT, world_water_start_ary, 1, MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(water_end,   1, MPI_INT, world_water_end_ary,   1, MPI_INT,MPI_COMM_WORLD);

    for(i=0; i<world_size; i++)
    {
        world_water_start[i] = world_water_start_ary[i];
        world_water_end[i]   = world_water_end_ary[i];
    }

    //printf("world_rank %3d num_waters %5d water_start %5d water_end %5d \n",world_rank,num_waters,*water_start,*water_end);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints information about the work load distribution                                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void print_work_load_stats_wat(int world_rank,int world_size,vector <int> &world_waters,vector <int> &world_water_start,vector <int> &world_water_end)
{
    int i = 0;

    if(world_rank == 0)
    {
        printf(" %10s %10s %20s %20s \n","Rank","Waters","Water_Start","Water_End");
        printf(" %10s-%10s-%20s-%20s \n","----------","----------","--------------------","--------------------");
        for(i=0; i<world_size; i++)
        {
            printf(" %10d %10d %20d %20d \n",i,world_waters[i],world_water_start[i],world_water_end[i]);
        }
        printf("\n");
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function breaks up the workload by waters                                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::parallelize_by_water(int num_waters_1)
{
    int i = 0;

    //count the number of waterss the mpi process is responsible for
    my_waters = compute_workload_wat(num_waters_1,world_rank,world_size);

    //create an array for collecting my_waters. then copy to the vector
    int world_waters_ary[world_size];
    MPI_Allgather(&my_waters, 1, MPI_INT, world_waters_ary, 1, MPI_INT,MPI_COMM_WORLD);

    //allocate memory for world_waters and copy world_waters_ary
    world_waters.resize(world_size,0);
    for(i=0; i<world_size; i++)
    {
        world_waters[i] = world_waters_ary[i];
    }

    //allocate memory for world_water_start and world_water_end
    world_water_start.resize(world_size,0);
    world_water_end.resize(world_size,0);

    //compute each mpi processes start and end water
    distribute_work_load_lip(num_waters_1,world_rank,world_size,my_waters,world_waters,&water_start,&water_end,world_water_start,world_water_end);

    //print_work_load_stats(s.world_rank,s.world_size,world_waters,world_water_start,world_water_end);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints the workload distribution for waters                                                 //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::workload_water()
{
    print_work_load_stats_wat(world_rank,world_size,world_waters,world_water_start,world_water_end);
}


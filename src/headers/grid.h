
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects grid data (double) from all ranks and adds it to rank 0 grid data.                 //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void collect_grid_d(int world_size,int world_rank,vector<vector <double>> &data_vec)
{
    int i = 0;                            //Standard variable used in loops
    int j = 0;                            //Standard variable used in loops
    int k = 0;                            //Standard variable used in loops
    MPI_Status status;                    //Used for mpi_recv
    int num_g_y = data_vec.size();
    int num_g_x = data_vec[0].size();
    double weights[num_g_y][num_g_x];     //An array to hold the grid being collected
    double data[num_g_y][num_g_x];

    //copy vector data to an array. 
    for(i=0; i<num_g_x; i++) //loop over x
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            data[j][i] = data_vec[j][i];
        }
    }

    //collect array data
    if(world_size > 1)
    {
        for(i=1; i<world_size; i++)
        {
            if(world_rank == 0)
            {
                MPI_Recv(weights, num_g_x*num_g_y, MPI_DOUBLE, i, 13, MPI_COMM_WORLD, &status);
                for(j=0; j<num_g_x; j++) //loop over x
                { 
                    for(k=0; k<num_g_y; k++) //loop over y
                    {
                        data_vec[k][j] = data_vec[k][j] + weights[k][j];
                    }
                }
            }
            else if(world_rank == i)
            {
                MPI_Send(data, num_g_x*num_g_y, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects grid data (int) from all ranks and adds it to rank 0 grid data.                    //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void collect_grid_i(int world_size,int world_rank,vector<vector <int>> &data_vec)
{
    int i = 0;                            //Standard variable used in loops
    int j = 0;                            //Standard variable used in loops
    int k = 0;                            //Standard variable used in loops
    MPI_Status status;                    //Used for mpi_recv
    int num_g_y = data_vec.size();
    int num_g_x = data_vec[0].size();
    int weights[num_g_y][num_g_x];        //An array to hold the grid being collected
    int data[num_g_y][num_g_x];

    //copy vector data to an array.
    for(i=0; i<num_g_x; i++) //loop over x
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            data[j][i] = data_vec[j][i];
        }
    }

    //collect array data
    if(world_size > 1)
    {
        for(i=1; i<world_size; i++)
        {
            if(world_rank == 0)
            {
                MPI_Recv(weights, num_g_x*num_g_y, MPI_INT, i, 13, MPI_COMM_WORLD, &status);
                for(j=0; j<num_g_x; j++) //loop over x
                {
                    for(k=0; k<num_g_y; k++) //loop over y
                    {
                        data_vec[k][j] = data_vec[k][j] + weights[k][j];
                    }
                }
            }
            else if(world_rank == i)
            {
                MPI_Send(data, num_g_x*num_g_y, MPI_INT, 0, 13, MPI_COMM_WORLD);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
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
// This function gathers data (int) for the grid points. Used when parallelization scheme is by grid point   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void gather_grid_i_gp(int world_size,int world_rank,int my_num_g_x,int num_g_x,int num_g_y,iv1d &world_num_g_x,iv2d &local_data,iv2d &global_data)
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
                int data_rec[num_g_y][world_num_g_x[i]];

                //receive data
                MPI_Recv(data_rec, world_num_g_x[i]*num_g_y, MPI_INT, i, 13, MPI_COMM_WORLD, &status);

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
                int data_snd[num_g_y][my_num_g_x];
                for(j=0; j<my_num_g_x; j++) //loop over x
                {
                    for(k=0; k<num_g_y; k++) //loop over y
                    {
                        data_snd[k][j] = local_data[k][j];
                    }
                }

                MPI_Send(data_snd, my_num_g_x*num_g_y, MPI_INT, 0, 13, MPI_COMM_WORLD);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function broadcasts grid data (double) to all the ranks.                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void broadcast_grid_d(int world_size,int world_rank,vector<vector <double>> &data_vec)
{
    int i = 0;                            //Standard variable used in loops
    int j = 0;                            //Standard variable used in loops
    int k = 0;                            //Standard variable used in loops
    MPI_Status status;                    //Used for mpi_recv
    int num_g_y = data_vec.size();
    int num_g_x = data_vec[0].size();
    double weights[num_g_x][num_g_y];     //An array to hold the grid being collected
    double data[num_g_y][num_g_x];

    //copy vector data to an array.
    for(i=0; i<num_g_x; i++) //loop over x
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            data[j][i] = data_vec[j][i];
        }
    }

    //now we communicate the average rmsd to all nodes
    if(world_size > 1)
    {
        //broadcast rmsd
        for(i=1; i<world_size; i++)
        {
            if(world_rank == 0)
            {
                MPI_Send(data, num_g_x*num_g_y, MPI_DOUBLE, i, 13, MPI_COMM_WORLD);
            }
            else if(world_rank == i)
            {
                MPI_Recv(weights, num_g_x*num_g_y, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD, &status);
                for(j=0; j<num_g_x; j++) //loop over x
                {
                    for(k=0; k<num_g_y; k++) //loop over y
                    {
                        data_vec[k][j] = weights[k][j];
                    }
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function broadcasts grid data (double) to all the ranks.                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void broadcast_grid_i(int world_size,int world_rank,vector<vector <int>> &data_vec)
{
    int i = 0;                            //Standard variable used in loops
    int j = 0;                            //Standard variable used in loops
    int k = 0;                            //Standard variable used in loops
    MPI_Status status;                    //Used for mpi_recv
    int num_g_y = data_vec.size();
    int num_g_x = data_vec[0].size();
    int weights[num_g_x][num_g_y];     //An array to hold the grid being collected
    int data[num_g_y][num_g_x];

    //copy vector data to an array.
    for(i=0; i<num_g_x; i++) //loop over x
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            data[j][i] = data_vec[j][i];
        }
    }

    //now we communicate the average rmsd to all nodes
    if(world_size > 1)
    {
        //broadcast rmsd
        for(i=1; i<world_size; i++)
        {
            if(world_rank == 0)
            {
                MPI_Send(data, num_g_x*num_g_y, MPI_INT, i, 13, MPI_COMM_WORLD);
            }
            else if(world_rank == i)
            {
                MPI_Recv(weights, num_g_x*num_g_y, MPI_INT, 0, 13, MPI_COMM_WORLD, &status);
                for(j=0; j<num_g_x; j++) //loop over x
                {
                    for(k=0; k<num_g_y; k++) //loop over y
                    {
                        data_vec[k][j] = weights[k][j];
                    }
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function excludes insignificant data.                                                                //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exclude_insignificant_data(double cutoff,double avg_rho,vector<vector <double>> &rho,vector<vector <int>> &nan)
{
    int i=0;
    int j=0;

    int num_g_y = rho.size();
    int num_g_x = rho[0].size();

    for(i=0; i<num_g_y; i++) //loop over y-dimension
    {
        for(j=0; j<num_g_x; j++) //loop over x-dimension
        {
            if(rho[i][j] > cutoff*avg_rho)
            {
                nan[i][j] = 0;
            }
            else 
            {
                nan[i][j] = 1;
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes the grid to an output file.                                                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_grid_to_file_d(int out_data_format,double cell_size,string my_file_name,vector<vector <int>> &nan,vector<vector <double>> &data)
{
    int i=0;
    int j=0;
    int num_g_y = data.size();
    int num_g_x = data[0].size();

    FILE *my_file = fopen(my_file_name.c_str(), "w");
    if(my_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",my_file_name.c_str());
    }

    if(out_data_format == 0) //matrix
    {
        for(i=0; i<num_g_y; i++) //loop over y-dimension
        {
            for(j=0; j<num_g_x; j++) //loop over x-dimension
            {
                if(nan[i][j] == 0)
                {
                    fprintf(my_file," %10.6f",data[i][j]);
                }
                else //data excluded
                {
                    fprintf(my_file," %10s ","NaN");
                }
            }
            fprintf(my_file,"\n");
        }
    }
    else if(out_data_format == 1) //vector format
    {
        for(i=0; i<num_g_x; i++) //loop over x-dimension
        {
            for(j=0; j<num_g_y; j++) //loop over y-dimension
            {
                double x = i*cell_size;
                double y = j*cell_size;

                if(nan[i][j] == 0)
                {   
                    fprintf(my_file," %10.3f  %10.3f  %10.3f \n",x,y,data[j][i]);
                }
                else //data excluded
                {
                    fprintf(my_file," %10.3f  %10.3f  %10s \n",x,y,"NaN");
                }
            }
            fprintf(my_file,"\n");
        }
    }

    fclose(my_file);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes the grid to an output file.                                                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_grid_to_file_i(int out_data_format,double cell_size,string my_file_name,vector<vector <int>> &nan,vector<vector <int>> &data)
{
    int i=0;
    int j=0;
    int num_g_y = data.size();
    int num_g_x = data[0].size();

    FILE *my_file = fopen(my_file_name.c_str(), "w");
    if(my_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",my_file_name.c_str());
    }

    if(out_data_format == 0) //matrix
    {
        for(i=0; i<num_g_y; i++) //loop over y-dimension
        {
            for(j=0; j<num_g_x; j++) //loop over x-dimension
            {
                if(nan[i][j] == 0)
                {
                    fprintf(my_file," %10d",data[i][j]);
                }
                else //data excluded
                {
                    fprintf(my_file," %10s ","NaN");
                }
            }
            fprintf(my_file,"\n");
        }
    }
    else if(out_data_format == 1) //vector format
    {
        for(i=0; i<num_g_x; i++) //loop over x-dimension
        {
            for(j=0; j<num_g_y; j++) //loop over y-dimension
            {
                double x = i*cell_size;
                double y = j*cell_size;

                if(nan[i][j] == 0)
                {
                    fprintf(my_file," %10.3f  %10.3f  %10d \n",x,y,data[j][i]);
                }
                else //data excluded
                {
                    fprintf(my_file," %10.3f  %10.3f  %10s \n",x,y,"NaN");
                }
            }
            fprintf(my_file,"\n");
        }
    }

    fclose(my_file);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes the a single frame of the grid to an output file.                                    //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_grid_to_file_sf(int out_data_format,double cell_size,string my_file_name,vector<vector <int>> &count,vector<vector <double>> &data)
{
    int i=0;
    int j=0;
    int num_g_y = data.size();
    int num_g_x = data[0].size();

    FILE *my_file = fopen(my_file_name.c_str(), "w");
    if(my_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",my_file_name.c_str());
    }

    if(out_data_format == 0) //matrix
    {
        for(i=0; i<num_g_y; i++) //loop over y-dimension
        {
            for(j=0; j<num_g_x; j++) //loop over x-dimension
            {
                if(count[i][j] > 0) //data present
                {
                    fprintf(my_file," %10.6f",data[i][j]);
                }
                else //no data at grid point
                {
                    fprintf(my_file," %10s ","NaN");
                }
            }
            fprintf(my_file,"\n");
        }
    } 
    else if(out_data_format == 1) //vector format
    {
        for(i=0; i<num_g_x; i++) //loop over x-dimension
        {
            for(j=0; j<num_g_y; j++) //loop over y-dimension
            {
                double x = i*cell_size;
                double y = j*cell_size;

                if(count[j][i] > 0) //data present
                {
                    fprintf(my_file," %10.3f  %10.3f  %10f \n",x,y,data[j][i]);
                }
                else //data excluded
                {
                    fprintf(my_file," %10.3f  %10.3f  %10s \n",x,y,"NaN");
                }
            }
            fprintf(my_file,"\n");
        }
    }

 
    fclose(my_file);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the average rho.                                                                    //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double get_average_rho(vector <vector <double>> &rho)
{
    int i = 0;
    int j = 0;
    double avg_rho = 0;

    int num_g_y = rho.size();
    int num_g_x = rho[0].size();

    for(i=0; i<num_g_x; i++) //loop over x-dimension
    {
        for(j=0; j<num_g_y; j++) //loop over y-dimension
        {
            avg_rho = avg_rho + rho[j][i];
        }
    }
    avg_rho = avg_rho/(double)(num_g_x*num_g_y);

    return avg_rho;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function divides each gridpoint by rho.                                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void normalize_the_grid(vector<vector <double>> &rho,vector<vector <double>> &data)
{
    int i = 0;
    int j = 0;
    int num_g_y = data.size();
    int num_g_x = data[0].size();

    for(i=0; i<num_g_x; i++) //loop over x-dimension
    {
        for(j=0; j<num_g_y; j++) //loop over y-dimension
        {
            if(rho[j][i] > 0)
            {
                data[j][i] = data[j][i]/rho[j][i];
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function determines the upper and lower bounds surrounding the mapping atom center.                  //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_bounds(double hx,double hy,double radius,double cell_size,int num_g_x,int num_g_y,int *lower_x,int *lower_y,int *upper_x,int *upper_y)
{   
    //calculate the upper and lower bounds
    *lower_x = floor((hx - radius)/cell_size);
    *lower_y = floor((hy - radius)/cell_size);
    *upper_x = ceil((hx + radius)/cell_size);
    *upper_y = ceil((hy + radius)/cell_size);

    //add a 1 point buffer around the paremeter
    *lower_x = *lower_x - 1;
    *lower_y = *lower_y - 1;
    *upper_x = *upper_x + 1;
    *upper_y = *upper_y + 1; 
    
    //make sure the bounds do not fall outside the grid
    //lower_x
    if(*lower_x < 0)
    {   
        *lower_x = 0;
    }
    else if(*lower_x > num_g_x)
    {   
        *lower_x = num_g_x - 1;
    }
    
    //lower_y
    if(*lower_y < 0)
    {   
        *lower_y = 0;
    }
    else if(*lower_y > num_g_y)
    {   
        *lower_y = num_g_y - 1;
    }
    
    //upper_x
    if(*upper_x > num_g_x)
    {   
        *upper_x = num_g_x;
    }
    else if(*upper_x < 0)
    {   
        *upper_x = 0;
    }
    
    //upper_y
    if(*upper_y > num_g_y)
    {   
        *upper_y = num_g_y;
    }
    else if(*upper_y < 0)
    {   
        *upper_y = 0;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function deposits a double to the grid.                                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void add_to_grid_d(double hx,double hy,double radius,double cell_size,vector<vector <double>> &grid,double data)
{
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int lower_x = 0;
    int lower_y = 0;
    int upper_x = 0;
    int upper_y = 0;
    int num_g_y = grid.size();
    int num_g_x = grid[0].size();

    get_bounds(hx,hy,radius,cell_size,num_g_x,num_g_y,&lower_x,&lower_y,&upper_x,&upper_y);

    for(i=lower_x; i<upper_x; i++) //loop over the grid x-axis
    {
        for(j=lower_y; j<upper_y; j++) //loop over the grid y-axis
        {
            double dist_x = i*cell_size - hx;
            double dist_y = j*cell_size - hy;

            double dist = sqrt(dist_x*dist_x + dist_y*dist_y);

            if(dist <= radius)
            {
                grid[j][i] = grid[j][i] + data;
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function deposits a double to the grid checking check first                                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void controlled_add_to_grid_d(double hx,double hy,double radius,double cell_size,vector<vector <double>> &grid,double data,vector<vector <int>> &check)
{
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int lower_x = 0;
    int lower_y = 0;
    int upper_x = 0;
    int upper_y = 0;
    int num_g_y = grid.size();
    int num_g_x = grid[0].size();

    get_bounds(hx,hy,radius,cell_size,num_g_x,num_g_y,&lower_x,&lower_y,&upper_x,&upper_y);

    for(i=lower_x; i<upper_x; i++) //loop over the grid x-axis
    {
        for(j=lower_y; j<upper_y; j++) //loop over the grid y-axis
        {
            if(check[j][i] == 0)
            {
                double dist_x = i*cell_size - hx;
                double dist_y = j*cell_size - hy;

                double dist = sqrt(dist_x*dist_x + dist_y*dist_y);

                if(dist <= radius)
                {
                    grid[j][i] = grid[j][i] + data;
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes a double to the grid.                                                                //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_to_grid_d(double hx,double hy,double radius,double cell_size,vector<vector <double>> &grid,double data)
{
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int lower_x = 0;
    int lower_y = 0;
    int upper_x = 0;
    int upper_y = 0;
    int num_g_y = grid.size();
    int num_g_x = grid[0].size();

    get_bounds(hx,hy,radius,cell_size,num_g_x,num_g_y,&lower_x,&lower_y,&upper_x,&upper_y);

    for(i=lower_x; i<upper_x; i++) //loop over the grid x-axis
    {
        for(j=lower_y; j<upper_y; j++) //loop over the grid y-axis
        {
            double dist_x = i*cell_size - hx;
            double dist_y = j*cell_size - hy;

            double dist = sqrt(dist_x*dist_x + dist_y*dist_y);

            if(dist <= radius)
            {
                grid[j][i] = data;
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function deposits an int to the grid.                                                                //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void add_to_grid_i(double hx,double hy,double radius,double cell_size,vector<vector <int>> &grid,int data)
{
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int lower_x = 0;
    int lower_y = 0;
    int upper_x = 0;
    int upper_y = 0;
    int num_g_y = grid.size();
    int num_g_x = grid[0].size();

    get_bounds(hx,hy,radius,cell_size,num_g_x,num_g_y,&lower_x,&lower_y,&upper_x,&upper_y);

    for(i=lower_x; i<upper_x; i++) //loop over the grid x-axis
    {
        for(j=lower_y; j<upper_y; j++) //loop over the grid y-axis
        {
            double dist_x = i*cell_size - hx;
            double dist_y = j*cell_size - hy;

            double dist = sqrt(dist_x*dist_x + dist_y*dist_y);

            if(dist <= radius)
            {
                grid[j][i] = grid[j][i] + data;
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function deposits an int to the grid checking check first.                                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void controlled_add_to_grid_i(double hx,double hy,double radius,double cell_size,vector<vector <int>> &grid,int data,vector<vector <int>> &check)
{
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int lower_x = 0;
    int lower_y = 0;
    int upper_x = 0;
    int upper_y = 0;
    int num_g_y = grid.size();
    int num_g_x = grid[0].size();

    get_bounds(hx,hy,radius,cell_size,num_g_x,num_g_y,&lower_x,&lower_y,&upper_x,&upper_y);

    for(i=lower_x; i<upper_x; i++) //loop over the grid x-axis
    {
        for(j=lower_y; j<upper_y; j++) //loop over the grid y-axis
        {
            if(check[j][i] == 0)
            {
                double dist_x = i*cell_size - hx;
                double dist_y = j*cell_size - hy;

                double dist = sqrt(dist_x*dist_x + dist_y*dist_y);

                if(dist <= radius)
                {
                    grid[j][i] = grid[j][i] + data;
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes an int to the grid.                                                                  //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_to_grid_i(double hx,double hy,double radius,double cell_size,vector<vector <int>> &grid,int data)
{
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int lower_x = 0;
    int lower_y = 0;
    int upper_x = 0;
    int upper_y = 0;
    int num_g_y = grid.size();
    int num_g_x = grid[0].size();

    get_bounds(hx,hy,radius,cell_size,num_g_x,num_g_y,&lower_x,&lower_y,&upper_x,&upper_y);

    for(i=lower_x; i<upper_x; i++) //loop over the grid x-axis
    {
        for(j=lower_y; j<upper_y; j++) //loop over the grid y-axis
        {
            double dist_x = i*cell_size - hx;
            double dist_y = j*cell_size - hy;

            double dist = sqrt(dist_x*dist_x + dist_y*dist_y);

            if(dist <= radius)
            {
                grid[j][i] = data;
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function deposits lipid density to the grid.                                                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void deposit_lipid_density(double hx,double hy,double radius,double cell_size,dv2d &grid,int frames,int num_atoms_deposit)
{
    int i        = 0;    //standard variable used in loops
    int j        = 0;    //standard variable used in loops
    int lower_x  = 0;    //lower index flanking atom in x
    int lower_y  = 0;    //lower index flanking atom in y
    int upper_x  = 0;    //upper index flanking atom in x
    int upper_y  = 0;    //upper index flanking atom in y
    double count = 0;    //how many grid points does atom add to

    int num_g_y = grid.size();
    int num_g_x = grid[0].size();

    get_bounds(hx,hy,radius,cell_size,num_g_x,num_g_y,&lower_x,&lower_y,&upper_x,&upper_y);

    //compute the normalizing factor
    for(i=lower_x; i<upper_x; i++) //loop over the grid x-axis
    {
        for(j=lower_y; j<upper_y; j++) //loop over the grid y-axis
        {
            double dist_x = i*cell_size - hx;
            double dist_y = j*cell_size - hy;

            double dist = sqrt(dist_x*dist_x + dist_y*dist_y);

            if(dist <= radius)
            {
                count = count + 1.0;
            }
        }
    }

    //deposit density to the grid
    for(i=lower_x; i<upper_x; i++) //loop over the grid x-axis
    {
        for(j=lower_y; j<upper_y; j++) //loop over the grid y-axis
        {
            double dist_x = i*cell_size - hx;
            double dist_y = j*cell_size - hy;

            double dist = sqrt(dist_x*dist_x + dist_y*dist_y);

            if(dist <= radius) 
            {       
                grid[j][i] = grid[j][i] + 1.0/(count*(double)frames*(double)num_atoms_deposit);
            }       
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints information about the grid.                                                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void print_grid_stats(double box_x,double box_y,matrix ibox,int num_g_x,int num_g_y,double cell_size,int world_rank)
{
    if(world_rank == 0)
    {
        if(box_x == 0 && box_y == 0) //user did not specify a box size
        {
            printf("Using box at frame zero to construct the grid as follows: \n");
            printf(" %10s %10s %10s %10s %10s \n","box_x","box_y","cell_size","num_g_x","num_g_y");
            printf("-%10s-%10s-%10s-%10s-%10s \n","----------","----------","----------","----------","----------");
            printf(" %10f %10f %10f %10d %10d \n",ibox[XX][XX],ibox[YY][YY],cell_size,num_g_x,num_g_y);
        }
        else
        {
            printf("Using user specified box to construct the grid as follows: \n");
            printf(" %10s %10s %10s %10s %10s \n","box_x","box_y","cell_size","num_g_x","num_g_y");
            printf("-%10s-%10s-%10s-%10s-%10s \n","----------","----------","----------","----------","----------");
            printf(" %10f %10f %10f %10d %10d \n",box_x,box_y,cell_size,num_g_x,num_g_y);
        }
        printf("\n");
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints information about the grid.                                                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void print_grid_stats_obj(double box_x,double box_y,int b_box_i,int num_g_x,int num_g_y,double cell_size,int world_rank)
{
    if(world_rank == 0)
    {
        if(b_box_i == 1) //user did not specify a box size
        {
            printf("Using box at frame zero to construct the grid as follows: \n");
        }
        else //user specified a box to use
        {
            printf("Using user specified box to construct the grid as follows: \n");
        }
        printf(" %10s %10s %10s %10s %10s \n","box_x","box_y","cell_size","num_g_x","num_g_y");
        printf("-%10s-%10s-%10s-%10s-%10s \n","----------","----------","----------","----------","----------");
        printf(" %10f %10f %10f %10d %10d \n",box_x,box_y,cell_size,num_g_x,num_g_y);

        printf("\n");
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function determines how many grid points in the x and y dimensions.                                  //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_grid_size(double box_x,double box_y,matrix ibox,int *num_g_x,int *num_g_y,double cell_size)
{
    if(box_x == 0 && box_y == 0) //user did not specify a box size
    {
        *num_g_x = (int)ceil(ibox[XX][XX]/cell_size)+1;
        *num_g_y = (int)ceil(ibox[YY][YY]/cell_size)+1;
    }
    else
    {
        *num_g_x = (int)ceil(box_x/cell_size)+1;
        *num_g_y = (int)ceil(box_y/cell_size)+1;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function determines how many grid points in the x and y dimensions.                                  //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_grid_size_obj(double my_box_x,double my_box_y,matrix ibox,int *num_g_x,int *num_g_y,double cell_size,double *box_x,double *box_y,int *b_box_i)
{
    if(my_box_x == 0 && my_box_y == 0) //user did not specify a box size
    {
        *num_g_x = (int)ceil(ibox[XX][XX]/cell_size)+1;
        *num_g_y = (int)ceil(ibox[YY][YY]/cell_size)+1;
        *box_x = ibox[YY][YY];
        *box_y = ibox[YY][YY];
        *b_box_i = 1;
    }
    else
    {
        *num_g_x = (int)ceil(my_box_x/cell_size)+1;
        *num_g_y = (int)ceil(my_box_y/cell_size)+1;
        *box_x = my_box_x;
        *box_y = my_box_y;
        *b_box_i = 0;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function sizes and initializes the grid.                                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void size_and_init_grid_d(int num_g_x,int num_g_y,vector<vector <double>> &grid)
{
    int i = 0;           //standard variable used in loops
    int j = 0;           //standard variable used in loops

    //set the size of the vectors
    grid.resize(num_g_y);

    for(i=0; i<num_g_y; i++) //loop over y
    {
        grid[i].resize(num_g_x);
    }

    //initialize the grid
    for(i=0; i<grid.size(); i++) //loop over lipid atoms
    {
        for(j=0; j<grid[i].size(); j++) //loop over grid for each atom
        {
            grid[i][j] = 0;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function sizes and initializes the grid.                                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void size_and_init_grid_i(int num_g_x,int num_g_y,vector<vector <int>> &grid)
{   
    int i = 0;           //standard variable used in loops
    int j = 0;           //standard variable used in loops
    
    //set the size of the vectors
    grid.resize(num_g_y);
    
    for(i=0; i<num_g_y; i++) //loop over y
    {   
        grid[i].resize(num_g_x);
    }
    
    //initialize the grid
    for(i=0; i<grid.size(); i++) //loop over lipid atoms
    {   
        for(j=0; j<grid[i].size(); j++) //loop over grid for each atom
        {   
            grid[i][j] = 0;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This is a class for working with the grid                                                                //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Grid
{
    private:
        string grid_file_name;                                                               //name of the file to write grid to
        string rho_file_name;                                                                //name of the file to write rho to
        int    i           = 0;                                                              //standard variable used in loops
        int    j           = 0;                                                              //standard variable used in loops
        int    world_size  = 0;                                                              //size of mpi world
        int    world_rank  = 0;                                                              //rank of mpi process
        int    odf         = 0;                                                              //format for output files
        double APS         = 0.005;                                                          //area of a grid square
        double cell_size   = sqrt(APS);                                                      //grid point spacing
        int    num_g_x     = 1;                                                              //number of grid points in x
        int    num_g_y     = 1;                                                              //number of grid points in y
        int    b_box_i     = 1;                                                              //is the initial box used?
        double size_grid_x = 0;                                                              //size of grid in x
        double size_grid_y = 0;                                                              //size of grid in y
        double box_x       = 0;                                                              //box x
        double box_y       = 0;                                                              //box y

    public:
        dv2d grid{};                                                                         //long term grid data       
        dv2d rho{};                                                                          //long term grid rho data
        dv2d frame_grid{};                                                                   //current frame grid 
        iv2d cell_count{};                                                                   //how many atoms added to current frame grid
        iv2d nan{};                                                                          //tells if the grid point should be excluded
        dv2d stdev{};                                                                        //standard deviation of grid data    
         
    public:
        void   get_dim(double box_x,double box_y,matrix ibox,double my_APS);                 //determine the grid dimensions
        void   print_dim();                                                                  //print the grid dimensions
        void   set_output(string my_grid_file_name,int odf);                                 //set the name and format of the output file
        void   clean_frame();                                                                //clean current frame grid
        void   clean_this_grid(dv2d this_grid);                                              //clears the provided grid 
        void   stamp(double hx,double hy,double radius,double data);                         //add a value to current frame grid
        void   direct_stamp(double hx,double hy,double radius,double data);                  //stamp directly to long term grid
        void   norm_frame();                                                                 //normalize the current frame
        void   add_frame();                                                                  //add the current frame to long term grid 
        void   add_frame_direct();                                                           //add the current frame to long term grid. add count to rho
        void   write_frame(int global_frame);                                                //write the current frame to out put file
        void   collect_grid();                                                               //collect grid/rho from mpi processes
        void   normalize();                                                                  //normalize the grid
        void   exclude_data(double cutoff,int report);                                       //exclude insignificant data
        void   write_grid();                                                                 //write the gid to file
        void   write_rho();                                                                  //write rho to file
        void   bcast_grid();                                                                 //broadcast the grid to all mpi processes
        void   get_stdev(int b_stdev,int b_clean,Trajectory traj);                           //compute the standard deviation
        void   match_grid(Grid model_grid);                                                  //set the grid dimensions to match another grid
        void   copy_rho(Grid model_rho);                                                     //copy rho from another grid
        int    num_x();                                                                      //return num_g_x
        int    num_y();                                                                      //return num_g_y
        double get_cell_size();                                                              //return cell_size
        double get_aps();                                                                    //return APS
        double get_box_x();                                                                  //return the box x dimension
        double get_box_y();                                                                  //return the box y dimension
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return determines the grid dimensions                                                      //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::get_dim(double my_box_x,double my_box_y,matrix ibox,double my_APS)
{
    //get the world size and assign ranks
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);      //get the world size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);      //get the process rank

    APS       = my_APS;
    cell_size = sqrt(APS);

    get_grid_size_obj(my_box_x,my_box_y,ibox,&num_g_x,&num_g_y,cell_size,&box_x,&box_y,&b_box_i);

    size_grid_x = num_g_x*cell_size;
    size_grid_y = num_g_y*cell_size;

    size_and_init_grid_d(num_g_x,num_g_y,grid);
    size_and_init_grid_d(num_g_x,num_g_y,rho);
    size_and_init_grid_d(num_g_x,num_g_y,frame_grid);
    size_and_init_grid_i(num_g_x,num_g_y,cell_count);
    size_and_init_grid_i(num_g_x,num_g_y,nan);
}  

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function prints the grid dimensions                                                                 //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::print_dim()
{
    print_grid_stats_obj(box_x,box_y,b_box_i,num_g_x,num_g_y,cell_size,world_rank);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function sets the output name and format                                                            //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::set_output(string my_grid_file_name,int my_odf)
{
    grid_file_name = my_grid_file_name;

    rho_file_name = add_tag(grid_file_name,"_rho");

    odf = my_odf;
}   

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return num_g_x                                                                             //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Grid::num_x()
{
    return num_g_x;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return num_g_y                                                                             //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Grid::num_y()
{
    return num_g_y;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns cell_size                                                                          //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Grid::get_cell_size()
{
    return cell_size;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return APS                                                                                 //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Grid::get_aps()
{
    return APS;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function deposits a quantity to the current frame grid                                              //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::stamp(double hx,double hy,double radius,double data)
{
   add_to_grid_d(hx,hy,radius,cell_size,frame_grid,data); 
   add_to_grid_i(hx,hy,radius,cell_size,cell_count,1);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function deposits a quantity to the long term grid                                                  //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::direct_stamp(double hx,double hy,double radius,double data)
{
   add_to_grid_d(hx,hy,radius,cell_size,grid,data);
   add_to_grid_d(hx,hy,radius,cell_size,rho,1.0);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function clears the current frame grid                                                              //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::clean_frame()
{
    for(i=0; i<num_g_x; i++) //loop over x
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            frame_grid[j][i] = 0;
            cell_count[j][i] = 0;
        }
    }
}  

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function clears the entire grid                                                                     //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::clean_this_grid(dv2d this_grid)
{
    for(i=0; i<num_g_x; i++) //loop over x
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            this_grid[j][i]       = 0.0;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function normalizes the current frame                                                               //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::norm_frame()
{
    for(i=0; i<num_g_x; i++) //loop over x
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            if(cell_count[j][i] > 0)
            {
                frame_grid[j][i] = frame_grid[j][i]/(double)cell_count[j][i];
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function adds the current frame grid to the long long term sum                                      //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::add_frame()
{
    for(i=0; i<num_g_x; i++) //loop over x
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            if(cell_count[j][i] > 0)
            {
                grid[j][i] = grid[j][i] + frame_grid[j][i];
                rho[j][i]  = rho[j][i]  + 1.0;
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function adds the current frame grid to the long long term sum. add count to rho                    //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::add_frame_direct()
{
    for(i=0; i<num_g_x; i++) //loop over x
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            if(cell_count[j][i] > 0)
            {
                grid[j][i] = grid[j][i] + frame_grid[j][i];
                rho[j][i]  = rho[j][i]  + cell_count[j][i];
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function writes the current frame to file                                                           //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::write_frame(int global_frame)
{
    string frame_file_name = add_tag(grid_file_name,to_string(global_frame));

    write_grid_to_file_sf(odf,cell_size,frame_file_name,cell_count,frame_grid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function collects the grid from mpi processes                                                       //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::collect_grid()
{
    collect_grid_d(world_size,world_rank,grid);
    collect_grid_d(world_size,world_rank,rho);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function normalizes the grid                                                                        //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::normalize()
{
    normalize_the_grid(rho,grid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function excludes insignificant data                                                                //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::exclude_data(double cutoff,int report)
{
    double avg_rho = get_average_rho(rho);

    if(report == 1)
    {
        printf("The average rho is %f. Excluding lattice points with fewer than %f samples. \n",avg_rho,cutoff*avg_rho);
    }
    exclude_insignificant_data(cutoff,avg_rho,rho,nan);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function writes the grid to output file                                                             //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::write_grid()                                                                                      
{
    write_grid_to_file_d(odf,cell_size,grid_file_name,nan,grid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function writes rho to output file                                                                  //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::write_rho()
{
    write_grid_to_file_d(odf,cell_size,rho_file_name,nan,rho);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function broadcasts the grid                                                                        //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::bcast_grid()
{
    broadcast_grid_d(world_size,world_rank,grid);
    broadcast_grid_d(world_size,world_rank,rho);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function computes the standard deviation                                                            //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::get_stdev(int b_stdev,int b_clean,Trajectory traj)
{
    if(b_stdev == 1)
    {
        int i = 0;                               //standard variable used in loops
        int j = 0;                               //standard variable used in loops
        int k = 0;                               //standard variable used in loops
        char my_string[200];                     //Used for reading data from single frame files

        size_and_init_grid_d(num_g_x,num_g_y,stdev);

        if(world_rank == 0)
        {
            printf("Computing the standard deviation \n");
        }

        //broadcast <data> to all ranks
        broadcast_grid_d(world_size,world_rank,grid);

        //now each rank computes the sum of squares
        for(traj.current_frame=0; traj.current_frame<traj.get_num_frames(); traj.current_frame++)
        {
            int global_frame = traj.get_frame_global();

            string single_file_name = add_tag(grid_file_name,to_string(global_frame));
 
            FILE *single_file = fopen(single_file_name.c_str(), "r");
            if(single_file == NULL)
            {
                printf("failure opening %s. Make sure the file exists. \n",single_file_name.c_str());
            }
            else
            {
                for(k=0; k<num_g_y; k++) //loop over y-dimension
                {
                    for(j=0; j<num_g_x; j++) //loop over x-dimension
                    {
                       int result = fscanf(single_file, "%s,", my_string);

                       if(strcmp(my_string, "NaN") != 0) //data present
                       {
                           double value  = atof(my_string);
                           double dif    = value - grid[k][j];
                           double square = dif*dif;
                           stdev[k][j]   = stdev[k][j] + square;
                        }
                    }
                }
                fclose(single_file);
            }
      
            //clean up single frame files
            if(b_clean == 1)
            {
                remove(single_file_name.c_str());
            }
        }

        //collect the sum_of_squares from all ranks
        collect_grid_d(world_size,world_rank,stdev);

        //now get the final standard deviation
        if(world_rank == 0)
        {
            for(k=0; k<num_g_y; k++) //loop over y-dimension
            {
                for(j=0; j<num_g_x; j++) //loop over x-dimension
                {
                    if(rho[k][j] > 0)
                    {
                        double norm        = stdev[k][j]/(rho[k][j]-1.0);
                        double stdeviation = sqrt(norm);

                        stdev[k][j] = stdeviation;
                    }
                }
            }

            //set the stdev file name
            string stdev_file_name = add_tag(grid_file_name,"_stdev");

            //now print the standard deviation to file
            write_grid_to_file_d(odf,cell_size,stdev_file_name,nan,stdev);

            //compute standard error
            for(k=0; k<num_g_y; k++) //loop over y-dimension
            {
                for(j=0; j<num_g_x; j++) //loop over x-dimension
                {
                    if(rho[k][j] > 0)
                    {
                        stdev[k][j] = stdev[k][j]/sqrt(rho[k][j]);
                    }
                }
            }

            //set the stem file name
            string stem_file_name = add_tag(grid_file_name,"_stem");

            //now print the standard error to file
            write_grid_to_file_d(odf,cell_size,stem_file_name,nan,stdev);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns the box x dimension                                                                //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Grid::get_box_x()
{
    return box_x;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns the box y dimension                                                                //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Grid::get_box_y()
{
    return box_y;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function copies the dimesion of a grid to another grid                                              //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::match_grid(Grid model_grid)
{
    //get the world size and assign ranks
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);      //get the world size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);      //get the process rank

    num_g_x   = model_grid.num_x();
    num_g_y   = model_grid.num_y();
    APS       = model_grid.get_aps();
    cell_size = model_grid.get_cell_size();

    box_x = get_box_x();
    box_y = get_box_y();

    size_grid_x = num_g_x*cell_size;
    size_grid_y = num_g_y*cell_size;

    size_and_init_grid_d(num_g_x,num_g_y,grid);
    size_and_init_grid_d(num_g_x,num_g_y,rho);
    size_and_init_grid_d(num_g_x,num_g_y,frame_grid);
    size_and_init_grid_i(num_g_x,num_g_y,cell_count);
    size_and_init_grid_i(num_g_x,num_g_y,nan);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function copies rho from another grid                                                               //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::copy_rho(Grid model_rho)
{
    for(i=0; i<num_g_x; i++) //loop over x
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            rho[j][i] = model_rho.rho[j][i];
        }
    }
}





//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This is a class for working with a single grid of ints                                                   //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Grid_i
{
    private:
        string grid_file_name;                                                                //name of the file to write grid to
        int    i           = 0;                                                               //standard variable used in loops
        int    j           = 0;                                                               //standard variable used in loops
        int    world_size  = 0;                                                               //size of mpi world
        int    world_rank  = 0;                                                               //rank of mpi process
        int    odf         = 0;                                                               //format for output files
        double APS         = 0.005;                                                           //area of a grid square
        double cell_size   = sqrt(APS);                                                       //grid point spacing
        int    num_g_x     = 1;                                                               //number of grid points in x
        int    num_g_y     = 1;                                                               //number of grid points in y
        double size_grid_x = 0;                                                               //size of grid in x
        double size_grid_y = 0;                                                               //size of grid in y

    public:
        iv2d grid{};                                                                          //long term grid data

    public:
        void   set_dim(double my_APS,double my_num_g_x,double my_num_g_y);                    //sets the grid dimensions
        void   set_output(string my_grid_file_name,int odf);                                  //set the name and format of the output file
        void   clean_grid();                                                                  //clean the grid
        void   stamp(double hx,double hy,double radius,int data);                             //add a value to the grid
        void   controlled_stamp(double hx,double hy,double radius,int data,Grid_i &check);    //add a value to the grid if check equals zero.
        void   paint(double hx,double hy,double radius,int data);                             //set the value for the grid points
        void   collect_grid();                                                                //collect grid/rho from mpi processes
        void   write_grid(Grid_i nan);                                                        //write the gid to file
        void   bcast_grid();                                                                  //broadcast the grid to all mpi processes
        void   copy_grid(Grid_i model_grid);                                                  //copy from another grid
        int    num_x();                                                                       //return num_g_x
        int    num_y();                                                                       //return num_g_y
        double get_cell_size();                                                               //return cell_size
        double get_aps();                                                                     //return APS
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return determines the grid dimensions                                                      //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_i::set_dim(double my_APS,double my_num_g_x,double my_num_g_y)
{
    //get the world size and assign ranks
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);      //get the world size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);      //get the process rank

    APS       = my_APS;
    cell_size = sqrt(APS);
    num_g_x   = my_num_g_x;
    num_g_y   = my_num_g_y;

    size_grid_x = num_g_x*cell_size;
    size_grid_y = num_g_y*cell_size;

    size_and_init_grid_i(num_g_x,num_g_y,grid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function sets the output name and format                                                            //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_i::set_output(string my_grid_file_name,int my_odf)
{
    grid_file_name = my_grid_file_name;

    odf = my_odf;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function clears the grid                                                                            //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_i::clean_grid()
{
    for(i=0; i<num_g_x; i++) //loop over x
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            grid[j][i] = 0;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function deposits a quantity to grid                                                                //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_i::stamp(double hx,double hy,double radius,int data)
{
   add_to_grid_i(hx,hy,radius,cell_size,grid,data);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function writes a quantity to grid                                                                  //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_i::paint(double hx,double hy,double radius,int data)
{
   write_to_grid_i(hx,hy,radius,cell_size,grid,data);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function collects the grid from mpi processes                                                       //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_i::collect_grid()
{
    collect_grid_i(world_size,world_rank,grid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function writes the grid to output file                                                             //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_i::write_grid(Grid_i nan)
{
    write_grid_to_file_i(odf,cell_size,grid_file_name,nan.grid,grid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function broadcasts the grid                                                                        //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_i::bcast_grid()
{
    broadcast_grid_i(world_size,world_rank,grid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function copies data from another grid                                                              //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_i::copy_grid(Grid_i model_grid)
{
    for(i=0; i<num_g_x; i++) //loop over x
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            grid[j][i] = model_grid.grid[j][i];
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return num_g_x                                                                             //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Grid_i::num_x()
{
    return num_g_x;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return num_g_y                                                                             //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Grid_i::num_y()
{
    return num_g_y;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns cell_size                                                                          //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Grid_i::get_cell_size()
{
    return cell_size;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return APS                                                                                 //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Grid_i::get_aps()
{
    return APS;
}







//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This is a class for working with a single grid of doubles                                                //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Grid_d
{
    private:
        string grid_file_name;                                                                //name of the file to write grid to
        int    i           = 0;                                                               //standard variable used in loops
        int    j           = 0;                                                               //standard variable used in loops
        int    world_size  = 0;                                                               //size of mpi world
        int    world_rank  = 0;                                                               //rank of mpi process
        int    odf         = 0;                                                               //format for output files
        double APS         = 0.005;                                                           //area of a grid square
        double cell_size   = sqrt(APS);                                                       //grid point spacing
        int    num_g_x     = 1;                                                               //number of grid points in x
        int    num_g_y     = 1;                                                               //number of grid points in y
        double size_grid_x = 0;                                                               //size of grid in x
        double size_grid_y = 0;                                                               //size of grid in y

    public:
        dv2d grid{};                                                                          //long term grid data

    public:
        void   set_dim(double my_APS,double my_num_g_x,double my_num_g_y);                    //sets the grid dimensions
        void   set_output(string my_grid_file_name,int odf);                                  //set the name and format of the output file
        void   clean_grid();                                                                  //clean the grid
        void   stamp(double hx,double hy,double radius,double data);                          //add a value to the grid
        void   controlled_stamp(double hx,double hy,double radius,double data,Grid_i &check); //add a value to the grid if check equals zero. 
        void   collect_grid();                                                                //collect grid/rho from mpi processes
        void   paint(double hx,double hy,double radius,double data);                          //set the value for the grid points
        void   normalize(Grid_d rho);                                                         //normalize the grid
        void   write_grid(Grid_i nan);                                                        //write the gid to file
        void   bcast_grid();                                                                  //broadcast the grid to all mpi processes
        void   copy_grid(Grid_d model_grid);                                                  //copy from another grid
        int    num_x();                                                                       //return num_g_x
        int    num_y();                                                                       //return num_g_y
        double get_cell_size();                                                               //return cell_size
        double get_aps();                                                                     //return APS
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return determines the grid dimensions                                                      //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_d::set_dim(double my_APS,double my_num_g_x,double my_num_g_y)
{
    //get the world size and assign ranks
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);      //get the world size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);      //get the process rank

    APS       = my_APS;
    cell_size = sqrt(APS);
    num_g_x   = my_num_g_x;
    num_g_y   = my_num_g_y;

    size_grid_x = num_g_x*cell_size;
    size_grid_y = num_g_y*cell_size;

    size_and_init_grid_d(num_g_x,num_g_y,grid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function sets the output name and format                                                            //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_d::set_output(string my_grid_file_name,int my_odf)
{
    grid_file_name = my_grid_file_name;

    odf = my_odf;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function clears the grid                                                                            //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_d::clean_grid()
{
    for(i=0; i<num_g_x; i++) //loop over x
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            grid[j][i] = 0;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function deposits a quantity to grid                                                                //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_d::stamp(double hx,double hy,double radius,double data)
{
   add_to_grid_d(hx,hy,radius,cell_size,grid,data);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function deposits a quantity to grid checking first check                                           //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_d::controlled_stamp(double hx,double hy,double radius,double data,Grid_i &check)
{
   controlled_add_to_grid_d(hx,hy,radius,cell_size,grid,data,check.grid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function writes a quantity to grid                                                                  //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_d::paint(double hx,double hy,double radius,double data)
{
   write_to_grid_d(hx,hy,radius,cell_size,grid,data);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function collects the grid from mpi processes                                                       //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_d::collect_grid()
{
    collect_grid_d(world_size,world_rank,grid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function normalizes the grid                                                                        //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_d::normalize(Grid_d rho)
{
    for(i=0; i<num_g_x; i++) //loop over x
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            if(rho.grid[j][i] > 0)
            {
                grid[j][i] = grid[j][i]/rho.grid[j][i];
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function writes the grid to output file                                                             //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_d::write_grid(Grid_i nan)
{
    write_grid_to_file_d(odf,cell_size,grid_file_name,nan.grid,grid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function broadcasts the grid                                                                        //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_d::bcast_grid()
{
    broadcast_grid_d(world_size,world_rank,grid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function copies rho from another grid                                                               //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_d::copy_grid(Grid_d model_grid)
{
    for(i=0; i<num_g_x; i++) //loop over x
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            grid[j][i] = model_grid.grid[j][i];
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return num_g_x                                                                             //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Grid_d::num_x()
{
    return num_g_x;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return num_g_y                                                                             //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Grid_d::num_y()
{
    return num_g_y;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns cell_size                                                                          //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Grid_d::get_cell_size()
{
    return cell_size;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return APS                                                                                 //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Grid_d::get_aps()
{
    return APS;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function excludes insignificant data                                                                //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exclude_data(Grid_d rho,Grid_i &nan,double cutoff,int report)
{
    double avg_rho = get_average_rho(rho.grid);

    if(report == 1)
    {
        printf("The average rho is %f. Excluding lattice points with fewer than %f samples. \n",avg_rho,cutoff*avg_rho);
    }
    exclude_insignificant_data(cutoff,avg_rho,rho.grid,nan.grid);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects grid data (double) from all ranks and adds it to rank 0 grid data.                 //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void collect_grid_d(int world_size,int world_rank,vector<vector<vector <double>>> &data_vec)
{
    int i = 0;                            //Standard variable used in loops
    int j = 0;                            //Standard variable used in loops
    int k = 0;                            //Standard variable used in loops
    int l = 0;                            //Standard variable used in loops
    MPI_Status status;                    //Used for mpi_recv
    int num_g_z = data_vec.size();
    int num_g_y = data_vec[0].size();
    int num_g_x = data_vec[0][0].size();
    double weights[num_g_z][num_g_y][num_g_x];     //An array to hold the grid being collected
    double data[num_g_z][num_g_y][num_g_x];

    //copy vector data to an array. 
    for(i=0; i<num_g_z; i++) //loop over z
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            for(k=0; k<num_g_x; k++) //loop over x
            {
                data[i][j][k] = data_vec[i][j][k];
            }
        }
    }

    //collect array data
    if(world_size > 1)
    {
        for(i=1; i<world_size; i++)
        {
            if(world_rank == 0)
            {
                MPI_Recv(weights, num_g_z*num_g_y*num_g_x, MPI_DOUBLE, i, 13, MPI_COMM_WORLD, &status);
                for(j=0; j<num_g_z; j++) //loop over z
                { 
                    for(k=0; k<num_g_y; k++) //loop over y
                    {
                        for(l=0; l<num_g_x; l++) //loop over x
                        {
                            data_vec[j][k][l] = data_vec[j][k][l] + weights[j][k][l];
                        }
                    }
                }
            }
            else if(world_rank == i)
            {
                MPI_Send(data, num_g_z*num_g_y*num_g_x, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD);
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
void collect_grid_i(int world_size,int world_rank,vector<vector<vector <int>>> &data_vec)
{
    int i = 0;                            //Standard variable used in loops
    int j = 0;                            //Standard variable used in loops
    int k = 0;                            //Standard variable used in loops
    int l = 0;                            //Standard variable used in loops
    MPI_Status status;                    //Used for mpi_recv
    int num_g_z = data_vec.size();
    int num_g_y = data_vec[0].size();
    int num_g_x = data_vec[0][0].size();
    int weights[num_g_z][num_g_y][num_g_x];        //An array to hold the grid being collected
    int data[num_g_z][num_g_y][num_g_x];

    //copy vector data to an array.
    for(i=0; i<num_g_z; i++) //loop over z
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            for(k=0; k<num_g_x; k++) //loop over x
            {
                data[i][j][k] = data_vec[i][j][k];
            }
        }
    }

    //collect array data
    if(world_size > 1)
    {
        for(i=1; i<world_size; i++)
        {
            if(world_rank == 0)
            {
                MPI_Recv(weights, num_g_z*num_g_y*num_g_x, MPI_INT, i, 13, MPI_COMM_WORLD, &status);

                for(j=0; j<num_g_z; j++) //loop over z
                {
                    for(k=0; k<num_g_y; k++) //loop over y
                    {
                        for(l=0; l<num_g_x; l++) //loop over x
                        {
                            data_vec[j][k][l] = data_vec[j][k][l] + weights[j][k][l];
                        }
                    }
                }
            }
            else if(world_rank == i)
            {
                MPI_Send(data, num_g_z*num_g_y*num_g_x, MPI_INT, 0, 13, MPI_COMM_WORLD);
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
void broadcast_grid_d(int world_size,int world_rank,vector<vector<vector <double>>> &data_vec)
{
    int i = 0;                            //Standard variable used in loops
    int j = 0;                            //Standard variable used in loops
    int k = 0;                            //Standard variable used in loops
    int l = 0;                            //Standard variable used in loops
    MPI_Status status;                    //Used for mpi_recv
    int num_g_z = data_vec.size();
    int num_g_y = data_vec[0].size();
    int num_g_x = data_vec[0][0].size();
    double weights[num_g_z][num_g_y][num_g_x];     //An array to hold the grid being collected
    double data[num_g_z][num_g_y][num_g_x];

    //copy vector data to an array.
    for(i=0; i<num_g_z; i++) //loop over z
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            for(k=0; k<num_g_x; k++) //loop over x
            {
                data[i][j][k] = data_vec[i][j][k];
            }
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
                MPI_Send(data, num_g_z*num_g_y*num_g_x, MPI_DOUBLE, i, 13, MPI_COMM_WORLD);
            }
            else if(world_rank == i)
            {
                MPI_Recv(weights, num_g_z*num_g_y*num_g_x, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD, &status);
                for(j=0; j<num_g_z; j++) //loop over z
                {
                    for(k=0; k<num_g_y; k++) //loop over y
                    {
                        for(l=0; l<num_g_x; l++) //loop over x
                        {
                            data_vec[j][k][l] = weights[j][k][l];
                        }
                    }
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function broadcasts grid data (int) to all the ranks.                                                //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void broadcast_grid_i(int world_size,int world_rank,vector<vector<vector <int>>> &data_vec)
{
    int i = 0;                            //Standard variable used in loops
    int j = 0;                            //Standard variable used in loops
    int k = 0;                            //Standard variable used in loops
    int l = 0;                            //Standard variable used in loops
    MPI_Status status;                    //Used for mpi_recv
    int num_g_z = data_vec.size();
    int num_g_y = data_vec[0].size();
    int num_g_x = data_vec[0][0].size();
    int weights[num_g_z][num_g_y][num_g_x];     //An array to hold the grid being collected
    int data[num_g_z][num_g_y][num_g_x];

    //copy vector data to an array.
    for(i=0; i<num_g_z; i++) //loop over z
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            for(k=0; k<num_g_x; k++) //loop over x
            {
                data[i][j][k] = data_vec[i][j][k];
            }
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
                MPI_Send(data, num_g_z*num_g_y*num_g_x, MPI_INT, i, 13, MPI_COMM_WORLD);
            }
            else if(world_rank == i)
            {
                MPI_Recv(weights, num_g_z*num_g_y*num_g_x, MPI_INT, 0, 13, MPI_COMM_WORLD, &status);
                for(j=0; j<num_g_z; j++) //loop over z
                {
                    for(k=0; k<num_g_y; k++) //loop over y
                    {
                        for(l=0; l<num_g_x; l++) //loop over x
                        {
                            data_vec[j][k][l] = weights[j][k][l];
                        }
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
void exclude_insignificant_data(double cutoff,double avg_rho,vector<vector<vector <double>>> &rho,vector<vector<vector <int>>> &nan)
{
    int i=0;
    int j=0;
    int k=0;

    int num_g_z = rho.size();
    int num_g_y = rho[0].size();
    int num_g_x = rho[0][0].size();

    for(i=0; i<num_g_z; i++) //loop over z-dimension
    {
        for(j=0; j<num_g_y; j++) //loop over y-dimension
        {
            for(k=0; k<num_g_x; k++) //loop over x-dimension
            {
                if(rho[i][j][k] > cutoff*avg_rho)
                {
                    nan[i][j][k] = 0;
                }
                else 
                {
                    nan[i][j][k] = 1;
                }  
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes the grid to an output file.                                                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_grid_to_file_d(double cell_size,string my_file_name,vector<vector<vector <int>>> &nan,vector<vector<vector <double>>> &data,double ex_val)
{
    int i=0;
    int j=0;
    int k=0;
    int num_g_z = data.size();
    int num_g_y = data[0].size();
    int num_g_x = data[0][0].size();

    FILE *my_file = fopen(my_file_name.c_str(), "w");
    if(my_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",my_file_name.c_str());
    }

    fprintf(my_file,"#data calclulated using the MosAT grid function \n");
    fprintf(my_file,"object 1 class gridpositions counts %d %d %d \n",num_g_x,num_g_y,num_g_z);
    fprintf(my_file,"origin 0.0 0.0 0.0 \n");
    fprintf(my_file,"delta %f 0 0 \n",10*cell_size);
    fprintf(my_file,"delta 0 %f 0 \n",10*cell_size);
    fprintf(my_file,"delta 0 0 %f \n",10*cell_size);
    fprintf(my_file,"object 2 class gridconnections counts %d %d %d \n",num_g_x,num_g_y,num_g_z);
    fprintf(my_file,"object 3 class array type double rank 0 items %d data follows \n",num_g_x*num_g_y*num_g_z);

    int count = 1;

    for(k=0; k<num_g_x; k++) //loop over x-dimension
    {
        for(j=0; j<num_g_y; j++) //loop over y-dimension
        {
            for(i=0; i<num_g_z; i++) //loop over z-dimension
            {
                if(nan[i][j][k] == 0)
                {
                    fprintf(my_file,"%e",data[i][j][k]);
                }
                else //data excluded
                {
                    fprintf(my_file,"%e",ex_val);
                }

                if(count%3 == 0 && count > 0)
                {
                    fprintf(my_file,"\n");
                }
                else 
                {
                    fprintf(my_file," ");
                }
                count++;
            }
        }
    }

    fclose(my_file);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes the grid to an output file.                                                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_grid_to_file_i(double cell_size,string my_file_name,vector<vector<vector <int>>> &nan,vector<vector<vector <int>>> &data,int ex_val)
{
    int i=0;
    int j=0;
    int k=0;
    int num_g_z = data.size();
    int num_g_y = data[0].size();
    int num_g_x = data[0][0].size();

    FILE *my_file = fopen(my_file_name.c_str(), "w");
    if(my_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",my_file_name.c_str());
    }

    fprintf(my_file,"#data calclulated using the MosAT grid function \n");
    fprintf(my_file,"object 1 class gridpositions counts %d %d %d \n",num_g_x,num_g_y,num_g_z);
    fprintf(my_file,"origin 0.0 0.0 0.0 \n");
    fprintf(my_file,"delta %f 0 0 \n",10*cell_size);
    fprintf(my_file,"delta 0 %f 0 \n",10*cell_size);
    fprintf(my_file,"delta 0 0 %f \n",10*cell_size);
    fprintf(my_file,"object 2 class gridconnections counts %d %d %d \n",num_g_x,num_g_y,num_g_z);
    fprintf(my_file,"object 3 class array type int rank 0 items %d data follows \n",num_g_x*num_g_y*num_g_z);

    int count = 1;

    for(k=0; k<num_g_x; k++) //loop over x-dimension
    {
        for(j=0; j<num_g_y; j++) //loop over y-dimension
        {
            for(i=0; i<num_g_z; i++) //loop over z-dimension
            {
                if(nan[i][j][k] == 0)
                {
                    fprintf(my_file,"%10d",data[i][j][k]);
                }
                else //data excluded
                {
                    fprintf(my_file,"%10d",ex_val);
                }

                if(count%3 == 0 && count > 0)
                {
                    fprintf(my_file,"\n");
                }
                else 
                {
                    fprintf(my_file," ");
                }
                count++;
            }
        }
    }

    fclose(my_file);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the average rho.                                                                    //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double get_average_rho(vector<vector<vector <double>>> &rho)
{
    int i = 0;
    int j = 0;
    int k = 0;
    double avg_rho = 0;

    int num_g_z = rho.size();
    int num_g_y = rho[0].size();
    int num_g_x = rho[0][0].size();

    for(i=0; i<num_g_z; i++) //loop over z-dimension
    {
        for(j=0; j<num_g_y; j++) //loop over y-dimension
        {
            for(k=0; k<num_g_x; k++) //loop over x-dimension
            {
                avg_rho = avg_rho + rho[i][j][k];
            }
        } 
    }
    avg_rho = avg_rho/(double)(num_g_x*num_g_y*num_g_z);

    return avg_rho;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function divides each gridpoint by rho.                                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void normalize_the_grid(vector<vector<vector <double>>> &rho,vector<vector<vector <double>>> &data)
{
    int i = 0;
    int j = 0;
    int k = 0;
    int num_g_z = data.size();
    int num_g_y = data[0].size();
    int num_g_x = data[0][0].size();

    for(i=0; i<num_g_z; i++) //loop over z-dimension
    {
        for(j=0; j<num_g_y; j++) //loop over y-dimension
        {
            for(k=0; k<num_g_x; k++) //loop over x-dimension
            {
                if(rho[i][j][k] > 0)
                {
                    data[i][j][k] = data[i][j][k]/rho[i][j][k];
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function determines the upper and lower bounds surrounding the mapping atom center.                  //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_bounds(double hx,double hy,double hz,double radius,double cell_size,int num_g_x,int num_g_y,int num_g_z,int *lower_x,int *lower_y,int *lower_z,int *upper_x,int *upper_y,int *upper_z)
{   
    //calculate the upper and lower bounds
    *lower_x = floor((hx - radius)/cell_size);
    *lower_y = floor((hy - radius)/cell_size);
    *lower_z = floor((hz - radius)/cell_size);
    *upper_x = ceil((hx + radius)/cell_size);
    *upper_y = ceil((hy + radius)/cell_size);
    *upper_z = ceil((hz + radius)/cell_size);

    //add a 1 point buffer around the paremeter
    *lower_x = *lower_x - 1;
    *lower_y = *lower_y - 1;
    *lower_z = *lower_z - 1;
    *upper_x = *upper_x + 1;
    *upper_y = *upper_y + 1; 
    *upper_z = *upper_z + 1;

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
   
    //lower_z
    if(*lower_z < 0)
    {
        *lower_z = 0;
    }
    else if(*lower_z > num_g_z)
    {
        *lower_z = num_g_z - 1;
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

    //upper_z
    if(*upper_z > num_g_z)
    {
        *upper_z = num_g_z;
    }
    else if(*upper_z < 0)
    {
        *upper_z = 0;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function deposits a double to the grid.                                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void add_to_grid_d(double hx,double hy,double hz,double radius,double cell_size,vector<vector<vector <double>>> &grid,double data)
{
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int lower_x = 0;
    int lower_y = 0;
    int lower_z = 0;
    int upper_x = 0;
    int upper_y = 0;
    int upper_z = 0;
    int num_g_z = grid.size();
    int num_g_y = grid[0].size();
    int num_g_x = grid[0][0].size();

    get_bounds(hx,hy,hz,radius,cell_size,num_g_x,num_g_y,num_g_z,&lower_x,&lower_y,&lower_z,&upper_x,&upper_y,&upper_z);

    for(i=lower_z; i<upper_z; i++) //loop over the grid z-axis
    {
        for(j=lower_y; j<upper_y; j++) //loop over the grid y-axis
        {
            for(k=lower_x; k<upper_x; k++) //loop over the grid x-axis
            {
                double dist_z = i*cell_size - hz;
                double dist_y = j*cell_size - hy;
                double dist_x = k*cell_size - hx;

                double dist = sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);

                if(dist <= radius)
                {
                    grid[i][j][k] = grid[i][j][k] + data;
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function deposits a double to the grid checking check first                                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void controlled_add_to_grid_d(double hx,double hy,double hz,double radius,double cell_size,vector<vector<vector <double>>> &grid,double data,vector<vector<vector <int>>> &check)
{
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int lower_x = 0;
    int lower_y = 0;
    int lower_z = 0;
    int upper_x = 0;
    int upper_y = 0;
    int upper_z = 0;
    int num_g_z = grid.size();
    int num_g_y = grid[0].size();
    int num_g_x = grid[0][0].size();

    get_bounds(hx,hy,hz,radius,cell_size,num_g_x,num_g_y,num_g_z,&lower_x,&lower_y,&lower_z,&upper_x,&upper_y,&upper_z);

    for(i=lower_z; i<upper_z; i++) //loop over the grid z-axis
    {
        for(j=lower_y; j<upper_y; j++) //loop over the grid y-axis
        {
            for(k=lower_x; k<upper_x; k++) //loop over the grid x-axis
            {
                if(check[i][j][k] == 0)
                {
                    double dist_z = i*cell_size - hz;
                    double dist_y = j*cell_size - hy;
                    double dist_x = k*cell_size - hx;

                    double dist = sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);

                    if(dist <= radius)
                    {
                        grid[i][j][k] = grid[i][j][k] + data;
                    }
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
void write_to_grid_d(double hx,double hy,double hz,double radius,double cell_size,vector<vector<vector <double>>> &grid,double data)
{
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int lower_x = 0;
    int lower_y = 0;
    int lower_z = 0;
    int upper_x = 0;
    int upper_y = 0;
    int upper_z = 0;
    int num_g_z = grid.size();
    int num_g_y = grid[0].size();
    int num_g_x = grid[0][0].size();

    get_bounds(hx,hy,hz,radius,cell_size,num_g_x,num_g_y,num_g_z,&lower_x,&lower_y,&lower_z,&upper_x,&upper_y,&upper_z);

    for(i=lower_z; i<upper_z; i++) //loop over the grid z-axis
    {
        for(j=lower_y; j<upper_y; j++) //loop over the grid y-axis
        {
            for(k=lower_x; k<upper_x; k++) //loop over the grid x-axis
            {
                double dist_z = i*cell_size - hz;
                double dist_y = j*cell_size - hy;
                double dist_x = k*cell_size - hx;

                double dist = sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);

                if(dist <= radius)
                {
                    grid[i][j][k] = data;
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function deposits an int to the grid.                                                                //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void add_to_grid_i(double hx,double hy,double hz,double radius,double cell_size,vector<vector<vector <int>>> &grid,int data)
{
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int lower_x = 0;
    int lower_y = 0;
    int lower_z = 0;
    int upper_x = 0;
    int upper_y = 0;
    int upper_z = 0;
    int num_g_z = grid.size();
    int num_g_y = grid[0].size();
    int num_g_x = grid[0][0].size();

    get_bounds(hx,hy,hz,radius,cell_size,num_g_x,num_g_y,num_g_z,&lower_x,&lower_y,&lower_z,&upper_x,&upper_y,&upper_z);

    for(i=lower_z; i<upper_z; i++) //loop over the grid z-axis
    {
        for(j=lower_y; j<upper_y; j++) //loop over the grid y-axis
        {
            for(k=lower_x; k<upper_x; k++) //loop over the grid x-axis
            {
                double dist_z = i*cell_size - hz;
                double dist_y = j*cell_size - hy;
                double dist_x = k*cell_size - hx;

                double dist = sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);

                if(dist <= radius)
                {
                    grid[i][j][k] = grid[i][j][k] + data;
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
void write_to_grid_i(double hx,double hy,double hz,double radius,double cell_size,vector<vector<vector <int>>> &grid,int data)
{
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int lower_x = 0;
    int lower_y = 0;
    int lower_z = 0;
    int upper_x = 0;
    int upper_y = 0;
    int upper_z = 0;
    int num_g_z = grid.size();
    int num_g_y = grid[0].size();
    int num_g_x = grid[0][0].size();

    get_bounds(hx,hy,hz,radius,cell_size,num_g_x,num_g_y,num_g_z,&lower_x,&lower_y,&lower_z,&upper_x,&upper_y,&upper_z);

    for(i=lower_z; i<upper_z; i++) //loop over the grid z-axis
    {
        for(j=lower_y; j<upper_y; j++) //loop over the grid y-axis
        {
            for(k=lower_x; k<upper_x; k++) //loop over the grid x-axis
            {
                double dist_z = i*cell_size - hz;
                double dist_y = j*cell_size - hy;
                double dist_x = k*cell_size - hx;

                double dist = sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);

                if(dist <= radius)
                {
                    grid[i][j][k] = data;
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints information about the grid.                                                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void print_grid_stats(double box_x,double box_y,double box_z,matrix ibox,int num_g_x,int num_g_y,int num_g_z,double cell_size,int world_rank)
{
    if(world_rank == 0)
    {
        if(box_x == 0 && box_y == 0) //user did not specify a box size
        {
            printf("Using box at frame zero to construct the grid as follows: \n");
            printf(" %10s %10s %10s %10s %10s %10s %10s \n","box_x","box_y","box_z","cell_size","num_g_x","num_g_y","num_g_z");
            printf("-%10s-%10s-%10s-%10s-%10s-%10s-%10s \n","----------","----------","----------","----------","----------","----------","----------");
            printf(" %10f %10f %10f %10f %10d %10d %10d \n",ibox[XX][XX],ibox[YY][YY],ibox[ZZ][ZZ],cell_size,num_g_x,num_g_y,num_g_z);
        }
        else
        {
            printf("Using user specified box to construct the grid as follows: \n");
            printf(" %10s %10s %10s %10s %10s %10s %10s \n","box_x","box_y","box_z","cell_size","num_g_x","num_g_y","num_g_z");
            printf("-%10s-%10s-%10s-%10s-%10s-%10s-%10s \n","----------","----------","----------","----------","----------","----------","----------");
            printf(" %10f %10f %10f %10f %10d %10d %10d \n",box_x,box_y,box_z,cell_size,num_g_x,num_g_y,num_g_z);
        }
        printf("\n");
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints information about the grid.                                                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void print_grid_stats_3d(double box_x,double box_y,double box_z,int b_box_i,int num_g_x,int num_g_y,int num_g_z,double cell_size,int world_rank)
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

        printf(" %10s %10s %10s %10s %10s %10s %10s \n","box_x","box_y","box_z","cell_size","num_g_x","num_g_y","num_g_z");
        printf("-%10s-%10s-%10s-%10s-%10s-%10s-%10s \n","----------","----------","----------","----------","----------","----------","----------");
        printf(" %10f %10f %10f %10f %10d %10d %10d \n",box_x,box_y,box_z,cell_size,num_g_x,num_g_y,num_g_z);

        printf("\n");
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function determines how many grid points in the x and y dimensions.                                  //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_grid_size(double box_x,double box_y,double box_z,matrix ibox,int *num_g_x,int *num_g_y,int *num_g_z,double cell_size)
{
    if(box_x == 0 && box_y == 0 && box_z == 0) //user did not specify a box size
    {
        *num_g_x = (int)ceil(ibox[XX][XX]/cell_size)+1;
        *num_g_y = (int)ceil(ibox[YY][YY]/cell_size)+1;
        *num_g_z = (int)ceil(ibox[ZZ][ZZ]/cell_size)+1;
    }
    else
    {
        *num_g_x = (int)ceil(box_x/cell_size)+1;
        *num_g_y = (int)ceil(box_y/cell_size)+1;
        *num_g_z = (int)ceil(box_z/cell_size)+1;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function determines how many grid points in the x and y dimensions.                                  //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_grid_size_obj(double my_box_x,double my_box_y,double my_box_z,matrix ibox,int *num_g_x,int *num_g_y,int *num_g_z,double cell_size,double *box_x,double *box_y,double *box_z,int *b_box_i)
{
    if(my_box_x == 0 && my_box_y == 0 && my_box_z == 0) //user did not specify a box size
    {
        *num_g_x = (int)ceil(ibox[XX][XX]/cell_size)+1;
        *num_g_y = (int)ceil(ibox[YY][YY]/cell_size)+1;
        *num_g_z = (int)ceil(ibox[ZZ][ZZ]/cell_size)+1;
        *box_x = ibox[YY][YY];
        *box_y = ibox[YY][YY];
        *box_z = ibox[ZZ][ZZ];
        *b_box_i = 1;
    }
    else
    {
        *num_g_x = (int)ceil(my_box_x/cell_size)+1;
        *num_g_y = (int)ceil(my_box_y/cell_size)+1;
        *num_g_z = (int)ceil(my_box_z/cell_size)+1;
        *box_x = my_box_x;
        *box_y = my_box_y;
        *box_z = my_box_z;
        *b_box_i = 0;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function sizes and initializes the grid.                                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void size_and_init_grid_d(int num_g_x,int num_g_y,int num_g_z,vector<vector<vector <double>>> &grid)
{
    int i = 0;           //standard variable used in loops
    int j = 0;           //standard variable used in loops
    int k = 0;           //standard variable used in loops

    //set the size of the vectors
    grid.resize(num_g_z);

    for(i=0; i<num_g_z; i++) //loop over z
    {
        grid[i].resize(num_g_y);
        
        for(j=0; j<num_g_y; j++) //loop over y
        {
            grid[i][j].resize(num_g_x);
        }
    }

    //initialize the grid
    for(i=0; i<grid.size(); i++) //loop over z
    {
        for(j=0; j<grid[i].size(); j++) //loop over y
        {
            for(k=0; k<grid[i][j].size(); k++) //loop over x
            {
                grid[i][j][k] = 0;
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function sizes and initializes the grid.                                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void size_and_init_grid_i(int num_g_x,int num_g_y,int num_g_z,vector<vector<vector <int>>> &grid)
{
    int i = 0;           //standard variable used in loops
    int j = 0;           //standard variable used in loops
    int k = 0;           //standard variable used in loops

    //set the size of the vectors
    grid.resize(num_g_z);

    for(i=0; i<num_g_z; i++) //loop over z
    {
        grid[i].resize(num_g_y);

        for(j=0; j<num_g_y; j++) //loop over y
        {
            grid[i][j].resize(num_g_x);
        }
    }

    //initialize the grid
    for(i=0; i<grid.size(); i++) //loop over z
    {
        for(j=0; j<grid[i].size(); j++) //loop over y
        {
            for(k=0; k<grid[i][j].size(); k++) //loop over x
            {
                grid[i][j][k] = 0;
            }
        }
    }
}








//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This is a class for working with the grid                                                                //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Grid_3d
{
    private:
        string grid_file_name;                                                               //name of the file to write grid to
        string rho_file_name;                                                                //name of the file to write rho to
        int    i           = 0;                                                              //standard variable used in loops
        int    j           = 0;                                                              //standard variable used in loops
        int    k           = 0;                                                              //standard variable used in loops
        int    world_size  = 0;                                                              //size of mpi world
        int    world_rank  = 0;                                                              //rank of mpi process
        double APS         = 0.005;                                                          //area of a grid square
        double cell_size   = sqrt(APS);                                                      //grid point spacing
        int    num_g_x     = 1;                                                              //number of grid points in x
        int    num_g_y     = 1;                                                              //number of grid points in y
        int    num_g_z     = 1;                                                               //number of grid points in z
        int    b_box_i     = 1;                                                              //is the initial box used?
        double size_grid_x = 0;                                                              //size of grid in x
        double size_grid_y = 0;                                                              //size of grid in y
        double size_grid_z = 0;                                                              //size of grid in z
        double box_x       = 0;                                                              //box x
        double box_y       = 0;                                                              //box y
        double box_z       = 0;                                                              //box z

    public:
        dv3d grid{};                                                                         //long term grid data       
        dv3d rho{};                                                                          //long term grid rho data
        dv3d frame_grid{};                                                                   //current frame grid 
        iv3d cell_count{};                                                                   //how many atoms added to current frame grid
        iv3d nan{};                                                                          //tells if the grid point should be excluded
        dv3d stdev{};                                                                        //standard deviation of grid data    

    public:
        void   get_dim(double box_x,double box_y,double box_z,matrix ibox,double my_APS);    //determine the grid dimensions
        void   print_dim();                                                                  //print the grid dimensions
        void   set_output(string my_grid_file_name);                                         //set the name and format of the output file
        void   clean_frame();                                                                //clean current frame grid
        void   clean_this_grid(dv3d this_grid);                                              //clears the provided grid 
        void   stamp(double hx,double hy,double hz,double radius,double data);               //add a value to current frame grid
        void   direct_stamp(double hx,double hy,double hz,double radius,double data);        //stamp directly to long term grid
        void   norm_frame();                                                                 //normalize the current frame
        void   add_frame();                                                                  //add the current frame to long term grid 
        void   add_frame_direct();                                                           //add the current frame to long term grid. add count to rho
        void   write_frame(int global_frame,double ex_val);                                  //write the current frame to out put file
        void   collect_grid();                                                               //collect grid/rho from mpi processes
        void   normalize();                                                                  //normalize the grid
        void   exclude_data(double cutoff,int report);                                       //exclude insignificant data
        void   write_grid(double ex_val);                                                    //write the gid to file
        void   write_rho(double ex_val);                                                     //write rho to file
        void   bcast_grid();                                                                 //broadcast the grid to all mpi processes
        void   get_stdev(int b_stdev,int b_clean,Trajectory traj,double ex_val);             //compute the standard deviation
        void   match_grid(Grid_3d model_grid);                                               //set the grid dimensions to match another grid
        void   copy_rho(Grid_3d model_rho);                                                  //copy rho from another grid
        int    num_x();                                                                      //return num_g_x
        int    num_y();                                                                      //return num_g_y
        int    num_z();                                                                      //return num_g_z
        double get_cell_size();                                                              //return cell_size
        double get_aps();                                                                    //return APS
        double get_box_x();                                                                  //return the box x dimension
        double get_box_y();                                                                  //return the box y dimension
        double get_box_z();                                                                  //return the box z dimension
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return determines the grid dimensions                                                      //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d::get_dim(double my_box_x,double my_box_y,double my_box_z,matrix ibox,double my_APS)
{
    //get the world size and assign ranks
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);      //get the world size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);      //get the process rank

    APS       = my_APS;
    cell_size = sqrt(APS);

    get_grid_size_obj(my_box_x,my_box_y,my_box_z,ibox,&num_g_x,&num_g_y,&num_g_z,cell_size,&box_x,&box_y,&box_z,&b_box_i);

    size_grid_x = num_g_x*cell_size;
    size_grid_y = num_g_y*cell_size;
    size_grid_z = num_g_z*cell_size;

    size_and_init_grid_d(num_g_x,num_g_y,num_g_z,grid);
    size_and_init_grid_d(num_g_x,num_g_y,num_g_z,rho);
    size_and_init_grid_d(num_g_x,num_g_y,num_g_z,frame_grid);
    size_and_init_grid_i(num_g_x,num_g_y,num_g_z,cell_count);
    size_and_init_grid_i(num_g_x,num_g_y,num_g_z,nan);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function prints the grid dimensions                                                                 //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d::print_dim()
{
    print_grid_stats_3d(box_x,box_y,box_z,b_box_i,num_g_x,num_g_y,num_g_z,cell_size,world_rank);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function sets the output name and format                                                            //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d::set_output(string my_grid_file_name)
{
    grid_file_name = my_grid_file_name;

    rho_file_name = add_tag(grid_file_name,"_rho");
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return num_g_x                                                                             //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Grid_3d::num_x()
{
    return num_g_x;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return num_g_y                                                                             //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Grid_3d::num_y()
{
    return num_g_y;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return num_g_z                                                                             //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Grid_3d::num_z()
{
    return num_g_z;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns cell_size                                                                          //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Grid_3d::get_cell_size()
{
    return cell_size;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return APS                                                                                 //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Grid_3d::get_aps()
{
    return APS;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function deposits a quantity to the current frame grid                                              //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d::stamp(double hx,double hy,double hz,double radius,double data)
{
   add_to_grid_d(hx,hy,hz,radius,cell_size,frame_grid,data);
   add_to_grid_i(hx,hy,hz,radius,cell_size,cell_count,1);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function deposits a quantity to the long term grid                                                  //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d::direct_stamp(double hx,double hy,double hz,double radius,double data)
{
   add_to_grid_d(hx,hy,hz,radius,cell_size,grid,data);
   add_to_grid_d(hx,hy,hz,radius,cell_size,rho,1.0);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function clears the current frame grid                                                              //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d::clean_frame()
{
    for(i=0; i<num_g_z; i++) //loop over z
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            for(k=0; k<num_g_x; k++) //loop over x
            {    
                frame_grid[i][j][k] = 0;
                cell_count[i][j][k] = 0;
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function clears the entire grid                                                                     //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d::clean_this_grid(dv3d this_grid)
{
    for(i=0; i<num_g_z; i++) //loop over z
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            for(k=0; k<num_g_x; k++) //loop over x
            {
                this_grid[i][j][k] = 0.0;
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function normalizes the current frame                                                               //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d::norm_frame()
{
    for(i=0; i<num_g_z; i++) //loop over z
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            for(k=0; k<num_g_x; k++) //loop over x
            {
                if(cell_count[i][j][k] > 0)
                {
                    frame_grid[i][j][k] = frame_grid[i][j][k]/(double)cell_count[i][j][k];
                }
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function adds the current frame grid to the long long term sum                                      //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d::add_frame()
{
    for(i=0; i<num_g_z; i++) //loop over z
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            for(k=0; k<num_g_x; k++) //loop over x
            {
                if(cell_count[i][j][k] > 0)
                {
                    grid[i][j][k] = grid[i][j][k] + frame_grid[i][j][k];
                    rho[i][j][k]  = rho[i][j][k]  + 1.0;
                }
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function adds the current frame grid to the long long term sum. add count to rho                    //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d::add_frame_direct()
{
    for(i=0; i<num_g_z; i++) //loop over z
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            for(k=0; k<num_g_x; k++) //loop over x
            {
                if(cell_count[i][j][k] > 0)
                {
                    grid[i][j][k] = grid[i][j][k] + frame_grid[i][j][k];
                    rho[i][j][k]  = rho[i][j][k]  + cell_count[i][j][k];
                }
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function writes the current frame to file                                                           //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d::write_frame(int global_frame,double ex_val)
{
    int i = 0;
    int j = 0;
    int k = 0;

    string frame_file_name = add_tag(grid_file_name,to_string(global_frame));

    iv3d this_nan(cell_count.size(),iv2d(cell_count[0].size(),iv1d(cell_count[0][0].size(),0)));
    for(i=0; i<cell_count.size(); i++)
    {
        for(j=0; j<cell_count[0].size(); j++)
        {
            for(k=0; k<cell_count[0][0].size(); k++)
            {
                if(cell_count[i][j][k] == 0) //exclude data
                {
                    this_nan[i][j][k] = 1;
                }
            }
        }
    }

    write_grid_to_file_d(cell_size,frame_file_name,this_nan,frame_grid,ex_val);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function collects the grid from mpi processes                                                       //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d::collect_grid()
{
    collect_grid_d(world_size,world_rank,grid);
    collect_grid_d(world_size,world_rank,rho);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function normalizes the grid                                                                        //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d::normalize()
{
    normalize_the_grid(rho,grid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function excludes insignificant data                                                                //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d::exclude_data(double cutoff,int report)
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
void Grid_3d::write_grid(double ex_val)
{
    write_grid_to_file_d(cell_size,grid_file_name,nan,grid,ex_val);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function writes rho to output file                                                                  //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d::write_rho(double ex_val)
{
    write_grid_to_file_d(cell_size,rho_file_name,nan,rho,ex_val);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function broadcasts the grid                                                                        //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d::bcast_grid()
{
    broadcast_grid_d(world_size,world_rank,grid);
    broadcast_grid_d(world_size,world_rank,rho);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function computes the standard deviation                                                            //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d::get_stdev(int b_stdev,int b_clean,Trajectory traj,double ex_val)
{
    if(b_stdev == 1)
    {
        int i = 0;                               //standard variable used in loops
        int j = 0;                               //standard variable used in loops
        int k = 0;                               //standard variable used in loops
        char my_string[200];                     //Used for reading data from single frame files
        char line[200];                          //holds the current line

        size_and_init_grid_d(num_g_x,num_g_y,num_g_z,stdev);

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
                //first 8 lines are headers
                fgets(line, sizeof line, single_file);
                fgets(line, sizeof line, single_file);
                fgets(line, sizeof line, single_file);
                fgets(line, sizeof line, single_file);
                fgets(line, sizeof line, single_file);
                fgets(line, sizeof line, single_file);
                fgets(line, sizeof line, single_file);
                fgets(line, sizeof line, single_file);

                for(i=0; i<num_g_x; i++) //loop over x-dimension
                {
                    for(j=0; j<num_g_y; j++) //loop over y-dimension
                    {
                        for(k=0; k<num_g_z; k++) //loop over z-dimension
                        {
                            int result    = fscanf(single_file, "%s,", my_string);
                            double value  = atof(my_string);

                            if(grid[k][j][i] != ex_val)
                            {
                                double dif     = value - grid[k][j][i];
                                double square  = dif*dif;
                                stdev[k][j][i] = stdev[k][j][i] + square;
                            }
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
            for(i=0; i<num_g_z; i++) //loop over z-dimension
            {
                for(j=0; j<num_g_y; j++) //loop over y-dimension
                {
                    for(k=0; k<num_g_x; k++) //loop over x-dimension
                    {
                        if(rho[i][j][k] > 0)
                        {
                            double norm        = stdev[i][j][k]/(rho[i][j][k]-1.0);
                            double stdeviation = sqrt(norm);

                            stdev[i][j][k] = stdeviation;
                        }
                    }
                }
            }

            //set the stdev file name
            string stdev_file_name = add_tag(grid_file_name,"_stdev");

            //now print the standard deviation to file
            write_grid_to_file_d(cell_size,stdev_file_name,nan,stdev,ex_val);

            for(i=0; i<num_g_z; i++) //loop over z-dimension
            {
                for(j=0; j<num_g_y; j++) //loop over y-dimension
                {
                    for(k=0; k<num_g_x; k++) //loop over x-dimension
                    {
                        if(rho[i][j][k] > 0)
                        {
                            stdev[i][j][k] = stdev[i][j][k]/sqrt(rho[i][j][k]);
                        }
                    }
                }
            }

            //set the stdev file name
            string stem_file_name = add_tag(grid_file_name,"_stem");

            //now print the standard deviation to file
            write_grid_to_file_d(cell_size,stem_file_name,nan,stdev,ex_val);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns the box x dimension                                                                //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Grid_3d::get_box_x()
{
    return box_x;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns the box y dimension                                                                //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Grid_3d::get_box_y()
{
    return box_y;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns the box z dimension                                                                //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Grid_3d::get_box_z()
{
    return box_z;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function copies the dimesion of a grid to another grid                                              //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d::match_grid(Grid_3d model_grid)
{
    //get the world size and assign ranks
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);      //get the world size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);      //get the process rank

    num_g_x   = model_grid.num_x();
    num_g_y   = model_grid.num_y();
    num_g_z   = model_grid.num_z();

    APS       = model_grid.get_aps();
    cell_size = model_grid.get_cell_size();

    box_x = get_box_x();
    box_y = get_box_y();
    box_z = get_box_z();

    size_grid_x = num_g_x*cell_size;
    size_grid_y = num_g_y*cell_size;
    size_grid_z = num_g_z*cell_size;

    size_and_init_grid_d(num_g_x,num_g_y,num_g_z,grid);
    size_and_init_grid_d(num_g_x,num_g_y,num_g_z,rho);
    size_and_init_grid_d(num_g_x,num_g_y,num_g_z,frame_grid);
    size_and_init_grid_i(num_g_x,num_g_y,num_g_z,cell_count);
    size_and_init_grid_i(num_g_x,num_g_y,num_g_z,nan);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function copies rho from another grid                                                               //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d::copy_rho(Grid_3d model_rho)
{
    for(i=0; i<num_g_z; i++) //loop over z
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            for(k=0; k<num_g_x; k++) //loop over x
            {
                rho[i][j][k] = model_rho.rho[i][j][k];
            }
        }
    }
}











//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This is a class for working with a single grid of ints                                                   //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Grid_3d_i
{
    private:
        string grid_file_name;                                                                //name of the file to write grid to
        int    i           = 0;                                                               //standard variable used in loops
        int    j           = 0;                                                               //standard variable used in loops
        int    k           = 0;                                                               //standard variable used in loops
        int    world_size  = 0;                                                               //size of mpi world
        int    world_rank  = 0;                                                               //rank of mpi process
        double APS         = 0.005;                                                           //area of a grid square
        double cell_size   = sqrt(APS);                                                       //grid point spacing
        int    num_g_x     = 1;                                                               //number of grid points in x
        int    num_g_y     = 1;                                                               //number of grid points in y
        int    num_g_z     = 1;                                                               //number of grid points in z
        int    b_box_i     = 1;                                                               //is the initial box used?
        double size_grid_x = 0;                                                               //size of grid in x
        double size_grid_y = 0;                                                               //size of grid in y
        double size_grid_z = 0;                                                               //size of grid in z
        double box_x       = 0;                                                               //box x
        double box_y       = 0;                                                               //box y
        double box_z       = 0;                                                               //box y

    public:
        iv3d grid{};                                                                          //long term grid data

    public:
        void   get_dim(double my_box_x,double my_box_y,double my_box_z,matrix ibox,double my_APS);//gets and sets the grid dimensions
        void   set_dim(double my_APS,double my_num_g_x,double my_num_g_y,double my_num_g_z);  //sets the grid dimensions
        void   set_output(string my_grid_file_name);                                          //set the name and format of the output file
        void   print_dim();                                                                   //print dimensions of the grid
        void   clean_grid();                                                                  //clean the grid
        void   stamp(double hx,double hy,double hz,double radius,int data);                   //add a value to the grid
        void   paint(double hx,double hy,double hz,double radius,int data);                   //set the value for the grid points
        void   collect_grid();                                                                //collect grid/rho from mpi processes
        void   write_grid(Grid_3d_i nan,int ex_val);                                          //write the gid to file
        void   bcast_grid();                                                                  //broadcast the grid to all mpi processes
        void   copy_grid(Grid_3d_i model_grid);                                               //copy from another grid
        int    num_x();                                                                       //return num_g_x
        int    num_y();                                                                       //return num_g_y
        int    num_z();                                                                       //return num_g_z
        double get_cell_size();                                                               //return cell_size
        double get_aps();                                                                     //return APS
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return determines the grid dimensions                                                      //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d_i::get_dim(double my_box_x,double my_box_y,double my_box_z,matrix ibox,double my_APS)
{
    //get the world size and assign ranks
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);      //get the world size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);      //get the process rank

    APS       = my_APS;
    cell_size = sqrt(APS);

    get_grid_size_obj(my_box_x,my_box_y,my_box_z,ibox,&num_g_x,&num_g_y,&num_g_z,cell_size,&box_x,&box_y,&box_z,&b_box_i);

    size_grid_x = num_g_x*cell_size;
    size_grid_y = num_g_y*cell_size;
    size_grid_z = num_g_z*cell_size;

    size_and_init_grid_i(num_g_x,num_g_y,num_g_z,grid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return determines the grid dimensions                                                      //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d_i::set_dim(double my_APS,double my_num_g_x,double my_num_g_y,double my_num_g_z)
{
    //get the world size and assign ranks
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);      //get the world size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);      //get the process rank

    APS       = my_APS;
    cell_size = sqrt(APS);
    num_g_x   = my_num_g_x;
    num_g_y   = my_num_g_y;
    num_g_z   = my_num_g_z;

    size_grid_x = num_g_x*cell_size;
    size_grid_y = num_g_y*cell_size;

    size_and_init_grid_i(num_g_x,num_g_y,num_g_z,grid);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function sets the output name and format                                                            //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d_i::set_output(string my_grid_file_name)
{
    grid_file_name = my_grid_file_name;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function prints the grid dimensions                                                                 //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d_i::print_dim()
{
    print_grid_stats_3d(box_x,box_y,box_z,b_box_i,num_g_x,num_g_y,num_g_z,cell_size,world_rank);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function clears the grid                                                                            //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d_i::clean_grid()
{
    for(i=0; i<num_g_z; i++) //loop over z
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            for(k=0; k<num_g_x; k++) //loop over x
            {
                grid[i][j][k] = 0;
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function deposits a quantity to grid                                                                //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d_i::stamp(double hx,double hy,double hz,double radius,int data)
{
   add_to_grid_i(hx,hy,hz,radius,cell_size,grid,data);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function writes a quantity to grid                                                                  //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d_i::paint(double hx,double hy,double hz,double radius,int data)
{
   write_to_grid_i(hx,hy,hz,radius,cell_size,grid,data);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function collects the grid from mpi processes                                                       //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d_i::collect_grid()
{
    collect_grid_i(world_size,world_rank,grid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function writes the grid to output file                                                             //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d_i::write_grid(Grid_3d_i nan,int ex_val)
{
    write_grid_to_file_i(cell_size,grid_file_name,nan.grid,grid,ex_val);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function broadcasts the grid                                                                        //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d_i::bcast_grid()
{
    broadcast_grid_i(world_size,world_rank,grid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function copies data from another grid                                                              //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d_i::copy_grid(Grid_3d_i model_grid)
{
    for(i=0; i<num_g_z; i++) //loop over z
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            for(k=0; k<num_g_x; k++) //loop over x
            {
                grid[i][j][k] = model_grid.grid[i][j][k];
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return num_g_x                                                                             //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Grid_3d_i::num_x()
{
    return num_g_x;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return num_g_y                                                                             //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Grid_3d_i::num_y()
{
    return num_g_y;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return num_g_z                                                                             //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Grid_3d_i::num_z()
{
    return num_g_z;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns cell_size                                                                          //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Grid_3d_i::get_cell_size()
{
    return cell_size;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return APS                                                                                 //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Grid_3d_i::get_aps()
{
    return APS;
}








//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This is a class for working with a single grid of doubles                                                //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Grid_3d_d
{
    private:
        string grid_file_name;                                                                //name of the file to write grid to
        int    i           = 0;                                                               //standard variable used in loops
        int    j           = 0;                                                               //standard variable used in loops
        int    k           = 0;                                                               //standard variable used in loops
        int    world_size  = 0;                                                               //size of mpi world
        int    world_rank  = 0;                                                               //rank of mpi process
        double APS         = 0.005;                                                           //area of a grid square
        double cell_size   = sqrt(APS);                                                       //grid point spacing
        int    num_g_x     = 1;                                                               //number of grid points in x
        int    num_g_y     = 1;                                                               //number of grid points in y
        int    num_g_z     = 1;                                                               //number of grid points in z
        int    b_box_i     = 1;                                                               //is the initial box used?
        double size_grid_x = 0;                                                               //size of grid in x
        double size_grid_y = 0;                                                               //size of grid in y
        double size_grid_z = 0;                                                               //size of grid in z
        double box_x       = 0;                                                               //box x
        double box_y       = 0;                                                               //box y
        double box_z       = 0;                                                               //box y

    public:
        dv3d grid{};                                                                          //long term grid data

    public:
        void   get_dim(double my_box_x,double my_box_y,double my_box_z,matrix ibox,double my_APS);         //gets and set the grid dimensions
        void   set_dim(double my_APS,double my_num_g_x,double my_num_g_y,double my_num_g_z);               //sets the grid dimensions
        void   set_output(string my_grid_file_name);                                                       //set the name and format of the output file
        void   print_dim();                                                                                //print dimensions of the grid
        void   clean_grid();                                                                               //clean the grid
        void   stamp(double hx,double hy,double hz,double radius,double data);                             //add a value to the grid
        void   controlled_stamp(double hx,double hy,double hz,double radius,double data,Grid_3d_i &check); //add a value to the grid if check equals zero. 
        void   collect_grid();                                                                             //collect grid/rho from mpi processes
        void   paint(double hx,double hy,double hz,double radius,double data);                             //set the value for the grid points
        void   normalize(Grid_3d_d rho);                                                                   //normalize the grid
        void   normalize_constant(double c);                                                               //normalize the grid using a constant normaizing factor
        void   write_grid(Grid_3d_i nan,double ex_val);                                                    //write the gid to file
        void   bcast_grid();                                                                               //broadcast the grid to all mpi processes
        void   copy_grid(Grid_3d_d model_grid);                                                            //copy from another grid
        int    num_x();                                                                                    //return num_g_x
        int    num_y();                                                                                    //return num_g_y
        int    num_z();                                                                                    //return num_g_z
        double get_cell_size();                                                                            //return cell_size
        double get_aps();                                                                                  //return APS
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return determines the grid dimensions                                                      //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d_d::get_dim(double my_box_x,double my_box_y,double my_box_z,matrix ibox,double my_APS)
{
    //get the world size and assign ranks
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);      //get the world size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);      //get the process rank

    APS       = my_APS;
    cell_size = sqrt(APS);

    get_grid_size_obj(my_box_x,my_box_y,my_box_z,ibox,&num_g_x,&num_g_y,&num_g_z,cell_size,&box_x,&box_y,&box_z,&b_box_i);

    size_grid_x = num_g_x*cell_size;
    size_grid_y = num_g_y*cell_size;
    size_grid_z = num_g_z*cell_size;

    size_and_init_grid_d(num_g_x,num_g_y,num_g_z,grid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return determines the grid dimensions                                                      //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d_d::set_dim(double my_APS,double my_num_g_x,double my_num_g_y,double my_num_g_z)
{
    //get the world size and assign ranks
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);      //get the world size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);      //get the process rank

    APS       = my_APS;
    cell_size = sqrt(APS);
    num_g_x   = my_num_g_x;
    num_g_y   = my_num_g_y;
    num_g_z   = my_num_g_z; 

    size_grid_x = num_g_x*cell_size;
    size_grid_y = num_g_y*cell_size;
    size_grid_z = num_g_z*cell_size; 

    size_and_init_grid_d(num_g_x,num_g_y,num_g_z,grid);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function sets the output name and format                                                            //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d_d::set_output(string my_grid_file_name)
{
    grid_file_name = my_grid_file_name;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function prints the grid dimensions                                                                 //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d_d::print_dim()
{
    print_grid_stats_3d(box_x,box_y,box_z,b_box_i,num_g_x,num_g_y,num_g_z,cell_size,world_rank);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function clears the grid                                                                            //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d_d::clean_grid()
{
    for(i=0; i<num_g_z; i++) //loop over z
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            for(k=0; k<num_g_x; k++) //loop over x
            {
                grid[i][j][k] = 0.0;
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function deposits a quantity to grid                                                                //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d_d::stamp(double hx,double hy,double hz,double radius,double data)
{
   add_to_grid_d(hx,hy,hz,radius,cell_size,grid,data);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function deposits a quantity to grid checking first check                                           //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d_d::controlled_stamp(double hx,double hy,double hz,double radius,double data,Grid_3d_i &check)
{
   controlled_add_to_grid_d(hx,hy,hz,radius,cell_size,grid,data,check.grid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function writes a quantity to grid                                                                  //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d_d::paint(double hx,double hy,double hz,double radius,double data)
{
   write_to_grid_d(hx,hy,hz,radius,cell_size,grid,data);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function collects the grid from mpi processes                                                       //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d_d::collect_grid()
{
    collect_grid_d(world_size,world_rank,grid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function normalizes the grid                                                                        //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d_d::normalize(Grid_3d_d rho)
{
    for(i=0; i<num_g_z; i++) //loop over z
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            for(k=0; k<num_g_x; k++) //loop over x
            {
                if(rho.grid[i][j][k] > 0.0)
                {
                    grid[i][j][k] = grid[i][j][k]/rho.grid[i][j][k];
                }
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function normalizes the grid                                                                        //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d_d::normalize_constant(double c)
{
    for(i=0; i<num_g_z; i++) //loop over z
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            for(k=0; k<num_g_x; k++) //loop over x
            {
                if(c != 0.0) //avoid division by zeor
                {
                    grid[i][j][k] = grid[i][j][k]/c;
                }
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function writes the grid to output file                                                             //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d_d::write_grid(Grid_3d_i nan,double ex_val)
{
    write_grid_to_file_d(cell_size,grid_file_name,nan.grid,grid,ex_val);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function broadcasts the grid                                                                        //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d_d::bcast_grid()
{
    broadcast_grid_d(world_size,world_rank,grid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function copies rho from another grid                                                               //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_3d_d::copy_grid(Grid_3d_d model_grid)
{
    for(i=0; i<num_g_z; i++) //loop over z
    {
        for(j=0; j<num_g_y; j++) //loop over y
        {
            for(k=0; k<num_g_x; k++) //loop over x
            {
                grid[i][j][k] = model_grid.grid[i][j][k];
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return num_g_x                                                                             //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Grid_3d_d::num_x()
{
    return num_g_x;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return num_g_y                                                                             //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Grid_3d_d::num_y()
{
    return num_g_y;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return num_g_z                                                                             //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Grid_3d_d::num_z()
{
    return num_g_z;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns cell_size                                                                          //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Grid_3d_d::get_cell_size()
{
    return cell_size;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return APS                                                                                 //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Grid_3d_d::get_aps()
{
    return APS;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function excludes insignificant data                                                                //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exclude_data(Grid_3d_d rho,Grid_3d_i &nan,double cutoff,int report)
{
    double avg_rho = get_average_rho(rho.grid);

    if(report == 1)
    {
        printf("The average rho is %f. Excluding lattice points with fewer than %f samples. \n",avg_rho,cutoff*avg_rho);
    }
    exclude_insignificant_data(cutoff,avg_rho,rho.grid,nan.grid);
}


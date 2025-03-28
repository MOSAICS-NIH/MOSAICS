
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function will collect a 1d array of doubles from the nodes and reconstruct them into a single array  //
// such that the time is continuous for the trajectory                                                       //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void collect_1d_double_array(int world_size,vector <int> world_frames,int world_rank,double world_ary[],double ary[])
{
   int i = 0;
   int j = 0;
   int frame_count = -1;
   MPI_Status status;

   for(i=0; i<world_size; i++) //loop over all nodes
    {
        for(j=0; j<world_frames[i]; j++) //loop ovrer each nodes frames
        {
            frame_count = frame_count + 1;
            if(i == 0 && world_rank == 0) //rank zero simply copies its data
            {
                world_ary[frame_count] = ary[j];
            }
            else // all non zero rank nodes must send data to rank 0 node
            {
                if(world_rank == i)
                {
                    double my_ary = ary[j];
                    MPI_Send(&my_ary, 1, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD);
                }
                else if(world_rank == 0)
                {
                    double my_ary = 0;
                    MPI_Recv(&my_ary, 1, MPI_DOUBLE, i, 13, MPI_COMM_WORLD, &status);
                    world_ary[frame_count] = my_ary;
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects a 1-d vector of doubles                                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void collect_dv1d(int world_size,int world_rank,dv1d &my_vec)
{
    int i = 0;
    int j = 0;

    if(world_size > 0)
    {
        for(i=1; i<world_size; i++)
        {
            if(world_rank == 0)
            {
                int size;                     //how many items to be received
                MPI_Recv(&size, 1, MPI_INT, i, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                double recv[size];            //array to hold received items
                MPI_Recv(recv, size, MPI_DOUBLE, i, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                //add the items to the vector
                for(j=0; j<size; j++)
                {
                    my_vec.push_back(recv[j]);
                }
            }
            else if(world_rank == i)
            {
                int size = my_vec.size();  //how many items to send

                MPI_Send(&size, 1, MPI_INT, 0, 13, MPI_COMM_WORLD);

                double snd[size];
                for(j=0; j<size; j++)
                {
                    snd[j] = my_vec[j];
                }
                MPI_Send(snd, size, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects a 1-d vector of ints                                                               //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void collect_iv1d(int world_size,int world_rank,iv1d &my_vec)
{
    int i = 0;
    int j = 0;

    if(world_size > 0)
    {
        for(i=1; i<world_size; i++)
        {
            if(world_rank == 0)
            {
                int size;                     //how many items to be received
                MPI_Recv(&size, 1, MPI_INT, i, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                int recv[size];            //array to hold received items
                MPI_Recv(recv, size, MPI_INT, i, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                //add the items to the vector
                for(j=0; j<size; j++)
                {
                    my_vec.push_back(recv[j]);
                }
            }
            else if(world_rank == i)
            {
                int size = my_vec.size();  //how many items to send

                MPI_Send(&size, 1, MPI_INT, 0, 13, MPI_COMM_WORLD);

                int snd[size];
                for(j=0; j<size; j++)
                {
                    snd[j] = my_vec[j];
                }
                MPI_Send(snd, size, MPI_INT, 0, 13, MPI_COMM_WORLD);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects a 1-d vector of doubles. returns new ve                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
dv1d collect_and_clone_dv1d(int world_size,int world_rank,dv1d &my_vec)
{
    int i = 0;
    int j = 0;

    dv1d temp_vec(0,0.0);

    if(world_rank == 0)
    {
        for(i=0; i<my_vec.size(); i++)
        {
            temp_vec.push_back(my_vec[i]);
        }
    }

    if(world_size > 0)
    {
        for(i=1; i<world_size; i++)
        {
            if(world_rank == 0)
            {
                int size;                     //how many items to be received
                MPI_Recv(&size, 1, MPI_INT, i, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                double recv[size];            //array to hold received items
                MPI_Recv(recv, size, MPI_DOUBLE, i, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                //add the items to the vector
                for(j=0; j<size; j++)
                {
                    temp_vec.push_back(recv[j]);
                }
            }
            else if(world_rank == i)
            {
                int size = my_vec.size();  //how many items to send

                MPI_Send(&size, 1, MPI_INT, 0, 13, MPI_COMM_WORLD);

                double snd[size];
                for(j=0; j<size; j++)
                {
                    snd[j] = my_vec[j];
                }
                MPI_Send(snd, size, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD);
            }
        }
    }
    return temp_vec;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects a 1-d vector of ints. returns a new vec                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
iv1d collect_and_clone_iv1d(int world_size,int world_rank,iv1d &my_vec)
{
    int i = 0;
    int j = 0;

    iv1d temp_vec(0,0);

    if(world_rank == 0)
    {
        for(i=0; i<my_vec.size(); i++)
        {
            temp_vec.push_back(my_vec[i]);
        }
    }

    if(world_size > 0)
    {
        for(i=1; i<world_size; i++)
        {
            if(world_rank == 0)
            {
                int size;                     //how many items to be received
                MPI_Recv(&size, 1, MPI_INT, i, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                int recv[size];            //array to hold received items
                MPI_Recv(recv, size, MPI_INT, i, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                //add the items to the vector
                for(j=0; j<size; j++)
                {
                    temp_vec.push_back(recv[j]);
                }
            }
            else if(world_rank == i)
            {
                int size = my_vec.size();  //how many items to send

                MPI_Send(&size, 1, MPI_INT, 0, 13, MPI_COMM_WORLD);

                int snd[size];
                for(j=0; j<size; j++)
                {
                    snd[j] = my_vec[j];
                }
                MPI_Send(snd, size, MPI_INT, 0, 13, MPI_COMM_WORLD);
            }
        }
    }
    return temp_vec;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects a 2-d vector of doubles                                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void collect_dv2d(int world_size,int world_rank,dv2d &my_vec)
{
    int i = 0;
    int j = 0;
    int k = 0;

    if(world_size > 0)
    {
        for(i=1; i<world_size; i++)
        {
            if(world_rank == 0)
            {
                int num_dv1d;                    //how many dv1d to be received
                MPI_Recv(&num_dv1d, 1, MPI_INT, i, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                for(j=0; j<num_dv1d; j++) //loop over y
                {
                    int size = 0;
 
                    MPI_Recv(&size, 1, MPI_INT, i, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                    double recv[size];            //array to hold received items
                    MPI_Recv(recv, size, MPI_DOUBLE, i, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                    //add the items to the vector
                    dv1d this_vec(0,0.0);
                    for(k=0; k<size; k++)
                    {
                        this_vec.push_back(recv[k]);
                    }

                    //add the vector to my_vec
                    my_vec.push_back(this_vec);
                }
            }
            else if(world_rank == i)
            {
                int num_dv1d = my_vec.size();  //how many items to send

                MPI_Send(&num_dv1d, 1, MPI_INT, 0, 13, MPI_COMM_WORLD);

                for(j=0; j<num_dv1d; j++) //loop over the 1d vectors
                {
                    int size = my_vec[j].size();

                    MPI_Send(&size, 1, MPI_INT, 0, 13, MPI_COMM_WORLD);

                    double snd[size];
                    for(k=0; k<size; k++)
                    {
                        snd[k] = my_vec[j][k];
                    }
                    MPI_Send(snd, size, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD);
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects a 1-d vector of doubles and summs the contents for each element                    //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void collect_and_sum_dv1d(int world_size,int world_rank,dv1d &my_vec)
{
    int i = 0;
    int j = 0;

    if(world_size > 0)
    {
        for(i=1; i<world_size; i++)
        {
            if(world_rank == 0)
            {
                int size;                     //how many items to be received
                MPI_Recv(&size, 1, MPI_INT, i, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                double recv[size];            //array to hold received items
                MPI_Recv(recv, size, MPI_DOUBLE, i, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                //add the items to the vector
                for(j=0; j<size; j++)
                {
                    my_vec[j] = my_vec[j] + recv[j];
                }
            }
            else if(world_rank == i)
            {
                int size = my_vec.size();  //how many items to send

                MPI_Send(&size, 1, MPI_INT, 0, 13, MPI_COMM_WORLD);

                double snd[size];
                for(j=0; j<size; j++)
                {
                    snd[j] = my_vec[j];
                }
                MPI_Send(snd, size, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function broadcasts a 1-d vector of doubles                                                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void broadcast_dv1d(int world_size,int world_rank,dv1d &my_vec)
{
    int i = 0;
    int j = 0;

    if(world_size > 0)
    {
        for(i=1; i<world_size; i++)
        {
            if(world_rank == 0)
            {
                int size = my_vec.size();  //how many items to send
          
                MPI_Send(&size, 1, MPI_INT, i, 13, MPI_COMM_WORLD);

                double snd[size];
                for(j=0; j<size; j++)
                {
                    snd[j] = my_vec[j];
                }
                MPI_Send(snd, size, MPI_DOUBLE, i, 13, MPI_COMM_WORLD);
            }
            else if(world_rank == i)
            {
                int size;                     //how many items to be received
                MPI_Recv(&size, 1, MPI_INT, 0, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                double recv[size];            //array to hold received items
                MPI_Recv(recv, size, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                for(j=0; j<size; j++)
                {
                    my_vec[j] = recv[j];
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function broadcasts a 1-d vector of ints                                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void broadcast_iv1d(int world_size,int world_rank,iv1d &my_vec)
{
    int i = 0;
    int j = 0;

    if(world_size > 0)
    {
        for(i=1; i<world_size; i++)
        {
            if(world_rank == 0)
            {
                int size = my_vec.size();  //how many items to send

                MPI_Send(&size, 1, MPI_INT, i, 13, MPI_COMM_WORLD);

                int snd[size];
                for(j=0; j<size; j++)
                {
                    snd[j] = my_vec[j];
                }
                MPI_Send(snd, size, MPI_INT, i, 13, MPI_COMM_WORLD);
            }
            else if(world_rank == i)
            {
                int size;                     //how many items to be received
                MPI_Recv(&size, 1, MPI_INT, 0, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                int recv[size];            //array to hold received items
                MPI_Recv(recv, size, MPI_INT, 0, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                for(j=0; j<size; j++)
                {
                    my_vec[j] = recv[j];
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects a 2D vectors of long ints to make a 3D vec                                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void collect_lv3d(int world_size,int world_rank,lv3d &my_vec)
{
    int i = 0;
    int j = 0;
    int k = 0;

    if(world_size > 0)
    {
        for(i=1; i<world_size; i++)
        {
            if(world_rank == 0)
            {
                int size_x = my_vec[0].size();    //how many items to receive in x
                int size_y = my_vec[0][0].size(); //how many items to receive in y

                long recv[size_x*size_y];     //array to hold received items
                MPI_Recv(recv, size_x*size_y, MPI_LONG, i, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                int counter = 0;

                //add the items to the vector
                for(j=0; j<size_x; j++)
                {
                    for(k=0; k<size_y; k++)
                    {
                        my_vec[i][j][k] = recv[counter];
                        counter++;
                    }
                }
            }
            else if(world_rank == i)
            {
                int size_x = my_vec[0].size();    //how many items to send in x
                int size_y = my_vec[0][0].size(); //how many items to send in y

                long snd[size_x*size_y];
                int counter = 0;
                for(j=0; j<size_x; j++)
                {
                    for(k=0; k<size_y; k++)
                    {
                        snd[counter] = my_vec[i][j][k];
                        counter++;
                    }
                }
                MPI_Send(snd, size_x*size_y, MPI_LONG, 0, 13, MPI_COMM_WORLD);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function broadcasts a 3D vector of long ints                                                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void broadcast_lv3d(int world_size,int world_rank,lv3d &my_vec)
{
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;

    for(i=1; i<world_size; i++)
    {
        if(world_rank == 0) //send data
        {
            int size_x = my_vec.size(); 
            int size_y = my_vec[0].size();
            int size_z = my_vec[0][0].size();

            long snd[size_x*size_y*size_z]; 

            int counter = 0;

            for(j=0; j<size_x; j++)
            {
                for(k=0; k<size_y; k++)
                {
                    for(l=0; l<size_z; l++)
                    {
                        snd[counter] = my_vec[j][k][l];
                        counter++;
                    }
                }
            }
            MPI_Send(snd, size_x*size_y*size_z, MPI_LONG, i, 13, MPI_COMM_WORLD);
        }
        else if(i==world_rank) //receive data
        {
            int size_x = my_vec.size();
            int size_y = my_vec[0].size();
            int size_z = my_vec[0][0].size();

            long recv[size_x*size_y*size_z];     //array to hold received items

            MPI_Recv(recv, size_x*size_y*size_z, MPI_LONG, 0, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            int counter = 0;
            for(j=0; j<size_x; j++)
            {
                for(k=0; k<size_y; k++)
                {
                    for(l=0; l<size_z; l++)
                    {
                       my_vec[j][k][l] = recv[counter];
                       counter++;
                    }
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects an integer and sums it                                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void collect_and_sum_int(int world_size,int world_rank,int *my_val)
{
    int i = 0;
    int j = 0;

    if(world_size > 0)
    {
        for(i=1; i<world_size; i++)
        {
            if(world_rank == 0)
            {
                int rcv_val = 0;
                MPI_Recv(&rcv_val, 1, MPI_INT, i, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                *my_val = *my_val + rcv_val;
            }
            else if(world_rank == i)
            {
                MPI_Send(my_val, 1, MPI_INT, 0, 13, MPI_COMM_WORLD);
            }
        }
    }
}

/*cannot compile some programs with this
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects a 2-d vector of doubles                                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void broadcast_rvec(int world_size,int world_rank,rvec *my_rvec,int size)
{
    int i = 0;
    int j = 0;

    if(world_size > 0)
    {
        for(i=1; i<world_size; i++)
        {
            if(world_rank == 0)
            {
                double snd[size];

                //send x
                for(j=0; j<size; j++)
                {
                    snd[j] = my_rvec[j][0];
                }
                MPI_Send(snd, size, MPI_DOUBLE, i, 13, MPI_COMM_WORLD);

                //send y
                for(j=0; j<size; j++)
                {
                    snd[j] = my_rvec[j][1];
                }
                MPI_Send(snd, size, MPI_DOUBLE, i, 13, MPI_COMM_WORLD);

                //send z
                for(j=0; j<size; j++)
                {
                    snd[j] = my_rvec[j][2];
                }
                MPI_Send(snd, size, MPI_DOUBLE, i, 13, MPI_COMM_WORLD);

            }
            else if(world_rank == i)
            {
                double recv[size];            //array to hold received items

                //receive x
                MPI_Recv(recv, size, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for(j=0; j<size; j++)
                {
                    my_rvec[j][0] = recv[j];
                }

                //receive y
                MPI_Recv(recv, size, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for(j=0; j<size; j++)
                {
                    my_rvec[j][1] = recv[j];
                }

                //receive z
                MPI_Recv(recv, size, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for(j=0; j<size; j++)
                {
                    my_rvec[j][2] = recv[j];
                }
            }
        }
    }
}
*/

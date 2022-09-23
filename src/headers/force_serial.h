
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function checks the parallelization scheme and adjusts my_frames accordingly.                       //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void check_serial(int world_rank,int world_size,enum Switch serial)
{
    if(serial == 1 && world_size > 1)
    {
        if(world_rank == 0)
        {
            printf("This program does not support parallelization. Please set the number of cores to 1 and try again. \n");
        }
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }
}


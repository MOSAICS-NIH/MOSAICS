
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function takes a filename and check that it contains the correct extension                          //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void check_extension_mpi(int world_rank,string arg,string filename,string tag)
{
    int i             = 0;                  //standard variable used in loops
    int bad_extension = 0;                  //tells if the correct extension was provided
    int filename_size = filename.length();  //how long is the filename
    int tag_size      = tag.length();       //how long is the tag

    //make sure tag is not longer than filename
    if(tag_size > filename_size)
    {
        bad_extension = 1;
    }
    else //check for proper extension
    {
        for(i=filename_size-tag_size; i<filename_size; i++)
        {
            if(filename[i] != tag[i + tag_size - filename_size])
            {
                bad_extension = 1;
                break;
            }
        }
    }

    //report error and terminate program
    if(bad_extension == 1)
    {
        if(world_rank == 0)
        {
            printf("The filename %s provided via the %s tag requires the %s extension. \n",filename.c_str(),arg.c_str(),tag.c_str());
        }
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }
}


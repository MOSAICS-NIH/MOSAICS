
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints info useful for debugging                                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void dbug(string tag,int world_rank)
{
    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(stderr,"world rank %5d tag %s \n",world_rank,tag.c_str());
    MPI_Barrier(MPI_COMM_WORLD);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This is a class for debugging code                                                                        //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Debug
{
    private:
        int tag        = 0;      //Tag to print


    public:
        void print_tag();   //print tag and increment
};

void Debug::print_tag()
{
    printf("tag %d \n",tag);
    tag = tag + 1;
}


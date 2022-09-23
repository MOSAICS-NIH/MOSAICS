
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This structure holds variables used throughout the code. Any variable that can be declared at the start   //
// of the program is stored here.                                                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct system_variables
{
    string  program_name;                         //Name of your analysis program
    int     world_size;                           //How many MPI ranks make the world
    int     world_rank;                           //Rank of a process
    int     counter;                              //How many times the "program run time" been displayed
    double  seconds;                              //Used to keep track of how long the program has been running
    clock_t t;                                    //Keeps the time for testing performance
    clock_t t0;                                   //Takes the time at the start
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function initializes the variables held in system_variables                                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void initialize_system_variables(system_variables *s)
{
    s->counter = 0;
    s->seconds = 0;
    s->t  = clock();
    s->t0 = clock();

    MPI_Comm_size(MPI_COMM_WORLD, &s->world_size);      //get the world size
    MPI_Comm_rank(MPI_COMM_WORLD, &s->world_rank);      //get the process rank
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints the credits for the program. Also prints a record of the input arguments.            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void name_and_record(int world_rank,int argc, const char *argv[],string program_name)
{
    int i = 0;
    if(world_rank == 0)
    {
        //set the length of the line
        int line_length = 130;

        //compute the padding needed
        int length = program_name.length();
        int lead_pad = ceil((line_length-4-length)/2);
        int trail_pad = line_length - 4 - length - lead_pad;

        cout << "\n";

        //line of slashes
        for(i=0; i<line_length; i++)
        {
            cout << "/";
        }
        cout << "\n";

        //2 slashes then spaces with 2 slashes at end
        cout << "//";
        for(i=0; i<line_length-4; i++)
        {
            cout << " ";
        }
        cout << "//";
        cout << "\n";

        //line with program name
        cout << "//";
        for(i=0; i<lead_pad; i++)
        {
            cout << " ";
        }
        cout << program_name;
        for(i=0; i<trail_pad; i++)
        {
            cout << " ";
        }
        cout << "//";
        cout << "\n";

        //2 slashes then spaces with 2 slashes at end
        cout << "//";
        for(i=0; i<line_length-4; i++)
        {
            cout << " ";
        }
        cout << "//";
        cout << "\n";

        //line of slashes
        for(i=0; i<line_length; i++)
        {
            cout << "/";
        }
        cout << "\n";

        //now we print the command line arguments for records
        printf("Input: ");
        for(i=0; i<argc; i++)
        {
            printf("%s ",argv[i]);
        }
        printf("\n\n");
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints final statistics for the program                                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void print_closing(int world_rank)
{
    MPI_Barrier(MPI_COMM_WORLD);

    if(world_rank == 0)
    {
        printf("\nAnalysis completed successfully.\n\n");
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function determines how much time has passed and gives an estimate of the time remaining.            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void time_stats(clock_t t,int *counter,int current_frame,int my_frames,int world_rank)
{
    if(world_rank == 0)
    {
        double seconds = 0;

        seconds = (clock() - t)/CLOCKS_PER_SEC;

        if((int)seconds/10 > *counter)
        {
            double percent_done = ((double)(current_frame+1)/(double)my_frames)*100;
            double estimated_total_time = 100*seconds/percent_done;
            double time_remaining = estimated_total_time - seconds;

            if(current_frame+1 == my_frames)
            {
                estimated_total_time = 0;
            }

            int phr = 0;
            int pmin = 0;
            int psec = 0;
            int lhr = 0;
            int lmin = 0;
            int lsec = 0;

            phr = (seconds)/(60*60);
            pmin = (seconds - (phr*60*60))/60;
            psec = seconds - (phr*60*60) - (pmin*60);

            lhr = ((int)time_remaining)/(60*60);
            lmin = (time_remaining - (lhr*60*60))/60;
            lsec = time_remaining - (lhr*60*60) - (lmin*60);

            if(percent_done != 0)
            {
                printf("Finished frame %7d with %5.1f percent done overall in %2d hr %2d min %2d sec. Estimated time to completion is %2d hr %2d min %2d sec. \n",current_frame+1,percent_done,phr,pmin,psec,lhr,lmin,lsec);
            }
            *counter = *counter + 1;
        }
    }
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This is a class for logging performance data throughout the program                                      //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Performance
{
    private:
        sv1d    title;
        dv1d    time;
        int     world_rank; 
        int     world_size;

    public:
        void log_time(double my_time,string my_title);
        void print_stats();
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function logs the time spent in a section                                                           //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Performance::log_time(double my_time,string my_title)
{
    time.push_back(my_time);
    title.push_back(my_title);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function prints the final performance statistics for each section                                   //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Performance::print_stats()
{
    int    i         = 0;
    int    hr        = 0;
    int    min       = 0;
    int    sec       = 0;
    int    total_sec = 0;
    int    hr_tot    = 0;
    int    min_tot   = 0;
    int    sec_tot   = 0;
    double percent   = 0;

    //get the world size and assign ranks
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);      //get the world size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);      //get the process rank

    if(world_rank == 0)
    {
        //print header information
        printf("\nPerformance stats:\n");
        printf("%20s %10s %10s %10s %10s \n","Routine","Hours","Minutes","Seconds","Percent");
        printf("%21s%11s%11s%11s%11s\n","---------------------","-----------","-----------","-----------","-----------");

        //compute the total time spent 
        for(i=0; i<title.size(); i++)
        {
            total_sec = total_sec + time[i];
        }

        //now compute hours, minutes, seconds etc. in each section
        for(i=0; i<title.size(); i++)
        {
            percent = 100*time[i]/total_sec;
            hr      = (time[i])/(60*60);
            min     = (time[i] - (hr*60*60))/60;
            sec     = (time[i]) - (hr*60*60) - (min*60);

            printf("%20s %10d %10d %10d %10.1f \n",title[i].c_str(),hr,min,sec,percent);
        }

        //compute the total time (hr,min,sec)
        hr_tot  = (total_sec)/(60*60);
        min_tot = (total_sec - (hr_tot*60*60))/60;
        sec_tot = (total_sec) - (hr_tot*60*60) - (min_tot*60);

        printf("%20s %10d %10d %10d %10.1f \n","Total",hr_tot,min_tot,sec_tot,100.0);
    }
}

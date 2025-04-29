
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function determines how much time has passed and gives an estimate of the time remaining.            //
// Takes in the current step starting at 1 (not zero).                                                       //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void be_time_stats(clock_t t,int *counter,int current_step,int my_steps,int world_rank,string my_string)
{
    if(world_rank == 0)
    {
        double seconds = 0;

        seconds = (clock() - t)/CLOCKS_PER_SEC;

        if((int)seconds/10 > *counter)
        {
            double percent_done         = ((double)(current_step)/(double)my_steps)*100.0;
            double estimated_total_time = 100.0*seconds/percent_done;
            double time_remaining       = estimated_total_time - seconds;

            if(current_step == my_steps)
            {
                estimated_total_time = 0.0;
            }

            int phr  = 0;
            int pmin = 0;
            int psec = 0;
            int lhr  = 0;
            int lmin = 0;
            int lsec = 0;

            phr  = (seconds)/(60*60);
            pmin = (seconds - (phr*60*60))/60;
            psec = seconds - (phr*60*60) - (pmin*60);

            lhr  = ((int)time_remaining)/(60*60);
            lmin = (time_remaining - (lhr*60*60))/60;
            lsec = time_remaining - (lhr*60*60) - (lmin*60);

            if(percent_done != 0)
            {
                printf("Finished %s %7d with %5.1f percent done overall in %2d hr %2d min %2d sec. Estimated time to completion is %2d hr %2d min %2d sec. \n",my_string.c_str(),current_step+1,percent_done,phr,pmin,psec,lhr,lmin,lsec);
            }
            *counter = *counter + 1;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function takes in binding events data and characterizes repeat bindings and makes a list of binding  //
// lipids (no repeats).                                                                                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void reduce_binding_list(iv1d &lipid_nr,iv1d &res_nr,sv1d &res_name,iv1d &bind_i,iv1d &bind_f,iv1d &dwell_t,
                         iv1d &repeats_res_nr,iv1d &repeats,iv1d &repeats_lipid_nr,dv1d &repeats_avg_dwell_time,
                         dv1d &repeats_avg_off_time,sv1d &repeats_res_name)
{
    int i = 0;                              //standard variable used in loops
    int j = 0;                              //standard variable used in loops
    iv1d prev_bind_f(iv1d(0,0));            //time lipid last left grid point

    //resize repeats vectors
    repeats.resize(0,0);
    repeats_res_nr.resize(0,0);
    repeats_lipid_nr.resize(0,0);
    repeats_res_name.resize(0);
    repeats_avg_dwell_time.resize(0,0.0);
    repeats_avg_off_time.resize(0,0.0);

    for(i=0; i<res_nr.size(); i++) //loop over binding events
    {
        int match = 0;
        for(j=0; j<repeats_res_nr.size(); j++) //loop over list of binders
        {
            if(repeats_res_nr[j] == res_nr[i]) //lipid already encountered.
            {
                repeats[j]                = repeats[j] + 1;
                repeats_avg_dwell_time[j] = repeats_avg_dwell_time[j] + (double)dwell_t[i];
                repeats_avg_off_time[j]   = repeats_avg_off_time[j] + (bind_i[i] - prev_bind_f[j]);
                prev_bind_f[j]            = bind_f[i];

                match = 1;
            }
        }

        if(match == 0) //first encounter
        {
            repeats.push_back(0);
            repeats_lipid_nr.push_back(lipid_nr[i]);
            repeats_res_nr.push_back(res_nr[i]);
            repeats_res_name.push_back(res_name[i]);
            repeats_avg_dwell_time.push_back((double)dwell_t[i]);
            repeats_avg_off_time.push_back(0.0);
            prev_bind_f.push_back(bind_f[i]);
        }
    }

    //compute average dwell time and average off times for each lipid
    for(i=0; i<repeats_res_nr.size(); i++) //loop over lipid binders
    {
        if(repeats[i] > 0) //avoid division by zero
        {
            repeats_avg_dwell_time[i] = repeats_avg_dwell_time[i]/(double)(repeats[i] + 1);
            repeats_avg_off_time[i]   = repeats_avg_off_time[i]/(double)(repeats[i]);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This is a class for working with a binding events file                                                    //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Binding_events
{
    private:

    public:
        int x_i          = 0;                                                                                                            //x grid point
        int y_i          = 0;                                                                                                            //y grid point
        int ef_frames    = 0;                                                                                                            //effective number of frames in 2d_kinetics (accounts for stride)
        int num_lipids   = 0;                                                                                                            //number of selected lipids in 2d_kinetics
        int num_g_x      = 0;                                                                                                            //number of grid points in x in 2d_kinetics
        int num_g_y      = 0;                                                                                                            //number of grid points in y in 2d_kinetics
        double ef_dt     = 0.0;                                                                                                          //effective dt in 2d_kinetics (accounts for stride)
        double APS       = 0;                                                                                                            //area per square used for grid in 2d_kinetics

        //lipid mixing stuff
        int lip_nr_1    = 0;
        int res_nr_1    = 0;
        int num_lip_1   = 0;
        int num_lip_2   = 0;
      
        //binding events
        iv1d lipid_nr{};                                                                                                                 //what is the lipid number (0-max), note the lipid number is not unique if combining binding events files 
        iv1d bind_i{};                                                                                                                   //what time did the lipid first bind to the grid point
        iv1d bind_f{};                                                                                                                   //what frame did the lipid leave the grid point
        iv1d dwell_t{};                                                                                                                  //dwell time for the binding event
        iv1d res_nr{};                                                                                                                   //residue number for the binding event
        sv1d res_name{};                                                                                                                 //residue name for the binding event

        //repeated visits
        iv1d repeats_res_nr{};                                                                                                           //residue number for repeated visits list
        iv1d repeats{};                                                                                                                  //number of repeat visits for each lipid
        iv1d repeats_lipid_nr{};                                                                                                         //lipid number for repeated visits list
        sv1d repeats_res_name{};                                                                                                         //name of residue for repeats
        dv1d repeats_avg_dwell_time{};                                                                                                   //average dwell time for repeated visits list
        dv1d repeats_avg_off_time{};                                                                                                     //average time between repeated visits

        //time line 
        iv2d bound_time_line{};                                                                                                          //stores if the lipid bound at each frame
        iv1d time_line_lipid_nr{};                                                                                                       //stores time line lipid_nr
        iv1d time_line_res_nr{};                                                                                                         //stores time line res_rn
        sv1d time_line_res_name{};                                                                                                       //stores time line res_name

        //lipid tessellations
        iv2d voro;                                                                                                                       //grid containing tessellations
        iv2d voro_nan;                                                                                                                   //grid containing tessellations nan info
        iv2d voro_frame;                                                                                                                 //current frame tessellation data
        iv2d voro_nan_frame;                                                                                                             //current frame tessellation nan data


        //binary be files
        ofstream be_file_o;                                                                                                              //file for writing binding events
        ifstream be_file_i;                                                                                                              //file for reading binding events

        //info file for gid data
        i64v2d info_pos;                                                                                                                 //hols the position of binding events data in the be file

    public:
        void    gen_lipid_nr();                                                                                                             //generates a unique lipid nr for each lipid
        void    organize_events(int mode);                                                                                                  //organize binding events
        void    reduce_list();                                                                                                              //characterize repeat visits
        void    size_timeline();                                                                                                            //allocate memory for the binding time line
        void    get_binding_timeline();                                                                                                     //converts binding events into a timeline
        void    add_to_binding_timeline();                                                                                                  //add binding events to time line 
        void    stamp_to_binding_timeline();                                                                                                //stamps binding events to existing timeline
        void    binding_events_from_timeline();                                                                                             //use binding timeline to create binding events
        void    write_time_line(string timeline_file_name);                                                                                 //write the binding timeline to file
        void    suppress_timeline_noise(int threshold);                                                                                     //mend fragmented binding events
        void    suppress_timeline_noise_be(int threshold);                                                                                  //mend fragmented binding events
	void    sweep_timeline_noise(int lipid_count);                                                                                      //constrain numbero of lipids bound in timeline for a given t
        void    get_complete_set(string in_file_name,iv1d &lipid_nr,iv1d &res_nr,sv1d &res_name);                                           //reads binding events file until a complete set of lipid_nr etc. is found
        dv1d    find_voro_center(int this_lipid);                                                                                           //returns the center of a lipid tessellation
        int     get_binding_events_grid(string in_file_name,i64v2d &this_be_pos,int this_x,int this_y);                                     //reads in a binding events file for a lattice point
        int     add_binding_events_grid(string in_file_name,i64v2d &this_be_pos,int this_x,int this_y);                                     //reads the binding events file for a lattice point and adds to current list
        void    write_binding_events_bin(string binding_events_file_name);                                                                  //writes out a binding events file in binary
        int     get_binding_events_bin(string binding_events_file_name);                                                                    //reads in a binary finding events file
        int64_t write_binding_events_tmp(ofstream &be_file_o,int64_t pos);                                                                  //writes a be file to the tmp file for an MPI core. 
        int64_t get_binding_events_tmp(ifstream &be_file_i,int64_t pos);                                                                    //reads binding events from the temporary be file
        int     get_info(string binding_events_file_name);                                                                                  //reads the info file for a binding evetns file (grid)
        int     get_binding_events_xy(string binding_events_file_name,int x,int y);                                                         //reads binding events data for a grid point 
        int     get_tessellations(string in_file_name,int my_gi,int my_gf,int my_num_g,int stride,int my_ef_frames,int world_rank);         //reads binding events file and extracts tessellations for each frame
        void    get_voro_frame(int this_frame,int world_size,int world_rank);                                                               //collects a single frame of the lipid tessellation
        void    write_voro_frame(string out_file_name);                                                                                     //writes the current frame tessellation data to file
        int     add_binding_events(string binding_events_file_name,int compound_lipid_count,int lipid_nr_offset);                           //reads the binding events file and adds to current list
        void    write_binding_events_legacy(string binding_events_file_name);                                                               //writes a binding events file in human readable form
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function creates a unique lipid number for each lipid                                                //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Binding_events::gen_lipid_nr()
{
    int i      = 0;          //standard variable used in loops
    int b_swap = 0;          //Used to tell if aray requires further sorting

    //organize by res_nr
    for(b_swap=1; b_swap > 0; )
    {
        b_swap = 0;
        for(i=0; i<lipid_nr.size()-1; i++)
        {
            if(res_nr[i] > res_nr[i+1])
            {
                int dwell_t_tmp     = dwell_t[i];
                int lipid_nr_tmp    = lipid_nr[i];
                int bind_i_tmp      = bind_i[i];
                int bind_f_tmp      = bind_f[i];
                int res_nr_tmp      = res_nr[i];
                string res_name_tmp = res_name[i];

                dwell_t[i]  = dwell_t[i+1];
                lipid_nr[i] = lipid_nr[i+1];
                bind_i[i]   = bind_i[i+1];
                bind_f[i]   = bind_f[i+1];
                res_nr[i]   = res_nr[i+1];
                res_name[i] = res_name[i+1];

                dwell_t[i+1]  = dwell_t_tmp;
                lipid_nr[i+1] = lipid_nr_tmp;
                bind_i[i+1]   = bind_i_tmp;
                bind_f[i+1]   = bind_f_tmp;
                res_nr[i+1]   = res_nr_tmp;
                res_name[i+1] = res_name_tmp;

                b_swap = 1;
            }
        }
    }

    //create unique lipid_nr for each lipid
    int prev_res_nr = -1;
    int lipid_count = -1;
    for(i=0; i<lipid_nr.size()-1; i++)
    {
        if(res_nr[i] != prev_res_nr)
        {
            prev_res_nr = res_nr[i];
            lipid_count++;
        }
        lipid_nr[i] = lipid_count;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function organizes the binding events                                                                //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Binding_events::organize_events(int mode)
{
    int i      = 0;          //standard variable used in loops
    int b_swap = 0;          //Used to tell if aray requires further sorting

    if(mode == 0) //order by initial binding time
    {
        for(b_swap=1; b_swap > 0; )
        {
            b_swap = 0;
            for(i=0; i<lipid_nr.size()-1; i++)
            {
                if(bind_i[i] > bind_i[i+1])
                {
                    int dwell_t_tmp     = dwell_t[i];
                    int lipid_nr_tmp    = lipid_nr[i];
                    int bind_i_tmp      = bind_i[i];
                    int bind_f_tmp      = bind_f[i];
                    int res_nr_tmp      = res_nr[i];
                    string res_name_tmp = res_name[i];

                    dwell_t[i]  = dwell_t[i+1];
                    lipid_nr[i] = lipid_nr[i+1];
                    bind_i[i]   = bind_i[i+1];
                    bind_f[i]   = bind_f[i+1];
                    res_nr[i]   = res_nr[i+1];
                    res_name[i] = res_name[i+1];

                    dwell_t[i+1]  = dwell_t_tmp;
                    lipid_nr[i+1] = lipid_nr_tmp;
                    bind_i[i+1]   = bind_i_tmp;
                    bind_f[i+1]   = bind_f_tmp;
                    res_nr[i+1]   = res_nr_tmp;
                    res_name[i+1] = res_name_tmp;

                    b_swap = 1;
                }
            }
        }
    }
    else if(mode == 1) //order by dwell time
    {
        //sort events with longest dwell time coming first
        for(b_swap=1; b_swap > 0; )
        {
            b_swap = 0;
            for(i=0; i<lipid_nr.size()-1; i++)
            {
                if(dwell_t[i] < dwell_t[i+1])
                {
                    int dwell_t_tmp     = dwell_t[i];
                    int lipid_nr_tmp    = lipid_nr[i];
                    int bind_i_tmp      = bind_i[i];
                    int bind_f_tmp      = bind_f[i];
                    int res_nr_tmp      = res_nr[i];
                    string res_name_tmp = res_name[i];

                    dwell_t[i]  = dwell_t[i+1];
                    lipid_nr[i] = lipid_nr[i+1];
                    bind_i[i]   = bind_i[i+1];
                    bind_f[i]   = bind_f[i+1];
                    res_nr[i]   = res_nr[i+1];
                    res_name[i] = res_name[i+1];

                    dwell_t[i+1]  = dwell_t_tmp;
                    lipid_nr[i+1] = lipid_nr_tmp;
                    bind_i[i+1]   = bind_i_tmp;
                    bind_f[i+1]   = bind_f_tmp;
                    res_nr[i+1]   = res_nr_tmp;
                    res_name[i+1] = res_name_tmp;

                    b_swap = 1;
                }
            }
        }
    }
    else if(mode == 2) //print by repeated visits
    {
        for(b_swap=1; b_swap > 0; )
        {
            b_swap = 0;
            for(i=0; i<repeats_res_nr.size()-1; i++)
            {
                if(repeats[i] < repeats[i+1])
                {
                    int repeats_res_nr_tmp             = repeats_res_nr[i];
                    int repeats_tmp                    = repeats[i];
                    int repeats_lipid_nr_tmp           = repeats_lipid_nr[i];
                    double repeats_avg_dwell_time_tmp  = repeats_avg_dwell_time[i];
                    double repeats_avg_off_time_tmp    = repeats_avg_off_time[i];
                    string repeats_res_name_tmp        = repeats_res_name[i];

                    repeats_res_nr[i]                  = repeats_res_nr[i+1];
                    repeats[i]                         = repeats[i+1];
                    repeats_lipid_nr[i]                = repeats_lipid_nr[i+1];
                    repeats_avg_dwell_time[i]          = repeats_avg_dwell_time[i+1];
                    repeats_avg_off_time[i]            = repeats_avg_off_time[i+1];
                    repeats_res_name[i]                = repeats_res_name[i+1];

                    repeats_res_nr[i+1]                = repeats_res_nr_tmp;
                    repeats[i+1]                       = repeats_tmp;
                    repeats_lipid_nr[i+1]              = repeats_lipid_nr_tmp;
                    repeats_avg_dwell_time[i+1]        = repeats_avg_dwell_time_tmp;
                    repeats_avg_off_time[i+1]          = repeats_avg_off_time_tmp;
                    repeats_res_name[i+1]              = repeats_res_name_tmp;

                    b_swap = 1;
                }
            }
        }
    }
    else if(mode == 3) //sort by res_nr
    {
        for(b_swap=1; b_swap > 0; )
        {
            b_swap = 0;
            for(i=0; i<repeats_res_nr.size()-1; i++)
            {
                if(repeats_res_nr[i] < repeats_res_nr[i+1])
                {
                    int repeats_res_nr_tmp          = repeats_res_nr[i];
                    int repeats_tmp                 = repeats[i];
                    int repeats_lipid_nr_tmp        = repeats_lipid_nr[i];
                    int repeats_avg_dwell_time_tmp  = repeats_avg_dwell_time[i];
                    double repeats_avg_off_time_tmp = repeats_avg_off_time[i];
                    string repeats_res_name_tmp     = repeats_res_name[i];

                    repeats_res_nr[i]               = repeats_res_nr[i+1];
                    repeats[i]                      = repeats[i+1];
                    repeats_lipid_nr[i]             = repeats_lipid_nr[i+1];
                    repeats_avg_dwell_time[i]       = repeats_avg_dwell_time[i+1];
                    repeats_avg_off_time[i]         = repeats_avg_off_time[i+1];
                    repeats_res_name[i]             = repeats_res_name[i+1];

                    repeats_res_nr[i+1]             = repeats_res_nr_tmp;
                    repeats[i+1]                    = repeats_tmp;
                    repeats_lipid_nr[i+1]           = repeats_lipid_nr_tmp;
                    repeats_avg_dwell_time[i+1]     = repeats_avg_dwell_time_tmp;
                    repeats_avg_off_time[i+1]       = repeats_avg_off_time_tmp;
                    repeats_res_name[i+1]           = repeats_res_name_tmp;

                    b_swap = 1;
                }
            }
        }
    }
    else if(mode == 4) //sort by final binding time
    {
        for(b_swap=1; b_swap > 0; )
        {
            b_swap = 0;
            for(i=0; i<lipid_nr.size()-1; i++)
            {
                if(bind_f[i] > bind_f[i+1])
                {
                    int dwell_t_tmp     = dwell_t[i];
                    int lipid_nr_tmp    = lipid_nr[i];
                    int bind_i_tmp      = bind_i[i];
                    int bind_f_tmp      = bind_f[i];
                    int res_nr_tmp      = res_nr[i];
                    string res_name_tmp = res_name[i];

                    dwell_t[i]  = dwell_t[i+1];
                    lipid_nr[i] = lipid_nr[i+1];
                    bind_i[i]   = bind_i[i+1];
                    bind_f[i]   = bind_f[i+1];
                    res_nr[i]   = res_nr[i+1];
                    res_name[i] = res_name[i+1];

                    dwell_t[i+1]  = dwell_t_tmp;
                    lipid_nr[i+1] = lipid_nr_tmp;
                    bind_i[i+1]   = bind_i_tmp;
                    bind_f[i+1]   = bind_f_tmp;
                    res_nr[i+1]   = res_nr_tmp;
                    res_name[i+1] = res_name_tmp;

                    b_swap = 1;
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function takes in binding evets vectors and adds to a binding timeline vector                        //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Binding_events::add_to_binding_timeline()
{
    int i = 0;
    int j = 0;

    for(i=0; i<res_nr.size(); i++) //loop over binding events
    {
        for(j=bind_i[i]; j<=bind_f[i]; j++) //loop over bound frames
        {
            time_line_lipid_nr[lipid_nr[i]] = lipid_nr[i];
            time_line_res_nr[lipid_nr[i]]   = res_nr[i];
            time_line_res_name[lipid_nr[i]] = res_name[i];

            bound_time_line[j][lipid_nr[i]] = 1;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function takes in binding evets vectors and stamps to a binding timeline vector                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Binding_events::stamp_to_binding_timeline()
{
    int i = 0;
    int j = 0;

    for(i=0; i<res_nr.size(); i++) //loop over binding events
    {
        for(j=bind_i[i]; j<=bind_f[i]; j++) //loop over bound frames
        {
            time_line_lipid_nr[lipid_nr[i]] = lipid_nr[i];
            time_line_res_nr[lipid_nr[i]]   = res_nr[i];
            time_line_res_name[lipid_nr[i]] = res_name[i];

            bound_time_line[j][lipid_nr[i]] = bound_time_line[j][lipid_nr[i]] + 1;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function sizes the timeline vectors based on num_lipids and ef_frames                                //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Binding_events::size_timeline()
{
    bound_time_line.resize(ef_frames,iv1d(num_lipids,0));
    time_line_lipid_nr.resize(num_lipids,0); 
    time_line_res_nr.resize(num_lipids,0);
    time_line_res_name.resize(num_lipids);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function uses information stored in the timeline to populate the binding events                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Binding_events::binding_events_from_timeline()
{
    int i = 0;
    int j = 0;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Resize binding events vectors                                                                             //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    lipid_nr.resize(0,0);
    res_nr.resize(0,0);
    res_name.resize(0);
    bind_i.resize(0,0);
    bind_f.resize(0,0);
    dwell_t.resize(0,0.0);

    for(i=0; i<num_lipids; i++) //loop over lipids
    {
        int bind_initial = -1;
        int bind_final   = -1;
        int prev_state   = -1;
        for(j=0; j<ef_frames; j++) //loop over time line frame
        {
            if(bound_time_line[j][i] == 1 && prev_state != 1) //lipid just bound
            {
                bind_initial = j;
            }
            else if( (bound_time_line[j][i] == 0 || j==ef_frames-1) && prev_state == 1) //lipid just left, write event to file
            {
                bind_final = j-1; //leaves last frame undound (careful if splicing)
                int dwell_time = bind_final + 1 - bind_initial;

                lipid_nr.push_back(time_line_lipid_nr[i]);
                res_nr.push_back(time_line_res_nr[i]);
                res_name.push_back(time_line_res_name[i]);
                bind_i.push_back(bind_initial);
                bind_f.push_back(bind_final);
                dwell_t.push_back(dwell_time);
            }
            prev_state = bound_time_line[j][i];
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes the binding timeline to file                                                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Binding_events::write_time_line(string timeline_file_name)
{
    int i = 0;
    int j = 0;

    FILE *timeline_file = fopen(timeline_file_name.c_str(), "w");
    if(timeline_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",timeline_file_name.c_str());
    }
    else
    {
        for(i=0; i<num_lipids; i++) //loop over lipids
        {
            for(j=0; j<ef_frames; j++) //loop over time line frame
            {
                fprintf(timeline_file," %1d ",bound_time_line[j][i]);
            }
            fprintf(timeline_file,"\n");
        }
        fclose(timeline_file);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function takes in binding events vectors and fills binding timeline vectors                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Binding_events::get_binding_timeline()
{
    int i = 0;
    int j = 0;

    //resize vectors
    bound_time_line.resize(0,iv1d(0,0));
    time_line_lipid_nr.resize(0,0);
    time_line_res_nr.resize(0,0);
    time_line_res_name.resize(0);

    //allocate memory for vectors
    bound_time_line.resize(ef_frames,iv1d(num_lipids,0));
    time_line_lipid_nr.resize(num_lipids,0);
    time_line_res_nr.resize(num_lipids,0);
    time_line_res_name.resize(num_lipids);

    //now add bound lipids to the time line
    for(i=0; i<res_nr.size(); i++) //loop over binding events
    {
        for(j=bind_i[i]; j<=bind_f[i]; j++) //loop over bound frames
        {
            time_line_lipid_nr[lipid_nr[i]] = lipid_nr[i];
            time_line_res_nr[lipid_nr[i]]   = res_nr[i];
            time_line_res_name[lipid_nr[i]] = res_name[i];

            bound_time_line[j][lipid_nr[i]] = 1;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function characterizes repeat bindings                                                               //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Binding_events::reduce_list()
{
    reduce_binding_list(lipid_nr,res_nr,res_name,bind_i,bind_f,dwell_t,repeats_res_nr,repeats,repeats_lipid_nr,
                        repeats_avg_dwell_time,repeats_avg_off_time,repeats_res_name);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function takes a binding timeline and fills in fragmentation gaps                                    //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Binding_events::suppress_timeline_noise(int threshold)
{
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int prev_state = 0;

    if(threshold > 0)
    {
        for(i=0; i<num_lipids; i++) //loop over time line lipids
        {
            for(j=0; j<ef_frames; j++) //loop over time line frames
            {
                if(bound_time_line[j][i] == 0 && prev_state == 1) //lipid just left
                {
                    //determine max number of frames to check
                    int max = j+threshold;
                    if(max > ef_frames) //cant read off end of array
                    {
                        max = ef_frames;
                    }

                    for(k=j; k<max; k++) //check if lipid rebinds
                    {
                        if(bound_time_line[k][i] == 1) //lipid rebinds
                        {
                            for(l=j; l<k; l++) //mend break
                            {
                                bound_time_line[l][i] = 1;
                            }
                            break;
                        }
                    }
                }
                prev_state = bound_time_line[j][i];
            }
        }
    }
}    

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function bypasses a time line while filling in fragmentation gaps                                    //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Binding_events::suppress_timeline_noise_be(int threshold)
{
    int i = 0;
    int j = 0;

    if(threshold > 0)
    {
        iv1d removed(lipid_nr.size(), 0);

        for(i=0; i<lipid_nr.size(); i++) //loop over leaving lipids
        {
            if(removed[i] == 0)
            {
                for(j=0; j<lipid_nr.size(); j++) //loop over incoming lipids
                {
                    if(lipid_nr[i] == lipid_nr[j]) //same lipid
                    {
                        if(bind_i[j] >= bind_i[i] && bind_f[j] <= bind_f[i] && i != j) //event falls inside another event
                        {
                            removed[j] = 1;
                        }

                        if(removed[j] == 0)
                        {
                            if(bind_i[j] > bind_f[i] && bind_i[j] <= (bind_f[i] + threshold)) 
                            {
                                bind_f[i]  = bind_f[j];
                                dwell_t[i] = bind_f[j] + 1 - bind_i[i];
                                removed[j] = 1;
                                i = i - 1;
                                goto next_iter;
                            }
                        }
                    }
                }
                next_iter:;
            }
        }

        //make structures to store new binding events info	
        iv1d this_lipid_nr{};                                                                                                                  
        iv1d this_bind_i{};                                                                                                                   
        iv1d this_bind_f{};                                                                                                                   
        iv1d this_dwell_t{};                                                                                                                  
        iv1d this_res_nr{};                                                                                                                   
        sv1d this_res_name{}; 

        //copy the remaining binding events
        for(i=0; i<lipid_nr.size(); i++) //loop over binding events
        {
            if(removed[i] == 0)
            {
                this_lipid_nr.push_back(lipid_nr[i]);
                this_bind_i.push_back(bind_i[i]);
                this_bind_f.push_back(bind_f[i]);
                this_dwell_t.push_back(dwell_t[i]);
                this_res_nr.push_back(res_nr[i]);
                this_res_name.push_back(res_name[i]);
            }
        }

        //replace old binding events with new ones
        lipid_nr = this_lipid_nr; 
        bind_i   = this_bind_i;
        bind_f   = this_bind_f;
        dwell_t  = this_dwell_t; 
        res_nr   = this_res_nr; 
        res_name = this_res_name; 
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function takes a binding timeline and filters out noise                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Binding_events::sweep_timeline_noise(int lipid_count)
{
    int i = 0;
    int j = 0;

    for(j=0; j<ef_frames; j++) //loop over time line frames
    {
        //store bound lipids (how many grid points was lipid bound)
        iv1d count(0,0);
        for(i=0; i<num_lipids; i++) //loop over time line lipids
        {
            if(bound_time_line[j][i] > 0) //lipid is bound
            {
                count.push_back(bound_time_line[j][i]);
            }
        }

        if(count.size() > 0)
        {
            //organize data in order of largest count
            int b_swap = 0;
            for(b_swap=1; b_swap > 0; )
            {
                b_swap = 0;
                for(i=0; i<count.size()-1; i++)
                {
                    if(count[i] < count[i+1])
                    {
                        int count_tmp = count[i];
                        count[i]      = count[i+1];
                        count[i+1]    = count_tmp;
                        b_swap        = 1;
                    }
                }
            }

            //determine the smallest count allowed
            int biggest_count = 0;
            if(lipid_count > count.size())
            {
                biggest_count = count[count.size()-1];
            }
            else
            {
                biggest_count = count[lipid_count-1];
            }

            //remove lipids from timeline (noise)
            for(i=0; i<num_lipids; i++) //loop over time line lipids
            {
                if(bound_time_line[j][i] >= biggest_count) //lipid is a top lipid
                {
                    bound_time_line[j][i] = 1;
                }
                else
                {
                    bound_time_line[j][i] = 0;
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads binding events files until a complete set of lipid_nr is acquired                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Binding_events::get_complete_set(string in_file_name,iv1d &lipid_nr,iv1d &res_nr,sv1d &res_name)
{
    int count       = 0;          //Keep track of position for adding data
    int string_size = 2000;       //How big of a string can we read?
    char my_string[string_size];  //String to hold read in data entries

    //string lipid_info_file_name = base_file_name_i + "_lipid_info.dat";
    string lipid_info_file_name = chop_and_add_tag(in_file_name,"_lipid_info.dat");

    FILE *lipid_info_file = fopen(lipid_info_file_name.c_str(), "r");
    if(lipid_info_file == NULL) //read binding events files to get a complete set
    {
        int i = 0; 
        int j = 0;
        int k = 0;
        int l = 0;

        printf("Could not find a lipid info file (%s). Will read binding events files until a complete lipid set is acquired. \n",lipid_info_file_name.c_str());

        for(i=0; i<=num_g_x; i++) //loop over x
        {
            printf("Working on column %d. Set contains %d lipids with %d needed for a complete set. \n",i,res_nr.size(),num_lipids);

            for(j=0; j<num_g_y; j++) //loop over y
            {
                Binding_events events;

		//string in_file_name = base_file_name_i + ".be";

                int result = events.get_info(in_file_name);

                if(result == 0)
                {
                    printf("Could not find binding events file (%s). This will most likely crash the program. \n",in_file_name.c_str());
                }
                else
                {
                    result = events.get_binding_events_xy(in_file_name,i,j);

                    if(result == 1 && events.lipid_nr.size() > 0) //binding events file present with data inside
                    {
                        for(k=0; k<events.lipid_nr.size(); k++) //loop over binding events
                        {
                            if(res_nr.size() < num_lipids) //not yet a complete set
                            {
                                int found = 0;

                                for(l=0; l<res_nr.size(); l++) //loop over residue names and numbers in set
                                {
                                    if(events.res_nr[k] == res_nr[l]) //lipid already in set
                                    {
                                        found = 1;
                                        break;
                                    }
                                }

                                if(found == 0) //add lipid to set
                                {
                                    res_nr.push_back(events.res_nr[k]);
                                    res_name.push_back(events.res_name[k]);
                                    lipid_nr.push_back(events.lipid_nr[k]);
                                }
                            }
                            else
                            {
                                goto end_loop;
                            }
                        }
                    }
                }
            }
        }
        end_loop:;

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Organize by lipid number                                                                                  //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        int b_swap = 0;          //Used to tell if aray requires further sorting

        for(b_swap=1; b_swap > 0; )
        {
            b_swap = 0;
            for(i=0; i<num_lipids-1; i++)
            {
                if(lipid_nr[i] > lipid_nr[i+1])
                {
                    int lipid_nr_tmp    = lipid_nr[i];
                    int res_nr_tmp      = res_nr[i];
                    string res_name_tmp = res_name[i];

                    lipid_nr[i] = lipid_nr[i+1];
                    res_nr[i]   = res_nr[i+1];
                    res_name[i] = res_name[i+1];

                    lipid_nr[i+1] = lipid_nr_tmp;
                    res_nr[i+1]   = res_nr_tmp;
                     res_name[i+1] = res_name_tmp;

                    b_swap = 1;
                }
            }
        }
        printf("A complete lipid set has been acquired. \n\n");
    }
    else //read in the data from the lipid_info.dat file
    {
        while(fscanf(lipid_info_file, "%s,", my_string) != EOF)
        {
           if(my_string[0] != '#')
           {
               if(count%3 == 0) //lipid_nr
               {
                   int lipid_number = atoi(my_string);
                   lipid_nr.push_back(lipid_number);
               }
               else if(count%3 == 1) //res_name
               {
                   string this_res_name(my_string);
                   res_name.push_back(this_res_name);
               }
               else if(count%3 == 2) //res_id
               {
                   int res_id = atoi(my_string);
                   res_nr.push_back(res_id);
               }
               count++;
           }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function analyzes lipid tessellation data and returns the center of a specified lipid                //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
dv1d Binding_events::find_voro_center(int this_lipid)
{
    int j     = 0;
    int k     = 0;
    int count = 0;
    dv1d center(2,0.0);

    for(j=0; j<num_g_y; j++) //loop over y-dimension
    {
        for(k=0; k<num_g_x; k++) //loop over x-dimension
        {
            if(voro_nan_frame[j][k] == 0) //check that grid point is not the protein
            {
                if(voro_frame[j][k] == this_lipid)
                {
                    center[0] = center[0] + (double)k;
                    center[1] = center[1] + (double)j;
                    count     = count + 1;
                }
            }
        }
    }

    if(count > 0)
    {
        center[0] = center[0]/(double)count;
        center[1] = center[1]/(double)count;
    }
    else
    {
        center[0] = -999999;
        center[1] = -999999;
    }
    
    return center;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads a binding events file for a particular lattice point                                  //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Binding_events::get_binding_events_grid(string in_file_name,i64v2d &this_be_pos,int this_x,int this_y)
{
    int i      = 0;  //standard variable used in loops
    int result = 0;  //tells if the file was read successfully or not

    int resize = 1;

    if(resize == 1)
    {
        lipid_nr.resize(0,0);
        res_nr.resize(0,0);
        res_name.resize(0);
        bind_i.resize(0,0);
        bind_f.resize(0,0);
        dwell_t.resize(0,0);
    }

    x_i         = 0;
    y_i         = 0;
    ef_frames   = 0;
    num_lipids  = 0;
    num_g_x     = 0;
    num_g_y     = 0;
    ef_dt       = 0.0;
    APS         = 0;

    ifstream be_file(in_file_name, ios::out | ios::binary);
    if(!be_file)
    {

    }
    else
    {
        //read header info
        be_file.read(reinterpret_cast<char *>(&ef_dt),       sizeof(double)); //ef_dt
        be_file.read(reinterpret_cast<char *>(&ef_frames),   sizeof(int));    //ef_frames
        be_file.read(reinterpret_cast<char *>(&num_lipids),  sizeof(int));    //num_lipids
        be_file.read(reinterpret_cast<char *>(&num_g_x),     sizeof(int));    //number of grid points in x
        be_file.read(reinterpret_cast<char *>(&num_g_y),     sizeof(int));    //number of grid points in y
        be_file.read(reinterpret_cast<char *>(&APS),         sizeof(double)); //APS

        //set the file position to the desired lattice point
        be_file.seekg(this_be_pos[this_y][this_x]);

        int num_events = 0;

        be_file.read(reinterpret_cast<char *>(&x_i),            sizeof(int));    //grid points index in x
        be_file.read(reinterpret_cast<char *>(&y_i),            sizeof(int));    //grid points index in y
        be_file.read(reinterpret_cast<char *>(&num_events),     sizeof(int));    //number of binding events

        for(i=0; i<num_events; i++) //loop over binding events
        {
            int this_lipid_nr       = 0;
            int this_res_nr         = 0;
            int this_bi             = 0;
            int this_bf             = 0;
            int this_time           = 0;
            string this_res_name;
            size_t this_res_name_size; 

            be_file.read(reinterpret_cast<char *>(&this_lipid_nr),      sizeof(int));                 //lipid number
            be_file.read(reinterpret_cast<char *>(&this_res_nr),        sizeof(int));                 //residue id
            be_file.read(reinterpret_cast<char *>(&this_res_name_size), sizeof(this_res_name_size));  //size of residue name

            this_res_name.resize(this_res_name_size);

            be_file.read(reinterpret_cast<char *>(&this_res_name[0]),   this_res_name_size);          //lipid number
	    be_file.read(reinterpret_cast<char *>(&this_bi),            sizeof(int));                 //bi
            be_file.read(reinterpret_cast<char *>(&this_bf),            sizeof(int));                 //bf
            be_file.read(reinterpret_cast<char *>(&this_time),          sizeof(int));                 //time

            lipid_nr.push_back(this_lipid_nr);
            res_nr.push_back(this_res_nr);
            res_name.push_back(this_res_name);
            bind_i.push_back(this_bi);
            bind_f.push_back(this_bf);
            dwell_t.push_back(this_time);
        }
        be_file.close();

	result = 1;
    }

    //lipid mixing stuff
    lip_nr_1  = x_i;
    res_nr_1  = y_i;
    num_lip_1 = num_g_x;
    num_lip_2 = num_lipids;

    return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads a binding events file for a grid point and adds them to existing binding events list  //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Binding_events::add_binding_events_grid(string in_file_name,i64v2d &this_be_pos,int this_x,int this_y)
{
    int i      = 0;  //standard variable used in loops
    int result = 0;  //tells if the file was read successfully or not

    ifstream be_file(in_file_name, ios::out | ios::binary);
    if(!be_file)
    {

    }
    else
    {
        //set the file position to the desired lattice point
        be_file.seekg(this_be_pos[this_y][this_x]);

        int num_events = 0;

        be_file.read(reinterpret_cast<char *>(&x_i),            sizeof(int));    //grid points index in x
        be_file.read(reinterpret_cast<char *>(&y_i),            sizeof(int));    //grid points index in y
        be_file.read(reinterpret_cast<char *>(&num_events),     sizeof(int));    //number of binding events

        for(i=0; i<num_events; i++) //loop over binding events
        {
            int this_lipid_nr       = 0;
            int this_res_nr         = 0;
            int this_bi             = 0;
            int this_bf             = 0;
            int this_time           = 0;
            string this_res_name;
            size_t this_res_name_size;

            be_file.read(reinterpret_cast<char *>(&this_lipid_nr),      sizeof(int));                 //lipid number
            be_file.read(reinterpret_cast<char *>(&this_res_nr),        sizeof(int));                 //residue id
            be_file.read(reinterpret_cast<char *>(&this_res_name_size), sizeof(this_res_name_size));  //size of residue name
            
            this_res_name.resize(this_res_name_size);

            be_file.read(reinterpret_cast<char *>(&this_res_name[0]),   this_res_name_size);          //lipid number
            be_file.read(reinterpret_cast<char *>(&this_bi),            sizeof(int));                 //bi
            be_file.read(reinterpret_cast<char *>(&this_bf),            sizeof(int));                 //bf
            be_file.read(reinterpret_cast<char *>(&this_time),          sizeof(int));                 //time

            lipid_nr.push_back(this_lipid_nr);
            res_nr.push_back(this_res_nr);
            res_name.push_back(this_res_name);
            bind_i.push_back(this_bi);
            bind_f.push_back(this_bf);
            dwell_t.push_back(this_time);
        }
        be_file.close();

        result = 1;
    }

    //lipid mixing stuff
    lip_nr_1  = x_i;
    res_nr_1  = y_i;
    num_lip_1 = num_g_x;
    num_lip_2 = num_lipids;

    return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes the binding events to file in binary                                                 //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Binding_events::write_binding_events_bin(string binding_events_file_name)
{
    int i          = 0;
    int num_events = lipid_nr.size();

    ofstream be_file(binding_events_file_name, ios::out | ios::binary);

    //write header info
    be_file.write(reinterpret_cast<const char *>(&ef_dt),      sizeof(double)); //ef_dt
    be_file.write(reinterpret_cast<const char *>(&ef_frames),  sizeof(int));    //ef_frames
    be_file.write(reinterpret_cast<const char *>(&num_lipids), sizeof(int));    //num_lipids
    be_file.write(reinterpret_cast<const char *>(&num_g_x),    sizeof(int));    //number of grid points in x
    be_file.write(reinterpret_cast<const char *>(&num_g_y),    sizeof(int));    //number of grid points in y
    be_file.write(reinterpret_cast<const char *>(&APS),        sizeof(double)); //APS
    be_file.write(reinterpret_cast<const char *>(&x_i),        sizeof(int));    //grid points index in x
    be_file.write(reinterpret_cast<const char *>(&y_i),        sizeof(int));    //grid points index in y
    be_file.write(reinterpret_cast<const char *>(&num_events), sizeof(int));    //number of binding events

    for(i=0; i<num_events; i++) //loop over binding events
    {
        size_t this_res_name_size = res_name[i].size();

        be_file.write(reinterpret_cast<const char *>(&lipid_nr[i]),         sizeof(int));                //lipid number
        be_file.write(reinterpret_cast<const char *>(&res_nr[i]),           sizeof(int));                //residue id
        be_file.write(reinterpret_cast<const char *>(&this_res_name_size),  sizeof(this_res_name_size)); //size of residue name
        be_file.write(reinterpret_cast<const char *>(res_name[i].c_str()),  this_res_name_size);         //lipid number
        be_file.write(reinterpret_cast<const char *>(&bind_i[i]),           sizeof(int));                //bi
        be_file.write(reinterpret_cast<const char *>(&bind_f[i]),           sizeof(int));                //bf
        be_file.write(reinterpret_cast<const char *>(&dwell_t[i]),          sizeof(int));                //time
    }
    be_file.close();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads in a binary binding events file                                                       //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Binding_events::get_binding_events_bin(string binding_events_file_name)
{
    int i          = 0;  //standard variable used in loops
    int num_events = 0;  //how many binding events
    int result     = 0;  //tells if the file was read successfully or not

    int resize = 1;

    if(resize == 1)
    {
        lipid_nr.resize(0,0);
        res_nr.resize(0,0);
        res_name.resize(0);
        bind_i.resize(0,0);
        bind_f.resize(0,0);
        dwell_t.resize(0,0);
    }

    x_i         = 0;
    y_i         = 0;
    ef_frames   = 0;
    num_lipids  = 0;
    num_g_x     = 0;
    num_g_y     = 0;
    ef_dt       = 0.0;
    APS         = 0;

    ifstream be_file(binding_events_file_name, ios::out | ios::binary);
    if(!be_file)
    {

    }
    else
    {
        //read header info
        be_file.read(reinterpret_cast<char *>(&ef_dt),      sizeof(double)); //ef_dt
        be_file.read(reinterpret_cast<char *>(&ef_frames),  sizeof(int));    //ef_frames
        be_file.read(reinterpret_cast<char *>(&num_lipids), sizeof(int));    //num_lipids
        be_file.read(reinterpret_cast<char *>(&num_g_x),    sizeof(int));    //number of grid points in x
        be_file.read(reinterpret_cast<char *>(&num_g_y),    sizeof(int));    //number of grid points in y
        be_file.read(reinterpret_cast<char *>(&APS),        sizeof(double)); //APS
        be_file.read(reinterpret_cast<char *>(&x_i),        sizeof(int));    //grid points index in x
        be_file.read(reinterpret_cast<char *>(&y_i),        sizeof(int));    //grid points index in y
        be_file.read(reinterpret_cast<char *>(&num_events), sizeof(int));    //number of binding events

        for(i=0; i<num_events; i++) //loop over binding events
        {
            int this_lipid_nr       = 0;
            int this_res_nr         = 0;
            int this_bi             = 0;
            int this_bf             = 0;
            int this_time           = 0;
            string this_res_name;
            size_t this_res_name_size;

            be_file.read(reinterpret_cast<char *>(&this_lipid_nr),      sizeof(int));                //lipid number
            be_file.read(reinterpret_cast<char *>(&this_res_nr),        sizeof(int));                //residue id
            be_file.read(reinterpret_cast<char *>(&this_res_name_size), sizeof(this_res_name_size)); //size of residue name

            this_res_name.resize(this_res_name_size);

            be_file.read(reinterpret_cast<char *>(&this_res_name[0]),   this_res_name_size);         //lipid number
            be_file.read(reinterpret_cast<char *>(&this_bi),            sizeof(int));                //bi
            be_file.read(reinterpret_cast<char *>(&this_bf),            sizeof(int));                //bf
            be_file.read(reinterpret_cast<char *>(&this_time),          sizeof(int));                //time

            lipid_nr.push_back(this_lipid_nr);
            res_nr.push_back(this_res_nr);
            res_name.push_back(this_res_name);
            bind_i.push_back(this_bi);
            bind_f.push_back(this_bf);
            dwell_t.push_back(this_time);
	}
        be_file.close();

	result = 1;
    }

    return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes the temporary binding events for lattice points for the MPI core                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int64_t Binding_events::write_binding_events_tmp(ofstream &be_file_o,int64_t pos)
{
    int i          = 0;
    int num_events = lipid_nr.size();
    int64_t new_pos;

    be_file_o.seekp(pos);

    //write header info
    be_file_o.write(reinterpret_cast<const char *>(&ef_dt),      sizeof(double)); //ef_dt
    be_file_o.write(reinterpret_cast<const char *>(&ef_frames),  sizeof(int));    //ef_frames
    be_file_o.write(reinterpret_cast<const char *>(&num_lipids), sizeof(int));    //num_lipids
    be_file_o.write(reinterpret_cast<const char *>(&num_g_x),    sizeof(int));    //number of grid points in x
    be_file_o.write(reinterpret_cast<const char *>(&num_g_y),    sizeof(int));    //number of grid points in y
    be_file_o.write(reinterpret_cast<const char *>(&APS),        sizeof(double)); //APS
    be_file_o.write(reinterpret_cast<const char *>(&x_i),        sizeof(int));    //grid points index in x
    be_file_o.write(reinterpret_cast<const char *>(&y_i),        sizeof(int));    //grid points index in y
    be_file_o.write(reinterpret_cast<const char *>(&num_events), sizeof(int));    //number of binding events

    for(i=0; i<num_events; i++) //loop over binding events
    {
        size_t this_res_name_size = res_name[i].size();

        be_file_o.write(reinterpret_cast<const char *>(&lipid_nr[i]),         sizeof(int));                //lipid number
        be_file_o.write(reinterpret_cast<const char *>(&res_nr[i]),           sizeof(int));                //residue id
        be_file_o.write(reinterpret_cast<const char *>(&this_res_name_size),  sizeof(this_res_name_size)); //size of residue name
        be_file_o.write(reinterpret_cast<const char *>(res_name[i].c_str()),  this_res_name_size);         //lipid number
        be_file_o.write(reinterpret_cast<const char *>(&bind_i[i]),           sizeof(int));                //bi
        be_file_o.write(reinterpret_cast<const char *>(&bind_f[i]),           sizeof(int));                //bf
        be_file_o.write(reinterpret_cast<const char *>(&dwell_t[i]),          sizeof(int));                //time
    }
    new_pos = be_file_o.tellp();

    return new_pos; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads in binding events from the temporary be file                                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int64_t Binding_events::get_binding_events_tmp(ifstream &be_file_i,int64_t pos)
{
    int i          = 0;  //standard variable used in loops
    int num_events = 0;  //how many binding events
    int64_t new_pos;

    int resize = 1;

    if(resize == 1)
    {
        lipid_nr.resize(0,0);
        res_nr.resize(0,0);
        res_name.resize(0);
        bind_i.resize(0,0);
        bind_f.resize(0,0);
        dwell_t.resize(0,0);
    }

    x_i         = 0;
    y_i         = 0;
    ef_frames   = 0;
    num_lipids  = 0;
    num_g_x     = 0;
    num_g_y     = 0;
    ef_dt       = 0.0;
    APS         = 0;

    be_file_i.seekg(pos);
 
    //read header info
    be_file_i.read(reinterpret_cast<char *>(&ef_dt),      sizeof(double)); //ef_dt
    be_file_i.read(reinterpret_cast<char *>(&ef_frames),  sizeof(int));    //ef_frames
    be_file_i.read(reinterpret_cast<char *>(&num_lipids), sizeof(int));    //num_lipids
    be_file_i.read(reinterpret_cast<char *>(&num_g_x),    sizeof(int));    //number of grid points in x
    be_file_i.read(reinterpret_cast<char *>(&num_g_y),    sizeof(int));    //number of grid points in y
    be_file_i.read(reinterpret_cast<char *>(&APS),        sizeof(double)); //APS
    be_file_i.read(reinterpret_cast<char *>(&x_i),        sizeof(int));    //grid points index in x
    be_file_i.read(reinterpret_cast<char *>(&y_i),        sizeof(int));    //grid points index in y
    be_file_i.read(reinterpret_cast<char *>(&num_events), sizeof(int));    //number of binding events

    for(i=0; i<num_events; i++) //loop over binding events
    {
        int this_lipid_nr       = 0;
        int this_res_nr         = 0;
        int this_bi             = 0;
        int this_bf             = 0;
        int this_time           = 0;
        string this_res_name;
        size_t this_res_name_size;

        be_file_i.read(reinterpret_cast<char *>(&this_lipid_nr),      sizeof(int));                //lipid number
	be_file_i.read(reinterpret_cast<char *>(&this_res_nr),        sizeof(int));                //residue id
	be_file_i.read(reinterpret_cast<char *>(&this_res_name_size), sizeof(this_res_name_size)); //size of residue name

        this_res_name.resize(this_res_name_size);

        be_file_i.read(reinterpret_cast<char *>(&this_res_name[0]),   this_res_name_size);         //lipid number
	be_file_i.read(reinterpret_cast<char *>(&this_bi),            sizeof(int));                //bi
	be_file_i.read(reinterpret_cast<char *>(&this_bf),            sizeof(int));                //bf
	be_file_i.read(reinterpret_cast<char *>(&this_time),          sizeof(int));                //time

        lipid_nr.push_back(this_lipid_nr);
        res_nr.push_back(this_res_nr);
        res_name.push_back(this_res_name);
        bind_i.push_back(this_bi);
        bind_f.push_back(this_bf);
        dwell_t.push_back(this_time);
    }

    new_pos = be_file_i.tellg();

    return new_pos;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads in the info file for the binding events file                                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Binding_events::get_info(string binding_events_file_name)
{
    int i           = 0;           //standard variable used in loops
    int j           = 0;           //standard variable used in loops
    int result      = 0;           //was the BE file opened successfully
    int num_events  = 0;           //how many binding events
    int string_size = 2000;        //size of string used for reading in data
    char my_string[string_size];   //used to read in data

    ifstream be_file(binding_events_file_name, ios::out | ios::binary);
    if(!be_file)
    {

    }
    else
    {
        //read header info
        be_file.read(reinterpret_cast<char *>(&ef_dt),      sizeof(double)); //ef_dt
        be_file.read(reinterpret_cast<char *>(&ef_frames),  sizeof(int));    //ef_frames
        be_file.read(reinterpret_cast<char *>(&num_lipids), sizeof(int));    //num_lipids
        be_file.read(reinterpret_cast<char *>(&num_g_x),    sizeof(int));    //number of grid points in x
        be_file.read(reinterpret_cast<char *>(&num_g_y),    sizeof(int));    //number of grid points in y
        be_file.read(reinterpret_cast<char *>(&APS),        sizeof(double)); //APS
        be_file.read(reinterpret_cast<char *>(&x_i),        sizeof(int));    //grid points index in x
        be_file.read(reinterpret_cast<char *>(&y_i),        sizeof(int));    //grid points index in y
        be_file.read(reinterpret_cast<char *>(&num_events), sizeof(int));    //number of binding events

        //allocate memory to hold file positions
        info_pos.resize(num_g_x);
        for(i=0; i<num_g_x; i++)
        {
            info_pos[i].resize(num_g_y);
        }

        //read in .info file for the binding events file
        string info_file_name = binding_events_file_name + ".info";
        FILE *this_file = fopen(info_file_name.c_str(), "r");
        if(this_file == NULL)
        {
            printf("Could not find an info file (%s). Will analyze the binding events file. \n",info_file_name.c_str());

            int64_t current_pos = 0;
            be_file.seekg(current_pos);

            for(i=0; i<num_g_x; i++) //loop over x
            {
                for(j=0; j<num_g_y; j++) //loop over y 
                {
                    info_pos[i][j] = current_pos;
                    current_pos = get_binding_events_tmp(be_file,current_pos);
                }
            }

            //write new info file
            int world_size = 0; 
            int world_rank = -1;
            MPI_Comm_size(MPI_COMM_WORLD, &world_size);      //get the world size
            MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);      //get the process rank

            if(world_rank == 0)
            {
                FILE *this_file = fopen(info_file_name.c_str(), "w");
                if(this_file == NULL)
                {
                    printf("Could not open file %s. \n",info_file_name.c_str());
                }
                else
                {
                    for(i=0; i<num_g_x; i++) //loop over x
                    {
                        for(j=0; j<num_g_y; j++) //loop over y 
                        {
                            fprintf(this_file," %ld ",info_pos[i][j]);
                        }
                        fprintf(this_file,"\n");
                    }
                    fclose(this_file);
                }
            }
        }
        else //read info file
        {
            for(i=0; i<num_g_x; i++) //loop over x
            {
                for(j=0; j<num_g_y; j++) //loop over y 
                {
                    fscanf(this_file, "%s,", my_string);
                    info_pos[i][j] = atol(my_string);
                }
                fprintf(this_file,"\n");
            }
            fclose(this_file);
        }

	//close binding events file
        be_file.close();

	result = 1;
    }
    return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads in a set of binding events for a lattice point                                        //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Binding_events::get_binding_events_xy(string binding_events_file_name,int x,int y)
{
    int i          = 0;  //standard variable used in loops
    int num_events = 0;  //how many binding events
    int result     = 0;  //was the BE file opened successfully

    int resize = 1;

    if(resize == 1)
    {
        lipid_nr.resize(0,0);
        res_nr.resize(0,0);
        res_name.resize(0);
        bind_i.resize(0,0);
        bind_f.resize(0,0);
        dwell_t.resize(0,0);
    }

    x_i         = 0;
    y_i         = 0;
    ef_frames   = 0;
    num_lipids  = 0;
    num_g_x     = 0;
    num_g_y     = 0;
    ef_dt       = 0.0;
    APS         = 0;

    ifstream be_file(binding_events_file_name, ios::out | ios::binary);
    if(!be_file)
    {

    }
    else
    {
        be_file.seekg(info_pos[x][y]);

        //read header info
        be_file.read(reinterpret_cast<char *>(&ef_dt),      sizeof(double)); //ef_dt
        be_file.read(reinterpret_cast<char *>(&ef_frames),  sizeof(int));    //ef_frames
        be_file.read(reinterpret_cast<char *>(&num_lipids), sizeof(int));    //num_lipids
        be_file.read(reinterpret_cast<char *>(&num_g_x),    sizeof(int));    //number of grid points in x
        be_file.read(reinterpret_cast<char *>(&num_g_y),    sizeof(int));    //number of grid points in y
        be_file.read(reinterpret_cast<char *>(&APS),        sizeof(double)); //APS
        be_file.read(reinterpret_cast<char *>(&x_i),        sizeof(int));    //grid points index in x
        be_file.read(reinterpret_cast<char *>(&y_i),        sizeof(int));    //grid points index in y
        be_file.read(reinterpret_cast<char *>(&num_events), sizeof(int));    //number of binding events

        for(i=0; i<num_events; i++) //loop over binding events
        {
            int this_lipid_nr       = 0;
            int this_res_nr         = 0;
            int this_bi             = 0;
            int this_bf             = 0;
            int this_time           = 0;
            string this_res_name;
            size_t this_res_name_size;

            be_file.read(reinterpret_cast<char *>(&this_lipid_nr),      sizeof(int));                //lipid number
            be_file.read(reinterpret_cast<char *>(&this_res_nr),        sizeof(int));                //residue id
            be_file.read(reinterpret_cast<char *>(&this_res_name_size), sizeof(this_res_name_size)); //size of residue name

            this_res_name.resize(this_res_name_size);

            be_file.read(reinterpret_cast<char *>(&this_res_name[0]),   this_res_name_size);         //lipid number
            be_file.read(reinterpret_cast<char *>(&this_bi),            sizeof(int));                //bi
            be_file.read(reinterpret_cast<char *>(&this_bf),            sizeof(int));                //bf
            be_file.read(reinterpret_cast<char *>(&this_time),          sizeof(int));                //time

            lipid_nr.push_back(this_lipid_nr);
            res_nr.push_back(this_res_nr);
            res_name.push_back(this_res_name);
            bind_i.push_back(this_bi);
            bind_f.push_back(this_bf);
            dwell_t.push_back(this_time);
        }
	be_file.close();

	result = 1;
    }

    return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads binding events files for the grid and builds the lipid tessellations for each frame   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Binding_events::get_tessellations(string in_file_name,int my_gi,int my_gf,int my_num_g,int stride,int my_ef_frames,int world_rank)
{
    int i          = 0;    //standard variable used in loops
    int j          = 0;    //standard variable used in loops 
    int k          = 0;    //standard variable used in loops
    int l          = 0;    //standard variable used in loops
    int grid_count = 0;    //How many lattice points have been looped over
    int counter    = 0;    //How many times the "program run time" been displayed
    clock_t t;             //Keeps the time for testing performance

    //take the initial time
    t = clock();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Report estimated memory requirement                                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    long   num_entries  = (long)my_ef_frames*(long)my_num_g;
    double mem_voro     = (double)num_entries*4.0/1000000.0;
    double mem_voro_nan = (double)num_entries*4.0/1000000.0;

    if(world_rank == 0)
    {
        printf("Collecting lipid tessellation data. Estimated memory required: %f MB. \n",mem_voro + mem_voro_nan);
        printf("-----------------------------------------------------------------------------------------------------------------------------------\n");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Allocate memory for tessellation data                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    voro.resize(my_ef_frames);
    voro_nan.resize(my_ef_frames);
    for(i=0; i<my_ef_frames; i++)
    {
        voro[i].resize(my_num_g,0);
        voro_nan[i].resize(my_num_g,1);
    }

    //read info file for binding events data 
    int result = get_info(in_file_name);

    if(result == 1)
    {
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Read binding events files and generate single frame tessellation data                                     //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        for(i=0; i<num_g_x; i++) //loop over x
        {
            for(j=0; j<num_g_y; j++) //loop over y
            {
                if(grid_count >= my_gi && grid_count <= my_gf)
                {
                    int ef_g = grid_count - my_gi;

                    result = get_binding_events_xy(in_file_name,i,j);

                    if(lipid_nr.size() > 0)
                    {
                        get_binding_timeline();
                        int ef_frame = 0;

                        for(k=0; k<ef_frames; k+=stride) //loop over frames
                        {
                            for(l=0; l<num_lipids; l++) //loop over lipids
                            {
                                if(bound_time_line[k][l] == 1)
                                {
                                    voro[ef_frame][ef_g]     = time_line_lipid_nr[l];
                                    voro_nan[ef_frame][ef_g] = 0;
                                    break;
                                }
                            }
                            ef_frame = ef_frame + 1;
                        }
                    }
                    //report progress and estimated time to completion
                    int current_step = ef_g + 1;
                    int my_steps     = my_gf - my_gi + 1;
                    be_time_stats(t,&counter,current_step,my_steps,world_rank,"lattice point");
                }
                grid_count++;
            }
        }
    }
    else
    {
        printf("failure opening %s. Make sure the file exists. \n",in_file_name.c_str());
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }

    return (clock() - t)/CLOCKS_PER_SEC;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects the lipid tessellation data for a single frame                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Binding_events::get_voro_frame(int this_frame,int world_size,int world_rank)
{
    int i          = 0;
    int j          = 0;
    int grid_count = 0;

    MPI_Barrier(MPI_COMM_WORLD);

    voro_frame.resize(num_g_y);
    voro_nan_frame.resize(num_g_y);
    for(i=0; i<num_g_y; i++)
    {
        voro_frame[i].resize(num_g_x,0);
        voro_nan_frame[i].resize(num_g_x,1);
    } 

    iv1d current_row_grid = collect_and_clone_iv1d(world_size,world_rank,voro[this_frame]);
    iv1d current_row_nan  = collect_and_clone_iv1d(world_size,world_rank,voro_nan[this_frame]);

    //put data back in the grid
    if(world_rank == 0)
    {	     
        for(i=0; i<num_g_x; i++)
        {
            for(j=0; j<num_g_y; j++)
            {
                voro_frame[j][i]     = current_row_grid[grid_count];
                voro_nan_frame[j][i] = current_row_nan[grid_count];
                grid_count++;
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes the current frame tessellation data to file                                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Binding_events::write_voro_frame(string out_file_name)
{
    int j = 0;
    int k = 0;

    FILE *out_file = fopen(out_file_name.c_str(),"w");

    if(out_file == NULL)
    {
        printf("failure writing %s. Make sure there is enough disc space available and that the target directory exists. \n",out_file_name.c_str());
        fflush(stdin);
    }
    else
    {
        for(j=0; j<num_g_y; j++) //loop over y-dimension
        {
            for(k=0; k<num_g_x; k++) //loop over x-dimension
            {
                if(voro_nan_frame[j][k] == 0)
                {
                    fprintf(out_file," %10d",voro_frame[j][k]);
                }
                else //data excluded
                {
                    fprintf(out_file," %10s ","NaN");
                }
            }
            fprintf(out_file,"\n");
        }
        fclose(out_file);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads in a binary binding events file and adds events to an existing set of data            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Binding_events::add_binding_events(string binding_events_file_name,int compound_lipid_count,int lipid_nr_offset)
{
    int i          = 0;  //standard variable used in loops
    int num_events = 0;  //how many binding events
    int result     = 0;  //tells if the file was read successfully or not

    int num_lipids_store = num_lipids;

    ifstream be_file(binding_events_file_name, ios::out | ios::binary);
    if(!be_file)
    {

    }
    else
    {
        //read header info
        be_file.read(reinterpret_cast<char *>(&ef_dt),      sizeof(double)); //ef_dt
        be_file.read(reinterpret_cast<char *>(&ef_frames),  sizeof(int));    //ef_frames
        be_file.read(reinterpret_cast<char *>(&num_lipids), sizeof(int));    //num_lipids
        be_file.read(reinterpret_cast<char *>(&num_g_x),    sizeof(int));    //number of grid points in x
        be_file.read(reinterpret_cast<char *>(&num_g_y),    sizeof(int));    //number of grid points in y
        be_file.read(reinterpret_cast<char *>(&APS),        sizeof(double)); //APS
        be_file.read(reinterpret_cast<char *>(&x_i),        sizeof(int));    //grid points index in x
        be_file.read(reinterpret_cast<char *>(&y_i),        sizeof(int));    //grid points index in y
        be_file.read(reinterpret_cast<char *>(&num_events), sizeof(int));    //number of binding events

        for(i=0; i<num_events; i++) //loop over binding events
        {
            int this_lipid_nr       = 0;
            int this_res_nr         = 0;
            int this_bi             = 0;
            int this_bf             = 0;
            int this_time           = 0;
            string this_res_name;
            size_t this_res_name_size;

            be_file.read(reinterpret_cast<char *>(&this_lipid_nr),      sizeof(int));                //lipid number
            be_file.read(reinterpret_cast<char *>(&this_res_nr),        sizeof(int));                //residue id
            be_file.read(reinterpret_cast<char *>(&this_res_name_size), sizeof(this_res_name_size)); //size of residue name

            this_res_name.resize(this_res_name_size);

            be_file.read(reinterpret_cast<char *>(&this_res_name[0]),   this_res_name_size);         //lipid number
            be_file.read(reinterpret_cast<char *>(&this_bi),            sizeof(int));                //bi
            be_file.read(reinterpret_cast<char *>(&this_bf),            sizeof(int));                //bf
            be_file.read(reinterpret_cast<char *>(&this_time),          sizeof(int));                //time

            lipid_nr.push_back(this_lipid_nr + lipid_nr_offset);
            res_nr.push_back(this_res_nr);
            res_name.push_back(this_res_name);
            bind_i.push_back(this_bi);
            bind_f.push_back(this_bf);
            dwell_t.push_back(this_time);
        }
        be_file.close();

        result = 1;
    }

    if(compound_lipid_count == 1)
    {
        num_lipids = num_lipids + num_lipids_store;
    }

    return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes the binding events to file                                                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Binding_events::write_binding_events_legacy(string binding_events_file_name)
{
    int i = 0;

    FILE *binding_events_file = fopen(binding_events_file_name.c_str(), "w");
    if(binding_events_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",binding_events_file_name.c_str());
    }
    else
    {
        fprintf(binding_events_file," x_i %10d y_i %10d ef_dt(ps) %10f ef_frames %d num_lipids %10d num_g_x %10d num_g_y %10d APS(nm^2) %10f \n\n",x_i,y_i,ef_dt,ef_frames,num_lipids,num_g_x,num_g_y,APS);
        fprintf(binding_events_file," %10s %10s %10s %15s %15s %20s \n","lipid","res_nr","res_name","bind_i(frame)","bind_f(frame)","dwell time(frames)");
        fprintf(binding_events_file," %10s-%10s-%10s-%15s-%15s-%20s \n","----------","----------","----------","---------------","---------------","--------------------");

        for(i=0; i<dwell_t.size(); i++)
        {
            fprintf(binding_events_file," %10d %10d %10s %15d %15d %20d \n",lipid_nr[i],res_nr[i],res_name[i].c_str(),bind_i[i],bind_f[i],dwell_t[i]);
        }
        fclose(binding_events_file);
    }
}

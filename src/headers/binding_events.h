
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads a binding events file                                                                 //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int read_binding_events(string in_file_name,iv1d &lipid_nr,iv1d &res_nr,sv1d &res_name,iv1d &bind_i,iv1d &bind_f,
                        iv1d &dwell_t,int *x_i,int *y_i,double *ef_dt,int *ef_frames,int *num_lipids,int *num_g_x,
                        int *num_g_y,double *APS,int resize,int lipid_nr_offset)
{
    int number_of_lines = 0;      //Number of lines in input files
    int items_per_line  = 0;      //How many item is a single line
    int k               = 0;      //Standard variable used in loops
    int l               = 0;      //Standard variable used in loops
    int m               = 0;      //Standard variable used in loops
    int outcome         = 0;      //Return whether the file was read succesfully or not
    char my_string[200];          //String to hold read in data entries

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Resize binding events vectors                                                                             //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(resize == 1)
    {
        lipid_nr.resize(0,0);
        res_nr.resize(0,0);
        res_name.resize(0);
        bind_i.resize(0,0);
        bind_f.resize(0,0);
        dwell_t.resize(0,0);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Here we open files for reading/writing                                                                    //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    FILE *in_file = fopen(in_file_name.c_str(), "r");
    if(in_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",in_file_name.c_str());
    }
    else
    {
        outcome = 1;

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Get information about the data files                                                                      //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        characterize_file(in_file_name,&number_of_lines,&items_per_line,4);

        //printf("in_file_name %20s number_of_lines %10d items_per_line %10d \n",in_file_name.c_str(),number_of_lines,items_per_line);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Read in binding events                                                                                    //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        for(k=0; k<number_of_lines+4; k++) //loop over rows
        {
            if(k==0) //first line containts x_i,y_i,ef_dt,ef_frames,num_lipids,num_g_x and num_g_y
            {
                for(l=0; l<16; l++) //loop over header items
                {
                    fscanf(in_file, "%s,", my_string);

                    if(l == 1) //x_i
                    {
                        *x_i = atoi(my_string);
                    }
                    if(l == 3) //y_i
                    {
                        *y_i = atoi(my_string);
                    }
                    if(l == 5) //ef_dt
                    {
                        *ef_dt = atof(my_string);
                    }
                    if(l == 7) //ef_frames
                    {
                        *ef_frames = atoi(my_string);
                    }
                    if(l == 9) //num_lipids
                    {
                        *num_lipids = atoi(my_string);
                    }
                    if(l == 11) //num_g_x
                    {
                        *num_g_x = atoi(my_string);
                    }
                    if(l == 13) //num_g_y
                    {
                        *num_g_y = atoi(my_string);
                    }
                    if(l == 15) //APS
                    {
                        *APS = atof(my_string);
                        next_line(in_file);
                    }
                }
            }
            else if(k>3) //first 4 lines are header information
            {
                for(l=0; l<items_per_line; l++) //loop over collumns
                {
                    fscanf(in_file, "%s,", my_string);

                    if(l == 0) //lipid number
                    {
                        lipid_nr.push_back(atoi(my_string) + lipid_nr_offset);
                    }
                    if(l == 1) //res number
                    {
                        res_nr.push_back(atoi(my_string));
                    }
                    if(l == 2) //res name
                    {
                        res_name.push_back(strdup(my_string));
                    }
                    if(l == 3) //bind_i
                    {
                        bind_i.push_back(atoi(my_string));
                    }
                    if(l == 4) //bind_f
                    {
                        bind_f.push_back(atoi(my_string));
                    }
                    if(l == 5) //dwell_t
                    {
                        dwell_t.push_back(atoi(my_string));
                    }
                }
            }
            else //read the line
            {
                next_line(in_file);
            }
        }
        fclose(in_file);
    }
    return outcome;
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

        //lipid blobs
        iv3d blob;                                                                                                                       //grid containing blobs
        iv3d blob_nan;                                                                                                                   //grid containing blobs nan info
        iv2d blob_frame;                                                                                                                 //current frame blob data
        iv2d blob_nan_frame;                                                                                                             //current frame blob nan data

    public:
        int  get_binding_events(string in_file_name);                                                                                    //reads the binding events file
        int  add_binding_events(string in_file_name,int compound_lipid_count,int lipid_nr_offset);                                       //reads the binding events file and adds to current list
        void gen_lipid_nr();                                                                                                             //generates a unique lipid nr for each lipid
        void organize_events(int mode);                                                                                                  //organize binding events
        void reduce_list();                                                                                                              //characterize repeat visits
        void size_timeline();                                                                                                            //allocate memory for the binding time line
        void get_binding_timeline();                                                                                                     //converts binding events into a timeline
        void add_to_binding_timeline();                                                                                                  //add binding events to time line 
        void stamp_to_binding_timeline();                                                                                                //stamps binding events to existing timeline
        void binding_events_from_timeline();                                                                                             //use binding timeline to create binding events
        void write_binding_events(string binding_events_file_name);                                                                      //write the binding events to file
        void write_time_line(string timeline_file_name);                                                                                 //write the binding timeline to file
        void suppress_timeline_noise(int threshold);                                                                                     //mend fragmented binding events
        void sweep_timeline_noise(int lipid_count);                                                                                      //constrain numbero of lipids bound in timeline for a given t
        void get_complete_set(string base_file_name_i,iv1d &lipid_nr,iv1d &res_nr,sv1d &res_name);                                       //reads binding events file until a complete set of lipid_nr etc. is found
        void get_blobs(string base_file_name_i,int my_xi,int my_xf,int my_num_g_x,int stride,int my_ef_frames,int world_rank);           //reads binding events files and extracts the lipid blobs for each frame
        void get_blobs_frame(int this_frame,int world_size,int world_rank);                                                              //collects a single frame of the lipid blobs
        void write_blobs_frame(string out_file_name);                                                                                    //writes the current frame blob data to file
        dv1d find_blobs_center(int this_lipid);                                                                                          //returns the center of a lipid blob
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads a binding events file                                                                 //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Binding_events::get_binding_events(string in_file_name)
{
    int resize = 1;

    x_i         = 0;                                      
    y_i         = 0;                                     
    ef_frames   = 0;                                     
    num_lipids  = 0;                                     
    num_g_x     = 0;                                     
    num_g_y     = 0;                                     
    ef_dt       = 0.0;                                  
    APS         = 0;                                      

    int result = read_binding_events(in_file_name,lipid_nr,res_nr,res_name,bind_i,bind_f,dwell_t,&x_i,&y_i,&ef_dt,&ef_frames,&num_lipids,&num_g_x,&num_g_y,&APS,resize,0);

    //lipid mixing stuff
    lip_nr_1  = x_i;
    res_nr_1  = y_i;
    num_lip_1 = num_g_x; 
    num_lip_2 = num_lipids;

    return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads a binding events file and add them to existing binding events list                    //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Binding_events::add_binding_events(string in_file_name,int compound_lipid_count,int lipid_nr_offset)
{   
    int num_lipids_store = num_lipids;

    int resize = 0;
    int result = read_binding_events(in_file_name,lipid_nr,res_nr,res_name,bind_i,bind_f,dwell_t,&x_i,&y_i,&ef_dt,&ef_frames,&num_lipids,&num_g_x,&num_g_y,&APS,resize,lipid_nr_offset);

    if(compound_lipid_count == 1)
    { 
        num_lipids = num_lipids + num_lipids_store;    
    }

    return result;
}

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
            else if(bound_time_line[j][i] == 0 && prev_state == 1) //lipid just left, write event to file
            {
                bind_final = j-1;
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
// This function writes the binding events to file                                                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Binding_events::write_binding_events(string binding_events_file_name)
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
void Binding_events::get_complete_set(string base_file_name_i,iv1d &lipid_nr,iv1d &res_nr,sv1d &res_name)
{
    int i = 0; 
    int j = 0;
    int k = 0;
    int l = 0;

    printf("Reading binding events until a complete lipid set is acquired. \n");

    for(i=0; i<=num_g_x; i++) //loop over x
    {
        printf("Working on column %d. Set contains %d lipids with %d needed for a complete set. \n",i,res_nr.size(),num_lipids);

        for(j=0; j<num_g_y; j++) //loop over y
        {
            Binding_events events;
            string in_file_name = base_file_name_i + "_" + to_string(i) + "_" + to_string(j) + ".be";
            int result          = events.get_binding_events(in_file_name);

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

    //for(i=0; i<num_lipids; i++)
    //{
    //    printf("lipid_nr %d res_nr %d res_name %s \n",lipid_nr[i],res_nr[i],res_name[i].c_str());
    //}

    printf("A complete lipid set has been acquired. \n\n");
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads binding events files for the grid and builds the lipid blobs for each frame           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Binding_events::get_blobs(string base_file_name_i,int my_xi,int my_xf,int my_num_g_x,int stride,int my_ef_frames,int world_rank)
{
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Back up header info in case a gridpoint (last one examined in get_blobs) has no be file this info could   //
    // be lost                                                                                                   //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int ef_frames_store  = ef_frames;
    int ef_dt_store      = ef_dt;
    int num_lipids_store = num_lipids;
    int num_g_x_store    = num_g_x;
    int num_g_y_store    = num_g_y;
    int APS_store        = APS;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Report estimated memory requirement                                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int num_entries     = my_ef_frames*num_g_y*my_num_g_x;
    double mem_blob     = num_entries*4.0/1000000.0;
    double mem_blob_nan = num_entries*4.0/1000000.0;

    if(world_rank == 0)
    {
        printf("Collecting lipid blob data. Estimated memory required: %f MB. \n",mem_blob + mem_blob_nan);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Allocate memory for blob data                                                                             //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    blob.resize(my_ef_frames);
    blob_nan.resize(my_ef_frames);
    for(i=0; i<my_ef_frames; i++)
    {
        blob[i].resize(num_g_y);
        blob_nan[i].resize(num_g_y);
        for(j=0; j<num_g_y; j++)
        {
            blob[i][j].resize(my_num_g_x,0);
            blob_nan[i][j].resize(my_num_g_x,1);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read binding events files and generate single frame blob data                                             //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(i=my_xi; i<=my_xf; i++) //loop over x
    {
        int ef_x = i - my_xi;

        for(j=0; j<num_g_y; j++) //loop over y
        {
            if(world_rank == 0)
            {
                printf("Working on x %d y %d \n",i,j);
            }

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                                                           //
            // Read in binding events                                                                                    //
            //                                                                                                           //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            string in_file_name = base_file_name_i + "_" + to_string(i) + "_" + to_string(j) + ".be";
            int result          = get_binding_events(in_file_name);

            if(result == 1)
            {
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
                                blob[ef_frame][j][ef_x] = time_line_lipid_nr[l];
                                blob_nan[ef_frame][j][ef_x] = 0;
                                break;
                            }
                        }
                        ef_frame = ef_frame + 1;
                    }
                }
            }
            else
            {
                num_g_y = num_g_y_store;
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Restore the header info                                                                                   //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ef_frames  = ef_frames_store;
    ef_dt      = ef_dt_store;
    num_lipids = num_lipids_store;
    num_g_x    = num_g_x_store;
    num_g_y    = num_g_y_store;
    APS        = APS_store;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects the lipid blobs data for a single frame                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Binding_events::get_blobs_frame(int this_frame,int world_size,int world_rank)
{
    int i = 0;
    int j = 0;

    MPI_Barrier(MPI_COMM_WORLD);

    blob_frame.resize(0);
    blob_nan_frame.resize(0);

    for(j=0; j<num_g_y; j++) //loop over y
    {
        iv1d current_row_grid = collect_and_clone_iv1d(world_size,world_rank,blob[this_frame][j]);
        iv1d current_row_nan  = collect_and_clone_iv1d(world_size,world_rank,blob_nan[this_frame][j]);

        blob_frame.push_back(current_row_grid);
        blob_nan_frame.push_back(current_row_nan);
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes the current frame blob data to file                                                  //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Binding_events::write_blobs_frame(string out_file_name)
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
                if(blob_nan_frame[j][k] == 0)
                {
                    fprintf(out_file," %10d",blob_frame[j][k]);
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
// This function analyzes lipid blob data and returns the center of a specified lipid                        //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
dv1d Binding_events::find_blobs_center(int this_lipid)
{
    int j     = 0;
    int k     = 0;
    int count = 0;
    dv1d center(2,0.0);

    for(j=0; j<num_g_y; j++) //loop over y-dimension
    {
        for(k=0; k<num_g_x; k++) //loop over x-dimension
        {
            if(blob_nan_frame[j][k] == 0) //check that grid point is not the protein
            {
                if(blob_frame[j][k] == this_lipid)
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


//define up to 4-dimensional vector of doubles
typedef vector<double>  dv1d;
typedef vector<dv1d>    dv2d;
typedef vector<dv2d>    dv3d;
typedef vector<dv3d>    dv4d;

//define up to 4-dimensional vector of ints
typedef vector<int>     iv1d;
typedef vector<iv1d>    iv2d;
typedef vector<iv2d>    iv3d;
typedef vector<iv3d>    iv4d;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This structure holds variables specific to the analysis program.. Any variable that can be declared at    //
// the start of the program is stored here.                                                                  //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct program_variables
{
    string program_description;                   //A breif description of the program

    string in_file_name;                          //Name of the input trajectory file
    string ref_file_name;                         //Name of the refernce file
    string lsq_index_file_name;                   //Name of the index file with atoms to do least squares fitting
    int stride;                                   //How many frames do we skip each time?
    int start_frame;                              //Dont read trajectory frames before this number
    int end_frame;                                //Dont read trajectory frames after this number
    int b_lsq;                                    //Do least squares fitting?
    int lsq_dim;                                  //Dimension of lsq fitting.
    int lsq_ref;                                  //Structure used for lsq fitting (0:ref 1:frame_0)

    //lipid mixing stuff
    string param_1_file_name;                     //Name of the parameter file 1
    string param_2_file_name;                     //Name of the parameter file 2
    string mix_file_name;                         //Name of the lipid mixing output file
    string lf_pdb_file_name;                      //Name of the lf_pdb_file
    string leaflet_finder_param_name;             //Name of the leaflet finder param file
    string map_1;                                 //Name of mapping atom 1
    string map_2;                                 //Name of mapping atom 2
    int b_lf_param;                               //Tells if the user included a lipid types parameter file
    int leaflet;                                  //Which leaflet should analysis be performed for (0:neither 1:upper 2:lower)
    int b_lf_pdb;                                 //Print pdb file with marked leaflets by beta value
    int num_lipids_1;                             //Number of lipids of the target type 1 in the target leaflet
    int num_lipids_2;                             //Number of lipids of the target type 2 in the target leaflet
    int mix_stride;                               //How often to print the mixing data to output file
    int out_data_format;                          //What format should the output data be in? (0:matrix 1:vector)
    int range;                                    //How many frames to compute running average over
    int b_report_binding_events;                  //Print to file the binding events as they are encountered
    int num_lipids;                               //How many lipids are in the target leaflet(s)
    int num_g_x;                                  //This is the number of grid points in x direction
    int num_g_y;                                  //This is the number of grid points in y direction
    int b_m1;                                     //m1 argument provided? 
    int b_m2;                                     //m2 argument provided?
    int b_APS;                                    //APS argument provided?
    int b_r;                                      //r argument provided?
    int b_cutoff;                                 //cutoff argument provided?
    int b_dist_proj;                              //do distance projection of solvation number?
    double cutoff;                                //Percentage of average rho used for excluding data
    double APS;                                   //This is the area of a grid square
    double radius;                                //Radius of the atom
    double box_x;                                 //user input box (x-dimension)
    double box_y;                                 //user input box (y-dimension) 
    double ef_dt;                                 //Effective delta t between frames
    double window_cutoff;                         //What percentage of window frames must a contact be present to be counted as a neighbor
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function initializes the variables held in program_variables                                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void initialize_program_variables(program_variables *p)
{
    p->stride      = 1;
    p->start_frame = 0;
    p->end_frame   = -1;
    p->b_lsq       = 0;
    p->lsq_dim     = 3;
    p->lsq_ref     = 0;

    //initialize program variables here
    p->mix_file_name           = "";
    p->lf_pdb_file_name        = "";
    p->leaflet                 = 0;
    p->b_lf_pdb                = 0;
    p->num_lipids_1            = 0;
    p->num_lipids_2            = 0;
    p->mix_stride              = 1000;
    p->ef_dt                   = 1;
    p->b_lf_param              = 0;
    p->cutoff                  = 0.0;
    p->APS                     = 0;
    p->radius                  = 0;
    p->box_x                   = 0;
    p->box_y                   = 0;
    p->range                   = 0;
    p->window_cutoff           = 0.8;
    p->b_report_binding_events = 0;
    p->num_lipids              = 0;
    p->num_g_x                 = 0;
    p->num_g_y                 = 0;
    p->out_data_format         = 0;
    p->b_m1                    = 0;
    p->b_m2                    = 0;
    p->b_APS                   = 0;
    p->b_r                     = 0;
    p->b_cutoff                = 0;
    p->b_dist_proj             = 0;

    //here we set the program description
    p->program_description = p->program_description + "Lipid Mixing is a program that measures the degree of lipid mixing as a simulation progresses. This is done by counting the percentage of lipids of a specific type to have been in a lipid's first shell. This percentage, which we call the mixing fraction, is computed for all lipids of a user-specified type and the average periodically reported.";
}


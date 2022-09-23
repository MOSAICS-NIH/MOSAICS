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
    string out_file_name;                         //Name of the output trajectory file
    string ref_file_name;                         //Name of the refernce file
    string lsq_index_file_name;                   //Name of the index file with atoms to do least squares fitting
    int stride;                                   //How many frames do we skip each time?
    int start_frame;                              //Dont read trajectory frames before this number
    int end_frame;                                //Dont read trajectory frames after this number
    int b_print;                                  //Print the output trajectory?
    int b_lsq;                                    //Do least squares fitting?
    int lsq_dim;                                  //Dimension of lsq fitting.
    int lsq_ref;                                  //Structure used for lsq fitting (0:ref 1:frame_0)

    //membrane thickness stuff
    string param_file_name;                       //Name of the parameter file
    string thk_file_name;                         //Name of the thickness output file
    string lf_pdb_file_name;                      //Name of the lf_pdb_file
    string leaflet_finder_param_name;             //Name of the leaflet finder param file
    int b_lf_param;                               //Tells if the user included a lipid types parameter file
    int leaflet;                                  //Which leaflet should analysis be performed for (0:neither 1:upper 2:lower)
    int out_data_format;                          //What format should the output data be in? (0:matrix 1:vector)
    int b_stdev;                                  //Compute the z-coord standard deviation?
    int b_clean;                                  //Delete single frame rmsd files after computing average?
    int b_lf_pdb;                                 //Print pdb file with marked leaflets by beta value
    int count_all;                                //Include all lipids in the opposing leaflet
    double upper_cutoff;                          //Dump pdb of frames with distance below this value 
    double lower_cutoff;                          //Dump pdb of frames with distance above this value
    double cutoff;                                //Percentage of average rho used for excluding data
    double APS;                                   //This is the area of a grid square
    double radius;                                //Radius of the atom
    double box_x;                                 //user input box (x-dimension)
    double box_y;                                 //user input box (y-dimension)
    double xy_cutoff;                             //How far in xy before counting lipids in min dist 
    double bin_width;                             //Width of bins
    double temp;                                  //Temperature of the simulation
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
    p->b_print     = 0;
    p->b_lsq       = 0;
    p->lsq_dim     = 3;
    p->lsq_ref     = 0;

    p->leaflet         = 1; 
    p->out_data_format = 0;
    p->cutoff          = 0.0;
    p->APS             = 0.005;
    p->radius          = 0.23;
    p->box_x           = 0;
    p->box_y           = 0;
    p->b_stdev         = 0;
    p->b_clean         = 0;
    p->b_lf_pdb        = 0;
    p->xy_cutoff       = 0.5;
    p->b_lf_param      = 0;
    p->upper_cutoff    = 0;
    p->lower_cutoff    = 0;
    p->bin_width       = 0.1;
    p->temp            = 310;
    p->count_all       = 0;

    //here we set the program description
    p->program_description = "Membrane Thickness is an analysis program designed to probe the membrane thickness and project this information onto the XY plane. The program works by finding pairs of lipid in opposing leaflets and measuring the delta z. ";  
}


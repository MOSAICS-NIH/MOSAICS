
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

    //nearest neighbors stuff
    string param_1_file_name;                     //Name of the parameter file 1
    string param_2_file_name;                     //Name of the parameter file 2
    string nbrs_file_name;                        //Name of the ld output file
    string lf_pdb_file_name;                      //Name of the output pdb file with leaflets indicated by B-factor
    string leaflet_finder_param_name;             //Name of the leaflet finder param file
    int b_lf_param;                               //Tells if the user included a lipid types parameter file
    int leaflet;                                  //Which leaflet? (0:both 1:upper 2:lower)
    int b_lf_pdb;                                 //Print pdb with leaflets indicated?
    int b_stdev;                                  //compute the fsd standard deviation?
    int b_clean;                                  //remove single frame data after stdev computation
    double APS;                                   //This is the area of a grid square
    double radius;                                //Radius of the atom
    double cutoff;                                //Percentage of average rho used for excluding data
    double box_x;                                 //Grid x dimension
    double box_y;                                 //Grid y dimension
    double box_z;                                 //Grid z dimension
    double local_rad;                             //Radius for counting lipids
    double ex_val;                                //Set excluded lattice points to this value
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

    //initialize program variables here
    p->leaflet         = 0;
    p->APS             = 0;
    p->radius          = 0;
    p->cutoff          = 0;
    p->box_x           = 0;
    p->box_y           = 0;
    p->box_z           = 0;
    p->local_rad       = 0;
    p->b_lf_pdb        = 0;
    p->b_lf_param      = 0;
    p->ex_val          = 0.0;

    //here we set the program description
    p->program_description = "Nearest Neighbors 3d is an analysis tool that probes the lipid packing density as a function of position in the XYZ. This is done by computing the number of lipid nearest-neighbors found within a user specified distance of each target lipid.";
 }


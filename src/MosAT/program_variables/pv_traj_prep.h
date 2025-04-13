
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

    //2d_kinetics stuff
    string pf_pdb_file_name;                      //Name of the output pdb file with protein indicated by B-factor
    string protein_finder_param_name;             //Name of the protein finder param file
    string selection_text_file_name;              //Name of the selection card with the selection text
    string param_file_name;                       //Name of selection card with list of operations
    string slim_index_file_name;                  //Name of index file with id's for atoms to keep
    int b_pf_param;                               //Tells if the user included a protein types parameter file
    int b_pf_pdb;                                 //Print the protein finder pdb?
    int b_test;                                   //Print info for testing molecule definitions
    int stop;                                     //Stop printing -test info after this many molecules
    int b_slim;                                   //Do we remove certain atoms from the trajectory? 
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
    p->b_pf_pdb    = 0;
    p->b_test      = 0;
    p->stop        = 0;
    p->b_slim      = 0;

    //here we set the program description
    p->program_description = "Traj Prep is a tool used to prepare trajectories for analysis. This tool helps the user perform basic operation such as fitting, centering, translating, and wrapping. To use the tool, the user provides a recipe telling which operations to perform in what order.";
}


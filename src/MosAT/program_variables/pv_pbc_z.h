
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
    string sf_pdb_file_name;                      //Name of the output pdb file with sol indicated by B-factor
    string solvent_finder_param_name;             //Name of the solvent finder param file
    string record_file_name;                      //Name of the output file with the jump record
    int stride;                                   //How many frames do we skip each time?
    int start_frame;                              //Dont read trajectory frames before this number
    int end_frame;                                //Dont read trajectory frames after this number
    int b_print;                                  //Print the output trajectory?
    int b_lsq;                                    //Do least squares fitting?
    int lsq_dim;                                  //Dimension of lsq fitting.
    int lsq_ref;                                  //Structure used for lsq fitting (0:ref 1:frame_0)
    int b_sf_param;                               //Tells if the user included a solvent types parameter file
    int b_sf_pdb;                                 //Print solvent finder pdb?
    int b_record;                                 //Report the jump records?
    int jump_count;                               //How many jumps detected?

    //add program variables here
    double cutoff;                                //Percentage of the box an atom must move by to be counted as a jump
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
    p->cutoff     = 0.50;
    p->b_sf_param = 0;
    p->b_sf_pdb   = 0;
    p->b_record   = 0;
    p->jump_count = 0;

    //here we set the program description
    p->program_description = "PBC_Z is a program designed to fix periodic boundary problems, specifically the jumping of lipids across the z-axis due to membrane bending/fluctuations. PBC_Z works by checking lipids in each frame for jumps by comparing the z-coordinate to that from the previous frame. If the difference is too big (compared to the box size) a jump is assumed.";
}


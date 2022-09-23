
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

    FILE *orientation_file;                       //File for printing the orientation data to
    string upper_file_name;                       //Name of the upper index file
    string lower_file_name;                       //Name of the lower index  file
    string orientation_file_name;                 //Name of the orientatio file
    int upper_capacity;                           //Number of items in the upper index file
    int lower_capacity;                           //Number of items in the lower index file
    int print_vec_pdb;
    int print_vec;                                //Print the orientation vector?
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
    p->upper_capacity = 0;
    p->lower_capacity = 0;

    //here we set the program description
    p->program_description = "Protein Orientation is a program the characterizes the orientation vector for a protein. The orientation vector is defined by choosing 2 sets of atoms. The center of mass of each set is then computed and the vector connecting the 2 centers gives the orientation vector. The angle of this vector is then computed relative to the z-axis. The vector is also projected onto the XY plane and the angle relative to the y-axis is computed. Output is a text file with the angles as a function of time.";
}


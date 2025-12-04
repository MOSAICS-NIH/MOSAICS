
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

    string hb_file_name;                          //Name of the hb output file
    string acc_file_name;                         //Name of the acceptor atoms file
    string don_file_name;                         //Name of the donor atoms file
    string bond_file_name;                        //Name of the bonds list file
    string n1_file_name;                          //Name of the group 1 index file
    string n2_file_name;                          //Name of the group 2 index  file
    string exclude_file_name;                     //Name of file with lids of excluded acceptor/donor pairs
    int b_test;                                   //Write info for checking hydrogen bonds?
    int b_exclude;                                //Did the user provide and exclude list
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
    p->b_test          = 0;
    p->b_exclude       = 0;

    //here we set the program description
    p->program_description = "H Bonds is a program designed for quantifying hydrogen bonds between 2 groups of atoms. The number of hydrogen bonds formed is reported as a function of simulation time."; 
}


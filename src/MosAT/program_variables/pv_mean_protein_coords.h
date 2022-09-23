
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

    FILE *mpc_file;                               //File for writing out the average coords
    FILE *pf_pdb_file;                            //File for writing the output pdb file with protein indicated by B-factor
    string mpc_file_name;                         //Name of the mpc output file
    string par_file_name;                         //Name of the par file
    string pf_pdb_file_name;                      //Name of the output pdb file with leaflets indicated by B-factor
    string protein_finder_param_name;             //Name of the protein finder param file
    int b_pf_param;                               //Tells if the user included a protein types parameter file
    int prot_size;                                //How many atoms in the protein
    int print_prot_pdb;                           //Print a pdb with the protein highlighted?
    int b_mean_dist;                              //Compute the average distance from the average coords?
    rvec *r_avg;                                  //Holds thd average coords
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

    p->prot_size      = 0;
    p->print_prot_pdb = 0;
    p->b_mean_dist    = 0;

    //here we set the program description
    p->program_description = "Mean Protein Coords is a program that computes the time average coordinates for the protein in a system.";
}


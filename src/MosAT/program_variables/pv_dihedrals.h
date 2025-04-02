
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

    string pf_pdb_file_name;                      //Name of the output pdb file with protein indicated by B-factor
    string protein_finder_param_name;             //Name of the protein finder param file
    string target_res_file_name;                  //Name of index file with target residues
    string dihedrals_file_name;                   //Name of index file with dihedral definitions
    string dih_file_name;                         //Name of the output file with dihedral angles
    string traj_2_file_name;                      //Input trajectory file with the target dihedral angles
    string target_res_cmp_file_name;              //List of dihedral angles that should be compared to the list in -res
    string exclude_file_name;                     //List of residues to be removed from target residues list
    int b_pf_param;                               //Tells if the user included a protein types parameter file
    int b_pf_pdb;                                 //Print the protein finder pdb?
    int type;                                     //Which dihedral angle to measure?
    int b_noisy;                                  //Tell the program to write info about the selected residues/atoms
    int b_traj_b;                                 //Did the user provide a second trajectory?
    int b_cmp;                                    //Compare 2 sets of angles?
    int b_test;                                   //Print Pymol commands to check angles
    int b_ex;                                     //Tells whether a new list should be generated or not
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function initializes the variables held in program_variables                                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void initialize_program_variables(program_variables *p)
{
    //initialize program variables here

    p->stride      = 1;
    p->start_frame = 0;
    p->end_frame   = -1;
    p->b_print     = 0;
    p->b_lsq       = 0;
    p->lsq_dim     = 3;
    p->lsq_ref     = 0;

    p->b_pf_pdb        = 0;
    p->b_pf_param      = 0;
    p->type            = 0;
    p->b_noisy         = 0;
    p->b_traj_b        = 0;
    p->b_cmp           = 0;
    p->b_ex            = 0;

    //here we set the program description
    p->program_description = p->program_description + "Dihedrals is an analysis tool used for measuring dihedral angels in a ";
    p->program_description = p->program_description + "protein simulation. To use the program the user provides a list of residues ";
    p->program_description = p->program_description + "for which the dihedral angles are to be measured and the type (phi vs psi vs chi1)";
    p->program_description = p->program_description + "The program will then compute the angles for each frame and display them to screen";
}



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

    double cutoff_z;                              //how far away (z) do we still count the atoms
    double cutoff_xy;                             //how far away (xy) do we still count the atoms

    FILE *pf_pdb_file;                            //File for writing the output pdb file with protein indicated by B-factor
    string pf_pdb_file_name;                      //Name of the output pdb file with protein indicated by B-factor
    string protein_finder_param_name;             //Name of the protein finder param file
    string atom_type;                             //Atom type when finding pairs
    int b_pf_param;                               //Tells if the user included a protein types parameter file
    int b_pf_pdb;                                 //Print pdb with protein indicated as beta factor
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

    p->cutoff_z   = 0.1;
    p->cutoff_xy  = 0.1;
    p->b_pf_pdb   = 0;
    p->b_pf_param = 0;

    //here we set the program description
    p->program_description = p->program_description + "Atoms in 2 Planes is a program used for selecting two sets of atoms such that when computing the center of mass of each group the vector v1 – v2 will align parallel with the z-axis. This is accomplished by selecting pairs of atoms that are far apart in the z-direction and close to each other in the x-y plane."; 
}


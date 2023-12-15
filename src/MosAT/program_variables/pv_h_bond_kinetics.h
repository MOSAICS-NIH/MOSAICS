
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

    string hbk_file_name;                         //Name of the h_bonds kinetics output file
    string pf_pdb_file_name;                      //Name of the output pdb file with protein indicated by B-factor
    string lip_a_file_name;                       //Name of the lipid acceptor atoms file
    string lip_d_file_name;                       //Name of the lipid donor atoms file
    string prot_a_file_name;                      //Name of the protein acceptor atoms file
    string prot_d_file_name;                      //Name of the protein donor atoms file
    string bond_file_name;                        //Name of the bonds list file
    string protein_finder_param_name;             //Name of the protein finder param file
    string selection_text_file_name;              //Name of the selection card with the selection text
    string report_file_name;                      //Name of the selection card holding the atom number for reporting h-bonds
    string target_res;                            //Target residue type
    string exclude_file_name;                     //Name of file with lids of excluded acceptor/donor pairs
    string be_file_name;                          //Name of the binding events file 
    int target_x;                                 //Target x lattice point
    int target_y;                                 //Target y lattice point
    int b_pf_param;                               //Tells if the user included a protein types parameter file
    int b_pf_pdb;                                 //Print the protein finder pdb?
    int b_test;                                   //Write info for checking hydrogen bonds?
    int b_report;                                 //Report frame number when a target h-bond is found
    int b_sel_text;                               //Refine the selection of protein atoms?
    int size_x;                                   //How many items in the x-direction of the contact matrix?
    int size_y;                                   //How many items in the y-direction of the contact matrix?
    int test;                                     //Print PyMol select commands for the contacts found
    int range;                                    //Half-width of the noise filter
    int dump;                                     //Dump residence time for all contacts on last trajectory frame
    int resi;                                     //The target residue 
    int b_target_res;                             //Tells if a target lipid type was provided or not
    int b_exclude;                                //Did the user provide and exclude list
    int b_be;                                     //Tells if a binding events file was provided or not
    int b_x;                                      //Tells if the target x was specified
    int b_y;                                      //Tells if the target y was specified

    int min;                                      //First atom of the bound lipid
    int max;                                      //Last atom of the bound lipid
    int resi_size;                                //How many atoms in the target lipid type
    int norm_factor;                              //Used to get the average
    double delta_t;                               //Effective time step between traj frame accounting for stride
    double bin_width;                             //The histogram bin width (ps)
    double max_dash_rad;                          //dash thickness for contacts with frequency of 1
    double min_dash_rad;                          //dash thickness for contacts with frequency of 0
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
    p->b_pf_pdb        = 0;
    p->b_test          = 0;
    p->b_sel_text      = 0;
    p->b_report        = 0;
    p->size_x          = 0;
    p->size_y          = 0;
    p->range           = 0;
    p->dump            = 0;
    p->bin_width       = 1000.0;
    p->max_dash_rad    = 0.3;
    p->min_dash_rad    = 0.0;
    p->min             = -1;
    p->max             = -1;
    p->norm_factor     = 0;
    p->resi_size       = 0;
    p->b_exclude       = 0;
    p->b_be            = 0;
    p->b_x             = 0;
    p->b_y             = 0;
    p->target_x        = -1;
    p->target_y        = -1;

    //here we set the program description
    p->program_description = "H-Bond Kinetics is an analysis tool used for characterizing hydrogen bond kinetics between one or more bound lipids and the protein. The bonds are identified by drawing lines in PyMOL thus showing all bonds made with the bound lipid."; 
}


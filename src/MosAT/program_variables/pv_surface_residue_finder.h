
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

    //add program variables here
    FILE *surface_file;                           //File for writing the final surface atom data
    FILE *surface_pdb_file;                       //File for writing the final surface atom data to pdb
    FILE *lf_pdb_file;                            //File for writing the output pdb file with leaflets indicated by B-factor
    FILE *pf_pdb_file;                            //File for writing the output pdb file with protein indicated by B-factor
    string selection_text_file_name;              //Name of the selection card with the selection text
    string histo_file_name;                       //Name of the histogram file
    string surface_pdb_file_name;                 //Name of the surface pdb file
    string lf_pdb_file_name;                      //Name of the output pdb file with leaflets indicated by B-factor
    string pf_pdb_file_name;                      //Name of the output pdb file with protein indicated by B-factor
    string param_file_name;                       //Name of the parameter file
    string leaflet_finder_param_name;             //Name of the leaflet finder param file
    string protein_finder_param_name;             //Name of the protein finder param file
    int b_pf_param;                               //Tells if the user included a protein types parameter file
    int b_lf_param;                               //Tells if the user included a lipid types parameter file
    int mem_size;                                 //how many atoms make the membrane
    int prot_size;                                //Number of atoms making the protein
    int leaflet;                                  //Which leaflet? (for this program always set to 0)
    int b_lf_pdb;                                 //Print pdb with leaflets indicated by beta factor
    int b_pf_pdb;                                 //Print pdb with protein indicated by beta factor
    int bin;                                      //Bin width for histogram
    int b_sel_text;                               //Refine the selection of protein atoms?
    double cutoff;                                //Above what percentage is significant?
    double contact_cutoff;                        //How far away before no longer counted as a contact
    double bulk;                                  //Lipid ratio in the bulk
    double screen_dist;                           //Screening distance for checking resi centers
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
    p->mem_size       = 0;
    p->prot_size      = 0;
    p->leaflet        = 0;
    p->prot_size      = 0;
    p->cutoff         = 0;
    p->contact_cutoff = 0.6;
    p->bulk           = 0;
    p->b_lf_pdb       = 0;
    p->b_pf_pdb       = 0;
    p->b_lf_param     = 0;
    p->bin            = 1;
    p->b_sel_text     = 0;
    p->screen_dist    = 4.0;

    //here we set the program description
    p->program_description = "Surface Residue Finder is a program that selects the surface atoms of a protein. This is done by counting the number of lipids that form one or more contacts with a given protein atom. Then, atopms with too few contacting lipids are excluded from the surface. The output is a pdb with surface residues indicated by the beta factor as well as PyMOL commands used for selecting the surace and core atoms.";
}


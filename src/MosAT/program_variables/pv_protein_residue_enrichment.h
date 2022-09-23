
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
    FILE *enrich_file;                            //File for writing the final enrichment data
    FILE *enrich_pdb_file;                        //File for writing the final enrichment data to pdb
    FILE *lf_pdb_file;                            //File for writing the output pdb file with leaflets indicated by B-factor
    FILE *pf_pdb_file;                            //File for writing the output pdb file with protein indicated by B-factor
    string enrich_file_name;                      //Name of the enrichment file
    string enrich_pdb_file_name;                  //Name of the enrichment pdb file
    string lf_pdb_file_name;                      //Name of the output pdb file with leaflets indicated by B-factor
    string pf_pdb_file_name;                      //Name of the output pdb file with protein indicated by B-factor
    string param_1_file_name;                     //Name of the param_1 file
    string param_2_file_name;                     //Name of the param_2 file
    string leaflet_finder_param_name;             //Name of the leaflet finder param file
    string protein_finder_param_name;             //Name of the protein finder param file
    int b_pf_param;                               //Tells if the user included a protein types parameter file
    int b_lf_param;                               //Tells if the user included a lipid types parameter file
    int lip_t1_capacity;                          //Number of atoms contained in the lipid_file_1
    int lip_t2_capacity;                          //Number of atoms contained in the lipid_file_2
    int target_atoms_capacity;                    //Number of atoms contained in the target atoms index file
    int mem_size;                                 //how many atoms make the membrane
    int prot_size;                                //Number of atoms making the protein
    int prot_res_count;                           //How many protein residues are there
    int lip_count_1;                              //How many lipids are there
    int lip_count_2;                              //How many lipids are there
    int leaflet;                                  //Which leaflet? (for this program always set to 0)
    int b_lf_pdb;                                 //Print pdb with leaflets indicated as beta factor
    int b_pf_pdb;                                 //Print pdb with protein indicated as beta factor
    double cutoff;                                //Above what percentage is significant?
    double contact_cutoff;                        //How far away before no longer counted as a contact
    double bulk;                                  //Lipid ratio in the bulk
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
    p->lip_t1_capacity       = 0;
    p->lip_t2_capacity       = 0;
    p->target_atoms_capacity = 0;
    p->mem_size              = 0;
    p->prot_size             = 0;
    p->prot_res_count        = 0;
    p->lip_count_1           = 0;
    p->lip_count_2           = 0;
    p->leaflet               = 0;
    p->prot_size             = 0;
    p->cutoff                = 0;
    p->contact_cutoff        = 0.6;
    p->bulk                  = 0;
    p->b_lf_pdb              = 0;
    p->b_pf_pdb              = 0;
    p->b_lf_param            = 0;

    //here we set the program description
    p->program_description = "Protein Residue Enrichment is a program that computes the lipid enrichment factor for each residue of the protein. That is, the number of lipids contacting the residue are counted for each protein residue. This is done for two groups of lipids A and B. Then the enrichment factor is computed using the number of contacting lipids and the ratio of lipids A and B in the bulk. Output is a pdb with the percent enrichment given as the beta factor.";
}


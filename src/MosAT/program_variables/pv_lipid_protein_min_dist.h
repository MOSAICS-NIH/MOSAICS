
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

    string lpmd_file_name;                        //Name of the lpmd output file
    string pf_pdb_file_name;                      //Name of the output pdb file with protein indicated by B-factor
    string protein_finder_param_name;             //Name of the protein finder param file
    string selection_text_file_name;              //Name of the file with the selection text
    string selection_card_file_name;              //Name of the selection card 
    string lf_pdb_file_name;                      //Name of the lf_pdb_file
    string leaflet_finder_param_name;             //Name of the leaflet finder param file
    string be_file_name;                          //Name of the binding events file 
    string target_res;                            //Target residue type
    int target_x;                                 //Target x lattice point
    int target_y;                                 //Target y lattice point
    int b_lf_param;                               //Tells if the user included a lipid types parameter file
    int leaflet;                                  //Which leaflet should analysis be performed for (0:neither 1:upper 2:lower)
    int b_lf_pdb;                                 //Print pdb file with marked leaflets by beta value
    int b_pf_param;                               //Tells if the user included a protein types parameter file
    int b_pf_pdb;                                 //Print the protein finder pdb?
    int b_target_res;                             //Tells if a target lipid type was provided or not
    int b_be;                                     //Tells if a binding events file was provided or not
    int b_x;                                      //Tells if the target x was specified
    int b_y;                                      //Tells if the target y was specified
    int size_x;                                   //How many items in the x-direction of the contact matrix?
    int size_y;                                   //How many items in the y-direction of the contact matrix?
    int frag_size;                                //How many frames to store contact matrix at any given time
    int clean;                                    //Remove contact matrices when analysis finishes
    int test;                                     //Print PyMol select commands for the contacts found
    int tier;                                     //Write the top clusters to a PDB file
    int resi;                                     //The target residue
    int resi_size;                                //The size of the target residue
    int min;                                      //First atom of the bound lipid
    int max;                                      //Last atom of the bound lipid
    int norm_factor;                              //Used to get the average
    double contact_cutoff;                        //Cutoff for counting contacts
    double cutoff;                                //Cutoff factor for counting members
    double screen_dist;                           //Screen lipids whose centers are within this disatnce (nm) of the protein                           
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
    p->contact_cutoff  = 0.6;
    p->size_x          = 0;
    p->size_y          = 0;
    p->cutoff          = 0.8;
    p->frag_size       = 100; 
    p->clean           = 0;
    p->test            = 0;
    p->tier            = 100;
    p->resi            = 0;
    p->resi_size       = 0;
    p->leaflet         = 0;
    p->b_lf_pdb        = 0;
    p->b_lf_param      = 0;
    p->screen_dist     = 4.0;
    p->min             = -1;
    p->max             = -1;
    p->norm_factor     = 0;
    p->b_be            = 0;
    p->b_x             = 0;
    p->b_y             = 0;
    p->target_x        = -1; 
    p->target_y        = -1;

    //here we set the program description
    p->program_description = "lipid Protein Min Dist is a program that characterizes for each atom of one or more bound lipids the mean minimum distance to the protein.";   
}



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

    //pdb editor stuff
    string selection_text_file_name;              //Name of file with atom selection text
    string lf_pdb_file_name;                      //Name of the output pdb file with leaflets indicated by B-factor
    string pf_pdb_file_name;                      //Name of the output pdb file with protein indicated by B-factor
    string sf_pdb_file_name;                      //Name of the output pdb file with sol indicated by B-factor
    string leaflet_finder_param_name;             //Name of the leaflet finder param file
    string protein_finder_param_name;             //Name of the protein finder param file
    string solvent_finder_param_name;             //Name of the solvent finder param file
    string chain_id_string;                       //Desired value for chain id
    string element;                               //Desired value for element
    string this_res_name;                         //Desired value for res_name
    string this_atom_name;                        //Desired value for atom_name
    char   chain_id;                              //Desired value for chain id
    int b_sf_param;                               //Tells if the user included a solvent types parameter file
    int b_pf_param;                               //Tells if the user included a protein types parameter file
    int b_lf_param;                               //Tells if the user included a lipid types parameter file
    int leaflet;                                  //Which leaflet? (0:both 1:upper 2:lower)
    int b_lf_pdb;                                 //Print the leaflet finder pdb?
    int b_pf_pdb;                                 //Print the protein finder pdb?
    int b_sf_pdb;                                 //Print solvent finder pdb?
    int b_b_factor;                               //Is the b_facor modified
    int b_occupancy;                              //Is the occupancy modified
    int b_chain_id;                               //Is the chain id modified
    int b_element;                                //Is the element modified
    int b_res_name;                               //Is the res_name modified
    int b_atom_name;                              //Is the atom_name modified
    double b_factor;                              //Desired value for B factor
    double occupancy;                             //Desired value for the occupancy
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

    p->leaflet     = 0;
    p->b_pf_pdb    = 0;
    p->b_sf_pdb    = 0;
    p->b_lf_pdb    = 0;
    p->b_lf_param  = 0;
    p->b_pf_param  = 0;
    p->b_sf_param  = 0;
    p->b_factor    = 0.0;
    p->occupancy   = 0.0;
    p->b_b_factor  = 0;
    p->b_occupancy = 0;
    p->b_chain_id  = 0;
    p->b_element   = 0;
    p->b_res_name  = 0;
    p->b_atom_name = 0;

    //here we set the program description
    p->program_description = p->program_description + "PDB Editor is an analysis tool used for setting items in a PDB file (B factor, occupancy, chain id, element) for a group of atoms using a selection text.";
}


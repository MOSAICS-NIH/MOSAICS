
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

    string lpsb_file_name;                        //Name of the lpsb output file
    string param_file_name;                       //Name of the lipid types parameter file
    string lf_pdb_file_name;                      //Name of the output pdb file with leaflets indicated by B-factor
    string pf_pdb_file_name;                      //Name of the output pdb file with protein indicated by B-factor
    string lip_p_file_name;                       //Name of the lipid positive atoms file
    string lip_n_file_name;                       //Name of the lipid negative atoms file
    string prot_p_file_name;                      //Name of the protein positive atoms file
    string prot_n_file_name;                      //Name of the protein negative atoms file
    string leaflet_finder_param_name;             //Name of the leaflet finder param file
    string protein_finder_param_name;             //Name of the protein finder param file
    string lipid_charge_file_name;                //Name of the selection card with the lipid charge groups
    string protein_charge_file_name;              //Name of the selection card with the protein charge groups
    int b_sf_param;                               //Tells if the user included a solvent types parameter file
    int b_pf_param;                               //Tells if the user included a protein types parameter file
    int b_lf_param;                               //Tells if the user included a lipid types parameter file
    int leaflet;                                  //Which leaflet? (0:both 1:upper 2:lower)
    int out_data_format;                          //What format should the output data be in? (0:matrix 1:vector)
    int b_lf_pdb;                                 //Print the leaflet finder pdb?
    int b_pf_pdb;                                 //Print the protein finder pdb?
    int b_sf_pdb;                                 //Print solvent finder pdb?
    int b_stdev;                                  //Compute the z-coord standard deviation?
    int b_clean;                                  //Delete single frame rmsd files after computing average?
    int b_test;                                   //Write info for checking salt bridges?
    double APS;                                   //This is the area of a grid square
    double radius;                                //Radius of the atom
    double cutoff;                                //Percentage of average rho used for excluding data
    double box_x;                                 //Grid x dimension
    double box_y;                                 //Grid y dimension
    double contact_cutoff;                        //Cutoff for counting contacts
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
    p->leaflet         = 0;
    p->out_data_format = 0;
    p->APS             = 0;
    p->radius          = 0;
    p->cutoff          = 0;
    p->box_x           = 0;
    p->box_y           = 0;
    p->b_pf_pdb        = 0;
    p->b_sf_pdb        = 0;
    p->contact_cutoff  = 0.6;
    p->b_stdev         = 0;
    p->b_clean         = 0;
    p->b_lf_pdb        = 0;
    p->b_lf_param      = 0;
    p->b_test          = 0;

    //here we set the program description
    p->program_description = "Lipid Salt Bridges is an analysis program desiged to compute the number of salt bridges formed between lipid molecules and the protein. This number is then projected onto the XY plane."; 
}


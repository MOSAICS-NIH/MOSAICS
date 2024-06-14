
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
//                                                                                                           d
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

    //greasy waters stuff
    string lim_file_name;                         //Name of the lipid immobilization output file
    string lf_pdb_file_name;                      //Name of the output pdb file with leaflets indicated by B-factor
    string crd_file_name;                         //Name of the selection card with lipid types
    string leaflet_finder_param_name;             //Name of the leaflet finder param file
    int b_lf_param;                               //Tells if the user included a lipid types parameter file
    int leaflet;                                  //Which leaflet? (0:both 1:upper 2:lower)
    int b_lf_pdb;                                 //Print pdb with indicated leaflets by beta factor?
    int num_g_x;                                  //Number of lattice points in the x-direction
    int num_g_y;                                  //Number of lattice points in the y-direction
    int odf;                                      //Grid data format
    double APS;                                   //This is the area of a grid square
    double radius;                                //Radius of the atom
    double box_x;                                 //Grid x dimension
    double box_y;                                 //Grid y dimension
    double cell_size;                             //Size of a lattice cell (nm)
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
    p->APS             = 0;
    p->radius          = 0;
    p->box_x           = 0;
    p->box_y           = 0;
    p->b_lf_pdb        = 0;
    p->b_lf_param      = 0;
    p->num_g_x         = 0;
    p->num_g_y         = 0;
    p->cell_size       = 0.0;
    p->odf             = 0;

    //here we set the program description
    p->program_description = "Lipid Immobilization is an analysis program that traces the path taken for each lipid molecule. Viewing this data can be used to identify lipid molecules that become bound within a simulation."; 
	   
}


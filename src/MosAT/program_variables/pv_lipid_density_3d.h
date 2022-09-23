
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

    //density stuff
    string rho_file_name;                         //Name of the lipid density output file
    string lf_pdb_file_name;                      //Name of the output pdb file with leaflets indicated by B-factor
    string pf_pdb_file_name;                      //Name of the output pdb file with protein indicated by B-factor
    string param_file_name;                       //Name of the param index file
    string leaflet_finder_param_name;             //Name of the leaflet finder param file
    string protein_finder_param_name;             //Name of the protein finder param file
    int b_lf_param;                               //Tells if the user included a lipid types parameter file
    int b_pf_param;                               //Tells if the user included a protein parameter file
    int leaflet;                                  //Which leaflet? (0:both 1:upper 2:lower)
    int b_lf_pdb;                                 //Print pdb with indicated leaflets by beta factor?
    int b_pf_pdb;                                 //Print pdb with indicated protein by beta factor?
    int b_stdev;                                  //Compute the z-coord standard deviation?
    int b_clean;                                  //Delete single frame rmsd files after computing average?
    int b_dist;                                   //Add density only local to the protein? 
    int ex_val;                                   //Set excluded lattice points to this value
    double dist_cutoff;                           //How close to the protien must the lipids be to add density?
    double APS;                                   //This is the area of a grid square
    double radius;                                //Radius of the atom
    double cutoff;                                //Percentage of average rho used for excluding data
    double box_x;                                 //Grid x dimension
    double box_y;                                 //Grid y dimension
    double box_z;                                 //Grid z dimension
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
    p->cutoff          = 0;
    p->box_x           = 0;
    p->box_y           = 0;
    p->box_z           = 0;
    p->b_lf_pdb        = 0;
    p->b_pf_pdb        = 0;
    p->b_stdev         = 0;
    p->b_clean         = 0;
    p->b_lf_param      = 0;
    p->b_pf_param      = 0;
    p->b_dist          = 0;
    p->dist_cutoff     = 0;
    p->ex_val          = 0;

    //here we set the program description
    p->program_description = "Lipid Density 3D is a program that deposits density to a 3-dimensional lattice around select lipid atoms.";
}


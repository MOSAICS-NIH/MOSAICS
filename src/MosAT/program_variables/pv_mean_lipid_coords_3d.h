
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

    FILE *param_file;                             //File for reading in the reference lipid data
    string rho_file_name;                         //Name of the rho output file 
    string rho_t_file_name;                       //Name of the rho output file
    string mlc_file_name;                         //Name of the mlc output file
    string param_file_name;                       //Name of the lipid types parameter file
    string lf_pdb_file_name;                      //Name of the output pdb file with leaflets indicated by B-factor
    string map_1;                                 //Name of mapping atom 1
    string map_2;                                 //Name of mapping atom 2
    string target_lip;                            //Name of target lipid
    string leaflet_finder_param_name;             //Name of the leaflet finder param file
    int b_lf_param;                               //Tells if the user included a lipid types parameter file
    int leaflet;                                  //Which leaflet? (0:both 1:upper 2:lower)
    int num_g_x;                                  //This is the number of grid points in x direction
    int num_g_y;                                  //This is the number of grid points in y direction
    int num_g_z;                                  //This is the number of grid points in z direction
    int b_lf_pdb;                                 //Print the leaflet finder pdb?
    int grid_stride;                              //How many grid points to skip when printing mean lipids to pdb 
    double APS;                                   //This is the area of a grid square
    double radius;                                //Radius of the atom
    double cell_size;                             //This is the lengh between grid points
    double cutoff;                                //Percentage of average rho used for excluding data
    double box_x;                                 //Grid x dimension
    double box_y;                                 //Grid y dimension
    double box_z;                                 //Grid z dimension
    double ex_val;                                //Set excluded lattice points to this value
    int size_x;                                   //How many grid points in x-direction to include in the pdb (after g_stride)    
    int size_y;                                   //How many grid points in y-direction to include in the pdb (after g_stride)
    int ef_size_grid;                             //Total number of grid points after accounting for g_stride
    int b_rho_t;                                  //Tells if rho_t should be used  for excluding insignificant data
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
    p->num_g_x         = 0;
    p->num_g_y         = 0;
    p->num_g_z         = 0;
    p->APS             = 0;
    p->radius          = 0;
    p->cell_size       = 0;
    p->cutoff          = 0;
    p->box_x           = 0;
    p->box_y           = 0;
    p->box_z           = 0;
    p->b_lf_pdb        = 0;
    p->grid_stride     = 1;
    p->b_lf_param      = 0;
    p->size_x          = 0;
    p->size_y          = 0;
    p->ef_size_grid    = 0;
    p->ex_val          = 0.0;
    p->b_rho_t         = 0;

    //here we set the program description
    p->program_description = p->program_description + "Mean Lipid Coords 3d is an analysis tool used to compute the time averaged lipid coordinates based on the lipid's position in XYZ. ";
}



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
    string param_1_file_name;                     //Name of the parameter file 1
    string param_2_file_name;                     //Name of the parameter file 2
    string lf_pdb_file_name;                      //Name of the output pdb file with leaflets indicated by B-factor
    string enrich_file_name;                      //Name of the enrich file
    string leaflet_finder_param_name;             //Name of the leaflet finder param file
    int b_lf_param;                               //Tells if the user included a lipid types parameter file
    int num_g_x;                                  //This is the number of grid points in x direction
    int num_g_y;                                  //This is the number of grid points in y direction
    int leaflet;                                  //Which leaflet? (0:both 1:upper 2:lower)
    int mem_size;                                 //Size of the target membrane
    int print_stride;                             //How frequeny to print single frame data?
    int b_lf_pdb;                                 //Print pdb file with marked leaflets by beta value
    int num_lip_t1;                               //How many lipids of type 1 are in the selected leaflet
    int num_lip_t2;                               //How many lipids of type 2 are in the selected leaflet
    int num_lipids;                               //Number of lipids
    double cutoff;                                //Percentage of average rho used for excluding data
    double APS;                                   //This is the area of a grid square
    double radius;                                //Radius of the atom
    double cell_size;                             //This is the lengh between grid points
    double bulk;                                  //Ratio of the lipids in the bulk
    double box_x;                                 //Grid x dimension
    double box_y;                                 //Grid y dimension
    double box_z;                                 //Grid z dimension
    double ex_val;                                //Set excluded lattice points to this value
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
    p->num_g_x         = 0;
    p->num_g_y         = 0;
    p->cutoff          = 0;
    p->APS             = 0;
    p->radius          = 0;
    p->cell_size       = 0;
    p->bulk            = 0;
    p->leaflet         = 0;
    p->mem_size        = 0;
    p->box_x           = 0;
    p->box_y           = 0;
    p->box_z           = 0;
    p->print_stride    = 0;
    p->b_lf_pdb        = 0;
    p->num_lip_t1      = 0;
    p->num_lip_t2      = 0;
    p->num_lipids      = 0;
    p->b_lf_param      = 0;
    p->ex_val          = 0.0;

    //here we set the program description
    p->program_description = p->program_description + "3D Enrichment is an analysis program that measures the enrichment factor for a given lipid type and ";
    p->program_description = p->program_description + "projects this data onto a 3-dimensional lattice. ";
}


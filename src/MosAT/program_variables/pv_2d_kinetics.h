
//define up to 4-dimensional vector of doubles
typedef vector<double>  dv1d;
typedef vector<dv1d>    dv2d;
typedef vector<dv2d>    dv3d;
typedef vector<dv3d>    dv4d;

//define up to 4-dimensional vector of ints
typedef vector<int>     iv1d;
typedef vector<iv1d>    iv2d;
typedef vector<iv2d>    iv3d;
typedef vector<iv3d>    iv4d;

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

    //2d_kinetics stuff
    string param_file_name;                       //Name of the parameter file with voronoi lipids and target atoms for density
    string k_file_name;                           //Name of the k_off output file
    string lf_pdb_file_name;                      //Name of the lf_pdb_file
    string leaflet_finder_param_name;             //Name of the leaflet finder param file
    string pf_pdb_file_name;                      //Name of the output pdb file with the protein indicated by B-factor
    string protein_finder_param_name;             //Name of the protein finder param file
    int b_lf_param;                               //Tells if the user included a lipid types parameter file
    int num_g_x;                                  //This is the number of grid points in x direction
    int num_g_y;                                  //This is the number of grid points in y direction
    int leaflet;                                  //Which leaflet should analysis be performed for (0:neither 1:upper 2:lower)
    int out_data_format;                          //What format should the output data be in? (0:matrix 1:vector)
    int b_lf_pdb;                                 //Print pdb file with marked leaflets by beta value
    int num_lipids;                               //How many lipids are in the target leaflet(s)
    int range;                                    //How many frames on each side of t should a running average be computed over?
    int b_pf_param;                               //Tells if the user included a protein finder parameter file
    int b_pf_pdb;                                 //Print pdb file with marked protein by beta value
    int v_prot;                                   //Include protein in voronoi diagram?
    int b_clean;                                  //Remove voronoi diagrams after computing dwell times
    int dump;                                     //Dump all bound lipids at the end of the trajectory? i.e. record how long they have been bound?
    double APS;                                   //This is the area of a grid square
    double radius;                                //Radius of the atom
    double cell_size;                             //This is the lengh between grid points
    double box_x;                                 //user input box (x-dimension)
    double box_y;                                 //user input box (y-dimension)
    double delta_t;                               //User specified time step (overwrites from trajectory)
    double voro_stamp_rad;                        //The stamping radius used to find candidates for the voronoi diagram
    double c_dist;                                //Distance cutoff for counting protein atoms in voronoi diagram
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
    p->num_g_x                 = 0;
    p->num_g_y                 = 0;
    p->leaflet                 = 0;
    p->out_data_format         = 0;
    p->APS                     = 0;
    p->radius                  = 0;
    p->cell_size               = 0;
    p->box_x                   = 0;
    p->box_y                   = 0;
    p->b_lf_pdb                = 0;
    p->num_lipids              = 0;
    p->range                   = 0;
    p->delta_t                 = 0; 
    p->b_lf_param              = 0;
    p->voro_stamp_rad          = 0.85;
    p->c_dist                  = 0.6;
    p->b_pf_param              = 0;
    p->b_pf_pdb                = 0;
    p->v_prot                  = 0;
    p->b_clean                 = 0;
    p->dump                    = 0;

    //here we set the program description
    p->program_description = "2D Kinetics is an analysis program designed for characterizing the lipid residence time and projecting this data onto ";
    p->program_description = p->program_description + "the XY plane. The main output from 2D Kinetics is a set of binding events files; one for each lattice point.";
}


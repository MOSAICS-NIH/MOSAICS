
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

    //apl stuff
    string apl_file_name;                         //Name of the area per lipid output file
    string param1_file_name;                      //Name of the parameter file with lipid types for apl projection
    string param2_file_name;                      //Name of the parameter file with lipid types for voronoi diagram
    string lf_pdb_file_name;                      //Name of the output pdb file with leaflets indicated by B-factor
    string leaflet_finder_param_name;             //Name of the leaflet finder param file
    string pf_pdb_file_name;                      //Name of the output pdb file with the protein indicated by B-factor
    string protein_finder_param_name;             //Name of the protein finder param file
    string traj_voro_file_name;                   //Name of the second trajectory file used for constructing voronoi diagrams
    int b_lf_param;                               //Tells if the user included a lipid types parameter file
    int b_pf_param;                               //Tells if the user included a protein finder parameter file
    int leaflet;                                  //Which leaflet? (0:both 1:upper 2:lower)
    int out_data_format;                          //What format should the output data be in? (0:matrix 1:vector)
    int b_stdev;                                  //Compute the z-coord standard deviation?
    int b_clean;                                  //Delete single frame rmsd files after computing average?
    int b_lf_pdb;                                 //Print pdb file with marked leaflets by beta value
    int b_pf_pdb;                                 //Print pdb file with marked protein by beta value
    int b_voronoi;                                //Write the single frame voronoi diagrams to file?
    int v_prot;                                   //Include protein in voronoi diagram?
    double APS;                                   //This is the area of a grid square
    double radius;                                //Radius of the atom
    double cutoff;                                //Percentage of average density used for excluding data
    double box_x;                                 //Grid x dimension
    double box_y;                                 //Grid y dimension
    double box_z;                                 //Grid z dimension
    double c_dist;                                //Distance cutoff for counting protein atoms in voronoi diagram
    double voro_stamp_rad;                        //The stamping radius used to find candidates for the voronoi diagram
    double bin_width;                             //Bin width for the histogram
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

    //epl stuff
    p->APS                      = 0;
    p->radius                   = 0;
    p->cutoff                   = 0;
    p->leaflet                  = 0;
    p->out_data_format          = 0;
    p->box_x                    = 0;
    p->box_y                    = 0;
    p->box_z                    = 0;
    p->b_stdev                  = 0;
    p->b_clean                  = 0;
    p->b_lf_pdb                 = 0;
    p->b_pf_pdb                 = 0;
    p->b_lf_param               = 0;
    p->b_pf_param               = 0;
    p->c_dist                   = 0;
    p->b_voronoi                = 0;
    p->voro_stamp_rad           = 0.85;
    p->v_prot                   = 1.0; 
    p->bin_width                = 0.1;
    p->ex_val                   = 0.0;

    //here we set the program description
    p->program_description = p->program_description + "APL 3d is an analysis tools used for computing the area per lipid based on the lipids position in XYZ. ";
}


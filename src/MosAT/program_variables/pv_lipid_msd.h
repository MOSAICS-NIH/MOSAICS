
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

    //lipid msd stuff
    string param_file_name;                       //Name of the parameter file
    string lf_pdb_file_name;                      //Name of the lf_pdb_file
    string leaflet_finder_param_name;             //Name of the leaflet finder param file
    string msd_file_name;                         //Name of the ouput file with MSD data
    int b_lf_param;                               //Tells if the user included a lipid types parameter file
    int leaflet;                                  //Which leaflet should analysis be performed for (0:neither 1:upper 2:lower)
    int b_lf_pdb;                                 //Print pdb file with marked leaflets by beta value
    int num_lipids;                               //Number of lipids of the target type in the target leaflet
    double ef_dt;                                 //The effective time step between frames
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

    p->leaflet                = 0;
    p->b_lf_pdb               = 0;
    p->b_lf_param             = 0;
    p->num_lipids             = 0;
    p->ef_dt                  = 0.0;

    //here we set the program description
    p->program_description = p->program_description + "Lipid MSD is an analysis tool designed to compute the mean square displacement for lipids in the xy plane as a function of time. This should yield a straight line where the slope is proportional to the diffusion coeficient."; 
}


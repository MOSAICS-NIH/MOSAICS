
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

    //Contact Analysis stuff
    string group_1_file_name;                     //Name of the group 1 atoms index file
    string group_2_file_name;                     //Name of the group 2 atoms index file
    string cont_file_name;                        //Name of the output file with contacts
    double cdist;                                 //Distance threshold (nm) for counting contacts
    double max_dash_rad;                          //dash thickness for contacts with frequency of 1
    double min_dash_rad;                          //dash thickness for contacts with frequency of 0
    int size_x;                                   //x-dimension of contact profile
    int size_y;                                   //y-dimenstion of contact profile
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function initializes the variables held in program_variables                                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void initialize_program_variables(program_variables *p)
{
    //initialize program variables here

    p->stride       = 1;
    p->start_frame  = 0;
    p->end_frame    = -1;
    p->b_print      = 0;
    p->b_lsq        = 0;
    p->lsq_dim      = 3;
    p->lsq_ref      = 0;
    p->cdist        = 1.0;
    p->size_x       = 0;
    p->size_y       = 0;
    p->max_dash_rad = 0.3;
    p->min_dash_rad = 0.0;

    //here we set the program description
    p->program_description = p->program_description + "Contact Analysis is an analysis tool used for characterizing contacts between two groups of atoms.";	    
}


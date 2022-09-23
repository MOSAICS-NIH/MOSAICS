
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

    //here we set the program description
    p->program_description = p->program_description + "MOSAICS Analysis Template (MosAT) is a versatile C++ program used for ";
    p->program_description = p->program_description + "reading trajectory data. The program is designed to ";
    p->program_description = p->program_description + "make easy the task of writing an analysis program without the sacrifice of performance. ";
    p->program_description = p->program_description + "Specifically, MosAT will read a trajectory without any additional programming required ";
    p->program_description = p->program_description + "by the user. MosAT also reads a reference file such as a .gro or .pdb file. Using the ";
    p->program_description = p->program_description + "reference file, MosAT generates the atom and residue names and numbers for all molecules ";
    p->program_description = p->program_description + "in the system. Additional data such as atomic masses, the initial coordinates and initial ";
    p->program_description = p->program_description + "box are also extracted. This makes possible the creation of data structures before the ";
    p->program_description = p->program_description + "main loop is executed at which point the trajectory is read. Furthermore, to reduce memory ";
    p->program_description = p->program_description + "constraints, only a single frame is ever read into memory at a given time. MosAT makes easy ";
    p->program_description = p->program_description + "the addition of command line arguments. This is done using the add_argument function. ";
    p->program_description = p->program_description + "Likewise, user defined atom selections are possible through use of selection cards. To ensure speed, ";
    p->program_description = p->program_description + "MosAT is pre-parallelized (using MPI) to take advantage of modern supercomputers. In ";
    p->program_description = p->program_description + "general, the block parallelization scheme used by MosAT splits the trajectory into chunks ";
    p->program_description = p->program_description + "where each core is assigned a section to read. In most cases, this removes any communication ";
    p->program_description = p->program_description + "requirements until after the trajectory has been read. For this reason, the scalability ";
    p->program_description = p->program_description + "with block parallelization is nearly linear. Alternatively, the parallelization scheme can ";
    p->program_description = p->program_description + "be set (by turning off block parallelization) so each rank reads the entire trajectory. ";
    p->program_description = p->program_description + "To ensure reliability and speed of an analysis program, MosAT makes full use of the trajectory file ";
    p->program_description = p->program_description + "readers developed by the Gromacs team and simplified for use with analysis tools by the MDTraj developers. ";
    p->program_description = p->program_description + "Additional functionality implemented in MosAT include least squares ";
    p->program_description = p->program_description + "fitting, a stride for skipping frames, and control over the begin and end frames read. ";
    p->program_description = p->program_description + "Compatible trajectory types include xtc, trr, gro and pdb with cross IO functionality supported.";
}


#ifndef PV_LIPID_NEIGHBOR_CONTENT_H
#define PV_LIPID_NEIGHBOR_CONTENT_H

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This structure holds the trajectory controls and Lipid Neighbor Content analysis settings.                //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct program_variables
{
    string program_description;                   //Description printed with the native help output
    string in_file_name;                          //Input coordinate trajectory
    string ref_file_name;                         //Reference structure for topology and leaflet assignment
    int    stride;                                //Read every stride-th original trajectory frame
    int    start_frame;                           //First permitted original trajectory frame
    int    end_frame;                             //Last permitted original trajectory frame; -1 means all

    string param_file_name;                       //Networked lipid-type, cutoff and atom-card selection
    string lnc_file_name;                         //Block-average output and base name for related files
    string lf_pdb_file_name;                      //Optional PDB showing native leaflet assignments
    string leaflet_finder_param_name;             //Optional additional leaflet-finder definitions
    int    b_lf_pdb;                              //Write the leaflet-assignment PDB
    int    b_lf_param;                             //Read additional leaflet-finder definitions
    int    leaflet;                               //Target: 0 full membrane, 1 upper or 2 lower
    int    block;                                 //Number of retained central frames in an averaging block
    int    num_lipids;                            //Number of selected lipids in the target leaflet(s)
    int    range;                                 //Nonnegative half-width of the centered binary noise filter
    int    window;                                //Complete centered window size, 2*range + 1
    int    read_first;                            //First effective frame to read, including left context
    int    read_last;                             //Last effective frame to read, including right context
    int    empty_count;                           //How many lipids had no neighbors
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function initializes the native trajectory defaults and describes the new analysis.                  //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void initialize_program_variables(program_variables *p)
{
    p->stride        = 1;
    p->start_frame   = 0;
    p->end_frame     = -1;

    p->b_lf_pdb      = 0;
    p->b_lf_param    = 0;
    p->leaflet       = 0;
    p->block         = 0;
    p->num_lipids    = 0;
    p->range         = 0;
    p->window        = 1;
    p->read_first    = -1;
    p->read_last     = -1;
    p->empty_count   = 0;

    p->program_description = "Lipid Neighbor Content measures filtered neighbor-type composition around each central lipid type. ";
    p->program_description += "The -crd card has three columns: residue name, central-type cutoff (nm), ";
    p->program_description += "and a linked one-column atom-name card. Use 1.3 nm for the default cutoff. ";
    p->program_description += "Linked paths are relative to the working directory. Each residue type ";
    p->program_description += "must occur once in the main card. Geometric centers use the selected atoms; ";
    p->program_description += "contacts use direct two-dimensional distances in XY strictly less than the cutoff, ";
    p->program_description += "with periodic minimum-image corrections. Prepare whole molecules first and avoid rotating the system. ";
    p->program_description += "The required -range specifies a nonnegative noise-filter half-width N in ";
    p->program_description += "effective frames after -b, -e and -stride. A raw pair contact is retained at ";
    p->program_description += "center t when present in at least half of the 2N+1 frames t-N through t+N. ";
    p->program_description += "Omega is fixed at 0.5; -range 0 reproduces the instantaneous contacts. ";
    p->program_description += "Only complete windows contribute: the first and last N effective centers ";
    p->program_description += "are omitted, and at least 2N+1 effective frames are required. Self pairs stay zero. ";
    p->program_description += "The program averages filtered n_ij, which tells if lipid j is a neighbor to i, ";
    p->program_description += "over lipids j thus giving the percentage of each lipid type in i's inner shell. ";
    p->program_description += "The data is then averaged over lipid i, to give the average composition for each ";
    p->program_description += "lipid type. The average composition for each lipid type is stored for each trajectory frame. ";
    p->program_description += "Rank 0 collects the averages and report raw data. The average composition data is ";
    p->program_description += "also averaged over blocks of the trajectory and this data is reported too. Last, the ";
    p->program_description += "change in composition between neighboring blocks is computed and summed over lipid type j ";
    p->program_description += "to give how much the composition changes for lipid type i between blocks. The program also ";
    p->program_description += "monitors for lipids that have no neighbors and statistics are reported for how many lipids ";
    p->program_description += "are found. Similalry, the average number of neighbors is reported for each lipid type along ";
    p->program_description += "with the variance over the total number of lipids analyzed. ";
}

#endif

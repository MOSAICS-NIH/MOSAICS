
#ifdef __APPLE__
#  define off64_t off_t
#  define fopen64 fopen
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This is a class for working with refference files (pdb and gro)                                          //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Ref
{
    private:
        string ref_file_name;                                                                //name of the ref file to read 
        int i             = 0;                                                               //standard variable used in loops
        int j             = 0;                                                               //standard variable used in loops
        int ref_f         = 0;                                                               //format of the ref file (0:gro, 1:pdb)
        int status        = 0;                                                               //tells if the ref file was found or not
        int world_size    = 0;                                                               //size of mpi world
        int world_rank    = 0;                                                               //rank of mpi process
        int box_dimension = 0;                                                               //how many dimensions in the box
        int frames        = 0;                                                               //the number of frames in the ref file
        int num_atoms     = 0;                                                               //how many atoms in the ref file
        int step          = 0;                                                               //the step
        gmx_bool bV       = 0;                                                               //does the ref file contain velocities
        gmx_bool bBox     = 0;                                                               //is a simulation box present
        real time         = 0.0;                                                             //the time stored in the ref

    public:
        iv1d atom_nr{};                                                                      //holds the atom numbers
        iv1d res_nr{};                                                                       //holds the res numbers
        sv1d atom_name{};                                                                    //holds the atom names 
        sv1d res_name{};                                                                     //holds the res names
        dv1d beta{};                                                                         //Beta factor for pdb files
        dv1d weight{};                                                                       //Weight used for pdb files
        sv1d element{};                                                                      //Element collumn in pdb file
        cv1d chain_id{};                                                                     //Chain id for pdb file
        char title[200];                                                                     //holds the title line
        rvec *r;                                                                             //holds the coordinates
        rvec *v;                                                                             //holds the velocities
        matrix box;                                                                          //holds the box dimensions

    public:
        void get_ref(string my_ref_file_name,string tag);                                    //read the refference file
        int  get_num_atoms();                                                                //returns the number of atoms
};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function reads the reference file (pdb or gro)                                                      //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Ref::get_ref(string my_ref_file_name,string tag)
{
    //get the world size and assign ranks
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);      //get the world size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);      //get the process rank

    ref_file_name = my_ref_file_name;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                          //
    // check the file extension is correct                                                                      //
    //                                                                                                          //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int length = ref_file_name.length();
    if(ref_file_name.at(length-4) == '.' && ref_file_name.at(length-3) == 'g' && ref_file_name.at(length-2) == 'r' && ref_file_name.at(length-1) == 'o') //gro
    {
        ref_f = 0;
    }
    else if(ref_file_name.at(length-4) == '.' && ref_file_name.at(length-3) == 'p' && ref_file_name.at(length-2) == 'd' && ref_file_name.at(length-1) == 'b') //pdb
    {
        ref_f = 1;
    }
    else
    {
        if(world_rank == 0)
        {
            printf("Supported reference file types are pdb and gro. Make sure your reference filename (%s) ends with one of these.\n",tag.c_str());
        }
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                          //
    // check the file exists                                                                                    //
    //                                                                                                          //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    FILE *test_file = fopen64(ref_file_name.c_str(), "r");
    if(test_file == NULL)
    {
        if(world_rank == 0)
        {
            printf("failure opening %s. Make sure the file exists. \n",ref_file_name.c_str());
        }
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }
    else
    {
        fclose(test_file);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                          //
    // open the ref file                                                                                        //
    //                                                                                                          //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    FILE *ref_file; 
    if(ref_f == 0 || ref_f == 1) //gro or pdb
    {
        ref_file = fopen64(ref_file_name.c_str(), "r");
        if(ref_file == NULL)
        {
            if(world_rank == 0)
            {
                printf("failure opening %s. Make sure the file exists. \n",ref_file_name.c_str());
            }
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }
    }

    if(ref_f == 0) //gro
    {
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                          //
        // analyze the ref file                                                                                     //
        //                                                                                                          //
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        analyze_gro_file_ref(&ref_file,&num_atoms,&frames,world_rank,box,&box_dimension);

        if(frames > 1)
        {
            if(world_rank == 0)
            {
                printf("Only single frame ref files are supported (%s). There were %d frames detected in %s. \n",tag.c_str(),frames,ref_file_name.c_str());
            }
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                          //
        // allocate memory to hold coords, etc.                                                                     //
        //                                                                                                          //
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        atom_nr.resize(num_atoms,0);
        res_nr.resize(num_atoms,0);
        atom_name.resize(num_atoms);
        res_name.resize(num_atoms);
        beta.resize(num_atoms,0);
        weight.resize(num_atoms,0);
        element.resize(num_atoms);
        chain_id.resize(num_atoms,'A');
        r = (rvec*)calloc(num_atoms , sizeof(rvec));
        v = (rvec*)calloc(num_atoms , sizeof(rvec));

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                          //
        // read the ref file                                                                                        //
        //                                                                                                          //
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        read_gro_frame_by_char(&ref_file,box,&num_atoms,atom_nr,res_nr,res_name,atom_name,r,v,title,world_rank,&time,&step,frames,&bV);
    }
    else if(ref_f == 1) //pdb
    {
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                          //
        // analyze the ref file                                                                                     //
        //                                                                                                          //
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        analyze_pdb_file_ref(&ref_file,&num_atoms,&frames,world_rank,ref_file_name);

        if(frames > 1)
        {
            if(world_rank == 0)
            {
                printf("Only single frame ref files are supported (%s). There were %d frames detected in %s. \n",tag.c_str(),frames,ref_file_name.c_str());
            }
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                          //
        // allocate memory to hold coords, etc.                                                                     //
        //                                                                                                          //
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        atom_nr.resize(num_atoms,0);
        res_nr.resize(num_atoms,0);
        atom_name.resize(num_atoms);
        res_name.resize(num_atoms);
        beta.resize(num_atoms,0);
        weight.resize(num_atoms,0);
        element.resize(num_atoms);
        chain_id.resize(num_atoms,'A');
        r = (rvec*)calloc(num_atoms , sizeof(rvec));
        v = (rvec*)calloc(num_atoms , sizeof(rvec));

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                          //
        // read the ref file                                                                                        //
        //                                                                                                          //
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        read_pdb_frame_by_char(&ref_file,box,atom_nr,res_nr,res_name,atom_name,r,title,world_rank,&time,&step,frames,beta,weight,element,chain_id,&bBox);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns the number of atoms                                                                //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Ref::get_num_atoms()
{
    return num_atoms;
}

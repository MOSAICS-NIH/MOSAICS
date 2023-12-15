
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This is a class for working with contact profiles between a lipid and the protein                        //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Contacts
{
    private:
        string base_file_name;                        //name used to generate other names for contact profiles
        string tmp_file_name;                         //name used to write out temporary contact profiles  
        string cont_file_name;                        //name of the complete contacts file
        int    world_rank;                            //rank of the core
        int    world_size;                            //size of mpi world
        int    prot_size;                             //how big is the protein
        int    size_y;                                //how big is the target resi
        vector <int64_t>  position{};                 //holds the starting position of each frame in the contacts profile
        FILE *cont_file;                              //File used to read the contacts file
        FILE *tmp_file;                               //File used to write temporary contacts files

    public:
        void init(string this_file_name,int this_size_y,int this_prot_size);   //initializes the contact matrix
        void prime_tmp_file();                                                 //opens the temporary file for writing
        void add_profile(iv2d &this_profile,int frame);                        //adds the contact profile for the current frame to temporary file
	void close_tmp_file();                                                 //closes the temporary contacts file
	void merge_profiles();                                                 //merges temporary contact files and analyzes the resulting file 
	iv2d get_profile(int current_frame);                                   //read a frame from the contact matrix
        iv2d get_profile_alt(int current_frame,iv2d &contact_profile);         //read a frame from the contact matrix. return record on contacts set
        void clear_profile(iv2d &history,iv2d &contact_profile);               //reset the current contact matrix
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function sets the base file name and and the matrix dimensions. Also sets up mpi variables.          //    
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Contacts::init(string this_file_name,int this_size_y,int this_prot_size)
{
    //get the world size and assign ranks
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);      //get the world size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);      //get the process rank

    base_file_name = this_file_name; 
    size_y         = this_size_y;
    prot_size      = this_prot_size; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function opens temporary contact files for writing                                                   //    
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Contacts::prime_tmp_file()
{
    string tag    = "_" + to_string(world_rank) + "_contacts.dat";
    tmp_file_name = chop_and_add_tag(base_file_name,tag);

    tmp_file             = fopen(tmp_file_name.c_str(), "w");
    if(tmp_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",tmp_file_name.c_str());
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function add the current profile to the tmp file                                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Contacts::add_profile(iv2d &this_profile,int frame)
{
    int i = 0;       //standard variable used in loops

    fprintf(tmp_file,"%s %d \n","frame",frame);

    for(i=0; i<this_profile.size(); i++) //loop over contacts
    {
        fprintf(tmp_file,"%d %d \n",this_profile[i][0],this_profile[i][1]);
    }	
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function closes temporary contact files after writing                                                //    
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Contacts::close_tmp_file()
{
    fclose(tmp_file);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function merges the temporary contact files and then gets the positon of each frame in the resulting //    
// file.                                                                                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Contacts::merge_profiles()
{
    int i = 0;               //standard variable used in loops
    int j = 0;               //standard variable used in loops

    string tag     = "_contacts.dat";
    cont_file_name = chop_and_add_tag(base_file_name,tag);

    if(world_rank == 0)
    {
        printf("\nSplicing together temporary contact profiles \n");

        //open file and copy tmp files into it
        cont_file = fopen(cont_file_name.c_str(), "w");
        if(cont_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",cont_file_name.c_str());
        }
        else
        {
            for(i=0; i<world_size; i++) //loop over mpi cores
            {
                //generate name of the temp contacts file
                string this_tag = "_" + to_string(i) + "_contacts.dat";
                tmp_file_name   = chop_and_add_tag(base_file_name,this_tag);

                //make an index and read in contacts data
                Index tmp_contacts;
                tmp_contacts.get_index(tmp_file_name);

                //remove temp contacts file
                remove(tmp_file_name.c_str());

                for(j=0; j<tmp_contacts.index_s.size(); j+=2) //loop over items in contacts data file
                {
                    fprintf(cont_file,"%s %s \n",tmp_contacts.index_s[j].c_str(),tmp_contacts.index_s[j+1].c_str());
                }
            }
            fclose(cont_file);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if(world_rank == 0)
    {
        printf("\nAnalyzing contacts file %s \n\n",cont_file_name.c_str());
    }

    //analyze contacts file to get the position of each frame
    cont_file = fopen(cont_file_name.c_str(), "r");
    if(cont_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",cont_file_name.c_str());
    }
    else
    {
        char line[200];                //used to read in a line of data
        char my_string_1[20];          //stores the current string
        char my_string_2[20];          //stores the current string
        int result      = 0;           //tells if reading the next string was successful

        //loop over trajectory
        while(fgets(line, sizeof line, cont_file) != NULL)
        {
            int line_offset = 0;                    //stores the position in the current line

            result = next_string(world_rank,200,line,my_string_1,20,&line_offset);
            result = next_string(world_rank,200,line,my_string_2,20,&line_offset);

            if(strcmp(my_string_1, "frame") == 0)
            {
                position.push_back(ftell(cont_file));
            }
        }
        fclose(cont_file);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads in a frame from the contacts profile data                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
iv2d Contacts::get_profile(int current_frame)
{
    char line[200];           //used to read in a line of data
    char my_string_1[20];     //stores the current string
    char my_string_2[20];     //stores the current string
    int result      = 0;      //tells if reading the next string was successful

    iv2d contact_profile(size_y, iv1d(prot_size,0));

    int64_t offset_r = position[current_frame];

    FILE *in_file;
    in_file = fopen(cont_file_name.c_str(), "r");
    if(in_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",cont_file_name.c_str());
    }
    else
    {
        //move to the desired frame
        fseek(in_file,offset_r,SEEK_SET);

        //read the current frame
        while(fgets(line, sizeof line, in_file) != NULL)
        {
            int line_offset = 0;                    //stores the position in the current line

            result = next_string(world_rank,200,line,my_string_1,20,&line_offset);
            result = next_string(world_rank,200,line,my_string_2,20,&line_offset);

            if(strcmp(my_string_1, "frame") == 0)
            {
                break;
            }
            else
            {
                contact_profile[atoi(my_string_1)][atoi(my_string_2)] = 1;
            }
        }
        fclose(in_file);
    }

    return contact_profile;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads in a frame from the contacts profile data and returns a record of the contacts        //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
iv2d Contacts::get_profile_alt(int current_frame,iv2d &contact_profile)
{
    char line[200];           //used to read in a line of data
    char my_string_1[20];     //stores the current string
    char my_string_2[20];     //stores the current string
    int result      = 0;      //tells if reading the next string was successful

    iv2d history(0, iv1d(2,0));      //store a record of which contacts are present

    int64_t offset_r = position[current_frame];

    FILE *in_file;
    in_file = fopen(cont_file_name.c_str(), "r");
    if(in_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",cont_file_name.c_str());
    }
    else
    {
        //move to the desired frame
        fseek(in_file,offset_r,SEEK_SET);

        //read the current frame
        while(fgets(line, sizeof line, in_file) != NULL)
        {
            int line_offset = 0;                    //stores the position in the current line

            result = next_string(world_rank,200,line,my_string_1,20,&line_offset);
            result = next_string(world_rank,200,line,my_string_2,20,&line_offset);

            if(strcmp(my_string_1, "frame") == 0)
            {
                break;
            }
            else
            {
                //get the atom number of the lipid and protein atoms
                int lipid_atom = atoi(my_string_1);
                int prot_atom  = atoi(my_string_2);

                //set the contact in the matrix
                contact_profile[lipid_atom][prot_atom] = 1;

                //store a record of the contact
                iv1d this_contact(2,0);    
                this_contact[0] = lipid_atom;
                this_contact[1] = prot_atom;
                history.push_back(this_contact); 
            }
        }
        fclose(in_file);
    }

    return history;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function clears the current contacs matrix                                                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Contacts::clear_profile(iv2d &history,iv2d &contact_profile)
{
    int i=0;    //standard variable used in loops

    for(i=0; i<history.size(); i++) //loop over contacts
    {
        int lipid_atom = history[i][0];
        int prot_atom  = history[i][1];

        contact_profile[lipid_atom][prot_atom] = 0;
    }
}

//define up to 4-dimensional vector of strings
typedef vector<string>  sv1d;
typedef vector<sv1d>    sv2d;
typedef vector<sv2d>    sv3d;
typedef vector<sv3d>    sv4d;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function determines which atoms belong to the protein and which do not.                             //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void protein_finder(vector <int> &b_protein,vector <string> &res_name,int num_atoms,int b_pf_param,string protein_finder_param_name)
{
    int i = 0;
    int j = 0;
    sv1d residue_types(0);

    residue_types.push_back("ALA"); 
    residue_types.push_back("ARG");
    residue_types.push_back("ASP");
    residue_types.push_back("ASN");
    residue_types.push_back("CYS");
    residue_types.push_back("GLU");
    residue_types.push_back("GLN");
    residue_types.push_back("GLY");
    residue_types.push_back("HIS");
    residue_types.push_back("ILE");
    residue_types.push_back("LEU");
    residue_types.push_back("LYS");
    residue_types.push_back("MET");
    residue_types.push_back("PHE");
    residue_types.push_back("PRO");
    residue_types.push_back("SER");
    residue_types.push_back("THR");
    residue_types.push_back("TRP");
    residue_types.push_back("TYR");
    residue_types.push_back("VAL");
    residue_types.push_back("HIH");
    residue_types.push_back("GLH");
    residue_types.push_back("ASH");
    residue_types.push_back("HSP");
    residue_types.push_back("HSE");
    residue_types.push_back("HSD");

    //check if user provided additional parameters
    if(b_pf_param == 1)
    {
        Index protein_finder_param;
        protein_finder_param.get_index(protein_finder_param_name);

        //add new parameters to lip_parameters
        for(i=0; i<protein_finder_param.index_s.size(); i++)
        {
            residue_types.push_back(protein_finder_param.index_s[i]);
        }
    }

    //identify protein residues
    for(i=0; i<num_atoms; i++)
    {
        b_protein[i] = 0;

        for(j=0; j<residue_types.size(); j++)
        {
            if(strcmp(res_name[i].c_str(), residue_types[j].c_str()) == 0) //atom belongs to the protein
            {
                b_protein[i] = 1;
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function counts how many atoms are in the protein.                                                  //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int count_protein_atoms(vector <int> &b_protein,int num_atoms)
{
    int i = 0;
    int count = 0;

    for(i=0; i<num_atoms; i++)
    {
        if(b_protein[i] == 1)
        {
            count++;
        }
    }
    return count;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function fills an index with the protein atoms.                                                     //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_prot(vector <int> &b_protein,int num_atoms,vector <int> &prot)
{
    int i = 0;
    int count = 0;

    for(i=0; i<num_atoms; i++)
    {
        if(b_protein[i] == 1)
        {
            prot[count] = i + 1;
            count++;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function sets the beta value according to whether an atom belongs to the protein or not.            //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void set_beta_pf(vector <double> &beta,vector <int> &b_protein,int num_atoms)
{
    int i = 0;

    for(i=0; i<num_atoms; i++)
    {
        if(b_protein[i] == 1) //protein atom
        {
            beta[i] = 1;
        }
        else //not protein atom
        {
            beta[i] = 0;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes a pdb with the protein indicated by beta value.                                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void dump_protein_pdb(int world_rank,int print_prot_pdb,string pf_pdb_file_name,vector <double> &beta,vector <int> &b_protein,int num_atoms,
                      matrix box_ref,vector <int> &atom_nr,vector <int> &res_nr,vector <string> &res_name,vector <string> &atom_name,rvec *r_ref,
                      char title[200],vector <double> &weight,vector <string> &element,vector <char> &chain_id)
{
    //print a pdb with distinguished protein
    if(print_prot_pdb == 1 && world_rank == 0)
    {
        FILE *pf_pdb_file = fopen(pf_pdb_file_name.c_str(), "w");
        if(pf_pdb_file == NULL)
        {
            printf("failure opening %s (pf). Make sure the file exists. \n",pf_pdb_file_name.c_str());
        }
        set_beta_pf(beta,b_protein,num_atoms);
        write_frame_pdb(box_ref,num_atoms,atom_nr,res_nr,res_name,atom_name,r_ref,title,world_rank,&pf_pdb_file,beta,weight,element,chain_id,0);
        fclose(pf_pdb_file);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints information about the protein (how many residues,how many atoms).                    //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void print_protein_info(vector <int> &prot,vector <int> &res_nr,vector <string> &res_name,int world_rank)
{
    if(world_rank == 0)
    {
        vector <string> res_types;

        int i = 0;
        int j = 0;
        int total_res = 0;
        int prev_res = -1;
        int count = 0;

        //count protein residues
        for(i=0; i<prot.size(); i++)
        {
            if(res_nr[prot[i]-1] != prev_res)
            {
                total_res++;
                prev_res = res_nr[prot[i]-1];
            }
        }

        //count residue types
        for(i=0; i<prot.size(); i++)
        {
            int recorded = 0;
            int size = res_types.size();

            for(j=0; j<size; j++)
            {
                if(strcmp(res_name[prot[i]-1].c_str(), res_types[j].c_str()) == 0) //res type already recorder
                {
                    recorded = 1;
                }
            }

            if(recorded == 0)
            {
                res_types.push_back(res_name[prot[i]-1]);
            }
        }

        printf("There are %6d atoms, %9d residues and  %2d residue types in the protein as summarized below: \n",prot.size(),total_res,res_types.size());
        printf(" %2s %10s %11s \n","index","type","count");
        printf("-%5s-%10s-%11s\n","-----","----------","-----------");

        //count residues for each type
        for(i=0; i<res_types.size(); i++) //loop over residue types
        {
            int prev_res = -1;
            int count = 0;

            for(j=0; j<prot.size(); j++) //loop over protein selection
            {
                if(strcmp(res_name[prot[j]-1].c_str(), res_types[i].c_str()) == 0) //residue is of correct type
                {
                    if(res_nr[prot[j]-1] != prev_res) //new residue
                    {
                        count++;
                        prev_res = res_nr[prot[j]-1];
                    }
                }
            }

            printf(" %5d %10s %11d \n",i,res_types[i].c_str(),count);
        }
        printf("\n");
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the number of residues making the protein                                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int count_protein_residues(vector <int> &res_nr,vector <int> &prot)
{
    int i = 0;
    int j = 0;
    int prev_res_nr = -1;
    int count = 0;

    for(i=0; i<prot.size(); i++) //loop over protein atoms
    {
        if(res_nr[prot[i]-1] != prev_res_nr)
        {
            prev_res_nr = res_nr[prot[i]-1];
            count++;
        }
    }
    return count;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function runs protein finder and makes a vector with the protein atoms                               //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::get_protein(string protein_finder_param_name,int b_pf_param)
{

    //resize b_protein
    b_protein.resize(num_atoms,0);

    //determine which atoms belong to the protein
    protein_finder(b_protein,res_name,num_atoms,b_pf_param,protein_finder_param_name);

    //count protein atoms
    int prot_size = count_protein_atoms(b_protein,num_atoms);

    //resize prot
    prot.resize(prot_size);

    //get the protein atoms
    get_prot(b_protein,num_atoms,prot);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes a pdb with the protein indicated by beta factor                                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::write_protein(string pf_pdb_file_name,int b_pf_pdb)
{
    dump_protein_pdb(world_rank,b_pf_pdb,pf_pdb_file_name,beta,b_protein,num_atoms,box_ref,atom_nr,res_nr,res_name,
                     atom_name,r_ref,title,weight,element,chain_id);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints information about the protein                                                        //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::get_prot_stats()
{
    print_protein_info(prot,res_nr,res_name,world_rank);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the number of protein residues                                                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::get_num_res_prot()
{
    return count_protein_residues(res_nr,prot);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the last atom of the current residue so after i++ we have the next residue          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::next_prot_res(int i)
{
    return jump_to_res_end(i,prot.size(),res_end,res_nr,prot,res_name,world_rank,res_start);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the first atom of the current residue                                               //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::p_res_start(int i)
{
    return res_start[res_nr[prot[i]-1]-1];
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the last atom of the current residue                                                //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::p_res_end(int i)
{
    return res_end[res_nr[prot[i]-1]-1];
}


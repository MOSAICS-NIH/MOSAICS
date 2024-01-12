
//define up to 4-dimensional vector of strings
typedef vector<string>  sv1d;
typedef vector<sv1d>    sv2d;
typedef vector<sv2d>    sv3d;
typedef vector<sv3d>    sv4d;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function takes an index holding atom types for a head and tail atom and determines wich lipids      //
// belong to the upper and lower leaflets.                                                                  //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void leaflet_finder(int num_atoms,vector <int> &atom_nr,vector <int> &res_nr,vector <string> &res_name,vector <string> &atom_name,rvec *r,
                    int world_rank,vector <int> &leaflets,vector <int> &res_start,vector <int> &res_end,int b_lf_param,string leaflet_finder_param_name)
{
    int i = 0;                           //standard variable used in loops
    int j = 0;                           //standard variable used in loops
    int k = 0;                           //standard variable used in loops
    int l = 0;                           //standard variable used in loops
    int m = 0;                           //standard variable used in loops
    int n = 0;                           //standard variable used in loops
    int current_res = -1;                //The current residue number
    double dz_1;                         //z-distance for tail 1
    double dz_2;                         //z-distance for tail 2
    double dz = 0;                       //z-distance which is finally used
    sv2d lip_parameters(0,sv1d(5));      //vector holding leaflet finder parameters

                              //type  |  h1  |  t1  |  h2  |  t2
    lip_parameters.push_back({ "POPE", "PO4", "C4A", "PO4", "C4B"  });        /*Martini POPE*/
    lip_parameters.push_back({ "POPG", "PO4", "C4A", "PO4", "C4B"  });        /*Martini POPG*/
    lip_parameters.push_back({ "POPC", "PO4", "C4A", "PO4", "C4B"  });        /*Martini POPC*/
    lip_parameters.push_back({ "DLPE", "PO4", "C3A", "PO4", "C3B"  });        /*Martini DLPE*/
    lip_parameters.push_back({ "DLPG", "PO4", "C3A", "PO4", "C3B"  });        /*Martini DLPG*/
    lip_parameters.push_back({ "DLPC", "PO4", "C3A", "PO4", "C3B"  });        /*Martini DLPC*/
    lip_parameters.push_back({ "DFF2", "PO4", "C4A", "PO4", "C4B"  });        /*Martini PIP2*/
    lip_parameters.push_back({ "PVP2", "PO4", "C4A", "PO4", "C4B"  });        /*Martini PIP2*/
    lip_parameters.push_back({ "POP2", "PO4", "C4A", "PO4", "C4B"  });        /*Martini PIP2*/
    lip_parameters.push_back({ "CHOL", "ROH", "C2",  "ROH", "C2"   });        /*Martini Cholesterol*/
    lip_parameters.push_back({ "POPC", "P"  , "C316","P"  , "C218" });        /*All atom POPC*/
    lip_parameters.push_back({  "DLPC", "P"  , "C312","P"  , "C212"});        /*All atom DLPC*/

    //check if user provided additional parameters
    if(b_lf_param == 1)
    {
        Index leaflet_finder_param;
        leaflet_finder_param.get_index(leaflet_finder_param_name);

        if(leaflet_finder_param.index_s.size()%5 != 0)
        {
            if(world_rank == 0)
            {
                printf("Leaflet finder parameter file does not contain a multiple of 5 items. Please check the formatting of %s. \n",leaflet_finder_param_name.c_str());
            }
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }
        else
        {
            int num_new_lipids = leaflet_finder_param.index_s.size()/5;

            //add new parameters to lip_parameters
            sv1d new_lipid(5);
            for(i=0; i<num_new_lipids; i++)
            {
                for(j=0; j<5; j++)
                {
                    new_lipid[j] = leaflet_finder_param.index_s[i*5 + j];

                    if(j == 4) //complete lipid stored in new_lipid
                    {
                        lip_parameters.push_back(new_lipid);
                    }
                }
            }
        }
    }

    //now we find the leaflets
    for(i=0; i<num_atoms; i++) //loop over system atoms
    {
        for(j=0; j</*bank_size*/ lip_parameters.size(); j++) //loop over lipid types
        {
            if(strcmp(res_name[i].c_str(), lip_parameters[j][0].c_str()) == 0) //lipid type is correct
            {
                int match_1 = 0;                       //did the first tail atoms match lip_parameters[j]
                int match_2 = 0;                       //did the second tail atoms match lip_parameters[j]
                double dz_1 = 0;                       //z-distance for tail 1
                double dz_2 = 0;                       //z-distance for tail 2
                double dz = 0;                         //z-distance which is finally used
                int min = res_start[res_nr[i]-1];      //first atom of this residue
                int max = res_end[res_nr[i]-1];        //last atom of this residue + 1

                //check the first tail atoms
                if(strcmp(atom_name[i].c_str(), lip_parameters[j][1].c_str()) == 0) //atom (i) is the first head atom
                {
                    for(k=min; k<=max; k++) //loop over residue atoms
                    {
                        if(strcmp(atom_name[k].c_str(), lip_parameters[j][2].c_str()) == 0) //atom (k) is the tail atom
                        {
                            dz_1 = r[i][2] - r[k][2];
                            match_1 = 1;
                        }
                    }

                    //check the second tail atoms
                    for(k=min; k<=max; k++) //loop over residue atoms
                    {
                        if(strcmp(atom_name[k].c_str(), lip_parameters[j][3].c_str()) == 0) //atom (k) is the second head atom
                        {
                            for(l=min; l<=max; l++) //loop over residue atoms
                            {
                                if(strcmp(atom_name[l].c_str(), lip_parameters[j][4].c_str()) == 0) //atom (l) is the second tail atom
                                {
                                    dz_2 = r[k][2] - r[l][2];
                                    match_2 = 1;
                                }
                            }
                        }
                    }

                    if(match_1 == 1 && match_2 == 1)
                    {
                        //should have now 2 distances
                        if(fabs(dz_1) > fabs(dz_2))
                        {
                            dz = dz_1;
                        }
                        else
                        {
                            dz = dz_2;
                        }

                        for(n=min; n<=max; n++) //loop over the residue atoms
                        {
                            if(dz > 0) //upper leaflet
                            {
                                leaflets[n] = 1;
                            }
                            else //lower leaflet
                            {
                                leaflets[n] = 2;
                            }
                        }
                    }
                }
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function sets the beta value according to whether an atom belongs to the upper leaflet, the lower   //
// leaflet or neither.                                                                                      //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void set_beta_leaflet_finder(vector <double> &beta,vector <int> &leaflets,int num_atoms)
{
    int i = 0;

    for(i=0; i<num_atoms; i++)
    {
        if(leaflets[i] == 0) //not a lipid atom
        {
            beta[i] = 0;
        }
        else if(leaflets[i] == 1) //upper leaflet
        {
            beta[i] = 1;
        }
        else if(leaflets[i] == 2) //lower leaflet
        {
            beta[i] = -1;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function counts how many atoms are in the selected leaflet (upper, lower or both).                  //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int count_leaflet_atoms(int leaflet,vector <int> &leaflets,int num_atoms)
{
    int i = 0;
    int count = 0;

    for(i=0; i<num_atoms; i++)
    {
        if(leaflet != 0 && leaflets[i] == leaflet)
        {
            count++;
        }
        else if(leaflet == 0 && (leaflets[i] == 1 || leaflets[i] == 2) )
        {
            count++;
        }
    }
    return count;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function fills an index with the selected lipid atoms (lower, upper or both).                       //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_membrane(int leaflet,vector <int> &leaflets,int num_atoms,vector <int> &mem)
{
    int i = 0;
    int count = 0;

    for(i=0; i<num_atoms; i++)
    {
        if(leaflet != 0 && leaflets[i] == leaflet)
        {
            mem[count] = i + 1;
            count++;
        }
        else if(leaflet == 0 && (leaflets[i] == 1 || leaflets[i] == 2) )
        {
            mem[count] = i + 1;
            count++;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function gets the opposing leaflet atoms.                                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int get_opposing_leaflet(int leaflet)
{
    int opposing = 0;

    if(leaflet == 1) //upper
    {
        opposing = 2;
    }
    else if(leaflet == 2) //lower
    {
        opposing = 1;
    }
    else if(leaflet == 0) //whole
    {
        opposing = 0;
    }

    return opposing;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints information about the target leaflet (lipid types and how many).                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_lipid_types(int mem_size,vector <int> &mem,vector <int> &res_nr,vector <string> &res_name,int world_rank,int leaflet)
{
    if(world_rank == 0)
    {
        vector <string> lipid_types;

        int i = 0;
        int j = 0;

        for(i=0; i<mem_size; i++)
        {
            int recorded = 0;
            int size = lipid_types.size();

            for(j=0; j<size; j++)
            {
                if(strcmp(res_name[mem[i]-1].c_str(), lipid_types[j].c_str()) == 0) //lipid type already recorder
                {
                    recorded = 1;
                }
            }

            if(recorded == 0)
            {
                lipid_types.push_back(res_name[mem[i]-1]);
            }
        }

        if(leaflet == 0)
        {
            printf("There are %2d lipid types in the selected leaflet (%s) as summarized below: \n",lipid_types.size(),"both");
        }
        else if(leaflet == 1)
        {
            printf("There are %2d lipid types in the selected leaflet (%s) as summarized below: \n",lipid_types.size(),"upper");
        }
        else if(leaflet == 2)
        {
            printf("There are %2d lipid types in the selected leaflet (%s) as summarized below: \n",lipid_types.size(),"lower");
        }
        printf(" %2s %10s %11s \n","index","type","count");
        printf("-%5s-%10s-%11s\n","-----","----------","-----------");

        //count each lipid type
        for(i=0; i<lipid_types.size(); i++) //loop over lipid types
        {
            int prev_res = -1;
            int count = 0;

            for(j=0; j<mem_size; j++) //loop over membrane selection
            {
                if(strcmp(res_name[mem[j]-1].c_str(), lipid_types[i].c_str()) == 0) //lipid is of correct type
                {
                    if(res_nr[mem[j]-1] != prev_res) //new lipid
                    {
                        count++;
                        prev_res = res_nr[mem[j]-1];
                    }
                }
            }

            printf(" %5d %10s %11d \n",i,lipid_types[i].c_str(),count);
        }
        printf("\n");
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints info about the user specified lipid selection (lipid types and how many).            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void lip_sel_stats(int mem_size,vector <int> &mem,vector <int> &res_nr,vector <string> &res_name,int world_rank,int lip_t_capacity,vector <string> &lip_t,string tag)
{
    if(world_rank == 0)
    {
        int i = 0;
        int j = 0;

        printf("You have selected (%s) the following lipid types: \n",tag.c_str());
        printf(" %5s %10s %11s \n","index","type","count");
        printf("-%5s-%10s-%11s\n","-----","----------","-----------");

        //count lipids for each type
        for(i=0; i<lip_t_capacity; i++) //loop over lipid types
        {
            int prev_res = -1;
            int count = 0;

            for(j=0; j<mem_size; j++) //loop over membrane selection
            {
                if(strcmp(res_name[mem[j]-1].c_str(), lip_t[i].c_str()) == 0) //lipid is of correct type
                {
                    if(res_nr[mem[j]-1] != prev_res) //new lipid
                    {
                        count++;
                        prev_res = res_nr[mem[j]-1];
                    }
                }
            }
            printf(" %5d %10s %11d \n",i,lip_t[i].c_str(),count);
        }
        printf("\n");
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes a pdb with the leaflets indicated by beta value.                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void dump_leaflets_pdb(int world_rank,int b_lf_pdb,string lf_pdb_file_name,vector <double> &beta,vector <int> &leaflets,int num_atoms,
                       matrix box_ref,vector <int> &atom_nr,vector <int> &res_nr,vector <string> &res_name,vector <string> &atom_name,rvec *r_ref,
                       char title[200],vector <double> &weight,vector <string> &element,vector <char> &chain_id)
{
    FILE *lf_pdb_file;

    if(world_rank == 0 && b_lf_pdb == 1)
    {
        lf_pdb_file = fopen(lf_pdb_file_name.c_str(), "w");
        if(lf_pdb_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",lf_pdb_file_name.c_str());
        }
        set_beta_leaflet_finder(beta,leaflets,num_atoms);
        write_frame_pdb(box_ref,num_atoms,atom_nr,res_nr,res_name,atom_name,r_ref,title,world_rank,&lf_pdb_file,beta,weight,element,chain_id,0);
        fclose(lf_pdb_file);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the last atom index in mem[] for the current residue.                               //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int jump_to_res_end(int i, int mem_size,vector <int> &res_end,vector <int> &res_nr,vector <int> &mem,vector <string> &res_name,int world_rank,vector <int> &res_start)
{
    int new_i = 0;       //index which is set to the last atom of the resideu and eventually returned.
    int j=0;             //standard variable used in loops

    for(j=i; j<mem_size; j++) //loop over mem index
    {
        if(mem[j]-1 == res_end[res_nr[mem[i]-1]-1]) //atom is the last atom for current res
        {
            //set the index
            new_i = j;
            break;
        }
    }

    if(new_i < i)
    {
        if(world_rank == 0)
        {
            printf("New atom encountered with an index smaller than the previous when attemping to move to the next residue. This is often caused by an error in the reference file (-ref) where a new residue type is listed  \n");
            printf("but the residue number is not changed from the previous. This causes errors because the leaflet, protein, and solvent finders are based on the residue name whereas the end of residues are determined by the residue number. \n"); 
	    printf("Thus, the error is encountered in programs that use leaflet, protein, or solvent finders. Most often, this is an error caused by the protein finder when the ends of the protein are capped. This can happen if the caps are \n"); 
	    printf("given a different name from the amino acid but the residue number is not increased. Please check your reference file for changes in residue name but not number. If this is the case, you can fix the error by giving each cap \n"); 
	    printf("a unique residue number (cannot be the same number as the residue coming befor or after the cap). Alternatively, you can add the missing residue, likely the cap, to the leaflet, protein, or solvent finder using the -lf_prm, \n");
	    printf("-pf_prm, or -sol_prm arguments. Note, the error occured on residue %d (internal numbering scheme). Printing info for problematic residue and terminating the program! \n\n",res_nr[i]);

            int this_start = res_start[res_nr[mem[i]-1]-1];
            int this_end   = res_end[res_nr[mem[i]-1]-1];

            printf("%9s %9s \n","res_id","res_name");
            printf("%9s-%9s \n","---------","---------");
            for(j=this_start; j<this_end; j++)
            {
                printf("%9d %9s \n",res_nr[j],res_name[j].c_str());
            }
	}

	MPI_Finalize();
	exit(EXIT_SUCCESS);
    }

    return new_i;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function assigns atoms to the leaflets                                                              //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::get_leaflets(int leaf,string leaflet_finder_param_name,int b_lf_param)
{
    //set the target leaflet
    leaflet = leaf;

    //resize leaflets vector to hold data
    leaflets.resize(num_atoms,0);

    //assign atoms to a leaflet or neither
    leaflet_finder(num_atoms,atom_nr,res_nr,res_name,atom_name,r_ref,world_rank,leaflets,res_start,res_end,b_lf_param,leaflet_finder_param_name);

    //get the target leaflet
    int mem_size = count_leaflet_atoms(leaflet,leaflets,num_atoms);
    target_leaflet.resize(mem_size);
    get_membrane(leaflet,leaflets,num_atoms,target_leaflet);

    //get the opposing leaflet
    op_leaflet = get_opposing_leaflet(leaflet);
    mem_size = count_leaflet_atoms(op_leaflet,leaflets,num_atoms);
    opposing_leaflet.resize(mem_size);
    get_membrane(op_leaflet,leaflets,num_atoms,opposing_leaflet);

    //get full membrane
    mem_size = count_leaflet_atoms(0,leaflets,num_atoms);
    full_membrane.resize(mem_size);
    get_membrane(0,leaflets,num_atoms,full_membrane);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function writes a pdb with the leaflets indicated by beta value                                     //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::write_leaflets(string lf_pdb_file_name,int b_lf_pdb)
{
    dump_leaflets_pdb(world_rank,b_lf_pdb,lf_pdb_file_name,beta,leaflets,num_atoms,box_ref,atom_nr,res_nr,res_name,atom_name,r_ref,title,weight,element,chain_id);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function sets the beta value according to the leaflets                                              //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::set_beta_lf()
{
    set_beta_leaflet_finder(beta,leaflets,num_atoms);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function prints info about the leaflet composition                                                  //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::get_leaflet_stats()
{
    get_lipid_types(target_leaflet.size(),target_leaflet,res_nr,res_name,world_rank,leaflet);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function prints info about the lipid selection                                                      //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::get_lipid_selection_stats(vector <string> lip_t,string title)
{
    lip_sel_stats(target_leaflet.size(),target_leaflet,res_nr,res_name,world_rank,lip_t.size(),lip_t,title);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function counts how many lipids there are in the target leaflet                                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::count_target_lipids()
{
    int i = 0;
    int lip_count = 0;
    int prev_lip = -1;

    for(i=0; i<target_leaflet.size(); i++) //loop membrane atoms
    {
        if(res_nr[target_leaflet[i]-1] != prev_lip)
        {
            lip_count++;
        }
        prev_lip = res_nr[target_leaflet[i]-1];
    }

    return lip_count;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function counts how many lipids there are in the opposing leaflet                                    //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::count_opposing_lipids()
{
    int i = 0;
    int lip_count = 0;
    int prev_lip = -1;

    for(i=0; i<opposing_leaflet.size(); i++) //loop membrane atoms
    {
        if(res_nr[opposing_leaflet[i]-1] != prev_lip)
        {
            lip_count++;
        }
        prev_lip = res_nr[opposing_leaflet[i]-1];
    }

    return lip_count;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function counts how many lipids there are in the target leaflet of a type                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::count_target_lipids_type(vector <string> lip_t)
{
    int i = 0;
    int j = 0;
    int lip_count = 0;
    int prev_lip = -1;

    for(i=0; i<target_leaflet.size(); i++) //loop membrane atoms
    {
        for(j=0; j<lip_t.size(); j++) //loop over lipid types
        {
            if(strcmp(res_name[target_leaflet[i]-1].c_str(), lip_t[j].c_str()) == 0) //lipid type is correct
            {
                if(res_nr[target_leaflet[i]-1] != prev_lip)
                {
                    lip_count++;
                }
                prev_lip = res_nr[target_leaflet[i]-1];
            }
        }
    }

    return lip_count;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function counts how many lipids there are in the opposing leaflet of a type                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::count_opposing_lipids_type(vector <string> lip_t)
{
    int i = 0;
    int j = 0;
    int lip_count = 0;
    int prev_lip = -1;

    for(i=0; i<opposing_leaflet.size(); i++) //loop membrane atoms
    {
        for(j=0; j<lip_t.size(); j++) //loop over lipid types
        {
            if(strcmp(res_name[opposing_leaflet[i]-1].c_str(), lip_t[j].c_str()) == 0) //lipid type is correct
            {
                if(res_nr[opposing_leaflet[i]-1] != prev_lip)
                {
                    lip_count++;
                }
                prev_lip = res_nr[opposing_leaflet[i]-1];
            }
        }
    }

    return lip_count;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the first atom of the current target lipid                                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::t_lip_start(int i)
{
    return res_start[res_nr[target_leaflet[i]-1]-1];
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the last atom of the current target lipid                                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::t_lip_end(int i)
{
    return res_end[res_nr[target_leaflet[i]-1]-1];
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the first atom of the current opposing lipid                                        //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::o_lip_start(int i)
{
    return res_start[res_nr[opposing_leaflet[i]-1]-1];
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the last atom of the current opposing lipid                                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::o_lip_end(int i)
{
    return res_end[res_nr[opposing_leaflet[i]-1]-1];
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the first atom of the current full mem lipid                                        //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::fm_lip_start(int i)
{
    return res_start[res_nr[full_membrane[i]-1]-1];
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the last atom of the current full mem lipid                                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::fm_lip_end(int i)
{
    return res_end[res_nr[full_membrane[i]-1]-1];
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the last atom of the current lipid so after i++ we have the next lipid              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::next_target_lipid(int i)
{
    return jump_to_res_end(i,target_leaflet.size(),res_end,res_nr,target_leaflet,res_name,world_rank,res_start);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the last atom of the current lipid so after i++ we have the next lipid              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::next_opposing_lipid(int i)
{
    return jump_to_res_end(i,opposing_leaflet.size(),res_end,res_nr,opposing_leaflet,res_name,world_rank,res_start);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the last atom of the current lipid so after i++ we have the next lipid              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::next_full_mem_lipid(int i)
{
    return jump_to_res_end(i,full_membrane.size(),res_end,res_nr,full_membrane,res_name,world_rank,res_start);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function detects jumps in the z direction                                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::check_lip_jumps_z()
{
    int i = 0;
    int j = 0;

    if(current_frame == 0) //store the coordinate
    {
        prev_z.resize(num_atoms,0.0);

        for(i=0; i<num_atoms; i++)
        {
            prev_z[i] = r[i][2];
        }
    }
    else if(current_frame >= 1)
    {
        for(i=0; i<full_membrane.size(); i++) //loop over system atoms
        {
            int min  = fm_lip_start(i);
            int max  = fm_lip_end(i);
            int jump = 0;

            //jump to next lipid
            i = next_full_mem_lipid(i);

            for(j=min; j<=max; j++)
            {
                if(r[j][2] - prev_z[j] > 0.5*box[ZZ][ZZ]) //jumped accross bottom of box
                {
                    jump = 1;
                }
                else if(r[j][2] - prev_z[j] < -1*0.5*box[ZZ][ZZ]) //jumped accross top of box
                {
                    jump = 1;
                }

                //update prev coord
                prev_z[j] = r[j][2];
            }
  
            if(jump == 1)
            {
                printf("Residue %d (%s) has moved more than half the box z-dimension %f between steps (%d:%d). This could indicate a fragmented leaflet. \n",res_nr[i],res_name[i].c_str(),box[ZZ][ZZ],get_frame_global()-1,get_frame_global());
            }
        }
    }
}


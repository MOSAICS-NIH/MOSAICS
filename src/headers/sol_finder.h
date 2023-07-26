
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function determines which atoms belong to water and which do not.                                   //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void w_finder(vector <int> &b_w,vector <string> &res_name,int num_atoms,int b_sf_param,string solvent_finder_param_name)
{
    const int num_res = 3;
    int i = 0;
    int j = 0;
    sv1d residue_types(0);

    residue_types.push_back("W");
    residue_types.push_back("TIP3");
    residue_types.push_back("WF");

    //check if user provided additional parameters
    if(b_sf_param == 1)
    {
        Index solvent_finder_param;
        solvent_finder_param.get_index(solvent_finder_param_name);

        //add new parameters to lip_parameters
        for(i=0; i<solvent_finder_param.index_s.size(); i++)
        {
            residue_types.push_back(solvent_finder_param.index_s[i]);
        }
    }

    for(i=0; i<num_atoms; i++)
    {
        b_w[i] = 0;

        for(j=0; j<residue_types.size(); j++)
        {
            if(strcmp(res_name[i].c_str(), residue_types[j].c_str()) == 0) //atom belongs to the solvent
            {
                b_w[i] = 1;
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function counts how many atoms are in the solvent.                                                  //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int count_sol_atoms(vector <int> &b_sol,int num_atoms)
{
    int i = 0;
    int count = 0;

    for(i=0; i<num_atoms; i++)
    {
        if(b_sol[i] == 1)
        {
            count++;
        }
    }
    return count;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function counts how many molecules are in the solvent.                                                  //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int count_sol_mol(vector <int> &b_sol,int num_atoms,vector <int> &res_nr)
{   
    int i = 0;
    int count = 0;
    int prev_res = -1;
   
    for(i=0; i<num_atoms; i++)
    {   
        if(b_sol[i] == 1 && res_nr[i] != prev_res)
        {   
            count++;
            prev_res = res_nr[i];
        }
    }
    return count;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function fills an index with the solvent atoms.                                                     //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_sol(vector <int> &b_sol,int num_atoms,vector <int> &sol)
{
    int i = 0;
    int count = 0;

    for(i=0; i<num_atoms; i++)
    {
        if(b_sol[i] == 1)
        {
            sol[count] = i + 1;
            count++;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function sets the beta value according to whether an atom belongs to the solvent or not.            //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void set_beta_sf(vector <double> &beta,vector <int> &b_sol,int num_atoms)
{
    int i = 0;

    for(i=0; i<num_atoms; i++)
    {
        if(b_sol[i] == 1) //sol atom
        {
            beta[i] = 1;
        }
        else //not sol atom
        {
            beta[i] = 0;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes a pdb with the solvent indicated by beta value.                                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void dump_sol_pdb(int world_rank,int b_sf_pdb,string sf_pdb_file_name,vector <double> &beta,vector <int> &b_sol,int num_atoms,
                  matrix box_ref,vector <int> &atom_nr,vector <int> &res_nr,vector <string> &res_name,vector <string> &atom_name,rvec *r_ref,
                  char title[200],vector <double> &weight,vector <string> &element,vector <char> &chain_id)
{
    FILE *sf_pdb_file;

    if(world_rank == 0 && b_sf_pdb == 1)
    {
        sf_pdb_file = fopen(sf_pdb_file_name.c_str(), "w");
        if(sf_pdb_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",sf_pdb_file_name.c_str());
        }
        set_beta_sf(beta,b_sol,num_atoms);
        write_frame_pdb(box_ref,num_atoms,atom_nr,res_nr,res_name,atom_name,r_ref,title,world_rank,&sf_pdb_file,beta,weight,element,chain_id,0);
        fclose(sf_pdb_file);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints information about the water molecuels.                                               //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void print_w_stats(vector <int> &sol,vector <int> &res_nr,vector <string> &res_name,int world_rank,int num_sol_mol)
{
    if(world_rank == 0)
    {
        vector <string> sol_types;

        int i = 0;
        int j = 0;
        int prev_res = -1;
        int count = 0;

        //count solvent types
        for(i=0; i<sol.size(); i++)
        {
            int recorded = 0;
            int size = sol_types.size();

            for(j=0; j<size; j++)
            {
                if(strcmp(res_name[sol[i]-1].c_str(), sol_types[j].c_str()) == 0) //res type already recorder
                {
                    recorded = 1;
                }
            }

            if(recorded == 0)
            {
                sol_types.push_back(res_name[sol[i]-1]);
            }
        }

        printf("There are %6d atoms, %9d molecules and %2d molecule types in the solvent as summarized below: \n",sol.size(),num_sol_mol,sol_types.size());
        printf(" %2s %10s %11s \n","index","type","count");
        printf("-%5s-%10s-%11s\n","-----","----------","-----------");

        //count molecules for each type
        for(i=0; i<sol_types.size(); i++) //loop over sol types
        {
            int prev_res = -1;
            int count = 0;

            for(j=0; j<sol.size(); j++) //loop over solvent selection
            {
                if(strcmp(res_name[sol[j]-1].c_str(), sol_types[i].c_str()) == 0) //residue is of correct type
                {
                    if(res_nr[sol[j]-1] != prev_res) //new residue
                    {
                        count++;
                        prev_res = res_nr[sol[j]-1];
                    }
                }
            }
            printf(" %5d %10s %11d \n",i,sol_types[i].c_str(),count);
        }
        printf("\n");
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function runs solvent finder and gets a vector with the solvent                                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::get_solvent(string solvent_finder_param_name,int b_sf_param)
{
    //resize b_sol
    b_sol.resize(num_atoms);

    //find solvent/water atoms
    w_finder(b_sol,res_name,num_atoms,b_sf_param,solvent_finder_param_name);

    //count the water atoms
    int sol_size = count_sol_atoms(b_sol,num_atoms);

    //resize sol
    sol.resize(sol_size,0);

    //get sol atoms
    get_sol(b_sol,num_atoms,sol);

    //count water molecules
    num_waters = count_sol_mol(b_sol,num_atoms,res_nr);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes a pdb with the solvent indicated by beta value                                       //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::write_sol(string sf_pdb_file_name,int b_sf_pdb)
{
    dump_sol_pdb(world_rank,b_sf_pdb,sf_pdb_file_name,beta,b_sol,num_atoms,box_ref,atom_nr,res_nr,res_name,atom_name,
                 r_ref,title,weight,element,chain_id);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints information about the solvent                                                        //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::get_sol_stats()
{
    //count solvent molecules
    int num_sol_mol = count_sol_mol(b_sol,num_atoms,res_nr);

    //print sol stats
    print_w_stats(sol,res_nr,res_name,world_rank,num_sol_mol);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the first atom of the current water                                                 //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::sol_start(int i)
{
    return res_start[res_nr[sol[i]-1]-1];
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the last atom of the current water                                                  //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::sol_end(int i)
{
    return res_end[res_nr[sol[i]-1]-1];
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the last atom of the current water so after i++ we have the next water              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::next_water(int i)
{
    return jump_to_res_end(i,sol.size(),res_end,res_nr,sol);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function counts how many solvent molecules there are of a type                                       //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::count_waters_type(vector <string> sol_t)
{
    int i         = 0;
    int j         = 0;
    int sol_count = 0;
    int prev_sol  = -1;

    for(i=0; i<sol.size(); i++) //loop solvent atoms
    {
        for(j=0; j<sol_t.size(); j++) //loop over sol types
        {
            if(strcmp(res_name[sol[i]-1].c_str(), sol_t[j].c_str()) == 0) //sol type is correct
            {
                if(res_nr[sol[i]-1] != prev_sol)
                {
                    sol_count++;
                }
                prev_sol = res_nr[sol[i]-1];
            }
        }
    }

    return sol_count;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects data for each residue of the protein                                               //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void collect_protein_res_data_i(int world_rank,int world_size,int prot_res_count,int res_data[])
{
    //now we collect the contacts from all nodes
    int dat = 0;                     //used for collecting data
    int i   = 0;                     //standard variable used in loops
    int j   = 0;                     //standard variable used in loops
    int world_dat[world_size];       //holds data collected from all cores

    for(i=0; i<prot_res_count; i++) //loop over protein residues
    {
        dat = res_data[i];

        MPI_Gather(&dat, 1, MPI_INT, world_dat, 1, MPI_INT, 0,MPI_COMM_WORLD);

        if(world_rank == 0)
        {
            dat = 0;

            for(j=0; j<world_size; j++) //loop over world
            {
                dat = dat + world_dat[j];
            }
            res_data[i] = dat;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes data to a pdb file as the beta factor                                                //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_data_pdb(matrix box_ref,int num_atoms,vector <int> &atom_nr,vector <int> &res_nr,vector <string> &res_name,vector <string> &atom_name,
                    rvec *r_ref,real time,int step,vector <double> &beta,vector <double> &weight,vector <char> &chain_id,
                    vector <string> &element,char title[200],int world_rank,string file_name)
{
    //open file and write data
    FILE *my_file = fopen(file_name.c_str(), "w");
    if(my_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",file_name.c_str());
    }
    else
    {
        //write data to pdb file
        write_frame_pdb(box_ref,num_atoms,atom_nr,res_nr,res_name,atom_name,r_ref,title,world_rank,&my_file,beta,weight,element,chain_id,0);

        //close file
        fclose(my_file);
    }
}


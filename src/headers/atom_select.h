
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function terminates the program                                                                     //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void terminate_program()
{
    MPI_Finalize();
    exit(EXIT_SUCCESS);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function analyzes a string and checks if it is an int (1) or a range of ints (1-10)                 //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int check_int_sel(string name)
{
    int i        = 0;      //standard variable used in loops
    int is_digit = 1;      //tells if the string is an int
    int range    = 0;      //counts the occurance of -

    //check for non numerical chars
    for(i=0; i<name.length(); i++) //loop over string
    {
        if(isdigit(name[i]) == 0) //a non numerical char was found
        {
            if(i > 0 && i << name.length()-1 && name[i] == '-') //check for a range indicator
            {
                range++;
            }
            else
            {
                is_digit = 0;
            }
        }
    }

    if(range > 1) 
    {
        is_digit = 0;
    }
    else if(is_digit == 1 && range == 1)
    {
        is_digit = 2; 
    }

    return is_digit;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function analyzes the current string and determines the instruction type                            //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int check_instruction_type(string this_string)
{
        int type = -1;

        //characterize the string as either a selection or an operator
        if(strcmp(this_string.c_str(), "id") == 0) //select by atom number
        {
            type = 0;
        }
        else if(strcmp(this_string.c_str(), "resi") == 0) //select by residue number
        {
            type = 1;
        }
        else if(strcmp(this_string.c_str(), "atom") == 0) //select by atom name
        {
            type = 2;
        }
        else if(strcmp(this_string.c_str(), "resn") == 0) //select by residue name
        {
            type = 3;
        }
        else if(strcmp(this_string.c_str(), "(") == 0) //open bracket
        {
            type = 6;
        }  
        else if(strcmp(this_string.c_str(), ")") == 0) //close bracket
        {
            type = 7;
        }
        else if(strcmp(this_string.c_str(), "and") == 0) //and
        {
            type = 8;
        }
        else if(strcmp(this_string.c_str(), "or") == 0) //compound
        {
            type = 9;
        }
        else if(strcmp(this_string.c_str(), "not") == 0) //invert selection
        {
            type = 10;
        }
        else if(strcmp(this_string.c_str(), "prot") == 0) //select protein
        {
            type = 11;
        }
        else if(strcmp(this_string.c_str(), "sol") == 0) //select solvent
        {
            type = 12;
        }
        else if(strcmp(this_string.c_str(), "upper") == 0) //select upper leaflet
        {
            type = 13;
        }
        else if(strcmp(this_string.c_str(), "lower") == 0) //select lower leaflet
        {
            type = 14;
        }
        else if(strcmp(this_string.c_str(), "hydrogen") == 0) //select hydrogens
        {
            type = 15;
        }
        else if(check_int_sel(this_string) == 1 || check_int_sel(this_string) == 2) //an atom/residue id was provided
        {
            type = 4;
        }
        else //assume an atom/residue name was provided
        {
            type = 5; 
        }

        return type;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function searches for the closing bracket                                                           //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int find_close_bracket(sv1d sel_text,int i)
{
    int j           = 0;   //standard variable used in loops
    int close_index = -1;  //loop index where bracket close is found
    int bracket     = 1;   //keep track of how many open brackets we have

    //search for bracket close
    for(j=i+1; j<sel_text.size(); j++) 
    {
        //analyze string type
        if(check_instruction_type(sel_text[j]) == 6)
        {
            bracket++;
        }
        else if(check_instruction_type(sel_text[j]) == 7)
        {
            bracket--;
        }

        //check for open ended brackets
        if(j==sel_text.size()-1 && bracket > 0)
        {
            printf("open ended bracket detected \n");
            terminate_program();
        }

        //bracket close found
        if(bracket == 0)
        {
            close_index = j;
            goto exit_loop;
        }
    }
    exit_loop:;
    
    return close_index;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function determines the range for the selection                                                     //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_sel_bounds(int *min,int *max,string sel_text)
{
    int i                 = 0; 
    int found_deliminator = 0; 

    string min_s; 
    string max_s; 

    for(i=0; i<sel_text.size(); i++) //loop over the text
    {
        if(sel_text[i] == '-')
        {
            found_deliminator = 1;
        }
        else 
        {
            if(found_deliminator == 0)
            {
                min_s.push_back(sel_text[i]);
            }
            else if(found_deliminator == 1)
            {
                max_s.push_back(sel_text[i]);
            }
        }
    } 
   
    *min = atoi(min_s.c_str());   

    if(max_s.size() == 0) //no deliminator detected
    {
        *max = atoi(min_s.c_str());
    }
    else //a deliminator was detected
    {
        *max = atoi(max_s.c_str());
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function compiles a list of names for the selection                                                 //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
sv1d get_sel_list(string sel_text)
{
    int    i     = 0;    //standard variable used in loops

    sv1d   list(0);      //make a list of names
    string this_string;  //the current item to add to the list

    if(sel_text[0] == '+' || sel_text[sel_text.size()-1] == '+')
    {
        printf("Selection lists should not begin or end with + (%s) \n",sel_text.c_str());
        terminate_program();
    }
    else
    {
        for(i=0; i<sel_text.size(); i++) //loop over the text
        {
            if(sel_text[i] == '+')
            {
                list.push_back(this_string);
                this_string = "";
            }
            else if(i == sel_text.size()-1)
            {
                this_string.push_back(sel_text[i]);
                list.push_back(this_string);
            }
            else
            {
                this_string.push_back(sel_text[i]);
            }
        }
    }

    return list; 
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function makes an atom selection by the atom id                                                     //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
iv1d select_by_atom_id(Trajectory traj,int min,int max)
{
    int i = 0;
 
    iv1d selection(0,0);

    for(i=0; i<traj.atoms(); i++) //loop over atoms
    {
        if(traj.atom_nr[i] >= min && traj.atom_nr[i] <= max) 
        {
            selection.push_back(traj.atom_nr[i]);
        }
    }

    return selection;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function makes an atom selection by the res id                                                      //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
iv1d select_by_res_id(Trajectory traj,int min,int max)
{
    int i = 0;

    iv1d selection(0,0);

    for(i=0; i<traj.atoms(); i++) //loop over atoms
    {
        if(traj.res_nr[i] >= min && traj.res_nr[i] <= max)
        {
            selection.push_back(traj.atom_nr[i]);
        }
    }

    return selection;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function makes an atom selection by the atom name                                                   //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
iv1d select_by_atom_name(Trajectory traj,sv1d list)
{
    int i = 0;
    int j = 0;

    iv1d selection(0,0);

    for(i=0; i<traj.atoms(); i++) //loop over atoms
    {
        for(j=0; j<list.size(); j++) //loop over list 
        {
            if(strcmp(traj.atom_name[i].c_str(), list[j].c_str()) == 0) 
            {
                selection.push_back(traj.atom_nr[i]);
                goto exit_loop_1;    
            }
        }
        exit_loop_1:;
    }

    return selection;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function makes an atom selection by the res name                                                    //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
iv1d select_by_res_name(Trajectory traj,sv1d list)
{
    int i = 0;
    int j = 0;

    iv1d selection(0,0);

    for(i=0; i<traj.atoms(); i++) //loop over atoms
    {
        for(j=0; j<list.size(); j++) //loop over list 
        {
            if(strcmp(traj.res_name[i].c_str(), list[j].c_str()) == 0)
            {
                selection.push_back(traj.atom_nr[i]);
                goto exit_loop_2;
            }
        }
        exit_loop_2:;
    }

    return selection;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function selects the protein atoms                                                                  //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
iv1d select_prot(Trajectory traj)
{
    int i = 0;

    iv1d selection(0,0);

    for(i=0; i<traj.prot.size(); i++) //loop over protein atoms
    {
        selection.push_back(traj.prot[i]);
    }

    return selection;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function selects the solvent atoms                                                                  //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
iv1d select_sol(Trajectory traj)
{
    int i = 0;

    iv1d selection(0,0);

    for(i=0; i<traj.sol.size(); i++) //loop over solvent atoms
    {
        selection.push_back(traj.sol[i]);
    }

    return selection;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function selects the upper leaflet                                                                  //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
iv1d select_upper(Trajectory traj)
{
    int i = 0;
    
    iv1d selection(0,0);

    for(i=0; i<traj.atoms(); i++) //loop over system atoms
    {
        if(traj.leaflets[i] == 1) //upper
        {
            selection.push_back(traj.atom_nr[i]);
        }
    }

    return selection;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function selects the lower leaflet                                                                  //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
iv1d select_lower(Trajectory traj)
{
    int i = 0;

    iv1d selection(0,0);

    for(i=0; i<traj.atoms(); i++) //loop over system atoms
    {
        if(traj.leaflets[i] == 2) //lower
        {
            selection.push_back(traj.atom_nr[i]);
        }
    }

    return selection;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function selects hydrogen atoms (assumed to be atoms beginning with H)                              //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
iv1d select_hydrogen(Trajectory traj)
{
    int i = 0;

    iv1d selection(0,0);

    for(i=0; i<traj.atoms(); i++) //loop over system atoms
    {
        if(traj.atom_name[i].at(0) == 'H') //atom is a hydrogen
        {
            selection.push_back(traj.atom_nr[i]);
        }
    }

    return selection;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function combines two selections                                                                    //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
iv1d prune(iv1d &sel_1,iv1d &sel_2)
{
    int i = 0;
    int j = 0;
    int k = 0;

    iv1d selection(0,0);

    for(i=0; i<sel_1.size(); i++) //loop over first selection
    {
        for(j=0; j<sel_2.size(); j++) //loop over second selection
        {
            if(sel_1[i] == sel_2[j])
            {
                int found = 0;
                for(k=0; k<selection.size(); k++) //loop over main selection atoms
                {
                    if(sel_1[i] == selection[k])
                    {
                        found = 1;
                    }
                }
  
                if(found == 0)
                {
                    selection.push_back(sel_1[i]);
                }
            }         
        }
    }

    return selection; 
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function combines two selections                                                                    //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
iv1d combine(iv1d &sel_1,iv1d &sel_2)
{
    int i = 0;
    int j = 0;

    iv1d selection(0,0);

    //check the first selection
    for(i=0; i<sel_1.size(); i++) //loop over first selection
    {
        int found = 0; 
        for(j=0; j<selection.size(); j++) //loop over combined selection
        {
            if(sel_1[i] == selection[j])
            {
                found = 1;
            }
        }
 
        if(found == 0)
        {
            selection.push_back(sel_1[i]);
        }
    }

    //check the second selection
    for(i=0; i<sel_2.size(); i++) //loop over second selection
    {
        int found = 0; 
        for(j=0; j<selection.size(); j++) //loop over combined selection
        {
            if(sel_2[i] == selection[j])
            {
                found = 1;
            }
        }
 
        if(found == 0)
        {
            selection.push_back(sel_2[i]);
        }
    }

    return selection;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function inverts the selection                                                                      //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
iv1d invert(Trajectory traj,iv1d current_sel)
{
    int i = 0;
    int j = 0;

    iv1d selection(0,0);

    for(i=0; i<traj.atoms(); i++) //loop over system atoms
    {
        int found = 0;

        for(j=0; j<current_sel.size(); j++) //loop over current selection
        {
            if(traj.atom_nr[i] == current_sel[j])
            {
                found = 1;
            }
        }

        if(found == 0)
        {
            selection.push_back(traj.atom_nr[i]);
        }
    }

    return selection;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function updates the selection                                                                      //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void update_selection(Trajectory traj,int *b_invert,int op,iv1d &current_sel,iv1d &this_selection)
{
    if(*b_invert == 1)
    {
        current_sel = invert(traj,current_sel);
    }

    if(op == 1)
    {
        this_selection = combine(this_selection,current_sel);
    }
    else if(op == -1)
    {
        this_selection = prune(this_selection,current_sel);
    }
    *b_invert = 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function makes an atom selection                                                                    //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
iv1d grab_atoms(Trajectory traj,sv1d sel_text)
{
    int i        = 0;          //standard variable used in loops
    int j        = 0;          //standard variable used in loops
    int op       = 1;          //tells how to handle new data and this_selection
    int b_invert = 0;          //invert the current selection?

    iv1d this_selection(0,0);  //holds the long term selection
    iv1d current_sel(0,0);     //holds the current selection
 
    //parse the selection text
    for(i=0; i<sel_text.size(); i++) //loop over items in text
    {
        int s_type = check_instruction_type(sel_text[i]);

        if(s_type == 0 || s_type == 1 || s_type == 2 || s_type == 3) //selection type [0:id, 1:res_nr, 2:atom_name, 3:res_name]
        {
            int next_s_type = -1;

            //check that the next entry completes the selection
            if(i+1 < sel_text.size())
            { 
                next_s_type = check_instruction_type(sel_text[i+1]);
     
                if(s_type == 0 && next_s_type == 4) //select by atom id
                {
                    int min = 0;
                    int max = 0;
                    get_sel_bounds(&min,&max,sel_text[i+1]);

                    current_sel = select_by_atom_id(traj,min,max);

                    update_selection(traj,&b_invert,op,current_sel,this_selection);
                }
                else if(s_type == 1 && next_s_type == 4) //select by res id
                {
                    int min = 0;
                    int max = 0;
                    get_sel_bounds(&min,&max,sel_text[i+1]);

                    current_sel = select_by_res_id(traj,min,max);

                    update_selection(traj,&b_invert,op,current_sel,this_selection);
                }
                else if(s_type == 2 && next_s_type == 5) //select by atom name
                {
                    sv1d list   = get_sel_list(sel_text[i+1]);
                    current_sel = select_by_atom_name(traj,list);

                    update_selection(traj,&b_invert,op,current_sel,this_selection);
                }
                else if(s_type == 3 && next_s_type == 5) //select by res name
                {
                    sv1d list    = get_sel_list(sel_text[i+1]);
                    current_sel = select_by_res_name(traj,list);

                    update_selection(traj,&b_invert,op,current_sel,this_selection);
                }
                else 
                {
                    printf("an incompatible selection pair was detected (%s) (%s) \n",sel_text[i].c_str(),sel_text[i+1].c_str());
                    terminate_program();
                }

                i++;
            }
            else
            {
                printf("a selection type was detected without a trailing indicator (%s) \n",sel_text[i].c_str());
                terminate_program();
            }
        }
        else if(s_type == 4) //atom/residue id
        {
            printf("an atom or residue id was detected without a preceding selection type (%s) \n",sel_text[i].c_str());
            terminate_program();
        }
        else if(s_type == 5) //atom/residue name
        {
            printf("an atom or residue name was detected without a preceding selection type (%s) \n",sel_text[i].c_str());
            terminate_program();
        }
        else if(s_type == 6) //open bracket
        {
            int start = i+1;
            int end   = find_close_bracket(sel_text,i) - 1;

            sv1d bracket_sel(0);
            for(j=start; j<=end; j++) //loop over text in brackets
            {
                bracket_sel.push_back(sel_text[j]);
            }
            current_sel = grab_atoms(traj,bracket_sel);

            update_selection(traj,&b_invert,op,current_sel,this_selection);

            i = end + 1;
        }
        else if(s_type == 7) //close bracket
        {
            printf("closed bracket detected without a preceding open bracket \n");
            terminate_program();
        }
        else if(s_type == 8) //and
        {
            op = -1; 
        }
        else if(s_type == 9) //or
        {
            op = +1;
        }
        else if(s_type == 10)
        {
            b_invert = 1; 
        }
        else if(s_type == 11) //select protein
        {
            current_sel = select_prot(traj); 

            update_selection(traj,&b_invert,op,current_sel,this_selection);
        }
        else if(s_type == 12) //select solvent
        {
            current_sel = select_sol(traj);

            update_selection(traj,&b_invert,op,current_sel,this_selection);
        }
        else if(s_type == 13) //select upper leaflet
        {
            current_sel = select_upper(traj);
    
            update_selection(traj,&b_invert,op,current_sel,this_selection);
        }
        else if(s_type == 14) //select solvent
        {
            current_sel = select_lower(traj);
    
            update_selection(traj,&b_invert,op,current_sel,this_selection);
        }
        else if(s_type == 15) //select hydrogens
        {
            current_sel = select_hydrogen(traj);

            update_selection(traj,&b_invert,op,current_sel,this_selection);
        }
    }
   
    return this_selection; 
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This is a class for making atom selections using a selection text                                        //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Selection
{
    public:
        iv1d sel;                                                         //This holds the atom numbers for the selection 

    public:
        void get_selection(Trajectory traj,sv1d sel_tex);                 //read in the index
        void highlight_sel(Trajectory traj,string pdb_filename);          //write a pdb file with the selection highlighted in the B-factor 
        iv1d tag_atoms(Trajectory traj);                                  //This function marks the selected atoms with a 1 and the rest with 0 
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function analyzes the selection text and the trajectory and makes the atom selection                //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Selection::get_selection(Trajectory traj,sv1d sel_text)
{
    sel = grab_atoms(traj,sel_text);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function highlights the selection using the B-factor of a pdb file                                  //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Selection::highlight_sel(Trajectory traj,string pdb_filename)
{
    int i = 0;
    int j = 0;
 
    dv1d    this_beta(traj.atoms(),0.0);                   //beta factor used in pdb file
  
    for(i=0; i<traj.atoms(); i++) //loop over system atoms
    {
        for(j=0; j<sel.size(); j++) //loop over selection
        {
            if(traj.atom_nr[i] == sel[j])
            {
                this_beta[i] = 1.0; 
            }
        }
    }

    //open pdb file for writing
    FILE *pdb_file;
    pdb_file = fopen(pdb_filename.c_str(), "w");
    if(pdb_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",pdb_filename.c_str());
    }
    else
    {
        write_frame_pdb(traj.ibox,traj.atoms(),traj.atom_nr,traj.res_nr,traj.res_name,traj.atom_name,traj.r_ref,traj.title,0,&pdb_file,this_beta,traj.weight,traj.element,traj.chain_id,1);
        fclose(pdb_file);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns a 1d vector with the selected atoms marked by 1 and other atoms by a 0             //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
iv1d Selection::tag_atoms(Trajectory traj)
{
    int i = 0;                           //standard variable used in loops

    iv1d tagged_atoms(traj.atoms(),0);   //holds an indicator for each atom

    //tag the selected atoms
    for(i=0; i<sel.size(); i++) //loop over selected atoms
    {
        tagged_atoms[sel[i]-1] = 1;
    }
    
    return tagged_atoms;
}

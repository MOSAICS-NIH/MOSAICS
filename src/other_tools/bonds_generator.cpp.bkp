
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <sstream>
#include <fstream>

using namespace std;

#include "../headers/multi_dim_vec.h"
#include "../headers/file_reader.h"
#include "../headers/grid_lt.h"
#include "../headers/common_routines.h"
#include "../headers/command_line_args.h"
#include "../headers/array.h"
#include "../headers/file_io_variable.h"
#include "../headers/file_naming.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function checks the filetype for acceptable extensions                                               //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
string check_extension_multi(string arg,string filename,sv1d &tag)
{
    int i             = 0;                  //standard variable used in loops
    int j             = 0;                  //standard variable used in loops
    int pass          = 0;                  //tells if the extension matches one of the tags
    string ext;                             //the extension of the file

    for(j=0; j<tag.size(); j++) //loop over acceptable extensions
    {
        int bad_extension = 0;                  //tells if the correct extension was provided
        int filename_size = filename.length();  //how long is the filename
        int tag_size      = tag[j].length();    //how long is the tag

        //make sure tag is not longer than filename
        if(tag_size > filename_size)
        {
            bad_extension = 1;
        }
        else //check for proper extension
        {
            for(i=filename_size-tag_size; i<filename_size; i++)
            {
                if(filename[i] != tag[j][i + tag_size - filename_size])
                {
                    bad_extension = 1;
                    break;
                }
            }
        }

        if(bad_extension == 0) //matches one of the acceptable extensions
        {
            pass = 1;
	    ext  = tag[j];
        }
    }

    //report error and terminate program
    if(pass == 0)
    {
        printf("The filename %s provided via the %s tag requires the %s extension. \n",filename.c_str(),arg.c_str(),tag[0].c_str());
        exit(EXIT_SUCCESS);
    }

    return ext;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function searches a data file for a tag and gives the position                                       //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int check_next_tag(sv2d &data_s,string main_tag,sv1d &alt_tag,int *pos_line,int *pos_items)
{
    int i      = 0;        //standard variable used in loops
    int j      = 0;        //standard variable used in loops
    int result = 0;        //tells if the tag was found

    for(i=*pos_line; i<data_s.size(); i++) //loop over lines
    {
        if(data_s[i].size() > 0) //prevent seg fault for empty lines
        {
            int test_1 = 0;
            int test_2 = 0;

            if(strcmp(data_s[i][0].c_str(), main_tag.c_str() ) == 0)
            {
                test_1 = 1;
            }
            else if(data_s[i].size() >= alt_tag.size())
            {
                int miss = 0;

                for(j=0; j<alt_tag.size(); j++) //loop over items in tag
                {
                    if(strcmp(data_s[i][j].c_str(), alt_tag[j].c_str()) != 0)
                    {
                        miss = 1;
                    }
                }

                if(miss == 0)
                {
                    test_2 = 1;
                }
            }

            if(test_1 == 1)
            {
                //record position of tag
                result     = 1;
                *pos_line  = i;
                *pos_items = 0;
		i          = data_s.size();
            }
            else if(test_2 == 1)
            {
                string full_tag;

                for(j=0; j<alt_tag.size(); j++)
                {
                     full_tag = full_tag + alt_tag[j];
                }

                //record position of tag
                result     = 1;
                *pos_line  = i;
                *pos_items = alt_tag.size()-1;
                i          = data_s.size();
            }
        }
    }
    return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads itp files looking for the correct molecule type                                       //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
sv3d check_itp_files(Data_file_vari &data,string working_dir)
{
    int i = 0;        //standard variable used in loops
    int j = 0;        //standard variable used in loops

    sv3d itp_files;   //holds data for the itp files

    //look for include statements
    for(i=0; i<data.data_s.size(); i++) //loop over lines
    {
        if(data.data_s[i].size() > 0) //prevent seg fault for empty lines
        {
            if(strcmp(data.data_s[i][0].c_str(), "#include" ) == 0 && data.data_s[i].size() > 1)
            {
                //extract file name
                string this_file_name;  //name of itp file to open

                for(j=0; j<data.data_s[i][1].length(); j++) //loop over characters in string
                {
                    if(data.data_s[i][1][j] != '"') //remove parentheses 
                    {
                        this_file_name.push_back(data.data_s[i][1][j]);
                    }
                }

                //create object to read in the current itp file
                Data_file_vari this_itp;
                int result = this_itp.get_data(working_dir + this_file_name);
                if(result == 0)
                {
                    printf("Could not open file %s. Will terminate program! \n",(working_dir + this_file_name).c_str());
                    exit(EXIT_SUCCESS);
		}
		itp_files.push_back(this_itp.data_s);
            }
        }
    }

    return itp_files;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function moves to the next line without a comment (;)                                                //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void skip_comments(sv3d &itp_files,int *pos_line,int j)
{
    int k = 0;    //standard variable used in loops

    for(k=0; k<1; )
    {
        if(itp_files[j][*pos_line][0][0] == ';')
        {
            *pos_line = *pos_line + 1;
        }
        else
        {
            k = 1;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This is the main function that executes other functions                                                   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, const char * argv[]) 
{
    //Here we define some variables used throughout
    string in_file_name;           //Name of input file data file
    string out_file_name;          //Name of output file
    string working_dir;            //The directory with the topology file
    int i             = 0;         //General variable used in loops
    int j             = 0;         //General variable used in loops
    int k             = 0;         //General variable used in loops
    int b_constraints = 0;         //Add bonds from [constraints] sections?
    sv1d cl_tags;                  //Holds a list of command line tags for the program

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Set program name/description and print info                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    string program_name = "Bonds Generator";

    print_credits(argc,argv,program_name);

    string program_description = "Bonds Generator is a tool that reads in a GROMACS topology file and generates a bonds list to be used with other MOSAICS tools.";

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Analyze the input arguments                                                                               //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    start_input_arguments(argc,argv,program_description);
    add_argument_s(argc,argv,"-d",   in_file_name,               "Input topology file (top)"                                , cl_tags, nullptr,      1);
    add_argument_s(argc,argv,"-o",   out_file_name,              "Input file with bonds list (crd)"                         , cl_tags, nullptr,      1);
    add_argument_s(argc,argv,"-wd",  working_dir,                "Directory containing the topology file"                   , cl_tags, nullptr,      1);
    add_argument_i(argc,argv,"-con", &b_constraints,             "Include bonds from [constraints] sections? (0:no, 1:yes)" , cl_tags, nullptr,      0);
    conclude_input_arguments(argc,argv,program_name,cl_tags);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check file extensions                                                                                     //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    check_extension("-o",out_file_name,".crd");

    //set acceptable extensions for topology file (.top or .psf)
    sv1d tags(2);
    tags[0] = ".top";
    tags[1] = ".psf";
    string this_ext = check_extension_multi("-d",in_file_name,tags);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Report bonds from a GROMACS topology file                                                                 //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(strcmp(this_ext.c_str(), ".top" ) == 0) //GROMACS topology 
    {
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Read molecules list                                                                                       //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //get filename of topology file
        string top_file_name = working_dir + in_file_name;

        //create object to store main topology data
        Data_file_vari data;
        int top_read_result = data.get_data(top_file_name);
        if(top_read_result == 0)
        {
           printf("Could not open file %s. Will terminate program! \n",top_file_name.c_str()); 
           exit(EXIT_SUCCESS); 
        } 

        //create tags for [molecules] section
        string main_tag = "[molecules]";    
        sv1d   alt_tag;
        alt_tag.push_back("[");
        alt_tag.push_back("molecules");
        alt_tag.push_back("]");

        int top_line_pos  = 0;   //hold the line position in topology file
        int top_items_pos = 0;   //hold the item position in topology file

        //check topology file for a [molecules] section
        int result = check_next_tag(data.data_s,main_tag,alt_tag,&top_line_pos,&top_items_pos);

        //read the itp files
        sv3d itp_files = check_itp_files(data,working_dir);

        int bonds_sum = 0;  //keeps track of how much should be added to the bonds data

        iv2d bonds(0,iv1d(2,0));    //stores the bonds

        if(result == 1) //molecules section found
        {
            for(i=top_line_pos+1; i<data.data_s.size(); i++) //loop over lines in file coming after tag. allowed since no [sections] are allowed after [molecules]
            {
                if(data.data_s[i].size() > 0) //empty line will crash the program
                {
                    if(strcmp(data.data_s[i][0].c_str(), ";" ) != 0) //ignore lines beginning with ;
                    {
                        if(data.data_s[i].size() > 1) //atleast 2 items needed (molt_type and count)
                        {
                            string this_mol   = data.data_s[i][0];                  //the current molecule type
                            int    this_count = atoi(data.data_s[i][1].c_str());    //number of copies of current molecule
                            int    num_atoms  = 0;                                  //how many atoms in the molecule

                            //create tags for [moleculetype] section
                            string moltype_main_tag = "[moleculetype]";  
                            sv1d   moltype_alt_tag;
                            moltype_alt_tag.push_back("[");
                            moltype_alt_tag.push_back("moleculetype");
                            moltype_alt_tag.push_back("]");

                            int type_found = 0;

                            iv2d these_bonds(0,iv1d(2,0));     //holds bonds for the current molecule type

                            //loop over itp files
                            for(j=0; j<itp_files.size(); j++) //loop over itp files
                            {
                                int pos_line  = 0;     //line index for the current itp file
                                int pos_items = 0;     //position in current line in current itp file
            		
                                //scan the itp files for a [moleculetype] section
                                while(check_next_tag(itp_files[j],moltype_main_tag,moltype_alt_tag,&pos_line,&pos_items))
                                {
                                    pos_line = pos_line + 1;
        
                                    //skip lines with ;
                                    skip_comments(itp_files,&pos_line,j);

                                    if(strcmp(itp_files[j][pos_line][0].c_str(), this_mol.c_str() ) == 0) //molecule type is correct
                                    {
                                        type_found = 1;

                                        //create a copy of the current itp file position used to look for additional [moleculetype] sections
                                        int next_type_line_pos  = pos_line;
                                        int next_type_items_pos = pos_items;

                                        //check for additional [moleculetype] sections in the current itp file
                                        int next_result = check_next_tag(itp_files[j],moltype_main_tag,moltype_alt_tag,&next_type_line_pos,&next_type_items_pos);

                                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                        //                                                                                                           //
                                        // Count atoms for the current molecule                                                                      //
                                        //                                                                                                           //
                                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                        //create tags for [atoms] section
                                        string atoms_main_tag = "[atoms]";
                                        sv1d   atoms_alt_tag;
                                        atoms_alt_tag.push_back("[");
                                        atoms_alt_tag.push_back("atoms");
                                        atoms_alt_tag.push_back("]");

                                        //set the position in the itp to match that for moleculetype. used to look for a [atoms] section
                                        int atoms_line_pos  = pos_line;
                                        int atoms_items_pos = pos_items;

                                        //check for a [atoms] section
                                        int atoms_result = check_next_tag(itp_files[j],atoms_main_tag,atoms_alt_tag,&atoms_line_pos,&atoms_items_pos);

                                        if(atoms_result == 1 && (next_result == 0 || atoms_line_pos < next_type_line_pos)) //[atoms] section found befor the next [moleculetype] section
                                        {
                                            atoms_line_pos++;

                                            //read in atoms
                                            for(k=0; k<1;)
                                            {
                                                if(atoms_line_pos < itp_files[j].size()) //line is inside file bounds
                                                {
                                                    if(itp_files[j][atoms_line_pos].size() >= 1) //line has at least 1 itmes
                                                    {
                                                        if(itp_files[j][atoms_line_pos][0][0] != '[') //still in current [atoms] section
                                                        {
                                                            if(itp_files[j][atoms_line_pos][0][0] != ';' && itp_files[j][atoms_line_pos][0][0] != '#')
                                                            {
                                                                //count atoms
                                                                num_atoms++;
                                                            }
                                                        }
                                                        else
                                                        {
                                                            k = 1;
                                                        }
                                                    }
                                                    atoms_line_pos++;
                                                }
                                                else
                                                {
                                                    k = 1;
                                                }
                                            }
                                        }

                                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                        //                                                                                                           //
                                        // Find bonds for the current molecule                                                                       //
                                        //                                                                                                           //
                                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                        //create tags for [bonds] section
                                        string bonds_main_tag = "[bonds]";   
                                        sv1d   bonds_alt_tag;
                                        bonds_alt_tag.push_back("[");
                                        bonds_alt_tag.push_back("bonds");
                                        bonds_alt_tag.push_back("]");

                                        //set the position in the itp to match that for moleculetype. used to look for a [bonds] section
                                        int bonds_line_pos  = pos_line;
                                        int bonds_items_pos = pos_items;

                                        //check for a [bonds] section
                                        int bonds_result = check_next_tag(itp_files[j],bonds_main_tag,bonds_alt_tag,&bonds_line_pos,&bonds_items_pos);

                                        if(bonds_result == 1 && (next_result == 0 || bonds_line_pos < next_type_line_pos)) //[bonds] section found befor the next [moleculetype] section
                                        {
                                            bonds_line_pos++;

                                            //read in bonds
                                            for(k=0; k<1;)
                                            {	 
                                                if(bonds_line_pos < itp_files[j].size()) //line is inside file bounds
                                                {
                                                    if(itp_files[j][bonds_line_pos].size() >= 2) //line has at least 2 itmes
                                                    {
                                                        if(itp_files[j][bonds_line_pos][0][0] != '[') //still in current [bonds] section
                                                        {
                                                            if(itp_files[j][bonds_line_pos][0][0] != ';' && itp_files[j][bonds_line_pos][0][0] != '#')
                                                            {
                                                                //store the bonds
                                                                iv1d this_bond(2,0);
                                                                this_bond[0] = atoi(itp_files[j][bonds_line_pos][0].c_str());
                                                                this_bond[1] = atoi(itp_files[j][bonds_line_pos][1].c_str());
 
                                                                these_bonds.push_back(this_bond);
                                                            }
                                                        }
                                                        else 
                                                        {
                                                            k = 1;
                                                        }
                                                    }
                                                    bonds_line_pos++;
                                                }
                                                else 
                                                {
                                                    k = 1;
                                                }
                                            }
                                        } 

                                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                        //                                                                                                           //
                                        // Find bonds from constraints for the current molecule                                                      //
                                        //                                                                                                           //
                                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                        if(b_constraints == 1)
                                        {                                 
                                            //create tags for [constraints] section
                                            string constraints_main_tag = "[constraints]";
                                            sv1d   constraints_alt_tag;
                                            constraints_alt_tag.push_back("[");
                                            constraints_alt_tag.push_back("constraints");
                                            constraints_alt_tag.push_back("]");

                                            //set the position in the itp to match that for moleculetype. used to look for a [constraints] section
                                            int constraints_line_pos  = pos_line;
                                            int constraints_items_pos = pos_items;

                                            //check for a [constraints] section
                                            int constraints_result = check_next_tag(itp_files[j],constraints_main_tag,constraints_alt_tag,&constraints_line_pos,&constraints_items_pos);

                                            if(constraints_result == 1 && (next_result == 0 || constraints_line_pos < next_type_line_pos)) //[constraints] section found befor the next [moleculetype] section
                                            {
                                                constraints_line_pos++;

                                                //read in constraints
                                                for(k=0; k<1;)
                                                {
                                                    if(constraints_line_pos < itp_files[j].size()) //line is inside file bounds
                                                    {
                                                        if(itp_files[j][constraints_line_pos].size() >= 2) //line has at least 2 itmes
                                                        {
                                                            if(itp_files[j][constraints_line_pos][0][0] != '[') //still in current [constraints] section
                                                            {
                                                                if(itp_files[j][constraints_line_pos][0][0] != ';' && itp_files[j][constraints_line_pos][0][0] != '#')
                                                                {
                                                                    //store the constraints
                                                                    iv1d this_bond(2,0);
                                                                    this_bond[0] = atoi(itp_files[j][constraints_line_pos][0].c_str());
                                                                    this_bond[1] = atoi(itp_files[j][constraints_line_pos][1].c_str());

                                                                    these_bonds.push_back(this_bond);
                                                                }
                                                            }
                                                            else
                                                            {
                                                                k = 1;
                                                            }
                                                        }
                                                        constraints_line_pos++;
                                                    }
                                                    else
                                                    {
                                                        k = 1;
                                                    }
                                                }
                                            }
                                        }
                                    }	
                                }
                            }

                            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                            //                                                                                                           //
                            // Add the new bonds                                                                                         //
                            //                                                                                                           //
                            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                            if(type_found == 1) //add bonds if any found
                            {
                                for(j=0; j<this_count; j++) //loop over copies of the current molecule
                                {
                                    for(k=0; k<these_bonds.size(); k++) //loop over the current molecules bonds
                                    {
                                        iv1d this_bond(2,0);
                                        this_bond[0] = these_bonds[k][0] + bonds_sum; 
                                        this_bond[1] = these_bonds[k][1] + bonds_sum; 
 
                                        bonds.push_back(this_bond);
                                    }
                                    bonds_sum = bonds_sum + num_atoms;
                                } 
                            }	
                            else
                            {
                                printf("Could not find [moleculetype] %s in itp files! Terminating program! \n",this_mol.c_str());
                                exit(EXIT_SUCCESS);
                            }	 
                        }
                    }
                }
                else //end of molecules section reached (empty line)
                {
                    i=data.data_s.size();
                }
            }
        }
        else
        {
            printf("could not find a [molecules] section in %s. Terminating program! \n",in_file_name.c_str());
            exit(EXIT_SUCCESS);
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Report bonds data                                                                                         //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        FILE *out_file;
        out_file = fopen(out_file_name.c_str(), "w");
        if(out_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",out_file_name.c_str());
        }
        else
        {
            for(i=0; i<bonds.size(); i++)
            {
                fprintf(out_file," %d %d \n",bonds[i][0],bonds[i][1]);
            }
            fclose(out_file);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Report bonds from a protein structure file (psf)                                                          //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(strcmp(this_ext.c_str(), ".psf" ) == 0) //psf 
    {
        iv2d bonds(0,iv1d(2,0));    //stores the bonds

        //get filename of topology file
        string top_file_name = working_dir + in_file_name;

        //create object to store main topology data
        Data_file_vari data;
        int top_read_result = data.get_data(top_file_name);
        if(top_read_result == 0)
        {
           printf("Could not open file %s. Will terminate program! \n",top_file_name.c_str());
           exit(EXIT_SUCCESS);
        }

        //check psf file for a !NBOND: section
        for(i=0; i<data.data_s.size(); i++) //loop over lines of psf file
        {
            if(data.data_s[i].size() > 1) //check that line has at least 2 entries
            {
                if(strcmp(data.data_s[i][1].c_str(), "!NBOND:" ) == 0) //bonds section found
                {
                    int num_bonds = atoi(data.data_s[i][0].c_str());
                    
                    i++; //move to the next line

                    for(j=0; j<num_bonds; j++)
                    {
                        int this_pos = (j%4)*2;

                        iv1d this_bond(2,0);
                        this_bond[0] = atoi(data.data_s[i][this_pos].c_str());
                        this_bond[1] = atoi(data.data_s[i][this_pos+1].c_str());
			bonds.push_back(this_bond);

                        if(this_pos == 6) //move to the next line
                        {
                            i++;
                        }
		    } 
		}
            }
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Report bonds data                                                                                         //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        FILE *out_file;
        out_file = fopen(out_file_name.c_str(), "w");
        if(out_file == NULL)
        {
            printf("failure opening %s. Make sure the file exists. \n",out_file_name.c_str());
        }
        else
        {
            for(i=0; i<bonds.size(); i++)
            {
                fprintf(out_file," %d %d \n",bonds[i][0],bonds[i][1]);
            }
            fclose(out_file);
        }
    }

    std::cout << "\nFormatting completed successfully" << "\n\n";

    return 0;
}


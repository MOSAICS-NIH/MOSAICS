
#include "gmx_lib/gmx_def.h"
#include "gmx_lib/gmx_bool.h"
#include "gmx_lib/gmx_util.h"
#include "gmx_lib/gmx_vec.h"
#include "gmx_lib/gmx_fit.h"
#include "gmx_lib/gmx_pdb.h"

#ifdef __APPLE__
#  define off64_t off_t
#  define fopen64 fopen
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function converts a string into an array of chars for use with c-based functions                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
char* cname(string file_name)
{
    int i = 0;      //standard variable used in loops
    char* cname;    //filename converted for use with c functions

    cname = (char *)malloc((file_name.length()+1) * sizeof(char));
    for(i=0; i<file_name.length(); i++)
    {
        cname[i] = file_name[i];
    }
    cname[file_name.length()] = '\0';
 
    return cname;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function analyzes a gro file. It gets the number of atoms in the system, determines how many frames  //
// are in the file, and stores the starting position of each.                                                // 
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void analyze_gro_file(FILE **in_file,int *num_atoms,int *frames,int world_rank,string in_file_name,
                      vector<int64_t> &pos,int64_t *filesize)
{
    if(world_rank == 0)
    {
        int  j = 0;               //standard variable used in loops
        char line[200];           //used to read in a line of data

        printf("Analyzing %s. \n",in_file_name.c_str());

        //store position of frame 1
        pos.push_back(ftell(*in_file));

        //loop over trajectory
        while(fgets(line, sizeof line, *in_file) != NULL)
        {
            if(j == 1)
            {
                *num_atoms = atoi(line);
            }
            if(j%(*num_atoms + 3) == (*num_atoms + 2) )
            {
                pos.push_back(ftell(*in_file));
            }
            j++;
        }
        //get the number of frames
        *frames = (j)/(*num_atoms + 3);

        //get the file size
        *filesize = ftell(*in_file);
 
        //rewind file for future reading
        rewind(*in_file);

        printf("Finished analyzing %s. \n",in_file_name.c_str());
        printf("Trajectory frames: %-20d \n\n",*frames);
    }

    //broadcast num_atoms, frames, and filesize
    MPI_Bcast(num_atoms, 1, MPI_INT,  0, MPI_COMM_WORLD);
    MPI_Bcast(frames,    1, MPI_INT,  0, MPI_COMM_WORLD);
    MPI_Bcast(filesize,  1, MPI_LONG, 0, MPI_COMM_WORLD);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function analyzes a reference gro file. It counts the number of atoms and the number of frames. It   //
// also extracts the box coordinates.                                                                        //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void analyze_gro_file_ref(FILE **ref_file,int *num_atoms,int *frames,int world_rank,matrix box_ref,int *box_dimension)
{
    int i               = 0;                    //standard variable used in loops 
    int j               = 0;                    //standard variable used in loops
    int result          = 0;                    //gives the rusult of reading the next string from a line 
    int line_offset     = 0;                    //stores the position in the current line
    int number_of_lines = 0;                    //how many lines are in the file
    char line[200];                             //stores a line
    char my_string[20];                         //stores the current string

    while(fgets(line, sizeof line, *ref_file) != NULL) 
    {
        if(number_of_lines == 1) //line contains the number of atoms
        {
            *num_atoms = atoi(line);
        }
        if(number_of_lines > 1 && number_of_lines == *num_atoms + 2) //line contains the box
        {
            int box_dimension_ref = 0;
            while(next_string(world_rank,200,line,my_string,20,&line_offset)) //extract the strings from the line
            {
                box_dimension_ref = box_dimension_ref + 1;
            }
            line_offset = 0;

            //store the box dimensions
            if(box_dimension_ref == 3)
            {
                result = next_string(world_rank,200,line,my_string,20,&line_offset);
                box_ref[XX][XX] = atof(my_string);

                result = next_string(world_rank,200,line,my_string,20,&line_offset);
                box_ref[YY][YY] = atof(my_string);

                result = next_string(world_rank,200,line,my_string,20,&line_offset);
                box_ref[ZZ][ZZ] = atof(my_string);
            }
            else if(box_dimension_ref == 9)
            {
                for(i=0; i<3; i++)
                {
                    for(j=0; j<3; j++)
                    {
                        result = next_string(world_rank,200,line,my_string,20,&line_offset);
                        box_ref[i][j] = atof(my_string);
                    }
                }
            }
        }
        number_of_lines = number_of_lines + 1;
    }
    //get the number of frames
    *frames = (number_of_lines)/(*num_atoms + 3);

    //rewind file for future reading
    rewind(*ref_file);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function analyzes a pdb file. It gets the number of atoms in the system, determines how many frames  //
// are in the file, and stores the starting position of each.                                                //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void analyze_pdb_file(FILE **in_file,int *num_atoms,int *frames,int world_rank,string in_file_name,
                      vector<int64_t> &pos,int64_t *filesize)
{
    if(world_rank == 0)
    {
        int result      = 0;            //gives the rusult of reading the next string from a line
        int line_offset = 0;            //stores the position in the current line
        char line[200];                 //stores a line
        char my_string[200];            //stores the current string

        printf("Analyzing %s. \n",in_file_name.c_str());

        //scan over the first frame and get the number of atoms
        while(fgets(line, sizeof line, *in_file) != NULL)
        {
            line_offset = 0;
            result = next_string(world_rank,200,line,my_string,200,&line_offset);

            if(strcmp(my_string, "ATOM") == 0) //atom
            {
                *num_atoms = *num_atoms + 1;
            }
            else if(strcmp(my_string, "ENDMDL") == 0) //end of frame
            {
                break;
            }
        }
        rewind(*in_file);

        //get the position of the first frame
        pos.push_back(ftell(*in_file)); 

        //count the frames and get the position of each 
        while(fgets(line, sizeof line, *in_file) != NULL)
        {
            line_offset = 0;
            result = next_string(world_rank,200,line,my_string,200,&line_offset);
            if(strcmp(my_string, "ENDMDL") == 0) //end of frame
            {
                pos.push_back(ftell(*in_file));
                *frames = *frames +1;
            }
        }

        //get the file size
        *filesize = ftell(*in_file);

        //rewind file for future reading
        rewind(*in_file);

        printf("Finished analyzing %s. \n",in_file_name.c_str());
        printf("Trajectory frames: %-20d \n\n",*frames);
    }

    //broadcast num_atoms, frames, and filesize
    MPI_Bcast(num_atoms, 1, MPI_INT,  0, MPI_COMM_WORLD);
    MPI_Bcast(frames,    1, MPI_INT,  0, MPI_COMM_WORLD);
    MPI_Bcast(filesize,  1, MPI_LONG, 0, MPI_COMM_WORLD);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function analyzes a reference pdb file. It counts the number of atoms and the number of frames. We   //
// get the box coordinates later when calling read_pdb_frame_by_char().                                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void analyze_pdb_file_ref(FILE **in_file,int *num_atoms,int *frames,int world_rank,string in_file_name)
{
    int result      = 0;                    //gives the rusult of reading the next string from a line
    int line_offset = 0;                    //stores the position in the current line
    char line[200];                         //stores a line
    char my_string[200];                    //stores the current string

    //scan over the first frame and get the number of atoms
    while(fgets(line, sizeof line, *in_file) != NULL)
    {
        line_offset = 0;
        result = next_string(world_rank,200,line,my_string,200,&line_offset);

        if(strcmp(my_string, "ATOM") == 0) //atom
        {
            *num_atoms = *num_atoms + 1;
        }
        else if(strcmp(my_string, "ENDMDL") == 0) //end of frame
        {
            break;
        }
    }
    rewind(*in_file);

    //count the frames 
    while(fgets(line, sizeof line, *in_file) != NULL)
    {
        line_offset = 0;
        result = next_string(world_rank,200,line,my_string,200,&line_offset);
        if(strcmp(my_string, "ENDMDL") == 0) //end of frame
        {
            *frames = *frames +1;
        }
    }

    //rewind file for future reading
    rewind(*in_file);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function analyzes a xtc file. It gets the number of atoms in the system, determines how many frames  //
// are in the file, and stores the starting position of each.                                                //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void analyze_xtc_file(int *num_atoms,int *frames,int world_rank,char* filename,vector<int64_t> &pos,
                      int *step,real *time,real *prec,int *status_xtc,int64_t *filesize)
{
    if(world_rank == 0)
    {
        int i = 0;          //standard variable used in loops
        rvec *x;            //holds the atomic coordinates for each frame
        matrix box_xtc;     //holds the box dimensions for each frame
        XDRFILE *xd;        //used for opening xtc files

        printf("Analyzing %s. \n",filename);

        //count number of atoms
        *status_xtc = read_xtc_natoms(filename,num_atoms);

        //allocate memory to hold coordinates
        x = (rvec *)calloc(*num_atoms , sizeof(*x));         

        //open xtc file
        xd = xdrfile_open(filename, "r");
        if(NULL == xd)
        {
            printf("Could not open trajectory file (%s). Check that the file exists. \n",filename);
        }
        else
        {
            //get position of fist frame
            pos.push_back(xdr_tell(xd));

            do 
            {
                *status_xtc = read_xtc(xd, *num_atoms, step, time, box_xtc, x, prec); //read a frame
                if (exdrENDOFFILE != *status_xtc) 
                {
                    //store position of frame
                    pos.push_back(xdr_tell(xd));

                    //count frame
                    (*frames)++;

                    //track the file size
                    *filesize = xdr_tell(xd);
                }
            } while (*status_xtc == exdrOK);

            xdrfile_close(xd);

            printf("Finished analyzing %s. \n",filename);
            printf("Trajectory frames: %-20d \n\n",*frames);
        }
        free(x);
    }

    //broadcast num_atoms, frames, and filesize
    MPI_Bcast(num_atoms, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(frames, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(filesize, 1, MPI_LONG, 0, MPI_COMM_WORLD);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function analyzes a trr file. It gets the number of atoms in the system, determines how many frames  //
// are in the file, and stores the starting position of each.                                                //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void analyze_trr_file(int *num_atoms,int *frames,int world_rank,char* filename,vector<int64_t> &pos,
                      int *step,real *time,int *status_trr,int64_t *filesize)
{
    if(world_rank == 0)
    {
        int i = 0;          //standard variable used in loops
        rvec *x;            //holds the atomic coordinates for each frame
        rvec *v;            //holds the velocities for each frame
        rvec *f;            //holds the forces for each frame
        matrix box_trr;     //holds the box dimensions for each frame
        XDRFILE *xd;        //used for opening trr files
        float lambda;       //holds trr lambda value

        printf("Analyzing %s. \n",filename);

        //count number of atoms
        *status_trr = read_trr_natoms(filename,num_atoms);

        //allocate memory to hold coordinates
        x = (rvec *)calloc(*num_atoms , sizeof(*x));
        v = (rvec *)calloc(*num_atoms , sizeof(*v));
        f = (rvec *)calloc(*num_atoms , sizeof(*f));

        //open trr file
        xd = xdrfile_open(filename, "r");
        if(NULL == xd)
        {
            printf("Could not open trajectory file (%s). Check that the file exists. \n",filename);
        }
        else
        {
            //get position of fist frame
            pos.push_back(xdr_tell(xd));

            do
            {
                int bBox_junk,bv_junk,bf_junk;  //needed to tell if a box, V, and F were present
                *status_trr = read_trr(xd, *num_atoms, step, time, &lambda, box_trr, x, v, f,&bBox_junk,&bv_junk,&bf_junk); //read frame
                if (exdrENDOFFILE != *status_trr)
                {
                    //store position of frame
                    pos.push_back(xdr_tell(xd));

                    //count frame
                    (*frames)++;

                    //track the file size
                    *filesize = xdr_tell(xd);
                }
            } while (*status_trr == exdrOK);

            xdrfile_close(xd);

            printf("Finished analyzing %s. \n",filename);
            printf("Trajectory frames: %-20d \n\n",*frames);
        }
        free(x);
        free(v);
        free(f);
    }

    //broadcast num_atoms, frames, and filesize
    MPI_Bcast(num_atoms, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(frames, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(filesize, 1, MPI_LONG, 0, MPI_COMM_WORLD);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function sets the file position for each mpi process when reading the trajectory. Set to the start   //
// of the next frame.                                                                                        //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void set_file_pos(int in_f,FILE *in_file,XDRFILE *xd_r,XDRFILE *trr_r,vector <int64_t> &pos_block,
                  int current_frame,int64_t *offset_r,vector <int64_t> &pos_stride,enum Switch block_parallel)
{
    //set the position of the current frame. This is done differently depending on the parallelization scheme. 
    if(block_parallel == on)
    {
        *offset_r = pos_block[current_frame];
    }
    else if(block_parallel == off)
    {
        *offset_r = pos_stride[current_frame];
    }

    //move to the current frame in the trajectory file
    if(in_f == 0) //gro
    {
        fseek(in_file,*offset_r,SEEK_SET);
    }
    else if(in_f == 1) //pdb
    {
        fseek(in_file,*offset_r,SEEK_SET);
    }
    else if(in_f == 2) //xtc
    {
        int seek_result = xdr_seek(xd_r, *offset_r, SEEK_SET);
    }
    else if(in_f == 3) //trr
    {
        int seek_result = xdr_seek(trr_r, *offset_r, SEEK_SET);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Returns the global frame for a mpi rank. Not the true global frame but the global frame considering only  //
// effective frames taking into acount stride and begin/end frame.                                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int get_global_frame(vector <int> &world_frames,int world_rank,int current_frame,int block_parallel)
{
    int i            = 0;       //standard variable used in loops
    int global_frame = 0;       //the global frame

    if(block_parallel == on)
    {
        for(i=0; i<world_rank; i++)
        {
            global_frame = global_frame + world_frames[i];
        }
        global_frame = global_frame + current_frame;
    }
    else
    {
        global_frame = current_frame;
    }

    return global_frame;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function extracts the current time and step from the title info                                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_step_and_time(char title[],int world_rank,real *time,int *step)
{
    char string200[200];     //store string extracted from a line
    int offset    = 0;       //offset in the current line
    int next_time = 0;       //tells if the next string is the time
    int next_step = 0;       //tells if the next string is the step
    int result    = 0;       //gives the result of reading a string from the line
    int i         = 0;       //standard variable used in loops

    while(i<1) //extract strings until i is set to 1 
    {
        init_carray(string200,200);
        result = next_string(world_rank,200,title,string200,200,&offset);

        //check if current string is the time or step
        if(next_time == 1)
        {
            *time = atof(string200);
        }
        else if(next_step == 1)
        {
            *step = atoi(string200);
            i = 1;
        }

        //check if the next string is the time or step
        if(strcmp(string200, "t=") == 0 || strcmp(string200, "time") == 0)
        {
            next_time = 1;
            next_step = 0;
        }
        else if(strcmp(string200, "step=") == 0 || strcmp(string200, "step") == 0)
        {
            next_time = 0;
            next_step = 1;
        }
        else
        {
            next_time = 0;
            next_step = 0;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads a frame from a gro file char by char.                                                 //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void read_gro_frame_by_char(FILE **in_file,matrix box,int *num_atoms,vector<int> &atom_nr,vector<int> &res_nr,vector<string> &res_name,
                            vector<string> &atom_name,rvec *r,rvec *v,char title[],int world_rank,real *time,int *step,int frames,gmx_bool *bV)
{
    int i               = 0;               //standard variable used in loops      
    int j               = 0;               //standard variable used in loops 
    int k               = 0;               //standard variable used in loops 
    int lines_per_frame = *num_atoms + 3;  //the lines per frame
    char string9[9];                       //store a string from current line
    char line[200];                        //holds the current line

    for(i=0; i<lines_per_frame; i++)
    {
        //The first line is unformatted text
        if(i == 0)
        {
            init_carray(line,200);
            init_carray(title,200);
            fgets(line, sizeof line, *in_file);

            for(k=0; k<200; k++)
            {
                title[k] = line[k];
            }
            if(frames > 1) //single frame gro files may not have time/step
            {
                get_step_and_time(title,world_rank,time,step);
            }
        }
        else if(i == 1) //Unformatted line containing number of atoms
        {
            init_carray(line,200);
            fgets(line, sizeof line, *in_file);
            *num_atoms = atoi(line);
        }
        else if(i == lines_per_frame-1) //Unformatted line containing box dimensions
        {
            init_carray(line,200);
            fgets(line, sizeof line, *in_file);

            double x1 = 0;  //holds one of the box coordinates
            double x2 = 0;  //holds one of the box coordinates
            double y1 = 0;  //holds one of the box coordinates
            double y2 = 0;  //holds one of the box coordinates
            double z1 = 0;  //holds one of the box coordinates
            double z2 = 0;  //holds one of the box coordinates

            sscanf(line, "%lf%lf%lf", &x1, &y1, &z1);
            box[XX][XX] = x1;
            box[YY][YY] = y1;
            box[ZZ][ZZ] = z1;

            if (sscanf (line, "%*f%*f%*f%lf%lf%lf%lf%lf%lf",
                        &x1, &y1, &z1, &x2, &y2, &z2) != 6)
            {
                x1 = y1 = z1 = x2 = y2 = z2 = 0.0;
            }
            box[XX][YY] = x1;
            box[XX][ZZ] = y1;
            box[YY][XX] = z1;
            box[YY][ZZ] = x2;
            box[ZZ][XX] = y2;
            box[ZZ][YY] = z2;

        }
        else //The remaining lines are formatted and contain the atoms/residues and their coords
        {
            int result = 0; 
            fgets(line, sizeof line, *in_file);
 
            int line_pos = 0;
            for(j=0; j<8; j++)
            {
                init_carray(string9,9);

                if(j==0) //Residue number
                {
                    for(k=0; k<5; k++)
                    {
                        string9[k] = line[line_pos]; 
                        line_pos++;
                    }
                    clean_string(string9,9);
                    res_nr[i-2] = atoi(string9);
                }
                else if(j==1) //Residue name
                {
                    for(k=0; k<5; k++)
                    {
                        string9[k] = line[line_pos];
                        line_pos++;
                    }
                    clean_string(string9,9);
                    res_name[i-2] = string9;
                }
                else if(j==2) //Atom name
                {
                    for(k=0; k<5; k++)
                    {
                        string9[k] = line[line_pos];
                        line_pos++;
                    }
                    clean_string(string9,9);
                    atom_name[i-2] = string9;
                }
                else if(j==3) //Atom number
                {
                    for(k=0; k<5; k++)
                    {
                        string9[k] = line[line_pos];
                        line_pos++;
                    }
                    clean_string(string9,9);
                    atom_nr[i-2] = atoi(string9);
                }
                else if(j==4) //x coord
                {
                    result = next_string(world_rank,200,line,string9,9,&line_pos);
                    r[i-2][0] = atof(string9);
                }
                else if(j==5) //y coord
                {
                    result = next_string(world_rank,200,line,string9,9,&line_pos);;
                    r[i-2][1] = atof(string9);
                }
                else if(j==6) //z coord
                {
                    result = next_string(world_rank,200,line,string9,9,&line_pos);
                    r[i-2][2] = atof(string9);
                }
                else if(j==7) //new line \n
                {
                    if(line[line_pos] == '\n') //end of line. no velocities found
                    {
                        result = next_string(world_rank,200,line,string9,9,&line_pos);
                    }
                    else //velocities present
                    {
                        result = next_string(world_rank,200,line,string9,9,&line_pos);
                        v[i-2][0] = atof(string9);

                        result = next_string(world_rank,200,line,string9,9,&line_pos);
                        v[i-2][1] = atof(string9);

                        result = next_string(world_rank,200,line,string9,9,&line_pos);
                        v[i-2][2] = atof(string9);
                        
                        *bV = 1;
                    }
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads a frame from the pdb char by char.                                                    //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void read_pdb_frame_by_char(FILE **in_file,matrix box,vector<int> &atom_nr,vector<int> &res_nr,vector<string> &res_name,vector<string> &atom_name,
                            rvec *r,char title[],int world_rank,real *time,int *step,int frames,vector<double> &beta,vector<double> &weight,
                            vector<string> &element,vector<char> &chain_id,gmx_bool *bBox)
{
    int i           = 0;            //standard variable used in loops
    int line_offset = 0;            //holds the current position in the line
    int result      = 0;            //tells the result of reading a string from the line
    int atom_count  = 0;            //used to count the atoms as they are encountered
    int num_chars   = 0;            //counts the characters in the current line
    char line[200];                 //holds the current line
    char string9[9];                //holds the current string to be stored i.e. the atom names etc.  
    char my_string[200];            //holds the first string of each row when searching for "ATOM" etc. 

    //we assume no box until the "CRYST1" line is found
    *bBox = 0;

    //read until "ENDMDL" is found. i.e. a single frame
    while(fgets(line, sizeof line, *in_file) != NULL)
    {
        line_offset = 0;
        result = next_string(world_rank,200,line,my_string,200,&line_offset);

        if(strcmp(my_string, "ATOM") == 0) //atom
        {
            //count the number of characters in the current line
            num_chars = 0;
            for(i=0; i<200; i++)
            {
                if(line[i] == '\n')
                {
                    break;
                }
                num_chars++;
            }

            //extract info from the line
            line_offset = 0;
            for(i=0; i<6; i++,line_offset++) //tag
            {
            }

            for(i=0; i<5; i++,line_offset++) //atom number
            {
                string9[i] = line[line_offset];
            }
            string9[5] = '\0';
            clean_string(string9,9);
            atom_nr[atom_count] = atoi(string9);
            line_offset++;

            for(i=0; i<4; i++,line_offset++) //atom name
            {
                string9[i] = line[line_offset];
            }
            string9[4] = '\0';
            clean_string(string9,9);
            atom_name[atom_count] = string9;
            line_offset++;

            for(i=0; i<4; i++,line_offset++) //res name
            {
                string9[i] = line[line_offset];
            }
            string9[4] = '\0';
            clean_string(string9,9);
            res_name[atom_count] = string9;

            for(i=0; i<1; i++,line_offset++) //chain identifier
            {
                string9[i] = line[line_offset];
            }
            string9[1] = '\0';
            chain_id[atom_count] = string9[0];

            for(i=0; i<4; i++,line_offset++) //res number
            {
                string9[i] = line[line_offset];
            }
            string9[4] = '\0';
            clean_string(string9,9);
            res_nr[atom_count] = atoi(string9);

            for(i=0; i<4; i++,line_offset++) //space
            {
            }

            for(i=0; i<8; i++,line_offset++) //x
            {
                string9[i] = line[line_offset];
            }
            string9[8] = '\0';
            clean_string(string9,9);
            r[atom_count][0] = 0.1*atof(string9);

            for(i=0; i<8; i++,line_offset++) //y
            {
                string9[i] = line[line_offset];
            }
            string9[8] = '\0'; 
            clean_string(string9,9);
            r[atom_count][1] = 0.1*atof(string9);
           
            for(i=0; i<8; i++,line_offset++) //z
            {
                string9[i] = line[line_offset];
            }
            string9[8] = '\0';
            clean_string(string9,9);
            r[atom_count][2] = 0.1*atof(string9);

            for(i=0; i<6; i++,line_offset++) //weight
            {
                string9[i] = line[line_offset];
            }
            string9[6] = '\0';
            clean_string(string9,9);
            weight[atom_count] = atof(string9);

            for(i=0; i<7; i++,line_offset++) //beta
            {
                string9[i] = line[line_offset];
            }
            string9[7] = '\0';
            clean_string(string9,9);
            beta[atom_count] = atof(string9);

            for(i=0; i<10; i++,line_offset++) //space
            {
            }
            
            if(num_chars < 79)
            {
                element[atom_count] = "  ";
            }
            else
            {
                for(i=0; i<2; i++,line_offset++) //element
                {              
                    string9[i] = line[line_offset];
                }
                string9[2] = '\0'; 
                clean_string(string9,9);
                element[atom_count] = string9;
            }
            atom_count++;
        }
        else if(strcmp(my_string, "CRYST1") == 0) //box
        {
            int ePbc;
            read_cryst1(line, &ePbc, box);
            *bBox = 1;
        }
        else if(strcmp(my_string, "TITLE") == 0) // title
        {
            line_offset = 10;
            i=0;
            while(line[line_offset] != '\n')
            {
                title[i] = line[line_offset];
                line_offset++;
                i++; 
            }
            title[i] = '\n';
            title[i+1] = '\0';
            //printf("%s",title);

            if(frames > 1) //single frame gro files may not have time/step
            {  
                get_step_and_time(title,world_rank,time,step);
            }
        }
        else if(strcmp(my_string, "ENDMDL") == 0) //end of frame
        {
            break;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function takes the atom and residue numbers which may not be continuous and makes them so. Atoms and //
// residues are indexed to begin with number 1.                                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_cont_indices(int num_atoms,vector<int> &atom_nr,vector<int> &res_nr)
{
    int i        = 0;    //standard variable used in loops
    int prev_res = -1;   //stores the previous residue number
    int count    = 0;    //counts the residues as they are encountered

    for(i=0; i<num_atoms; i++)
    {
        //make the atom numbers continuous
        atom_nr[i] = i + 1;

        //make the residue numbers continuous
        if(res_nr[i] != prev_res) //new residue
        {
            prev_res = res_nr[i];
            count++;
            res_nr[i] = count;
        }
        else
        {
            res_nr[i] = count;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads the next frame from the trajectory. To enable cross io certain variables must be      //
// determined regardless the file io types. For example bBox, bF, bV, title, chain_id, weight, beta, and     //
// element  must be determined.                                                                              // 
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void read_frame(FILE **in_file,matrix box,int *num_atoms,vector<int> &atom_nr,vector<int> &res_nr,vector<string> &res_name,vector<string> &atom_name,
                rvec *r,char title[],int world_rank,real *time,int *step,int frames,int in_f,XDRFILE *xd_r,real *prec,int *status_xtc,
                int current_frame,XDRFILE *trr_r,rvec *v,rvec *f,int *status_trr,int *box_dimension,vector<double> &beta,vector<double> &weight,
                vector<string> &element,vector<char> &chain_id,vector <int> &world_frames,int *global_frame,gmx_bool *bBox,gmx_bool *bV,
                gmx_bool *bF,int64_t offset_r,Switch block_parallel)
{
    int i = 0;    //standard variable used in loops

    if(in_f == 0) //gro
    {
        //read in the trajectory frame (overrides ref file resid,atomid etc.)
        read_gro_frame_by_char(in_file,box,num_atoms,atom_nr,res_nr,res_name,
                               atom_name,r,v,title,world_rank,time,step,frames,bV);

        //make the atoms and residue id continuous
        get_cont_indices(*num_atoms,atom_nr,res_nr);

        //set values needed to write data to pdb format
        for(i=0; i<(*num_atoms); i++)
        {
            chain_id[i] = ' ';
            weight[i]   = 0.0;
            beta[i]     = 1.0;
            element[i]  = "  "; 
        } 

        //get the global frame number
        *global_frame = get_global_frame(world_frames,world_rank,current_frame,block_parallel);

        //check if some items were included like box, velocities and forces
        *bBox = 1; //gro file always has a box line
        *bF = 0;   //gro file does not contain forces
    }
    else if(in_f == 1) //pdb
    {
        //read in the trajectory frame
        read_pdb_frame_by_char(in_file,box,atom_nr,res_nr,res_name,
                               atom_name,r,title,world_rank,time,step,frames,
                               beta,weight,element,chain_id,bBox);

        //make the atoms and residue id continuous (overrides ref file resid,atomid etc.)
        get_cont_indices(*num_atoms,atom_nr,res_nr);

        //get the global frame number
        *global_frame = get_global_frame(world_frames,world_rank,current_frame,block_parallel);

        //check if some items were included like box, velocities and forces
        *bV = 0;   //pdbs do not contain velocities
        *bF = 0;   //pdbs do not contain forces
    }
    else if(in_f == 2) //xtc
    {
        //set box to 0.0
        clear_mat(box);

        //read in the trajectory frame
        *status_xtc = read_xtc(xd_r, *num_atoms, step, time, box, r, prec);

        if(*status_xtc != 0) //could not read the frame
        {
            //note: we could get a cleaner exit here if we were willing to add a communication step. 
            printf("An MPI core with world rank %d was unable to read trajectory frame %d. Perhaps there is something wrong with the trajectory or .info file? Attempting to terminate analysis. \n",world_rank,get_global_frame(world_frames,world_rank,current_frame,block_parallel));
            exit(EXIT_SUCCESS);
        }

        //set values needed to write data to pdb format
        for(i=0; i<(*num_atoms); i++)
        {   
            chain_id[i] = ' ';
            weight[i]   = 0.0;
            beta[i]     = 1.0;  
            element[i]  = "  ";
        } 

        //get the global frame number
        *global_frame = get_global_frame(world_frames,world_rank,current_frame,block_parallel);

        //generate a "title"
        string title_s = "step "; 
        title_s = title_s + to_string(*step) + "    time " + to_string(*time) + "\n";
        strcpy(title, title_s.c_str());

        //check if some items were included like box, velocities and forces
        if(box[XX][XX] == 0 && box[XX][YY] == 0 && box[XX][ZZ] == 0 &&
           box[YY][XX] == 0 && box[YY][YY] == 0 && box[YY][ZZ] == 0 &&
           box[ZZ][XX] == 0 && box[ZZ][YY] == 0 && box[ZZ][ZZ] == 0)
        {
            *bBox = 0;
        }
        else
        {
            *bBox = 1;
        }
        *bV = 0;   //xtc file does not contain velcocities
        *bF = 0;   //xtc file does not contain forces
    }   
    else if(in_f == 3) //trr
    {
        float lambda;  //needed to read trr files

        //read in the trajectory frame
        *status_trr = read_trr(trr_r, *num_atoms, step, time, &lambda, box, r, v, f,bBox,bV,bF);

        if(*status_trr != 0) //could not read the frame
        {
            //note: we could get a cleaner exit here if we were willing to add a communication step. 
            printf("An MPI core with world rank %d was unable to read trajectory frame %d. Perhaps there is something wrong with the trajectory or .info file? Attempting to terminate analysis. \n",world_rank,get_global_frame(world_frames,world_rank,current_frame,block_parallel));
            exit(EXIT_SUCCESS);
        }

        //set values needed to write data to pdb format
        for(i=0; i<(*num_atoms); i++)
        {
            chain_id[i] = ' ';
            weight[i]   = 0.0;
            beta[i]     = 1.0;
            element[i]  = "  ";
        }

        //get the global frame number
        *global_frame  = get_global_frame(world_frames,world_rank,current_frame,block_parallel);

        //generate a "title"
        string title_s = "step ";
        title_s = title_s + to_string(*step) + "    time " + to_string(*time) + "\n";
        strcpy(title, title_s.c_str());
    } 

    //get the box dimension
    if(box[XX][YY] == 0.0 && box[XX][ZZ] == 0.0 && box[YY][XX] == 0.0 && box[YY][ZZ] == 0.0 && box[ZZ][XX] == 0.0 && box[ZZ][YY] == 0.0)
    {
        *box_dimension = 3;
    }
    else 
    {
        *box_dimension = 9;
    }

    //set the precision (needed if out_f is xtc)
    *prec = 1000.0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints the current gro frame data to the screen. Does not account for the max atom/res      //
// numbers (99999/99999) allowed for a gro file.                                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void print_gro_frame(matrix box,int num_atoms,vector<int> &atom_nr,vector<int> &res_nr,vector<string> &res_name,vector<string> &atom_name,
                     rvec *r,char *title,int world_rank,int box_dimension)
{
    int i               = 0;              //standard variable used in loops
    int j               = 0;              //standard variable used in loops
    int k               = 0;              //standard variable used in loops
    int lines_per_frame = num_atoms + 3;  //standard variable used in loops

    for(i=0; i<lines_per_frame; i++)
    {
        //The first line is text, the second gives number of atoms and the last gives the box
        if(i == 0) //title
        {
            printf("%s",title);
        }
        else if(i == 1) //atoms count
        {
            printf("%5d \n",num_atoms);
        }
        else if(i == lines_per_frame-1) //box
        {
            if(box_dimension == 3)
            {
                printf("   %7.5f   %7.5f   %7.5f\n",box[XX][XX],box[YY][YY],box[ZZ][ZZ]);
            }
            else if(box_dimension == 9)
            {
                for(j=0; j<3; j++)
                {
                    for(k=0; k<3; k++)
                    {
                        printf("   %7.5f",box[j][k]);
                    }
                }
                printf("\n");
            }
        }
        else //atom
        {
            printf("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",res_nr[i-2],res_name[i-2].c_str(),atom_name[i-2].c_str(),atom_nr[i-2],r[i-2][0],r[i-2][1],r[i-2][2]);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes out a trajectory frame in gro fromat. This function adjusts the atom and res id to   // 
// account for the maximum value allowed (99999/99999) for these numbers in a gro file.                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_frame_gro(matrix box,int num_atoms,vector<int> &atom_nr,vector<int> &res_nr,vector<string> &res_name,vector<string> &atom_name,
                     rvec *r,char *title,int world_rank,FILE **out_file,int box_dimension,gmx_bool bV,rvec *v)
{
    int i                = 0;              //standard variable used in loops
    int j                = 0;              //standard variable used in loops
    int k                = 0;              //standard variable used in loops
    int adjusted_atom_nr = 0;              //adjusted atom id (99999 max)
    int adjusted_res_nr  = 0;              //adjusted res id (99999 max)

    for(i=0; i<num_atoms+3; i++) //loop over the current frame
    {
        if(i >= 2)
        {
            //account for the maxium allowed resid and atom id
            adjusted_atom_nr = atom_nr[i-2]%100000;
            adjusted_res_nr  = res_nr[i-2]%100000;
        }

        //The first line is text, the second gives number of atoms and the last gives the box
        if(i == 0) //title
        {
            fprintf(*out_file,"%s",title);
        }
        else if(i == 1) //atom count 
        {
            fprintf(*out_file,"%d\n",num_atoms);
        }
        else if(i == num_atoms+2) //box
        {
            if(box_dimension == 3)
            {
                fprintf(*out_file,"  %8.5f  %8.5f  %8.5f\n",box[XX][XX],box[YY][YY],box[ZZ][ZZ]);
            }
            else if(box_dimension == 9)
            {
                fprintf(*out_file,"  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n",box[XX][XX],box[XX][YY],box[XX][ZZ],box[YY][XX],box[YY][YY],box[YY][ZZ],box[ZZ][XX],box[ZZ][YY],box[ZZ][ZZ]);
            }
        }
        else //atom
        {
            if(bV == 0) //contains no velocities
            {
                fprintf(*out_file,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",adjusted_res_nr,res_name[i-2].c_str(),atom_name[i-2].c_str(),adjusted_atom_nr,r[i-2][0],r[i-2][1],r[i-2][2]);
            }
            else //velocities present
            {
                fprintf(*out_file,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",adjusted_res_nr,res_name[i-2].c_str(),atom_name[i-2].c_str(),adjusted_atom_nr,r[i-2][0],r[i-2][1],r[i-2][2],v[i-2][0],v[i-2][1],v[i-2][2]);
            }
        }
    }
    fflush(*out_file);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes out a trajectory frame in pdb fromat. The x,y,z units are converted to angstrom here //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_frame_pdb(matrix box,int num_atoms,vector<int> &atom_nr,vector<int> &res_nr,vector<string> &res_name,vector<string> &atom_name,
                     rvec *r,char *title,int world_rank,FILE **out_file,vector<double> &beta,vector<double> &weight,
                     vector<string> &element,vector<char> &chain_id,int global_frame)
{
    int i = 0;  //standared variable used in loops

    //print the title
    fprintf(*out_file,"%-6s    %s","TITLE",title);

    //print the box
    gmx_write_pdb_box(*out_file,-1,box);

    //print the model number
    fprintf(*out_file,"MODEL%9d\n",global_frame);    

    //print the atom lines
    for(i=0; i<num_atoms; i++)
    {
        fprintf_atomline_pdb(out_file,"ATOM",atom_nr[i],atom_name[i].c_str(),' ',res_name[i].c_str(),chain_id[i],res_nr[i],' ',10.0*r[i][0],10.0*r[i][1],10.0*r[i][2],weight[i],beta[i],element[i].c_str());
    }

    //print the TER and ENDMDL tags
    fprintf(*out_file,"TER\n");
    fprintf(*out_file,"ENDMDL\n");
    fflush(*out_file);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes out a trajectory frame in xtc fromat.                                                //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_frame_xtc(int world_rank,XDRFILE *xd_w,int step,real time,matrix box, rvec *r, real prec,int num_atoms)
{
    int test = write_xtc(xd_w,num_atoms,step,time,box,r,prec);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes out a trajectory frame in trr fromat.                                                //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_frame_trr(int world_rank,XDRFILE *trr_w,int step,real time,matrix box, rvec *r,int num_atoms,
                     rvec *v,rvec *f,gmx_bool *bBox,gmx_bool *bV,gmx_bool *bF)
{   
    float lambda = 0.0; //needed when writing trr files

    write_trr(trr_w,num_atoms,step,time,lambda,
                       *bBox ? box : nullptr,
                       r,
                       *bV     ? v : nullptr,
                       *bF     ? f : nullptr,
                       bBox,bV,bF);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes the current trajectory frame to the output file. The function determines the output  //
// type and calls the appropriate helper function.                                                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_frame(matrix box,int num_atoms,vector<int> &atom_nr,vector<int> &res_nr,vector<string> &res_name,vector<string> &atom_name,rvec *r,
                 char *title,int world_rank,FILE **out_file,int out_f,XDRFILE *xd_w,int step,real time,
                 real prec,XDRFILE *trr_w,rvec *v,rvec *f,int box_dimension,vector<double> &beta,
                 vector<double> &weight,vector<string> &element,vector<char> &chain_id,int *global_frame,gmx_bool *bBox,gmx_bool *bV,
                 gmx_bool *bF,int b_print)
{
    if(b_print == 1)
    {
        if(out_f == 0) //gro
        {
            write_frame_gro(box,num_atoms,atom_nr,res_nr,res_name,atom_name,r,title,
                            world_rank,out_file,box_dimension,*bV,v);
        }
        else if(out_f == 1) //pdb
        {
            write_frame_pdb(box,num_atoms,atom_nr,res_nr,res_name,atom_name,r,title,
                      world_rank,out_file,beta,weight,element,chain_id,*global_frame);
        }
        else if(out_f == 2) //xtc
        {
            write_frame_xtc(world_rank,xd_w,step,time,box,r,prec,num_atoms);
        }
        else if(out_f == 3) //trr
        {
            write_frame_trr(world_rank,trr_w,step,time,box,r,num_atoms,v,f,bBox,bV,bF);   
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads in temporary trajectory files produced by each mpi rank and writes them to a single   //
// trajectory file. Also generates an info file for the new trajectory.                                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void finalize_traj(int world_rank,XDRFILE *xd_r,XDRFILE *xd_w,string out_file_name,int world_size,int *step,
                   real *time,real *prec,matrix box,int *status_xtc,rvec *r,int *num_atoms,XDRFILE *trr_r,
                   XDRFILE *trr_w,rvec *v,rvec *f,FILE **out_file,int out_f,int b_print,int *bBox,int *bV,int *bF,
                   vector <int> &world_frames)
{
    MPI_Barrier(MPI_COMM_WORLD); //make certain all ranks have finished reading traj before splicing it together

    if(world_rank == 0 && b_print == 1) 
    {  
        int i             = 0;    //standard variable used in loop   
        int world_size_ef = 0;    //how many mpi processes actually wrote a tmp traj file

        //create items needed to create an info file for the new trajectory
        vector <int64_t>  out_pos{};  //store the position of each frame written to .info file
        int64_t out_filesize = 0;     //This is the size of the output trajectory file in bits
        int     out_frames   = 0;     //This is the number of trajectory frames in the output file
  
        //figure out how many ranks actually read data
        for(i=0; i<world_size; i++)
        {
            if(world_frames[i] > 0)
            {
                world_size_ef++;
            }
            out_frames = out_frames + world_frames[i];
        }
   
        if(out_f == 0) //gro
        {
            printf("\nFinalizing trajectory. \n");

            char line[200];     //used to read in a line of the gro file
            char my_string[20]; //holds a string from the current line

            //open output file for writing 
            *out_file = fopen64(out_file_name.c_str(), "w");
            if(*out_file == NULL)
            {
                printf("failure opening %s. Make sure the file exists. \n",out_file_name.c_str());
            }

            //store first frame position
            out_pos.push_back(ftell(*out_file));

            for(i=0; i<world_size_ef; i++) //loop over temporary files
            {
                //get tmp file name
                string out_file_name_tmp = out_file_name;
                out_file_name_tmp.pop_back();
                out_file_name_tmp.pop_back();
                out_file_name_tmp.pop_back();
                out_file_name_tmp.pop_back();
                out_file_name_tmp = out_file_name_tmp + "_" + to_string(i) + ".gro";

                //report progress
                printf("  Working on %s \n",out_file_name_tmp.c_str());

                //open tmp file for reading
                FILE *out_file_tmp = fopen64(out_file_name_tmp.c_str(), "r");
                if(out_file_tmp == NULL)
                {
                    printf("failure opening %s. Make sure the file exists. \n",out_file_name_tmp.c_str());
                }

                int j = 0;
                while(fgets(line, sizeof line, out_file_tmp) != NULL) //read in the frames
                {
                    fprintf(*out_file,"%s",line); //write current frame to output file

                    if(j%(*num_atoms + 3) == (*num_atoms + 2) )
                    {
                        out_pos.push_back(ftell(*out_file));
                        out_filesize = ftell(*out_file);
                    }
                    j++;
                }
                fclose(out_file_tmp);

                //remove tmp file
                remove(out_file_name_tmp.c_str());
            }
            fclose(*out_file);
        }
        else if(out_f == 1) //pdb
        {
            printf("\nFinalizing trajectory. \n");

            char line[200];     //used to read in a line of the pdb file
            char my_string[20]; //holds a string from the current line

            //open output file for writing 
            *out_file = fopen64(out_file_name.c_str(), "w");
            if(*out_file == NULL)
            {
                printf("failure opening %s. Make sure the file exists. \n",out_file_name.c_str());
            }

            //store first frame position
            out_pos.push_back(ftell(*out_file));

            for(i=0; i<world_size_ef; i++) //loop over temporary files
            {
                //get tmp file name
                string out_file_name_tmp = out_file_name;
                out_file_name_tmp.pop_back();
                out_file_name_tmp.pop_back();
                out_file_name_tmp.pop_back();
                out_file_name_tmp.pop_back();
                out_file_name_tmp = out_file_name_tmp + "_" + to_string(i) + ".pdb";

                //report progress
                printf("  Working on %s \n",out_file_name_tmp.c_str());

                //open tmp file for reading
                FILE *out_file_tmp = fopen64(out_file_name_tmp.c_str(), "r");
                if(out_file_tmp == NULL)
                {
                    printf("failure opening %s. Make sure the file exists. \n",out_file_name_tmp.c_str());
                }

                int j = 0;
                while(fgets(line, sizeof line, out_file_tmp) != NULL) //read in the frames
                {
                    fprintf(*out_file,"%s",line); //write current frame to output file

                    int line_offset = 0;
                    int result      = next_string(world_rank,20,line,my_string,20,&line_offset);
                    if(strcmp(my_string, "ENDMDL") == 0) //end of frame
                    {
                        out_pos.push_back(ftell(*out_file));
                        out_filesize = ftell(*out_file);
                    }
                }
                fclose(out_file_tmp);
 
                //remove tmp file
                remove(out_file_name_tmp.c_str());
            }     
            fclose(*out_file);  
        }
        else if(out_f == 2) //xtc
        {
            printf("\nFinalizing trajectory. \n");

            //open output file for writing 
            xd_w = xdrfile_open(cname(out_file_name), "w");

            for(i=0; i<world_size_ef; i++)  //loop over temporary files
            {
                //get tmp file name
                string out_file_name_tmp = out_file_name;
                out_file_name_tmp.pop_back();
                out_file_name_tmp.pop_back();
                out_file_name_tmp.pop_back();
                out_file_name_tmp.pop_back();
                out_file_name_tmp = out_file_name_tmp + "_" + to_string(i) + ".xtc";

                //report progress
                printf("  Working on %s \n",out_file_name_tmp.c_str());

                //open temp file for reading
                xd_r = xdrfile_open(cname(out_file_name_tmp), "r");

                //allocate memory to hold coords
                rvec *tmp_r;
                tmp_r = (rvec *)calloc(*num_atoms , sizeof(*tmp_r));

                while(read_xtc(xd_r, *num_atoms, step, time, box, tmp_r, prec) == exdrOK) //read in the frames
                {
                    out_pos.push_back(xdr_tell(xd_w));

                    *status_xtc = write_xtc(xd_w,*num_atoms,*step,*time,box,tmp_r,*prec); //write current frame to output file

                    out_filesize = xdr_tell(xd_w);
                }

                //close temp file
                xdrfile_close(xd_r);

                //remove temporary xtc file when finished with it
                remove(out_file_name_tmp.c_str());

                //free up memory;
                free(tmp_r);
            }

            //close output file
            xdrfile_close(xd_w);
        } 
        else if(out_f == 3) //trr
        {
            printf("\nFinalizing trajectory. \n");

            float lambda;   //needed to read/write trr files
            rvec *tmp_r;    //holds the temporary coords
            rvec *tmp_v;    //holds the temporary vel 
            rvec *tmp_f;    //holds the temporary forces 

            //allocate memory to hold coords etc.   
            tmp_r = (rvec *)calloc(*num_atoms , sizeof(*tmp_r));
            if(*bV == 1)
            {
                tmp_v = (rvec *)calloc(*num_atoms , sizeof(*tmp_v));
            }
            if(*bF == 1)
            {
                tmp_f = (rvec *)calloc(*num_atoms , sizeof(*tmp_f));
            }

            //open output file for writing 
            trr_w = xdrfile_open(cname(out_file_name), "w");

            for(i=0; i<world_size_ef; i++) //loop over temporary files
            {  
                //get tmp file name 
                string out_file_name_tmp = out_file_name;
                out_file_name_tmp.pop_back();
                out_file_name_tmp.pop_back();
                out_file_name_tmp.pop_back();
                out_file_name_tmp.pop_back();
                out_file_name_tmp = out_file_name_tmp + "_" + to_string(i) + ".trr";

                //report progress
                printf("  Working on %s \n",out_file_name_tmp.c_str());

                //open temp file for reading
                trr_r = xdrfile_open(cname(out_file_name_tmp), "r");
                
                while(read_trr(trr_r, *num_atoms, step, time, &lambda, box, tmp_r, tmp_v, tmp_f,bBox,bV,bF) == exdrOK) //read in the frames
                {
                    out_pos.push_back(xdr_tell(trr_w));

                    //write current frame to output file
                    write_trr(trr_w,*num_atoms,*step,*time,lambda,
                                       *bBox ? box : nullptr,
                                       tmp_r,
                                       *bV     ? tmp_v : nullptr,
                                       *bF     ? tmp_f : nullptr,
                                       bBox,bV,bF);

                    out_filesize = xdr_tell(trr_w); 
                }


                //close temp file
                xdrfile_close(trr_r);

                //remove temporary trr file when finished with it
                remove(out_file_name_tmp.c_str());
            }
 
            //close output file
            xdrfile_close(trr_w);
           
            //free memory 
            free(tmp_r);
            if(*bV == 1)
            {
                free(tmp_v);
            }
            if(*bF == 1)
            {
                free(tmp_f);
            }
        }

        //write info file for the new trajectory
        string info_file_name = out_file_name + ".info";
        printf("Writing info file %s \n\n",info_file_name.c_str());
        FILE *info_file = fopen(info_file_name.c_str(),"w");
        fprintf(info_file," %ld \n",out_filesize);
        fprintf(info_file," %d \n",*num_atoms);
        fprintf(info_file," %d \n",out_frames);
        for(i=0; i<out_frames; i++)
        {
            fprintf(info_file," %ld \n",out_pos[i]);
        }
        fclose(info_file);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function closes trajectory/ref files before ending the program                                       //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void close_files(int in_f,XDRFILE *xd_r,XDRFILE *xd_w,FILE *in_file,FILE *out_file,XDRFILE *trr_r,XDRFILE *trr_w,
                 int out_f,int b_print,int world_rank,int my_frames)
{
    if(b_print == 1 && my_frames > 0)
    {
        //close output trajectory files
        if(out_f == 0 || out_f == 1) //gro,pdb
        {
            fclose(out_file);
        }
        else if(out_f == 2) //xtc
        {   
            xdrfile_close(xd_w);
        }
        else if(out_f == 3) //trr
        {
            xdrfile_close(trr_w);
        }
    }

    //close input trajectory files
    if(in_f == 0 || in_f == 1) //xtc,pdb
    {
        fclose(in_file);
    }
    else if(in_f == 2) //xtc
    {
        xdrfile_close(xd_r);
    }
    else if(in_f == 3) //trr
    {
        xdrfile_close(trr_r);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function finds the first and last atom of each residue.                                             //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_first_last(int num_atoms,vector<int> &res_nr,vector <int> &res_start,vector <int> &res_end)
{
    int i           = 0;   //standard variable used in loops
    int current_res = -1;  //counts the residues as they are encountered

    for(i=0; i<num_atoms; i++)
    {
        if(res_nr[i] != current_res)
        {
            res_start[res_nr[i]-1] = i;
            res_end[res_nr[i]-1]   = i;
            current_res            = res_nr[i];
        }
        else
        {
            res_end[res_nr[i]-1] = i;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function fits the atom selection in lsq_index to the reference structure.                           //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void lsq_fit(int dimension,int b_lsq,int num_atoms,vector <int> lsq_index_vec,real mass_lsq[],
             real mass[],rvec *r_ref,rvec *r,int current_frame,vector <int> &atom_nr,rvec *lsq_shift,int world_rank,int world_size)
{
    if(b_lsq == 1)
    {
        //notes about the procedure: 
        //1. reset_x uses the index numbers to reference position in the mass array.
        //   since atom numbers start at 1 we must shift the index down by 1.
        //2. do_fit uses the masses to determine which atoms to include in the fit.
        //   any atom with mass 0 is excluded from the fit.
        //3. reset_x and do_fit take an array. We must therefore copy masses from our vec to an array.
        //4. to prevent a build up of error we do not shift the reference structure. 
        //   we instead make a copy of the reference and shift that.
        //5. we only take note of how much the structure is shifted on frame 0. 
        //   the trajectory is thus shifted to the same position each time after the rotation.
        //   this is equivalant to tracking the shift each frame if the protein was centered around the selection
        //   prior to running routine.

        int i=0;                  //Standard variable used in loops
        int j=0;                  //Standard variable used in loops
        double pre_rmsd = 0;      //RMSD before lsq fit
        double post_rmsd = 0;     //RMSD after lsq fit

        //make a copy of the ref structure so no changes accumulate
        rvec *r_ref_copy;
        r_ref_copy = (rvec *)calloc(num_atoms , sizeof(*r_ref_copy));
        for(i=0; i<num_atoms; i++)
        {
            r_ref_copy[i][0] = r_ref[i][0];
            r_ref_copy[i][1] = r_ref[i][1];
            r_ref_copy[i][2] = r_ref[i][2];
        }

        //gromacs functions take an index as an array. So we copy our vector to an array
        //we also shift the atom number down by 1 here. 
        int lsq_index[lsq_index_vec.size()];
        for(i=0; i<lsq_index_vec.size(); i++)
        {
            lsq_index[i] = lsq_index_vec[i] - 1;
        }

        //first the mass is set for atoms in lsq_index. atoms with mass zero are not fit
        if(current_frame == 0)
        {
            for(i=0; i<num_atoms; i++) //loop over system atoms
            {
                for(j=0; j<lsq_index_vec.size(); j++) //loop over index atoms
                {
                    if(lsq_index[j] + 1 == atom_nr[i]) //shifted down by 1. see notes above. 
                    {
                        mass_lsq[i] = mass[i];
                    }
                }
            }
        }

        //store a pre-shifted coord so the frame/ref can be restored
        if(current_frame == 0 && world_rank == 0)
        {
            lsq_shift[0][0] = r[0][0]; 
            lsq_shift[0][1] = r[0][1];
            lsq_shift[0][2] = r[0][2];
        }

        //shift the current frame so com is at the origin
        reset_x_ndim(dimension, lsq_index_vec.size(), lsq_index, num_atoms, nullptr,r, mass_lsq);

        //shift the reference structure so com is at the origin
        reset_x_ndim(dimension, lsq_index_vec.size(), lsq_index, num_atoms, nullptr,r_ref_copy, mass_lsq);

        //determine how much the system was shifted by
        if(current_frame == 0 && world_rank == 0)
        {
            lsq_shift[0][0] = lsq_shift[0][0] - r[0][0];
            lsq_shift[0][1] = lsq_shift[0][1] - r[0][1];
            lsq_shift[0][2] = lsq_shift[0][2] - r[0][2];
        }

	if(current_frame == 0) //note this will cause trouble if one of the cores has no frames to read! BAD PROGRAMMER BAD!
        {
            MPI_Bcast(lsq_shift, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }

        //calculate rmsd before fit
        pre_rmsd = rmsdev_ind(lsq_index_vec.size(), lsq_index,mass,r,r_ref_copy);

        //do the fitting
        do_fit_ndim(dimension, num_atoms, mass_lsq, r_ref_copy, r);

        //calculate rmsd after fit
        post_rmsd = rmsdev_ind(lsq_index_vec.size(), lsq_index,mass,r,r_ref_copy);
        //printf("pre_rmsd %f post_rmsd %f \n",pre_rmsd,post_rmsd);

        //move the center of mass back to original position
        for(i=0; i<num_atoms; i++)
        {
            r[i][0] = r[i][0] + lsq_shift[0][0];
            r[i][1] = r[i][1] + lsq_shift[0][1];
            r[i][2] = r[i][2] + lsq_shift[0][2];
        }
 
        //free up memory
        free(r_ref_copy);
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function fits the atom selection in lsq_index and leaves the center at the origin                   //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void this_lsq_fit(int dimension,int b_lsq,int num_atoms,Index &this_lsq_index,real mass_lsq[],
                  real mass[],rvec *r_ref,rvec *r,int current_frame,vector <int> &atom_nr,int world_rank)
{
    if(b_lsq == 1)
    {
        int i=0;     //Standard variable used in loops
        int j=0;     //Standard variable used in loops

        //make a copy of the ref structure so no changes accumulate
        rvec *r_ref_copy;
        r_ref_copy = (rvec *)calloc(num_atoms , sizeof(*r_ref_copy));
        for(i=0; i<num_atoms; i++)
        {
            r_ref_copy[i][0] = r_ref[i][0];
            r_ref_copy[i][1] = r_ref[i][1];
            r_ref_copy[i][2] = r_ref[i][2];
        }

        //gromacs functions take an index as an array. So we copy our vector to an array
        //we also shift the atom number down by 1 here. 
        int lsq_index[this_lsq_index.index_s.size()];
        for(i=0; i<this_lsq_index.index_s.size(); i++)
        {
            lsq_index[i] = this_lsq_index.index_i[i] - 1;
        }

        //first the mass is set for atoms in lsq_index. atoms with mass zero are not fit
        if(current_frame == 0)
        {
            for(i=0; i<num_atoms; i++) //loop over system atoms
            {
                for(j=0; j<this_lsq_index.index_s.size(); j++) //loop over index atoms
                {
                    if(lsq_index[j] + 1 == atom_nr[i]) //shifted down by 1. see notes above. 
                    {
                        mass_lsq[i] = mass[i];
                    }
                }
            }
        }

        //shift the current frame so com is at the origin
        reset_x_ndim(dimension, this_lsq_index.index_s.size(), lsq_index, num_atoms, nullptr,r, mass_lsq);

        //shift the reference structure so com is at the origin
        reset_x_ndim(dimension, this_lsq_index.index_s.size(), lsq_index, num_atoms, nullptr,r_ref_copy, mass_lsq);

        //do the fitting
        do_fit_ndim(dimension, num_atoms, mass_lsq, r_ref_copy, r);

        //free up memory
        free(r_ref_copy);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function fits moves the reference structure to the origin and leaves it there                       //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void reset_ref(int dimension,int b_lsq,int num_atoms,Index &this_lsq_index,real mass_lsq[],
               real mass[],rvec *r_ref,rvec *r,int current_frame,vector <int> &atom_nr,int world_rank)
{
    if(b_lsq == 1)
    {
        int i=0;     //Standard variable used in loops
        int j=0;     //Standard variable used in loops

        //gromacs functions take an index as an array. So we copy our vector to an array
        //we also shift the atom number down by 1 here. 
        int lsq_index[this_lsq_index.index_s.size()];
        for(i=0; i<this_lsq_index.index_s.size(); i++)
        {
            lsq_index[i] = this_lsq_index.index_i[i] - 1;
        }

        //first the mass is set for atoms in lsq_index. atoms with mass zero are not fit
        if(current_frame == 0)
        {
            for(i=0; i<num_atoms; i++) //loop over system atoms
            {
                for(j=0; j<this_lsq_index.index_s.size(); j++) //loop over index atoms
                {
                    if(lsq_index[j] + 1 == atom_nr[i]) //shifted down by 1. see notes above. 
                    {
                        mass_lsq[i] = mass[i];
                    }
                }
            }
        }

        //shift the reference structure (not a copy) so com is at the origin
        reset_x_ndim(dimension, this_lsq_index.index_s.size(), lsq_index, num_atoms, nullptr,r_ref, mass_lsq);

    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This is a class for reading/writing trajectory data                                                      //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Trajectory
{
    private:
        enum Switch block_parallel;                   //Controlls the parallelization scheme
        FILE *in_file;                                //File for working with the input trajectory file
        FILE *out_file;                               //File for writing the output trajectory file
        FILE *ref_file;                               //File for reading the reference file(if gro or pdb) 
        string out_file_name;                         //Name of the output trajectory file
        string out_file_name_tmp;                     //Name of the temporary out trajectory file
        string ref_file_name;                         //Name of the refernce file
        string lsq_index_file_name;                   //Name of the index file with atoms to do least squares fitting
        string traj_file_name;                        //Name of the input trajectory file
        int argc                 = 0;                 //How many command line arguments are there
        int world_size           = 0;                 //How many mpi process in the world
        int world_rank           = 0;                 //The rank for the mpi process
        int frames               = 0;                 //Total number of frames in the trajectory file
        int num_atoms            = 0;                 //Number of atoms in the system
        int num_atoms_ref        = 0;                 //Number of atoms from reference file
        int in_f                 = 0;                 //What is the format of the input file? (gro:0 pdb:1 xtc:2 trr:3)
        int out_f                = 0;                 //What is the format of the output file? (gro:0 pdb:1 xtc:2 trr:3)
        int ref_f                = 0;                 //What is the format of the reference file? (gro:0 pdb:1)
        //int box_dimension        = 9;                 //How many entries for the box? (3 or 9)
        int my_frames            = 0;                 //Number of frames the mpi process is responsible for (after -b, -e and -stride)
        int effective_frames     = 0;                 //Number of frames to be analyzed in total (summing over each mpi process and accounting for -b, -e and stride)
        int stride               = 1;                 //How many frames do we skip each time?
        int start_frame          = 0;                 //Dont read trajectory frames before this number
        int end_frame            = -1;                //Dont read trajectory frames after this number
        int b_print              = 0;                 //Print the output trajectory?
        int b_lsq                = 0;                 //Do least squares fitting?
        int lsq_dim              = 3;                 //Dimension of lsq fitting (2 or 3).
        int lsq_ref              = 0;                 //Structure used for lsq fitting (0:ref 1:frame_0)
        int num_residues         = 0;                 //How many residues in the system
        int64_t offset_r         = 0;                 //Possition in trajectory file where the mpi process should move to for reading
        int64_t filesize         = 0;                 //This is the size of the trajectory file in bits
        rvec *lsq_shift;                              //How much is the system translated by when doing least squares fitting?

        //xtc stuff
        XDRFILE *xd_r;                                //Gromacs structure for reading xtc files
        XDRFILE *xd_w;                                //Gromacs structure for writing xtc files
        real prec;                                    //Trajectory precision
        int status_xtc;                               //Used for reading a frame from an xtc file

        //trr stuf
        XDRFILE *trr_r;                               //Gromacs structure for reading trr files
        XDRFILE *trr_w;                               //Gromacs structure for writing trr files
        int status_trr;                               //Used for reading a frame from a trr file

        //vectors
        vector <int64_t>  pos{};                      //Stores start pos of each frame in the complete traj file
        vector <int64_t>  pos_stride{};               //Stores start pos of all frames to be analyzed after accounting for -stride, -b and -e
        vector <int64_t>  pos_block{};                //Stores start pos of frames for the current mpi process (after accounting for -stride, -b and -e and block parallel)
        vector <int>      world_frames{};             //Array holding how many frames each mpi process is responsible for (holds my_frames for each mpi process)
        vector <int>      world_frame_i{};            //stores the starting frame for each mpi process
        vector <int>      world_frame_f{};            //stores the final frame for each mpi process
        vector <int>      lsq_index{};                //Index with atoms to include in lsq fit
 
    public:
        int current_frame     = 0;                    //Current frame the mpi process is working on
        rvec *r;                                      //Holds coordinates
        rvec *v;                                      //Holds velocities
        rvec *f;                                      //Holds forces
        rvec *r_ref;                                  //Holds coordinates from reference file
        rvec *v_ref;                                  //Holds velocities from reference file
        rvec *r0;                                     //Holds the coordinates from first frame in traj
        rvec *r_store;                                //A vector in which the user can store the current coords for later use
        real time            = 0;                     //Holds the time of the current frame
        int step             = 0;                     //Holds the step of the current frame
        matrix box;                                   //The system box updated for each frame
        matrix ibox;                                  //The box from the first frame of the trajectory
        matrix box_ref;                               //The box from the refference file
        real *mass;                                   //Stores the mass of each atom in the system
        real *mass_lsq;                               //Stores masses used for least squares fitting
        char title[200];                              //System title. The first line of a gro file. Also added to pdb files
        vector <int>    atom_nr{};                    //Stores the atom number
        vector <int>    res_nr{};                     //Stores the residue number
        vector <string> res_name{};                   //Stores the reisdue name
        vector <string> atom_name{};                  //Stores the atom name
        vector <double> beta{};                       //Beta factor for pdb files
        vector <double> weight{};                     //Weight used for pdb files
        vector <string> element{};                    //Element collumn in pdb file
        vector <char>   chain_id{};                   //Chain id for pdb file
        vector <int>    res_start{};                  //First atom of the residue
        vector <int>    res_end{};                    //Last atom of the residue
        gmx_bool bBox            = 0;                 //Does the trajectory file contain a box?
        gmx_bool bF              = 0;                 //Does the trajectory file contain forces?
        gmx_bool bV              = 0;                 //Does the trajectory file contain velocities?
        int global_frame         = 0;                 //The frame number from the whole trajectory
        int box_dimension        = 9;                 //How many entries for the box? (3 or 9)

        //leaflet finder
        int             leaflet;                      //Which leaflet is the target leaflet?
        int             op_leaflet;                   //Which leaflet is the opposing leaflet?
        vector <int>    leaflets{};                   //Tells wheter each atom belongs to the upper, lower or neither leaflet
        vector <int>    target_leaflet{};             //The target leaflet 
        vector <int>    opposing_leaflet{};           //The opposing leaflet
        vector <int>    full_membrane{};              //The complete lipid bilayer
        vector <double> prev_z{};                     //The previous step coords (useful for detecting pbc jumps)

        //protein finder
        vector <int>    prot{};                       //The protein atoms
        vector <int>    b_protein{};                  //Tells if each atoms belongs to the protein 

        //solvent finder
        vector <int>    sol{};                        //Holds the solvent atoms
        vector <int>    b_sol{};                      //Tells if each atom belongs to the solvent
        int             num_waters;                   //How many solvent molecules are there?

        //parallelization schemes
        int             my_num_g_x = 0;               //How many grid points in x is the mpi process responsible for
        int             my_xi      = 0;               //The mpi processes first grid point 
        int             my_xf      = 0;               //The mpi processes last grid point
        vector <int>    world_num_g_x{};              //How many grid points in x is the mpi process responsible for (world)
        vector <int>    world_xi{};                   //The start grid point for each mpi process
        vector <int>    world_xf{};                   //The end grid point for each mpi process
        int             my_lipids   = 0;              //How many lipids is the mpi process responsible for
        int             lipid_start = 0;              //First lipid the mpi process is responsible for
        int             lipid_end   = 0;              //Last lipid the mpi process is responsible for
        vector <int>    world_lipids{};               //How many lipids is the mpi process responsible for (world)
        vector <int>    world_lipid_start{};          //The first lipid the mpi process is reponsible for (world)
        vector <int>    world_lipid_end{};            //The last lipid the mpi process is responsible for (world)
        int             my_prots     = 0;             //How many protein atoms is the mpi process responsible for
        int             prot_start   = 0;             //First protein atom the mpi process is responsible for
        int             prot_end     = 0;             //Last protein atom the mpi process is responsible for
        vector <int>    world_prots{};                //How many protein atoms is the mpi process responsible for (world)
        vector <int>    world_prot_start{};           //The first protein atom the mpi process is reponsible for (world)
        vector <int>    world_prot_end{};             //The last protein atom the mpi process is responsible for (world)
        int             my_aminos     = 0;            //How many protein residues is the mpi process responsible for
        int             amino_start   = 0;            //First protein residue the mpi process is responsible for
        int             amino_end     = 0;            //Last protein residue the mpi process is responsible for
        vector <int>    world_aminos{};               //How many protein residues is the mpi process responsible for (world)
        vector <int>    world_amino_start{};          //The first protein residue the mpi process is reponsible for (world)
        vector <int>    world_amino_end{};            //The last protein residue the mpi process is responsible for (world)
        int             my_waters   = 0;              //How many waters is the mpi process responsible for
        int             water_start = 0;              //First water the mpi process is responsible for
        int             water_end   = 0;              //Last water the mpi process is responsible for
        vector <int>    world_waters{};               //How many waters is the mpi process responsible for (world)
        vector <int>    world_water_start{};          //The first water the mpi process is reponsible for (world)
        vector <int>    world_water_end{};            //The last water the mpi process is responsible for (world)
        int             my_atoms   = 0;               //How many atoms is the mpi process responsible for
        int             atom_start = 0;               //First atom the mpi process is responsible for
        int             atom_end   = 0;               //Last atom the mpi process is responsible for
        vector <int>    world_atoms{};                //How many atoms is the mpi process responsible for (world)
        vector <int>    world_atom_start{};           //The first atom the mpi process is reponsible for (world)
        vector <int>    world_atom_end{};             //The last atom the mpi process is responsible for (world)
        int             my_num_g   = 0;               //How many grid points in xy is the mpi process responsible for
        int             my_gi      = 0;               //The mpi processes first grid point 
        int             my_gf      = 0;               //The mpi processes last grid point
        vector <int>    world_num_g{};                //How many grid points in xy is the mpi process responsible for (world)
        vector <int>    world_gi{};                   //The start grid point for each mpi process
        vector <int>    world_gf{};                   //The end grid point for each mpi process

    public:
        void        set_input_arguments(int argc, const char * argv[]);                                 //give the traj the command line args
        void        set_block_parallel(Switch setting);                                                 //set the parallelization scheme
        void        set_traj(string in_file_name);                                                      //set the trajectory file name
        void        set_ref(string my_ref_file_name);                                                   //set the reference file name
        void        set_traj_w(string my_out_file_name,int my_b_print);                                 //set the output traj file name
        void        set_lsq(string my_lsq_index_file_name,int my_b_lsq,int my_lsq_dim,int my_lsq_ref);  //set up the lsq fitting
        void        set_res(int my_stride,int my_start_frame,int my_end_frame);                         //set the stride and begin/end frames
        double      build();                                                                            //analyze the traj/ref files
        void        workload();                                                                         //print info about work load distribution
        int         get_frames();                                                                       //return the number of frames in the trajectory
        int         get_start_frame();                                                                  //return start_frame
        int         get_end_frame();                                                                    //return end_frame
        int         get_stride();                                                                       //return stride
        int         get_num_frames();                                                                   //return the number of frames when looping over traj
        vector<int> get_num_frames_world();                                                             //returns number of frames each mpi process is responsible for (world_frames)
        void        read_traj_frame();                                                                  //reads a frame from the traj
        void        do_fit();                                                                           //do least squares fitting on current frame
        void        do_this_fit(Index &this_lsq_index,int this_lsq_dim, int this_lsq_ref);              //do least squares fitting on current frame. leave center at the origin
        void        place_ref_at_origin(Index &this_lsq_index,int this_lsq_dim,int this_lsq_ref);       //moves com of ref structure to the origin and leaves it there
        void        write_traj_frame();                                                                 //write the current frame to output traj
        double      finalize_trajectory();                                                              //splice together temporary traj files
        int         atoms();                                                                            //returns the num_atoms in the trajectory
        int         get_ef_frames();                                                                    //return the effective num frames after -b, -e and -stride
        real        get_prec();                                                                         //return the precision of the trajectory
        int         get_frame_global();                                                                 //returns the global frame   
        int         get_frame_full();                                                                   //returns the frame index for the full trajectory (no -b -e -stride or block parallel)
        int         get_b_lsq();                                                                        //returns b_lsq
        int         get_box_dim();                                                                      //returns the box dimension
        int         get_num_residues();                                                                 //returns the number of residues
        void        report_progress();                                                                  //reports that analysis is beginning
        dv1d        center(sv1d &target_atoms,int min,int max);                                         //computes the center of a molecule
        dv1d        center_store(sv1d &target_atoms,int min,int max);                                   //computes the center of a molecule using stored coords
        dv1d        center_i(iv1d &group);                                                              //computes the center of a group of atoms from an index
        dv1d        com(sv1d &target_atoms,int min,int max);                                            //computes the center of mass of a group of atoms from an index
        dv1d        com_i(iv1d group);                                                                  //computes the center of mass of a group of atoms from an index
        double      gyrate(sv1d &target_atoms,int min,int max);                                         //computes the radius of gyration a group of atoms from an index
        double      dihedral_angle(int i,int j,int k,int l);                                            //computes the dihedral anle for an atom selection
        void        store_coords();                                                                     //copies the current coords to r_store
        int         get_res_start(int i);                                                               //returns the first atom of the current residue
        int         get_res_end(int i);                                                                 //returns the last atom of the current residue
        int         next_residue(int i);                                                                //primes loop index i (traj.atom_name[i]) for the next residue
        void        check_broken_molecule(int flag_beta);                                               //checks current frame for broken residues
        double      get_dist(dv1d &vec_a,dv1d &vec_b);                                                  //returns the distance between 2 vectors
        void        expand_traj(Trajectory traj_ref);                                                   //generates surrounding 8 periodic images for the reference traj
        dv2d        get_centers_target_lf();                                                            //returns the centers of the lipids in the target leaflet
        dv2d        get_centers_opposing_lf();                                                          //returns the centers of the lipids in the opposing leaflet
        dv2d        get_centers_full_mem();                                                             //returns the centers of the lipids in the full membrane
        dv2d        get_centers_prot();                                                                 //returns the center of the residues in the protein 
        dv2d        get_centers_sol();                                                                  //returns the center of the residues in the solvent
        dv2d        get_centers_system();                                                               //returns the center of the residues in the molecular system
        int         get_global_frame_i();                                                               //returns the global index of the first frame to be read by the core
        int         get_global_frame_f();                                                               //returns the global index of the last frame to be read by the core

        //leaflet finder 
        void        get_leaflets(int leaf,string leaflet_finder_param_name,int b_lf_param);             //assign atoms to the leaflets
        void        write_leaflets(string lf_pdb_file_name,int b_lf_pdb);                               //write pdb with leaflets indicated by beta factor
        void        set_beta_lf();                                                                      //set beta based on leaflets
        void        get_leaflet_stats();                                                                //print info about the leaflet and lipids
        void        get_lipid_selection_stats(vector <string> lip_t,string title);                      //print info about the selected lipids 
        int         count_target_lipids();                                                              //count the number of lipids in the target leaflet
        int         count_target_lipids_type(vector <string> lip_t);                                    //count the number of lipids of a type in the target leaflet
        int         count_opposing_lipids();                                                            //count the number of lipids in the opposing leaflet
        int         count_opposing_lipids_type(vector <string> lip_t);                                  //count the number of lipids of a type in the opposing leaflet
        int         t_lip_start(int i);                                                                 //returns the first atom (index) of the target lipid
        int         t_lip_end(int i);                                                                   //returns the last atom (index) of the target lipid
        int         o_lip_start(int i);                                                                 //returns the first atom (index) of the opposing lipid
        int         o_lip_end(int i);                                                                   //returns the last atom (index) of the opposing lipid
        int         fm_lip_start(int i);                                                                //returns the first atom (index) of the full membrane lipid
        int         fm_lip_end(int i);                                                                  //returns the last atom (index) of the full membrane lipid
        int         next_target_lipid(int i);                                                           //primes loop index (target_leaflet) for next lipid
        int         next_opposing_lipid(int i);                                                         //primes loop index (opposing_leaflet) for next lipid
        int         next_full_mem_lipid(int i);                                                         //primes loop index (full_membrane) for next lipid
        void        check_lip_jumps_z();                                                                //checks lipids for periodic boundary jumps in z-direction 

        //protein finder
        void        get_protein(string protein_finder_param_name,int b_pf_param);                       //get protein atoms
        void        write_protein(string pf_pdb_file_name,int b_pf_pdb);                                //print pdb with protein highlighted by beta value 
        void        get_prot_stats();                                                                   //print info about the protein
        int         get_num_res_prot();                                                                 //return the number of protein residues
        int         next_prot_res(int i);                                                               //primes loop index (prot) for the nex residue
        int         p_res_start(int i);                                                                 //returns the first atom (index) of the residue
        int         p_res_end(int i);                                                                   //returns the last atom (index) of the residue

        //solvent finder
        void        get_solvent(string solvent_finder_param_name,int b_sf_param);                       //get the solvent atoms
        void        write_sol(string sf_pdb_file_name,int b_sf_pdb);                                    //print pdb with solven atoms indicated by beta factor
        void        get_sol_stats();                                                                    //print info about the solvent
        int         sol_start(int i);                                                                   //returns the first atom (index) of the water   
        int         sol_end(int i);                                                                     //returns the last atom (index) of the water
        int         next_water(int i);                                                                  //primes loop index (sol) for next water
        int         count_waters_type(vector <string> sol_t);                                           //counts how many solvent molecules of a given type

        //parallelization schemes
        void        parallelize_by_grid(int num_g_x);                                                   //distribute workload by grid points (num_g_x)
        void        workload_grid();                                                                    //prints the workload distribution for grid points
        void        parallelize_by_lipid(int num_lipids_1);                                             //distribute workload by lipids
        void        workload_lipid();                                                                   //prints the workload distribution for lipids
        void        parallelize_by_prot(int num_prots_1);                                               //distribute workload by protein atoms
        void        workload_prot();                                                                    //prints the workload distribution for protein atoms
        void        parallelize_by_amino(int num_aminos_1);                                             //distribute workload by protein residues
        void        workload_amino();                                                                   //prints the workload distribution for protein residues
        void        parallelize_by_water(int num_waters_1);                                             //distribute workload by waters
        void        workload_water();                                                                   //prints the workload distribution for waters
        void        parallelize_by_selection(int num_atoms_1);                                          //distribute workload by a custom atom selection
        void        workload_sel();                                                                     //prints the workload distribution for custom atom selections
        void        parallelize_by_grid_alt(int num_g_x,int num_g_y);                                   //distribute the workload by the grid points (num_g_x and num_g_y)
        void        workload_grid_alt();                                                                //prints the workload distribution for grid point (num_g_x and num_g_y)

        const char **argv;                                                                              //The command line arguments
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function sets the command line args for the trajectory                                              //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::set_input_arguments(int my_argc, const char * my_argv[])
{
    int i = 0;    //standard variable used in loops

    argc = my_argc;

    for(i=0; i<argc; i++)
    {
        argv[i] = my_argv[i]; 
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns the global frame for the current frame                                             //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::get_frame_global()
{
    return get_global_frame(world_frames,world_rank,current_frame,block_parallel);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns the precision of the trajectory                                                    //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
real Trajectory::get_prec()
{
    return prec;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns the number of frames in the trajectory                                             //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::get_frames()
{
    return frames;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns the first trajectory frame to be analyzed                                          //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::get_start_frame()
{
    return start_frame;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns the last trajectory frame to be analyzed                                           //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::get_end_frame()
{
    return end_frame;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns the stride used for the analysis                                                   //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::get_stride()
{
    return stride;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns the effective number of frames. i.e. the total number of frames to be analyzed     //
// counting all mpi processes and accounting for -b, -e, and -stride.                                       //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::get_ef_frames()
{
    return effective_frames;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns the frame index for the full trajectory (no -b -e stride or block parallel)        //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::get_frame_full()
{
    int i                  = 0;     //standard variable used in loops
    int frame_number_full  = 0;     //frame index in the full trajectory
    int global_frame_count = 0;     //counts how many frames analyzed which is related to frame_number_full

    for(i=0; i<get_frames(); i++)
    {
        if(i >= get_start_frame() && i <= get_end_frame() && i%get_stride() == 0) //i is one of the frames analyzed
        {
            if(global_frame_count == get_frame_global())
            {
                frame_number_full = i;
            }
            global_frame_count++;
        }
    }

    return frame_number_full; 
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns b_lsq and is used to tell if the least squares fitting was used.                   //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::get_b_lsq()
{
    return b_lsq;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns the box dimension                                                                  //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::get_box_dim()
{
    return box_dimension;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns the number of residues                                                             //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::get_num_residues()
{
    return num_residues;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function sets the block parallelization scheme for the trajectory                                   //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::set_block_parallel(Switch setting)
{
    block_parallel = setting;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints to the screen the number of frames each mpi process is responsible for as well as    //
// the firs and last frame for each.                                                                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::workload()
{
    if(world_rank == 0)
    {
        if(block_parallel == on)
        {
            int i = 0;
            printf("Distributing the work load across %d mpi processes. \n",world_size);

            printf("%6s %12s %12s %12s \n","Rank","Frames","First_frame","Last_frame");
            printf("---------------------------------------------\n");
            for(i=0; i<world_size; i++)
            {
                printf("%6d %12d %12d %12d \n",i,world_frames[i],world_frame_i[i],world_frame_f[i]);
            }
            printf("\n");
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns the number of frames to be read by a given mpi process (used to loop over the traj)//
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::get_num_frames()
{
    return my_frames;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function reads in a frame from the trajectory                                                       //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::read_traj_frame()
{
    //adjust the file position
    set_file_pos(in_f,in_file,xd_r,trr_r,pos_block,current_frame,&offset_r,pos_stride,block_parallel);

    //read in the data
    read_frame(&in_file,box,&num_atoms,atom_nr,res_nr,res_name,atom_name,r,title,world_rank,&time,&step,frames,
               in_f,xd_r,&prec,&status_xtc,current_frame,trr_r,v,f,&status_trr,
               &box_dimension,beta,weight,element,chain_id,world_frames,&global_frame,&bBox,&bV,&bF,offset_r,block_parallel);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function does least squares fitting of the current frame                                            //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::do_fit()
{
    lsq_fit(lsq_dim,b_lsq,num_atoms,lsq_index,mass_lsq,mass,(lsq_ref == 0) ? r_ref : r0,r,current_frame,atom_nr,lsq_shift,world_rank,world_size);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function does least squares fitting of the current frame. leaves center at the origin.              //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::do_this_fit(Index &this_lsq_index,int this_lsq_dim, int this_lsq_ref)
{
    this_lsq_fit(this_lsq_dim,1,num_atoms,this_lsq_index,mass_lsq,mass,(this_lsq_ref == 0) ? r_ref : r0,r,current_frame,atom_nr,world_rank);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function moves the com of the reference structure to the origin. leaves center at the origin.       //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::place_ref_at_origin(Index &this_lsq_index,int this_lsq_dim,int this_lsq_ref)
{
    reset_ref(this_lsq_dim,1,num_atoms,this_lsq_index,mass_lsq,mass,(this_lsq_ref == 0) ? r_ref : r0,r,current_frame,atom_nr,world_rank);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function writes the current trajectory frame to the temporary output file                           //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::write_traj_frame()
{
    write_frame(box,num_atoms,atom_nr,res_nr,res_name,atom_name,r,title,world_rank,&out_file,out_f,xd_w,step,time,
                prec,trr_w,v,f,box_dimension,beta,weight,element,chain_id,&global_frame,&bBox,&bV,&bF,
                b_print);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function reads the temporary traj files and splices them together                                   //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Trajectory::finalize_trajectory()

{
    //take the time when finalize_traj started
    clock_t t_i = clock();

    //close the temporary output files that were being written to as well as the input trajectory files that were read 
    close_files(in_f,xd_r,xd_w,in_file,out_file,trr_r,trr_w,out_f,b_print,world_rank,my_frames);

    //read back the temp files and copy data to a single output trajectory file
    finalize_traj(world_rank,xd_r,xd_w,out_file_name,world_size,&step,&time,&prec,box,&status_xtc,r,
                  &num_atoms,trr_r,trr_w,v,f,&out_file,out_f,b_print,
                  &bBox,&bV,&bF,world_frames);

    //get the total time spent splicing files
    double time_fin = (clock() - t_i)/CLOCKS_PER_SEC;

    return time_fin;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function sets the trajectory file name and type                                                     //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::set_traj(string in_file_name)
{
    //get the world size and assign ranks
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);      //get the world size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);      //get the process rank

    //set the input trajectory file name
    traj_file_name = in_file_name;

    //here we check the input file format (gro:0 pdb:1 xtc:2)
    int length = traj_file_name.length();
    if(traj_file_name.at(length-4) == '.' && traj_file_name.at(length-3) == 'g' && traj_file_name.at(length-2) == 'r' && traj_file_name.at(length-1) == 'o') //gro
    {
        in_f = 0;
    }
    else if(traj_file_name.at(length-4) == '.' && traj_file_name.at(length-3) == 'p' && traj_file_name.at(length-2) == 'd' && traj_file_name.at(length-1) == 'b') //pdb
    {
        in_f = 1;
    }
    else if(traj_file_name.at(length-4) == '.' && traj_file_name.at(length-3) == 'x' && traj_file_name.at(length-2) == 't' && traj_file_name.at(length-1) == 'c') //xtc
    {
        in_f = 2;
    }
    else if(traj_file_name.at(length-4) == '.' && traj_file_name.at(length-3) == 't' && traj_file_name.at(length-2) == 'r' && traj_file_name.at(length-1) == 'r') //trr
    {
        in_f = 3;
    }
    else
    {
        if(world_rank == 0)
        {
            printf("Supported trajectory types are xtc, trr, gro and pdb. Make sure your input file name ends with one of these.\n");
        }
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }

    //check that file exists
    FILE *test_file = fopen64(traj_file_name.c_str(), "r");
    if(test_file == NULL)
    {
        if(world_rank == 0)
        {
            printf("failure opening %s. Make sure the file exists. \n",traj_file_name.c_str());
        }
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }
    else 
    {
        fclose(test_file);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function sets the reference file name and type                                                      //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::set_ref(string my_ref_file_name)
{
    //get the world size and assign ranks
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);      //get the world size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);      //get the process rank

    //set the reference file name
    ref_file_name = my_ref_file_name; 

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
            printf("Supported reference file types are pdb and gro. Make sure your reference file name ends with one of these.\n");
        }
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }

    //check that file exists
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
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function sets the output trajectory name and type and b_print                                       //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::set_traj_w(string my_out_file_name,int my_b_print)
{
    //get the world size and assign ranks
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);      //get the world size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);      //get the process rank

    //set the output trajectory file name and b_print
    out_file_name  = my_out_file_name;
    b_print = my_b_print;

    if(b_print == 1)
    {
        int length = out_file_name.length();
        if(out_file_name.at(length-4) == '.' && out_file_name.at(length-3) == 'g' && out_file_name.at(length-2) == 'r' && out_file_name.at(length-1) == 'o') //gro
        {
            out_f = 0;

            //set temporary output trajectory name
            out_file_name_tmp = out_file_name;
            out_file_name_tmp.pop_back();
            out_file_name_tmp.pop_back();
            out_file_name_tmp.pop_back();
            out_file_name_tmp.pop_back();
            out_file_name_tmp = out_file_name_tmp + "_" + to_string(world_rank) + ".gro";
        }
        else if(out_file_name.at(length-4) == '.' && out_file_name.at(length-3) == 'p' && out_file_name.at(length-2) == 'd' && out_file_name.at(length-1) == 'b') //pdb
        {
            out_f = 1;

            //set temporary output trajectory name
            out_file_name_tmp = out_file_name;
            out_file_name_tmp.pop_back();
            out_file_name_tmp.pop_back();
            out_file_name_tmp.pop_back();
            out_file_name_tmp.pop_back();
            out_file_name_tmp = out_file_name_tmp + "_" + to_string(world_rank) + ".pdb";
        }
        else if(out_file_name.at(length-4) == '.' && out_file_name.at(length-3) == 'x' && out_file_name.at(length-2) == 't' && out_file_name.at(length-1) == 'c') //xtc
        {
            out_f = 2;

            //set temporary output trajectory name
            out_file_name_tmp = out_file_name;
            out_file_name_tmp.pop_back();
            out_file_name_tmp.pop_back();
            out_file_name_tmp.pop_back();
            out_file_name_tmp.pop_back();
            out_file_name_tmp = out_file_name_tmp + "_" + to_string(world_rank) + ".xtc";
        }
        else if(out_file_name.at(length-4) == '.' && out_file_name.at(length-3) == 't' && out_file_name.at(length-2) == 'r' && out_file_name.at(length-1) == 'r') //trr
        {
            out_f = 3;

            //set temporary output trajectory name
            out_file_name_tmp = out_file_name;
            out_file_name_tmp.pop_back();
            out_file_name_tmp.pop_back();
            out_file_name_tmp.pop_back();
            out_file_name_tmp.pop_back();
            out_file_name_tmp = out_file_name_tmp + "_" + to_string(world_rank) + ".trr";
        }
        else
        {
            if(world_rank == 0)
            {
                printf("Supported file types are xtc, trr, gro and pdb. Make sure your output file name ends with one of these.\n");
            }
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }

        //check if the file exists and back it up if so. 
        if(world_rank == 0)
        {
            backup(out_file_name);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function sets the lsq fitting paramters and reads in the lsq index file data                        //
//                                                                                                          //
////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
void Trajectory::set_lsq(string my_lsq_index_file_name,int my_b_lsq,int my_lsq_dim,int my_lsq_ref)
{
    //get the world size and assign ranks
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);      //get the world size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);      //get the process rank

    //set the lsq parameters
    lsq_index_file_name = my_lsq_index_file_name;
    b_lsq               = my_b_lsq;                                    
    lsq_dim             = my_lsq_dim;                                  
    lsq_ref             = my_lsq_ref;                                  

    //allocate memory to hold first shift record
    //lsq_shift = (rvec *)calloc(num_atoms , sizeof(*lsq_shift));
    lsq_shift = (rvec *)calloc(        1 , sizeof(*lsq_shift));

    //initialize the lsq index
    lsq_index.resize(0,0);

    //read in the lsq index
    if(b_lsq == 1)
    {
        //check the index file extenstion 
        int length = lsq_index_file_name.length();
        if(lsq_index_file_name.at(length-4) == '.' && lsq_index_file_name.at(length-3) == 'n' && lsq_index_file_name.at(length-2) == 'd' && lsq_index_file_name.at(length-1) == 'x') //ndx
        {
        }
        else //wrong extension
        {
            if(world_rank == 0)
            {
                printf("Supported file types for least squares fitting are .ndx. Make sure your index file ends with this.\n");
            }
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }

        //create data for working with index file
        char my_string[20];    //used for reading in strings from the lsq index file
        FILE *index_file;      //file used for reading in data from the lsq index file

        //open index file
        index_file = fopen64(lsq_index_file_name.c_str(), "r");
        if(index_file == NULL)
        {
            if(world_rank == 0)
            {
                printf("failure opening least squares fitting index file (%s). Make sure the file exists. \n",lsq_index_file_name.c_str());
            }
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }
        else 
        {
            //read data from file
            while(fscanf(index_file, "%s,", my_string) == 1)
            {
                lsq_index.push_back(atoi(my_string));
            }

            //close index file
            fclose(index_file);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function sets the trajectory resolution parameters, i.e., the stride, and begin/end frame           //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::set_res(int my_stride,int my_start_frame,int my_end_frame)
{
    //get the world size and assign ranks
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);      //get the world size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);      //get the process rank

    //set the resolution parameters
    stride      = my_stride;
    start_frame = my_start_frame;
    end_frame   = my_end_frame;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function readies the trajectory for reading/writing.                                                //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Trajectory::build()
{
    //Major objectives of this function include
    //analyze trajectory file to get position of each frame
    //analyze the reference file
    //analyze first trajectory frame
    //allocate memory for select data structures
    //find the first and last atom of each residue
    //distribute the workload over mpi processes
    //open trajectory files for reading/writing
    
    int i = 0;  //standard variable used in loops
    int j = 0;  //standard variable used in loops

    //take the time when build started
    clock_t t_i = clock();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                          //
    // Shape vectors                                                                                            //
    //                                                                                                          //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //get the world size and assign ranks
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);      //get the world size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);      //get the process rank

    //size vectors needed to analyze traj
    world_frames.resize(world_size,0);
    world_frame_i.resize(world_size,0);
    world_frame_f.resize(world_size,0);
    pos.resize(0,0);

    //clear boxes
    clear_mat(box);
    clear_mat(ibox);

    //clear title
    init_carray(title,200);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                          //
    // Open input trajectory and ref files for types gro and pdb (xtc and trr files are opened using GROMACS    //
    // functions)                                                                                               //
    //                                                                                                          //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //input trajectory
    if(in_f == 0 || in_f == 1) //gro or pdb
    {
        in_file = fopen64(traj_file_name.c_str(), "r");
        if(in_file == NULL)
        {
            if(world_rank == 0)
            {
                printf("failure opening %s. Make sure the file exists. \n",traj_file_name.c_str());
            }
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }
    }
    //reference file
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
    MPI_Barrier(MPI_COMM_WORLD);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                          //
    // Analyze trajectory to get the number of atoms and frames and the position of each. Check for an info     //
    // file. if present, read it instead. otherwise, analyze the trajectory and write an info file              //    
    //                                                                                                          //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //check if trajectory info file exists
    string info_file_name = traj_file_name + ".info";
    FILE *info_file = fopen(info_file_name.c_str(),"r");
    if(info_file == NULL)
    {
        if(world_rank == 0)
        {
            printf("Could not find an info file for %s. Will analyze the trajectory.\n\n",traj_file_name.c_str());
        }

        if(in_f == 0) //gro
        {
            analyze_gro_file(&in_file,&num_atoms,&frames,world_rank,
                             traj_file_name,pos,&filesize);
        }
        else if(in_f == 1) //pdb
        {
            analyze_pdb_file(&in_file,&num_atoms,&frames,world_rank,
                             traj_file_name,pos,&filesize);
        }
        else if(in_f == 2) //xtc
        {
            analyze_xtc_file(&num_atoms,&frames,world_rank,cname(traj_file_name),
                             pos,&step,&time,&prec,&status_xtc,&filesize);
        }
        else if(in_f == 3) //trr
        {
            analyze_trr_file(&num_atoms,&frames,world_rank,cname(traj_file_name),
                             pos,&step,&time,&status_trr,&filesize);
        }

        //wite info file
        if(world_rank == 0)
        {
            printf("Writing info file %s \n\n",info_file_name.c_str());

            FILE *info_file = fopen(info_file_name.c_str(),"w");
            fprintf(info_file," %ld \n",filesize);
            fprintf(info_file," %d \n",num_atoms);
            fprintf(info_file," %d \n",frames);
            for(i=0; i<frames; i++)
            {
                fprintf(info_file," %ld \n",pos[i]);
            }
            fclose(info_file);
        }
    }
    else //read the .info file 
    {
        char my_string[20];   //used to read in each item in the info file
        int line_count = 0;   //count the lines as they are read

        if(world_rank == 0)
        {
            printf("Analyzing %s. \n",info_file_name.c_str());
        }

        while(fscanf(info_file, "%s,", my_string) == 1)
        {
            if(line_count == 0) //file size
            {
                filesize = atol(my_string);
            }
            else if(line_count == 1) //num_atoms
            {
                num_atoms = atoi(my_string);
            }
            else if(line_count == 2) //frames
            {
                frames = atoi(my_string);
            }
            else //position of each frame
            {
                pos.push_back(atol(my_string));
            }
            line_count++;
        }
        fclose(info_file);

        if(world_rank == 0)
        {
            printf("Finished analyzing %s. \n",info_file_name.c_str());
            printf("Trajectory frames: %-20d \n\n",frames);
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                          //
    // Resize pos for world_rank not equal zero since only rank zero analyzed the trajectory                    //
    //                                                                                                          //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(world_rank != 0)
    {
        pos.resize(frames,0);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                          //
    // Broadcast pos data. We copy data to an array and then broadcast the array. Could proably broadcast the   //
    // vector directly but not certain how.                                                                     //
    //                                                                                                          //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int64_t pos_ary[frames];      //array to hold frame position data

    //copy trajectory frame position data to array
    if(world_rank == 0)
    {
        for(i=0; i<frames; i++)
        {
            pos_ary[i] = pos[i];
        }
    }

    //broadcast trajectory frame position data
    MPI_Bcast(pos_ary, frames, MPI_LONG, 0, MPI_COMM_WORLD);

    //copy trajectory frame position back to the vector
    for(i=0; i<frames; i++)
    {
        pos[i] = pos_ary[i];
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                          //
    // Initialize end frame number                                                                              //
    //                                                                                                          //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //set the default end frame value here since the number of frames is now known
    if(end_frame == -1) //was not specified by user
    {
        end_frame = frames - 1;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                          //
    // Allocate memory to hold trajectory data                                                                  //
    //                                                                                                          //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //atom and residue names and numbers
    atom_nr.resize(num_atoms,0);
    res_nr.resize(num_atoms,0);
    res_name.resize(num_atoms);
    atom_name.resize(num_atoms);

    //infomration for pdb format
    beta.resize(num_atoms,0);
    weight.resize(num_atoms,0);
    element.resize(num_atoms);
    chain_id.resize(num_atoms,'A');

    //allocate membory for r/v/f
    r       = (rvec *)calloc(num_atoms , sizeof(*r));
    v       = (rvec *)calloc(num_atoms , sizeof(*v));
    f       = (rvec *)calloc(num_atoms , sizeof(*f));
    r_ref   = (rvec *)calloc(num_atoms , sizeof(*r_ref));
    r_store = (rvec *)calloc(num_atoms , sizeof(*r_store));
    v_ref   = (rvec *)calloc(num_atoms , sizeof(*v_ref));
    r0      = (rvec *)calloc(num_atoms , sizeof(*r0));

    //allocate memory to hold atomic masses
    mass     = (real *)calloc(num_atoms , sizeof(real));
    mass_lsq = (real *)calloc(num_atoms , sizeof(real));

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                          //
    // Get initial coords                                                                                       //
    //                                                                                                          //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(in_f == 0) //gro
    {
        int dummy_num_atoms = num_atoms; //no need to read num_atoms from frame 1 here. Already done in analyze_traj()

        //read frist frame from trajectory
        read_gro_frame_by_char(&in_file,ibox,&dummy_num_atoms,atom_nr,res_nr,res_name,atom_name,r0,v,
                               title,world_rank,&time,&step,frames,&bV);

        //rewind file for future reading
        rewind(in_file);
    }
    else if(in_f == 1) //pdb
    {
        //read frist frame from trajectory
        read_pdb_frame_by_char(&in_file,ibox,atom_nr,res_nr,res_name,atom_name,r0,
                               title,world_rank,&time,&step,frames,beta,weight,element,
                               chain_id,&bBox);

        //rewind file for future reading
        rewind(in_file);
    }
    else if(in_f == 2) //xtc
    {
        int i = 0;          //standard variable used in loops
        rvec *x;            //holds the atomic coordinates for each frame
        XDRFILE *xd;        //used for opening xtc files

        //allocate memory to hold coordinates
        x = (rvec *)calloc(num_atoms , sizeof(*x));

        //open xtc file for reading
        xd = xdrfile_open(cname(traj_file_name), "r");
        if(NULL == xd)
        {
            printf("Could not open trajectory file (%s). Check that the file exists. \n",cname(traj_file_name));
        }
        else
        {
            status_xtc = read_xtc(xd, num_atoms, &step, &time, ibox, x, &prec); //read first frame
            xdrfile_close(xd);
        }

        //copy x to r0
        for(i=0; i<num_atoms; i++) //loop over system atoms
        {
            for(j=0; j<3; j++) //loop over x,y,z coords
            {
                r0[i][j] = x[i][j];
            }
        }
 
        //free memory
        free(x);
    }
    else if(in_f == 3) //trr
    {
        int i = 0;          //standard variable used in loops
        rvec *x;            //holds the atomic coordinates for each frame
        rvec *v;            //holds the velocities for each frame
        rvec *f;            //holds the forces for each frame
        XDRFILE *xd;        //used for opening trr files
        float lambda;       //needed for working with trr files

        //allocate memory to hold coordinates
        x = (rvec *)calloc(num_atoms , sizeof(*x));
        v = (rvec *)calloc(num_atoms , sizeof(*v));
        f = (rvec *)calloc(num_atoms , sizeof(*f));

        //open trr file
        xd = xdrfile_open(cname(traj_file_name), "r");
        if(NULL == xd)
        {
            printf("Could not open trajectory file (%s). Check that the file exists. \n",cname(traj_file_name));
        }
        else
        {
            int bBox_junk,bv_junk,bf_junk;    //tells if the trajectory contains box, vel, and force data
            status_trr = read_trr(xd, num_atoms, &step, &time, &lambda, ibox, x, v, f,&bBox_junk,&bv_junk,&bf_junk); //read first frame
            xdrfile_close(xd);
        }

        //copy coords to r0
        for(i=0; i<num_atoms; i++)  //loop over system atoms
        {
            for(j=0; j<3; j++) //loop over x,y,z coords
            {
                r0[i][j] = x[i][j];
            }
        }

        //free memory
        free(x);
        free(v);
        free(f);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                          //
    // Analyze reference                                                                                        //
    //                                                                                                          //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    gmx_bool bV_ref = 0;     //tells whether a box is present or not
    int ref_frames  = 0;     //we count the frames (should be 1) but dont really need them

    if(world_rank == 0)
    {
        printf("Analyzing reference file (%s). \n",ref_file_name.c_str());
    }

    //now read the reference file
    if(ref_f == 0) //gro
    {
        //analyze the reference gro file to get the number of atoms and frames as well as the box dimensions
        analyze_gro_file_ref(&ref_file,&num_atoms_ref,&ref_frames,world_rank,box_ref,&box_dimension);

        //check that ref and traj files are compatible
        if(num_atoms_ref != num_atoms)
        {
            if(world_rank == 0)
            {
                printf("Reference and trajectory files are not compatible. Trajectory has %d atoms while the reference file contains %d atoms. \n",num_atoms,num_atoms_ref);
            }
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }

        //read a single frame of the reference gro file to get the atom names and numbers etc. 
        read_gro_frame_by_char(&ref_file,box_ref,&num_atoms_ref,atom_nr,res_nr,res_name,atom_name,r_ref,v_ref,
                               title,world_rank,&time,&step,ref_frames,&bV_ref);

        //adjust the atom and residue number to be continuous and start at 1 
        get_cont_indices(num_atoms,atom_nr,res_nr);

        //set masses to 1 (masses are unknown for ref gro file)
        for(i=0; i<num_atoms; i++)
        {
            mass[i] = 1.0;
        }
    }
    else if(ref_f == 1) //pdb
    {
        //analyze the reference pdb file to get the number of atoms and frames
        analyze_pdb_file_ref(&ref_file,&num_atoms_ref,&ref_frames,world_rank,ref_file_name);

        //check that ref and traj files are compatible
        if(num_atoms_ref != num_atoms)
        {
            if(world_rank == 0)
            {
                printf("Reference and trajectory files are not compatible. Trajectory has %d atoms while the reference file contains %d atoms. \n",num_atoms,num_atoms_ref);
            }
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }

        //read a single frame of the reference pdb file to get the atom names and numbers etc. The box dimensions are also acquired here.
        read_pdb_frame_by_char(&ref_file,box_ref,atom_nr,res_nr,res_name,atom_name,r_ref,title,world_rank,&time,&step,ref_frames,beta,weight,element,chain_id,&bBox);

        //adjust the atom and residue number to be continuous and start at 1 
        get_cont_indices(num_atoms,atom_nr,res_nr);

        //check if a box was present
        if(bBox == 0 && world_rank == 0)
        {
            printf("Could not find a box in the reference file (%s). \n",ref_file_name.c_str());
        }

        //set masses to the beta factor for ref pdb files
        for(i=0; i<num_atoms; i++)
        {
            mass[i] = beta[i];
        }
    }

    if(world_rank == 0)
    {
        printf("Finished analyzing %s. \n\n",ref_file_name.c_str());
    }

    //check that the data was read correctly
    if(world_rank == 0)
    {
        //print_gro_frame(ibox,num_atoms,atom_nr,res_nr,res_name,atom_name,r,title,world_rank,box_dimension);
        //print_gro_frame(ibox,num_atoms,atom_nr,res_nr,res_name,atom_name,r_ref,title,world_rank,box_dimension);
        //printf("box_ref_x %f box_ref_y %f box_ref_z %f \n",box_ref[XX][XX],box_ref[YY][YY],box_ref[ZZ][ZZ]);
        //printf("world_rank %d num_atoms_ref %d \n",world_rank,num_atoms_ref);
        //printf("atom_nr[%d] %d atom_nr[%d] %d res_nr[%d] %d res_nr[%d] %d \n",0,atom_nr[0],num_atoms-1,atom_nr[num_atoms-1],0,res_nr[0],num_atoms-1,res_nr[num_atoms-1]);
    }

    //close ref file
    if(ref_f == 0 || ref_f == 1) //gro or pdb
    {
        fclose(ref_file);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                          //
    // Get residue start end. That is we get the first and last atom of each residue                            //
    //                                                                                                          //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    num_residues = res_nr[num_atoms-1];            //How many residues are there in the system

    //allocate memory for vectors
    res_start.resize(num_residues,0);
    res_end.resize(num_residues,0);

    //find the first and last atom of each residue 
    get_first_last(num_atoms,res_nr,res_start,res_end);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                          //
    // Count frames for each mpi process (my_frames) as well as the total after -b, -e, and -strid (ef_frames)  // 
    //                                                                                                          //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int count = 0;         //Used to keep count of how many frames we are actually working with after -b, -e, and -stride

    //here we determine how many frames each node is responsible for
    for(i=0; i<frames; i++) //loop over all traj frames
    {
        if(i >= start_frame && i <= end_frame) //account for -b -e
        {
            if(i%stride == 0) //account for -stride
            {
                if(count%world_size == world_rank) 
                {
                    my_frames = my_frames + 1;
                }
                count++;
            }
        }
    }

    //create an array to do communication. then copy data back to the vector
    int world_frames_ary[world_size];  //temporary array to hold my_frames from each mpi process

    //Collect my_frames and distribute to the mpi world
    MPI_Allgather(&my_frames, 1,MPI_INT,world_frames_ary, 1, MPI_INT, MPI_COMM_WORLD );

    for(i=0; i<world_size; i++)
    {
        world_frames[i] = world_frames_ary[i];
    }

    //set the effective frame count
    effective_frames = count;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                          //
    // allocate memory to hold each mpi processes frames                                                        //
    //                                                                                                          //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    pos_stride.resize(effective_frames,0);
    pos_block.resize(my_frames,0);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                          //
    // Get starting position for each frame accounting for strid -b and -e                                      //
    //                                                                                                          //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    count = 0;      //Used to keep count of how many frames we are actually working with after -b, -e, and -stride

    for(i=0; i<frames; i++) //loop over frames in traj
    {
        if(i >= start_frame && i <= end_frame) //account for -b -e
        {
            if(i%stride == 0) //account for -stride
            {
                pos_stride[count] = pos[i];
                count++;
            }
        }
    }
    count = 0;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                          //
    // Get starting position for each mpi processes frames accounting for stride, -b, -e, and block parallel    //
    //                                                                                                          //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int rank = 0;       //start with rank zero and increment over the world

    for(i=0; i<effective_frames; i++) //loop over the effective frames accounting for -b, -e, and -stride
    {
        if(rank == world_rank) //add the frame position
        {
            pos_block[count] = pos_stride[i];
        }
        count++;

        if(count == world_frames[rank]) //finished frames. move to next rank
        {
            rank++;
            count = 0;
        }
    }
    count = 0;
    rank  = 0;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                          //
    // Get the first and last frame each mpi process is responsible for                                         //
    //                                                                                                          //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int first      = 1;                      //keeps track of whether it is the mpi processes first frame
    int my_frame_i = 0;                      //first frame the mpi process is responsible for
    int my_frame_f = 0;                      //last frame the mpi process is responsible for

    for(i=0; i<effective_frames; i++) //loop over the effective frames accounting for -b, -e, and -stride
    {
        if(rank == world_rank)
        {
            if(first == 1)
            {
                my_frame_i = i;
                first = 0;
            }
            my_frame_f = i;
        }

        count++;
        if(count == world_frames[rank]) //finished frames. move to next rank
        {
            rank++;
            count = 0;
        }
    }

    //create arrays so we can collect frame_i and frame_f
    int world_frame_i_ary[world_size];                     //temporary array used to collect data
    int world_frame_f_ary[world_size];                     //temporary array used to collect data

    //collect my_frame_i and my_frame_f
    MPI_Allgather(&my_frame_i, 1,MPI_INT,world_frame_i_ary, 1, MPI_INT, MPI_COMM_WORLD );
    MPI_Allgather(&my_frame_f, 1,MPI_INT,world_frame_f_ary, 1, MPI_INT, MPI_COMM_WORLD );

    //copy world_frame_i and world_frame_f arrays to the vectors
    for(i=0; i<world_size; i++)
    {
        world_frame_i[i] = world_frame_i_ary[i];
        world_frame_f[i] = world_frame_f_ary[i];
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                          //
    // Adjust frames when block parallel is off (each process reads all frames exlcuding -b, -e, -stride)       //
    //                                                                                                          //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(block_parallel == off)
    {
        my_frames = effective_frames;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                          //
    // Open traj files for reading/writing                                                                      //
    //                                                                                                          //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //input trajectory files
    if(in_f == 2) //xtc
    {
        xd_r = xdrfile_open(cname(traj_file_name), "r");
    }
    else if(in_f == 3) //trr
    {
        trr_r = xdrfile_open(cname(traj_file_name), "r");
    }

    //temporary output trajectory files
    if(my_frames > 0)
    {
        if((out_f == 0 || out_f == 1) && b_print == 1) //gro or pdb
        {
            out_file = fopen64(out_file_name_tmp.c_str(), "w");
            if(out_file == NULL)
            {
                if(world_rank == 0)
                {
                    printf("failure opening %s. Make sure the file exists. \n",out_file_name_tmp.c_str());
                }
                MPI_Finalize();
                exit(EXIT_SUCCESS);
            }
        }
        else if(out_f == 2 && b_print == 1) //xtc
        {
            xd_w = xdrfile_open(cname(out_file_name_tmp), "w");
        }
        else if(out_f == 3 && b_print == 1) //trr
        {
            trr_w = xdrfile_open(cname(out_file_name_tmp), "w");
        }
    }

    //get the total time spent in build
    double time_build = (clock() - t_i)/CLOCKS_PER_SEC;

    return time_build;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// Returns the number of atoms in the system                                                                //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::atoms()
{
    return num_atoms;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// Returns the number of frames each mpi process is responsible for reading                                 //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<int> Trajectory::get_num_frames_world()
{
    return world_frames;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// Prints to screen that the trajectory will now be read                                                    //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::report_progress()
{
    if(world_rank == 0)
    {
        printf("Reading the file %s and preparing for analysis. \n",traj_file_name.c_str());
        printf("%10s\n","-----------------------------------------------------------------------------------------------------------------------------------");
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the geometric center of a group of atoms                                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
dv1d Trajectory::center(sv1d &target_atoms,int min,int max)
{
    int    j         = 0;     //standard variable used in loops
    int    k         = 0;     //standard variable used in loops
    int    num_atoms = 0;     //how many atoms contributing to the center
    double cx        = 0;     //used to add up the x component from each atom
    double cy        = 0;     //used to add up the y component from each atom
    double cz        = 0;     //used to add up the z component from each atom

    dv1d   r_center(3,0.0);   //the geometric center

    for(j=min; j<=max; j++) //loop over molecule atoms
    {
        for(k=0; k<target_atoms.size(); k++) //loop over target atoms
        {
            if(strcmp(atom_name[j].c_str(), target_atoms[k].c_str()) == 0) //atom type is correct
            {
                cx = cx + r[j][0];
                cy = cy + r[j][1];
                cz = cz + r[j][2];
                num_atoms++;
                break;
            }
        }
    }
    r_center[0] = cx/(double)num_atoms;
    r_center[1] = cy/(double)num_atoms;
    r_center[2] = cz/(double)num_atoms;

    return r_center;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the center of mass for a group of atoms                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
dv1d Trajectory::com(sv1d &target_atoms,int min,int max)
{
    int    j         = 0;     //standard variable used in loops
    int    k         = 0;     //standard variable used in loops
    double cx        = 0;     //used to add up the x component from each atom times the mass
    double cy        = 0;     //used to add up the y component from each atom times the mass
    double cz        = 0;     //used to add up the z component from each atom times the mass
    double M         = 0;     //total mass of all atoms

    dv1d   r_com(3,0.0);      //the center of mass

    for(j=min; j<=max; j++) //loop over molecule atoms
    {
        for(k=0; k<target_atoms.size(); k++) //loop over target atoms
        {
            if(strcmp(atom_name[j].c_str(), target_atoms[k].c_str()) == 0) //atom type is correct
            {
                cx = cx + mass[j]*r[j][0];
                cy = cy + mass[j]*r[j][1];
                cz = cz + mass[j]*r[j][2];
                M  = M + mass[j];
                break;
            }
        }
    }
    if(M==0)
    {
        printf("Sum of masses was zero when computing the center of mass. Check your mass data in the reference file. \n");
    }
    else
    {
        r_com[0] = cx/M;
        r_com[1] = cy/M;
        r_com[2] = cz/M;
    }

    return r_com;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the center of mass for a group (index) of atoms                                    //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
dv1d Trajectory::com_i(iv1d group)
{
    int    i         = 0;     //standard variable used in loops
    double cx        = 0.0;   //used to add up the x component from each atom times the mass
    double cy        = 0.0;   //used to add up the y component from each atom times the mass
    double cz        = 0.0;   //used to add up the z component from each atom times the mass
    double M         = 0.0;   //total mass of all atoms

    dv1d   r_com(3,0.0);      //the center of mass

    for(i=0; i<group.size(); i++) //loop over group of atoms
    {
        cx = cx + mass[group[i]-1]*r[group[i]-1][0];
        cy = cy + mass[group[i]-1]*r[group[i]-1][1];
        cz = cz + mass[group[i]-1]*r[group[i]-1][2];
        M  = M  + mass[group[i]-1];
    }
  
    if(M==0)
    {
        printf("Sum of masses was zero when computing the center of mass. Check your mass data in the reference file. \n");
    }
    else 
    {
        r_com[0] = cx/M;
        r_com[1] = cy/M;
        r_com[2] = cz/M;
    }

    return r_com;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the geometric center of a group (index) of atoms                                   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
dv1d Trajectory::center_i(iv1d &group)
{
    int    i         = 0;     //standard variable used in loops
    double cx        = 0;     //used to add up the x component from each atom
    double cy        = 0;     //used to add up the y component from each atom
    double cz        = 0;     //used to add up the z component from each atom

    dv1d   r_center(3,0.0); //the geometric center

    for(i=0; i<group.size(); i++) //loop over group of atoms
    {
        cx = cx + r[group[i]-1][0];
        cy = cy + r[group[i]-1][1];
        cz = cz + r[group[i]-1][2];
    }
    r_center[0] = cx/(double)group.size();
    r_center[1] = cy/(double)group.size();
    r_center[2] = cz/(double)group.size();

    return r_center;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the geometric center of a group of atoms from stored coordinates                   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
dv1d Trajectory::center_store(sv1d &target_atoms,int min,int max)
{
    int    j         = 0;     //standard variable used in loops
    int    k         = 0;     //standard variable used in loops
    int    num_atoms = 0;     //how many atoms contributing to the center
    double cx        = 0;     //used to add up the x component from each atom
    double cy        = 0;     //used to add up the y component from each atom
    double cz        = 0;     //used to add up the z component from each atom

    dv1d   r_center(3,0.0);   //the geometric center

    for(j=min; j<=max; j++) //loop over molecule atoms
    {
        for(k=0; k<target_atoms.size(); k++) //loop over target atoms
        {
            if(strcmp(atom_name[j].c_str(), target_atoms[k].c_str()) == 0) //atom type is correct
            {
                cx = cx + r_store[j][0];
                cy = cy + r_store[j][1];
                cz = cz + r_store[j][2];
                num_atoms++;
                break;
            }
        }
    }
    r_center[0] = cx/(double)num_atoms;
    r_center[1] = cy/(double)num_atoms;
    r_center[2] = cz/(double)num_atoms;

    return r_center;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the radius of gyration for a group of atoms                                        //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Trajectory::gyrate(sv1d &target_atoms,int min,int max)
{
    int    j         = 0;     //standard variable used in loops
    int    k         = 0;     //standard variable used in loops
    double M         = 0;     //total mass of all atoms
    double gy        = 0;     //the radius of gyration

    dv1d   dif(3,0.0);        //dif between atomic coordinates and the center of mass

    dv1d R = com(target_atoms,min,max);  //the center of mass  

    for(j=min; j<=max; j++) //loop over molecule atoms
    {
        for(k=0; k<target_atoms.size(); k++) //loop over target atoms
        {
            if(strcmp(atom_name[j].c_str(), target_atoms[k].c_str()) == 0) //atom type is correct
            {
                dif[0] = r[j][0] - R[0];
                dif[1] = r[j][1] - R[1];
                dif[2] = r[j][2] - R[2];

                double square = dif[0]*dif[0] + dif[1]*dif[1] + dif[2]*dif[2];

                gy = gy + mass[j]*square;
                M  = M + mass[j];
                break;
            }
        }
    }
    gy = sqrt(gy/M);

    return gy;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function takes in 4 atoms and returns the dihedral angle (radians)                                   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Trajectory::dihedral_angle(int i,int j,int k,int l)
{  
    //Notes: The dihedral angle between 4 atoms (i,j,k,l) is defined as the angle between the two planes made by 
    //       atoms i,j,k and atoms j,k,l. To compute the angle we compute a vector that is normal to each plane (m and n).
    //       Then, the angle between the 2 vectors is the same as the angle between the 2 planes. The angle between the vectors 
    //       is computed with standard trigonometry relations. 
  
    rvec difm_ij;   //the vector connecting atoms i and j
    rvec difm_kj;   //the vector connecting atoms k and j
    rvec difm_kl;   //the vector connecting atoms k and l
    rvec m;         //the vector normal to the plane made by atoms i,j,k
    rvec n;         //the vector normal to the plane made by atoms k,j,l

    //rij
    difm_ij[0] = r[i][0] - r[j][0];
    difm_ij[1] = r[i][1] - r[j][1];
    difm_ij[2] = r[i][2] - r[j][2];
    
    //rkj
    difm_kj[0] = r[k][0] - r[j][0];
    difm_kj[1] = r[k][1] - r[j][1];
    difm_kj[2] = r[k][2] - r[j][2];
    
    //rkl
    difm_kl[0] = r[k][0] - r[l][0];
    difm_kl[1] = r[k][1] - r[l][1];
    difm_kl[2] = r[k][2] - r[l][2];
    
    //now compute the vectors normal to the planes
    cprod(difm_ij,difm_kj,m);
    cprod(difm_kj,difm_kl,n);

    //compute the dihedral angle
    double angle = gmx_angle(m,n);

    //now determine the sign of the angle
    double ipr   = iprod(difm_ij, n);
    int    sign  = (ipr < 0.0) ? -1.0 : 1.0;
    
    return (sign)*angle;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function takes the current coordinates and stores them in r_store                                    //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::store_coords()
{
    int i = 0;       //standard variable used in loops 

    for(i=0; i<num_atoms; i++)
    {
        r_store[i][0] = r[i][0];
        r_store[i][1] = r[i][1];
        r_store[i][2] = r[i][2];
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function takes an atom index i (traj.atom_name[i]) and returns the first atom of the current residue //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::get_res_start(int i)
{
    return res_start[res_nr[i]-1];
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function takes an atom index i (traj.atom_name[i]) and returns the last atom of the current residue  //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::get_res_end(int i)
{
    return res_end[res_nr[i]-1];
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the last atom of the current residue so after i++ we have the next residue          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::next_residue(int i)
{
    return res_end[res_nr[i]-1];
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function checks for broken residues.                                                                 //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::check_broken_molecule(int flag_beta)
{
    int i = 0;     //standard variable used in loops
    int j = 0;     //standard variable used in loops
    int k = 0;     //standard variable used in loops
    int l = 0;     //standard variable used in loops

    for(i=0; i<atoms(); i++) //loop over system atoms
    {
        int min = get_res_start(i);  //get first atom of residue
        int max = get_res_end(i);    //get last atom of residue

        i = next_residue(i);         //prime i for next residue

        int broken = 0;              //tells if the residue is broken

        for(j=min; j<=max; j++) //loop over residue atoms
        {
            beta[j] = 0.0;

            for(k=min; k<=max; k++) //loop over residue atoms
            {
                double dx = r[j][0] - r[k][0];
                double dy = r[j][1] - r[k][1];
                double dz = r[j][2] - r[k][2];

                double distance = sqrt(dx*dx + dy*dy + dz*dz);

                if(distance > 0.5*box[XX][XX] || distance > 0.5*box[YY][YY] || distance > 0.5*box[ZZ][ZZ])
                {
                    broken = 1;

                    if(flag_beta == 1)
                    {
                        for(l=min; l<=max; l++) //loop over residue atoms
                        {
                            beta[l] = 1.0;
                        }
                    }

                    goto end_loop;
                }
            }
        }
        end_loop:;

        if(broken == 1 && world_rank == 0)
        {
            printf("Frame %d Residue %d %s is probably broken. ",current_frame,res_nr[min],res_name[min].c_str());
            for(j=min; j<=max; j++) //loop over residue atoms
            {
                printf("atom %d %s ",j,atom_name[j].c_str());
            }
            printf("\n");
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function measures the distance between two atoms                                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Trajectory::get_dist(dv1d &vec_a,dv1d &vec_b)
{
    double dx = vec_a[0] - vec_b[0];
    double dy = vec_a[1] - vec_b[1];
    double dz = vec_a[2] - vec_b[2];

    double dist = sqrt(dx*dx + dy*dy + dz*dz);

    return dist;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function takes a trajectory and generates the surrounding 8 periodic images                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trajectory::expand_traj(Trajectory traj_ref)
{
    init_carray(title,200);

    //size the trajectory structures
    atom_nr.resize(9*traj_ref.atoms(),0);
    res_nr.resize(9*traj_ref.atoms(),0);
    res_name.resize(9*traj_ref.atoms());
    atom_name.resize(9*traj_ref.atoms());

    //information for pdb format
    beta.resize(9*traj_ref.atoms(),0);
    weight.resize(9*traj_ref.atoms(),0);
    element.resize(9*traj_ref.atoms());
    chain_id.resize(9*traj_ref.atoms(),'A');

    //allocate membory for r/v/f
    r       = (rvec *)calloc(9*traj_ref.atoms() , sizeof(*r));
    v       = (rvec *)calloc(9*traj_ref.atoms() , sizeof(*v));
    f       = (rvec *)calloc(9*traj_ref.atoms() , sizeof(*f));
    r_ref   = (rvec *)calloc(9*traj_ref.atoms() , sizeof(*r_ref));
    r_store = (rvec *)calloc(9*traj_ref.atoms() , sizeof(*r_store));
    v_ref   = (rvec *)calloc(9*traj_ref.atoms() , sizeof(*v_ref));
    r0      = (rvec *)calloc(9*traj_ref.atoms() , sizeof(*r0));

    //allocate memory to hold atomic masses
    mass = (real *)calloc(9*traj_ref.atoms() , sizeof(real));
    if(b_lsq == 1)
    {
        mass_lsq = (real *)calloc(9*traj_ref.atoms() , sizeof(real));
    }

    res_start.resize(9*traj_ref.get_num_residues(),0);
    res_end.resize(9*traj_ref.get_num_residues(),0);

    //get the world size and assign ranks
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);      //get the world size
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);      //get the process rank

    world_frames.resize(world_size,0);

    //adjust variables
    frames               = traj_ref.get_frames();              
    num_atoms            = 9*traj_ref.atoms();                 
    num_atoms_ref        = 9*traj_ref.atoms();                 
    box_dimension        = traj_ref.get_box_dim();                 
    my_frames            = traj_ref.get_num_frames();                 
    effective_frames     = traj_ref.get_ef_frames();                 
    stride               = traj_ref.get_stride();                 
    start_frame          = traj_ref.get_start_frame();                 
    end_frame            = traj_ref.get_end_frame();
    prec                 = traj_ref.get_prec();

    world_frames = traj_ref.get_num_frames_world();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                          //
    // Open traj files for reading/writing                                                                      //
    //                                                                                                          //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //input trajectory
    if(in_f == 0 || in_f == 1) //gro or pdb
    {
        in_file = fopen64(traj_file_name.c_str(), "r");
        if(in_file == NULL)
        {
            if(world_rank == 0)
            {
                printf("failure opening %s. Make sure the file exists. \n",traj_file_name.c_str());
            }
            MPI_Finalize();
            exit(EXIT_SUCCESS);
        }
    }

    //input trajectory files
    if(in_f == 2) //xtc
    {
        xd_r = xdrfile_open(cname(traj_file_name), "r");
    }
    else if(in_f == 3) //trr
    {
        trr_r = xdrfile_open(cname(traj_file_name), "r");
    }

    //temporary output trajectory files
    if(my_frames > 0)
    {
        if((out_f == 0 || out_f == 1) && b_print == 1) //gro or pdb
        {
            out_file = fopen64(out_file_name_tmp.c_str(), "w");
            if(out_file == NULL)
            {
                if(world_rank == 0)
                {
                    printf("failure opening %s. Make sure the file exists. \n",out_file_name_tmp.c_str());
                }
                MPI_Finalize();
                exit(EXIT_SUCCESS);
            }
        }
        else if(out_f == 2 && b_print == 1) //xtc
        {
            xd_w = xdrfile_open(cname(out_file_name_tmp), "w");
        }
        else if(out_f == 3 && b_print == 1) //trr
        {
            trr_w = xdrfile_open(cname(out_file_name_tmp), "w");
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the centers of the lipids in the target leaflet                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
dv2d Trajectory::get_centers_target_lf()
{
    int i = 0;                                //standard variable used in loops
    int j = 0;                                //standard variable used in loops

    dv2d lipid_centers(0,dv1d(3,0.0));        //hold the center of each lipid molecule in the target leaflet

    for(i=0; i<target_leaflet.size(); i++) //loop over target leaflet
    {
        //get the first and last atom of the current lipid
        int min = t_lip_start(i);
        int max = t_lip_end(i);

        //jump to the next lipid
        i = next_target_lipid(i);

        iv1d target_atoms(0);

        for(j=min; j<=max; j++) //loop over current lipid
        {
            if(atom_name[j].at(0) != 'H') //atom is not a hydrogen
            {
                target_atoms.push_back(atom_nr[j]);
            }
        }

        dv1d this_center = center_i(target_atoms);
        lipid_centers.push_back(this_center);
    }

    return lipid_centers; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the centers of the lipids in the opposing leaflet                                   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
dv2d Trajectory::get_centers_opposing_lf()
{
    int i = 0;                                //standard variable used in loops
    int j = 0;                                //standard variable used in loops

    dv2d opposing_centers(0,dv1d(3,0.0));     //hold the center of each lipid molecule in the opposing leaflet

    for(i=0; i<opposing_leaflet.size(); i++) //loop over opposing leaflet
    {
        //get the first and last atom of the current lipid
        int min = o_lip_start(i);
        int max = o_lip_end(i);

        //jump to the next lipid
        i = next_opposing_lipid(i);

        iv1d target_atoms(0);

        for(j=min; j<=max; j++) //loop over current lipid
        {
            if(atom_name[j].at(0) != 'H') //atom is not a hydrogen
            {
                target_atoms.push_back(atom_nr[j]);
            }
        }

        dv1d this_center = center_i(target_atoms);
        opposing_centers.push_back(this_center);
    }

    return opposing_centers;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the centers of the lipids in the full leaflet                                       //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
dv2d Trajectory::get_centers_full_mem()
{
    int i = 0;                                //standard variable used in loops
    int j = 0;                                //standard variable used in loops

    dv2d full_mem_centers(0,dv1d(3,0.0));     //hold the center of each lipid molecule in the full mem

    for(i=0; i<full_membrane.size(); i++) //loop over full membrane 
    {
        //get the first and last atom of the current lipid
        int min = fm_lip_start(i);
        int max = fm_lip_end(i);

        //jump to the next lipid
        i = next_full_mem_lipid(i);

        iv1d target_atoms(0);

        for(j=min; j<=max; j++) //loop over current lipid
        {
            if(atom_name[j].at(0) != 'H') //atom is not a hydrogen
            {
                target_atoms.push_back(atom_nr[j]);
            }
        }

        dv1d this_center = center_i(target_atoms);
        full_mem_centers.push_back(this_center);
    }

    return full_mem_centers;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the centers of the residues in the protein                                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
dv2d Trajectory::get_centers_prot()
{
    int i = 0;                                //standard variable used in loops
    int j = 0;                                //standard variable used in loops

    dv2d protein_centers(0,dv1d(3,0.0));      //hold the center of each residue in the protein

    for(i=0; i<prot.size(); i++)  //loop over protein atoms
    {
        //get the first and last atom of the current residue
        int min = p_res_start(i);
        int max = p_res_end(i);

        //jump to the next residue
        i = next_prot_res(i);

        iv1d target_atoms(0);

        for(j=min; j<=max; j++) //loop over current residue
        {
            if(atom_name[j].at(0) != 'H') //atom is not a hydrogen
            {
                target_atoms.push_back(atom_nr[j]);
            }
        }

        dv1d this_center = center_i(target_atoms);
        protein_centers.push_back(this_center);
    }

    return protein_centers; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the centers of the residues in the solvent                                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
dv2d Trajectory::get_centers_sol()
{
    int i = 0;                                //standard variable used in loops
    int j = 0;                                //standard variable used in loops

    dv2d sol_centers(0,dv1d(3,0.0));          //hold the center of each residue in the solvent

        for(i=0; i<sol.size(); i++)  //loop over solvent atoms
        {
            //get the first and last atom of the current water
            int min = sol_start(i);
            int max = sol_end(i);

            //jump to the next water
            i = next_water(i);

            iv1d target_atoms(0);

            for(j=min; j<=max; j++) //loop over current water
            {
                if(atom_name[j].at(0) != 'H') //atom is not a hydrogen
                {
                    target_atoms.push_back(atom_nr[j]);
                }
            }

            dv1d this_center = center_i(target_atoms);
            sol_centers.push_back(this_center);
        }

    return sol_centers;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the centers of the residues in the system                                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
dv2d Trajectory::get_centers_system()
{
    int i = 0;                                //standard variable used in loops
    int j = 0;                                //standard variable used in loops

    dv2d system_centers(0,dv1d(3,0.0));       //hold the center of each residue in the system

    for(i=0; i<atoms(); i++)  //loop over system atoms
    {
        //get the first and last atom of the current residue
        int min = get_res_start(i);
        int max = get_res_end(i);

        //jump to the next residue
        i = next_residue(i);

        iv1d target_atoms(0);

        for(j=min; j<=max; j++) //loop over current residue
        {
            if(atom_name[j].at(0) != 'H') //atom is not a hydrogen
            {
                target_atoms.push_back(atom_nr[j]);
            }
        }

        dv1d this_center = center_i(target_atoms);
        system_centers.push_back(this_center);
    }

    return system_centers;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the first frame to be read for the core in the global index                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::get_global_frame_i()
{
    return world_frame_i[world_rank]; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the last frame to be read for the core in the global index                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Trajectory::get_global_frame_f()
{
    return world_frame_f[world_rank];
}

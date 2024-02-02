
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function counts the number of lines in a text file.                                                  //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int count_lines(string text_file_name)
{
    int num_lines = 0;         //number of lines to be returned by function
    int c = 0;                 //used for reading in a character with fgetc

    //open the text file
    FILE *text_file = fopen(text_file_name.c_str(), "r");
    if(text_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",text_file_name.c_str());
    }
    else
    {
        //count lines in text file
        do
        {
            c = fgetc(text_file);
            if (c == '\n')
            {
                num_lines = num_lines + 1;
            }
        }
        while (c != EOF);

        //close the text file
        fclose(text_file);
    }

    return num_lines;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function will read in a line from the buffer. Works similar to fgets but for a buffer instead of a   //
// file.                                                                                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void read_line_buf(int world_rank,int nchars,char buf[],char my_line[],int *buf_offset)
{
    int i = 0;
    int count = 0;
    for(i=0; i<200; i++)
    {
        my_line[i] = ' ';
    }

    for(i=(*buf_offset); i<nchars; i++)
    {
        if(buf[i] != '\n')
        {
            my_line[count] = buf[i];
            count++;
        }
        else
        {
            my_line[count] = '\n';
            my_line[count+1] = '\0';
            *buf_offset = i + 1;
            break;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function will read in the next string from the buffer. Works similar to fscanf but for a buffer      //
// instead of a file.                                                                                        //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int next_string(int world_rank,int nchars,char buf[],char my_string[],int string_size,int *buf_offset)
{
    int result = 0;

    int i = 0;
    int count = 0;
    int start = 0;
    int end = 0;

    for(i=0; i<string_size; i++)
    {
        my_string[i] = ' ';
    }

    for(i=(*buf_offset); i<nchars; i++)
    {
        if(buf[i] != ' ' && buf[i] != '\n')
        {
            start = 1;
            if(end == 0)
            {
                my_string[count] = buf[i];
                count++;
            }
            else if(end == 1)
            {
                //next string reached
                *buf_offset = i;
                break;
            }

            result = 1;
        }
        else if(buf[i] == ' ' && start == 1)
        {
            my_string[count] = '\0';
            end = 1;
        }
        else if(buf[i] == ' ' && start == 0)
        {
            //Do nothing. Read in next character.
        }
        else if(buf[i] == '\n')
        {
            //end of line reached
            *buf_offset = i;
            my_string[count] = '\0';
            break;
        }
    }
    return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function takes an array of chars and removes blank spaces and add \0 to the end                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void clean_string(char my_string[],int size)
{
    int i = 0;
    int count = 0;
    char my_string_copy[size];

    //make a copy of the string and clear the original
    for(i=0; i<size; i++)
    {
        my_string_copy[i] = my_string[i];
        my_string[i] = ' ';
    }

    //now remove blank spaces and add \0 to end
    for(i=0; i<size; i++)
    {
        if(my_string_copy[i] != ' ')
        {
            my_string[count] = my_string_copy[i];
            count++;
        }
    }
    my_string[count] = '\0';
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This functionc reads to the next line                                                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void next_line(FILE *in_file)
{
    int i = 0;
    int c = 0;

    do
    {
      c = fgetc(in_file);
      if (c == '\n')
      {
          i++;
      }
    }
    while (i < 1);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function counts the number of lines and items per line in a file                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void characterize_file(string my_file_name,int *number_of_lines,int *items_per_line,int header_lines)
{
    FILE *my_file;                //File for reading in data
    int string_size = 2000;       //How big of a string can we read?
    char my_string[string_size];  //String to hold read in data entries
    int capacity  = 0;            //Number of items in file
    int c         = 0;            //Used to read in chars using fgetc
    int num_lines = 0;            //Number of lines in the file
    int i         = 0;            //Standard variable used in loops

    my_file = fopen(my_file_name.c_str(), "r");
    if(my_file == NULL) //file does not exist
    {
        printf("failure opening %s. Make sure the file exists. \n",my_file_name.c_str());
    }
    else
    {
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Count the number of lines in the data file                                                                //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        do
        {
          c = fgetc(my_file);
          if (c == '\n')
          {
              num_lines++;
          }
        }
        while (c != EOF);
        rewind(my_file);

        num_lines = num_lines - header_lines;

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Count the items per line                                                                                  //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //skip first header_lines lines in file (header)
        for(i=0; i<header_lines; )
        {
            c = fgetc(my_file);
            if (c == '\n')
            {
                i++;
            }
        }

        while(fscanf(my_file, "%s,", my_string) != EOF/*== 1*/)
        {
            capacity = capacity + 1;
        }
        rewind(my_file);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                                           //
        // Set the number_of_lines and items_per_line                                                                //
        //                                                                                                           //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        *items_per_line = capacity/num_lines;
        *number_of_lines = num_lines;
        fclose(my_file);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads to the next line                                                                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void read_line(FILE *in_file,string &line)
{
    int i = 0;
    int j = 0;
    int c = 0;
    int char_line_size = 2000;
    char char_line[char_line_size];

    do
    {
      c = fgetc(in_file);
      char_line[j] = c;

      if (c == '\n')
      {
          char_line[j+1] = '\0';
          i++;
      }
      j++;
    }
    while (i < 1);

    line = strdup(char_line);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function checks if a file exisits and backs it up if so. (don't over write files)                    //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void backup(string in_file_name)
{
    int i = 0;
    int j = 1;
    string bkp_file_name;

    //try opening the file to see if it exists
    FILE *in_file = fopen(in_file_name.c_str(), "r");
    if(in_file == NULL)
    {
    }
    else //file exists. back it up!
    {
        //create a file name for the backup
        for(i=0; i<1; )
        {
            bkp_file_name = in_file_name + ".bkp." + to_string(j);  
            FILE *test_file = fopen(bkp_file_name.c_str(), "r");
            if(test_file == NULL)
            {
                i=1;
            }
            j++;
        }

        //backup the file
        printf("Backing up file %s to %s \n",in_file_name.c_str(),bkp_file_name.c_str());
        ifstream source(in_file_name.c_str(), ios::binary);
        ofstream dest(bkp_file_name.c_str(), ios::binary);

        istreambuf_iterator<char> begin_source(source);
        istreambuf_iterator<char> end_source;
        ostreambuf_iterator<char> begin_dest(dest); 
        copy(begin_source, end_source, begin_dest);

        source.close();
        dest.close();
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function parses a c++ string (a line) to give the individual strings in the line                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
sv1d parse_string(string this_string)
{
    int i      = 0;           //standard variable used in loops
    int primed = 0;           //tells if a non space character has been encountered yet
    sv1d parsed_items;        //holds the list of string in the current line
    string current_string;    //holds characters for the current string

    for(i=0; i<this_string.length(); i++) //loop over characters in string
    {
        if( (this_string[i] == ' ' || this_string[i] == '\t' || i==this_string.length()-1) && primed == 1) //store string
        {
            parsed_items.push_back(current_string);
            string new_string;
            current_string = new_string;
            primed = 0;
        }
        else if( (this_string[i] == ' ' || this_string[i] == '\t' ) && primed == 0) //waiting to begin new string
        {
           //do nothing! 
        }
        else //character encountered 
        {
            current_string.push_back(this_string[i]);
            primed = 1;
        }
    }

    return parsed_items;
}


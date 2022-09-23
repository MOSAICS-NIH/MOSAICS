
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This is a class for working with the text files.                                                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Data_file
{
    private:
        string file_name;             //Name of the input file
        int headers         = 0;      //Number of header lines
        int number_of_lines = 0;      //Number of lines in input files
        int items_per_line  = 0;      //Number of items per line in the input file

    public:
        sv1d header{};                  //vector holding the header lines (string)
        sv2d data_s{};                  //vector holding the data (string)
        dv2d data_d{};                  //vector holding the data (double)
        iv2d data_i{};                  //vector holding the data (int)

    public:
        void get_data(string in_file_name);                                                     //read in the data
        int  size_x();                                                                          //return the size in x
        int  size_y();                                                                          //return the size in y
        void set_headers(int num_headers);                                                      //set the number of header lines
        void initialize(int x_dim, int y_dim);                                                  //allocate memory for vectors
        void push_back_row(sv1d &this_row);                                                     //adds a row to the vectors
        void write_data(string out_file_name,int b_header);                                     //write the data to output file
        void show_data();                                                                       //print data to screen
        void double_to_string();                                                                //comvert double to string
        void string_to_all();                                                                   //convert strings to int and double
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function sets the format of the grid                                                                 //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Data_file::get_data(string in_file_name)
{
    int i = 0;
    int j = 0;
    int k = 0;
    int string_size = 20000;
    file_name = in_file_name;
    char my_string[string_size];

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Characterize the data file                                                                                //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    characterize_file(file_name,&number_of_lines,&items_per_line,headers);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read in the data                                                                                          //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    data_s.resize(number_of_lines,sv1d(items_per_line));
    data_d.resize(number_of_lines,dv1d(items_per_line,0.0));
    data_i.resize(number_of_lines,iv1d(items_per_line,0));

    FILE *in_file = fopen(file_name.c_str(), "r");
    if(in_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",file_name.c_str());
    }
    else
    {
        for(i=0; i<headers; i++) //loop over header lines
        {
            string this_line; 
            read_line(in_file,this_line);
            header.push_back(this_line);
        }

        for(i=0; i<number_of_lines; i++) //loop over lines
        {
            for(j=0; j<items_per_line; j++) //loop over items per line
            {
                //for(k=0; k<string_size; k++)
                //{
                //    my_string[k] = ' ';
                //}
                fscanf(in_file, "%s,", my_string);
           
                data_s[i][j] = strdup(my_string); 
                data_d[i][j] = atof(my_string);
                data_i[i][j] = atoi(my_string);
            }
        }

        fclose(in_file);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function allocates memory for the vectors                                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Data_file::initialize(int x_dim, int y_dim)
{
    data_s.resize(y_dim,sv1d(x_dim));
    data_d.resize(y_dim,dv1d(x_dim,0.0));
    data_i.resize(y_dim,iv1d(x_dim,0));
    header.resize(headers);

    items_per_line=x_dim;
    number_of_lines=y_dim; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function adds a row to the vectors                                                                   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Data_file::push_back_row(sv1d &this_row)
{
    //check that a proper row is being added
    if(this_row.size() != items_per_line)
    {
        printf("Attempting to add a row with a different number of collumns (%d vs %d). This will cause problems. \n",items_per_line,this_row.size());
    }

    //grow the string vector
    data_s.push_back(this_row);

    //make dummy vectors so that we grow the integer and double vectors too
    dv1d dummy_d(this_row.size(),0.0);
    iv1d dummy_i(this_row.size(),0);
    
    //grow int and double vectors
    data_d.push_back(dummy_d);
    data_i.push_back(dummy_i);

    //increase the number of lines
    number_of_lines = number_of_lines + 1;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Write the data to output file                                                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Data_file::write_data(string out_file_name,int b_header)
{
    int i = 0;
    int j = 0;

    FILE *out_file = fopen(out_file_name.c_str(), "w");
    if(out_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",out_file_name.c_str());
    }   
    else
    {
        if(b_header == 1) //print header lines
        {
            for(i=0; i<headers; i++)
            {
                fprintf(out_file,"%s",header[i].c_str());
            }
        }
        
        for(i=0; i<number_of_lines; i++) //loop over y
        {
            for(j=0; j<items_per_line; j++) //loop over x
            {
                fprintf(out_file,"%15s  ",data_s[i][j].c_str());
            }
            fprintf(out_file,"\n");
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the size in x                                                                       //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Data_file::size_x()
{
    return items_per_line;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the size in y                                                                       //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Data_file::size_y()
{
    return number_of_lines;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function sets the number of header lines                                                             //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Data_file::set_headers(int num_headers)
{
    headers = num_headers; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints a data file for testing                                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Data_file::show_data()
{
    int i = 0;
    int j = 0;

    for(i=0; i<data_s.size(); i++)
    {
        for(j=0; j<data_s[i].size(); j++)
        {
            printf(" %10s ",data_s[i][j].c_str());
        }
        printf("\n");
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function converts the doubles to string                                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Data_file::double_to_string()
{
    int i = 0;
    int j = 0;

    for(i=0; i<number_of_lines; i++) //loop over y
    {
        for(j=0; j<items_per_line; j++) //loop over x
        {
            data_s[i][j] = to_string(data_d[i][j]);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function converts the doubles to string                                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Data_file::string_to_all()                                                                           
{
    int i = 0;                                                                                               
    int j = 0;

    for(i=0; i<number_of_lines; i++) //loop over y
    {
        for(j=0; j<items_per_line; j++) //loop over x
        {
            data_i[i][j] = atoi(data_s[i][j].c_str());
            data_d[i][j] = atof(data_s[i][j].c_str());
        }
    }   
}


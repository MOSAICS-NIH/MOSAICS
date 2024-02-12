
#ifdef __APPLE__
#  define off64_t off_t
#  define fopen64 fopen
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This is a class for reading index files                                                                  //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Index
{
    private:
        FILE   *index_file;                          //file for reading index       
        string index_file_name;                      //name of the index file
        char   my_string[200];                        //c-string for reading in data
        int    i = 0;                                //standard variable used in loops

    public:
        vector <string>    index_s{};                //Index of string
        vector <int>       index_i{};                //Index of int
        vector <double>    index_d{};                //Index of double
 
    public:
        void           get_index(string my_index_file_name);                   //read in the index
        void           display_index();                                        //display contents of an index file
        vector<string> get_column_s(int num_collumns,int collumn);             //extract a collumn from the index (string)
        vector<int>    get_column_i(int num_collumns,int collumn);             //extract a collumn from the index (int)
        vector<double> get_column_d(int num_collumns,int collumn);             //extract a collumn from the index (double)
        vector<string> get_row_s(int num_collumns,int row);                    //extract a row from the index (string)
        int check_next_tag(int *pos,string &tag);
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function reads in an index file                                                                     //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Index::get_index(string my_index_file_name)
{
    index_file_name = my_index_file_name;

    //open index file
    index_file = fopen64(index_file_name.c_str(), "r");
    if(index_file == NULL)
    {
        printf("failure opening index file (%s). Make sure the file exists. \n",index_file_name.c_str());
        exit(EXIT_SUCCESS);
    }
    else
    {
        //read data from file
        while(fscanf(index_file, "%s,", my_string) == 1)
        {
            string entry = strdup(my_string);

            if(entry[0] != '#')
            {
                index_s.push_back(strdup(my_string));
                index_i.push_back(atoi(my_string));
                index_d.push_back(atof(my_string));
            }
        }

        //close index file
        fclose(index_file);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function extracts a column from the index (string)                                                  //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<string> Index::get_column_s(int num_columns,int column)
{
    int i = 0;

    int vec_size = index_s.size()/num_columns;

    vector <string> extract{};
    extract.resize(vec_size);

    for(i=column; i<index_s.size(); i+=num_columns)
    {
        extract[i/num_columns] = index_s[i];
    }

    return extract;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function extracts a column from the index (int)                                                     //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<int> Index::get_column_i(int num_columns,int column)
{
    int i = 0;

    int vec_size = index_i.size()/num_columns;

    vector <int> extract{};
    extract.resize(vec_size);

    for(i=column; i<index_i.size(); i+=num_columns)
    {
        extract[i/num_columns] = index_i[i];
    }

    return extract;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function extracts a column from the index (double)                                                  //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<double> Index::get_column_d(int num_columns,int column)
{
    int i = 0;

    int vec_size = index_d.size()/num_columns;

    vector <double> extract{};
    extract.resize(vec_size);

    for(i=column; i<index_d.size(); i+=num_columns)
    {
        extract[i/num_columns] = index_d[i];
    }

    return extract;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function prints the contents of an index file                                                       //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Index::display_index()
{
    int i = 0;
 
    for(i=0; i<index_s.size(); i++)
    {
        printf(" %s \n",index_s[i].c_str());
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function extracts a row from the index (string)                                                     //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<string> Index::get_row_s(int num_columns,int row)
{
    int i     = 0;
    int count = 0;

    vector <string> extract{};
    extract.resize(num_columns);

    for(i=row*num_columns; i<(row+1)*num_columns; i++)
    {
        extract[count] = index_s[i];
        count++;
    }

    return extract;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function searches a data file for a tag and gives the position                                       //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Index::check_next_tag(int *pos,string &tag)
{
    int i      = 0;        //standard variable used in loops
    int j      = 0;        //standard variable used in loops
    int result = 0;        //returns if a tag was found
    int start  = *pos;     //first item in loop

    *pos = index_s.size(); //set position to end of file if no tag is found

    for(i=start; i<index_s.size(); i++) //loop items in index
    {
        if(index_s[i][0] == '[' && index_s[i][index_s[i].length()-1] == ']')
        {
            string this_tag;
            for(j=1; j<index_s[i].length()-1; j++) //loop over characters in string
            { 
		this_tag = this_tag + index_s[i][j];
            }
            tag = this_tag;

	    result = 1;
            *pos   = i; 
            i      = index_s.size();
        }
    }
    return result;
}


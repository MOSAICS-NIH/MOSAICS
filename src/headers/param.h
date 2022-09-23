
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This is a class for reading complex parameter files                                                      //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Param
{
    private:
        string param_file_name;                      //name of the index file
       
    public:
        Index index;                                 //Index for reading main file
        int num_lip_t;                               //The number of lipid types in file
        int num_col;                                 //How many columns in main index file (how many items per row?)
        int num_row;                                 //How many rows are in the main index file (how many items per col?)
        int file_col;                                //Which column holds the file names?
        int num_col_secondary;                       //Number of columns in the secondary files
        sv2d param_main_s;                           //This holds the data (string) in the primary parameter file 
        dv2d param_main_d;                           //This holds the data (double) in the primary parameter file 
        iv2d param_main_i;                           //This holds the data (int) in the primary parameter file 
        sv3d param_sec_s;                            //This holds the data (string) in the secondary parameter files
        dv3d param_sec_d;                            //This holds the data (double) in the secondary parameter files
        iv3d param_sec_i;                            //This holds the data (int) in the secondary parameter files

    public:
        void get_param(string my_param_file_name,int my_num_col,int my_file_col,int my_num_col_secondary);  //read in the index
        void print_param();                                                                                 //this function reports the complex data  
        vector<string> get_column_s(int column);                                                            //returns a column from the main data file
        int main_size_x();                                                                                  //returns the size of param_main in x
        int main_size_y();                                                                                  //returns the size of param_main in y
        int sec_size_x();                                                                                   //returns the size of param_sec in x
        int sec_size_y(int z);                                                                              //returns the size of param_sec in y
        int sec_size_z();                                                                                   //returns the size of param_sec in z
        int check_file();                                                                                   //checks the integrity of the data files
        sv1d get_column_sec_s(int target_file,int target_col);                                              //returns a column from the secondary file
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function reads in a complex parameter file                                                          //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Param::get_param(string my_param_file_name,int my_num_col,int my_file_col,int my_num_col_secondary)
{
    int i = 0;          //standard variable used in loops
    int j = 0;          //standard variable used in loops
    int k = 0;          //standard variable used in loops

    //set the file information
    num_col           = my_num_col;
    file_col          = my_file_col; 
    num_col_secondary = my_num_col_secondary; 
    param_file_name   = my_param_file_name; 
   
    //read the index files
    index.get_index(param_file_name);

    //determine the number of rows
    num_row = index.index_s.size()/num_col;

    //get the number of lipid types 
    num_lip_t = index.get_column_s(num_col,file_col).size();

    //store data from main parameter file
    for(i=0; i<num_row; i++) //loop over rows
    {
        param_main_s.push_back(index.get_row_s(num_col,i));
    }

    //read in files contained in parameter file
    for(i=0; i<num_lip_t; i++) //loop over files
    {
        sv1d files = index.get_column_s(num_col,file_col);
        Index current_lip;
        current_lip.get_index(files[i]);

        int num_row_secondary = current_lip.index_s.size()/num_col_secondary;

        sv2d current_data(0);
        for(j=0; j<num_row_secondary; j++) //loop over the rows
        {
            current_data.push_back(current_lip.get_row_s(num_col_secondary,j));
        }

        param_sec_s.push_back(current_data);
    }

    //resize main structures holding ints and doubles
    param_main_i.resize(num_row);
    param_main_d.resize(num_row);
    for(i=0; i<num_row; i++) 
    {
        param_main_i[i].resize(num_col);
        param_main_d[i].resize(num_col);
    } 

    //resize secondary structures holding ints and doubles
    param_sec_i.resize(param_sec_s.size());
    param_sec_d.resize(param_sec_s.size());
    for(i=0; i<param_sec_s.size(); i++)
    {
        param_sec_i[i].resize(param_sec_s[i].size());
        param_sec_d[i].resize(param_sec_s[i].size());
        for(j=0; j<param_sec_s[i].size(); j++)
        {
            param_sec_i[i][j].resize(param_sec_s[i][j].size());
            param_sec_d[i][j].resize(param_sec_s[i][j].size());
        }
    }

    //convert main data into int and double
    for(i=0; i<num_row; i++) //loop over the rows
    {
        for(j=0; j<num_col; j++) //loop over the columns
        {
            param_main_i[i][j] = atoi(param_main_s[i][j].c_str());
            param_main_d[i][j] = atof(param_main_s[i][j].c_str());
        }
    }

    //convert secondary data into int and double
    for(i=0; i<num_lip_t; i++) //loop over files
    {
        for(j=0; j<param_sec_s[i].size(); j++) //loop over columns 
        {
            for(k=0; k<param_sec_s[i][j].size(); k++) //loop over items in row
            {
                param_sec_i[i][j][k] = atoi(param_sec_s[i][j][k].c_str());
                param_sec_d[i][j][k] = atof(param_sec_s[i][j][k].c_str());
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function reports the data taken from a complex parameter file                                       //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Param::print_param()
{
    int i = 0;          //standard variable used in loops
    int j = 0;          //standard variable used in loops
    int k = 0;          //standard variable used in loops
    
    //report the priary data
    for(i=0; i<num_row; i++) //loop over columns
    {
        for(j=0; j<num_col; j++) //loop over row 
        {
            printf(" %5s ",param_main_s[i][j].c_str());
        }
        printf("\n");
    }
    printf("\n");

    //report the secondary data
    for(i=0; i<num_lip_t; i++) //loop over files
    {
        for(j=0; j<param_sec_s[i].size(); j++) //loop over columns 
        {
            for(k=0; k<param_sec_s[i][j].size(); k++) //loop over items in row
            {
                printf(" %5s ",param_sec_s[i][j][k].c_str());
            }
            printf("\n");
        }
        printf("\n");
    }
    printf("\n");
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns a column from the main data file                                                   //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<string> Param::get_column_s(int column)
{
    int i = 0;
 
    sv1d my_vec(0);

    for(i=0; i<param_main_s.size(); i++)
    {
        my_vec.push_back(param_main_s[i][column]);
    }

    return my_vec;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns the size of param_main in x                                                        //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Param::main_size_x()
{
    return num_col; 
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns the size of param_main in y                                                        //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Param::main_size_y()
{
    return num_row; 
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns the size of param_secondary in x                                                   //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Param::sec_size_x()
{
    return num_col_secondary;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns the size of param_secondary in y                                                   //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Param::sec_size_y(int z)
{
    return param_sec_s[z].size();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function returns the size of param_secondary in z                                                   //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Param::sec_size_z()
{
    return param_sec_s.size();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function checks that the files contain the right amount of data                                     //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Param::check_file()
{
    int pass = 1;
    int i    = 0;

    //create an index to read the main file
    Index index_main;
    
    //read the index files
    index_main.get_index(param_file_name);

    if(index_main.index_s.size()%num_col != 0) 
    {
        printf("%s appears to be currupt. The number of items in this file should be a multiple of %d. \n",param_file_name.c_str(),num_col);
        pass = 0;
    }     

    //examine files contained in parameter file
    for(i=0; i<num_lip_t; i++) //loop over files
    {
        sv1d files = index.get_column_s(num_col,file_col);
        Index current_lip;
        current_lip.get_index(files[i]);
 
        if(current_lip.index_s.size()%num_col_secondary != 0)
        {
            printf("%s appears to be currupt. The number of items in this file should be a multiple of %d. \n",files[i].c_str(),num_col_secondary);
            pass = 0;
        }
    }
 
    return pass;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function checks that the files contain the right amount of data                                     //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
sv1d Param::get_column_sec_s(int target_file,int target_col)
{
    int i = 0;

    sv1d my_vec(0);

    for(i=0; i<sec_size_y(target_file); i++) //loop over y
    {
        my_vec.push_back(param_sec_s[target_file][i][target_col]);
    }
    return my_vec;
}



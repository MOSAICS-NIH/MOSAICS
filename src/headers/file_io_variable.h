
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This is a class for working with the text files with variable line sizes                                  //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Data_file_vari
{
    private:
        string file_name;             //Name of the input file

    public:
        sv2d data_s{};                  //vector holding the data (string)

    public:
        int get_data(string in_file_name);                                                     //read in the data
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function sets the format of the grid                                                                 //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Data_file_vari::get_data(string in_file_name)
{
    int i           = 0;                //standard variable used in loops
    int j           = 0;                //standard variable used in loops
    int k           = 0;                //standard variable used in loops
    int string_size = 20000;            //size of the c string
    int result      = 0;                //tells if the file was read or not
    char my_string[string_size];        //c string
         
    file_name  = in_file_name;         //name of file to be read

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Read in the data                                                                                          //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int num_lines = count_lines(file_name);

    FILE *in_file = fopen(file_name.c_str(), "r");
    if(in_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",file_name.c_str());
    }
    else
    {
        string this_line;      //store the current line

        for(i=0; i<num_lines; i++) //loop over lines
        {
            //read current line  
            read_line(in_file,this_line);

            //parse items in the line
            sv1d items = parse_string(this_line);

            //add items to data_s
            data_s.push_back(items);

	    result = 1;
	}
    }
 
    return result;     
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //  
// This function will characterize a grid in matrix format                                                   //  
//                                                                                                           //  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void characterize_grid_matrix(string in_file_name,int *y_size,int *capacity,int *x_size)
{
    int number_of_lines = 0;
    int cap             = 0;
    int items_per_line  = 0;
    int c               = 0;            //Used in fgetc
    char my_string[200];                //String to hold read in data entries
 
    //open file for reading
    FILE *in_file = fopen(in_file_name.c_str(), "r");
    if(in_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",in_file_name.c_str());
       exit(EXIT_SUCCESS);
    }
    else
    {
        //here we get the number of lines in the data file
        do
        {
          c = fgetc(in_file);
          if (c == '\n')
          {
              number_of_lines = number_of_lines + 1;
          }
        }
        while (c != EOF);
        rewind(in_file);

        //here we get the number of items in the file 
        while(fscanf(in_file, "%s,", my_string) == 1)
        {
            cap = cap + 1;
        }
        rewind(in_file);
        items_per_line = (cap)/(number_of_lines);

        //close the data file
        fclose(in_file);

        //set the dimensions of the grid
        *x_size   = items_per_line;
        *y_size   = number_of_lines;
        *capacity = cap;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function will characterize a grid in vector format                                                   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void characterize_grid_vector(string in_file_name,int *capacity,int *y_size,int *x_size)
{
    int number_of_lines = 0;
    int cap             = 0;
    int lines_per_block = 0;
    int blocks          = 0;
    int lines_data      = 0;                //How many lines of data in the vector file
    int empty_lines     = 0;                //How many empty lines in the vector file
    int c               = 0;                //Used in fgetc
    char my_string[200];                    //String to hold read in data entries

    //open file for reading
    FILE *in_file = fopen(in_file_name.c_str(), "r");
    if(in_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",in_file_name.c_str());
       exit(EXIT_SUCCESS);
    }
    else
    {
        //here we get the number of lines in the data file 
        do
        {
          c = fgetc(in_file);
          if (c == '\n')
          {
              number_of_lines = number_of_lines + 1;
          }
        }
        while (c != EOF);
        rewind(in_file);

        //count how many items in the file
        while(fscanf(in_file, "%s,", my_string) == 1)
        {
            cap = cap + 1;
        }
        rewind(in_file);

        //count the number of blocks and how many lines per block
        lines_data = (cap)/3;   
        empty_lines = number_of_lines - lines_data;  
        lines_per_block = (lines_data)/(empty_lines);     
        blocks = ((cap)/3)/(lines_per_block);  

        //set the grid dimensions
        *x_size   = blocks;
        *y_size   = lines_per_block; 
        *capacity = cap;       
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads in grid data for the matrix format                                                    //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void read_grid_matrix(int y_size,int x_size,string in_file_name,dv2d &grid,int nan_val,iv2d &nan)
{
    int i = 0;                          //Standard variable used in loops
    int j = 0;                          //Standard variable used in loops
    int k = 0;                          //Standard variable used in loops
    char my_string[200];                //String to hold read in data entries

    //open file for reading
    FILE *in_file = fopen(in_file_name.c_str(), "r");
    if(in_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",in_file_name.c_str());
        exit(EXIT_SUCCESS);
    }
    else
    {
        for(i=0; i<y_size; i++) //loop over y
        {
            for(j=0; j<x_size; j++) //loop over x
            {
                //for(k=0; k<200; k++)
                //{
                //    my_string[k] = ' ';
                //}
                fscanf(in_file, "%s,", my_string);
                if(strcmp(my_string, "NaN") == 0)
                {
                    grid[j][i] = nan_val;
                    nan[j][i]  = 1;
                }
                else
                {
                    grid[j][i] = atof(my_string);
                    nan[j][i]  = 0;
                }
            }
        }
        fclose(in_file);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads in grid data for the vector format                                                    //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void read_grid_vector(int x_size,int y_size,string in_file_name,dv2d &grid,dv2d &grid_x,dv2d &grid_y,int nan_val,iv2d &nan)
{
    int i = 0;                          //Standard variable used in loops
    int j = 0;                          //Standard variable used in loops
    int k = 0;                          //Standard variable used in loops
    int l = 0;                          //Standard variable used in loops
    char my_string[200];                //String to hold read in data entries

    //open file for reading
    FILE *in_file = fopen(in_file_name.c_str(), "r");
    if(in_file == NULL)
    {   
        printf("failure opening %s. Make sure the file exists. \n",in_file_name.c_str());
       exit(EXIT_SUCCESS);
    }
    else
    {   
        for(i=0; i<x_size; i++) //loop over each block (x changes)
        {
            for(j=0; j<y_size; j++) //loop over lines in the block (y changes)
            {
                for(k=0; k<3; k++) //loop over items in the line (x,y,z)
                {
                    //for(l=0; l<200; l++)
                    //{
                    //    my_string[l] = ' ';
                    //}
                    fscanf(in_file, "%s,", my_string);

                    if(k == 0) //x
                    {
                        grid_x[i][j] = atof(my_string);
                    }
                    else if(k == 1) //y
                    {
                        grid_y[i][j] = atof(my_string);
                    }
                    else //z
                    {
                        if(strcmp(my_string, "NaN") == 0) //NaN detected
                        {
                            grid[i][j] = nan_val;
                            nan[i][j]  = 1;
                        }
                        else
                        {
                            grid[i][j] = atof(my_string);
                            nan[i][j]  = 0;
                        }
                    }
                }
            }
        }
        fclose(in_file);
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes the grid to an output file for matrix format.                                        //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_grid_to_file_matrix(int y_size,int x_size,string out_file_name,dv2d &grid,iv2d &nan)
{   
    int i=0;            //standard variable used in loops
    int j=0;            //standard variable used in loops
    
    FILE *out_file = fopen(out_file_name.c_str(), "w");
    if(out_file == NULL)
    {   
        printf("failure opening %s. Make sure the file exists. \n",out_file_name.c_str());
    }
    else
    { 
        for(i=0; i<y_size; i++) //loop over y-dimension
        {   
            for(j=0; j<x_size; j++) //loop over x-dimension
            {   
                if(nan[j][i] == 0)
                {   
                    fprintf(out_file," %10.6f",grid[j][i]);
                }
                else //data excluded
                {   
                    fprintf(out_file," %10s ","NaN");
                }
            }
            fprintf(out_file,"\n");
        }
        fclose(out_file);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes the grid to an output file for vector format.                                        //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_grid_to_file_vector(int x_size,int y_size,string my_file_name,dv2d &grid,dv2d &grid_x,dv2d &grid_y,iv2d &nan)
{
    int i=0;            //standard variable used in loops
    int j=0;            //standard variable used in loops

    FILE *my_file = fopen(my_file_name.c_str(), "w");
    if(my_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",my_file_name.c_str());
    }
    else
    {   
        for(i=0; i<x_size; i++) //loop over each block (x changes)
        {   
            for(j=0; j<y_size; j++) //loop over lines in the block (y changes)
            { 
                fprintf(my_file," %10.3f ",grid_x[i][j]); 
                fprintf(my_file," %10.3f ",grid_y[i][j]);      
 
                if(nan[i][j] == 0) //passes cutoff
                {
                    fprintf(my_file," %10.3f  \n",grid[i][j]);   
                }
                else //data excluded
                {
                    fprintf(my_file," %10s \n","NaN");
                }
            }
            fprintf(my_file,"\n");
        }
        fclose(my_file);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function excluded insignificant data for the matrix format.                                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exclude_grid_data_vector(int x_size,int y_size,double cutoff,double avg_rho,dv2d &rho,iv2d &nan)
{
    int i=0;            //standard variable used in loops
    int j=0;            //standard variable used in loops
    int k=0;            //standard variable used in loops

    for(i=0; i<x_size; i++) //loop over each block (x changes)
    {   
        for(j=0; j<y_size; j++) //loop over lines in the block (y changes)
        {   
           if(rho[i][j] > cutoff*avg_rho)
           {
               nan[i][j] = 0;
           }
           else //data excluded
           {
               nan[i][j] = 1;
           }
        }  
    }   
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Main function checks if the grid point pass the rectangular selection critereon (vector)                  //
// No much different from the matrix version but changed i,j to x,y and from int to double.                  //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int check_rectangular_pass_vector(int invert,int target_x,int target_y,int range_x,int range_y,double x,double y)
{   
    int rectangular_selection_pass = 0;

    int in_rectangle  = y <= target_y + range_y && y >= target_y - range_y  && x <= target_x + range_x && x >= target_x - range_x;
 
    if(invert == 0) 
    {   
        if(in_rectangle == 1)
        {   
            rectangular_selection_pass = 1;
        }
        else
        {   
            rectangular_selection_pass = 0;
        }
    }
    else if(invert == 1)
    {   
        if(in_rectangle == 1)
        {   
            rectangular_selection_pass = 0;
        }
        else
        {   
            rectangular_selection_pass = 1;
        }
    }
    
    return rectangular_selection_pass;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function takes a masking file for the protein and makes a mask a distance from the protein (vector)  //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void distance_projection_vector(int x_size,int y_size,int invert,int target_x,int target_y,int range_x,int range_y,
                                double cell_size,int current_iteration,double range,double res,dv2d &init_mask,dv2d &mask,
				dv2d &mask_x,dv2d &mask_y,iv2d &mask_nan,int cumulative)
{
    int i = 0;       //standard variable used in loops
    int j = 0;       //standard variable used in loops
    int k = 0;       //standard variable used in loops
    int l = 0;       //standard variable used in loops

    for(i=0; i<x_size; i++) //loop over x
    {
        for(j=0; j<y_size; j++) //loop over y
        {
            //check rectangular selection critereon
            int rectangular_selection_pass = check_rectangular_pass_vector(invert,target_x,target_y,range_x,range_y,i,j);

            if(rectangular_selection_pass == 1) //passes rectangular selection critereon
            {
                double min_dist = 999999999;                   //The minimum distance of the grid point to the initial mask selection
                double max = current_iteration*res + range;    //The Largest distance from the initial mask selection allowed for the current iteration
                double min = current_iteration*res - range;    //The smallest distance from the initial mask selection allowed for the current iteration
                if(cumulative == 1)
                {
                    min = 0;
                }

                //compute the minimum distance of the current grid point (i,j) from the initial masking selection
                for(k=0; k<x_size; k++) //loop over x
                {
                    for(l=0; l<y_size; l++) //loop over y
                    {
                        if(init_mask[k][l] == 1) //check initial mask
                        {
                            double dx = (double)(i-k)*cell_size;   //distance in x from initial masking selection
                            double dy = (double)(j-l)*cell_size;   //distance in y from initial masking selection

                            double dist = sqrt(dx*dx + dy*dy);     //distance in x,y from initial masking selection

                            if(dist < min_dist) //update min_dist
                            {
                                min_dist = dist;
                            }
                        }
                    }
                }
                //printf("min_dist %f min %f max %f range %f \n",min_dist,min,max,range);

                //check if the grid point (i,j) falls in the acceptable range to be included in the grid selection for the current iteration
                if(min_dist < max && min_dist > min && min_dist != 0) //min_dist is in the acceptable range
                {
                    mask[i][j] = 1;
                }
                else //min_dist is not in the acceptable range
                {
                    mask[i][j] = 0;
                }
            }
            else //grid point falls outside the rectangular selection
            {
                mask[i][j] = 0;
            }

            //set x and y
            mask_x[i][j] = i*cell_size; 
            mask_y[i][j] = j*cell_size;

            //set the nan flag
            mask_nan[i][j] = 0;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This is a class for working with the grid.                                                                //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Grid_io
{       
    private:
        string file_name;        //Name of the input file
        int capacity   = 0;      //Number of items in the input files
        int x_size     = 0;      //Number of lines in matrix file
        int y_size     = 0;      //Number of items per line in matrix file
        int odf        = 0;      //Data file format
        double nan_val = 0.0;    //Value added to grid when NaN is encountered

    public:
        dv2d grid{};                  //vector holding grid data 
        dv2d grid_x{};                //vector holding grid data x
        dv2d grid_y{};                //vector holding grid data y
        iv2d nan{};                   //vector holding grid data NaN

    public:
        void set_format(int format);                                                                              //set the grid format
        void get_grid(string grid_file_name);                                                                     //read in grid data
        void write_grid(string out_file_name);                                                                    //write out grid data
        int size_x();                                                                                             //return the size in x
        int size_y();                                                                                             //return the size in y
        int size();                                                                                               //return how many items in file
        void size_grid(int size_x,int size_y);                                                                    //allocates memory for the grid
        void exclude_grid_data(double cutoff,double avg_rho,dv2d &rho);                                           //exclude insignificant data  
        void init_grid(double val);                                                                               //set grid point equal to val
        void set_nan(double val);                                                                                 //set nan tags equal to val
        double get_average();                                                                                     //return the average over the grid
        double get_sum();                                                                                         //return the sum over the grid
        void distance_projection(int invert,int target_x,int target_y,int range_x,int range_y,double cell_size,   //computes a mask around another mask
                                 int current_iteration,double range,double res,dv2d &init_mask,int cumulative);
        void print_dim(int b_header);                                                                             //prints the dimensions of grid
        void grow_mask(dv2d &data,int target_x,int target_y,double val);                                          //grows a mask around a point 
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function sets the format of the grid                                                                 //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_io::set_format(int format)
{
    odf = format; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads in the grid data.                                                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_io::get_grid(string grid_file_name)
{
    file_name = grid_file_name; 

    if(odf == 0) //matrix format 
    {
        characterize_grid_matrix(grid_file_name,&y_size,&capacity,&x_size);

        grid.resize(x_size,dv1d(y_size,0.0));
        nan.resize(x_size,iv1d(y_size,0));

        read_grid_matrix(y_size,x_size,grid_file_name,grid,nan_val,nan);
    }
    else if(odf == 1) //vector format
    {
        characterize_grid_vector(grid_file_name,&capacity,&y_size,&x_size);
    
	grid.resize(x_size,dv1d(y_size,0.0));
        grid_x.resize(x_size,dv1d(y_size,0.0));
        grid_y.resize(x_size,dv1d(y_size,0.0));
        nan.resize(x_size,iv1d(y_size,0));
    
        read_grid_vector(x_size,y_size,grid_file_name,grid,grid_x,grid_y,nan_val,nan);
    }
}   

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes the grid to an output file                                                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_io::write_grid(string out_file_name)
{
    if(odf == 0) //matrix format
    {
        write_grid_to_file_matrix(y_size,x_size,out_file_name,grid,nan);
    }
    else if(odf == 1) //vector format
    {
        write_grid_to_file_vector(x_size,y_size,out_file_name,grid,grid_x,grid_y,nan);
    }
}   

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the size in x                                                                       //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Grid_io::size_x()
{
    return x_size;
}   

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the size in y                                                                       //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Grid_io::size_y()
{
    return y_size;
}   

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns how many items in file                                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Grid_io::size()
{
    return capacity;
}   

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function allocates memory for the grid                                                               //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_io::size_grid(int size_x,int size_y)
{
    x_size   = size_x; 
    y_size   = size_y; 
    capacity = size_x*size_y; 

    grid.resize(x_size,dv1d(y_size,0.0));
    nan.resize(x_size,iv1d(y_size,0.0));

    if(odf == 1) //vector format
    {
        grid_x.resize(x_size,dv1d(y_size,0.0));
        grid_y.resize(x_size,dv1d(y_size,0.0));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the average over the grid                                                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Grid_io::get_average()
{
    int i      = 0;
    int j      = 0;
    double avg = 0;

    for(i=0; i<x_size; i++) //loop over each block (x changes)
    {
        for(j=0; j<y_size; j++) //loop over lines in the block (y changes)
        {
            avg = avg + grid[i][j];
        }
    }
    avg = avg/((double)x_size*(double)y_size);

    return avg;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the sum over the grid                                                               //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Grid_io::get_sum()
{
    int i      = 0;
    int j      = 0;
    double sum = 0;

    for(i=0; i<x_size; i++) //loop over each block (x changes)
    {
        for(j=0; j<y_size; j++) //loop over lines in the block (y changes)
        {
            sum = sum + grid[i][j];
        }
    }

    return sum;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function excluded insignificant data                                                                 //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_io::exclude_grid_data(double cutoff,double avg_rho,dv2d &rho)
{
    exclude_grid_data_vector(x_size,y_size,cutoff,avg_rho,rho,nan);
}   

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function initializes the z value of the grid                                                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_io::init_grid(double val)
{
    int i = 0;
    int j = 0;
    
    for(i=0; i<x_size; i++) //loop over each block (x changes)
    {
        for(j=0; j<y_size; j++) //loop over lines in the block (y changes)
        {
            grid[i][j] = val; 
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function initializes the nan value of the grid                                                       //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_io::set_nan(double val)
{
    int i = 0;
    int j = 0;
 
    for(i=0; i<x_size; i++) //loop over each block (x changes)
    {
        for(j=0; j<y_size; j++) //loop over lines in the block (y changes)
        {
            grid[i][j] = val;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes a mask a certain distance from the initial mask                                    //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_io::distance_projection(int invert,int target_x,int target_y,int range_x,int range_y,double cell_size,int current_iteration,double range,
                                      double res,dv2d &init_mask,int cumulative)
{
    distance_projection_vector(x_size,y_size,invert,target_x,target_y,range_x,range_y,cell_size,current_iteration,range,res,
                               init_mask,grid,grid_x,grid_y,nan,cumulative);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints info about the grid dimensions                                                       //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_io::print_dim(int b_header)
{
    if(b_header == 1)
    {
        printf("%30s %20s %20s %20s \n","","size_x","size_y","size");
        printf("%30s-%20s-%20s-%20s \n","------------------------------","--------------------","--------------------","--------------------");
    }
    printf("%30s %20d %20d %20d \n",file_name.c_str(),x_size,y_size,capacity);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function grows a mask around a point                                                                 //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid_io::grow_mask(dv2d &data,int target_x,int target_y,double val)
{
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int mask_grew = 0;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Plant seed                                                                                                //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(data[target_x][target_y] == val)
    {
        grid[target_x][target_y] = 1;
        mask_grew = 1;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Grow the mask                                                                                             //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    while(mask_grew == 1)
    {
        mask_grew = 0;

        for(i=0; i<x_size; i++) //loop over x
        {
            for(j=0; j<y_size; j++) //loop over y
            {
                if(grid[i][j] == 1.0) ////grid point has a 1
                {
                    for(k=i-1; k<=i+1; k++) //check neighboring x
                    {
                        for(l=j-1; l<=j+1; l++) //check neighboring y
                        {
                            if(data[k][l] == val)
                            {
                                if(grid[k][l] == 0) //add the grid point
                                {
                                    grid[k][l] = 1;
                                    mask_grew = 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function compares 2 grids                                                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int comp_grid(Grid_io g1,Grid_io g2)
{
    int same = 1;

    if( g1.size_x() !=  g2.size_x())
    {
        same = 0;
    }
    if(g1.size_y() != g2.size_y())
    {
        same = 0;
    }
    if(g1.size() != g2.size())
    {
        same = 0;
    }
   
    return same;  
}


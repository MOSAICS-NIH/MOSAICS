
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This is a class for binning data of ints                                                                 //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Histogram_i
{
    private: 
        int bin_width;                                       //The bin width
        int size;                                            //How many data points making sample
        int start;                                           //Starting value for loop variable
        int end;                                             //End value for loop variable
        int num_bins;                                        //How many bins
        int b_force_range = 0;                               //Forces the range covered by histogram
        int forced_smallest;                                 //Smallest value to be covered (forced)
        int forced_largest;                                  //Largest value to be covered (forced)

    public:
        iv1d bins;                                           //Holds the number of hits for each bin

    public:
        void bin_data(iv1d &data,int bin_width);             //Analyze data and make a histogram
        void write_histo(string out_file_name,string label); //Write the histogram to file
        void set_range(int min,int max);                     //Manually set the smallest and largest values to include in histo
        double get_avg(iv1d &data);                          //Computes average over the data set
        double get_stdev(iv1d &data);                        //Computes standard deviation over the data set
        int    get_bin_width();                              //returns the bin width
        int get_size();                                      //returns the number of data points analyzed
        int get_start();                                     //returns starting value for loop value
        int get_end();                                       //returns end value for loop variable
        int get_num_bins();                                  //returns the number of bins
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function bins the data                                                                               //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Histogram_i::bin_data(iv1d &data,int my_bin_width)
{
    int i = 0;
    int j = 0;

    bin_width = my_bin_width; 
    size      = data.size();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Find largest value                                                                                        //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int largest_value = 0;
    for(i=0; i<size; i++) //loop over data
    {
        if(data[i] > largest_value)
        {
            largest_value = data[i];
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Find smallest value                                                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int smallest_value = 0;
    for(i=0; i<size; i++) //loop over data
    {
        if(data[i] < smallest_value)
        {
            smallest_value = data[i];
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Bin data                                                                                                  //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    num_bins  = (int)ceil((largest_value-smallest_value)/(bin_width));
    start     = floor(smallest_value/bin_width);
    end       = start + num_bins;
    int count = 0;

    bins.resize(num_bins,0);  

    for(i=0; i<size; i++) //loop over data
    {
        for(j=start; j<end; j++) //loop over num bins
        {
            if(data[i] >= j*(bin_width) && data[i] < (j+1)*(bin_width))
            {
                bins[count] = bins[count] + 1;
            }
            count++;
        }
        count = 0;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //  
// This function writes a histogram to output file                                                           //  
//                                                                                                           //  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Histogram_i::write_histo(string out_file_name,string label)
{
    int i     = 0;
    int count = 0;

    FILE *out_file = fopen(out_file_name.c_str(), "w");
    if(out_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",out_file_name.c_str());
    }
    else
    {
        double sum_prob    = 0;

        fprintf(out_file," %10s #%15s \n","#Column_1",label.c_str());
        fprintf(out_file," %10s #%15s \n","#Column_2","Count");
        fprintf(out_file," %10s #%15s \n","#Column_3","Probability");

        for(i=start; i<end; i++)
        {
            double probability = (double)bins[count]/(double)size;
            fprintf(out_file,"%10.4f %10d %10.8f \n",(double)i*bin_width,bins[count],probability);
            sum_prob = sum_prob + probability;
            count++;
        }
        fclose(out_file);
        printf("Summing probability over bins. sum_prob = %f \n",sum_prob);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function manually sets the range covered in the histogram                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Histogram_i::set_range(int min,int max)
{
    b_force_range   = 1;
    forced_smallest = min;
    forced_largest  = max;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the average over the data set                                                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Histogram_i::get_avg(iv1d &data)
{
    int i      = 0;
    double avg = 0;

    for(i=0; i<data.size(); i++) //loop over data
    {
        avg = avg + (double)data[i];
    }
    avg = avg/(double)data.size();

    return avg;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the standard deviation over the data set                                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Histogram_i::get_stdev(iv1d &data)
{
    int i = 0;

    //compute the average
    double avg = 0;
    for(i=0; i<data.size(); i++) //loop over data
    {
        avg = avg + (double)data[i];
    }
    avg = avg/(double)data.size();

    //compute standard deviation 
    double stdev = 0;
    for(i=0; i<data.size(); i++) //loop over the data
    {
        stdev = stdev + pow(((double)data[i] - avg),2);
    }
    stdev = stdev/(double)(data.size()-1);
    stdev = sqrt(stdev);

    return stdev;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the bin width                                                                       //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Histogram_i::get_bin_width()
{
    return bin_width; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the number of data points analyzed                                                  //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Histogram_i::get_size()
{
    return size;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the starting value for loop value                                                   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Histogram_i::get_start()
{
    return start;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the end value for loop value                                                        //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Histogram_i::get_end()
{
    return end;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the number of bins                                                                  //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Histogram_i::get_num_bins()
{
    return num_bins;
}






//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This is a class for binning data of doubles                                                              //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Histogram_d
{
    private:
        int     size;                                        //How many data points making sample
        int     start;                                       //Starting value for loop variable
        int     end;                                         //End value for loop variable
        int     num_bins;                                    //How many bins
        int     b_force_range = 0;                           //Forces the range covered by histogram
        double  bin_width;                                   //The bin width
        double  forced_smallest;                             //Smallest value to be covered (forced)
        double  forced_largest;                              //Largest value to be covered (forced)

    public:
        dv1d bins;                                           //Holds the number of hits for each bin

    public:
        void bin_data(dv1d &data,double bin_width);          //Analyze data and make a histogram
        void write_histo(string out_file_name,string label); //Write the histogram to file
        void set_range(double min,double max);               //Manually set the smallest and largest values to include in histo
        double get_avg(dv1d &data);                          //Computes average over the data set
        double get_stdev(dv1d &data);                        //Computes standard deviation over the data set
        double get_bin_width();                              //returns the bin width
        int get_size();                                      //returns the number of data points analyzed
        int get_start();                                     //returns starting value for loop value
        int get_end();                                       //returns end value for loop variable
        int get_num_bins();                                  //returns the number of bins
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function bins the data                                                                               //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Histogram_d::bin_data(dv1d &data,double my_bin_width)
{
    int i = 0;
    int j = 0;

    bin_width = my_bin_width;
    size      = data.size();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Find largest value                                                                                        //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double largest_value = 0;
    for(i=0; i<size; i++) //loop over data
    {
        if(data[i] > largest_value)
        {
            largest_value = data[i];
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Find smallest value                                                                                       //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double smallest_value = 0;
    for(i=0; i<size; i++) //loop over data
    {
        if(data[i] < smallest_value)
        {
            smallest_value = data[i];
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Check if the range is set manually                                                                        //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(b_force_range == 1)
    {
        smallest_value = forced_smallest;
        largest_value  = forced_largest;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                           //
    // Bin data                                                                                                  //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    num_bins  = (int)ceil((largest_value-smallest_value)/(bin_width));
    start     = floor(smallest_value/bin_width);
    end       = start + num_bins;
    int count = 0;

    bins.resize(num_bins,0.0);

    for(i=0; i<size; i++) //loop over data
    {
        for(j=start; j<end; j++) //loop over num bins
        {
            if(data[i] >= (double)j*(bin_width) && data[i] < (double)(j+1)*(bin_width))
            {
                bins[count] = bins[count] + 1.0;
            }
            count++;
        }
        count = 0;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes a histogram to output file                                                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Histogram_d::write_histo(string out_file_name,string label)
{
    int i     = 0;
    int count = 0;

    FILE *out_file = fopen(out_file_name.c_str(), "w");
    if(out_file == NULL)
    {
        printf("failure opening %s. Make sure the file exists. \n",out_file_name.c_str());
    }
    else
    {
        double sum_prob    = 0;

        fprintf(out_file," %10s #%15s \n","#Column_1",label.c_str());
        fprintf(out_file," %10s #%15s \n","#Column_2","Count");
        fprintf(out_file," %10s #%15s \n","#Column_3","Probability");

        for(i=start; i<end; i++)
        {
            double probability = bins[count]/(double)size;
            fprintf(out_file,"%10.4f %10d %10.8f \n",(double)i*bin_width,(int)bins[count],probability);
            sum_prob = sum_prob + probability;
            count++;
        }
        fclose(out_file);
        printf("Summing probability over bins. sum_prob = %f \n",sum_prob);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function manually sets the range covered in the histogram                                            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Histogram_d::set_range(double min,double max) 
{
    b_force_range   = 1;
    forced_smallest = min; 
    forced_largest  = max;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the average over the data set                                                      //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Histogram_d::get_avg(dv1d &data)
{
    int i      = 0;
    double avg = 0;

    for(i=0; i<data.size(); i++) //loop over data
    {
        avg = avg + data[i];
    }
    avg = avg/(double)data.size();

    return avg;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function computes the standard deviation over the data set                                           //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Histogram_d::get_stdev(dv1d &data)
{
    int i = 0;
 
    //compute the average
    double avg = 0;
    for(i=0; i<data.size(); i++) //loop over data
    {
        avg = avg + data[i];
    }
    avg = avg/(double)data.size();

    //compute standard deviation 
    double stdev = 0;
    for(i=0; i<data.size(); i++) //loop over the data
    {
        stdev = stdev + pow((data[i] - avg),2);
    }
    stdev = stdev/(double)(data.size()-1);
    stdev = sqrt(stdev);

    return stdev;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the bin width                                                                       //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Histogram_d::get_bin_width()
{
    return bin_width;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the number of data points analyzed                                                  //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Histogram_d::get_size()
{
    return size;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the starting value for loop value                                                   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Histogram_d::get_start()
{
    return start;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the end value for loop value                                                        //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Histogram_d::get_end()
{
    return end;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function returns the number of bins                                                                  //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Histogram_d::get_num_bins()
{
    return num_bins;
}


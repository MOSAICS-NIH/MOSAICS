
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This is a class for binning data of ints                                                                 //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Histogram_i
{
    private: 
        int  bin_width;                                      //The bin width
        int  size;                                           //How many data points making sample
        int  start;                                          //Starting value for loop variable
        int  end;                                            //End value for loop variable
        int  num_bins;                                       //How many bins
    public:
        iv2d bins;                                           //Holds the value, number of hits and probability for each bin

    public:
        void bin_data(iv1d &data,int bin_width);             //Analyze data and make a histogram
        void write_histo(string out_file_name,string label); //Write the histogram to file
        double get_avg(iv1d &data);                          //Computes average over the data set
        double get_stdev(iv1d &data);                        //Computes standard deviation over the data set
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

    bins.resize(num_bins);  
    for(i=0; i<num_bins; i++)
    {
        bins[i].resize(2,0);
    } 

    for(i=0; i<size; i++) //loop over data
    {
        for(j=start; j<end; j++) //loop over num bins
        {
            if(data[i] >= j*(bin_width) && data[i] < (j+1)*(bin_width))
            {
                bins[count][1] = bins[count][1] + 1;
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
            double probability = (double)bins[count][1]/(double)size;
            fprintf(out_file,"%10.4f %10d %10.8f \n",(double)i*bin_width,bins[count][1],probability);
            sum_prob = sum_prob + probability;
            count++;
        }
        fclose(out_file);
        printf("Summing probability over bins. sum_prob = %f \n",sum_prob);
    }
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
        double  bin_width;                                   //The bin width

    public:
        dv2d bins;                                           //Holds the value, number of hits and probability for each bin

    public:
        void bin_data(dv1d &data,double bin_width);          //Analyze data and make a histogram
        void write_histo(string out_file_name,string label); //Write the histogram to file
        double get_avg(dv1d &data);                          //Computes average over the data set
        double get_stdev(dv1d &data);                        //Computes standard deviation over the data set
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
    // Bin data                                                                                                  //
    //                                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    num_bins  = (int)ceil((largest_value-smallest_value)/(bin_width));
    start     = floor(smallest_value/bin_width);
    end       = start + num_bins;
    int count = 0;

    bins.resize(num_bins);
    for(i=0; i<num_bins; i++)
    {
        bins[i].resize(2,0.0);
    }

    for(i=0; i<size; i++) //loop over data
    {
        for(j=start; j<end; j++) //loop over num bins
        {
            if(data[i] >= (double)j*(bin_width) && data[i] < (double)(j+1)*(bin_width))
            {
                bins[count][1] = bins[count][1] + 1.0;
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
            double probability = bins[count][1]/(double)size;
            fprintf(out_file,"%10.4f %10d %10.8f \n",(double)i*bin_width,(int)bins[count][1],probability);
            sum_prob = sum_prob + probability;
            count++;
        }
        fclose(out_file);
        printf("Summing probability over bins. sum_prob = %f \n",sum_prob);
    }
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

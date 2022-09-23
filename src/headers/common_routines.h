
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// Print the command line arguments for records                                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void print_arguments(int argc, const char * argv[])
{
    int i = 0;

    printf("Input: ");
    for(i=0; i<argc; i++)
    {
        printf("%s ",argv[i]);
    }
    printf("\n\n");
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints the credits for the program. Also prints a record of the input arguments.            //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void print_credits(int argc, const char *argv[],string program_name)
{
    int i = 0;
    //set the length of the line
    int line_length = 130;

    //compute the padding needed
    int length = program_name.length();
    int lead_pad = ceil((line_length-4-length)/2);
    int trail_pad = line_length - 4 - length - lead_pad;

    cout << "\n";

    //line of slashes
    for(i=0; i<line_length; i++)
    {
        cout << "/";
    }
    cout << "\n";

    //2 slashes then spaces with 2 slashes at end
    cout << "//";
    for(i=0; i<line_length-4; i++)
    {
        cout << " ";
    }
    cout << "//";
    cout << "\n";

    //line with program name
    cout << "//";
    for(i=0; i<lead_pad; i++)
    {
        cout << " ";
    }
    cout << program_name;
    for(i=0; i<trail_pad; i++)
    {
        cout << " ";
    }
    cout << "//";
    cout << "\n";

    //2 slashes then spaces with 2 slashes at end
    cout << "//";
    for(i=0; i<line_length-4; i++)
    {
        cout << " ";
    }
    cout << "//";
    cout << "\n";

    //line of slashes
    for(i=0; i<line_length; i++)
    {
        cout << "/";
    }
    cout << "\n";

    //now we print the command line arguments for records
    printf("Input: ");
    for(i=0; i<argc; i++)
    {
        printf("%s ",argv[i]);
    }
    printf("\n\n");
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function checks a string and determines if it is an integer or not                                   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int check_int(string name)
{
    int i        = 0;      //standard variable used in loops
    int is_digit = 1;      //tells if the string is an int
    int min      = 0;      //set the loop lower limit when looking for non numberical chars

    //first character could be the sign
    if(name[0] == '+' || name[0] == '-') //first digit is the sign
    {
        min = 1;
    }

    //check for non numerical chars
    for(i=min; i<name.length(); i++) //loop over string
    {
        if(isdigit(name[i]) == 0) //a non numerical char was found
        {
            is_digit = 0;
        }
    }

    return is_digit;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function checks a string and determines if it is an floating point number or not                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int check_float(string name)
{
    int i        = 0;      //standard variable used in loops
    int is_float = 1;      //tells if the string is a float
    int decimal  = 0;      //count decimals found     
    int numbers  = 0;      //count numbers found
    int chars    = 0;      //count non numberical chars found other than decimal
    int sign     = 0;      //was a sign provided
    int min      = 0;      //set the loop lower limit when looking for non numberical chars

    //first and last characters should be a number
    if(isdigit(name[0]) == 0) //first digit is not a number
    {
        //check for + or -
        if( (name[0] == '+' || name[0] == '-') && isdigit(name[1]) != 0)
        {
            sign = 1;
        }
        else
        {
            is_float = 0;
        }
    }
    if(isdigit(name[name.length()-1]) == 0)
    {
        is_float = 0;
    }

    //count decimals
    //count numbers
    //check for non numerical chars
    if(sign == 0) //check first char
    {
        min = 0;
    }
    else //ignore first char (it is + or -)
    {
        min = 1;
    }
    for(i=min; i<name.length(); i++) //loop over string
    {
        if(isdigit(name[i]) != 0) //a numerical char was found
        {
            numbers++;
        }
        else if(name[i] == '.') //found the decimal
        {
            decimal++;
        }
        else
        {
            chars++;
        }
    }

    //should have at least 2 numbers and a single decimal
    if(decimal != 1 || numbers < 2 || chars > 0)
    {
        is_float = 0;
    }

    return is_float;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints the program description and the header for the help options                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int check_help(int argc, const char *argv[])
{
    int i    = 0;
    int help = 0;

    for(i=0; i<argc; i++)
    {
        if(strcmp(argv[i], "-h") == 0)
        {
            help = 1;
        }
    }
    return help;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints the program description and the header for the help options                          //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void start_input_arguments(int argc, const char *argv[],string description)
{
    int i = 0;

    stringstream ss (description);
    string item;
    string first = "Description:";
    int length = 0;
    int line_length = 0;

    //print the program description
    for(i=0; i<argc; i++)
    {
        if(strcmp(argv[i], "-h") == 0)
        {
            length = first.length();
            line_length = line_length + length + 1;
            printf("%s ",first.c_str());

            while (getline (ss, item, ' '))
            {
                length = item.length();
                line_length = line_length + length + 1;

                if(line_length > /*81*/ 130)
                {
                    printf("\n");
                    printf("%s ",item.c_str());
                    line_length = length + 1;
                }
                else
                {
                    printf("%s ",item.c_str());
                }
            }
            printf("\n\n");
        }
    }

    //print the help options heerder
    for(i=0; i<argc; i++)
    {
        if(strcmp(argv[i], "-h") == 0)
        {
            printf("        To use this program please enter the following command line options. \n");
            printf("\n");
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints the final help information. Concludes the analysis of input arguments.               //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void conclude_input_arguments(int argc, const char *argv[],string program_name,sv1d &cl_tags)
{
    int i = 0;
    int j = 0;

    for(i=0; i<argc; i++)
    {
        if(strcmp(argv[i], "-h") == 0)
        {
            printf("        -h      O    S    print help options                           \n");
            printf("                                                             \n");
            printf("        Thank you for using %s.",program_name.c_str());
            printf("\n");
            printf("\n");
            exit(EXIT_SUCCESS);
        }
    }

    //check for provided tags that are not recognized by the program
    int add_line = 0;
    for(i=0; i<argc; i++)
    {
        if(argv[i][0] == '-')
        {
            int found = 0;

            for(j=0; j<cl_tags.size(); j++)
            {
                if(strcmp(argv[i], cl_tags[j].c_str()) == 0)
                {
                    found = 1;
                }
            }

            if(found == 0)
            {
                printf("Warning! Command line argument not recognized (%8s). Perhaps this is a typo? You should check the list of acceptable arguments using the -h tag. \n",argv[i]);
                add_line = 1;
            }
        }
    }
    if(add_line == 1) //print a new line so that the warning messeges stand out. 
    {
        printf("\n");
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function enables the user to add a new input argument (string). Takes a tag, a string and a          //
// description of the arguement. The function sets the sting value from the command line and prints the      //
// description if -h is called.                                                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void add_argument_s(int argc, const char *argv[],string tag,string &argument,string info,sv1d &cl_tags,int *found_tag,int b_require)
{
    int i         = 0;        //standard variable used in loops
    int tag_found = 0;        //tells if the user provided the tag
    int h_found   = 0;        //tells if the user provided -h
    string require_tag;       //tag indicating if the tag is required by the user or optional

    //determine if the tag is required by the user and set required_tag
    if(b_require == 0)
    {
        require_tag = "O";
    }
    else if(b_require == 1)
    {
        require_tag = "R";
    }

    //return whether the tag is found
    if(found_tag != nullptr)
    {
        *found_tag = 0;
    }

    //loop over aguments. look for the tag. look for -h. set argument if tag is found
    for(i=0; i<argc; i++)
    {
        if(strcmp(argv[i], tag.c_str()) == 0)
        {
            tag_found = 1;

            argument = strdup(argv[i+1]);
            if(found_tag != nullptr)
            {
                *found_tag = 1;
            }
        }
        if(strcmp(argv[i], "-h") == 0)
        {
            h_found = 1;

            printf("        %-8s%-5s%-5s%-s\n",tag.c_str(),require_tag.c_str(),"S",info.c_str());
        }
    }

    //check that a required tag was actually found
    if(b_require == 1 && tag_found == 0 && h_found == 0)
    {
        printf("To use this program you must include the %8s command line argument. \n\n",tag.c_str());
        exit(EXIT_SUCCESS);
    }

    //add tag to list of acceptable tags
    cl_tags.push_back(tag);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function enables the user to add a new input argument (int). Takes a tag, a int and an               //
// description of the arguement. The function sets the int value from the command line and prints the        //
// description if -h is called.                                                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void add_argument_i(int argc, const char *argv[],string tag,int *argument,string info,sv1d &cl_tags,int *found_tag,int b_require)
{
    int i         = 0;        //standard variable used in loops
    int tag_found = 0;        //tells if the user provided the tag
    int h_found   = 0;        //tells if the user provided -h
    string require_tag;       //tag indicating if the tag is required by the user or optional

    //determine if the tag is required by the user and set required_tag
    if(b_require == 0)
    {
        require_tag = "O";
    }
    else if(b_require == 1)
    {
        require_tag = "R";
    }

    //return whether the tag is found
    if(found_tag != nullptr)
    {
        *found_tag = 0;
    }

    //loop over aguments. look for the tag. look for -h. set argument if tag is found
    for(i=0; i<argc; i++)
    {
        if(strcmp(argv[i], tag.c_str()) == 0)
        {
            tag_found = 1;

            //check that the user provided an int
            if(check_int(argv[i+1]) == 1 || check_help(argc,argv) == 1) //an int was provided
            {
                *argument = atoi(argv[i+1]);
                if(found_tag != nullptr)
                {
                    *found_tag = 1;
                }
            }
            else //int was not provided
            {
                printf("The %s command line argument requires an integer value while %s was detected. \n\n",tag.c_str(),argv[i+1]);
                exit(EXIT_SUCCESS);
            }
        }
        if(strcmp(argv[i], "-h") == 0)
        {
            h_found = 1;

            printf("        %-8s%-5s%-5s%-s\n",tag.c_str(),require_tag.c_str(),"I",info.c_str());
        }
    }

    //check that a required tag was actually found
    if(b_require == 1 && tag_found == 0 && h_found == 0)
    {
        printf("To use this program you must include the %8s command line argument. \n\n",tag.c_str());
        exit(EXIT_SUCCESS);
    }

    //add tag to list of acceptable tags
    cl_tags.push_back(tag);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function enables the user to add a new input argument (float). Takes a tag, a float and a            //
// description of the arguement. The function sets the float value from the command line and prints the      //
// description if -h is called.                                                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void add_argument_f(int argc, const char *argv[],string tag,float *argument,string info,sv1d &cl_tags,int *found_tag,int b_require)
{
    int i         = 0;        //standard variable used in loops
    int tag_found = 0;        //tells if the user provided the tag
    int h_found   = 0;        //tells if the user provided -h
    string require_tag;       //tag indicating if the tag is required by the user or optional

    //determine if the tag is required by the user and set required_tag
    if(b_require == 0)
    {
        require_tag = "O";
    }
    else if(b_require == 1)
    {
        require_tag = "R";
    }

    //return whether the tag is found
    if(found_tag != nullptr)
    {
        *found_tag = 0;
    }

    //loop over aguments. look for the tag. look for -h. set argument if tag is found
    for(i=0; i<argc; i++)
    {
        if(strcmp(argv[i], tag.c_str()) == 0)
        {
            tag_found = 1;

            //check that the user provided a float
            if(check_float(argv[i+1]) == 1 || check_help(argc,argv) == 1) //a float was provided
            {
                *argument = atof(argv[i+1]);
                if(found_tag != nullptr)
                {
                    *found_tag = 1;
                }
            }
            else //float was not provided 
            {
                printf("The %s command line argument requires a floating point number while %s was detected. \n\n",tag.c_str(),argv[i+1]);
                exit(EXIT_SUCCESS);
            }
        }
        if(strcmp(argv[i], "-h") == 0)
        {
            h_found = 1;

            printf("        %-8s%-5s%-5s%-s\n",tag.c_str(),require_tag.c_str(),"F",info.c_str());
        }
    }

    //check that a required tag was actually found
    if(b_require == 1 && tag_found == 0 && h_found == 0)
    {
        printf("To use this program you must include the %8s command line argument. \n\n",tag.c_str());
        exit(EXIT_SUCCESS);
    }

    //add tag to list of acceptable tags
    cl_tags.push_back(tag);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function enables the user to add a new input argument (double). Takes a tag, a double and a          //
// description of the arguement. The function sets the double value from the command line and prints the     //
// description if -h is called.                                                                              //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void add_argument_d(int argc, const char *argv[],string tag,double *argument,string info,sv1d &cl_tags,int *found_tag,int b_require)
{
    int i         = 0;        //standard variable used in loops
    int tag_found = 0;        //tells if the user provided the tag
    int h_found   = 0;        //tells if the user provided -h
    string require_tag;       //tag indicating if the tag is required by the user or optional

    //determine if the tag is required by the user and set required_tag
    if(b_require == 0)
    {
        require_tag = "O";
    }
    else if(b_require == 1)
    {
        require_tag = "R";
    }

    //return whether the tag is found
    if(found_tag != nullptr)
    {
        *found_tag = 0;
    }

    //loop over aguments. look for the tag. look for -h. set argument if tag is found
    for(i=0; i<argc; i++)
    {
        if(strcmp(argv[i], tag.c_str()) == 0)
        {
            tag_found = 1;

            //check that the user provided a float
            if(check_float(argv[i+1]) == 1 || check_help(argc,argv) == 1) //a float was provided
            {
                *argument = atof(argv[i+1]);
                if(found_tag != nullptr)
                {
                    *found_tag = 1;
                }
            }
            else //float was not provided 
            {
                printf("The %s command line argument requires a floating point number while %s was detected. \n\n",tag.c_str(),argv[i+1]);
                exit(EXIT_SUCCESS);
            }
        }
        if(strcmp(argv[i], "-h") == 0)
        {
            h_found = 1;

            printf("        %-8s%-5s%-5s%-s\n",tag.c_str(),require_tag.c_str(),"F",info.c_str());
        }
    }

    //check that a required tag was actually found
    if(b_require == 1 && tag_found == 0 && h_found == 0)
    {
        printf("To use this program you must include the %8s command line argument. \n\n",tag.c_str());
        exit(EXIT_SUCCESS);
    }

    //add tag to list of acceptable tags
    cl_tags.push_back(tag);
}


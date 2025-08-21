//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function takes a string and removes the .whatever from the end, adds a tag and puts back .whatever. //
// If no . is found it simply adds the tag to the end. Is also safe when using ../file.dat or ../file       //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
string add_tag(string in_file_name,string tag)
{
    string out_file_name = in_file_name;
    string crop          = "";
    int size             = out_file_name.size();
    int count            = 0;
    int i                = 0;
    int pair             = 0;
    int last_dot         = -1;

    for(i=0; i<size; i++) //loop over string
    {
        if(out_file_name[i] == '.')
        {
            pair = 0;

            if(i>0)
            {
                if(out_file_name[i-1] == '.')
                {
                    pair = 1;
                }
            }
            if(i<size-1)
            {
                if(out_file_name[i+1] == '.')
                {
                    pair = 1;
                }
            }

            if(pair == 0)
            {
                last_dot = i;
            }
        }
    }

    if(last_dot > -1) //a dot was found
    {
        while(out_file_name[size-1] != '.')
        {
            crop = out_file_name[size-1] + crop;
            out_file_name.pop_back();
            size = out_file_name.size();
        }
        crop = out_file_name[size-1] + crop;
        out_file_name.pop_back();

        out_file_name = out_file_name + tag + crop;

    }
    else //no dot was found
    {
        out_file_name = out_file_name + tag;
    }

    return out_file_name;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function takes a string and removes the .whatever from the end, adds a tag.                         //
// If no . is found it simply adds the tag to the end. Is also safe when using ../file.dat or ../file       //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
string chop_and_add_tag(string in_file_name,string tag)
{
    string out_file_name = in_file_name;
    string crop          = "";
    int size             = out_file_name.size();
    int count            = 0;
    int i                = 0;
    int pair             = 0;
    int last_dot         = -1;

    for(i=0; i<size; i++) //loop over string
    {
        if(out_file_name[i] == '.')
        {
            pair = 0;

            if(i>0)
            {
                if(out_file_name[i-1] == '.')
                {
                    pair = 1;
                }
            }
            if(i<size-1)
            {
                if(out_file_name[i+1] == '.')
                {
                    pair = 1;
                }
            }

            if(pair == 0)
            {
                last_dot = i;
            }
        }
    }

    if(last_dot > -1) //a dot was found
    {
        while(out_file_name[size-1] != '.')
        {
            crop = out_file_name[size-1] + crop;
            out_file_name.pop_back();
            size = out_file_name.size();
        }
        crop = out_file_name[size-1] + crop;
        out_file_name.pop_back();

        out_file_name = out_file_name + tag;

    }
    else //no dot was found
    {
        out_file_name = out_file_name + tag;
    }

    return out_file_name;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function takes a string and removes the filename leaving only the directory                         //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
string chop_filename(string in_file_name)
{
    string out_file_name = "";
    int i                = 0;
    int last_slash       = -1;

    for(i=0; i<in_file_name.size(); i++) //loop over string
    {
        if(in_file_name[i] == '/')
        {
            last_slash = i;
        }
    }

    if(last_slash > -1) //a slash was found
    {
        for(i=0; i <= last_slash; i++)
        {
            out_file_name.push_back(in_file_name[i]);
        }
    }
    else //no slash was found
    {
        out_file_name = "";
    }

    return out_file_name;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function takes a filename and check that it contains the correct extension                          //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void check_extension(string arg,string filename,string tag)
{
    int i             = 0;                  //standard variable used in loops
    int bad_extension = 0;                  //tells if the correct extension was provided
    int filename_size = filename.length();  //how long is the filename
    int tag_size      = tag.length();       //how long is the tag

    //make sure tag is not longer than filename
    if(tag_size > filename_size)
    {
        bad_extension = 1;
    }
    else //check for proper extension
    {
        for(i=filename_size-tag_size; i<filename_size; i++)
        {
            if(filename[i] != tag[i + tag_size - filename_size])
            {
                bad_extension = 1;
                break;
            }
        }
    }

    //report error and terminate program
    if(bad_extension == 1)
    {
        printf("The filename %s provided via the %s tag requires the %s extension. \n",filename.c_str(),arg.c_str(),tag.c_str());
        exit(EXIT_SUCCESS);
    }
}


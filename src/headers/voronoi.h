
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return a voronoi diagram                                                                   //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
Grid_i voronoi_diagram(Trajectory &traj,double aps,int num_g_x,int num_g_y,Param &param,double c_dist,double radius,int b_prot,int b_dynamic,int b_pbc)
{
    //when a dynamic box is used then the number of lattice points is chosen to best match the area of the box.
    //otherwise, the number of lattice points in each direction is specified by the user which might overestimate
    //the size of the box. The difference matters when the voronoi cell is used to estimate the area per lipid. 
    //in this case, the area will be overestimated for lipids near the boundaries if the lattice area overshoots 
    //the area of the box. in other instances, it is not important if the lattice dimensions perfectly match the box
    //and we generally try to make the grid a little bigger. 
    int    i         = 0;                      //standard variable used in loops
    int    j         = 0;                      //standard variable used in loops
    int    k         = 0;                      //standard variable used in loops
    int    l         = 0;                      //standard variable used in loops
    int    m         = 0;                      //standard variable used in loops
    int    n         = 0;                      //standard variable used in loops
    int    x         = 0;                      //used for pbc x loops
    int    y         = 0;                      //used for pbc y loops
    int x_min        = 0;                      //lower bound for pbc loop in x
    int x_max        = 0;                      //upper bound for pbc loop in x
    int y_min        = 0;                      //lower bound for pbc loop in y
    int y_max        = 0;                      //upper bound for pbc loop in y
    double cell_size = sqrt(aps);              //This disatance (nm) between lattice points

    //create a grid to hold the voronoi diagram 
    Grid_i voronoi;

    //check pbc conditions
    if(b_pbc == 1)
    {
        x_min = -1;
        x_max =  1;
        y_min = -1;
        y_max =  1;
    }

    //get grid dimensions if a dynamic box is used
    if(b_dynamic == 1)
    {
        //over write the num_g_x/y. this is okay since a copy is used in this function 
        num_g_x = (int)ceil(traj.box[XX][XX]/cell_size);
        num_g_y = (int)ceil(traj.box[YY][YY]/cell_size);
    }

    //set the grid dimensions
    voronoi.set_dim(aps,num_g_x,num_g_y);

    //printf("using dynamic grid num_g_x %d num_g_y %d x_min %d x_max %d y_mix %d y_max %d traj.box[XX][XX] %f traj.box[YY][YY] %f \n",num_g_x,num_g_y,x_min,x_max,y_min,y_max,traj.box[XX][XX],traj.box[YY][YY]);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                          //
    // select protein atoms in the same plane as target lipid atoms                                             //
    //                                                                                                          //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //create vector to hold protein atoms in same plane as lipids
    iv1d protein_atoms(0,0);

    if(b_prot == 1)
    {
        for(i=0; i<traj.target_leaflet.size(); i++) //loop over the target leaflet atoms
        {
            //get the first and last atom of the current lipid
            int min = traj.t_lip_start(i);
            int max = traj.t_lip_end(i);

            //jump to the next lipid
            i = traj.next_target_lipid(i);

            for(j=0; j<param.main_size_y(); j++) //loop over lipid types
            {
                if(strcmp(traj.res_name[min].c_str(), param.param_main_s[j][0].c_str() ) == 0) //lipid type is correct
                {
                    for(k=min; k<=max; k++) //loop over current residue atoms
                    {
                        for(l=0; l<param.sec_size_y(j); l++) //loop over target atoms
                        {
                            if(strcmp(traj.atom_name[k].c_str(), param.param_sec_s[j][l][0].c_str() ) == 0) //atom is a target atom
                            {
                                for(m=0; m<traj.prot.size(); m++) //loop over protein atoms
                                {
                                    double dx = traj.r[k][0] - traj.r[traj.prot[m]-1][0];
                                    double dy = traj.r[k][1] - traj.r[traj.prot[m]-1][1];
                                    double dz = traj.r[k][2] - traj.r[traj.prot[m]-1][2];

                                    double distance = sqrt(dx*dx + dy*dy + dz*dz);
                                    double dist_z   = sqrt(dz*dz);

                                    if(distance < c_dist && dist_z < 0.5*c_dist) //add the atom to voronoi diagram
                                    {
                                        int duplicate = 0;

                                        //check that atoms has not been counted already
                                        for(n=0; n<protein_atoms.size(); n++) //loop over protein atoms
                                        {
                                            if(protein_atoms[n] == k) //atom already added
                                            {
                                                duplicate = 1;
                                                break;
                                            }
                                        }

                                        //add atom to the list
                                        if(duplicate == 0)
                                        {
                                            protein_atoms.push_back(traj.prot[m]-1);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                          //
    // stamp target lipid atoms to grid, i.e., make a list of candidates                                        //
    //                                                                                                          //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    //create grid to hold stamped atom id list 
    iv3d candidates(num_g_y,iv2d(num_g_x,iv1d(0,0)));
    iv3d res_id(num_g_y,iv2d(num_g_x,iv1d(0,0)));

    //check lipid atoms for shortest distance
    for(i=0; i<traj.target_leaflet.size(); i++) //loop over the target leaflet atoms
    {
        //get the first and last atom of the current lipid
        int min = traj.t_lip_start(i);
        int max = traj.t_lip_end(i);

        //jump to the next lipid
        i = traj.next_target_lipid(i);

        for(j=0; j<param.main_size_y(); j++) //loop over lipid types
        {
            if(strcmp(traj.res_name[min].c_str(), param.param_main_s[j][0].c_str() ) == 0) //lipid type is correct
            {
                for(k=min; k<=max; k++) //loop over current residue atoms
                {
                    for(l=0; l<param.sec_size_y(j); l++) //loop over target atoms
                    {
                        if(strcmp(traj.atom_name[k].c_str(), param.param_sec_s[j][l][0].c_str() ) == 0) //atom is a target atom
                        {
                            for(x=x_min; x<=x_max; x++) //shift in x direction
                            {
                                for(y=y_min; y<=y_max; y++) //shift in y direction
                                {
                                    double hx = traj.r[k][0] + (double)x*traj.box[XX][XX];
                                    double hy = traj.r[k][1] + (double)y*traj.box[YY][YY];

                                    int lower_x = 0;
                                    int lower_y = 0;
                                    int upper_x = 0;
                                    int upper_y = 0;

                                    get_bounds(hx,hy,radius,voronoi.get_cell_size(),num_g_x,num_g_y,&lower_x,&lower_y,&upper_x,&upper_y);

                                    for(m=lower_x; m<upper_x; m++) //loop over the grid x-axis
                                    {
                                        for(n=lower_y; n<upper_y; n++) //loop over the grid y-axis
                                        {
                                            double dist_x = m*voronoi.get_cell_size() - hx;
                                            double dist_y = n*voronoi.get_cell_size() - hy;

                                            double dist = sqrt(dist_x*dist_x + dist_y*dist_y);

                                            if(dist <= radius)
                                            {
                                                candidates[n][m].push_back(k);
                                                res_id[n][m].push_back(traj.res_nr[k]);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                          //
    // add target protein atoms to list of candidates                                                           //
    //                                                                                                          //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    for(i=0; i<protein_atoms.size(); i++) //loop over target protein atoms
    {
        double hx = traj.r[protein_atoms[i]][0];
        double hy = traj.r[protein_atoms[i]][1];
  
        int lower_x = 0;
        int lower_y = 0;
        int upper_x = 0;
        int upper_y = 0;

        get_bounds(hx,hy,radius,voronoi.get_cell_size(),num_g_x,num_g_y,&lower_x,&lower_y,&upper_x,&upper_y);

        for(j=lower_x; j<upper_x; j++) //loop over the grid x-axis
        {
            for(k=lower_y; k<upper_y; k++) //loop over the grid y-axis
            {
                double dist_x = j*voronoi.get_cell_size() - hx;
                double dist_y = k*voronoi.get_cell_size() - hy;

                double dist = sqrt(dist_x*dist_x + dist_y*dist_y);

                if(dist <= radius)
                {
                    candidates[k][j].push_back(protein_atoms[i]);
                    res_id[k][j].push_back(-1);
                }
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                          //
    // construct a voronoi diagram                                                                              //
    //                                                                                                          //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    for(i=0; i<voronoi.num_x(); i++) //loop over x
    {
        for(j=0; j<voronoi.num_y(); j++) //loop over y
        {
            double min_dist = 999999.9;
            int    winner   = -1;
             
            if(candidates[j][i].size() > 0) //check candidates
            {
                for(k=0; k<candidates[j][i].size(); k++) //loop over candidates
                {
                    for(x=x_min; x<=x_max; x++) //shift in x direction
                    {
                        for(y=y_min; y<=y_max; y++) //shift in y direction
                        {
                            double dx = traj.r[candidates[j][i][k]][0] + (double)x*traj.box[XX][XX] - i*voronoi.get_cell_size();
                            double dy = traj.r[candidates[j][i][k]][1] + (double)y*traj.box[YY][YY] - j*voronoi.get_cell_size();

                            double distance = sqrt(dx*dx + dy*dy);
 
                            if(distance < min_dist)
                            {
                                min_dist = distance;
                                winner   = res_id[j][i][k];
                            }
                        }
                    }
                }
            }
            else //check all lipids and protein atoms 
            {
                //printf("Could not find a candidate. Computing all distances. x %d y %d \n",i,j);

                //check lipid atoms for shortest distance
                for(k=0; k<traj.target_leaflet.size(); k++) //loop over the target leaflet atoms
                {
                    //get the first and last atom of the current lipid
                    int min = traj.t_lip_start(k);
                    int max = traj.t_lip_end(k);

                    //jump to the next lipid
                    k = traj.next_target_lipid(k);

                    for(l=0; l<param.main_size_y(); l++) //loop over lipid types
                    {
                        if(strcmp(traj.res_name[min].c_str(), param.param_main_s[l][0].c_str() ) == 0) //lipid type is correct
                        {
                            for(m=min; m<=max; m++) //loop over current residue atoms
                            {
                                for(n=0; n<param.sec_size_y(l); n++) //loop over target atoms
                                {
                                    if(strcmp(traj.atom_name[m].c_str(), param.param_sec_s[l][n][0].c_str() ) == 0) //atom is a target atom
                                    {
                                        for(x=x_min; x<=x_max; x++) //shift in x direction
                                        {
                                            for(y=y_min; y<=y_max; y++) //shift in y direction
                                            {
                                                double dx = traj.r[m][0] + (double)x*traj.box[XX][XX] - i*voronoi.get_cell_size();
                                                double dy = traj.r[m][1] + (double)y*traj.box[YY][YY] - j*voronoi.get_cell_size();

                                                double distance = sqrt(dx*dx + dy*dy);

                                                if(distance < min_dist)
                                                {
                                                    min_dist = distance;
                                                    winner   = traj.res_nr[min];
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                //check protein atoms for shortest distance
                for(k=0; k<protein_atoms.size(); k++) //loop over protein atoms
                {
                    double dx = traj.r[protein_atoms[k]][0] - i*voronoi.get_cell_size();
                    double dy = traj.r[protein_atoms[k]][1] - j*voronoi.get_cell_size();

                    double distance = sqrt(dx*dx + dy*dy);

                    if(distance < min_dist)
                    {
                        min_dist = distance;
                        winner   = -1;
                    }
                }
            }

            //set the grid point
            voronoi.grid[j][i] = winner;
        }
    }

    return voronoi;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                          //
// This function return a voronoi diagram using the center of mass of the atoms selection                   //
//                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
Grid_i voronoi_diagram_com(Trajectory &traj,double aps,int num_g_x,int num_g_y,Param &param,double c_dist,double radius,int b_prot,int b_dynamic,int b_pbc)
{
    //when a dynamic box is used then the number of lattice points is chosen to best match the area of the box.
    //otherwise, the number of lattice points in each direction is specified by the user which might overestimate
    //the size of the box. The difference matters when the voronoi cell is used to estimate the area per lipid. 
    //in this case, the area will be overestimated for lipids near the boundaries if the lattice area overshoots 
    //the area of the box. in other instances, it is not important if the lattice dimensions perfectly match the box
    //and we generally try to make the grid a little bigger. 
    int    i          = 0;                      //standard variable used in loops
    int    j          = 0;                      //standard variable used in loops
    int    k          = 0;                      //standard variable used in loops
    int    l          = 0;                      //standard variable used in loops
    int    m          = 0;                      //standard variable used in loops
    int    n          = 0;                      //standard variable used in loops
    int    x          = 0;                      //used for pbc x loops
    int    y          = 0;                      //used for pbc y loops
    int x_min         = 0;                      //lower bound for pbc loop in x
    int x_max         = 0;                      //upper bound for pbc loop in x
    int y_min         = 0;                      //lower bound for pbc loop in y
    int y_max         = 0;                      //upper bound for pbc loop in y
    int found_centers = 0;                      //Have the centers of mass been computed yet?
    double cell_size  = sqrt(aps);              //This disatance (nm) between lattice points

    dv2d com_all(0,dv1d(2,0.0));                //Store the center of mass for each lipid

    //create a grid to hold the voronoi diagram 
    Grid_i voronoi;

    //check pbc conditions
    if(b_pbc == 1)
    {
        x_min = -1;
        x_max =  1;
        y_min = -1;
        y_max =  1;
    }

    //get grid dimensions if a dynamic box is used
    if(b_dynamic == 1)
    {
        //over write the num_g_x/y. this is okay since a copy is used in this function 
        num_g_x = (int)ceil(traj.box[XX][XX]/cell_size);
        num_g_y = (int)ceil(traj.box[YY][YY]/cell_size);
    }

    //set the grid dimensions
    voronoi.set_dim(aps,num_g_x,num_g_y);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                          //
    // stamp target lipid atoms to grid, i.e., make a list of candidates                                        //
    //                                                                                                          //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    //create grid to hold stamped atom id list 
    iv3d candidates(num_g_y,iv2d(num_g_x,iv1d(0,0)));
    dv3d com_x(num_g_y,dv2d(num_g_x,dv1d(0,0.0)));
    dv3d com_y(num_g_y,dv2d(num_g_x,dv1d(0,0.0)));

    //check lipid atoms for shortest distance
    for(i=0; i<traj.target_leaflet.size(); i++) //loop over the target leaflet atoms
    {
        //get the first and last atom of the current lipid
        int min = traj.t_lip_start(i);
        int max = traj.t_lip_end(i);

        //jump to the next lipid
        i = traj.next_target_lipid(i);

        for(j=0; j<param.main_size_y(); j++) //loop over lipid types
        {
            if(strcmp(traj.res_name[min].c_str(), param.param_main_s[j][0].c_str() ) == 0) //lipid type is correct
            {
                sv1d targets = param.get_column_sec_s(j,0);
                dv1d com     = traj.com(targets,min,max);

                com_all.push_back(com);

                for(x=x_min; x<=x_max; x++) //shift in x direction
                {
                    for(y=y_min; y<=y_max; y++) //shift in y direction
                    {
                        double hx = com[0] + (double)x*traj.box[XX][XX];
                        double hy = com[1] + (double)y*traj.box[YY][YY];

                        int lower_x = 0;
                        int lower_y = 0;
                        int upper_x = 0;
                        int upper_y = 0;

                        get_bounds(hx,hy,radius,voronoi.get_cell_size(),num_g_x,num_g_y,&lower_x,&lower_y,&upper_x,&upper_y);

                        for(m=lower_x; m<upper_x; m++) //loop over the grid x-axis
                        {
                            for(n=lower_y; n<upper_y; n++) //loop over the grid y-axis
                            {
                                double dist_x = m*voronoi.get_cell_size() - hx;
                                double dist_y = n*voronoi.get_cell_size() - hy;

                                double dist = sqrt(dist_x*dist_x + dist_y*dist_y);

                                if(dist <= radius)
                                {
                                    candidates[n][m].push_back(traj.res_nr[min]);
                                    com_x[n][m].push_back(com[0]);
                                    com_y[n][m].push_back(com[1]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                                          //
    // construct a voronoi diagram                                                                              //
    //                                                                                                          //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    for(i=0; i<voronoi.num_x(); i++) //loop over x
    {
        for(j=0; j<voronoi.num_y(); j++) //loop over y
        {
            double min_dist = 999999.9;
            int    winner   = -1;

            if(candidates[j][i].size() > 0) //check candidates
            {
                for(k=0; k<candidates[j][i].size(); k++) //loop over candidates
                {
                    for(x=x_min; x<=x_max; x++) //shift in x direction
                    {
                        for(y=y_min; y<=y_max; y++) //shift in y direction
                        {
                            double dx = com_x[j][i][k] + (double)x*traj.box[XX][XX] - i*voronoi.get_cell_size();
                            double dy = com_y[j][i][k] + (double)y*traj.box[YY][YY] - j*voronoi.get_cell_size();

                            double distance = sqrt(dx*dx + dy*dy);

                            if(distance < min_dist)
                            {
                                min_dist = distance;
                                winner   = candidates[j][i][k];
                            }
                        }
                    }
                }
            }
            else //check all lipids center of mass 
            {
                int counter = -1; //count lipids as they are encountered

                //check lipid atoms for shortest distance
                for(k=0; k<traj.target_leaflet.size(); k++) //loop over the target leaflet atoms
                {
                    //get the first and last atom of the current lipid
                    int min = traj.t_lip_start(k);
                    int max = traj.t_lip_end(k);

                    //jump to the next lipid
                    k = traj.next_target_lipid(k);

                    for(l=0; l<param.main_size_y(); l++) //loop over lipid types
                    {
                        if(strcmp(traj.res_name[min].c_str(), param.param_main_s[l][0].c_str() ) == 0) //lipid type is correct
                        {
                            counter++;

                            for(x=x_min; x<=x_max; x++) //shift in x direction
                            {
                                for(y=y_min; y<=y_max; y++) //shift in y direction
                                {
                                    double dx = com_all[counter][0] + (double)x*traj.box[XX][XX] - i*voronoi.get_cell_size();
                                    double dy = com_all[counter][1] + (double)y*traj.box[YY][YY] - j*voronoi.get_cell_size();

                                    double distance = sqrt(dx*dx + dy*dy);

                                    if(distance < min_dist)
                                    {
                                        min_dist = distance;
                                        winner   = traj.res_nr[min];
                                    }
                                }
                            }
                        }
                    }
                }
            }

            //set the grid point
            voronoi.grid[j][i] = winner;
        }
    }

    return voronoi;
}


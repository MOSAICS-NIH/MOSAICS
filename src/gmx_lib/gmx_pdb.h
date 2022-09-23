/*
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 */

enum
{
    epbcXYZ, epbcNONE, epbcXY, epbcSCREW, epbcNR
};

//! Set to true if warning has been printed
static gmx_bool bWarnedGuess = FALSE;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function prints an atom line to a pdb file and follows the approach taken in gromacs pdbio.cp        //
// function gmx_fprintf_pdb_atomline                                                                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void print_atomline_pdb(string record,int atom_seq_number,const char *atom_name,char alternate_location,
                        const char *res_name,char chain_id,int res_seq_number,char res_insertion_code,real x,real y,
                        real z,real occupancy,real b_factor,const char *element)
{
    char     tmp_atomname[6], tmp_resname[6];
    gmx_bool start_name_in_col13;

    //Format atom name
    if (atom_name != nullptr)
    {
        //If the atom name is an element name with two chars, it should start already in column 13.
        //Otherwise it should start in column 14, unless the name length is 4 chars.
        if ( (element != nullptr) && (strlen(element) >= 2) && (gmx_strncasecmp(atom_name, element, 2) == 0) )
        {
            start_name_in_col13 = TRUE;
        }
        else
        {
            start_name_in_col13 = (strlen(atom_name) >= 4);
        }
        snprintf(tmp_atomname, sizeof(tmp_atomname), start_name_in_col13 ? "" : " ");
        strncat(tmp_atomname, atom_name, 4);
        tmp_atomname[5] = '\0';
    }
    else
    {
        tmp_atomname[0] = '\0';
    }

    //Format residue name 
    strncpy(tmp_resname, (res_name != nullptr) ? res_name : "", 4);

    //Make sure the string is terminated if strlen was > 4 
    tmp_resname[4] = '\0';

    //String is properly terminated, so now we can use strcat. By adding a
    //space we can write it right-justified, and if the original name was
    //three characters or less there will be a space added on the right side.    
    strcat(tmp_resname, " ");

    //Truncate integers so they fit 
    atom_seq_number = atom_seq_number % 100000;
    res_seq_number  = res_seq_number % 10000;

    printf("%-6s%5d %-4.4s%c%4.4s%c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n",
           record.c_str(),
           atom_seq_number,
           tmp_atomname,
           alternate_location,
           tmp_resname,
           chain_id,
           res_seq_number,
           res_insertion_code,
           x, y, z,
           occupancy,
           b_factor,
           (element != nullptr) ? element : "");
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function reads the box line from a pdb file.                                                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
static void read_cryst1(char *line, int *ePBC, matrix box)
{
    #define SG_SIZE 11
    char   sa[12], sb[12], sc[12], sg[SG_SIZE+1], ident;
    double fa, fb, fc, alpha, beta, gamma, cosa, cosb, cosg, sing;
    int    syma, symb, symc;
    int    ePBC_file;

    enum {epbcXYZ, epbcNONE, epbcXY, epbcSCREW, epbcNR};
    #define DEG2RAD          (M_PI/180.0) 

    sscanf(line, "%*s%s%s%s%lf%lf%lf", sa, sb, sc, &alpha, &beta, &gamma);

    ePBC_file = -1;
    if (strlen(line) >= 55)
    {
        strncpy(sg, line+55, SG_SIZE);
        sg[SG_SIZE] = '\0';
        ident       = ' ';
        syma        = 0;
        symb        = 0;
        symc        = 0;
        sscanf(sg, "%c %d %d %d", &ident, &syma, &symb, &symc);
        if (ident == 'P' && syma ==  1 && symb <= 1 && symc <= 1)
        {
            fc        = strtod(sc, nullptr)*0.1;
            ePBC_file = (fc > 0 ? epbcXYZ : epbcXY);
        }
        if (ident == 'P' && syma == 21 && symb == 1 && symc == 1)
        {
            ePBC_file = epbcSCREW;
        }
    }
    if (ePBC)
    {
        *ePBC = ePBC_file;
    }

    if (box)
    {
        fa = strtod(sa, nullptr)*0.1;
        fb = strtod(sb, nullptr)*0.1;
        fc = strtod(sc, nullptr)*0.1;
        if (ePBC_file == epbcSCREW)
        {
            fa *= 0.5;
        }
        clear_mat(box);
        box[XX][XX] = fa;
        if ((alpha != 90.0) || (beta != 90.0) || (gamma != 90.0))
        {
            if (alpha != 90.0)
            {
                cosa = cos(alpha*DEG2RAD);
            }
            else
            {
                cosa = 0;
            }
            if (beta != 90.0)
            {
                cosb = cos(beta*DEG2RAD);
            }
            else
            {
                cosb = 0;
            }
            if (gamma != 90.0)
            {
                cosg = cos(gamma*DEG2RAD);
                sing = sin(gamma*DEG2RAD);
            }
            else
            {
                cosg = 0;
                sing = 1;
            }
            box[YY][XX] = fb*cosg;
            box[YY][YY] = fb*sing;
            box[ZZ][XX] = fc*cosb;
            box[ZZ][YY] = fc*(cosa - cosb*cosg)/sing;
            box[ZZ][ZZ] = sqrt(fc*fc
                                    - box[ZZ][XX]*box[ZZ][XX] - box[ZZ][YY]*box[ZZ][YY]);
        }
        else
        {
            box[YY][YY] = fb;
            box[ZZ][ZZ] = fc;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function is the same as printf_atomline_pdb but uses fprintf instead of printf                       //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void fprintf_atomline_pdb(FILE **out_file,string record,int atom_seq_number,const char *atom_name,char alternate_location,
                        const char *res_name,char chain_id,int res_seq_number,char res_insertion_code,real x,real y,
                        real z,real occupancy,real b_factor,const char *element)
{
    char     tmp_atomname[6], tmp_resname[6];
    gmx_bool start_name_in_col13;

    //Format atom name
    if (atom_name != nullptr)
    {
        //If the atom name is an element name with two chars, it should start already in column 13.
        //Otherwise it should start in column 14, unless the name length is 4 chars.
        if ( (element != nullptr) && (strlen(element) >= 2) && (gmx_strncasecmp(atom_name, element, 2) == 0) )
        {
            start_name_in_col13 = TRUE;
        }
        else
        {
            start_name_in_col13 = (strlen(atom_name) >= 4);
        }
        snprintf(tmp_atomname, sizeof(tmp_atomname), start_name_in_col13 ? "" : " ");
        strncat(tmp_atomname, atom_name, 4);
        tmp_atomname[5] = '\0';
    }
    else
    {
        tmp_atomname[0] = '\0';
    }

    //Format residue name
    strncpy(tmp_resname, (res_name != nullptr) ? res_name : "", 4);
    //Make sure the string is terminated if strlen was > 4
    tmp_resname[4] = '\0';
    //String is properly terminated, so now we can use strcat. By adding a
    //space we can write it right-justified, and if the original name was
    //three characters or less there will be a space added on the right side.

    strcat(tmp_resname, " ");

    //Truncate integers so they fit
    atom_seq_number = atom_seq_number % 100000;
    res_seq_number  = res_seq_number % 10000;

    fprintf(*out_file,"%-6s%5d %-4.4s%c%4.4s%c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n",
           record.c_str(),
           atom_seq_number,
           tmp_atomname,
           alternate_location,
           tmp_resname,
           chain_id,
           res_seq_number,
           res_insertion_code,
           x, y, z,
           occupancy,
           b_factor,
           (element != nullptr) ? element : "");
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function guesses the pbc based on box dimensions                                                     //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int guess_ePBC(const matrix box)
{
    int ePBC;

    if (box[XX][XX] > 0 && box[YY][YY] > 0 && box[ZZ][ZZ] > 0)
    {
        ePBC = epbcXYZ;
    }
    else if (box[XX][XX] > 0 && box[YY][YY] > 0 && box[ZZ][ZZ] == 0)
    {
        ePBC = epbcXY;
    }
    else if (box[XX][XX] == 0 && box[YY][YY] == 0 && box[ZZ][ZZ] == 0)
    {
        ePBC = epbcNONE;
    }
    else
    {
        if (!bWarnedGuess)
        {
            fprintf(stderr, "WARNING: Unsupported box diagonal %f %f %f, "
                    "will not use periodic boundary conditions\n\n",
                    box[XX][XX], box[YY][YY], box[ZZ][ZZ]);
            bWarnedGuess = TRUE;
        }
        ePBC = epbcNONE;
    }

    return ePBC;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function writes the system box to a pdb file                                                         //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void gmx_write_pdb_box(FILE *out, int ePBC, const matrix box)
{
    real alpha, beta, gamma;

    if (ePBC == -1)
    {
        ePBC = guess_ePBC(box);
    }

    if (ePBC == epbcNONE)
    {
        return;
    }

    if (norm2(box[YY])*norm2(box[ZZ]) != 0)
    {
        alpha = RAD2DEG*gmx_angle(box[YY], box[ZZ]);
    }
    else
    {
        alpha = 90;
    }
    if (norm2(box[XX])*norm2(box[ZZ]) != 0)
    {
        beta  = RAD2DEG*gmx_angle(box[XX], box[ZZ]);
    }
    else
    {
        beta  = 90;
    }
    if (norm2(box[XX])*norm2(box[YY]) != 0)
    {
        gamma = RAD2DEG*gmx_angle(box[XX], box[YY]);
    }
    else
    {
        gamma = 90;
    }
    fprintf(out, "REMARK    THIS IS A SIMULATION BOX\n");
    if (ePBC != epbcSCREW)
    {
        fprintf(out, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n",
                10*norm(box[XX]), 10*norm(box[YY]), 10*norm(box[ZZ]),
                alpha, beta, gamma, "P 1", 1);
    }
    else
    {
        /* Double the a-vector length and write the correct space group */
        fprintf(out, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n",
                20*norm(box[XX]), 10*norm(box[YY]), 10*norm(box[ZZ]),
                alpha, beta, gamma, "P 21 1 1", 1);

    }
}


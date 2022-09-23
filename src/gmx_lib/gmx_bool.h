/*
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

/*! \brief
 * Boolean type for use in \Gromacs C code.
 *
 * There is no standard size for 'bool' in C++, so when
 * we previously defined it to int for C code the data types
 * (and structs) would have different size depending on your compiler,
 * both at \Gromacs build time and when you use the library.
 * The only way around this is to NOT assume anything about the C++ type,
 * so we cannot use the name 'bool' in our C code anymore.
 */
typedef int gmx_bool;

#ifndef FALSE
/** False value for ::gmx_bool. */
#  define FALSE   0
#endif
#ifndef TRUE
/** True value for ::gmx_bool. */
#  define TRUE    1
#endif
/** Number of gmx_bool values. */
#define BOOL_NR 2

/*
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016, by the GROMACS development team, led by
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

#define XX      0 /* Defines for indexing in */
#define YY      1 /* vectors                 */
#define ZZ      2
#define DIM     3 /* Dimension of vectors    */

static inline void clear_mat(matrix a)
{
    const real nul = 0.0;

    a[XX][XX] = a[XX][YY] = a[XX][ZZ] = nul;
    a[YY][XX] = a[YY][YY] = a[YY][ZZ] = nul;
    a[ZZ][XX] = a[ZZ][YY] = a[ZZ][ZZ] = nul;
}

static inline real iprod(const rvec a, const rvec b)
{
    return (a[XX]*b[XX]+a[YY]*b[YY]+a[ZZ]*b[ZZ]);
}

static inline real norm2(const rvec a)
{
    return a[XX]*a[XX]+a[YY]*a[YY]+a[ZZ]*a[ZZ];
}

static inline void cprod(const rvec a, const rvec b, rvec c)
{
    c[XX] = a[YY]*b[ZZ]-a[ZZ]*b[YY];
    c[YY] = a[ZZ]*b[XX]-a[XX]*b[ZZ];
    c[ZZ] = a[XX]*b[YY]-a[YY]*b[XX];
}

/* WARNING:
 * As norm() uses sqrt() (which is slow) _only_ use it if you are sure you
 * don't need 1/norm(), otherwise use norm2()*invnorm(). */
static inline real norm(const rvec a)
{
    return sqrt(iprod(a, a));
}

/* This routine calculates the angle between a & b without any loss of accuracy close to 0/PI.
 * If you only need cos(theta), use the cos_angle() routines to save a few cycles.
 * This routine is faster than it might appear, since atan2 is accelerated on many CPUs (e.g. x86).
 */
static inline real gmx_angle(const rvec a, const rvec b)
{
    rvec w;
    real wlen, s;

    cprod(a, b, w);

    wlen  = norm(w);
    s     = iprod(a, b);

    return atan2(wlen, s);
}

static inline void clear_rvec(rvec a)
{
    /* The ibm compiler has problems with inlining this
     * when we use a const real variable
     */
    a[XX] = 0.0;
    a[YY] = 0.0;
    a[ZZ] = 0.0;
}

static inline void rvec_dec(rvec a, const rvec b)
{
    real x, y, z;

    x = a[XX]-b[XX];
    y = a[YY]-b[YY];
    z = a[ZZ]-b[ZZ];

    a[XX] = x;
    a[YY] = y;
    a[ZZ] = z;
}



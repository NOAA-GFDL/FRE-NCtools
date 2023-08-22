/***********************************************************************
 *                   GNU Lesser General Public License
 *
 * This file is part of the GFDL FRE NetCDF tools package (FRE-NCTools).
 *
 * FRE-NCtools is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * FRE-NCtools is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with FRE-NCTools.  If not, see
 * <http://www.gnu.org/licenses/>.
 **********************************************************************/

#ifndef NCTOOLS_CONSTANT_H
#define NCTOOLS_CONSTANT_H

//TODO: redefine constants without using macros.

#define RADIUS        (6371000.)
#define STRING        255

#ifndef M_PI
#define M_PI		(3.14159265358979323846)
#endif
#ifndef M_PI_2
#define M_PI_2      (1.57079632679489661923)
#endif


//TODO: change # defines to constants, e.g. 
//consexpression double  R2D {180/M_PI};
#define R2D (180/M_PI)
#define D2R (M_PI/180)
#define TPI (2.0*M_PI)

#endif


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
#ifndef GET_CONTACT_
#define GET_CONTACT_
/**********************************************************************
                       get_contact.h
    This header file is used to compute aligned-contact between tiles
**********************************************************************/

int get_align_contact(int tile1, int tile2, int nx1, int ny1, int nx2, int ny2, 
                      const double *x1, const double *y1, const double *x2, 
                      const double *y2, double periodx, double periody,
		      int *istart1, int *iend1, int *jstart1, int *jend1, 
                      int *istart2, int *iend2, int *jstart2, int *jend2);
int get_overlap_contact( int tile1, int tile2, int nx1, int ny1, int nx2, int ny2, 
                      const double *x1, const double *y1, const double *x2, 
                      const double *y2, int *istart1, int *iend1, int *jstart1, int *jend1, 
                      int *istart2, int *iend2, int *jstart2, int *jend2);
#endif

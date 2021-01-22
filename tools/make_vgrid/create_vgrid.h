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
 * License along with FRE-NCTools (LICENSE.md).  If not, see
 * <http://www.gnu.org/licenses/>.
 **********************************************************************/
/***********************************************************************
                       creater_vgrid.h
    This header file contains interface to create vertical grid on supergrid. 
    refinement =2 is assumed in this routine.
    contact: Zhi.Liang@noaa.gov
************************************************************************/

#ifndef CREATE_VGRID_H_
#define CREATE_VGRID_H_
void create_vgrid(int nbnds, double *bnds, int *nz, double *dbnds, double stretch, 
                  int use_legacy, double *zeta, int *np, const char *center);
#endif

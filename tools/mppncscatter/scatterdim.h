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
#ifndef SCATTERDIM_H
#define SCATTERDIM_H

#include "common.h"

typedef enum {NOSCATTER, SCATTERX, SCATTERY} scatter_t;

typedef struct ScatterDimStruct {
  int id;
  size_t len;
  char name[NC_MAX_NAME];
  
  scatter_t scatter_type;
  size_t* scatter_start;
  size_t* scatter_end;
  size_t* scatter_len;
  size_t scatter_ndiv; /* Number of divisions. */
} ScatterDim;

void ScatterDim_free(ScatterDim* p);

ScatterDim* ScatterDim_new(int id, size_t len, const char * name, scatter_t scatter_type, int ndiv);

#endif /* SCATTERDIM_H */
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
#include "scatterdim.h"

ScatterDim* ScatterDim_new(int id, size_t len, const char * name, scatter_t scatter_type, int ndiv) {
  ScatterDim* p = XMALLOC(ScatterDim,1);
  if (p==NULL) return p;
  
  p->id = id;
  p->len = len;
  strcpy(p->name, name);
  p->scatter_type = scatter_type;
  p->scatter_ndiv = ndiv;

  p->scatter_start = XMALLOC(size_t, ndiv);
  p->scatter_end = XMALLOC(size_t, ndiv);
  p->scatter_len = XMALLOC(size_t, ndiv);
  
  return p;
}
/*-------------------------------------------------------------------*/
void ScatterDim_free(ScatterDim* p) {
  if (p == NULL) return;
  
  if (p->scatter_start) {
    XFREE(p->scatter_start);
  }

  if (p->scatter_end) {
    XFREE(p->scatter_end);
  }

  if (p->scatter_len) {
    XFREE(p->scatter_len);
  }  
}

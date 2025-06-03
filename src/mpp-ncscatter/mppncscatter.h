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
#ifndef MPPNCSCATTER_H
#define MPPNCSCATTER_H

#include "strlist.h"
#include "opt.h"
#include "scatterdim.h"
#include "util.h"

/*-------------------------------------------------------------------*/
#define handle_error(status) {                      \
	if (status != NC_NOERR) {                         \
		fprintf(stderr, "%s\n", nc_strerror(status));   \
		exit(-1);                                       \
	}                                                 \
}

void def_dim(int nc, int *ncids, int ndims, ScatterDim *scatterdims[], mnsopts *opts);

void def_var(int nc, int *ncids, int nvars, int ndims, ScatterDim *scatterdims[], mnsopts *opts);

void free_scatter_dims(ScatterDim* dims[], int ndims);

int get_num_files(mnsopts* opts);

void get_num_divs(mnsopts* opts, int* nx, int* ny);

void get_scatter_dims(int nc, int ndims, int nvars, ScatterDim* scatterdims[], mnsopts *opt);

void get_scatter_extents(ScatterDim* scatterdims[], int ndims);

void get_scatter_extents_iolayout(ScatterDim* scatterdims[], int ndims, int nxio, int nyio);

void hyperslabcopy(nc_type type, size_t *dimlen, int *dimids, size_t *start, size_t *count, int ndim, char *ti, short *si, int *ii, float *fi, double *di, char *t, short *s, int *i, float *f, double *d);

void mpp_compute_extent(size_t isg, size_t ieg, size_t ndivs, size_t* start, size_t* end);

int mppncscatter(mnsopts* popts);

void print_scatter_dims(ScatterDim* scatterdims[], int ndims);

void printsizetarray(size_t *a, int n);

void put_var(int nc, int *ncids, int ndims, int nvars, ScatterDim *scatterdims[], mnsopts *opt);

void scatter_dims(int nc, int ndims, int nvars, ScatterDim* scatterdims[], mnsopts* opt);

/*-------------------------------------------------------------------*/

#endif /* MPPNCSCATTER_H */
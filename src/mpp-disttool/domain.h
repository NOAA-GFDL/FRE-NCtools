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
#ifndef MPPDISTTOOL_DOMAIN_H
#define MPPDISTTOOL_DOMAIN_H

#include <stddef.h>
#include <netcdf.h>
#include "strlist.h"

typedef enum { NOSCATTER = 0, SCATTERX, SCATTERY } scatter_t;

typedef struct {
    int      id;
    size_t   len;
    char     name[NC_MAX_NAME + 1];
    scatter_t scatter_type;
    size_t  *scatter_start;
    size_t  *scatter_end;
    size_t  *scatter_len;
    size_t   scatter_ndiv;
} ScatterDim;

/* Options for the MPP scatter operation. */
typedef struct {
    char    dryrun;
    char   *filein;
    int     nx;
    int     ny;
    int     nxio;
    int     nyio;
    char   *prefix;
    int     start;
    int     verbose;
    int     width;
    char  **xdims;
    int     xdims_len;
    char  **ydims;
    int     ydims_len;
} scatter_opts_t;

/* ScatterDim lifecycle. */
ScatterDim *ScatterDim_new(int id, size_t len, const char *name,
                           scatter_t scatter_type, int ndiv);
void        ScatterDim_free(ScatterDim *p);
void        scatter_dims_free(ScatterDim *dims[], int ndims);

/* Symmetric MPP tiling (from mppncscatter).
   Fills start[0..ndivs-1] and end[0..ndivs-1] with 0-based indices. */
void mpp_compute_extent(size_t isg, size_t ieg, size_t ndivs,
                        size_t *start, size_t *end);

/* Copy a hyperslab of multi-dimensional data (any nc_type). */
void hyperslabcopy(nc_type type, size_t *dimlen, int *dimids,
                   size_t *start, size_t *count, int ndim,
                   char *ti, short *si, int *ii, float *fi, double *di,
                   char *t,  short *s,  int *i,  float *f,  double *d);

/* Populate scatterdims[] from the open NetCDF file.
   Also applies io_layout compression if nxio/nyio are set. */
void scatter_dims(int nc, int ndims, int nvars,
                  ScatterDim *scatterdims[], scatter_opts_t *opt);

/* Output file count from scatter_opts_t. */
int scatter_get_num_files(scatter_opts_t *opts);

/* Define dimensions and variables in output file set. */
void scatter_def_dim(int nc, int *ncids, int ndims,
                     ScatterDim *scatterdims[], scatter_opts_t *opts);
void scatter_def_var(int nc, int *ncids, int nvars, int ndims,
                     ScatterDim *scatterdims[], scatter_opts_t *opts);

/* Write variable data into scattered output files. */
void scatter_put_var(int nc, int *ncids, int ndims, int nvars,
                     ScatterDim *scatterdims[], scatter_opts_t *opt);

#endif /* MPPDISTTOOL_DOMAIN_H */

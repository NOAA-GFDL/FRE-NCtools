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
#ifndef MPPDISTTOOL_COMPRESS_H
#define MPPDISTTOOL_COMPRESS_H

#include <stddef.h>
#include <stdbool.h>
#include <netcdf.h>

/* Mirrors the Fortran diminfo_type from nfu_compress.F90.
   For each dimension of a compressed variable:
     idx    – mapping from compressed position to flat uncompressed index
              (NULL for non-compressed dims; values are 0-based; -1 = invalid)
     length – size of this dimension in the compressed (input) array
     stride – stride of this dimension in the uncompressed output array */
typedef struct {
    int    *idx;
    int     length;
    int     stride;
} diminfo_t;

/* Check whether a dimension is CF "compressed by gathering".
   Returns NC_NOERR if the dimension variable has a "compress" attribute.
   On success, fills ndims_out, dimids_out[0..ndims-1], dimlens_out[0..ndims-1]
   with the expanded (uncompressed) dimension info.
   Pass NULL for outputs you don't need. */
int compress_inq_dim(int ncid, int dimid,
                     int    *ndims_out,
                     int    *dimids_out,    /* [NC_MAX_DIMS] */
                     size_t *dimlens_out);  /* [NC_MAX_DIMS] */

/* Initialise a diminfo_t array for a compressed variable.
   Fills diminfo[0..ndims-1] based on the variable's dimension list.
   Allocates diminfo[i].idx for each compressed dimension.
   Call compress_free_diminfo() when done. */
int compress_build_diminfo(int ncid, int varid,
                           diminfo_t *diminfo, int *ndims_out,
                           size_t *varsize_out);

/* Free allocated idx arrays from compress_build_diminfo(). */
void compress_free_diminfo(diminfo_t *diminfo, int ndims);

/* Read a compressed variable into a pre-allocated uncompressed array.
   'data' must hold at least uncompressed_size doubles (initialised by caller).
   'mask' (optional) is set to true at positions that were filled.
   Pass NULL for mask if not needed. */
int compress_get_var_double(int ncid, int varid,
                            double *data, bool *mask,
                            size_t uncompressed_size);

/* Variant that reads only a sub-range along the first dimension.
   start[0..ndims-1] and count[0..ndims-1] follow nc_get_vara conventions
   but must have count[i]==1 for i>0.  Pass NULL to read the whole variable. */
int compress_get_var_double_range(int ncid, int varid,
                                  const size_t *start, const size_t *count,
                                  double *data, bool *mask,
                                  size_t uncompressed_size);

/* Write an uncompressed array into a compressed variable. */
int compress_put_var_double(int ncid, int varid,
                            const double *src, size_t uncompressed_size);

/* Partial write: nc_put_vara_double for the compressed variable. */
int compress_put_vara_double(int ncid, int varid,
                             const size_t *start, const size_t *count,
                             const double *buf);

/* ---- Merge-sort helpers (ported from combine-ncc.F90) ---- */

/* Rank array x in ascending order: idx[0..n-1] receives indices of x
   in non-decreasing order of x values (0-based). */
void rank_ascending(const int *x, int *idx, int n);

/* Check if any variable in the file is compressed by gathering.
   Returns true if at least one dimension variable has a "compress" attribute. */
bool file_is_compressed(int ncid);

#endif /* MPPDISTTOOL_COMPRESS_H */

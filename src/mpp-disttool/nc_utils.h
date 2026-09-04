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
#ifndef MPPDISTTOOL_NC_UTILS_H
#define MPPDISTTOOL_NC_UTILS_H

#include <stddef.h>
#include <stdbool.h>
#include <netcdf.h>

/* Clone a dimension from incid to oncid (creates if not present).
   If dimension already exists in oncid, verifies size matches. */
int ncu_clone_dim(int incid, int dimid, int oncid);

/* Clone a variable definition (and all its attributes) from incid to oncid.
   Dimensions in the variable are looked up by name in oncid. */
int ncu_clone_var(int incid, int varid, int oncid);

/* Copy all global attributes from incid to oncid. */
int ncu_copy_global_atts(int incid, int oncid);

/* Get integer variable data by name. 'data' must be pre-allocated. */
int ncu_get_var_int(int ncid, const char *name, int *data);

/* Query variable by id: optionally return name, xtype, ndims, dimids,
   dimlens (per-dim lengths), natts, is_dim flag, has_records flag,
   varsize (total elements), recsize (elements per record), nrec.
   Pass NULL for outputs you don't need. */
int ncu_inq_var(int ncid, int varid,
                char *name_out,        /* NC_MAX_NAME+1 or NULL */
                nc_type *xtype_out,
                int *ndims_out,
                int *dimids_out,       /* [NC_MAX_VAR_DIMS] or NULL */
                size_t *dimlens_out,   /* [NC_MAX_VAR_DIMS] or NULL */
                int *natts_out,
                bool *is_dim_out,
                bool *has_records_out,
                size_t *varsize_out,
                size_t *recsize_out,
                int *nrec_out);

/* Query dimension by id. Pass NULL for outputs not needed. */
int ncu_inq_dim(int ncid, int dimid,
                char *name_out,      /* NC_MAX_NAME+1 or NULL */
                size_t *len_out,
                bool *is_unlim_out);

/* Return the nc_create flags (format) appropriate for creating an output
   file that matches the format of incid.  If incid is classic, returns
   NC_64BIT_OFFSET (preserving mppnccombine behaviour). */
int ncu_format_cmode(int incid);

/* Bytes per netCDF element type. */
size_t ncu_type_size(nc_type xtype);

/* Read a single record of a variable (record dimension = rec, 0-based).
   'buf' must be pre-allocated to recsize elements of the variable type.
   Works for any nc_type. */
int ncu_get_rec(int ncid, int varid, int rec, void *buf);

/* Write a single record (rec, 0-based) to a variable. */
int ncu_put_rec(int ncid, int varid, int rec, const void *buf);

/* Generic vara read/write for any nc_type. */
int ncu_get_vara(int ncid, int varid, const size_t *start,
                 const size_t *count, void *buf);
int ncu_put_vara(int ncid, int varid, const size_t *start,
                 const size_t *count, const void *buf);

#endif /* MPPDISTTOOL_NC_UTILS_H */

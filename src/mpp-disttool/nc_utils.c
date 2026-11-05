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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include "nc_utils.h"

#define NC_CHECK(call) do { \
    int _err = (call); \
    if (_err != NC_NOERR) return _err; \
} while (0)

size_t
ncu_type_size(nc_type xtype)
{
    switch (xtype) {
    case NC_BYTE: case NC_CHAR: return 1;
    case NC_SHORT:              return 2;
    case NC_INT:                return 4;
    case NC_FLOAT:              return 4;
    case NC_DOUBLE:             return 8;
    default:                    return 0;
    }
}

int
ncu_clone_dim(int incid, int dimid, int oncid)
{
    char name[NC_MAX_NAME + 1];
    size_t len;
    int unlimdimid, newdimid;

    NC_CHECK(nc_inq_dim(incid, dimid, name, &len));
    NC_CHECK(nc_inq_unlimdim(incid, &unlimdimid));
    if (dimid == unlimdimid)
        len = NC_UNLIMITED;

    /* If already present in oncid, accept as-is. */
    if (nc_inq_dimid(oncid, name, &newdimid) == NC_NOERR)
        return NC_NOERR;

    NC_CHECK(nc_def_dim(oncid, name, len, &newdimid));
    return NC_NOERR;
}

int
ncu_clone_var(int incid, int varid, int oncid)
{
    char name[NC_MAX_NAME + 1];
    char dname[NC_MAX_NAME + 1];
    char attname[NC_MAX_NAME + 1];
    nc_type xtype;
    int ndims, dimids[NC_MAX_VAR_DIMS], natts, ovarid;

    NC_CHECK(nc_inq_var(incid, varid, name, &xtype, &ndims, dimids, &natts));

    /* Remap dimension IDs to the output file's dimension IDs. */
    for (int i = 0; i < ndims; i++) {
        NC_CHECK(nc_inq_dimname(incid, dimids[i], dname));
        NC_CHECK(nc_inq_dimid(oncid, dname, &dimids[i]));
    }

    NC_CHECK(nc_def_var(oncid, name, xtype, ndims, dimids, &ovarid));

    for (int i = 0; i < natts; i++) {
        NC_CHECK(nc_inq_attname(incid, varid, i, attname));
        NC_CHECK(nc_copy_att(incid, varid, attname, oncid, ovarid));
    }
    return NC_NOERR;
}

int
ncu_copy_global_atts(int incid, int oncid)
{
    int ngatts;
    char attname[NC_MAX_NAME + 1];

    NC_CHECK(nc_inq_natts(incid, &ngatts));
    for (int i = 0; i < ngatts; i++) {
        NC_CHECK(nc_inq_attname(incid, NC_GLOBAL, i, attname));
        NC_CHECK(nc_copy_att(incid, NC_GLOBAL, attname, oncid, NC_GLOBAL));
    }
    return NC_NOERR;
}

int
ncu_get_var_int(int ncid, const char *name, int *data)
{
    int varid;
    NC_CHECK(nc_inq_varid(ncid, name, &varid));
    return nc_get_var_int(ncid, varid, data);
}

int
ncu_inq_var(int ncid, int varid,
            char *name_out,
            nc_type *xtype_out,
            int *ndims_out,
            int *dimids_out,
            size_t *dimlens_out,
            int *natts_out,
            bool *is_dim_out,
            bool *has_records_out,
            size_t *varsize_out,
            size_t *recsize_out,
            int *nrec_out)
{
    char name[NC_MAX_NAME + 1];
    nc_type xtype;
    int ndims, dimids[NC_MAX_VAR_DIMS], natts;
    int unlimdim;

    NC_CHECK(nc_inq_var(ncid, varid, name, &xtype, &ndims, dimids, &natts));
    NC_CHECK(nc_inq_unlimdim(ncid, &unlimdim));

    if (name_out)  strcpy(name_out, name);
    if (xtype_out) *xtype_out = xtype;
    if (ndims_out) *ndims_out = ndims;
    if (natts_out) *natts_out = natts;
    if (dimids_out)
        memcpy(dimids_out, dimids, ndims * sizeof(int));

    if (is_dim_out) {
        int dummy;
        *is_dim_out = (nc_inq_dimid(ncid, name, &dummy) == NC_NOERR);
    }
    if (has_records_out) {
        *has_records_out = false;
        for (int i = 0; i < ndims; i++)
            if (dimids[i] == unlimdim) { *has_records_out = true; break; }
    }

    if (dimlens_out || varsize_out || recsize_out || nrec_out) {
        size_t vsize = 1, rsize = 1;
        for (int i = 0; i < ndims; i++) {
            size_t dlen;
            NC_CHECK(nc_inq_dimlen(ncid, dimids[i], &dlen));
            if (dimlens_out) dimlens_out[i] = dlen;
            vsize *= dlen;
            if (dimids[i] != unlimdim) rsize *= dlen;
        }
        if (varsize_out) *varsize_out = vsize;
        if (recsize_out) *recsize_out = rsize;
        if (nrec_out) {
            *nrec_out = 1;
            if (unlimdim != -1) {
                bool has_rec = false;
                for (int i = 0; i < ndims; i++)
                    if (dimids[i] == unlimdim) { has_rec = true; break; }
                if (has_rec) {
                    size_t nrec_sz;
                    NC_CHECK(nc_inq_dimlen(ncid, unlimdim, &nrec_sz));
                    *nrec_out = (int)nrec_sz;
                }
            }
        }
    }
    return NC_NOERR;
}

int
ncu_inq_dim(int ncid, int dimid,
            char *name_out,
            size_t *len_out,
            bool *is_unlim_out)
{
    char name[NC_MAX_NAME + 1];
    size_t len;
    int unlimdim;

    NC_CHECK(nc_inq_dim(ncid, dimid, name, &len));
    if (name_out) strcpy(name_out, name);
    if (len_out)  *len_out = len;
    if (is_unlim_out) {
        NC_CHECK(nc_inq_unlimdim(ncid, &unlimdim));
        *is_unlim_out = (dimid == unlimdim);
    }
    return NC_NOERR;
}

int
ncu_format_cmode(int incid)
{
    int fmt;
    if (nc_inq_format(incid, &fmt) != NC_NOERR)
        return NC_64BIT_OFFSET;
    switch (fmt) {
    case NC_FORMAT_NETCDF4:
        return NC_NETCDF4;
    case NC_FORMAT_NETCDF4_CLASSIC:
        return NC_NETCDF4 | NC_CLASSIC_MODEL;
    case NC_FORMAT_64BIT_DATA:
        return NC_64BIT_DATA;
    /* classic and 64-bit offset both produce 64-bit offset output,
       matching mppnccombine behaviour */
    default:
        return NC_64BIT_OFFSET;
    }
}

int
ncu_get_vara(int ncid, int varid, const size_t *start,
             const size_t *count, void *buf)
{
    nc_type xtype;
    NC_CHECK(nc_inq_vartype(ncid, varid, &xtype));
    switch (xtype) {
    case NC_BYTE: case NC_CHAR:
        return nc_get_vara_text(ncid, varid, start, count, (char *)buf);
    case NC_SHORT:
        return nc_get_vara_short(ncid, varid, start, count, (short *)buf);
    case NC_INT:
        return nc_get_vara_int(ncid, varid, start, count, (int *)buf);
    case NC_FLOAT:
        return nc_get_vara_float(ncid, varid, start, count, (float *)buf);
    case NC_DOUBLE:
        return nc_get_vara_double(ncid, varid, start, count, (double *)buf);
    default:
        return NC_EBADTYPE;
    }
}

int
ncu_put_vara(int ncid, int varid, const size_t *start,
             const size_t *count, const void *buf)
{
    nc_type xtype;
    NC_CHECK(nc_inq_vartype(ncid, varid, &xtype));
    switch (xtype) {
    case NC_BYTE: case NC_CHAR:
        return nc_put_vara_text(ncid, varid, start, count, (const char *)buf);
    case NC_SHORT:
        return nc_put_vara_short(ncid, varid, start, count, (const short *)buf);
    case NC_INT:
        return nc_put_vara_int(ncid, varid, start, count, (const int *)buf);
    case NC_FLOAT:
        return nc_put_vara_float(ncid, varid, start, count, (const float *)buf);
    case NC_DOUBLE:
        return nc_put_vara_double(ncid, varid, start, count, (const double *)buf);
    default:
        return NC_EBADTYPE;
    }
}

int
ncu_get_rec(int ncid, int varid, int rec, void *buf)
{
    int ndims, dimids[NC_MAX_VAR_DIMS], unlimdim;
    size_t start[NC_MAX_VAR_DIMS], count[NC_MAX_VAR_DIMS];

    NC_CHECK(nc_inq_unlimdim(ncid, &unlimdim));
    NC_CHECK(nc_inq_varndims(ncid, varid, &ndims));
    NC_CHECK(nc_inq_vardimid(ncid, varid, dimids));

    for (int i = 0; i < ndims; i++) {
        if (dimids[i] == unlimdim) {
            start[i] = (size_t)rec;
            count[i] = 1;
        } else {
            start[i] = 0;
            NC_CHECK(nc_inq_dimlen(ncid, dimids[i], &count[i]));
        }
    }
    return ncu_get_vara(ncid, varid, start, count, buf);
}

int
ncu_put_rec(int ncid, int varid, int rec, const void *buf)
{
    int ndims, dimids[NC_MAX_VAR_DIMS], unlimdim;
    size_t start[NC_MAX_VAR_DIMS], count[NC_MAX_VAR_DIMS];

    NC_CHECK(nc_inq_unlimdim(ncid, &unlimdim));
    NC_CHECK(nc_inq_varndims(ncid, varid, &ndims));
    NC_CHECK(nc_inq_vardimid(ncid, varid, dimids));

    for (int i = 0; i < ndims; i++) {
        if (dimids[i] == unlimdim) {
            start[i] = (size_t)rec;
            count[i] = 1;
        } else {
            start[i] = 0;
            NC_CHECK(nc_inq_dimlen(ncid, dimids[i], &count[i]));
        }
    }
    return ncu_put_vara(ncid, varid, start, count, buf);
}

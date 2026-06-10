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
/*
 * C port of lib/libnfu/nfu_compress.F90 and the rank_ascending /
 * mergerank / merge routines from src/land_utils/combine-ncc.F90.
 */
#define _POSIX_C_SOURCE 200809L
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <netcdf.h>
#include "compress.h"
#include "nc_utils.h"
#include "xmalloc.h"

#define NC_CHECK(call) do { \
    int _err = (call); \
    if (_err != NC_NOERR) return _err; \
} while (0)

/* ------------------------------------------------------------------ */
/* compress_inq_dim – check if a dimension is compressed by gathering  */
/* ------------------------------------------------------------------ */
int
compress_inq_dim(int ncid, int dimid,
                 int    *ndims_out,
                 int    *dimids_out,
                 size_t *dimlens_out)
{
    char dimname[NC_MAX_NAME + 1];
    char compress_att[2048];
    int varid;
    size_t att_len;
    nc_type att_type;

    NC_CHECK(nc_inq_dimname(ncid, dimid, dimname));

    /* Does the dimension variable exist and have a "compress" attribute? */
    if (nc_inq_varid(ncid, dimname, &varid) != NC_NOERR)
        return NC_ENOTATT;
    if (nc_inq_att(ncid, varid, "compress", &att_type, &att_len) != NC_NOERR)
        return NC_ENOTATT;
    if (att_type != NC_CHAR || att_len == 0)
        return NC_ENOTATT;
    if (att_len >= sizeof(compress_att))
        att_len = sizeof(compress_att) - 1;
    NC_CHECK(nc_get_att_text(ncid, varid, "compress", compress_att));
    compress_att[att_len] = '\0';

    /* Parse space-separated dimension names in compress_att.
       The Fortran code scans from the end, producing them in reverse order.
       We do forward scanning, which gives the same set. */
    int n = 0;
    char *tok = strtok(compress_att, " \t");
    while (tok) {
        if (ndims_out || dimids_out || dimlens_out) {
            if (dimids_out || dimlens_out) {
                int did;
                NC_CHECK(nc_inq_dimid(ncid, tok, &did));
                if (dimids_out)  dimids_out[n]  = did;
                if (dimlens_out) {
                    size_t dlen;
                    NC_CHECK(nc_inq_dimlen(ncid, did, &dlen));
                    dimlens_out[n] = dlen;
                }
            }
            n++;
        }
        tok = strtok(NULL, " \t");
    }
    if (ndims_out) *ndims_out = n;
    return NC_NOERR;
}

/* ------------------------------------------------------------------ */
/* compress_build_diminfo                                               */
/* ------------------------------------------------------------------ */
int
compress_build_diminfo(int ncid, int varid,
                       diminfo_t *diminfo, int *ndims_out,
                       size_t *varsize_out)
{
    int ndims, dimids[NC_MAX_VAR_DIMS];
    NC_CHECK(nc_inq_varndims(ncid, varid, &ndims));
    NC_CHECK(nc_inq_vardimid(ncid, varid, dimids));

    if (ndims_out)  *ndims_out  = ndims;

    int stride = 1;
    for (int i = 0; i < ndims; i++) {
        size_t dimlen;
        int cdims, cdimids[NC_MAX_DIMS];
        size_t cdimlens[NC_MAX_DIMS];

        NC_CHECK(nc_inq_dimlen(ncid, dimids[i], &dimlen));
        diminfo[i].length = (int)dimlen;
        diminfo[i].idx    = NULL;

        int length_for_stride;
        if (compress_inq_dim(ncid, dimids[i],
                              &cdims, cdimids, cdimlens) == NC_NOERR) {
            /* Compressed dimension: load index array. */
            diminfo[i].idx = XMALLOC(int, dimlen);
            char cname[NC_MAX_NAME + 1];
            NC_CHECK(nc_inq_dimname(ncid, dimids[i], cname));
            NC_CHECK(ncu_get_var_int(ncid, cname, diminfo[i].idx));
            /* Stride contribution = product of uncompressed dim sizes. */
            length_for_stride = 1;
            for (int k = 0; k < cdims; k++)
                length_for_stride *= (int)cdimlens[k];
        } else {
            length_for_stride = diminfo[i].length;
        }
        diminfo[i].stride = stride;
        stride *= length_for_stride;
    }
    if (varsize_out) {
        /* Total compressed size = product of compressed lengths. */
        size_t vsz = 1;
        for (int i = 0; i < ndims; i++)
            vsz *= (size_t)diminfo[i].length;
        *varsize_out = vsz;
    }
    return NC_NOERR;
}

void
compress_free_diminfo(diminfo_t *diminfo, int ndims)
{
    for (int i = 0; i < ndims; i++) {
        if (diminfo[i].idx) {
            free(diminfo[i].idx);
            diminfo[i].idx = NULL;
        }
    }
}

/* ------------------------------------------------------------------ */
/* Internal: scatter compressed buffer into uncompressed output array  */
/* ------------------------------------------------------------------ */
static void
scatter_to_output(const double *buffer, size_t bufsize,
                  double *data, bool *mask,
                  const diminfo_t *diminfo, int ndims)
{
    int local_idx[NC_MAX_VAR_DIMS];
    memset(local_idx, 0, ndims * sizeof(int));

    for (size_t i = 0; i < bufsize; i++) {
        /* Compute 0-based flat output index. */
        int ii = 0;
        bool valid = true;
        for (int n = 0; n < ndims; n++) {
            if (diminfo[n].idx != NULL) {
                int v = diminfo[n].idx[local_idx[n]];
                if (v < 0) { valid = false; break; }
                ii += v * diminfo[n].stride;
            } else {
                ii += local_idx[n] * diminfo[n].stride;
            }
        }
        if (valid) {
            data[ii] = buffer[i];
            if (mask) mask[ii] = true;
        }
        /* Increment multi-dimensional counter. */
        for (int n = 0; n < ndims; n++) {
            local_idx[n]++;
            if (local_idx[n] < diminfo[n].length) break;
            local_idx[n] = 0;
        }
    }
}

/* ------------------------------------------------------------------ */
/* compress_get_var_double                                              */
/* ------------------------------------------------------------------ */
int
compress_get_var_double(int ncid, int varid,
                        double *data, bool *mask,
                        size_t uncompressed_size)
{
    return compress_get_var_double_range(ncid, varid, NULL, NULL,
                                         data, mask, uncompressed_size);
}

int
compress_get_var_double_range(int ncid, int varid,
                               const size_t *start, const size_t *count,
                               double *data, bool *mask,
                               size_t uncompressed_size)
{
    diminfo_t diminfo[NC_MAX_VAR_DIMS];
    int ndims;
    size_t varsize;

    NC_CHECK(compress_build_diminfo(ncid, varid, diminfo, &ndims, &varsize));

    size_t bufsize = (start && count) ? count[0] : varsize;
    double *buffer = XMALLOC(double, bufsize);

    int err;
    if (start && count) {
        err = nc_get_vara_double(ncid, varid, start, count, buffer);
        /* For sub-range reads, treat as 1D. */
        diminfo_t di1[1];
        di1[0].idx    = diminfo[0].idx;
        di1[0].length = (int)count[0];
        di1[0].stride = diminfo[0].stride;
        if (err == NC_NOERR)
            scatter_to_output(buffer, bufsize, data, mask, di1, 1);
        /* Don't free di1[0].idx since it's borrowed from diminfo[0]. */
    } else {
        err = nc_get_var_double(ncid, varid, buffer);
        if (err == NC_NOERR)
            scatter_to_output(buffer, bufsize, data, mask, diminfo, ndims);
    }

    free(buffer);
    compress_free_diminfo(diminfo, ndims);
    return err;
}

/* ------------------------------------------------------------------ */
/* compress_put_var_double                                              */
/* ------------------------------------------------------------------ */
int
compress_put_var_double(int ncid, int varid,
                        const double *src, size_t uncompressed_size)
{
    diminfo_t diminfo[NC_MAX_VAR_DIMS];
    int ndims;
    size_t varsize;

    NC_CHECK(compress_build_diminfo(ncid, varid, diminfo, &ndims, &varsize));

    double *buffer = XMALLOC(double, varsize);
    int local_idx[NC_MAX_VAR_DIMS];
    memset(local_idx, 0, ndims * sizeof(int));

    for (size_t i = 0; i < varsize; i++) {
        int ii = 0;
        for (int n = 0; n < ndims; n++) {
            if (diminfo[n].idx != NULL)
                ii += diminfo[n].idx[local_idx[n]] * diminfo[n].stride;
            else
                ii += local_idx[n] * diminfo[n].stride;
        }
        buffer[i] = src[ii];
        for (int n = 0; n < ndims; n++) {
            local_idx[n]++;
            if (local_idx[n] < diminfo[n].length) break;
            local_idx[n] = 0;
        }
    }

    int err = nc_put_var_double(ncid, varid, buffer);
    free(buffer);
    compress_free_diminfo(diminfo, ndims);
    return err;
}

int
compress_put_vara_double(int ncid, int varid,
                         const size_t *start, const size_t *count,
                         const double *buf)
{
    return nc_put_vara_double(ncid, varid, start, count, buf);
}

/* ------------------------------------------------------------------ */
/* Merge sort – ported from combine-ncc.F90:358-430                    */
/* ------------------------------------------------------------------ */
static void
merge_arrays(const int *x,
             const int *a, int na,
             const int *b, int nb,
             int *c, int nc_size)
{
    int i = 0, j = 0, k = 0;
    while (i < na && j < nb) {
        if (x[a[i]] <= x[b[j]])
            c[k++] = a[i++];
        else
            c[k++] = b[j++];
    }
    while (i < na)
        c[k++] = a[i++];
    (void)nc_size;
}

static void
mergerank(const int *x, int *a, int n, int *t)
{
    if (n < 2) return;
    if (n == 2) {
        if (x[a[0]] > x[a[1]]) { int v = a[0]; a[0] = a[1]; a[1] = v; }
        return;
    }
    int na = (n + 1) / 2;
    int nb = n - na;
    mergerank(x, a,      na, t);
    mergerank(x, a + na, nb, t);
    if (x[a[na - 1]] > x[a[na]]) {
        memcpy(t, a, na * sizeof(int));
        merge_arrays(x, t, na, a + na, nb, a, n);
    }
}

void
rank_ascending(const int *x, int *idx, int n)
{
    for (int i = 0; i < n; i++) idx[i] = i;
    int *t = XMALLOC(int, (n + 1) / 2);
    mergerank(x, idx, n, t);
    free(t);
}

/* ------------------------------------------------------------------ */
/* file_is_compressed                                                   */
/* ------------------------------------------------------------------ */
bool
file_is_compressed(int ncid)
{
    int ndims;
    if (nc_inq_ndims(ncid, &ndims) != NC_NOERR) return false;
    for (int dimid = 0; dimid < ndims; dimid++) {
        char name[NC_MAX_NAME + 1];
        int varid, attid;
        if (nc_inq_dimname(ncid, dimid, name) != NC_NOERR) continue;
        if (nc_inq_varid(ncid, name, &varid) == NC_NOERR &&
            nc_inq_attid(ncid, varid, "compress", &attid) == NC_NOERR)
            return true;
    }
    return false;
}

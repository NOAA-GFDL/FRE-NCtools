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
 * Domain decomposition helpers extracted from src/mpp-ncscatter/mppncscatter.c
 * and src/mpp-ncscatter/scatterdim.c.
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <netcdf.h>
#include "domain.h"
#include "xmalloc.h"

/* ------------------------------------------------------------------ */
/* ScatterDim                                                           */
/* ------------------------------------------------------------------ */
ScatterDim *
ScatterDim_new(int id, size_t len, const char *name,
               scatter_t scatter_type, int ndiv)
{
    ScatterDim *p = XMALLOC(ScatterDim, 1);
    p->id           = id;
    p->len          = len;
    /* Safe string copy with explicit null termination */
    size_t name_len = strlen(name);
    if (name_len > NC_MAX_NAME) name_len = NC_MAX_NAME;
    memcpy(p->name, name, name_len);
    p->name[name_len] = '\0';
    p->scatter_type = scatter_type;
    p->scatter_ndiv = (size_t)ndiv;
    p->scatter_start = XMALLOC(size_t, ndiv > 0 ? ndiv : 1);
    p->scatter_end   = XMALLOC(size_t, ndiv > 0 ? ndiv : 1);
    p->scatter_len   = XMALLOC(size_t, ndiv > 0 ? ndiv : 1);
    return p;
}

void
ScatterDim_free(ScatterDim *p)
{
    if (!p) return;
    XFREE(p->scatter_start);
    XFREE(p->scatter_end);
    XFREE(p->scatter_len);
    free(p);
}

void
scatter_dims_free(ScatterDim *dims[], int ndims)
{
    for (int i = 0; i < ndims; i++) {
        ScatterDim_free(dims[i]);
        dims[i] = NULL;
    }
}

/* ------------------------------------------------------------------ */
/* mpp_compute_extent (verbatim port from mppncscatter.c)              */
/* ------------------------------------------------------------------ */
#define EVEN(n) (!((n) & 1))
#define ODD(n)  (((n) & 1))
#define MAX(a,b) ((a)>(b)?(a):(b))

void
mpp_compute_extent(size_t isg, size_t ieg, size_t ndivs,
                   size_t *start, size_t *end)
{
    size_t n    = ieg - isg + 1;
    size_t iss  = isg;
    size_t ndiv, imax = ieg, ndmax = ndivs, ie, ndmirror;
    char   symmetrize;

    for (ndiv = 0; ndiv < ndivs; ++ndiv) {
        symmetrize = (EVEN(ndivs) && EVEN(n)) ||
                     (ODD(ndivs) && ODD(n)) ||
                     (ODD(ndivs) && EVEN(n) && (ndivs < (n / 2)));

        if (ndiv == 0) { imax = ieg; ndmax = ndivs; }

        if (ndiv < ((ndivs - 1) / 2 + 1)) {
            ie = iss + (size_t)ceil(((float)(imax - iss + 1.0)) / (ndmax - ndiv)) - 1;
            ndmirror = (ndivs - 1) - ndiv;
            if ((ndmirror > ndiv) && symmetrize) {
                start[ndmirror] = MAX(isg + ieg - ie, ie + 1);
                end[ndmirror]   = MAX(isg + ieg - iss, ie + 1);
                imax = start[ndmirror] - 1;
                ndmax = ndmax - 1;
            }
        } else {
            if (symmetrize) {
                iss = start[ndiv];
                ie  = end[ndiv];
            } else {
                ie = iss + (size_t)ceil(((float)(imax - iss + 1.0)) / (ndmax - ndiv)) - 1;
            }
        }
        start[ndiv] = iss;
        end[ndiv]   = ie;
        iss = ie + 1;
    }
}

/* ------------------------------------------------------------------ */
/* hyperslabcopy (verbatim port from mppncscatter.c)                   */
/* ------------------------------------------------------------------ */
void
hyperslabcopy(nc_type type, size_t *dimlen, int *dimids,
              size_t *start, size_t *count, int ndim,
              char *ti, short *si, int *ii, float *fi, double *di,
              char *t,  short *s,  int *iv, float *f,  double *d)
{
    size_t k, j, offset = 0, i0, stridek, stridej;
    switch (ndim) {
    case 0:
        switch (type) {
        case NC_BYTE: case NC_CHAR: *t = *ti; break;
        case NC_SHORT: *s = *si; break;
        case NC_INT:   *iv = *ii; break;
        case NC_FLOAT: *f = *fi; break;
        case NC_DOUBLE:*d = *di; break;
        default: break;
        }
        break;
    case 1:
        switch (type) {
        case NC_BYTE: case NC_CHAR: memcpy(t,  ti + start[0], count[0] * sizeof(char)); break;
        case NC_SHORT: memcpy(s,  si + start[0], count[0] * sizeof(short)); break;
        case NC_INT:   memcpy(iv, ii + start[0], count[0] * sizeof(int)); break;
        case NC_FLOAT: memcpy(f,  fi + start[0], count[0] * sizeof(float)); break;
        case NC_DOUBLE:memcpy(d,  di + start[0], count[0] * sizeof(double)); break;
        default: break;
        }
        break;
    case 2:
        i0      = start[0] * dimlen[dimids[1]];
        stridek = dimlen[dimids[1]];
        switch (type) {
        case NC_BYTE: case NC_CHAR:
            for (k = 0; k < count[0]; ++k) { memcpy(t+offset, ti+i0+k*stridek+start[1], count[1]*sizeof(char)); offset += count[1]; } break;
        case NC_SHORT:
            for (k = 0; k < count[0]; ++k) { memcpy(s+offset, si+i0+k*stridek+start[1], count[1]*sizeof(short)); offset += count[1]; } break;
        case NC_INT:
            for (k = 0; k < count[0]; ++k) { memcpy(iv+offset, ii+i0+k*stridek+start[1], count[1]*sizeof(int)); offset += count[1]; } break;
        case NC_FLOAT:
            for (k = 0; k < count[0]; ++k) { memcpy(f+offset, fi+i0+k*stridek+start[1], count[1]*sizeof(float)); offset += count[1]; } break;
        case NC_DOUBLE:
            for (k = 0; k < count[0]; ++k) { memcpy(d+offset, di+i0+k*stridek+start[1], count[1]*sizeof(double)); offset += count[1]; } break;
        default: break;
        }
        break;
    case 3:
        i0      = start[0]*dimlen[dimids[1]]*dimlen[dimids[2]] + start[1]*dimlen[dimids[2]] + start[2];
        stridek = dimlen[dimids[1]]*dimlen[dimids[2]];
        stridej = dimlen[dimids[2]];
        switch (type) {
        case NC_BYTE: case NC_CHAR:
            for (k = 0; k < count[0]; ++k) { for (j = 0; j < count[1]; ++j) { memcpy(t+offset, ti+i0+k*stridek+j*stridej, count[2]*sizeof(char)); offset += count[2]; } } break;
        case NC_SHORT:
            for (k = 0; k < count[0]; ++k) { for (j = 0; j < count[1]; ++j) { memcpy(s+offset, si+i0+k*stridek+j*stridej, count[2]*sizeof(short)); offset += count[2]; } } break;
        case NC_INT:
            for (k = 0; k < count[0]; ++k) { for (j = 0; j < count[1]; ++j) { memcpy(iv+offset, ii+i0+k*stridek+j*stridej, count[2]*sizeof(int)); offset += count[2]; } } break;
        case NC_FLOAT:
            for (k = 0; k < count[0]; ++k) { for (j = 0; j < count[1]; ++j) { memcpy(f+offset, fi+i0+k*stridek+j*stridej, count[2]*sizeof(float)); offset += count[2]; } } break;
        case NC_DOUBLE:
            for (k = 0; k < count[0]; ++k) { for (j = 0; j < count[1]; ++j) { memcpy(d+offset, di+i0+k*stridek+j*stridej, count[2]*sizeof(double)); offset += count[2]; } } break;
        default: break;
        }
        break;
    case 4:
        i0      = start[1]*dimlen[dimids[2]]*dimlen[dimids[3]] + start[2]*dimlen[dimids[3]] + start[3];
        stridek = dimlen[dimids[2]]*dimlen[dimids[3]];
        stridej = dimlen[dimids[3]];
        switch (type) {
        case NC_BYTE: case NC_CHAR:
            for (k = 0; k < count[1]; ++k) { for (j = 0; j < count[2]; ++j) { memcpy(t+offset, ti+i0+k*stridek+j*stridej, count[3]*sizeof(char)); offset += count[3]; } } break;
        case NC_SHORT:
            for (k = 0; k < count[0]; ++k) { for (j = 0; j < count[1]; ++j) { memcpy(s+offset, si+i0+k*stridek+j*stridej, count[2]*sizeof(short)); offset += count[2]; } } break;
        case NC_INT:
            for (k = 0; k < count[1]; ++k) { for (j = 0; j < count[2]; ++j) { memcpy(iv+offset, ii+i0+k*stridek+j*stridej, count[3]*sizeof(int)); offset += count[3]; } } break;
        case NC_FLOAT:
            for (k = 0; k < count[1]; ++k) { for (j = 0; j < count[2]; ++j) { memcpy(f+offset, fi+i0+k*stridek+j*stridej, count[3]*sizeof(float)); offset += count[3]; } } break;
        case NC_DOUBLE:
            for (k = 0; k < count[1]; ++k) { for (j = 0; j < count[2]; ++j) { memcpy(d+offset, di+i0+k*stridek+j*stridej, count[3]*sizeof(double)); offset += count[3]; } } break;
        default: break;
        }
        break;
    }
}

/* ------------------------------------------------------------------ */
/* Internal scatter helpers                                             */
/* ------------------------------------------------------------------ */
static void
get_num_divs(scatter_opts_t *opts, int *nx, int *ny)
{
    if (!opts) { *nx = *ny = 0; return; }
    if (opts->nxio && opts->nyio) { *nx = opts->nxio; *ny = opts->nyio; }
    else { *nx = opts->nx; *ny = opts->ny; }
}

int
scatter_get_num_files(scatter_opts_t *opts)
{
    int nx, ny;
    if (!opts) return 0;
    get_num_divs(opts, &nx, &ny);
    return nx * ny;
}

static void
get_scatter_dims_from_file(int nc, int ndims, int nvars,
                            ScatterDim *scatterdims[], scatter_opts_t *opt)
{
    int recid;
    if (nc_inq_unlimdim(nc, &recid) != NC_NOERR) recid = -1;

    for (int dimid = 0; dimid < ndims; dimid++) {
        scatter_t stype = NOSCATTER;
        size_t dimlen = 0;
        int ndiv = 0;
        char name[NC_MAX_NAME + 1];
        int varid;
        nc_type att_type;
        char att[256];

        if (nc_inq_dimname(nc, dimid, name) != NC_NOERR) goto done;

        if (dimid == recid) { dimlen = NC_UNLIMITED; goto done; }
        if (nc_inq_dimlen(nc, dimid, &dimlen) != NC_NOERR) goto done;

        if (instringlist(opt->xdims, name, opt->xdims_len)) { stype = SCATTERX; goto done; }
        if (instringlist(opt->ydims, name, opt->ydims_len)) { stype = SCATTERY; goto done; }

        if (nc_inq_varid(nc, name, &varid) != NC_NOERR) goto done;

        /* Check units attribute for degrees_east / degrees_north */
        int found = 0;
        if (nc_inq_atttype(nc, varid, "units", &att_type) == NC_NOERR &&
            att_type == NC_CHAR) {
            memset(att, 0, sizeof(att));
            if (nc_get_att_text(nc, varid, "units", att) == NC_NOERR) {
                if (!strncasecmp(att, "degrees_east",12) || !strncasecmp(att,"degree_east",11) ||
                    !strncasecmp(att, "degrees_e",9)     || !strncasecmp(att,"degree_e",8)     ||
                    !strncasecmp(att, "degreee",7)        || !strncasecmp(att,"degreese",8)) {
                    stype = SCATTERX; found = 1;
                } else if (!strncasecmp(att,"degrees_north",13) || !strncasecmp(att,"degree_north",12) ||
                           !strncasecmp(att,"degrees_n",9)       || !strncasecmp(att,"degree_n",8)     ||
                           !strncasecmp(att,"degreesn",8)         || !strncasecmp(att,"degreen",7)) {
                    stype = SCATTERY; found = 1;
                }
            }
        }
        if (!found) {
            if (nc_inq_atttype(nc, varid, "cartesian_axis", &att_type) == NC_NOERR &&
                att_type == NC_CHAR) {
                memset(att, 0, sizeof(att));
                if (nc_get_att_text(nc, varid, "cartesian_axis", att) == NC_NOERR) {
                    if (!strncasecmp(att, "x", 1)) stype = SCATTERX;
                    else if (!strncasecmp(att, "y", 1)) stype = SCATTERY;
                }
            }
        }

    done:
        ndiv = (stype == SCATTERX) ? opt->nx : (stype == SCATTERY) ? opt->ny : 0;
        scatterdims[dimid] = ScatterDim_new(dimid, dimlen, name, stype, ndiv);
    }
}

static void
get_scatter_extents(ScatterDim *scatterdims[], int ndims)
{
    for (int i = 0; i < ndims; i++) {
        ScatterDim *pd = scatterdims[i];
        if (!pd || pd->scatter_type == NOSCATTER) continue;
        mpp_compute_extent(0, pd->len - 1, pd->scatter_ndiv,
                           pd->scatter_start, pd->scatter_end);
        for (size_t j = 0; j < pd->scatter_ndiv; j++)
            pd->scatter_len[j] = pd->scatter_end[j] - pd->scatter_start[j] + 1;
    }
}

static void
get_scatter_extents_iolayout(ScatterDim *scatterdims[], int ndims,
                              int nxio, int nyio)
{
    for (int d = 0; d < ndims; d++) {
        ScatterDim *pd = scatterdims[d];
        if (!pd || pd->scatter_type == NOSCATTER) continue;
        int ndivio = (pd->scatter_type == SCATTERX) ? nxio : nyio;
        int step = (int)pd->scatter_ndiv / ndivio;
        size_t *nstart = XMALLOC(size_t, ndivio);
        size_t *nend   = XMALLOC(size_t, ndivio);
        size_t *nlen   = XMALLOC(size_t, ndivio);
        int i = 0, k = 0;
        while (i < (int)pd->scatter_ndiv) {
            nstart[k] = pd->scatter_start[i];
            nend[k]   = pd->scatter_end[i + step - 1];
            nlen[k]   = nend[k] - nstart[k] + 1;
            k++; i += step;
        }
        pd->scatter_ndiv = (size_t)ndivio;
        XFREE(pd->scatter_start); pd->scatter_start = nstart;
        XFREE(pd->scatter_end);   pd->scatter_end   = nend;
        XFREE(pd->scatter_len);   pd->scatter_len   = nlen;
    }
}

void
scatter_dims(int nc, int ndims, int nvars,
             ScatterDim *scatterdims[], scatter_opts_t *opt)
{
    get_scatter_dims_from_file(nc, ndims, nvars, scatterdims, opt);
    get_scatter_extents(scatterdims, ndims);
    if (opt->nxio && opt->nyio) {
        if (opt->nx % opt->nxio) {
            fprintf(stderr, "Error: x divisions not wholly divisible by io-layout x.\n");
            exit(1);
        }
        if (opt->ny % opt->nyio) {
            fprintf(stderr, "Error: y divisions not wholly divisible by io-layout y.\n");
            exit(1);
        }
        get_scatter_extents_iolayout(scatterdims, ndims, opt->nxio, opt->nyio);
    }
}

/* ------------------------------------------------------------------ */
/* scatter_def_dim                                                      */
/* ------------------------------------------------------------------ */
void
scatter_def_dim(int nc, int *ncids, int ndims,
                ScatterDim *scatterdims[], scatter_opts_t *opts)
{
    int nx, ny, dummy;
    get_num_divs(opts, &nx, &ny);

    for (int dimid = 0; dimid < ndims; dimid++) {
        ScatterDim *sd = scatterdims[dimid];
        if (!sd) continue;
        for (int yi = 0; yi < ny; yi++) {
            for (int xi = 0; xi < nx; xi++) {
                int i = yi * nx + xi;
                if (opts->dryrun) continue;
                switch (sd->scatter_type) {
                case NOSCATTER:
                    nc_def_dim(ncids[i], sd->name, sd->len, &dummy); break;
                case SCATTERX:
                    nc_def_dim(ncids[i], sd->name, sd->scatter_len[xi], &dummy); break;
                case SCATTERY:
                    nc_def_dim(ncids[i], sd->name, sd->scatter_len[yi], &dummy); break;
                }
            }
        }
    }
}

/* ------------------------------------------------------------------ */
/* scatter_def_var                                                      */
/* ------------------------------------------------------------------ */
void
scatter_def_var(int nc, int *ncids, int nvars, int ndims,
                ScatterDim *scatterdims[], scatter_opts_t *opts)
{
    int nx, ny;
    get_num_divs(opts, &nx, &ny);
    int scatatt[4];
    scatatt[0] = 1;

    for (int varid = 0; varid < nvars; varid++) {
        char varname[NC_MAX_NAME + 1], attname[NC_MAX_NAME + 1];
        nc_type type;
        int ndimvar, dimids[NC_MAX_DIMS], natt, newvarid;

        if (nc_inq_var(nc, varid, varname, &type, &ndimvar, dimids, &natt) != NC_NOERR) continue;

        for (int yi = 0; yi < ny; yi++) {
            for (int xi = 0; xi < nx; xi++) {
                int ifile = yi * nx + xi;
                if (!opts->dryrun) {
                    if (nc_def_var(ncids[ifile], varname, type, ndimvar, dimids, &newvarid) != NC_NOERR) {
                        fprintf(stderr, "Error: failed to define var '%s' in file %d.\n", varname, ifile);
                        exit(-1);
                    }
                    for (int j = 0; j < natt; j++) {
                        nc_inq_attname(nc, varid, j, attname);
                        nc_copy_att(nc, varid, attname, ncids[ifile], newvarid);
                    }
                }
                /* Add domain_decomposition attribute for 1D dimension vars. */
                if (ndimvar == 1) {
                    ScatterDim *sd = scatterdims[dimids[0]];
                    if (!sd) continue;
                    if (sd->scatter_type == SCATTERX && strcmp(varname, sd->name) == 0) {
                        scatatt[1] = (int)sd->len;
                        scatatt[2] = (int)sd->scatter_start[xi] + 1;
                        scatatt[3] = (int)(sd->scatter_start[xi] + sd->scatter_len[xi]);
                        if (!opts->dryrun)
                            nc_put_att_int(ncids[ifile], varid, "domain_decomposition", NC_INT, 4, scatatt);
                    } else if (sd->scatter_type == SCATTERY && strcmp(varname, sd->name) == 0) {
                        scatatt[1] = (int)sd->len;
                        scatatt[2] = (int)sd->scatter_start[yi] + 1;
                        scatatt[3] = (int)(sd->scatter_start[yi] + sd->scatter_len[yi]);
                        if (!opts->dryrun)
                            nc_put_att_int(ncids[ifile], varid, "domain_decomposition", NC_INT, 4, scatatt);
                    }
                }
            }
        }
    }
}

/* ------------------------------------------------------------------ */
/* scatter_put_var                                                      */
/* ------------------------------------------------------------------ */
void
scatter_put_var(int nc, int *ncids, int ndims, int nvars,
                ScatterDim *scatterdims[], scatter_opts_t *opt)
{
    if (opt->dryrun) return;

    int nx, ny;
    get_num_divs(opt, &nx, &ny);

    int recid;
    size_t nrec = 0;
    nc_inq_unlimdim(nc, &recid);
    if (recid != -1) nc_inq_dimlen(nc, recid, &nrec);

    size_t dimlen[NC_MAX_DIMS];
    for (int i = 0; i < ndims; i++) nc_inq_dimlen(nc, i, &dimlen[i]);

    /* Find maximum buffer sizes needed. */
    size_t maxchar = 0, maxshort = 0, maxint = 0, maxfloat = 0, maxdouble = 0;
    for (int varid = 0; varid < nvars; varid++) {
        nc_type type;
        int ndimvar, dimids[NC_MAX_DIMS], natt;
        nc_inq_var(nc, varid, NULL, &type, &ndimvar, dimids, &natt);
        size_t sz = 1;
        for (int i = 0; i < ndimvar; i++)
            if (dimids[i] != recid) sz *= dimlen[dimids[i]];
        switch (type) {
        case NC_BYTE: case NC_CHAR: if (sz > maxchar)   maxchar   = sz; break;
        case NC_SHORT:              if (sz > maxshort)  maxshort  = sz; break;
        case NC_INT:                if (sz > maxint)    maxint    = sz; break;
        case NC_FLOAT:              if (sz > maxfloat)  maxfloat  = sz; break;
        case NC_DOUBLE:             if (sz > maxdouble) maxdouble = sz; break;
        default: break;
        }
    }

    char   *tp = maxchar   ? (char   *)malloc(maxchar   * sizeof(char))   : NULL;
    char   *otp= maxchar   ? (char   *)malloc(maxchar   * sizeof(char))   : NULL;
    short  *sp = maxshort  ? (short  *)malloc(maxshort  * sizeof(short))  : NULL;
    short  *osp= maxshort  ? (short  *)malloc(maxshort  * sizeof(short))  : NULL;
    int    *ip = maxint    ? (int    *)malloc(maxint    * sizeof(int))     : NULL;
    int    *oip= maxint    ? (int    *)malloc(maxint    * sizeof(int))     : NULL;
    float  *fp = maxfloat  ? (float  *)malloc(maxfloat  * sizeof(float))  : NULL;
    float  *ofp= maxfloat  ? (float  *)malloc(maxfloat  * sizeof(float))  : NULL;
    double *dp = maxdouble ? (double *)malloc(maxdouble * sizeof(double)) : NULL;
    double *odp= maxdouble ? (double *)malloc(maxdouble * sizeof(double)) : NULL;

    size_t outstart[4] = {0, 0, 0, 0};
    size_t instart[4], count[4];

    /* Process non-record variables. */
    for (int varid = 0; varid < nvars; varid++) {
        char varname[NC_MAX_NAME+1];
        nc_type type;
        int ndimvar, dimids[NC_MAX_DIMS], natt;
        nc_inq_var(nc, varid, varname, &type, &ndimvar, dimids, &natt);
        if (recid != -1 && ndimvar > 0 && dimids[0] == recid) continue;

        switch (type) {
        case NC_BYTE: case NC_CHAR: nc_get_var_text(nc, varid, tp);   break;
        case NC_SHORT:              nc_get_var_short(nc, varid, sp);   break;
        case NC_INT:                nc_get_var_int(nc, varid, ip);     break;
        case NC_FLOAT:              nc_get_var_float(nc, varid, fp);   break;
        case NC_DOUBLE:             nc_get_var_double(nc, varid, dp);  break;
        default: break;
        }

        for (int yi = 0; yi < ny; yi++) {
            for (int xi = 0; xi < nx; xi++) {
                int i = yi * nx + xi;
                if (ndimvar > 0) {
                    for (int dimi = 0; dimi < ndimvar; dimi++) {
                        ScatterDim *sd = scatterdims[dimids[dimi]];
                        if (!sd) continue;
                        switch (sd->scatter_type) {
                        case NOSCATTER: instart[dimi]=0; count[dimi]=sd->len; break;
                        case SCATTERX:  instart[dimi]=sd->scatter_start[xi]; count[dimi]=sd->scatter_len[xi]; break;
                        case SCATTERY:  instart[dimi]=sd->scatter_start[yi]; count[dimi]=sd->scatter_len[yi]; break;
                        }
                    }
                    hyperslabcopy(type, dimlen, dimids, instart, count, ndimvar,
                                  tp, sp, ip, fp, dp, otp, osp, oip, ofp, odp);
                    switch (type) {
                    case NC_BYTE: case NC_CHAR: nc_put_vara_text(ncids[i], varid, outstart, count, otp); break;
                    case NC_SHORT:              nc_put_vara_short(ncids[i], varid, outstart, count, osp); break;
                    case NC_INT:                nc_put_vara_int(ncids[i], varid, outstart, count, oip); break;
                    case NC_FLOAT:              nc_put_vara_float(ncids[i], varid, outstart, count, ofp); break;
                    case NC_DOUBLE:             nc_put_vara_double(ncids[i], varid, outstart, count, odp); break;
                    default: break;
                    }
                } else {
                    switch (type) {
                    case NC_BYTE: case NC_CHAR: nc_put_var_text(ncids[i], varid, tp); break;
                    case NC_SHORT:              nc_put_var_short(ncids[i], varid, sp); break;
                    case NC_INT:                nc_put_var_int(ncids[i], varid, ip); break;
                    case NC_FLOAT:              nc_put_var_float(ncids[i], varid, fp); break;
                    case NC_DOUBLE:             nc_put_var_double(ncids[i], varid, dp); break;
                    default: break;
                    }
                }
            }
        }
    }

    /* Process record variables. */
    for (size_t reci = 0; reci < nrec; reci++) {
        instart[0] = reci; outstart[0] = reci; count[0] = 1;
        for (int varid = 0; varid < nvars; varid++) {
            char varname[NC_MAX_NAME+1];
            nc_type type;
            int ndimvar, dimids[NC_MAX_DIMS], natt, dimids2[NC_MAX_DIMS];
            nc_inq_var(nc, varid, varname, &type, &ndimvar, dimids, &natt);
            if (dimids[0] != recid) continue;
            int ndimvar2 = ndimvar - 1;

            for (int dimi = 1; dimi < ndimvar; dimi++) {
                count[dimi] = dimlen[dimids[dimi]];
                dimids2[dimi-1] = dimids[dimi];
            }

            switch (type) {
            case NC_BYTE: case NC_CHAR: nc_get_vara_text(nc,varid,instart,count,tp); break;
            case NC_SHORT:              nc_get_vara_short(nc,varid,instart,count,sp); break;
            case NC_INT:                nc_get_vara_int(nc,varid,instart,count,ip); break;
            case NC_FLOAT:              nc_get_vara_float(nc,varid,instart,count,fp); break;
            case NC_DOUBLE:             nc_get_vara_double(nc,varid,instart,count,dp); break;
            default: break;
            }

            for (int yi = 0; yi < ny; yi++) {
                for (int xi = 0; xi < nx; xi++) {
                    int i = yi * nx + xi;
                    size_t instart2[4]={0,0,0,0}, count2[4]={1,1,1,1};
                    for (int dimi = 1; dimi < ndimvar; dimi++) {
                        ScatterDim *sd = scatterdims[dimids[dimi]];
                        if (!sd) continue;
                        switch (sd->scatter_type) {
                        case NOSCATTER: instart[dimi]=0; count[dimi]=sd->len; break;
                        case SCATTERX:  instart[dimi]=sd->scatter_start[xi]; count[dimi]=sd->scatter_len[xi]; break;
                        case SCATTERY:  instart[dimi]=sd->scatter_start[yi]; count[dimi]=sd->scatter_len[yi]; break;
                        }
                        count2[dimi-1] = count[dimi];
                        instart2[dimi-1] = instart[dimi];
                    }
                    hyperslabcopy(type, dimlen, dimids2, instart2, count2, ndimvar2,
                                  tp, sp, ip, fp, dp, otp, osp, oip, ofp, odp);
                    switch (type) {
                    case NC_BYTE: case NC_CHAR: nc_put_vara_text(ncids[i],varid,outstart,count,otp); break;
                    case NC_SHORT:              nc_put_vara_short(ncids[i],varid,outstart,count,osp); break;
                    case NC_INT:                nc_put_vara_int(ncids[i],varid,outstart,count,oip); break;
                    case NC_FLOAT:              nc_put_vara_float(ncids[i],varid,outstart,count,ofp); break;
                    case NC_DOUBLE:             nc_put_vara_double(ncids[i],varid,outstart,count,odp); break;
                    default: break;
                    }
                }
            }
        }
    }

    free(tp); free(otp); free(sp); free(osp);
    free(ip); free(oip); free(fp); free(ofp);
    free(dp); free(odp);
}

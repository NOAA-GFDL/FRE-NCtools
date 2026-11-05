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
 * Three combine paths:
 *   combine_mpp      – port of mppnccombine.c (modernised API)
 *   combine_land     – port of combine-ncc.F90
 *   combine_iceberg  – port of iceberg_comb.sh (no ncrcat/ncatted)
 */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>
#include <netcdf.h>
#include "cmd_combine.h"
#include "cmd_check.h"
#include "nc_utils.h"
#include "compress.h"
#include "xmalloc.h"

#define NC_CHECK(call) do { \
    int _err = (call); \
    if (_err != NC_NOERR) { \
        fprintf(stderr, "NetCDF error: %s\n", nc_strerror(_err)); \
        return 1; \
    } \
} while (0)

#ifndef MAX_BF
#define MAX_BF 100
#endif
#ifndef DEFAULT_BF
#define DEFAULT_BF 1
#endif
#ifndef NC_BLKSZ
#define NC_BLKSZ 65536
#endif

/* ================================================================== */
/* combine_mpp – modern port of mppnccombine.c                        */
/* ================================================================== */

/* Per-file metadata for mppnccombine algorithm. */
typedef struct {
    int     ncid;
    int     ndims, nvars, ngatts, recdim;
    char    varname[NC_MAX_VARS][NC_MAX_NAME + 1];
    nc_type datatype[NC_MAX_VARS];
    int     varndims[NC_MAX_VARS];
    int     vardim[NC_MAX_VARS][NC_MAX_DIMS];
    int     natts[NC_MAX_VARS];
    unsigned char vardecomp[NC_MAX_VARS];
    char    dimname[NC_MAX_DIMS][NC_MAX_NAME + 1];
    size_t  dimsize[NC_MAX_DIMS];
    size_t  dimfullsize[NC_MAX_DIMS];
    size_t  dimstart[NC_MAX_DIMS]; /* 0-based */
    size_t  dimend[NC_MAX_DIMS];   /* 0-based; SIZE_MAX = not decomposed */
    unsigned char varmiss[NC_MAX_VARS];
    unsigned char varmissval[NC_MAX_VARS][8];
} fileinfo_t;

static void   ***varbuf = NULL;
static size_t    mem_allocated = 0;
static size_t    estimated_maxrss = 0;
static int       g_mem_dry_run = 0;
static int       g_print_mem  = 0;

static int g_min(int a, int b) { return a < b ? a : b; }

static size_t nc_type_sz(nc_type t) {
    switch (t) {
    case NC_BYTE: case NC_CHAR: return 1;
    case NC_SHORT:  return 2;
    case NC_INT:    return 4;
    case NC_FLOAT:  return 4;
    case NC_DOUBLE: return 8;
    default:        return 0;
    }
}

/* Generic vara read/write by nc_type. */
static int get_vara_typed(int ncid, int varid, const size_t *st,
                           const size_t *cnt, void *buf, nc_type t) {
    switch (t) {
    case NC_BYTE: case NC_CHAR: return nc_get_vara_text(ncid,varid,st,cnt,(char*)buf);
    case NC_SHORT:  return nc_get_vara_short(ncid,varid,st,cnt,(short*)buf);
    case NC_INT:    return nc_get_vara_int(ncid,varid,st,cnt,(int*)buf);
    case NC_FLOAT:  return nc_get_vara_float(ncid,varid,st,cnt,(float*)buf);
    case NC_DOUBLE: return nc_get_vara_double(ncid,varid,st,cnt,(double*)buf);
    default: return NC_EBADTYPE;
    }
}

static int put_vara_typed(int ncid, int varid, const size_t *st,
                           const size_t *cnt, const void *buf, nc_type t) {
    switch (t) {
    case NC_BYTE: case NC_CHAR: return nc_put_vara_text(ncid,varid,st,cnt,(const char*)buf);
    case NC_SHORT:  return nc_put_vara_short(ncid,varid,st,cnt,(const short*)buf);
    case NC_INT:    return nc_put_vara_int(ncid,varid,st,cnt,(const int*)buf);
    case NC_FLOAT:  return nc_put_vara_float(ncid,varid,st,cnt,(const float*)buf);
    case NC_DOUBLE: return nc_put_vara_double(ncid,varid,st,cnt,(const double*)buf);
    default: return NC_EBADTYPE;
    }
}

static int mpp_flush_decomp(fileinfo_t *outf, int nfiles, int r, int bf, int verbose) {
    for (int v = 0; v < outf->nvars; v++) {
        if (!outf->vardecomp[v]) continue;
        size_t outstart[NC_MAX_DIMS] = {0};
        size_t count[NC_MAX_DIMS];
        int varrecdim = -1;
        for (int d = 0; d < outf->varndims[v]; d++) {
            outstart[d] = 0;
            if (outf->vardim[v][d] == outf->recdim) {
                count[d] = 1; varrecdim = d;
            } else {
                count[d] = outf->dimfullsize[outf->vardim[v][d]];
            }
        }
        if (varrecdim >= 0) outstart[varrecdim] = (size_t)r;
        if (varrecdim < 0 && r > 0) continue;
        if (!varbuf || !varbuf[r % bf] || !varbuf[r % bf][v]) continue;
        if (put_vara_typed(outf->ncid, v, outstart, count,
                           varbuf[r % bf][v], outf->datatype[v]) != NC_NOERR) {
            fprintf(stderr, "Error: cannot write variable %d.\n", v);
            return 1;
        }
    }
    return 0;
}

static int mpp_process_vars(fileinfo_t *inf, fileinfo_t *outf,
                             int appendnc, int *nrecs, int *nblocks, int *bf,
                             int r, int nfiles, int f, int verbose, int missing) {
    int recdimsize = (inf->recdim < 0) ? 1 : (int)inf->dimsize[inf->recdim];

    if (*nrecs == 1) {
        *nrecs = recdimsize;
        if (*bf >= 1) {
            if (*bf > *nrecs) *bf = *nrecs;
        } else {
            *bf = g_min(MAX_BF, *nrecs);
        }
        *nblocks = ((*nrecs % *bf) != 0) ? (*nrecs / *bf) + 1 : (*nrecs / *bf);
    } else if (recdimsize != *nrecs) {
        fprintf(stderr, "Error: different number of records than first input.\n");
        return 1;
    }

    /* Allocate varbuf on first call. */
    if (!varbuf) {
        /* Validate buffer factor to prevent overflow */
        if (*bf <= 0 || *bf > MAX_BF) {
            fprintf(stderr, "Invalid buffer factor: %d\n", *bf);
            return 1;
        }
        /* Check for overflow: bf * sizeof(void**) */
        if ((size_t)(*bf) > SIZE_MAX / sizeof(void**)) {
            fprintf(stderr, "Buffer allocation size overflow\n");
            return 1;
        }
        size_t ptr_array_size = (size_t)(*bf) * sizeof(void**);
        /* Check for overflow: bf * NC_MAX_VARS * sizeof(void*) */
        if ((size_t)(*bf) > SIZE_MAX / NC_MAX_VARS ||
            (size_t)(*bf) * NC_MAX_VARS > SIZE_MAX / sizeof(void*)) {
            fprintf(stderr, "Buffer allocation size overflow\n");
            return 1;
        }
        size_t data_array_size = (size_t)(*bf) * NC_MAX_VARS * sizeof(void*);
        /* Check for overflow when adding the two sizes */
        if (ptr_array_size > SIZE_MAX - data_array_size) {
            fprintf(stderr, "Buffer allocation size overflow\n");
            return 1;
        }
        size_t nbytes = ptr_array_size + data_array_size;
        if (g_mem_dry_run) estimated_maxrss += nbytes;
        varbuf = (void***)calloc(nbytes, 1);
        if (!varbuf) { fprintf(stderr,"Cannot allocate varbuf.\n"); exit(1); }
        for (int z = 0; z < *bf; z++)
            varbuf[z] = (void**)((size_t)varbuf + (*bf)*sizeof(void**) + z*NC_MAX_VARS*sizeof(void*));
    }

    for (int v = 0; v < inf->nvars; v++) {
        size_t recsize = 1, recfullsize = 1;
        int varrecdim = -1;
        size_t instart[NC_MAX_DIMS] = {0};
        size_t outstart[NC_MAX_DIMS] = {0};
        size_t count[NC_MAX_DIMS];

        for (int d = 0; d < inf->varndims[v]; d++) {
            if (inf->vardim[v][d] == inf->recdim) {
                count[d] = 1; varrecdim = d;
            } else {
                count[d] = inf->dimsize[inf->vardim[v][d]];
                recsize *= count[d];
                instart[d] = 0;
                outstart[d] = (inf->dimend[inf->vardim[v][d]] == (size_t)(-1)) ? 0 : inf->dimstart[inf->vardim[v][d]];
                recfullsize *= inf->dimfullsize[inf->vardim[v][d]];
            }
        }

        /* Skip if not the right record / already written. */
        if (r > 0) {
            int dimid_check;
            bool is_dim = (nc_inq_dimid(inf->ncid, inf->varname[v], &dimid_check) == NC_NOERR);
            if (is_dim) {
                if (dimid_check == inf->recdim) { if (f != 0) continue; }
                else continue;
            } else {
                if (varrecdim < 0) continue;
                if (!inf->vardecomp[v] && f > 0) continue;
            }
        } else {
            if (!inf->vardecomp[v] && appendnc) continue;
        }

        void *values = malloc(nc_type_sz(inf->datatype[v]) * recsize);
        if (!values) {
            fprintf(stderr, "Cannot allocate values buffer.\n");
            return 1;
        }

        if (varrecdim >= 0) instart[varrecdim] = outstart[varrecdim] = (size_t)r;

        if (get_vara_typed(inf->ncid, v, instart, count, values, inf->datatype[v]) != NC_NOERR) {
            fprintf(stderr,"Error reading variable '%s'.\n", inf->varname[v]);
            free(values); return 1;
        }

        if (!inf->vardecomp[v] && !g_mem_dry_run) {
            if (put_vara_typed(outf->ncid, v, outstart, count, values, inf->datatype[v]) != NC_NOERR) {
                fprintf(stderr,"Error writing variable '%s'.\n", inf->varname[v]);
                free(values); return 1;
            }
        } else {
            /* Buffer decomposed var. */
            if (!varbuf[r % (*bf)][v]) {
                size_t vsz = nc_type_sz(inf->datatype[v]) * recfullsize;
                if (g_mem_dry_run) { estimated_maxrss += vsz; varbuf[r%(*bf)][v] = (void*)"x"; free(values); continue; }
                varbuf[r % (*bf)][v] = calloc(vsz, 1);
                if (!varbuf[r % (*bf)][v]) { fprintf(stderr,"Cannot allocate varbuf entry.\n"); free(values); return 1; }
                mem_allocated += vsz;
                if (missing && outf->varmiss[v]) {
                    /* Fill with missing_value. */
                    size_t esz = nc_type_sz(inf->datatype[v]);
                    for (size_t s = 0; s < recfullsize; s++)
                        memcpy((char*)varbuf[r%(*bf)][v] + s*esz, outf->varmissval[v], esz);
                }
            }

            /* Scatter decomposed data into full buffer. */
            int ndv = inf->varndims[v];
            if (ndv == 0) {
                memcpy(varbuf[r%(*bf)][v], values, nc_type_sz(inf->datatype[v]));
                free(values); continue;
            }

            /* Build stride/offset for scattering. */
            size_t imax, jmax=1, kmax=1, lmax=1;
            size_t imaxfull, jmaxfull=1, kmaxfull=1;

            imax = inf->dimsize[inf->vardim[v][ndv-1]];
            imaxfull = inf->dimfullsize[inf->vardim[v][ndv-1]];
            if (ndv > 1 && inf->vardim[v][ndv-2] != inf->recdim) {
                jmax = inf->dimsize[inf->vardim[v][ndv-2]];
                jmaxfull = inf->dimfullsize[inf->vardim[v][ndv-2]];
            }
            if (ndv > 2 && inf->vardim[v][ndv-3] != inf->recdim) {
                kmax = inf->dimsize[inf->vardim[v][ndv-3]];
                kmaxfull = inf->dimfullsize[inf->vardim[v][ndv-3]];
            }
            if (ndv > 3 && inf->vardim[v][ndv-4] != inf->recdim) {
                lmax = inf->dimsize[inf->vardim[v][ndv-4]];
                (void)inf->dimfullsize[inf->vardim[v][ndv-4]]; /* lmaxfull not used in scatter */
            }

            size_t ioff = (outstart[ndv-1]);
            size_t joff = (ndv > 1) ? outstart[ndv-2] : 0;
            size_t koff = (ndv > 2) ? outstart[ndv-3] : 0;
            size_t loff = (ndv > 3) ? outstart[ndv-4] : 0;
            if (varrecdim >= 0) {
                if (ndv==1) ioff=0; else if (ndv==2) joff=0;
                else if (ndv==3) koff=0; else loff=0;
            }

            size_t ijfull  = imaxfull * jmaxfull;
            size_t ijkfull = ijfull   * kmaxfull;
            size_t esz     = nc_type_sz(inf->datatype[v]);
            size_t b = 0;
            for (size_t l = 0; l < lmax; l++)
                for (size_t k = 0; k < kmax; k++)
                    for (size_t j = 0; j < jmax; j++)
                        for (size_t i = 0; i < imax; i++) {
                            size_t off = (i+ioff) + (j+joff)*imaxfull + (k+koff)*ijfull + (l+loff)*ijkfull;
                            memcpy((char*)varbuf[r%(*bf)][v] + off*esz, (char*)values + b*esz, esz);
                            b++;
                        }
        }
        free(values);
    }
    return 0;
}

static int mpp_process_file(const char *ncname, int appendnc,
                             fileinfo_t *outf, const char *outncname,
                             int *nfiles, int *nrecs, int *nblocks, int *bf,
                             int block, int f, int headerpad, int verbose,
                             int missing, int deflate, int deflation, int shuffle) {
    fileinfo_t *inf = XMALLOC(fileinfo_t, 1);
    size_t blksz = NC_BLKSZ;

    if (nc_open(ncname, NC_NOWRITE, &inf->ncid) != NC_NOERR) {
        fprintf(stderr,"Error: cannot open '%s'\n", ncname);
        free(inf); return 1;
    }

    int nfiles2 = -1;
    nc_get_att_int(inf->ncid, NC_GLOBAL, "NumFilesInSet", &nfiles2);
    if (nfiles2 > 0) *nfiles = nfiles2;

    NC_CHECK(nc_inq(inf->ncid, &inf->ndims, &inf->nvars, &inf->ngatts, &inf->recdim));

    for (int d = 0; d < inf->ndims; d++) {
        NC_CHECK(nc_inq_dim(inf->ncid, d, inf->dimname[d], &inf->dimsize[d]));
        inf->dimfullsize[d] = inf->dimsize[d];
        inf->dimstart[d] = 0;
        inf->dimend[d]   = (size_t)(-1); /* not decomposed */
    }

    if (block == 0 && !g_mem_dry_run) {
        outf->nvars  = inf->nvars;
        outf->recdim = inf->recdim;
    }

    for (int v = 0; v < inf->nvars; v++) {
        NC_CHECK(nc_inq_var(inf->ncid, v, inf->varname[v], &inf->datatype[v],
                            &inf->varndims[v], inf->vardim[v], &inf->natts[v]));

        int dimid_check;
        if (nc_inq_dimid(inf->ncid, inf->varname[v], &dimid_check) == NC_NOERR) {
            int decomp[4];
            if (nc_get_att_int(inf->ncid, v, "domain_decomposition", decomp) == NC_NOERR) {
                inf->dimfullsize[dimid_check] = (size_t)(decomp[1] - decomp[0] + 1);
                inf->dimstart[dimid_check]    = (size_t)(decomp[2] - decomp[0]); /* 0-based */
                inf->dimend[dimid_check]      = (size_t)(decomp[3] - decomp[0]); /* 0-based */
            }
        }
    }

    for (int v = 0; v < inf->nvars; v++) {
        inf->vardecomp[v] = 0;
        for (int d = 0; d < inf->varndims[v]; d++) {
            if (inf->dimend[inf->vardim[v][d]] != (size_t)(-1)) {
                inf->vardecomp[v] = 1; break;
            }
        }
        if (!appendnc && !g_mem_dry_run) {
            if (block == 0) {
                outf->varndims[v] = inf->varndims[v];
                for (int d = 0; d < inf->ndims; d++)
                    outf->dimfullsize[d] = inf->dimfullsize[d];
                for (int d = 0; d < inf->varndims[v]; d++)
                    outf->vardim[v][d] = inf->vardim[v][d];
                outf->vardecomp[v] = inf->vardecomp[v];
                outf->datatype[v]  = inf->datatype[v];
                strcpy(outf->varname[v], inf->varname[v]);
                outf->varmiss[v] = 0;
            }
        }
    }

    if (!appendnc && !g_mem_dry_run) {
        /* Determine output format. */
        int ncinformat, ncoutformat;
        nc_inq_format(inf->ncid,  &ncinformat);
        nc_inq_format(outf->ncid, &ncoutformat);

        if (ncoutformat == NC_FORMAT_CLASSIC) {
            int newmode;
            switch (ncinformat) {
            case NC_FORMAT_CLASSIC: case NC_FORMAT_64BIT_OFFSET:
                newmode = NC_CLOBBER | NC_64BIT_OFFSET; break;
            case NC_FORMAT_NETCDF4:
                newmode = NC_CLOBBER | NC_NETCDF4; break;
            case NC_FORMAT_NETCDF4_CLASSIC:
                newmode = NC_CLOBBER | NC_NETCDF4 | NC_CLASSIC_MODEL; break;
            default:
                newmode = NC_CLOBBER | NC_64BIT_OFFSET; break;
            }
            nc_close(outf->ncid);
            if (nc__create(outncname, newmode, 0, &blksz, &outf->ncid) != NC_NOERR) {
                fprintf(stderr,"Error: cannot recreate output '%s'.\n", outncname);
                free(inf); return 1;
            }
            nc_set_fill(outf->ncid, NC_NOFILL, NULL);
        }

        /* Define dimensions in output. */
        for (int d = 0; d < inf->ndims; d++) {
            int dummy;
            size_t dlen = (d == inf->recdim) ? NC_UNLIMITED : inf->dimfullsize[d];
            nc_def_dim(outf->ncid, inf->dimname[d], dlen, &dummy);
        }

        /* Define variables + copy attributes. */
        for (int v = 0; v < inf->nvars; v++) {
            int newv;
            nc_def_var(outf->ncid, inf->varname[v], inf->datatype[v],
                       inf->varndims[v], inf->vardim[v], &newv);
            char attname[NC_MAX_NAME + 1];
            for (int n = 0; n < inf->natts[v]; n++) {
                nc_inq_attname(inf->ncid, v, n, attname);
                if (missing && strcmp(attname, "missing_value") == 0) {
                    outf->varmiss[v] = 1;
                    size_t esz = nc_type_sz(inf->datatype[v]);
                    nc_get_att(inf->ncid, v, "missing_value", outf->varmissval[v]);
                    (void)esz;
                }
                if (strcmp(attname, "domain_decomposition") == 0) continue;
                nc_copy_att(inf->ncid, v, attname, outf->ncid, newv);
            }
        }

        /* Copy global attributes (skip NumFilesInSet). */
        char attname[NC_MAX_NAME + 1];
        for (int n = 0; n < inf->ngatts; n++) {
            nc_inq_attname(inf->ncid, NC_GLOBAL, n, attname);
            if (strcmp(attname, "NumFilesInSet") == 0) continue;
            if (strcmp(attname, "filename") == 0) {
                nc_put_att_text(outf->ncid, NC_GLOBAL, attname,
                                strlen(outncname), outncname);
                continue;
            }
            nc_copy_att(inf->ncid, NC_GLOBAL, attname, outf->ncid, NC_GLOBAL);
        }

        if (deflate && deflation > 0)
            for (int v = 0; v < inf->nvars; v++)
                nc_def_var_deflate(outf->ncid, v, shuffle, deflate, deflation);

        nc__enddef(outf->ncid, (size_t)headerpad, 4, 0, 4);
    }

    /* Read / buffer variables. */
    int r = block * (*bf);
    int err = 0;
    do {
        err += mpp_process_vars(inf, outf, appendnc, nrecs, nblocks, bf,
                                r, *nfiles, f, verbose, missing);
        r++;
        appendnc = 1;
    } while (r < g_min((block + 1) * (*bf), *nrecs));

    nc_close(inf->ncid);
    free(inf);
    return err ? 1 : 0;
}

int
combine_mpp(combine_opts_t *opts)
{
    g_mem_dry_run    = opts->mem_dry_run;
    g_print_mem      = opts->print_mem;
    mem_allocated    = 0;
    estimated_maxrss = 0;
    varbuf           = NULL;

    int bf      = opts->blocking > 0 ? opts->blocking : DEFAULT_BF;
    int nstart  = opts->nstart;
    int nend    = opts->nend;
    int force   = opts->force;
    int verbose = opts->verbose;
    int missing = opts->missing;
    int removein = opts->removein;
    int headerpad = opts->headerpad > 0 ? opts->headerpad : 16384;
    int deflate = opts->deflate;
    int deflation = opts->deflation;
    int shuffle = opts->shuffle;

    int nfiles = -1, nrecs = 1, nblocks = 1;

    const char *outncname = opts->outfile;
    size_t blksz = NC_BLKSZ;

    fileinfo_t *outf = XMALLOC(fileinfo_t, 1);
    memset(outf, 0, sizeof(*outf));

    struct stat statbuf;
    int create_mode = NC_NOCLOBBER;
    if (opts->format_64) create_mode |= NC_64BIT_OFFSET;
    if (opts->format_n4) create_mode = NC_NOCLOBBER | NC_NETCDF4 | NC_CLASSIC_MODEL;

    if (!opts->appendnc) {
        if (stat(outncname, &statbuf) == 0) {
            fprintf(stderr,"Error: output file '%s' already exists.\n", outncname);
            free(outf); return 1;
        }
        if (nc__create(outncname, create_mode, 0, &blksz, &outf->ncid) != NC_NOERR) {
            fprintf(stderr,"Error: cannot create output '%s'.\n", outncname);
            free(outf); return 1;
        }
        nc_set_fill(outf->ncid, NC_NOFILL, NULL);
    } else {
        if (nc_open(outncname, NC_WRITE, &outf->ncid) != NC_NOERR) {
            fprintf(stderr,"Error: cannot open '%s' for appending.\n", outncname);
            free(outf); return 1;
        }
    }

    int appendnc = opts->appendnc;
    int infile_errors = 0, outfile_errors = 0;
    int peWidth = -1;
    char infilename[2048];
    (void)0;

    /* --- Explicit input files on command line --- */
    if (opts->ninfiles > 0) {
        for (int block = 0; block < nblocks; block++) {
            int f = 0;
            for (int a = 0; a < opts->ninfiles; a++) {
                if (stat(opts->infiles[a], &statbuf) != 0) {
                    if (!force) {
                        fprintf(stderr,"ERROR: missing file '%s'.\n", opts->infiles[a]);
                        unlink(outncname); free(outf); return 9;
                    }
                    infile_errors = 1;
                }
                if (verbose) printf("  processing '%s'\n", opts->infiles[a]);
                if (mpp_process_file(opts->infiles[a], appendnc, outf, outncname,
                                     &nfiles, &nrecs, &nblocks, &bf, block, f,
                                     headerpad, verbose, missing, deflate, deflation, shuffle))
                    infile_errors = 1;
                appendnc = 1; f++;
                if (f == nfiles || a == opts->ninfiles - 1) {
                    if (g_mem_dry_run) {
                        printf("%.0f\n", ceil((float)estimated_maxrss/(1024*1024)));
                        goto done;
                    }
                    for (int r = block * bf; r < g_min((block+1)*bf, nrecs); r++)
                        outfile_errors += mpp_flush_decomp(outf, nfiles, r, bf, verbose);
                    f = 0; continue;
                }
            }
            for (int k = 0; k < bf; k++)
                for (int v = 0; v < NC_MAX_VARS; v++)
                    if (varbuf && varbuf[k][v]) { free(varbuf[k][v]); varbuf[k][v] = NULL; }
        }
    } else {
        /* --- Auto-discover files by numeric extension --- */
        if (nend < 0) nend = nstart + 1;
        for (int block = 0; block < nblocks; block++) {
            int f = 0;
            for (int a = nstart; a <= nend; a++) {
                if (peWidth < 0) {
                    sprintf(infilename, "%s.%04d", outncname, a);
                    if (stat(infilename, &statbuf) == 0) peWidth = 4;
                    else {
                        sprintf(infilename, "%s.%06d", outncname, a);
                        if (stat(infilename, &statbuf) == 0) peWidth = 6;
                        else continue;
                    }
                }
                sprintf(infilename, "%s.%0*d", outncname, peWidth, a);
                if (stat(infilename, &statbuf) != 0) {
                    if (!force) {
                        fprintf(stderr,"ERROR: missing '%s'.\n", infilename);
                        unlink(outncname); free(outf); return 9;
                    }
                    infile_errors = 1;
                }
                if (mpp_process_file(infilename, appendnc, outf, outncname,
                                     &nfiles, &nrecs, &nblocks, &bf, block, f,
                                     headerpad, verbose, missing, deflate, deflation, shuffle))
                    infile_errors = 1;
                if (a == nstart && nfiles > 0 && opts->nend < 0)
                    nend = nstart + nfiles;
                appendnc = 1; f++;
                if (f == nfiles || a == nend) {
                    if (g_mem_dry_run) {
                        printf("%.0f\n", ceil((float)estimated_maxrss/(1024*1024)));
                        goto done;
                    }
                    for (int r = block*bf; r < g_min((block+1)*bf, nrecs); r++)
                        outfile_errors += mpp_flush_decomp(outf, nfiles, r, bf, verbose);
                    f = 0; break;
                }
            }
            for (int k = 0; k < bf; k++)
                for (int v = 0; v < NC_MAX_VARS; v++)
                    if (varbuf && varbuf[k][v]) { free(varbuf[k][v]); varbuf[k][v] = NULL; }
        }
    }

done:
    if (nc_sync(outf->ncid) != NC_NOERR) outfile_errors++;
    if (nc_close(outf->ncid) != NC_NOERR) outfile_errors++;
    free(outf);

    if (!infile_errors && !outfile_errors) {
        if (removein) {
            if (opts->ninfiles > 0) {
                for (int a = 0; a < opts->ninfiles; a++) {
                    if (stat(opts->infiles[a], &statbuf) == 0) unlink(opts->infiles[a]);
                }
            } else {
                for (int a = nstart; a <= nend; a++) {
                    sprintf(infilename, "%s.%0*d", outncname, peWidth > 0 ? peWidth : 4, a);
                    if (stat(infilename, &statbuf) == 0) unlink(infilename);
                }
            }
        }
        return 0;
    }
    fprintf(stderr,"Warning: output file may be incomplete.\n");
    return 1;
}

/* ================================================================== */
/* combine_land – port of combine-ncc.F90                             */
/* ================================================================== */
int
combine_land(combine_opts_t *opts)
{
    if (opts->ninfiles < 1) {
        fprintf(stderr,"combine --land: no input files specified.\n");
        return 1;
    }

    const char *outfile = opts->outfile;
    int nfiles = opts->ninfiles;
    int *input = XMALLOC(int, nfiles);
    int ncid;

    /* Open all input files. */
    int in_format = NC_FORMAT_CLASSIC;
    for (int i = 0; i < nfiles; i++) {
        if (nc_open(opts->infiles[i], NC_NOWRITE, &input[i]) != NC_NOERR) {
            fprintf(stderr,"combine --land: cannot open '%s'.\n", opts->infiles[i]);
            free(input); return 1;
        }
        nc_inq_format(input[i], &in_format);
    }

    /* Create output file with matching format. */
    int cmode;
    switch (in_format) {
    case NC_FORMAT_NETCDF4:          cmode = NC_NETCDF4; break;
    case NC_FORMAT_NETCDF4_CLASSIC:  cmode = NC_NETCDF4 | NC_CLASSIC_MODEL; break;
    case NC_FORMAT_64BIT_OFFSET:     cmode = NC_CLOBBER | NC_64BIT_OFFSET; break;
    default:                         cmode = NC_CLOBBER | NC_64BIT_OFFSET; break;
    }

    size_t blksz = 65536;
    if (nc__create(outfile, cmode, 0, &blksz, &ncid) != NC_NOERR) {
        fprintf(stderr,"combine --land: cannot create '%s'.\n", outfile);
        free(input); return 1;
    }

    int ref = input[nfiles - 1]; /* Use last file as template. */

    /* --- Define dimensions --- */
    int ndims;
    nc_inq_ndims(ref, &ndims);

    /* For each dimension determine output length. */
    int *dim_is_compressed  = XMALLOC(int, ndims);
    int *dim_out_len        = XMALLOC(int, ndims);

    for (int dimid = 0; dimid < ndims; dimid++) {
        char dname[NC_MAX_NAME + 1];
        size_t dlen;
        bool is_unlim;
        ncu_inq_dim(ref, dimid, dname, &dlen, &is_unlim);

        int varid, attid;
        dim_is_compressed[dimid] = 0;
        dim_out_len[dimid] = (int)dlen;

        if (nc_inq_varid(ref, dname, &varid) == NC_NOERR &&
            nc_inq_attid(ref, varid, "compress", &attid) == NC_NOERR) {
            dim_is_compressed[dimid] = 1;
            /* Count valid (non-negative) indices across all files. */
            int total_len = 0;
            for (int fi = 0; fi < nfiles; fi++) {
                size_t slen;
                if (ncu_inq_dim(input[fi], dimid, NULL, &slen, NULL) != NC_NOERR) continue;
                int *buf = XMALLOC(int, slen);
                if (ncu_get_var_int(input[fi], dname, buf) == NC_NOERR)
                    for (size_t k = 0; k < slen; k++)
                        if (buf[k] >= 0) total_len++;
                free(buf);
            }
            dim_out_len[dimid] = total_len > 0 ? total_len : 1;
        }

        int dummy;
        nc_def_dim(ncid, dname, is_unlim ? NC_UNLIMITED : (size_t)dim_out_len[dimid], &dummy);
    }

    /* --- Define variables --- */
    int nvars;
    nc_inq_nvars(ref, &nvars);
    for (int v = 0; v < nvars; v++)
        ncu_clone_var(ref, v, ncid);

    /* --- Copy global attributes --- */
    ncu_copy_global_atts(ref, ncid);

    nc__enddef(ncid, 16384, 4, 0, 4);

    /* --- Copy uncompressed variables --- */
    for (int v = 0; v < nvars; v++) {
        char vname[NC_MAX_NAME + 1];
        nc_type xtype;
        int vndims, dimids[NC_MAX_DIMS];
        size_t dimlens[NC_MAX_DIMS];
        int nrec, recsize;
        bool has_records;

        ncu_inq_var(ref, v, vname, &xtype, &vndims, dimids, dimlens, NULL,
                    NULL, &has_records, NULL, (size_t*)&recsize, &nrec);

        /* Does this variable have any compressed dimension? */
        int has_compressed = 0;
        for (int d = 0; d < vndims; d++)
            if (dim_is_compressed[dimids[d]]) { has_compressed++; break; }

        if (has_compressed) continue; /* handled later */

        /* Copy uncompressed variable. */
        if (xtype == NC_CHAR) {
            size_t total = (size_t)recsize * (size_t)nrec;
            char *text = XMALLOC(char, total);
            nc_get_var_text(ref, v, text);
            int ov; nc_inq_varid(ncid, vname, &ov);
            nc_put_var_text(ncid, ov, text);
            free(text);
        } else {
            double *buf = XMALLOC(double, (size_t)recsize);
            for (int rec = 0; rec < nrec; rec++) {
                ncu_get_rec(ref, v, rec, buf);
                int ov; nc_inq_varid(ncid, vname, &ov);
                ncu_put_rec(ncid, ov, rec, buf);
            }
            free(buf);
        }
    }

    /* --- Copy compressed variables --- */
    for (int dimid = 0; dimid < ndims; dimid++) {
        if (!dim_is_compressed[dimid]) continue;

        char dname[NC_MAX_NAME + 1];
        nc_inq_dimname(ref, dimid, dname);

        /* Get per-file compressed dimension sizes. */
        int *sizes = XMALLOC(int, nfiles);
        int buflen = 0;
        for (int fi = 0; fi < nfiles; fi++) {
            size_t slen;
            if (ncu_inq_dim(input[fi], dimid, NULL, &slen, NULL) == NC_NOERR)
                sizes[fi] = (int)slen;
            else
                sizes[fi] = 0;
            buflen += sizes[fi];
        }

        if (dim_out_len[dimid] == 0) { free(sizes); continue; }

        /* Build reordering index. */
        int *rank_buf = XMALLOC(int, buflen);
        int k = 0;
        for (int fi = 0; fi < nfiles; fi++) {
            if (sizes[fi] == 0) continue;
            if (ncu_get_var_int(input[fi], dname, rank_buf + k) != NC_NOERR)
                memset(rank_buf + k, -1, sizes[fi] * sizeof(int));
            k += sizes[fi];
        }

        int *rank = XMALLOC(int, buflen);
        rank_ascending(rank_buf, rank, buflen);

        /* Skip leading negatives. */
        int skip = 0;
        while (skip < buflen && rank_buf[rank[skip]] < 0) skip++;
        for (int i = 0; i < buflen - skip; i++) rank[i] = rank[i + skip];
        for (int i = buflen - skip; i < buflen; i++) rank[i] = -1;

        /* Process variables depending on this compressed dimension. */
        for (int v = 0; v < nvars; v++) {
            char vname[NC_MAX_NAME + 1];
            nc_type xtype;
            int vndims, dimids[NC_MAX_DIMS];
            size_t dimlens[NC_MAX_DIMS];
            int nrec;
            bool has_records;
            ncu_inq_var(ref, v, vname, &xtype, &vndims, dimids, dimlens,
                        NULL, NULL, &has_records, NULL, NULL, &nrec);

            int cdim = -1;
            for (int d = 0; d < vndims; d++)
                if (dimids[d] == dimid) { cdim = d; break; }
            if (cdim < 0) continue;

            int ovarid;
            if (nc_inq_varid(ncid, vname, &ovarid) != NC_NOERR) continue;

            /* Product of non-compressed dimensions. */
            size_t slice_count = 1;
            for (int d = 0; d < vndims; d++)
                if (d != cdim) slice_count *= dimlens[d];

            double *ibuf  = XMALLOC(double, (size_t)buflen);
            double *obuf  = XMALLOC(double, (size_t)dim_out_len[dimid]);
            size_t start[NC_MAX_DIMS], cnt[NC_MAX_DIMS];

            for (size_t si = 0; si < slice_count; si++) {
                /* Determine start indices for this slice. */
                size_t ii = si;
                for (int d = 0; d < vndims; d++) {
                    if (d == cdim) { start[d] = 0; }
                    else { start[d] = ii % dimlens[d]; ii /= dimlens[d]; }
                }
                for (int d = 0; d < vndims; d++) cnt[d] = 1;

                /* Read from all input files. */
                k = 0;
                for (int fi = 0; fi < nfiles; fi++) {
                    if (sizes[fi] == 0) continue;
                    int ivarid;
                    if (nc_inq_varid(input[fi], vname, &ivarid) != NC_NOERR) {
                        k += sizes[fi]; continue;
                    }
                    cnt[cdim] = (size_t)sizes[fi];
                    nc_get_vara_double(input[fi], ivarid, start, cnt, ibuf + k);
                    k += sizes[fi];
                }

                /* Reorder. */
                for (int i = 0; i < dim_out_len[dimid]; i++)
                    obuf[i] = ibuf[rank[i]];

                /* Write to output. */
                cnt[cdim] = (size_t)dim_out_len[dimid];
                nc_put_vara_double(ncid, ovarid, start, cnt, obuf);
            }
            free(ibuf); free(obuf);
        }
        free(rank_buf); free(rank); free(sizes);
    }

    nc_sync(ncid); nc_close(ncid);
    for (int i = 0; i < nfiles; i++) nc_close(input[i]);
    free(input); free(dim_is_compressed); free(dim_out_len);

    if (opts->removein)
        for (int i = 0; i < nfiles; i++) unlink(opts->infiles[i]);
    return 0;
}

/* ================================================================== */
/* combine_iceberg – port of iceberg_comb.sh (no ncrcat/ncatted)      */
/* ================================================================== */
int
combine_iceberg(combine_opts_t *opts)
{
    if (opts->ninfiles < 2) {
        fprintf(stderr,"combine --iceberg: need at least 2 input files.\n");
        return 1;
    }
    const char *outfile = opts->outfile;

    /* Validate: all files must be iceberg restart files. */
    for (int i = 0; i < opts->ninfiles; i++) {
        if (!check_is_iceberg(opts->infiles[i])) {
            fprintf(stderr,"combine --iceberg: '%s' is not an iceberg restart.\n",
                    opts->infiles[i]);
            return 1;
        }
    }

    /* Filter to files that actually have icebergs (non-empty UNLIMITED dim). */
    char **combine_files = XMALLOC(char *, opts->ninfiles);
    int ncombine = 0;
    for (int i = 0; i < opts->ninfiles; i++)
        if (check_iceberg_has_data(opts->infiles[i]))
            combine_files[ncombine++] = opts->infiles[i];

    if (ncombine == 0) {
        /* No icebergs: copy first input as output. */
        FILE *src = fopen(opts->infiles[0], "rb");
        FILE *dst = fopen(outfile, "wb");
        if (!src || !dst) {
            fprintf(stderr,"combine --iceberg: cannot copy file.\n");
            if (src) fclose(src); if (dst) fclose(dst);
            free(combine_files); return 1;
        }
        char buf[65536]; size_t n;
        while ((n = fread(buf, 1, sizeof(buf), src)) > 0)
            fwrite(buf, 1, n, dst);
        fclose(src); fclose(dst);
        /* Remove NumFilesInSet attribute. */
        int ncid;
        if (nc_open(outfile, NC_WRITE, &ncid) == NC_NOERR) {
            nc_del_att(ncid, NC_GLOBAL, "NumFilesInSet");
            nc_close(ncid);
        }
        free(combine_files);
        return 0;
    }

    /* Open first non-empty file as template. */
    int ref_ncid;
    if (nc_open(combine_files[0], NC_NOWRITE, &ref_ncid) != NC_NOERR) {
        fprintf(stderr,"combine --iceberg: cannot open '%s'.\n", combine_files[0]);
        free(combine_files); return 1;
    }

    /* Determine output format. */
    int in_fmt = NC_FORMAT_CLASSIC;
    nc_inq_format(ref_ncid, &in_fmt);
    int cmode;
    switch (in_fmt) {
    case NC_FORMAT_NETCDF4:         cmode = NC_NETCDF4; break;
    case NC_FORMAT_NETCDF4_CLASSIC: cmode = NC_NETCDF4 | NC_CLASSIC_MODEL; break;
    case NC_FORMAT_64BIT_OFFSET:    cmode = NC_CLOBBER | NC_64BIT_OFFSET; break;
    default:                        cmode = NC_CLOBBER | NC_64BIT_OFFSET; break;
    }

    int out_ncid;
    size_t blksz = 65536;
    if (nc__create(outfile, cmode, 0, &blksz, &out_ncid) != NC_NOERR) {
        fprintf(stderr,"combine --iceberg: cannot create '%s'.\n", outfile);
        nc_close(ref_ncid); free(combine_files); return 1;
    }
    nc_set_fill(out_ncid, NC_NOFILL, NULL);

    /* Query unlimited dimension name and total combined size. */
    int unlimdimid;
    char unlimited_name[NC_MAX_NAME + 1] = "";
    /* total_recs not needed; NC_UNLIMITED handles growth automatically. */

    if (nc_inq_unlimdim(ref_ncid, &unlimdimid) == NC_NOERR && unlimdimid >= 0)
        nc_inq_dimname(ref_ncid, unlimdimid, unlimited_name);

    for (int i = 0; i < ncombine; i++) {
        int tmpid;
        if (nc_open(combine_files[i], NC_NOWRITE, &tmpid) == NC_NOERR) {
            int uid;
            if (nc_inq_unlimdim(tmpid, &uid) == NC_NOERR && uid >= 0) {
                size_t dlen; nc_inq_dimlen(tmpid, uid, &dlen); (void)dlen;
            }
            nc_close(tmpid);
        }
    }

    /* Define dimensions in output. */
    int ndims;
    nc_inq_ndims(ref_ncid, &ndims);
    for (int d = 0; d < ndims; d++) {
        char dname[NC_MAX_NAME + 1]; size_t dlen;
        nc_inq_dim(ref_ncid, d, dname, &dlen);
        int dummy;
        if (d == unlimdimid)
            nc_def_dim(out_ncid, dname, NC_UNLIMITED, &dummy);
        else
            nc_def_dim(out_ncid, dname, dlen, &dummy);
    }

    /* Define variables. */
    int nvars;
    nc_inq_nvars(ref_ncid, &nvars);
    for (int v = 0; v < nvars; v++) ncu_clone_var(ref_ncid, v, out_ncid);

    /* Copy global attributes, skip NumFilesInSet. */
    int ngatts;
    nc_inq_natts(ref_ncid, &ngatts);
    for (int i = 0; i < ngatts; i++) {
        char aname[NC_MAX_NAME + 1];
        nc_inq_attname(ref_ncid, NC_GLOBAL, i, aname);
        if (strcmp(aname, "NumFilesInSet") == 0) continue;
        nc_copy_att(ref_ncid, NC_GLOBAL, aname, out_ncid, NC_GLOBAL);
    }

    nc_enddef(out_ncid);

    /* Copy non-record variables from first file. */
    for (int v = 0; v < nvars; v++) {
        char vname[NC_MAX_NAME + 1];
        nc_type xtype;
        int vndims, dimids[NC_MAX_DIMS], natts;
        nc_inq_var(ref_ncid, v, vname, &xtype, &vndims, dimids, &natts);

        bool is_rec = false;
        for (int d = 0; d < vndims; d++)
            if (dimids[d] == unlimdimid) { is_rec = true; break; }
        if (is_rec) continue;

        size_t total = 1;
        for (int d = 0; d < vndims; d++) { size_t dl; nc_inq_dimlen(ref_ncid,dimids[d],&dl); total*=dl; }
        size_t esz = ncu_type_size(xtype);
        void *buf = malloc(total * esz);
        ncu_get_vara(ref_ncid, v, NULL, NULL, buf);
        int ov; nc_inq_varid(out_ncid, vname, &ov);
        /* Build start/count for put. */
        size_t st[NC_MAX_DIMS], cnt[NC_MAX_DIMS];
        for (int d = 0; d < vndims; d++) { st[d]=0; nc_inq_dimlen(ref_ncid,dimids[d],&cnt[d]); }
        ncu_put_vara(out_ncid, ov, st, cnt, buf);
        free(buf);
    }

    /* Concatenate record variables along UNLIMITED dim. */
    size_t out_rec_offset = 0;
    for (int fi = 0; fi < ncombine; fi++) {
        int in_ncid;
        if (nc_open(combine_files[fi], NC_NOWRITE, &in_ncid) != NC_NOERR) continue;
        int uid; size_t nrecs_in = 0;
        if (nc_inq_unlimdim(in_ncid, &uid) == NC_NOERR && uid >= 0)
            nc_inq_dimlen(in_ncid, uid, &nrecs_in);

        for (size_t rec = 0; rec < nrecs_in; rec++) {
            for (int v = 0; v < nvars; v++) {
                char vname[NC_MAX_NAME + 1];
                nc_type xtype;
                int vndims, dimids[NC_MAX_DIMS], natts;
                nc_inq_var(ref_ncid, v, vname, &xtype, &vndims, dimids, &natts);

                int rec_dim = -1;
                for (int d = 0; d < vndims; d++)
                    if (dimids[d] == unlimdimid) { rec_dim = d; break; }
                if (rec_dim < 0) continue;

                /* Get var from input file. */
                int iv;
                if (nc_inq_varid(in_ncid, vname, &iv) != NC_NOERR) continue;

                size_t st_in[NC_MAX_DIMS], cnt[NC_MAX_DIMS];
                size_t st_out[NC_MAX_DIMS];
                for (int d = 0; d < vndims; d++) {
                    if (d == rec_dim) {
                        st_in[d] = rec; cnt[d] = 1;
                        st_out[d] = out_rec_offset + rec;
                    } else {
                        st_in[d] = 0; st_out[d] = 0;
                        nc_inq_dimlen(in_ncid, dimids[d], &cnt[d]);
                    }
                }
                size_t total = 1;
                for (int d = 0; d < vndims; d++) total *= cnt[d];
                size_t esz = ncu_type_size(xtype);
                void *buf = malloc(total * esz);
                ncu_get_vara(in_ncid, iv, st_in, cnt, buf);

                int ov; nc_inq_varid(out_ncid, vname, &ov);
                ncu_put_vara(out_ncid, ov, st_out, cnt, buf);
                free(buf);
            }
        }
        out_rec_offset += nrecs_in;
        nc_close(in_ncid);
    }
    nc_close(ref_ncid);

    /* Remove NumFilesInSet (already excluded above, but ensure it's absent). */
    nc_redef(out_ncid);
    nc_del_att(out_ncid, NC_GLOBAL, "NumFilesInSet");
    nc_enddef(out_ncid);

    nc_sync(out_ncid);
    nc_close(out_ncid);
    free(combine_files);

    if (opts->removein)
        for (int i = 0; i < opts->ninfiles; i++) unlink(opts->infiles[i]);
    return 0;
}

/* ================================================================== */
/* cmd_combine – argument parsing + dispatch                           */
/* ================================================================== */
static void usage_combine(void) {
    printf("Usage: mppdisttool combine [-o outfile] [-v] [-r] [--land|--iceberg|--mpp]\n"
           "                           [-h N] [-m] [-k N] [-64|-n4] [-a] [-f] infile...\n");
}

int
cmd_combine(int argc, char **argv)
{
    combine_opts_t opts;
    memset(&opts, 0, sizeof(opts));
    opts.headerpad = 16384;
    opts.blocking  = DEFAULT_BF;
    opts.nend      = -1;

    char **infiles = XMALLOC(char *, argc);
    int ninfiles = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i],"--land") == 0)        opts.land_mode    = 1;
        else if (strcmp(argv[i],"--iceberg") == 0) opts.iceberg_mode = 1;
        else if (strcmp(argv[i],"--mpp") == 0)     opts.mpp_mode     = 1;
        else if (strcmp(argv[i],"-v") == 0)         opts.verbose++;
        else if (strcmp(argv[i],"-r") == 0)         opts.removein = 1;
        else if (strcmp(argv[i],"-a") == 0)         opts.appendnc = 1;
        else if (strcmp(argv[i],"-f") == 0)         opts.force    = 1;
        else if (strcmp(argv[i],"-m") == 0)         opts.missing  = 1;
        else if (strcmp(argv[i],"-64") == 0)        opts.format_64 = 1;
        else if (strcmp(argv[i],"-n4") == 0)        opts.format_n4 = 1;
        else if (strcmp(argv[i],"-M") == 0)         opts.print_mem = 1;
        else if (strcmp(argv[i],"-x") == 0)         opts.mem_dry_run = 1;
        else if ((strcmp(argv[i],"-o") == 0) && i+1 < argc) {
            opts.outfile = argv[++i];
        } else if ((strcmp(argv[i],"-h") == 0) && i+1 < argc) {
            opts.headerpad = atoi(argv[++i]);
        } else if ((strcmp(argv[i],"-k") == 0) && i+1 < argc) {
            opts.blocking = atoi(argv[++i]);
        } else if ((strcmp(argv[i],"-n") == 0) && i+1 < argc) {
            opts.nstart = atoi(argv[++i]);
        } else if ((strcmp(argv[i],"-e") == 0) && i+1 < argc) {
            opts.nend = atoi(argv[++i]);
        } else if ((strcmp(argv[i],"-d") == 0) && i+1 < argc) {
            opts.deflate = 1; opts.deflation = atoi(argv[++i]);
        } else if (strcmp(argv[i],"-s") == 0) {
            opts.shuffle = 1;
        } else if (strcmp(argv[i],"--help") == 0 || strcmp(argv[i],"-?") == 0) {
            usage_combine(); free(infiles); return 0;
        } else if (argv[i][0] != '-') {
            infiles[ninfiles++] = argv[i];
        }
    }

    /* First non-option arg is output if --mpp mode with no -o. */
    if (!opts.outfile && ninfiles > 0 && (opts.mpp_mode || (!opts.land_mode && !opts.iceberg_mode))) {
        opts.outfile = infiles[0];
        opts.infiles = infiles + 1;
        opts.ninfiles = ninfiles - 1;
    } else {
        opts.infiles  = infiles;
        opts.ninfiles = ninfiles;
    }

    /* If no mode explicitly given, auto-detect from first input file. */
    if (!opts.land_mode && !opts.iceberg_mode && !opts.mpp_mode) {
        if (opts.ninfiles > 0) {
            const char *probe = opts.infiles[0];
            int probeid;
            if (nc_open(probe, NC_NOWRITE, &probeid) == NC_NOERR) {
                if (file_is_compressed(probeid)) opts.land_mode = 1;
                else if (check_is_iceberg(probe)) opts.iceberg_mode = 1;
                else opts.mpp_mode = 1;
                nc_close(probeid);
            } else opts.mpp_mode = 1;
        } else opts.mpp_mode = 1;
    }

    int ret;
    if (opts.land_mode)        ret = combine_land(&opts);
    else if (opts.iceberg_mode) ret = combine_iceberg(&opts);
    else                        ret = combine_mpp(&opts);

    free(infiles);
    return ret;
}

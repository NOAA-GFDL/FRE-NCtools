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
 * Two scatter paths:
 *   scatter_mpp  – domain decomposition (from mppncscatter.c + domain.c)
 *   scatter_land – CF compressed-by-gathering scatter (from scatter-ncc.F90)
 */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <netcdf.h>
#include "cmd_scatter.h"
#include "domain.h"
#include "nc_utils.h"
#include "compress.h"
#include "xmalloc.h"
#include "strlist.h"

/* ================================================================== */
/* scatter_mpp – wraps domain.c logic (ported from mppncscatter)      */
/* ================================================================== */
static int
scatter_mpp(scatter_opts_t *opts)
{
    int nc;
    if (nc_open(opts->filein, NC_NOWRITE, &nc) != NC_NOERR) {
        fprintf(stderr,"scatter: cannot open '%s'\n", opts->filein);
        return -1;
    }

    int ndims, nvars, ngatts, unlimdimid;
    nc_inq(nc, &ndims, &nvars, &ngatts, &unlimdimid);

    int format_flag = NC_FORMAT_CLASSIC;
    nc_inq_format(nc, &format_flag);
    int create_flags;
    switch (format_flag) {
    case NC_FORMAT_64BIT_OFFSET:    create_flags = NC_64BIT_OFFSET;  break;
    case NC_FORMAT_NETCDF4:         create_flags = NC_NETCDF4;        break;
    case NC_FORMAT_NETCDF4_CLASSIC: create_flags = NC_NETCDF4 | NC_CLASSIC_MODEL; break;
    default:                        create_flags = 0;                 break;
    }

    /* Strip directory prefix from filein for output names. */
    const char *prefix = opts->filein;
    const char *p = strrchr(prefix, '/');
    if (p) prefix = p + 1;

    /* Build output filename format string. */
    char outnameformat[512];
    if (opts->prefix && strlen(opts->prefix) > 0)
        snprintf(outnameformat, sizeof(outnameformat), "%s/%%s.%%0%dd",
                 opts->prefix, opts->width);
    else
        snprintf(outnameformat, sizeof(outnameformat), "%%s.%%0%dd", opts->width);

    ScatterDim *scatterdims[NC_MAX_DIMS];
    scatter_dims(nc, ndims, nvars, scatterdims, opts);

    int nfiles = scatter_get_num_files(opts);
    int *ncids = XMALLOC(int, nfiles);

    for (int i = 0; i < nfiles; i++) {
        char output[512];
        snprintf(output, sizeof(output), outnameformat, prefix, i + opts->start);

        if (opts->verbose) printf("Info: Creating '%s'\n", output);

        if (!opts->dryrun) {
            int nc_err = nc_create(output, NC_CLOBBER | create_flags, &ncids[i]);
            if (nc_err != NC_NOERR) {
                fprintf(stderr, "scatter: cannot create '%s': %s\n",
                        output, nc_strerror(nc_err));
                /* Clean up already-created files */
                for (int j = 0; j < i; j++) {
                    nc_close(ncids[j]);
                }
                free(ncids);
                nc_close(nc);
                return -1;
            }
            int dummy;
            nc_set_fill(ncids[i], NC_NOFILL, &dummy);
            nc_put_att_int(ncids[i], NC_GLOBAL, "NumFilesInSet", NC_INT, 1, &nfiles);
        }
    }

    /* Copy global attributes. */
    if (!opts->dryrun) {
        char attname[NC_MAX_NAME + 1];
        for (int j = 0; j < ngatts; j++) {
            nc_inq_attname(nc, NC_GLOBAL, j, attname);
            for (int i = 0; i < nfiles; i++)
                nc_copy_att(nc, NC_GLOBAL, attname, ncids[i], NC_GLOBAL);
        }
    }

    scatter_def_dim(nc, ncids, ndims, scatterdims, opts);
    scatter_def_var(nc, ncids, nvars, ndims, scatterdims, opts);

    if (!opts->dryrun)
        for (int i = 0; i < nfiles; i++) nc_enddef(ncids[i]);

    scatter_put_var(nc, ncids, ndims, nvars, scatterdims, opts);

    nc_close(nc);
    if (!opts->dryrun)
        for (int i = 0; i < nfiles; i++) { nc_sync(ncids[i]); nc_close(ncids[i]); }

    scatter_dims_free(scatterdims, ndims);
    free(ncids);
    return 0;
}

/* ================================================================== */
/* scatter_land – port of scatter-ncc.F90                             */
/* ================================================================== */
static int
scatter_land(const char *infile, int npex, int npey)
{
    int in_ncid;
    if (nc_open(infile, NC_NOWRITE, &in_ncid) != NC_NOERR) {
        fprintf(stderr,"scatter --land: cannot open '%s'\n", infile);
        return 1;
    }

    int in_fmt;
    nc_inq_format(in_ncid, &in_fmt);
    int cmode;
    switch (in_fmt) {
    case NC_FORMAT_NETCDF4:         cmode = NC_NETCDF4; break;
    case NC_FORMAT_NETCDF4_CLASSIC: cmode = NC_NETCDF4 | NC_CLASSIC_MODEL; break;
    case NC_FORMAT_64BIT_OFFSET:    cmode = NC_CLOBBER | NC_64BIT_OFFSET; break;
    default:                        cmode = NC_CLOBBER | NC_64BIT_OFFSET; break;
    }

    int nfiles_out = npex * npey;

    /* Find lon/lat/zfull dimensions in input. */
    int ndims;
    nc_inq_ndims(in_ncid, &ndims);

    int nlon = 0, nlat = 0, nz = 1;
    for (int d = 0; d < ndims; d++) {
        char dname[NC_MAX_NAME + 1]; size_t dlen;
        nc_inq_dim(in_ncid, d, dname, &dlen);
        if (strcmp(dname,"lon")   == 0) nlon = (int)dlen;
        if (strcmp(dname,"lat")   == 0) nlat = (int)dlen;
        if (strcmp(dname,"zfull") == 0) nz   = (int)dlen;
    }

    if (nlon == 0 || nlat == 0) {
        fprintf(stderr,"scatter --land: input lacks 'lon' or 'lat' dimension.\n");
        nc_close(in_ncid); return 1;
    }
    if (nlon % npex != 0 || nlat % npey != 0) {
        fprintf(stderr,"scatter --land: domain not evenly divisible by npex/npey.\n");
        nc_close(in_ncid); return 1;
    }
    int nlon_local = nlon / npex;
    int nlat_local = nlat / npey;

    /* Create output files. */
    int *out_ncids = XMALLOC(int, nfiles_out);
    int *dimlen_list = XMALLOC(int, nfiles_out);
    size_t blksz = 65536;
    for (int n = 0; n < nfiles_out; n++) {
        char outname[2048];
        snprintf(outname, sizeof(outname), "%s.%04d", infile, n);
        int nc_err = nc__create(outname, cmode, 0, &blksz, &out_ncids[n]);
        if (nc_err != NC_NOERR) {
            fprintf(stderr, "scatter --land: cannot create '%s': %s\n",
                    outname, nc_strerror(nc_err));
            /* Clean up already-created files */
            for (int j = 0; j < n; j++) {
                nc_close(out_ncids[j]);
            }
            free(out_ncids);
            free(dimlen_list);
            nc_close(in_ncid);
            return 1;
        }
    }

    /* Count compressed dimension sizes per output tile. */

    for (int d = 0; d < ndims; d++) {
        char dname[NC_MAX_NAME + 1]; size_t dlen; bool is_unlim;
        ncu_inq_dim(in_ncid, d, dname, &dlen, &is_unlim);

        int is_compressed = 0;
        int varid, attid;
        if (nc_inq_varid(in_ncid, dname, &varid) == NC_NOERR &&
            nc_inq_attid(in_ncid, varid, "compress", &attid) == NC_NOERR)
            is_compressed = 1;

        memset(dimlen_list, 0, nfiles_out * sizeof(int));

        if (is_compressed) {
            /* Use compress_get_var_double to get full uncompressed data with mask. */
            int cndims; size_t cdimlens[NC_MAX_DIMS];
            compress_inq_dim(in_ncid, d, &cndims, NULL, cdimlens);
            size_t full_size = 1;
            for (int k = 0; k < cndims; k++) full_size *= cdimlens[k];

            double *buf  = XMALLOC(double, full_size);
            bool   *mask = XMALLOC(bool,   full_size);
            memset(mask, 0, full_size * sizeof(bool));
            compress_get_var_double(in_ncid, varid, buf, mask, full_size);

            int npts = nlon * nlat;
            for (size_t l = 0; l < full_size; l++) {
                if (!mask[l]) continue;
                int ij = (int)(l % (size_t)npts);
                int i  = ij % nlon;
                int j  = ij / nlon;
                int ii = i / nlon_local;
                int jj = j / nlat_local;
                int nn = jj * npex + ii;
                if (nn < nfiles_out) dimlen_list[nn]++;
            }
            free(buf); free(mask);
        }

        size_t out_dlen = (size_t)dlen;
        for (int n = 0; n < nfiles_out; n++) {
            int dummy;
            if (is_compressed) out_dlen = (size_t)(dimlen_list[n] > 0 ? dimlen_list[n] : 1);
            nc_def_dim(out_ncids[n], dname, is_unlim ? NC_UNLIMITED : out_dlen, &dummy);
        }
    }

    /* Clone variable definitions and global attributes. */
    int nvars;
    nc_inq_nvars(in_ncid, &nvars);
    for (int v = 0; v < nvars; v++)
        for (int n = 0; n < nfiles_out; n++)
            ncu_clone_var(in_ncid, v, out_ncids[n]);

    int ngatts;
    nc_inq_natts(in_ncid, &ngatts);
    for (int i = 0; i < ngatts; i++) {
        char aname[NC_MAX_NAME + 1];
        nc_inq_attname(in_ncid, NC_GLOBAL, i, aname);
        for (int n = 0; n < nfiles_out; n++)
            nc_copy_att(in_ncid, NC_GLOBAL, aname, out_ncids[n], NC_GLOBAL);
    }

    for (int n = 0; n < nfiles_out; n++)
        nc__enddef(out_ncids[n], 16384, 4, 0, 4);

    /* Query number of records. */
    int unlimdimid, nrec = 1;
    if (nc_inq_unlimdim(in_ncid, &unlimdimid) == NC_NOERR && unlimdimid >= 0) {
        size_t nrec_sz;
        nc_inq_dimlen(in_ncid, unlimdimid, &nrec_sz);
        nrec = (int)nrec_sz;
    }

    int npts      = nlon * nlat;
    int npts_local = nlon_local * nlat_local;

    /* For each variable and each record, scatter data. */
    for (int tlev = 0; tlev < nrec; tlev++) {
        for (int v = 0; v < nvars; v++) {
            char vname[NC_MAX_NAME + 1];
            nc_type xtype;
            int vndims, dimids[NC_MAX_DIMS];
            size_t dimlens[NC_MAX_DIMS];
            bool has_records, is_compressed_var;
            ncu_inq_var(in_ncid, v, vname, &xtype, &vndims, dimids, dimlens,
                        NULL, NULL, &has_records, NULL, NULL, NULL);

            if (!has_records && tlev > 0) continue;

            /* Check if first dim is compressed. */
            is_compressed_var = false;
            if (vndims > 0) {
                char d0name[NC_MAX_NAME+1];
                int varid0, attid0;
                nc_inq_dimname(in_ncid, dimids[0], d0name);
                if (nc_inq_varid(in_ncid,d0name,&varid0)==NC_NOERR &&
                    nc_inq_attid(in_ncid,varid0,"compress",&attid0)==NC_NOERR)
                    is_compressed_var = true;
            }

            /* Get uncompressed variable size. */
            int cndims_v; size_t cdimlens_v[NC_MAX_DIMS];
            compress_inq_dim(in_ncid, dimids[0], &cndims_v, NULL, cdimlens_v);

            size_t vsize;
            if (is_compressed_var) {
                /* Use first compressed dim's uncompressed dims + rest. */
                vsize = 1;
                for (int k = 0; k < cndims_v; k++) vsize *= cdimlens_v[k];
                for (int d = 1; d < vndims; d++) if (dimids[d] != unlimdimid) vsize *= dimlens[d];
            } else {
                vsize = 1;
                for (int d = 0; d < vndims; d++) if (dimids[d] != unlimdimid) vsize *= dimlens[d];
            }

            double *buf  = XMALLOC(double, vsize > 0 ? vsize : 1);
            bool   *mask = XMALLOC(bool,   vsize > 0 ? vsize : 1);
            memset(mask, 0, vsize * sizeof(bool));

            /* Read variable from input. */
            if (is_compressed_var) {
                /* Determine level dimension. */
                int k_dim = -1;
                for (int d = 0; d < vndims; d++) {
                    char ddn[NC_MAX_NAME+1];
                    nc_inq_dimname(in_ncid, dimids[d], ddn);
                    if (strcmp(ddn,"zfull")==0) k_dim = d;
                }
                int nz_var = (k_dim >= 0) ? nz : 1;
                (void)nz_var; /* recsize not used in this path */

                size_t st[NC_MAX_DIMS], cnt[NC_MAX_DIMS];
                st[0] = (size_t)tlev; cnt[0] = 1;
                for (int d = 1; d < vndims; d++) { st[d]=1; cnt[d]=1; }
                /* Read all at once. */
                compress_get_var_double(in_ncid, v, buf, mask, vsize);
            } else {
                if (has_records) {
                    /* 1D or scalar record variable. */
                    size_t st1[1] = {(size_t)tlev}, cnt1[1] = {1};
                    nc_get_vara_double(in_ncid, v, st1, cnt1, buf);
                } else {
                    nc_get_var_double(in_ncid, v, buf);
                }
            }

            /* Scatter to output tiles. */
            for (int n = 0; n < nfiles_out; n++) {
                int ii_tile = n % npex;
                int jj_tile = n / npex;

                if (is_compressed_var) {
                    /* Collect points belonging to this tile. */
                    double *tbuf = XMALLOC(double, (size_t)npts_local);
                    bool   *tmask = XMALLOC(bool, (size_t)npts_local);
                    memset(tbuf, 0, npts_local * sizeof(double));
                    memset(tmask, 0, npts_local * sizeof(bool));
                    int tcount = 0;

                    for (int l = 0; l < (int)vsize; l++) {
                        if (!mask[l]) continue;
                        int ij = l % npts;
                        int gi = ij % nlon;
                        int gj = ij / nlon;
                        int li_tile = gi / nlon_local;
                        int lj_tile = gj / nlat_local;
                        if (li_tile != ii_tile || lj_tile != jj_tile) continue;
                        int li_local = gi % nlon_local;
                        int lj_local = gj % nlat_local;
                        int ll = lj_local * nlon_local + li_local;
                        if (ll < npts_local) { tbuf[ll] = buf[l]; tmask[ll] = true; tcount++; }
                    }
                    if (tcount > 0) {
                        /* Write packed data to output tile. */
                        double *packed = XMALLOC(double, tcount);
                        int pk = 0;
                        for (int l = 0; l < npts_local; l++)
                            if (tmask[l]) packed[pk++] = tbuf[l];

                        int ov;
                        if (nc_inq_varid(out_ncids[n], vname, &ov) == NC_NOERR) {
                            size_t st[4]={0,0,0,0}, cnt[4]={0,1,1,1};
                            if (has_records) st[unlimdimid >= 0 ? 0 : 0] = (size_t)tlev;
                            cnt[0] = (size_t)tcount;
                            nc_put_vara_double(out_ncids[n], ov, st, cnt, packed);
                        }
                        free(packed);
                    }
                    free(tbuf); free(tmask);
                } else {
                    /* Non-compressed variable: write same value to all tiles. */
                    int ov;
                    if (nc_inq_varid(out_ncids[n], vname, &ov) == NC_NOERR) {
                        size_t st[1] = {(size_t)tlev};
                        size_t cnt[1] = {1};
                        if (has_records)
                            nc_put_vara_double(out_ncids[n], ov, st, cnt, buf);
                        else
                            nc_put_var_double(out_ncids[n], ov, buf);
                    }
                }
            }
            free(buf); free(mask);
        }
    }

    nc_close(in_ncid);
    for (int n = 0; n < nfiles_out; n++) { nc_sync(out_ncids[n]); nc_close(out_ncids[n]); }
    free(out_ncids); free(dimlen_list);
    return 0;
}

/* ================================================================== */
/* cmd_scatter – argument parsing + dispatch                           */
/* ================================================================== */
static void usage_scatter(void) {
    printf("Usage: mppdisttool scatter -x N -y N [-i N] [-j N] [-p prefix]\n"
           "                           [-s N] [-w N] [-X dims] [-Y dims] [-n] [-v] infile\n"
           "       mppdisttool scatter --land -i ndiv_x -j ndiv_y infile\n");
}

int
cmd_scatter(int argc, char **argv)
{
    scatter_opts_t opts;
    memset(&opts, 0, sizeof(opts));
    opts.width = 4;

    int land_mode = 0;
    int npex_land = 0, npey_land = 0;

    int dummy_n;
    newstringlist(&opts.xdims, &dummy_n, NC_MAX_DIMS);
    newstringlist(&opts.ydims, &dummy_n, NC_MAX_DIMS);

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i],"--land") == 0)                          land_mode = 1;
        else if (strcmp(argv[i],"-n") == 0)                         opts.dryrun = 1;
        else if (strcmp(argv[i],"-v") == 0 || strcmp(argv[i],"-V")==0) opts.verbose = 1;
        else if ((strcmp(argv[i],"-x") == 0) && i+1 < argc)        opts.nx = atoi(argv[++i]);
        else if ((strcmp(argv[i],"-y") == 0) && i+1 < argc)        opts.ny = atoi(argv[++i]);
        else if ((strcmp(argv[i],"-i") == 0) && i+1 < argc) {
            int v = atoi(argv[++i]);
            if (land_mode) npex_land = v; else opts.nxio = v;
        } else if ((strcmp(argv[i],"-j") == 0) && i+1 < argc) {
            int v = atoi(argv[++i]);
            if (land_mode) npey_land = v; else opts.nyio = v;
        } else if ((strcmp(argv[i],"-p") == 0) && i+1 < argc) {
            opts.prefix = argv[++i];
        } else if ((strcmp(argv[i],"-s") == 0) && i+1 < argc) {
            opts.start = atoi(argv[++i]);
        } else if ((strcmp(argv[i],"-w") == 0) && i+1 < argc) {
            opts.width = atoi(argv[++i]);
        } else if ((strcmp(argv[i],"-X") == 0) && i+1 < argc) {
            getstringlist(argv[++i], &opts.xdims, &opts.xdims_len);
        } else if ((strcmp(argv[i],"-Y") == 0) && i+1 < argc) {
            getstringlist(argv[++i], &opts.ydims, &opts.ydims_len);
        } else if (strcmp(argv[i],"--help") == 0 || strcmp(argv[i],"-h") == 0) {
            usage_scatter(); return 0;
        } else if (argv[i][0] != '-') {
            opts.filein = argv[i];
        }
    }

    if (!opts.filein) {
        fprintf(stderr,"scatter: missing input file.\n");
        usage_scatter(); return 1;
    }

    int ret;
    if (land_mode) {
        if (npex_land <= 0 || npey_land <= 0) {
            fprintf(stderr,"scatter --land: -i ndiv_x and -j ndiv_y required.\n");
            return 1;
        }
        ret = scatter_land(opts.filein, npex_land, npey_land);
    } else {
        if (opts.nx <= 0 || opts.ny <= 0) {
            fprintf(stderr,"scatter: -x N and -y N required.\n");
            usage_scatter(); return 1;
        }
        ret = scatter_mpp(&opts);
    }

    freestringlist(&opts.xdims, NC_MAX_DIMS);
    freestringlist(&opts.ydims, NC_MAX_DIMS);
    return ret;
}

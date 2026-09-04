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
 * Port of src/land_utils/decompress-ncc.F90.
 * Reads one or more compressed-by-gathering files and writes a single
 * decompressed (full-grid) NetCDF file.
 */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <netcdf.h>
#include "cmd_decompress.h"
#include "nc_utils.h"
#include "compress.h"
#include "xmalloc.h"

int
cmd_decompress(int argc, char **argv)
{
    int debug = 0;
    bool add_missing = false;
    const char **infiles = NULL;
    int ninfiles = 0;
    const char *outfile = NULL;

    /* Parse arguments. */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i],"-v")==0 || strcmp(argv[i],"-D")==0 ||
            strcmp(argv[i],"--debug-level")==0 || strcmp(argv[i],"--verbosity-level")==0) {
            if (i+1 < argc) debug = atoi(argv[++i]); else debug++;
        } else if (strcmp(argv[i],"-m")==0 || strcmp(argv[i],"--add-missing-value")==0) {
            add_missing = true;
        } else if (strcmp(argv[i],"-h")==0 || strcmp(argv[i],"--help")==0) {
            printf("Usage: mppdisttool decompress [-v level] [-m] infile... outfile\n");
            return 0;
        } else if (argv[i][0] != '-') {
            /* Collect positional args; last one is outfile. */
            if (!infiles) infiles = XMALLOC(const char *, argc);
            infiles[ninfiles++] = argv[i];
        }
    }

    if (ninfiles < 2) {
        fprintf(stderr,"decompress: at least one input and one output file required.\n");
        return 1;
    }
    outfile = infiles[--ninfiles]; /* Last arg is output. */

    /* Open input files. */
    int *input = XMALLOC(int, ninfiles);
    int in_format = NC_FORMAT_CLASSIC;
    for (int i = 0; i < ninfiles; i++) {
        if (nc_open(infiles[i], NC_NOWRITE, &input[i]) != NC_NOERR) {
            fprintf(stderr,"decompress: cannot open '%s'.\n", infiles[i]);
            free(input); return 1;
        }
        nc_inq_format(input[i], &in_format);
    }

    int ref = input[ninfiles - 1]; /* Template. */

    /* Determine output format. */
    int cmode;
    switch (in_format) {
    case NC_FORMAT_NETCDF4:         cmode = NC_NETCDF4; break;
    case NC_FORMAT_NETCDF4_CLASSIC: cmode = NC_NETCDF4 | NC_CLASSIC_MODEL; break;
    case NC_FORMAT_64BIT_OFFSET:    cmode = NC_CLOBBER | NC_64BIT_OFFSET; break;
    default:                        cmode = NC_CLOBBER | NC_64BIT_OFFSET; break;
    }

    size_t blksz = 65536;
    int ncid;
    if (nc__create(outfile, cmode, 0, &blksz, &ncid) != NC_NOERR) {
        fprintf(stderr,"decompress: cannot create '%s'.\n", outfile);
        free(input); return 1;
    }

    /* ---- Define dimensions (skip compressed, include expanded ones) ---- */
    int ndims;
    nc_inq_ndims(ref, &ndims);
    for (int d = 0; d < ndims; d++) {
        char dname[NC_MAX_NAME+1]; size_t dlen; bool is_unlim;
        ncu_inq_dim(ref, d, dname, &dlen, &is_unlim);

        /* Skip compressed dimensions. */
        int varid, attid;
        if (nc_inq_varid(ref, dname, &varid) == NC_NOERR &&
            nc_inq_attid(ref, varid, "compress", &attid) == NC_NOERR)
            continue;

        int dummy;
        nc_def_dim(ncid, dname, is_unlim ? NC_UNLIMITED : dlen, &dummy);
    }

    /* ---- Define variables, replacing compressed dims with expanded dims ---- */
    int nvars;
    nc_inq_nvars(ref, &nvars);
    for (int v = 0; v < nvars; v++) {
        char vname[NC_MAX_NAME+1];
        nc_type xtype;
        int vndims, dimids[NC_MAX_VARS], natts;
        nc_inq_var(ref, v, vname, &xtype, &vndims, dimids, &natts);
        if (debug) printf("decompress: defining '%s'\n", vname);

        /* Is this a compressed dim variable? */
        bool is_dim = false, is_comp = false;
        int dummy_did;
        if (nc_inq_dimid(ref, vname, &dummy_did) == NC_NOERR) is_dim = true;
        int attid;
        if (is_dim && nc_inq_varid(ref,vname,&dummy_did)==NC_NOERR &&
            nc_inq_attid(ref,dummy_did,"compress",&attid)==NC_NOERR) is_comp = true;
        if (is_dim && is_comp) continue; /* Don't define compressed dim vars. */

        /* Expand any compressed dims in the variable's dim list. */
        int new_ndims = 0;
        int new_dimids[NC_MAX_DIMS];
        for (int d = 0; d < vndims; d++) {
            char dname[NC_MAX_NAME+1];
            nc_inq_dimname(ref, dimids[d], dname);
            /* Is this dim compressed? */
            int cdims; int cdimids[NC_MAX_DIMS]; size_t cdimlens[NC_MAX_DIMS];
            if (compress_inq_dim(ref, dimids[d], &cdims, cdimids, cdimlens) == NC_NOERR) {
                /* Replace with expanded dims. */
                for (int k = 0; k < cdims; k++) {
                    char cdname[NC_MAX_NAME+1];
                    nc_inq_dimname(ref, cdimids[k], cdname);
                    int new_did;
                    if (nc_inq_dimid(ncid, cdname, &new_did) == NC_NOERR)
                        new_dimids[new_ndims++] = new_did;
                }
            } else {
                int new_did;
                if (nc_inq_dimid(ncid, dname, &new_did) == NC_NOERR)
                    new_dimids[new_ndims++] = new_did;
            }
        }

        int new_varid;
        nc_def_var(ncid, vname, xtype, new_ndims, new_dimids, &new_varid);
        for (int n = 0; n < natts; n++) {
            char attname[NC_MAX_NAME+1];
            nc_inq_attname(ref, v, n, attname);
            nc_copy_att(ref, v, attname, ncid, new_varid);
        }

        /* Optionally add missing_value if compressed and not present. */
        if (add_missing && is_comp) {
            int attid2;
            if (nc_inq_attid(ref, v, "missing_value", &attid2) != NC_NOERR &&
                nc_inq_attid(ref, v, "_FillValue",    &attid2) != NC_NOERR) {
                double mv;
                switch (xtype) {
                case NC_DOUBLE: mv = NC_FILL_DOUBLE; break;
                case NC_FLOAT:  mv = NC_FILL_FLOAT;  break;
                default:        mv = NC_FILL_INT;    break;
                }
                nc_put_att_double(ncid, new_varid, "missing_value", xtype, 1, &mv);
            }
        }
    }

    /* Copy global attributes. */
    ncu_copy_global_atts(ref, ncid);
    nc__enddef(ncid, 16384, 4, 0, 4);

    /* Extend record dimension if needed. */
    for (int v = 0; v < nvars; v++) {
        size_t cdimlens[NC_MAX_DIMS];
        char vname[NC_MAX_NAME+1];
        nc_type xtype;
        int vndims, dimids[NC_MAX_DIMS], natts;
        ncu_inq_var(ref, v, vname, &xtype, &vndims, dimids, cdimlens, &natts, NULL, NULL, NULL, NULL, NULL);
        /* If has record dim, write a dummy int to extend. */
        int unlimdim_out;
        nc_inq_unlimdim(ncid, &unlimdim_out);
        if (unlimdim_out < 0) break;
        char unlim_name[NC_MAX_NAME+1];
        nc_inq_dimname(ncid, unlimdim_out, unlim_name);
        int ov;
        if (nc_inq_varid(ncid, vname, &ov) != NC_NOERR) continue;
        int ov_ndims, ov_dimids[NC_MAX_DIMS];
        nc_inq_varndims(ncid, ov, &ov_ndims);
        nc_inq_vardimid(ncid, ov, ov_dimids);
        bool ov_has_rec = false;
        for (int d = 0; d < ov_ndims; d++)
            if (ov_dimids[d] == unlimdim_out) { ov_has_rec = true; break; }
        if (!ov_has_rec) continue;

        /* Get max dimlens of this var. */
        size_t end_idx[NC_MAX_DIMS];
        for (int d = 0; d < ov_ndims; d++) {
            size_t dlen; nc_inq_dimlen(ref, dimids[d < vndims ? d : 0], &dlen);
            end_idx[d] = dlen > 0 ? dlen - 1 : 0;
        }
        nc_put_var1_int(ncid, ov, end_idx, (const int[]){0});
        break;
    }

    /* ---- Gather compressed data and write to output ---- */
    int out_nvars;
    nc_inq_nvars(ncid, &out_nvars);
    for (int v = 0; v < out_nvars; v++) {
        char vname[NC_MAX_NAME+1];
        nc_type xtype;
        size_t vsize;
        ncu_inq_var(ncid, v, vname, &xtype, NULL, NULL, NULL, NULL, NULL, NULL, &vsize, NULL, NULL);
        if (debug) printf("decompress: processing '%s'\n", vname);

        double *buf  = XMALLOC(double, vsize > 0 ? vsize : 1);
        bool   *mask = XMALLOC(bool,   vsize > 0 ? vsize : 1);

        /* Determine missing value. */
        double missing = NC_FILL_DOUBLE;
        ncu_inq_var(ncid, v, NULL, &xtype, NULL, NULL, NULL, NULL, NULL, NULL, &vsize, NULL, NULL);
        double ocean_val = 0.0;
        bool do_ocean = (nc_get_att_double(ncid, v, "ocean_fillvalue", &ocean_val) == NC_NOERR);
        if (nc_get_att_double(ncid, v, "missing_value", &missing) != NC_NOERR)
            if (nc_get_att_double(ncid, v, "_FillValue",  &missing) != NC_NOERR) {
                switch (xtype) {
                case NC_DOUBLE: missing = NC_FILL_DOUBLE; break;
                case NC_FLOAT:  missing = NC_FILL_FLOAT; break;
                default:        missing = NC_FILL_INT; break;
                }
            }
        for (size_t k = 0; k < vsize; k++) buf[k] = missing;
        memset(mask, 0, vsize * sizeof(bool));

        /* Read from all input files. */
        for (int fi = 0; fi < ninfiles; fi++) {
            int iv;
            if (nc_inq_varid(input[fi], vname, &iv) != NC_NOERR) continue;
            compress_get_var_double(input[fi], iv, buf, mask, vsize);
        }

        if (do_ocean)
            for (size_t k = 0; k < vsize; k++)
                if (!mask[k]) buf[k] = ocean_val;

        /* Write to output. */
        nc_put_var_double(ncid, v, buf);
        free(buf); free(mask);
    }

    nc_sync(ncid); nc_close(ncid);
    for (int i = 0; i < ninfiles; i++) nc_close(input[i]);
    free(input);
    if (infiles) free((void*)infiles);
    return 0;
}

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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <netcdf.h>
#include "cmd_check.h"
#include "compress.h"

/* Icebergs file format version thresholds (matches iceberg_comb.sh). */
#define ICEBERG_MAJOR_MIN 0
#define ICEBERG_MINOR_MIN 1

bool
check_is_compressed(const char *filename)
{
    int ncid;
    if (nc_open(filename, NC_NOWRITE, &ncid) != NC_NOERR) return false;
    bool result = file_is_compressed(ncid);
    nc_close(ncid);
    return result;
}

bool
check_is_iceberg(const char *filename)
{
    int ncid;
    if (nc_open(filename, NC_NOWRITE, &ncid) != NC_NOERR) return false;

    int major = -1, minor = -1, nfiles = -1;
    nc_get_att_int(ncid, NC_GLOBAL, "file_format_major_version", &major);
    nc_get_att_int(ncid, NC_GLOBAL, "file_format_minor_version", &minor);
    nc_get_att_int(ncid, NC_GLOBAL, "NumFilesInSet",             &nfiles);
    nc_close(ncid);

    if (major < 0 || minor < 0 || nfiles <= 0) return false;
    if (major > ICEBERG_MAJOR_MIN) return true;
    if (major == ICEBERG_MAJOR_MIN && minor >= ICEBERG_MINOR_MIN) return true;
    return false;
}

bool
check_iceberg_has_data(const char *filename)
{
    int ncid, unlimdimid;
    if (nc_open(filename, NC_NOWRITE, &ncid) != NC_NOERR) return false;
    if (nc_inq_unlimdim(ncid, &unlimdimid) != NC_NOERR || unlimdimid < 0) {
        nc_close(ncid);
        return false;
    }
    size_t dimlen = 0;
    nc_inq_dimlen(ncid, unlimdimid, &dimlen);
    nc_close(ncid);
    return dimlen > 0;
}

int
cmd_check(int argc, char **argv)
{
    bool do_compressed = false;
    bool do_iceberg    = false;
    bool quiet         = false;
    const char *filename = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--compressed") == 0)  { do_compressed = true; }
        else if (strcmp(argv[i], "--iceberg") == 0) { do_iceberg = true; }
        else if (strcmp(argv[i], "--quiet") == 0)   { quiet = true; }
        else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printf("Usage: mppdisttool check [--compressed|--iceberg] [--quiet] FILE\n");
            return 0;
        } else if (argv[i][0] != '-') {
            filename = argv[i];
        }
    }

    if (!filename) {
        fprintf(stderr, "check: missing file operand\n");
        return 1;
    }
    if (!do_compressed && !do_iceberg) {
        /* Default: behave like is-compressed */
        do_compressed = true;
    }

    if (do_compressed) {
        bool r = check_is_compressed(filename);
        if (!quiet && r)
            printf("%s is compressed\n", filename);
        return r ? 0 : -1;
    }
    if (do_iceberg) {
        bool r = check_is_iceberg(filename);
        if (!quiet && r)
            printf("%s is an iceberg restart file\n", filename);
        return r ? 0 : 255;
    }
    return 1;
}

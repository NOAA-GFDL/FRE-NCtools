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
#ifndef MPPDISTTOOL_CMD_COMBINE_H
#define MPPDISTTOOL_CMD_COMBINE_H

/* Options shared by all combine paths. */
typedef struct {
    char        *outfile;       /* output file name */
    char       **infiles;       /* input file names */
    int          ninfiles;
    int          verbose;
    int          removein;      /* remove input files on success */
    int          appendnc;      /* append to existing output */
    int          headerpad;     /* extra header padding bytes */
    int          missing;       /* init decomposed vars with missing_value */
    int          force;         /* combine even if files missing */
    int          format_64;     /* force 64-bit offset output */
    int          format_n4;     /* force NETCDF4_CLASSIC output */
    int          deflate;       /* deflate (NETCDF4 only) */
    int          deflation;     /* deflation level */
    int          shuffle;       /* shuffle filter */
    int          nstart;        /* first PE extension (no explicit infiles) */
    int          nend;          /* last PE extension (-1 = auto) */
    int          blocking;      /* record blocking factor */
    int          mem_dry_run;   /* print estimated memory and exit */
    int          print_mem;     /* print memory usage statistics */
    /* land path */
    int          land_mode;
    /* iceberg path */
    int          iceberg_mode;
    /* MPP (default) */
    int          mpp_mode;
} combine_opts_t;

/* Entry point: "mppdisttool combine [opts] [-o out] infiles..." */
int cmd_combine(int argc, char **argv);

/* Individual path entry points (used by cmd_auto). */
int combine_land(combine_opts_t *opts);
int combine_iceberg(combine_opts_t *opts);
int combine_mpp(combine_opts_t *opts);

#endif /* MPPDISTTOOL_CMD_COMBINE_H */

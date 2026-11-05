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
 * Port of src/combine-restarts/combine_restarts.in.
 * Scans CWD for files matching *{res,nc}.[0-9]{4}, groups by base name,
 * then routes each group through the appropriate combine path.
 */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <dirent.h>
#include <regex.h>
#include <unistd.h>
#include <netcdf.h>
#include "cmd_auto.h"
#include "cmd_combine.h"
#include "cmd_check.h"
#include "compress.h"
#include "xmalloc.h"

/* ------------------------------------------------------------------ */
/* Simple dynamic string array                                          */
/* ------------------------------------------------------------------ */
typedef struct { char **items; int n, cap; } strarray_t;

static void sa_init(strarray_t *a) { a->items=NULL; a->n=0; a->cap=0; }
static void sa_push(strarray_t *a, const char *s) {
    if (a->n >= a->cap) {
        a->cap = a->cap ? a->cap*2 : 16;
        a->items = (char**)realloc(a->items, (size_t)a->cap * sizeof(char*));
    }
    a->items[a->n++] = strdup(s);
}
static void sa_free(strarray_t *a) {
    for (int i=0;i<a->n;i++) free(a->items[i]);
    free(a->items); a->n=0; a->cap=0;
}
static int sa_contains(strarray_t *a, const char *s) {
    for (int i=0;i<a->n;i++) if (strcmp(a->items[i],s)==0) return 1;
    return 0;
}
static int sa_cmp(const void *a, const void *b) {
    return strcmp(*(const char**)a, *(const char**)b);
}

/* ------------------------------------------------------------------ */
/* File grouping                                                        */
/* ------------------------------------------------------------------ */
/* Match files in CWD whose names contain "res" or "nc" as a word and
   end in .[0-9]{4}.  Returns base names (suffix stripped) de-duplicated. */
static strarray_t
find_base_names(int verbose)
{
    strarray_t bases;
    sa_init(&bases);

    regex_t re_tail, re_word;
    /* .[0-9]{4} at end of name */
    if (regcomp(&re_tail, "\\.[0-9][0-9][0-9][0-9]$", REG_EXTENDED|REG_NOSUB) != 0) return bases;
    /* word "res" or "nc" (bounded by non-alphanumeric or start/end) */
    if (regcomp(&re_word, "(^|[^a-zA-Z0-9_])(res|nc)([^a-zA-Z0-9_]|$)", REG_EXTENDED|REG_NOSUB) != 0) {
        regfree(&re_tail); return bases;
    }

    DIR *dir = opendir(".");
    if (!dir) { regfree(&re_tail); regfree(&re_word); return bases; }

    struct dirent *ent;
    while ((ent = readdir(dir)) != NULL) {
        const char *name = ent->d_name;
        if (regexec(&re_tail, name, 0, NULL, 0) != 0) continue;
        if (regexec(&re_word, name, 0, NULL, 0) != 0) continue;

        /* Strip trailing .[0-9]{4} to get base name. */
        char base[4096];
        size_t blen = strlen(name);
        if (blen < 5) continue;  /* Need at least ".NNNN" */
        if (blen - 5 >= sizeof(base)) {
            fprintf(stderr, "Warning: filename too long, skipping: %s\n", name);
            continue;
        }
        memcpy(base, name, blen - 5);
        base[blen - 5] = '\0'; /* remove ".NNNN" */

        if (!sa_contains(&bases, base)) {
            if (verbose > 1) fprintf(stderr,"debug2: found base '%s'\n", base);
            sa_push(&bases, base);
        }
    }
    closedir(dir);
    regfree(&re_tail); regfree(&re_word);

    /* Sort for deterministic order. */
    if (bases.n > 0)
        qsort(bases.items, (size_t)bases.n, sizeof(char*), sa_cmp);
    return bases;
}

/* Collect input files for a given base name. */
static strarray_t
collect_inputs(const char *base)
{
    strarray_t files;
    sa_init(&files);

    char pattern[4096];
    snprintf(pattern, sizeof(pattern), "\\.[0-9][0-9][0-9][0-9]$");
    regex_t re;
    if (regcomp(&re, pattern, REG_EXTENDED|REG_NOSUB) != 0) return files;

    DIR *dir = opendir(".");
    if (!dir) { regfree(&re); return files; }

    struct dirent *ent;
    size_t blen = strlen(base);
    while ((ent = readdir(dir)) != NULL) {
        const char *name = ent->d_name;
        if (strncmp(name, base, blen) != 0) continue;
        if (strlen(name) != blen + 5) continue; /* base + .NNNN */
        if (regexec(&re, name, 0, NULL, 0) != 0) continue;
        sa_push(&files, name);
    }
    closedir(dir);
    regfree(&re);

    if (files.n > 0)
        qsort(files.items, (size_t)files.n, sizeof(char*), sa_cmp);
    return files;
}

/* ------------------------------------------------------------------ */
/* cmd_auto                                                             */
/* ------------------------------------------------------------------ */
int
cmd_auto(int argc, char **argv)
{
    int verbose   = 0;
    int dryrun    = 0;
    int removein  = 0;
    const char *combineNccOpts = NULL;
    const char *mppCombOpts    = "-h 16384 -m";

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i],"-v")==0)        verbose++;
        else if (strcmp(argv[i],"-n")==0)   dryrun = 1;
        else if (strcmp(argv[i],"-r")==0)   removein = 1;
        else if (strcmp(argv[i],"-C")==0 && i+1<argc) combineNccOpts = argv[++i];
        else if (strcmp(argv[i],"-M")==0 && i+1<argc) mppCombOpts    = argv[++i];
        else if (strcmp(argv[i],"-h")==0 || strcmp(argv[i],"--help")==0) {
            printf("Usage: mppdisttool auto [-C combine-ncc-opts] [-M mppnccombine-opts]\n"
                   "                        [-n] [-r] [-v]\n");
            return 0;
        }
    }

    strarray_t bases = find_base_names(verbose);
    if (bases.n == 0) {
        fprintf(stderr,"auto: no restart files found in current directory.\n");
        return 0;
    }

    int ret = 0;
    for (int b = 0; b < bases.n; b++) {
        const char *base = bases.items[b];
        if (verbose) fprintf(stderr,"debug1: processing '%s'\n", base);

        strarray_t inputs = collect_inputs(base);
        if (inputs.n == 0) { sa_free(&inputs); continue; }

        /* Build combine_opts_t. */
        combine_opts_t opts;
        memset(&opts, 0, sizeof(opts));
        opts.outfile   = (char*)base;
        opts.infiles   = inputs.items;
        opts.ninfiles  = inputs.n;
        opts.verbose   = verbose;
        opts.removein  = removein;
        opts.headerpad = 16384;
        opts.blocking  = 1;
        opts.nend      = -1;

        /* Parse extra options passed via -C / -M. */
        if (combineNccOpts) {
            /* Just ignore for now — all options go through combine_opts_t. */
        }
        if (mppCombOpts) {
            /* Parse simple known flags from mppCombOpts string. */
            char buf[1024];
            strncpy(buf, mppCombOpts, sizeof(buf)-1);
            char *tok = strtok(buf, " \t");
            while (tok) {
                if (strcmp(tok,"-h")==0 && (tok=strtok(NULL," \t"))) opts.headerpad = atoi(tok);
                else if (strcmp(tok,"-k")==0 && (tok=strtok(NULL," \t"))) opts.blocking = atoi(tok);
                else if (strcmp(tok,"-m")==0) opts.missing = 1;
                else if (strcmp(tok,"-64")==0) opts.format_64 = 1;
                else if (strcmp(tok,"-n4")==0) opts.format_n4 = 1;
                tok = strtok(NULL, " \t");
            }
        }

        /* Probe first file for type. */
        int probeid = -1;
        nc_open(inputs.items[0], NC_NOWRITE, &probeid);

        int is_compressed_file = 0, is_iceberg_file = 0;
        if (probeid >= 0) {
            is_compressed_file = file_is_compressed(probeid) ? 1 : 0;
            nc_close(probeid);
        }
        if (!is_compressed_file)
            is_iceberg_file = check_is_iceberg(inputs.items[0]) ? 1 : 0;

        if (verbose > 1) {
            fprintf(stderr,"debug2: '%s': compressed=%d iceberg=%d\n",
                    base, is_compressed_file, is_iceberg_file);
        }

        int r = 0;
        if (dryrun) {
            /* Print what would be executed. */
            if (is_compressed_file) {
                printf("combine-ncc");
                if (combineNccOpts) printf(" %s", combineNccOpts);
            } else if (is_iceberg_file) {
                printf("iceberg_comb.sh");
            } else {
                printf("mppnccombine");
                if (mppCombOpts) printf(" %s", mppCombOpts);
                printf(" %s", base);
            }
            for (int f = 0; f < inputs.n; f++) printf(" %s", inputs.items[f]);
            if (is_compressed_file) printf(" %s", base);
            printf("\n");
            if (removein) {
                printf("rm -f");
                for (int f = 0; f < inputs.n; f++) printf(" %s", inputs.items[f]);
                printf("\n");
            }
        } else {
            if (is_compressed_file) {
                opts.land_mode = 1;
                r = combine_land(&opts);
            } else if (is_iceberg_file) {
                opts.iceberg_mode = 1;
                r = combine_iceberg(&opts);
            } else {
                opts.mpp_mode = 1;
                r = combine_mpp(&opts);
            }
            if (r != 0) {
                fprintf(stderr,"auto: error combining '%s'\n", base);
                ret = 1;
            }
        }
        sa_free(&inputs);
    }
    sa_free(&bases);
    return ret;
}

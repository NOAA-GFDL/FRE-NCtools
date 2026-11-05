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
 * mppdisttool – single binary replacing seven MPP distribution tools.
 * Supports subcommand dispatch and argv[0] backward-compatibility symlinks.
 */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>   /* basename() */
#include "cmd_check.h"
#include "cmd_combine.h"
#include "cmd_scatter.h"
#include "cmd_decompress.h"
#include "cmd_auto.h"

/* print_version is defined in src/print_version.c and linked via libver.a. */
extern void print_version(const char *command_name);

/* ------------------------------------------------------------------ */
/* Subcommand dispatch table                                            */
/* ------------------------------------------------------------------ */
typedef struct { const char *name; int (*fn)(int, char **); } subcmd_t;

static const subcmd_t commands[] = {
    { "check",      cmd_check      },
    { "combine",    cmd_combine    },
    { "scatter",    cmd_scatter    },
    { "decompress", cmd_decompress },
    { "auto",       cmd_auto       },
    { NULL, NULL }
};

static void usage(void) {
    printf("mppdisttool – unified MPP distribution tool\n\n"
           "Usage: mppdisttool <subcommand> [options]\n\n"
           "Subcommands:\n"
           "  check      Test if a file is compressed or an iceberg restart\n"
           "  combine    Combine decomposed NetCDF files into one\n"
           "  scatter    Decompose a NetCDF file into many\n"
           "  decompress Expand compressed-by-gathering NetCDF files\n"
           "  auto       Auto-detect and combine restart files in CWD\n"
           "\nBackward-compatibility symlinks (argv[0] dispatch):\n"
           "  mppncscatter    -> scatter\n"
           "  mppnccombine    -> combine --mpp\n"
           "  combine-ncc     -> combine --land\n"
           "  scatter-ncc     -> scatter --land\n"
           "  decompress-ncc  -> decompress\n"
           "  is-compressed   -> check --compressed\n"
           "\nRun 'mppdisttool <subcommand> --help' for subcommand options.\n");
}

/* ------------------------------------------------------------------ */
/* argv[0] backward-compatibility mapping                              */
/* ------------------------------------------------------------------ */
/* Returns extra argv words to prepend, or NULL if no match. */
static const char *compat_argv0_extra(const char *progname) {
    if (strcmp(progname, "mppncscatter")  == 0) return "scatter";
    if (strcmp(progname, "mppnccombine")  == 0) return "combine\0--mpp";
    if (strcmp(progname, "combine-ncc")   == 0) return "combine\0--land";
    if (strcmp(progname, "scatter-ncc")   == 0) return "scatter\0--land";
    if (strcmp(progname, "decompress-ncc")== 0) return "decompress";
    if (strcmp(progname, "is-compressed") == 0) return "check\0--compressed";
    return NULL;
}

/* Build a new argv with extra words prepended after argv[0]. */
static char **build_compat_argv(int argc, char **argv,
                                 const char *extra, int *new_argc_out) {
    /* Count extra words (NUL-terminated, double-NUL ends sequence). */
    int nextra = 0;
    const char *p = extra;
    while (*p) { nextra++; p += strlen(p) + 1; }

    int new_argc = argc + nextra;
    char **new_argv = (char **)malloc((size_t)(new_argc + 1) * sizeof(char *));
    new_argv[0] = argv[0];
    p = extra;
    for (int i = 1; i <= nextra; i++) {
        new_argv[i] = (char *)p;
        p += strlen(p) + 1;
    }
    for (int i = 1; i < argc; i++)
        new_argv[nextra + i] = argv[i];
    new_argv[new_argc] = NULL;
    *new_argc_out = new_argc;
    return new_argv;
}

/* ------------------------------------------------------------------ */
/* main                                                                 */
/* ------------------------------------------------------------------ */
int
main(int argc, char *argv[])
{
    /* Determine invocation name for argv[0] compat. */
    char *prog_copy = strdup(argv[0]);
    const char *progname = basename(prog_copy);

    /* Check for argv[0] backward-compatibility mapping. */
    const char *extra = compat_argv0_extra(progname);
    if (extra) {
        int new_argc;
        char **new_argv = build_compat_argv(argc, argv, extra, &new_argc);
        free(prog_copy);
        /* Dispatch through subcommand. */
        const char *subcmd = new_argv[1];
        for (const subcmd_t *c = commands; c->name; c++) {
            if (strcmp(c->name, subcmd) == 0) {
                int ret = c->fn(new_argc - 1, new_argv + 1);
                free(new_argv);
                return ret;
            }
        }
        free(new_argv);
        fprintf(stderr,"mppdisttool: unknown subcommand '%s'\n", subcmd);
        return 1;
    }
    free(prog_copy);

    /* Normal subcommand dispatch. */
    if (argc < 2) {
        usage();
        return 1;
    }

    const char *subcmd = argv[1];

    if (strcmp(subcmd,"-V")==0 || strcmp(subcmd,"--version")==0) {
        print_version("mppdisttool");
        return 0;
    }
    if (strcmp(subcmd,"-h")==0 || strcmp(subcmd,"--help")==0) {
        usage();
        return 0;
    }

    for (const subcmd_t *c = commands; c->name; c++) {
        if (strcmp(c->name, subcmd) == 0)
            return c->fn(argc - 1, argv + 1);
    }

    fprintf(stderr,"mppdisttool: unknown subcommand '%s'\n"
            "Run 'mppdisttool --help' for usage.\n", subcmd);
    return 1;
}

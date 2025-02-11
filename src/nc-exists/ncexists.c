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
  Copyright (C) 2011 NOAA Geophysical Fluid Dynamics Lab, Princeton, NJ
*/
#include <getopt.h>
#include <ctype.h>
#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "version.h"

char *program_name = "ncexists";

void
usage (int status)
{
  fprintf (stdout,
           "Usage: %s -f filename -g global_attribute\n"
           "   or: %s -f filename -v variable [-a attribute]\n",
           program_name, program_name);
  fputs ("Prints 1 to stdout if variable or attribute is found, "
         "0 if not found.\n",
         stdout);
  fputs ("\n"
         "  -f, --filename=FILENAME    netCDF file\n"
         "  -g, --global=ATTRIBUTE     global attribute\n"
         "  -v, --variable=VARIABLE    variable\n"
         "  -a, --attribute=ATTRIBUTE  variable attribute\n"
         "  -h, --help                 display this help and exit\n"
         "  -V, --version              output version information and exit\n",
         stdout);
  exit (status);
}

void
handle_error (int status)
{
  if (status != NC_NOERR)
    {
      fprintf (stderr, "%s: %s\n",
               program_name, nc_strerror (status));
      exit (EXIT_FAILURE);
    }
}

int
main (int argc, char **argv)
{
  int status;
  char *filename = NULL;
  char *attr = NULL;
  char *gattr = NULL;
  int ncid;
  char *var_name = NULL;
  size_t vr_len;
  int var_id;
  size_t t_len;
  int global = 0;
  int c;
  nc_type vr_type, t_type;

  static struct option const long_options[] = {
    {"filename", required_argument, 0, 'f'},
    {"variable", required_argument, 0, 'v'},
    {"attribute", required_argument, 0, 'a'},
    {"global", required_argument, 0, 'g'},
    {"help", no_argument, 0, 'h'},
    {"version", no_argument, 0, 'V'},
    {0, 0, 0, 0}
  };

  if (argc == 1)
    {
      fprintf (stderr, "%s: No arguments specified\n", program_name);
      usage (EXIT_FAILURE);
    }
  else
    {
      while ((c = getopt_long (argc, argv,
                               "hVf:v:g:a:",
                               long_options, NULL)) != -1)
        {
          switch (c)
            {
            case 'f':
              filename = optarg;
              break;

            case 'g':
              gattr = optarg;
              global = 1;
              break;

            case 'v':
              var_name = optarg;
              break;

            case 'a':
              attr = optarg;
              break;

            case 'h':
              usage (EXIT_SUCCESS);
              break;

            case 'V':
              print_version (program_name);
              exit (EXIT_SUCCESS);
              break;

            case '?':
              usage(EXIT_FAILURE);
              break;

            default:
              usage(EXIT_FAILURE);
              break;
            }
        }

      // Ensure the user has specified a filename
      if (filename == NULL)
        {
          fprintf (stderr, "%s: No filename specified\n",
                   program_name);
          usage (EXIT_FAILURE);
        }

      // Ensure the user has specified a variable or global attribute
      if (var_name == NULL && gattr == NULL)
        {
          fprintf (stderr, "%s: No variable or global attribute specified\n",
                   program_name);
          exit (EXIT_FAILURE);
        }

      // check if the filename exists and is readable
      int access_status = access (filename, R_OK);
      if (access_status != 0)
        {
          fprintf (stderr, "%s: %s: %s\n",
                   program_name, filename, strerror (errno));
          exit (EXIT_FAILURE);
        }

      // Open the netCDF file
      handle_error (nc_open (filename, NC_NOWRITE, &ncid));

      if (global == 1)
        {
          status = nc_inq_att (ncid, NC_GLOBAL, gattr, &t_type, &t_len);
          if (status == NC_NOERR)
            {
              printf ("1\n");
            }
          else
            {
              printf ("0\n");
            }
        }
      else
        {
          status = nc_inq_varid (ncid, var_name, &var_id);
          if (status == NC_NOERR)
            {
              if (attr == NULL)
                {
                  printf ("1\n");
                }
              else
                {
                  status = nc_inq_att (ncid, var_id, attr, &vr_type, &vr_len);
                  if (status == NC_NOERR)
                    {
                      printf ("1\n");
                    }
                  else
                    {
                      printf ("0\n");
                    }
                }
            }
          else
            {
              printf ("0\n");
            }
        }

      status = nc_close (ncid); /* close netCDF dataset */
      if (status != NC_NOERR)
        handle_error (status);
    }

  exit (EXIT_SUCCESS);
}

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
#include <config.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <netcdf.h>

#include "version.h"

#define CHECK_NC_ERRSTAT(ierr) check_error(ierr,__FILE__,__LINE__)

#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

void
usage (const char* program_name, int status)
{
   if (status != EXIT_SUCCESS)
      fprintf (stderr, "Try `%s --help' for more information.\n",
               program_name);
   else
      {
         printf ("Usage: %s [OPTION]... FILE\n", program_name);
         fputs("\
Test if any variable in FILE is compressed-by-gathering.\n\
Exit status is 0 if compressed, -1 if not compressed.\n\
\n\
  -d, --debug              add verbosity to the output\n\
  -h, --help               display this help and exit\n\
  -V, --version            output version information and exit\n",
               stdout);
         printf("\
\n\
Examples:\n\
  %s -d cana.res.nc\n",
          program_name);
      }
   exit (status);
}

void check_error(int ierr, const char* file, int line)
{
   if (ierr!=NC_NOERR) {
      fprintf(stderr,"ERROR :: %s\nfile: %s, line: %d\n",
         nc_strerror(ierr), file, line);
      exit(ierr);
   }
}

int main(int argc, char* argv[])
{
   static struct option const long_options[] =
   {
      {"debug", no_argument, NULL, 'd'},
      {"help", no_argument, NULL, 'h'},
      {"version", no_argument, NULL, 'V'}
   };

   int c;
   int ncid; // id of netcdf file
   int dimid; // dimension id
   int varid; // variable id
   int attid; // attribute id
   int ndims; // number of dimensions
   int verbose=0; // level of verbosity
   char name[NC_MAX_NAME+1];

   if(argc<2) {
      fprintf (stderr,
               "%s: missing file operand\n",
               argv[0]);
      usage(argv[0], EXIT_FAILURE);
   }

   // parse command-line arguments
   while ((c = getopt_long (argc, argv, "dhV", long_options, NULL)) != -1)
      {
         switch (c)
            {
            case 'd':
               verbose++;
               break;

            case 'V':
               print_version(argv[0]);
               exit (EXIT_SUCCESS);
               break;

            case 'h':
               usage(argv[0], EXIT_SUCCESS);
               break;

            default:
               usage(argv[0], EXIT_FAILURE);
            }
      }

   int optind = argc - 1;

   CHECK_NC_ERRSTAT(nc_open(argv[optind], NC_NOWRITE, &ncid));
   CHECK_NC_ERRSTAT(nc_inq_ndims(ncid, &ndims));
   if (verbose) fprintf(stderr, "Found %d dimensions\n",ndims);
   for (dimid=0; dimid<ndims; dimid++)
      {
         CHECK_NC_ERRSTAT(nc_inq_dimname(ncid, dimid, name));
         if (verbose) fprintf(stderr, "Examining \"%s\"\n",name);
         if (nc_inq_varid(ncid, name, &varid)==NC_NOERR &&
             nc_inq_attid(ncid, varid, "compress", &attid)==NC_NOERR)
            {
               return EXIT_SUCCESS;
            }
      }
   CHECK_NC_ERRSTAT(nc_close(ncid));

   return -1;
}

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
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <netcdf.h>

#define CHECK_NC_ERRSTAT(ierr) check_error(ierr,__FILE__,__LINE__)

void usage(const char* name)
{
   printf("==============================================================================\n");
   printf("%s -- tests for compressed-by-gathering variables\n",name);
   printf("==============================================================================\n");
   printf("This utility tests if any of the variables in the given netcdf file are \n");
   printf("compressed-by-gathering. It returns 0 if yes, and -1 otherwise. \n");

   printf("\nFor information of compression, see:");
   printf("http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.5/cf-conventions.html#compression-by-gathering\n");

   printf("\nUsage:\n");
   printf("%s [-d] [-h] file-name\n", name);
   printf("  -d -- adds verbosity to the output\n");
   printf("  -h -- prints this help message and exits\n");

   printf("\nFor example:\n");
   printf("%s -d cana.res.nc.0000\n", name);

   printf("==============================================================================\n");
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
   int c;
   int ncid; // id of netcdf file
   int dimid; // dimension id
   int varid; // variable id
   int attid; // attribute id
   int ndims; // number of dimensions
   int verbose=0; // level of verbosity
   char name[NC_MAX_NAME+1];
   
   if(argc<2) {
      usage(argv[0]);
      return 1;
   }
   
   // parse command-line arguments
   while ((c = getopt(argc, argv, "dh")) != -1)
      switch (c) {
      case 'd':
	 verbose++;
	 break;
      case 'h':
      default:
	 usage(argv[0]);
	 return 1;
      }

   int optind = argc - 1;
   
   CHECK_NC_ERRSTAT(nc_open(argv[optind],NC_NOWRITE,&ncid));
   CHECK_NC_ERRSTAT(nc_inq_ndims(ncid,&ndims));
   if (verbose) printf("Found %d dimensions\n",ndims);
   for (dimid=0; dimid<ndims; dimid++) {
      CHECK_NC_ERRSTAT(nc_inq_dimname(ncid,dimid,name));
      if (verbose) printf("Examining \"%s\"\n",name);
      if (nc_inq_varid(ncid,name,&varid)==NC_NOERR &&
          nc_inq_attid(ncid,varid,"compress",&attid)==NC_NOERR) {
          return 0;
      }
   }
   CHECK_NC_ERRSTAT(nc_sync(ncid));
   CHECK_NC_ERRSTAT(nc_close(ncid));

   return -1;
}

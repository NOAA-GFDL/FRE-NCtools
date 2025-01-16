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
#include "opt.h"
#include "mppncscatter.h"

/*-------------------------------------------------------------------*/
int main(int argc, char** argv)
{
	int status;
	mnsopts opts;

	status = 0;

	initmnsopts(&opts);
	
	/* parse command-line args.  */
	status = getmnsopts(argc, argv, &opts);
	
	status = (status ? status : mppncscatter(&opts));               
	
	freemnsopts(&opts);
	
 	if (opts.verbose) {
		fprintf(stdout, "Info: Done.\n");
		fflush(stdout);
	}

	return status;
} 
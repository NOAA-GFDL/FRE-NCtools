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
#include <math.h>
#include <string.h>
#include "fregrid_util.h"
#include "mpp.h"
#include "mpp_io.h"
#include "globals_acc.h"

/*******************************************************************************
void set_remap_file( )
*******************************************************************************/
void set_remap_file_acc( int ntiles, const char *mosaic_file, const char *remap_file, Interp_config_acc *interp_acc,
                         unsigned int *opcode, int save_weight_only)
{
  int    i, len, m, fid, vid;
  size_t start[4], nread[4];
  char str1[STRING], tilename[STRING];
  int file_exist;

  if(!remap_file) return;

  for(i=0; i<4; i++) {
    start[i] = 0; nread[i] = 1;
  }
  nread[1] = STRING;

  len = strlen(remap_file);
  if(len >= STRING) mpp_error("setoutput_remap_file(fregrid_util): length of remap_file should be less than STRING");
  if( strcmp(remap_file+len-3, ".nc")==0 ) {
    strncpy(str1, remap_file, len-3);
    str1[len-3] = 0;
  }
  else
    strcpy(str1, remap_file);

  (*opcode) |= WRITE;

  if(ntiles>1) {
     fid = mpp_open(mosaic_file, MPP_READ);
     vid   = mpp_get_varid(fid, "gridtiles");
  }

  for(m=0; m<ntiles; m++) {
    interp_acc[m].file_exist = 0;
    if(ntiles > 1) {
      start[0] = m;
      mpp_get_var_value_block(fid, vid, start, nread, tilename);
      if(strlen(str1) + strlen(tilename) > STRING -5) mpp_error("set_output_remap_file(fregrid_util): length of str1 + "
                                                                "length of tilename should be no greater than STRING-5");
      sprintf(interp_acc[m].remap_file, "%s.%s.nc", str1, tilename);
    }
    else
      sprintf(interp_acc[m].remap_file, "%s.nc", str1);
    /* check interp_acc file to be read (=1) or write ( = 2) */
    if(!save_weight_only && mpp_file_exist(interp_acc[m].remap_file)) {
      (*opcode) |= READ;
      interp_acc[m].file_exist = 1;
    }

  }

  if(ntiles>1) mpp_close(fid);

};/* set_remap_file */

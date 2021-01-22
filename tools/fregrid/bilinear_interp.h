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
 * License along with FRE-NCTools (LICENSE.md).  If not, see
 * <http://www.gnu.org/licenses/>.
 **********************************************************************/
#ifndef BILINEAR_INTERP_H_
#define BILINEAR_INTERP_H_
void setup_bilinear_interp(int ntiles_in, const Grid_config *grid_in, int ntiles_out, const Grid_config *grid_out, 
                      Interp_config *interp, unsigned int opcode, double dlon_in, double dlat_in, double lonbegin, double latbegin);
void do_scalar_bilinear_interp(const Interp_config *interp, int vid, int ntiles_in, const Grid_config *grid_in, const Grid_config *grid_out,
                          const Field_config *field_in, Field_config *field_out, int finer_step, int fill_missing);
void do_vector_bilinear_interp(Interp_config *interp, int vid, int ntiles_in, const Grid_config *grid_in, int ntiles_out, 
			  const Grid_config *grid_out, const Field_config *u_in,  const Field_config *v_in,
			  Field_config *u_out, Field_config *v_out, int finer_step, int fill_missing);
#endif

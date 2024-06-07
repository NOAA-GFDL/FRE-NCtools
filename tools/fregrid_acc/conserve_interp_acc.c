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
#include <string.h>
#include <netcdf.h>
#include <math.h>
#include <time.h>
#include "globals_acc.h"
#include "conserve_interp_acc.h"
#include "interp_utils_acc.h"
#include "create_xgrid_acc.h"
#include "create_xgrid_utils_acc.h"
#include "general_utils_acc.h"
#include "mpp.h"
#include "mpp_io.h"
#include "read_mosaic.h"

/*******************************************************************************
  void setup_conserve_interp
  Setup the interpolation weight for conservative interpolation
*******************************************************************************/
void setup_conserve_interp_acc(int ntiles_input_grid, Grid_config *input_grid, int ntiles_output_grid,
			   Grid_config *output_grid, Xgrid_config *xgrid, unsigned int opcode)
{
  int nlon_input_cells, nlat_input_cells, ncells_input;
  int nlon_output_cells, nlat_output_cells, tile;

  Grid_cells_struct_config *input_grid_cells=
    (Grid_cells_struct_config *)malloc(ntiles_input_grid * sizeof(Grid_cells_struct_config));
  Grid_cells_struct_config output_grid_cells;

  if( opcode & READ) {
    read_remap_file_acc(ntiles_input_grid, ntiles_output_grid, output_grid, input_grid, xgrid, opcode);
    copy_xgrid_to_device_acc(ntiles_input_grid, ntiles_output_grid, xgrid, opcode);
    exit(0);
    return;
  }

  for(int n=0; n<ntiles_output_grid; n++) {

    int nlon_output_cells = output_grid[n].nxc;
    int nlat_output_cells = output_grid[n].nyc;
    int npts_output_grid = (nlon_output_cells+1)*(nlat_output_cells+1);

    xgrid[n].nxcells = 0;

    copy_grid_to_device_acc(npts_output_grid, output_grid[n].latc, output_grid[n].lonc);

    get_grid_cells_struct_acc( nlon_output_cells, nlat_output_cells,
                               output_grid[n].lonc, output_grid[n].latc, &output_grid_cells );

    for(int m=0; m<ntiles_input_grid; m++){

      int nlon_input_cells = input_grid[m].nx;
      int nlat_input_cells = input_grid[m].ny;
      int ncells_input = nlon_input_cells * nlat_input_cells;
      int npts_input_grid = (nlon_input_cells+1)*(nlat_input_cells+1);
      int jlat_overlap_starts=0, jlat_overlap_ends=0, nxcells=0, upbound_nxcells=0;
      int *approx_nxcells_per_ij1, *ij2_start, *ij2_end;

      copy_grid_to_device_acc(npts_input_grid, input_grid[m].latc, input_grid[m].lonc);

      get_skip_cells_acc(nlon_input_cells*nlat_input_cells, &(input_grid_cells[m].skip_cells));

      //get the input grid portion (bounding indices) that overlaps with the output grid in the latitudonal direction.
      get_bounding_indices_acc(nlon_output_cells, nlat_output_cells, nlon_input_cells, nlat_input_cells,
                               output_grid[n].latc, input_grid[m].latc, &jlat_overlap_starts, &jlat_overlap_ends);

      create_upbound_nxcells_arrays_on_device_acc( ncells_input, &approx_nxcells_per_ij1, &ij2_start, &ij2_end);

      upbound_nxcells = get_upbound_nxcells_2dx2d_acc( nlon_input_cells, nlat_input_cells,
                                                       nlon_output_cells, nlat_output_cells,
                                                       jlat_overlap_starts, jlat_overlap_ends,
                                                       input_grid[m].lonc, input_grid[m].latc,
                                                       output_grid[n].lonc, output_grid[n].latc,
                                                       input_grid_cells[m].skip_cells,
                                                       &output_grid_cells,
                                                       approx_nxcells_per_ij1, ij2_start, ij2_end);

      if(opcode & GREAT_CIRCLE) {
        /*nxcells = create_xgrid_great_circle_acc(&nlon_input_cells, &nlat_input_cells,
          &nlon_output_cells, &nlat_output_cells,
          input_grid[m].lonc, input_grid[m].latc,
          output_grid[n].lonc, output_grid[n].latc,
          input_grid_cells[m].skip_cells,
          input_parent_lon_indices, input_parent_lat_indices,
          output_parent_lon_indices, output_parent_lat_indices,
          xcell_area, xcell_centroid_lon, xcell_centroid_lat);
        */
      }
      else {

        if(opcode & CONSERVE_ORDER1) {
          nxcells= create_xgrid_2dx2d_order1_acc( nlon_input_cells, nlat_input_cells,
                                                  nlon_output_cells, nlat_output_cells,
                                                  jlat_overlap_starts, jlat_overlap_ends,
                                                  input_grid[m].lonc, input_grid[m].latc,
                                                  output_grid[n].lonc, output_grid[n].latc,
                                                  upbound_nxcells,
                                                  input_grid_cells[m].skip_cells,
                                                  &output_grid_cells,
                                                  approx_nxcells_per_ij1, ij2_start, ij2_end,
                                                  xgrid[n].per_intile+m);
          xgrid[n].nxcells+=nxcells;

          printf("HERE %d %d\n", upbound_nxcells, nxcells);
          fflush(stdout);
        }
        else if(opcode & CONSERVE_ORDER2) {
          clock_t time_start, time_end;
          time_start = clock();
          nxcells= create_xgrid_2dx2d_order2_acc( nlon_input_cells, nlat_input_cells,
                                                  nlon_output_cells, nlat_output_cells,
                                                  jlat_overlap_starts, jlat_overlap_ends,
                                                  input_grid[m].lonc, input_grid[m].latc,
                                                  output_grid[n].lonc, output_grid[n].latc,
                                                  upbound_nxcells,
                                                  input_grid_cells[m].skip_cells,
                                                  &output_grid_cells,
                                                  approx_nxcells_per_ij1, ij2_start, ij2_end,
                                                  xgrid[n].per_intile+m, input_grid[m].cell_area);
          time_end = clock();
          xgrid[n].nxcells+=nxcells;
          printf("HERE %d %d %f\n", upbound_nxcells, nxcells, (float)(time_end - time_start)/CLOCKS_PER_SEC);
          fflush(stdout);
        }
        else {
          mpp_error("conserve_interp: interp_method should be CONSERVE_ORDER1 or CONSERVE_ORDER2");
        }
      } //conserve_order methods

      free_upbound_nxcells_array_from_all_acc(ncells_input, approx_nxcells_per_ij1, ij2_start, ij2_end);
      free_skip_cells_on_all_acc( ncells_input, input_grid_cells[n].skip_cells);
      delete_grid_from_device_acc(npts_input_grid, input_grid[m].lonc, input_grid[m].latc);

    } //input tile

      /* free the memory */
    free(input_grid_cells);

  }//output tile

  if( opcode & WRITE) write_remap_file(ntiles_output_grid, ntiles_input_grid, output_grid, input_grid, xgrid, opcode);
  if(opcode & CHECK_CONSERVE) check_area_conservation(ntiles_output_grid, ntiles_input_grid, output_grid, xgrid);

  if(mpp_pe() == mpp_root_pe())printf("NOTE: done calculating index and weight for conservative interpolation\n");
  exit(0);
}; /* setup_conserve_interp */


/*******************************************************************************
 void read_remap_file
 Reads in the weight/remap file if provided and copies the data to the device
*******************************************************************************/
void read_remap_file_acc(int ntiles_input_grid, int ntiles_output_grid, Grid_config *output_grid, Grid_config *input_grid,
                         Xgrid_config *xgrid, unsigned int opcode)
{

  int *input_lon=NULL, *input_lat=NULL, *output_lon=NULL, *output_lat=NULL;
  double *xcell_area=NULL, *xcell_centroid_lon=NULL, *xcell_centroid_lat=NULL;

  int nxcells_acc;

  for(int n=0; n<ntiles_output_grid; n++) {
    if( xgrid[n].file_exist ) {
      int *t_in, *ind_acc;
      int fid, vid;

      int nlon_output = output_grid[n].nxc;
      int nxcells = read_mosaic_xgrid_size(xgrid[n].remap_file);
      xgrid[n].nxcells = nxcells;

      t_in  = (int *)malloc(nxcells*sizeof(int));
      input_lon = (int *)malloc(nxcells*sizeof(int));
      input_lat = (int *)malloc(nxcells*sizeof(int));
      output_lat = (int *)malloc(nxcells*sizeof(int));
      output_lon = (int *)malloc(nxcells*sizeof(int));
      xcell_area = (double *)malloc(nxcells*sizeof(double));
      if( opcode & CONSERVE_ORDER2) {
        xcell_centroid_lon = (double *)malloc(nxcells*sizeof(double));
        xcell_centroid_lat = (double *)malloc(nxcells*sizeof(double));
      }

      if(opcode & CONSERVE_ORDER1)
        read_mosaic_xgrid_order1(xgrid[n].remap_file, input_lon, input_lat, output_lon, output_lat, xcell_area);
      else
        read_mosaic_xgrid_order2(xgrid[n].remap_file, input_lon, input_lat, output_lon, output_lat, xcell_area,
                                 xcell_centroid_lon, xcell_centroid_lat);

      //rescale the xgrid area
      for(int i=0; i<xgrid[n].nxcells; i++) xcell_area[i] *= GAREA;

      //read in tile number of input parent cells
      fid = mpp_open(xgrid[n].remap_file, MPP_READ);
      vid = mpp_get_varid(fid, "tile1");
      mpp_get_var_value(fid, vid, t_in);
      mpp_close(fid);

      //tile number starts from 1, not 0, in the weight files
      for(int i=0 ; i<nxcells ; i++) t_in[i]--;

      //get number of nxcells per input tile
      for(int m=0 ; m<ntiles_input_grid ; m++) xgrid[n].per_intile[m].nxcells=0;
      for(int i=0 ; i<nxcells ; i++) xgrid[n].per_intile[ t_in[i] ].nxcells++;

      for(int m=0 ; m<ntiles_input_grid ; m++) {
        nxcells_acc = xgrid[n].per_intile[m].nxcells;
        xgrid[n].per_intile[m].input_parent_cell_indices  = (int *)malloc(nxcells_acc*sizeof(int));
        xgrid[n].per_intile[m].output_parent_cell_indices = (int *)malloc(nxcells_acc*sizeof(int));
        xgrid[n].per_intile[m].xcell_area  = (double *)malloc(nxcells_acc*sizeof(double));
        if(opcode & CONSERVE_ORDER2) {
          xgrid[n].per_intile[m].dcentroid_lon = (double *)malloc(nxcells_acc*sizeof(double));
          xgrid[n].per_intile[m].dcentroid_lat = (double *)malloc(nxcells_acc*sizeof(double));
        }
      }

      ind_acc = (int *)calloc(ntiles_input_grid, sizeof(int));
      for(int i=0 ; i<nxcells ; i++) {
        int itile, ii, nlon_input;
        itile = t_in[i];
        ii=ind_acc[itile];
        nlon_input = input_grid[itile].nxc;
        xgrid[n].per_intile[itile].input_parent_cell_indices[ii] = input_lat[i]*nlon_input + input_lon[i];
        xgrid[n].per_intile[itile].output_parent_cell_indices[ii] = output_lat[i]*nlon_output + output_lon[i];
        xgrid[n].per_intile[itile].xcell_area[ii] = xcell_area[i];
        if( opcode & CONSERVE_ORDER2) {
          xgrid[n].per_intile[itile].dcentroid_lon[ii] = xcell_centroid_lon[i];
          xgrid[n].per_intile[itile].dcentroid_lat[ii] = xcell_centroid_lat[i];
        }
        ind_acc[itile]++;
      }

      free(t_in) ; free(ind_acc);
      free(input_lon) ; free(input_lat); free(output_lon); free(output_lat);
      free(xcell_area); free(xcell_centroid_lon); free(xcell_centroid_lat);

    }//if file exists
  } //ntiles out

  printf("NOTE: Finish reading index and weight for conservative interpolation from file.\n");

} //end read_remap_file

/*******************************************************************************
 void write_remap_file
 write out the traditional remap file.
*******************************************************************************/
void write_remap_file(const int ntiles_out, const int ntiles_in, Grid_config *output_grid,
                      Grid_config *input_grid, Xgrid_config *xgrid, unsigned int opcode)
{

  //copy and pasted from the original code with start and write

  for(int n=0 ; n<ntiles_out ; n++) {

    Xgrid_config *p_xgrid = xgrid+n;
    int nxcells=p_xgrid->nxcells;
    int nlon_cells, ii;

    size_t start[4] = {0,0,0,0}, nwrite[4] = {1, 1, 1, 1};
    int *data_int;
    double *data_double;

    int fid = mpp_open( xgrid[n].remap_file, MPP_WRITE);
    int dim_string = mpp_def_dim(fid, "string", STRING);
    int dim_ncells = mpp_def_dim(fid, "ncells", nxcells);
    int dim_two    = mpp_def_dim(fid, "two", 2);
    int dims[4] = {dim_ncells, dim_two, 0, 0};

    int id_tile1 = mpp_def_var(fid, "tile1", NC_INT, 1, &dim_ncells, 1,
                               "standard_name", "tile_number_in_mosaic1");
    int id_tile1_cell = mpp_def_var(fid, "tile1_cell", NC_INT, 2, dims, 1,
                                    "standard_name", "parent_cell_indices_in_mosaic1");
    int id_tile2_cell = mpp_def_var(fid, "tile2_cell", NC_INT, 2, dims, 1,
                                    "standard_name", "parent_cell_indices_in_mosaic2");
    int id_xgrid_area = mpp_def_var(fid, "xgrid_area", NC_DOUBLE, 1, &dim_ncells, 2,
                                    "standard_name", "exchange_grid_area", "units", "m2");
    int id_tile1_dist = (opcode & CONSERVE_ORDER2) ?
      mpp_def_var(fid, "tile1_distance", NC_DOUBLE, 2, dims, 1,
                  "standard_name", "distance_from_parent1_cell_centroid") : 0 ;

    nwrite[0] = nxcells;
    mpp_end_def(fid);

    data_int = (int *)malloc(nxcells*sizeof(int));

    //update data on host
    for(int m=0 ; m<ntiles_in ; m++) {
      Xinfo_per_input_tile *p_xgrid_for_tile_m = p_xgrid->per_intile+m;
      int m_nxcells = p_xgrid_for_tile_m->nxcells;
#pragma acc update host( p_xgrid_for_tile_m->input_parent_cell_indices[:m_nxcells], \
                         p_xgrid_for_tile_m->output_parent_cell_indices[:m_nxcells], \
                         p_xgrid_for_tile_m->xcell_area[:m_nxcells],    \
                         p_xgrid_for_tile_m->dcentroid_lon[:m_nxcells], \
                         p_xgrid_for_tile_m->dcentroid_lat[:m_nxcells])
    }

    //input tile
    ii = 0;
    for( int m=0 ; m<ntiles_in ; m++ ) {
      Xinfo_per_input_tile *p_xgrid_for_tile_m = p_xgrid->per_intile+m;
      int m_nxcells = p_xgrid_for_tile_m->nxcells;
      for( int i=0 ; i<m_nxcells ; i++ ){
        data_int[ii] = m+1;
        ii++;
      }
    }
    mpp_put_var_value(fid, id_tile1, data_int);

    // i (x, lon) indices of input parent
    ii=0;
    for( int m=0 ; m<ntiles_in ; m++){
      Xinfo_per_input_tile *p_xgrid_for_tile_m = p_xgrid->per_intile+m;
      int m_nxcells = p_xgrid_for_tile_m->nxcells;
      nlon_cells = input_grid[m].nxc;
      for( int i=0 ; i<m_nxcells; i++) {
        data_int[ii] = p_xgrid_for_tile_m->input_parent_cell_indices[i]%nlon_cells+1;
        ii++;
      }
    }
    mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, data_int);

    // i (x, lon) indices of output parent
    ii=0;
    nlon_cells = output_grid[n].nxc;
    for( int m=0 ; m<ntiles_in ; m++){
      Xinfo_per_input_tile *p_xgrid_for_tile_m = p_xgrid->per_intile+m;
      int m_nxcells = p_xgrid_for_tile_m->nxcells;
      for( int i=0 ; i<m_nxcells; i++) {
        data_int[ii] = p_xgrid_for_tile_m->output_parent_cell_indices[i]%nlon_cells+1;
        ii++;
      }
    }
    mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, data_int);

    start[1]=1;

    // j (y, lat) indices of input parent
    ii=0;
    for( int m=0 ; m<ntiles_in ; m++){
      Xinfo_per_input_tile *p_xgrid_for_tile_m = p_xgrid->per_intile+m;
      int m_nxcells = p_xgrid_for_tile_m->nxcells;
      nlon_cells = input_grid[m].nxc;
      for( int i=0 ; i<m_nxcells; i++) {
        data_int[ii] = p_xgrid_for_tile_m->input_parent_cell_indices[i]/nlon_cells+1;
        ii++;
      }
    }
    mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, data_int);

    // j (y, lat) indices of output parent
    ii=0;
    nlon_cells = output_grid[n].nxc;
    for( int m=0 ; m<ntiles_in ; m++){
      Xinfo_per_input_tile *p_xgrid_for_tile_m = p_xgrid->per_intile+m;
      int m_nxcells = p_xgrid_for_tile_m->nxcells;
      for( int i=0 ; i<m_nxcells; i++) {
        data_int[ii] = p_xgrid_for_tile_m->output_parent_cell_indices[i]/nlon_cells+1;
        ii++;
      }
    }
    mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, data_int);

    free(data_int);
    data_double = (double *)malloc(nxcells*sizeof(double));

    // exchange cell area
    ii=0;
    for( int m=0 ; m<ntiles_in ; m++){
      Xinfo_per_input_tile *p_xgrid_for_tile_m = p_xgrid->per_intile+m;
      int m_nxcells = p_xgrid_for_tile_m->nxcells;
      for( int i=0 ; i<m_nxcells; i++) {
        data_double[ii] = p_xgrid_for_tile_m->xcell_area[i];
        ii++;
      }
    }
    mpp_put_var_value(fid, id_xgrid_area, data_double);

    if(opcode & CONSERVE_ORDER2) {
      ii=0; start[1] = 0 ;
      for( int m=0 ; m<ntiles_in ; m++) {
        Xinfo_per_input_tile *p_xgrid_for_tile_m = p_xgrid->per_intile+m;
        int m_nxcells = p_xgrid_for_tile_m->nxcells;
        for( int i=0 ; i<m_nxcells ; i++) {
          data_double[ii] = p_xgrid_for_tile_m->dcentroid_lon[i];
          ii++;
        }
      }
      mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, data_double);

      ii=0; start[1] = 1 ;
      for( int m=0 ; m<ntiles_in ; m++) {
        Xinfo_per_input_tile *p_xgrid_for_tile_m = p_xgrid->per_intile+m;
        int m_nxcells = p_xgrid_for_tile_m->nxcells;
        for( int i=0 ; i<m_nxcells ; i++) {
          data_double[ii] = p_xgrid_for_tile_m->dcentroid_lat[i];
          ii++;
        }
      }
      mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, data_double);
    }
    free(data_double);

    mpp_close(fid);

  }//ntiles_out

}

void check_area_conservation(const int ntiles_output_grid, const int ntiles_input_grid, Grid_config *output_grid,
                             Xgrid_config *xgrid)
{

  for(int n=0; n<ntiles_output_grid; n++) {

    int nlon_cells = output_grid[n].nxc;
    int nlat_cells = output_grid[n].nyc;
    int ncells = nlon_cells*nlat_cells;

    int max_ij=0;
    double max_ratio=0.0, ratio_change=0.0, max_area=0.0;
    double *recomputed_output_area = (double *)calloc(ncells,sizeof(double));

    /* sum over exchange grid to get the area of output grid cells*/
    for(int m=0; m<ntiles_input_grid; m++) {
      Xinfo_per_input_tile *m_xinfo = xgrid[n].per_intile+m;
      int nxcells = m_xinfo->nxcells;
      for(int i=0; i<nxcells; i++) {
        int ii = m_xinfo->output_parent_cell_indices[i];
        if(ii<0) printf("SOMETHING WRONG, %d\n", ii);
        if(ii>ncells) printf("SOMETHING WRONG, %d\n", ii);
        recomputed_output_area[ii] += m_xinfo->xcell_area[i];
      }
    }

    /* compare actual area and recomputed_output_area */
    for(int ij=0 ; ij<ncells ; ij++) {
      double actual_cell_area=output_grid[n].cell_area[ij];
      ratio_change = fabs(actual_cell_area-recomputed_output_area[ij])/actual_cell_area;
      if(ratio_change > max_ratio) {
        max_ratio = ratio_change;
        max_area = recomputed_output_area[ij];
        max_ij = ij;
      }
      if( ratio_change > 1.e-4 ) {
        printf("(i,j)=(%d,%d), change = %g, area1=%g, recomputed_output_area=%g\n",
               ij%nlon_cells, ij/nlon_cells, ratio_change, output_grid[n].cell_area[ij],recomputed_output_area[ij]);
      }
    }
    printf("The maximum ratio change at (%d,%d) = %g, area1=%g, recomputed_output_area=%g\n",
           max_ij%nlon_cells, max_ij/nlon_cells, max_ratio, output_grid[n].cell_area[max_ij], max_area);

    free(recomputed_output_area);
  }//for each output tile

}

/*******************************************************************************
 void do_scalar_conserve_interp( )
 doing conservative interpolation
*******************************************************************************/
void do_scalar_conserve_interp_acc(Xgrid_config *xgrid, int varid, int ntiles_input_grid, const Grid_config *input_grid,
                                   int ntiles_output_grid, const Grid_config *output_grid, const Field_config *field_in,
                                   Field_config *field_out, unsigned int opcode, int nz)
{
  int nx1, ny1, nx2, ny2, i1, j1, i2, j2, tile, n, m, i, j, n1, n2;
  int k, n0;
  int has_missing, halo, interp_method;
  int weight_exist;
  int cell_measures, cell_methods;
  double area, missing, di, dj, area_missing;
  double *out_area;
  int    *out_miss;
  double gsum_out;
  int monotonic;
  int target_grid;
  Monotone_config *monotone_data;

  gsum_out = 0;
  interp_method = field_in->var[varid].interp_method;
  halo = 0;
  monotonic = 0;
  if(interp_method == CONSERVE_ORDER2) {
    halo = 1;
    monotonic = opcode & MONOTONIC;
  }

  area_missing = field_in->var[varid].area_missing;
  has_missing = field_in->var[varid].has_missing;
  weight_exist = input_grid[0].weight_exist;
  cell_measures = field_in->var[varid].cell_measures;
  cell_methods = field_in->var[varid].cell_methods;
  target_grid = opcode & TARGET;
  if( field_in->var[varid].use_volume ) target_grid = 0;

  missing = -MAXVAL;
  if(has_missing) missing = field_in->var[varid].missing;

  if( nz>1 && has_missing ) mpp_error("conserve_interp: has_missing should be false when nz > 1");
  if( nz>1 && cell_measures ) mpp_error("conserve_interp: cell_measures should be false when nz > 1");
  if( nz>1 && cell_methods == CELL_METHODS_SUM ) mpp_error("conserve_interp: cell_methods should not be sum when nz > 1");
  /*  if( nz>1 && monotonic ) mpp_error("conserve_interp: monotonic should be false when nz > 1"); */

  if(monotonic) monotone_data = (Monotone_config *)malloc(ntiles_input_grid*sizeof(Monotone_config));

  for(m=0; m<ntiles_output_grid; m++) {
    nx2 = output_grid[m].nxc;
    ny2 = output_grid[m].nyc;
    out_area = (double *)malloc(nx2*ny2*nz*sizeof(double));
    out_miss = (int *)malloc(nx2*ny2*nz*sizeof(int));
    for(i=0; i<nx2*ny2*nz; i++) {
      field_out[m].data[i] = 0.0;
      out_area[i] = 0.0;
      out_miss[i] = 0;
    }
    if(interp_method == CONSERVE_ORDER1) {
      if(has_missing) {
	for(n=0; n<xgrid[m].nxcells; n++) {
	  i2   = xgrid[m].output_parent_lon_indices[n];
	  j2   = xgrid[m].output_parent_lat_indices[n];
	  i1   = xgrid[m].input_parent_lon_indices [n];
	  j1   = xgrid[m].input_parent_lat_indices [n];
	  tile = xgrid[m].input_parent_tile [n];
	  area = xgrid[m].xcell_area [n];
	  nx1  = input_grid[tile].nx;
	  ny1  = input_grid[tile].ny;
          if(weight_exist) area *= input_grid[tile].weight[j1*nx1+i1];
	  n1 = j1*nx1+i1;
	  n0 = j2*nx2+i2;

	  if( field_in[tile].data[n1] != missing ) {
	    if( cell_methods == CELL_METHODS_SUM )
	      area /= input_grid[tile].cell_area[n1];
            else if( cell_measures ) {
	      if(field_in[tile].area[n1] == area_missing) {
	         printf("name=%s,tile=%d,i1,j1=%d,%d,i2,j2=%d,%d\n",field_in->var[varid].name,tile,i1,j1,i2,j2);
	         mpp_error("conserve_interp: data is not missing but area is missing");
	      }
	      area *= (field_in[tile].area[n1]/input_grid[tile].cell_area[n1]);
	    }
  	    field_out[m].data[n0] += (field_in[tile].data[n1]*area);
            out_area[n0] += area;
	    out_miss[n0] = 1;
          }
        }
      }
      else {
	for(n=0; n<xgrid[m].nxcells; n++) {
	  i2   = xgrid[m].output_parent_lon_indices[n];
	  j2   = xgrid[m].output_parent_lat_indices[n];
	  i1   = xgrid[m].input_parent_lon_indices [n];
	  j1   = xgrid[m].input_parent_lat_indices [n];
	  tile = xgrid[m].input_parent_tile [n];
	  area = xgrid[m].xcell_area [n];
	  nx1  = input_grid[tile].nx;
	  ny1  = input_grid[tile].ny;
	  if(weight_exist) area *= input_grid[tile].weight[j1*nx1+i1];
	  for(k=0; k<nz; k++) {
	    n1 = k*nx1*ny1 + j1*nx1+i1;
	    n0 = k*nx2*ny2 + j2*nx2+i2;
	    if(  cell_methods == CELL_METHODS_SUM )
	      area /= input_grid[tile].cell_area[n1];
            else if( cell_measures )
	      area *= (field_in[tile].area[n1]/input_grid[tile].cell_area[n1]);
  	    field_out[m].data[n0] += (field_in[tile].data[n1]*area);
	    out_area[n0] += area;
	    out_miss[n0] = 1;
	  }
	}
      }
    }
    else if(monotonic) {
      int ii, jj;
      double f_bar;
      double *xdata;
      for(n=0; n<ntiles_input_grid; n++) {
	nx1 =  input_grid[n].nx;
	ny1 =  input_grid[n].ny;
	monotone_data[n].f_bar_max = (double *)malloc(nx1*ny1*sizeof(double));
	monotone_data[n].f_bar_min = (double *)malloc(nx1*ny1*sizeof(double));
	monotone_data[n].f_max     = (double *)malloc(nx1*ny1*sizeof(double));
	monotone_data[n].f_min     = (double *)malloc(nx1*ny1*sizeof(double));
	for(j=0; j<ny1; j++) for(i=0; i<nx1; i++) {
	  n1 = j*nx1+i;

	  monotone_data[n].f_bar_max[n1] = -MAXVAL;
	  monotone_data[n].f_bar_min[n1] = MAXVAL;
	  monotone_data[n].f_max[n1]     = -MAXVAL;
	  monotone_data[n].f_min[n1]     = MAXVAL;
	  n1 = j*nx1+i;
	  for(jj=j-1; jj<=j+1; jj++) for(ii=i-1; ii<=i+1; ii++) {
	    n2 = (jj+1)*(nx1+2)+ii+1;
	    if( field_in[n].data[n2] != missing ) {
	      if( field_in[n].data[n2] > monotone_data[n].f_bar_max[n1] ) monotone_data[n].f_bar_max[n1] = field_in[n].data[n2];
	      if( field_in[n].data[n2] < monotone_data[n].f_bar_min[n1] ) monotone_data[n].f_bar_min[n1] = field_in[n].data[n2];
	    }
	  }
	}
      }

      xdata = (double *)malloc(xgrid[m].nxcells*sizeof(double));
      for(n=0; n<xgrid[m].nxcells; n++) {
	i1   = xgrid[m].input_parent_lon_indices [n];
	j1   = xgrid[m].input_parent_lat_indices [n];
	di   = xgrid[m].dcentroid_lon[n];
	dj   = xgrid[m].dcentroid_lat[n];
	tile = xgrid[m].input_parent_tile [n];
	n1 = j1*nx1+i1;
        n2 = (j1+1)*(nx1+2)+i1+1;
	if( field_in[tile].data[n2] != missing ) {
	  if( field_in[tile].grad_mask[n1] ) { /* use zero gradient */
	    xdata[n] = field_in[tile].data[n2];
	  }
	  else {
	    xdata[n] = field_in[tile].data[n2]+field_in[tile].grad_x[n1]*di+field_in[tile].grad_y[n1]*dj;
	  }
	  if(monotonic) {
	    if( xdata[n] > monotone_data[tile].f_max[n1]) monotone_data[tile].f_max[n1] = xdata[n];
	    if( xdata[n] < monotone_data[tile].f_min[n1]) monotone_data[tile].f_min[n1] = xdata[n];
	  }
	}
	else
	  xdata[n] = missing;
      }

      /* get the global f_max and f_min */
      if(mpp_npes() >1) {
	for(n=0; n<ntiles_input_grid; n++) {
	  mpp_min_double(nx1*ny1, monotone_data[n].f_min);
	  mpp_max_double(nx1*ny1, monotone_data[n].f_max);
	}
      }

      /* adjust the exchange grid cell data to make it monotonic */
      for(n=0; n<xgrid[m].nxcells; n++) {
	i1   = xgrid[m].input_parent_lon_indices [n];
	j1   = xgrid[m].input_parent_lat_indices [n];
	tile = xgrid[m].input_parent_tile [n];
	n1 = j1*nx1+i1;
	n2 = (j1+1)*(nx1+2)+i1+1;
	f_bar = field_in[tile].data[n2];
	if(xdata[n] == missing) continue;

	if( monotone_data[tile].f_max[n1] > monotone_data[tile].f_bar_max[n1] ) {
	  /* z1l: Due to truncation error, we might get xdata[n] > f_bar_max[n1]. So
	     we allow some tolerance. What is the suitable tolerance? */
	  xdata[n] = f_bar + ((xdata[n]-f_bar)/(monotone_data[tile].f_max[n1]-f_bar))
	    * (monotone_data[tile].f_bar_max[n1]-f_bar);
	  if( xdata[n] > monotone_data[tile].f_bar_max[n1]) {
	    if(xdata[n] - monotone_data[tile].f_bar_max[n1] < TOLERANCE ) xdata[n] = monotone_data[tile].f_bar_max[n1];
	    if( xdata[n] > monotone_data[tile].f_bar_max[n1]) {
	      printf(" n = %d, n1 = %d, xdata = %f, f_bar_max=%f\n", n, n1, xdata[n], monotone_data[tile].f_bar_max[n1]);
	      mpp_error(" xdata is greater than f_bar_max ");
	    }
	  }
	}
	else if( monotone_data[tile].f_min[n1] < monotone_data[tile].f_bar_min[n1] ) {
	  /* z1l: Due to truncation error, we might get xdata[n] < f_bar_min[n1]. So
	     we allow some tolerance. What is the suitable tolerance? */
	  xdata[n] = f_bar + ((xdata[n]-f_bar)/(monotone_data[tile].f_min[n1]-f_bar)) * (monotone_data[tile].f_bar_min[n1]-f_bar);
	  if( xdata[n] < monotone_data[tile].f_bar_min[n1]) {
	    if(monotone_data[tile].f_bar_min[n1] - xdata[n]< TOLERANCE ) xdata[n] = monotone_data[tile].f_bar_min[n1];
	    if( xdata[n] < monotone_data[tile].f_bar_min[n1]) {
	      printf(" n = %d, n1 = %d, xdata = %f, f_bar_min=%f\n", n, n1, xdata[n], monotone_data[tile].f_bar_min[n1]);
	      mpp_error(" xdata is less than f_bar_min ");
	    }
	  }
	}
      }
      for(n=0; n<ntiles_input_grid; n++) {
	free(monotone_data[n].f_bar_max);
	free(monotone_data[n].f_bar_min);
	free(monotone_data[n].f_max);
	free(monotone_data[n].f_min);
      }

      /* remap onto destination grid */
      for(n=0; n<xgrid[m].nxcells; n++) {
	i2   = xgrid[m].output_parent_lon_indices[n];
	j2   = xgrid[m].output_parent_lat_indices[n];
	i1   = xgrid[m].input_parent_lon_indices [n];
	j1   = xgrid[m].input_parent_lat_indices [n];
	tile = xgrid[m].input_parent_tile [n];
	area = xgrid[m].xcell_area [n];
	if(xdata[n] == missing) continue;
	if(weight_exist) area *= input_grid[tile].weight[j1*nx1+i1];
	n1 = j1*nx1+i1;
	n0 = j2*nx2+i2;
	if( cell_methods == CELL_METHODS_SUM )
	  area /= input_grid[tile].cell_area[n1];
	else if( cell_measures )
	  area *= (field_in[tile].area[n1]/input_grid[tile].cell_area[n1]);
	field_out[m].data[n0] += xdata[n]*area;
	out_area[n0] += area;
      }
      free(xdata);
    }
    else {
      if(has_missing) {
	for(n=0; n<xgrid[m].nxcells; n++) {
	  i2   = xgrid[m].output_parent_lon_indices[n];
	  j2   = xgrid[m].output_parent_lat_indices[n];
	  i1   = xgrid[m].input_parent_lon_indices [n];
	  j1   = xgrid[m].input_parent_lat_indices [n];
	  di   = xgrid[m].dcentroid_lon[n];
	  dj   = xgrid[m].dcentroid_lat[n];
	  tile = xgrid[m].input_parent_tile [n];
	  area = xgrid[m].xcell_area [n];
	  nx1  = input_grid[tile].nx;
	  ny1  = input_grid[tile].ny;

	  if(weight_exist) area *= input_grid[tile].weight[j1*nx1+i1];
	  n2 = (j1+1)*(nx1+2)+i1+1;
	  n0 = j2*nx2+i2;
	  if( field_in[tile].data[n2] != missing ) {
	    n1 = j1*nx1+i1;
	    n0 = j2*nx2+i2;
            if( cell_methods == CELL_METHODS_SUM )
	      area /= input_grid[tile].cell_area[n1];
            else if( cell_measures ) {
	      if(field_in[tile].area[n1] == area_missing) {
                printf("name=%s,tile=%d,i1,j1=%d,%d,i2,j2=%d,%d\n",field_in->var[varid].name,tile,i1,j1,i2,j2);
	        mpp_error("conserve_interp: data is not missing but area is missing");
              }
	      area *= (field_in[tile].area[n1]/input_grid[tile].cell_area[n1]);
            }
	    if(field_in[tile].grad_mask[n1]) { /* use zero gradient */
	      field_out[m].data[n0] += field_in[tile].data[n2]*area;
	    }
	    else {
	      field_out[m].data[n0] += (field_in[tile].data[n2]+field_in[tile].grad_x[n1]*di
					+field_in[tile].grad_y[n1]*dj)*area;
	    }
	    out_area[n0] += area;
	    out_miss[n0] = 1;
	  }
	}
      }
      else {
	for(n=0; n<xgrid[m].nxcells; n++) {
	  i2   = xgrid[m].output_parent_lon_indices[n];
	  j2   = xgrid[m].output_parent_lat_indices[n];
	  i1   = xgrid[m].input_parent_lon_indices [n];
	  j1   = xgrid[m].input_parent_lat_indices [n];
	  di   = xgrid[m].dcentroid_lon[n];
	  dj   = xgrid[m].dcentroid_lat[n];
	  tile = xgrid[m].input_parent_tile [n];
	  area = xgrid[m].xcell_area [n];

	  nx1  = input_grid[tile].nx;
	  ny1  = input_grid[tile].ny;
	  if(weight_exist) area *= input_grid[tile].weight[j1*nx1+i1];
	  for(k=0; k<nz; k++) {
	    n0 = k*nx2*ny2 + j2*nx2+i2;
	    n1 = k*nx1*ny1+j1*nx1+i1;
	    n2 = k*(nx1+2)*(ny1+2)+(j1+1)*(nx1+2)+i1+1;
	    if( cell_methods == CELL_METHODS_SUM )
	      area /= input_grid[tile].cell_area[n1];
	    else if( cell_measures )
	      area *= (field_in[tile].area[n1]/input_grid[tile].cell_area[n1]);
	    field_out[m].data[n0] += (field_in[tile].data[n2]+field_in[tile].grad_x[n1]*di
				      +field_in[tile].grad_y[n1]*dj)*area;
	    out_area[n0] += area;
            out_miss[n0] = 1;
	  }
	}
      }
    }

    if(opcode & CHECK_CONSERVE) {
      for(i=0; i<nx2*ny2*nz; i++) {
	if(out_area[i] > 0) gsum_out += field_out[m].data[i];
      }
    }

    if ( cell_methods == CELL_METHODS_SUM ) {
      for(i=0; i<nx2*ny2*nz; i++) {
        if(out_area[i] == 0) {
          if(out_miss[i] == 0)
            for(k=0; k<nz; k++) field_out[m].data[k*nx2*ny2+i] = missing;
          else
            for(k=0; k<nz; k++) field_out[m].data[k*nx2*ny2+i] = 0.0;
        }
      }
    }
    else {
      for(i=0; i<nx2*ny2*nz; i++) {
	if(out_area[i] > 0)
	  field_out[m].data[i] /= out_area[i];
	else if(out_miss[i] == 1)
	  field_out[m].data[i] = 0.0;
	else
	  field_out[m].data[i] = missing;
      }

      if( (target_grid) ) {
	for(i=0; i<nx2*ny2; i++) out_area[i] = 0.0;
	for(n=0; n<xgrid[m].nxcells; n++) {
	  i2   = xgrid[m].output_parent_lon_indices[n];
	  j2   = xgrid[m].output_parent_lat_indices[n];
	  i1   = xgrid[m].input_parent_lon_indices [n];
	  j1   = xgrid[m].input_parent_lat_indices [n];
	  tile = xgrid[m].input_parent_tile [n];
	  area = xgrid[m].xcell_area [n];
	  nx1  = input_grid[tile].nx;
	  ny1  = input_grid[tile].ny;
	  n0 = j2*nx2+i2;
	  n1 = j1*nx1+i1;
	  if(cell_measures )
	    out_area[n0] += (area*field_in[tile].area[n1]/input_grid[tile].cell_area[n1]);
	  else
	    out_area[n0] += area;
	}
	for(i=0; i<nx2*ny2*nz; i++) {
	  if(field_out[m].data[i] != missing) {
	    i2 = i%(nx2*ny2);
	    field_out[m].data[i] *=  (out_area[i2]/output_grid[m].cell_area[i2]);
	  }
	}
      }
    }

    free(out_area);
    free(out_miss);
  }


  /* conservation check if needed */
  if(opcode & CHECK_CONSERVE) {
    double gsum_in, dd;
    gsum_in = 0;
    for(n=0; n<ntiles_input_grid; n++) {

      nx1  = input_grid[n].nx;
      ny1  = input_grid[n].ny;


      if( cell_measures ) {
        for(j=0; j<ny1; j++) for(i=0; i<nx1; i++) {
  	  dd = field_in[n].data[(j+halo)*(nx1+2*halo)+i+halo];
	  if(dd != missing) gsum_in += dd*field_in[n].area[j*nx1+i];
        }
      }
      else if ( cell_methods == CELL_METHODS_SUM ) {
        for(j=0; j<ny1; j++) for(i=0; i<nx1; i++) {
	  dd = field_in[n].data[(j+halo)*(nx1+2*halo)+i+halo];
	  if(dd != missing) gsum_in += dd;
        }
      }
      else {
        for(k=0; k<nz; k++) for(j=0; j<ny1; j++) for(i=0; i<nx1; i++) {
	  dd = field_in[n].data[k*(nx1+2*halo)*(ny1+2*halo)+(j+halo)*(nx1+2*halo)+i+halo];
	  if(dd != missing) gsum_in += dd*input_grid[n].cell_area[j*nx1+i];
        }
      }
    }
    mpp_sum_double(1, &gsum_out);

    if(mpp_pe() == mpp_root_pe()) printf("the flux(data*area) sum of %s: input = %g, output = %g, diff = %g. \n",
					 field_in->var[varid].name, gsum_in, gsum_out, gsum_out-gsum_in);

  }


}; /* do_scalar_conserve_interp */


/*******************************************************************************
 void do_vector_conserve_interp( )
 doing conservative interpolation
*******************************************************************************/
void do_vector_conserve_interp_acc(Xgrid_config *xgrid, int varid, int ntiles_input_grid, const Grid_config *input_grid, int ntiles_output_grid,
                               const Grid_config *output_grid, const Field_config *u_in,  const Field_config *v_in,
                               Field_config *u_out, Field_config *v_out, unsigned int opcode)
{
  int          nx1, ny1, nx2, ny2, i1, j1, i2, j2, tile, n, m, i;
  double       area, missing, tmp_x, tmp_y;
  double       *out_area;

  missing = u_in->var[varid].missing;
  /* first rotate input data */
  for(n = 0; n < ntiles_input_grid; n++) {
    if(input_grid[n].rotate) {
      nx1 = input_grid[n].nx;
      ny1 = input_grid[n].ny;
      for(i=0; i<nx1*ny1; i++) {
	tmp_x = u_in[n].data[i];
	tmp_y = v_in[n].data[i];
	if( tmp_x != missing && tmp_y != missing) {
	  u_in[n].data[i] = tmp_x * input_grid[n].cosrot[i] - tmp_y * input_grid[n].sinrot[i];
	  v_in[n].data[i] = tmp_x * input_grid[n].sinrot[i] + tmp_y * input_grid[n].cosrot[i];
	}
      }
    }
  }

  for(m=0; m<ntiles_output_grid; m++) {
    nx2 = output_grid[m].nxc;
    ny2 = output_grid[m].nyc;
    out_area = (double *)malloc(nx2*ny2*sizeof(double));

    for(i=0; i<nx2*ny2; i++) {
      u_out[m].data[i] = 0.0;
      v_out[m].data[i] = 0.0;
    }
    for(i=0; i<nx2*ny2; i++) out_area[i] = 0.0;

    for(n=0; n<xgrid[m].nxcells; n++) {
      i2   = xgrid[m].output_parent_lon_indices[n];
      j2   = xgrid[m].output_parent_lat_indices[n];
      i1   = xgrid[m].input_parent_lon_indices [n];
      j1   = xgrid[m].input_parent_lat_indices [n];
      tile = xgrid[m].input_parent_tile [n];
      area = xgrid[m].xcell_area [n];
      nx1  = input_grid[tile].nx;
      ny1  = input_grid[tile].ny;
      tmp_x = u_in[tile].data[j1*nx1+i1];
      tmp_y = v_in[tile].data[j1*nx1+i1];
      if( tmp_x != missing && tmp_y != missing ) {
	u_out[m].data[j2*nx2+i2] += u_in[tile].data[j1*nx1+i1]*area;
	v_out[m].data[j2*nx2+i2] += v_in[tile].data[j1*nx1+i1]*area;
	out_area[j2*nx2+i2] += area;
      }
    }
    if(opcode & TARGET) {
      for(i=0; i<nx2*ny2; i++) {
	if(out_area[i] > 0) {
	  u_out[m].data[i] /= output_grid[m].area[i];
	  v_out[m].data[i] /= output_grid[m].area[i];
	}
	else {
	  u_out[m].data[i] = missing;
	  v_out[m].data[i] = missing;
	}
      }
    }
    else {
      for(i=0; i<nx2*ny2; i++) {
	if(out_area[i] > 0) {
	  u_out[m].data[i] /= out_area[i];
	  v_out[m].data[i] /= out_area[i];
	}
	else {
	  u_out[m].data[i] = missing;
	  v_out[m].data[i] = missing;
	}
      }
    }
    /* rotate the data if needed */
    if(output_grid[m].rotate) {
      for(i=0; i<nx2*ny2; i++) {
	tmp_x = u_out[m].data[i];
	tmp_y = v_out[m].data[i];
	if( tmp_x != missing && tmp_y != missing) {
	  u_out[m].data[i] =  tmp_x * output_grid[m].cosrot[i] + tmp_y * output_grid[m].sinrot[i];
	  v_out[m].data[i] = -tmp_x * output_grid[m].sinrot[i] + tmp_y * output_grid[m].cosrot[i];
	}
      }
    }
    free(out_area);
  }

}; /* do_vector_conserve_interp */

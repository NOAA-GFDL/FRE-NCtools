TESTS_ENVIRONMENT = top_srcdir=$(abs_top_srcdir); export top_srcdir;\
    for d in $$( find $(abs_top_builddir)/postprocessing $(abs_top_builddir)/tools -mindepth 1 -not -path '*/\.*' -type d ); do \
      PATH="$$d"'$(PATH_SEPARATOR)'"$$PATH"; \
    done;


LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
             $(top_srcdir)/tap-driver.sh

TESTS = 01-mppnccombine.sh \
        check_programs_exist.sh \
        test_grid_coupled_model.sh \
        test_grid_coupled_nest.sh \
        test_remap_c48_regular.sh \
        test_regrid_extrap.sh \
        test_runoff_regrid.sh \
        test_check_mask.sh \
        test_check_mask_MOM.sh 

EXTRA_DIST = $(TESTS) \
             mppnccombine/mppnccombine.ncl.0000 \
             mppnccombine/mppnccombine.ncl.0001 \
             mppnccombine/mppnccombine.ncl.0002 \
             mppnccombine/mppnccombine.ncl.0003 \
             grid_coupled_model/OCCAM_p5degree.nc \
             grid_coupled_model/C48_mosaic.nc \
             grid_coupled_model/C48_grid.tile1.nc \
             grid_coupled_model/C48_grid.tile2.nc \
             grid_coupled_model/C48_grid.tile3.nc \
             grid_coupled_model/C48_grid.tile4.nc \
             grid_coupled_model/C48_grid.tile5.nc \
             grid_coupled_model/C48_grid.tile6.nc \
             grid_coupled_model/19800101.atmos_daily.tile1.nc \
             grid_coupled_model/19800101.atmos_daily.tile2.nc \
             grid_coupled_model/19800101.atmos_daily.tile3.nc \
             grid_coupled_model/19800101.atmos_daily.tile4.nc \
             grid_coupled_model/19800101.atmos_daily.tile5.nc \
             grid_coupled_model/19800101.atmos_daily.tile6.nc \
             regrid_extrap/ocean_hgrid.nc \
             regrid_extrap/ocean_mosaic.nc \
             regrid_extrap/ocean_vgrid.nc \
             regrid_extrap/WOA09_ann_theta.nc \
             runoff_regrid \
             river_regrid \
             checkmask_baltic1

clean-local:
	$(RM) *.nc

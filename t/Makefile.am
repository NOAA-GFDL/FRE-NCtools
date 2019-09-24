TESTS_ENVIRONMENT = top_srcdir=$(abs_top_srcdir); export top_srcdir;\
    for d in $$( find $(abs_top_builddir)/postprocessing $(abs_top_builddir)/tools -mindepth 1 -not -path '*/\.*' -type d ); do \
      PATH="$$d"'$(PATH_SEPARATOR)'"$$PATH"; \
    done;


LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
             $(top_srcdir)/tap-driver.sh

TESTS = 01-mppnccombine.sh \
        check_programs_exist.sh \
        test_grid_coupled_model.sh \
        test_grid_coupled_nest.sh 

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
             grid_coupled_model/C48_grid.tile6.nc

clean-local:
	$(RM) *.nc

TESTS_ENVIRONMENT = top_srcdir=$(abs_top_srcdir); export top_srcdir;\
    for d in $$( find $(abs_top_builddir)/postprocessing $(abs_top_builddir)/tools -mindepth 1 -not -path '*/\.*' -type d ); do \
      PATH="$$d"'$(PATH_SEPARATOR)'"$$PATH"; \
    done;


LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
             $(top_srcdir)/tap-driver.sh

TESTS = Test02-mppnccombine.sh \
        check_programs_exist.sh \
        Test03-grid_coupled_model.sh \
        Test04-grid_coupled_nest.sh \
        Test05-remap_c48_regular.sh \
        Test06-regrid_extrap.sh \
        Test07-runoff_regrid.sh \
        Test08-check_mask.sh \
        Test09-check_mask_MOM.sh \
        Test10-remap_land_res.sh \
        Test11-make_regional_mosaic.sh \
        Test12-mppnscatter.sh \
        test_make_quick_mosaic.sh \
        test_remap_ocean.sh \
        test_regrid_land.sh

EXTRA_DIST = $(TESTS) \
             Test02-input\
             Test03-input \
             Test06-input \
             Test07-input \
             river_regrid \
             Test08-input \
             Test10-input \
             Test12-input \
             regrid_with_cell_measures

clean-local:
	$(RM) *.nc

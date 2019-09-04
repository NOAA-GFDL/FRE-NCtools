TESTS_ENVIRONMENT = \
    for d in $$( find $(abs_top_builddir)/postprocessing $(abs_top_builddir)/tools -mindepth 1 -not -path '*/\.*' -type d ); do \
      PATH="$$d"'$(PATH_SEPARATOR)'"$$PATH"; \
    done;


TEST_LOG_DRIVER = $(PROVE)

TESTS = 01-make_hgrid.sh \
        01-make_vgrid.sh

EXTRA_DIST = $(TESTS)

clean-local:
	$(RM) *.nc

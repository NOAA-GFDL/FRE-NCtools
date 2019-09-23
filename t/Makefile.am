TESTS_ENVIRONMENT = top_srcdir=$(abs_top_srcdir); export top_srcdir;\
    for d in $$( find $(abs_top_builddir)/postprocessing $(abs_top_builddir)/tools -mindepth 1 -not -path '*/\.*' -type d ); do \
      PATH="$$d"'$(PATH_SEPARATOR)'"$$PATH"; \
    done;


LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
             $(top_srcdir)/tap-driver.sh

TESTS = 01-mppnccombine.sh \
        check_programs_exist.sh 

EXTRA_DIST = $(TESTS) \
             mppnccombine/mppnccombine.ncl.0000 \
             mppnccombine/mppnccombine.ncl.0001 \
             mppnccombine/mppnccombine.ncl.0002 \
             mppnccombine/mppnccombine.ncl.0003


clean-local:
	$(RM) *.nc

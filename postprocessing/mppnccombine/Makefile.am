bin_PROGRAMS = mppnccombine

mppnccombine_SOURCES = mppnccombine.c
mppnccombine_CFLAGS = $(NETCDF_CFLAGS) $(AM_CFLAGS)
mppnccombine_LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS) -lm

clean-local:
	$(RM) *.mod

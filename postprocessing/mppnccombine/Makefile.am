bin_PROGRAMS = mppnccombine

mppnccombine_SOURCES = mppnccombine.c
mppnccombine_CFLAGS = $(BASE_CFLAGS) $(AM_CFLAGS)
mppnccombine_LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

clean-local:
	$(RM) *.mod

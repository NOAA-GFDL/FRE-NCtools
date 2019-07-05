bin_PROGRAMS = combine_blobs

combine_blobs_SOURCES = combine_blobs.c
combine_blobs_CFLAGS = $(BASE_CFLAGS) $(AM_CFLAGS)
combine_blobs_LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

clean-local:
	$(RM) *.mod

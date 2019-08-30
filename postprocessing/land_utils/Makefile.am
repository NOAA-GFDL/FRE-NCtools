bin_PROGRAMS = combine-ncc decompress-ncc is-compressed
bin_SCRIPTS = combine-lumip-vars.py

#combine_ncc_FC = $(FC)
combine_ncc_SOURCES = nfu.F90 nfu_compress.F90 combine-ncc.F90
combine_ncc_FCFLAGS = $(NETCDF_FCFLAGS) $(BASE_FCFLAGS) $(R8_FCFLAGS) $(AM_FCFLAGS)
combine_ncc_LDADD = $(NETCDF_FCLDFLAGS) $(NETCDF_FCLIBS)

decompress_ncc_SOURCES = nfu.F90 nfu_compress.F90 decompress-ncc.F90
decompress_ncc_FCFLAGS = $(NETCDF_FCFLAGS) $(BASE_FCFLAGS) $(R8_FCFLAGS) $(AM_FCFLAGS)
decompress_ncc_LDADD = $(NETCDF_FCLDFLAGS) $(NETCDF_FCLIBS)

is_compressed_SOURCES = is-compressed.c
is_compressed_CFLAGS = $(NETCDF_CFLAGS) $(BASE_CFLAGS) $(AM_CFLAGS)
is_compressed_LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

clean-local:
	$(RM) *.mod

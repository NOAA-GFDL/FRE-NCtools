bin_PROGRAMS = combine-ncc decompress-ncc is-compressed scatter-ncc
dist_bin_SCRIPTS = combine-lumip-vars.py

AM_FCFLAGS =  $(NETCDF_FCFLAGS)

combine_ncc_SOURCES = nfu.F90 nfu_compress.F90 combine-ncc.F90
combine_ncc_FCFLAGS = $(FC_REAL_KIND8_FLAG) $(AM_FCFLAGS)
combine_ncc_LDADD = $(NETCDF_FCLDFLAGS) $(NETCDF_FCLIBS)

decompress_ncc_SOURCES = nfu.F90 nfu_compress.F90 decompress-ncc.F90
decompress_ncc_FCFLAGS = $(FC_REAL_KIND8_FLAG) $(AM_FCFLAGS)
decompress_ncc_LDADD = $(NETCDF_FCLDFLAGS) $(NETCDF_FCLIBS)

scatter_ncc_SOURCES = nfu.F90 nfu_compress.F90 scatter-ncc.F90
scatter_ncc_FCFLAGS = $(FC_REAL_KIND8_FLAG) $(AM_FCFLAGS)
scatter_ncc_LDADD = $(NETCDF_FCLDFLAGS) $(NETCDF_FCLIBS)

is_compressed_SOURCES = is-compressed.c
is_compressed_CFLAGS = $(NETCDF_CFLAGS) $(AM_CFLAGS)
is_compressed_LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

nfu.$(FC_MODEXT): nfu.$(OBJEXT)
nfu_compress.$(FC_MODEXT): nfu_compress.$(OBJEXT)

nfu_compress.$(OBJEXT): nfu.$(FC_MODEXT)
combine_ncc.$(OBJEXT): nfu.$(FC_MODEXT) nfu_compress.$(FC_MODEXT)
decompress_ncc.$(OBJEXT): nfu.$(FC_MODEXT) nfu_compress.$(FC_MODEXT)
scatter_ncc.$(OBJEXT): nfu.$(FC_MODEXT) nfu_compress.$(FC_MODEXT)

clean-local:
	$(RM) *.$(FC_MODEXT)

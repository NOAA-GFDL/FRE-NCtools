bin_PROGRAMS = cr_lake_files

AM_CFLAGS = $(NETCDF_CFLAGS) -I$(top_srcdir)/tools/libfrencutils
LDADD = $(NETCDF_LDFLAGS) -lnetcdff $(NETCDF_LIBS) -lhdf5_hl -lhdf5 -lz -ldl 
AM_FCFLAGS = -I$(top_srcdir)/tools/libfrencutils $(FC_DEFAULT_REAL_KIND8_FLAG) $(NETCDF_FCFLAGS)

noinst_LIBRARIES = libcrlake.a
libcrlake_a_SOURCES = affinity.c \
                      memuse.c \
                      gradient_c2l.c \
                      create_xgrid.c \
                      interp.c \
                      threadloc.c \
                      nsclock.c \
                      mosaic_util.c \
                      read_mosaic.c

cr_lake_files_SOURCES = platform.F90 mpp_parameter.F90 mpp_data.F90 mpp.F90 mpp_memutils.F90 mpp_efp.F90 mpp_pset.F90 mpp_domains.F90 mpp_io.F90 fms_io.F90 memutils.F90 constants.F90 fms.F90 horiz_interp_type.F90 horiz_interp_conserve.F90 horiz_interp_spherical.F90 horiz_interp_bicubic.F90 horiz_interp_bilinear.F90 horiz_interp.F90 cr_lake_files.f90
cr_lake_files_LDADD = libcrlake.a $(LDADD)

mpp_data.$(OBJEXT): mpp_data.F90 mpp_parameter.$(OBJEXT)
mpp.$(OBJEXT): mpp.F90 mpp_data.$(OBJEXT) mpp_parameter.$(OBJEXT)
mpp_memutils.$(OBJEXT): mpp_memutils.F90 mpp.$(OBJEXT)
mpp_efp.$(OBJEXT): mpp_efp.F90 mpp.$(OBJEXT)
mpp_pset.$(OBJEXT): mpp_pset.F90 mpp.$(OBJEXT)
mpp_domains.$(OBJEXT): mpp_domains.F90 mpp.$(OBJEXT) mpp_parameter.$(OBJEXT) mpp_data.$(OBJEXT) mpp_memutils.$(OBJEXT) mpp_pset.$(OBJEXT) mpp_efp.$(OBJEXT)
mpp_io.$(OBJEXT): mpp_io.F90 mpp.$(OBJEXT) mpp_parameter.$(OBJEXT) mpp_domains.$(OBJEXT)
fms_io.$(OBJEXT): fms_io.F90 mpp_io.$(OBJEXT) mpp_domains.$(OBJEXT) mpp.$(OBJEXT) mpp_parameter.$(OBJEXT) platform.$(OBJEXT)
memutils.$(OBJEXT): memutils.F90 mpp_io.$(OBJEXT) mpp.$(OBJEXT)
fms.$(OBJEXT): fms.F90 mpp.$(OBJEXT) mpp_domains.$(OBJEXT) mpp_io.$(OBJEXT) fms_io.$(OBJEXT) memutils.$(OBJEXT) constants.$(OBJEXT)
horiz_interp_type.$(OBJEXT): horiz_interp_type.F90 mpp.$(OBJEXT)
horiz_interp_conserve.$(OBJEXT): horiz_interp_conserve.F90 mpp.$(OBJEXT) fms.$(OBJEXT) fms_io.$(OBJEXT) constants.$(OBJEXT) horiz_interp_type.$(OBJEXT)
horiz_interp_spherical.$(OBJEXT): horiz_interp_spherical.F90 mpp.$(OBJEXT) fms.$(OBJEXT) constants.$(OBJEXT) horiz_interp_type.$(OBJEXT)
horiz_interp_bicubic.$(OBJEXT): horiz_interp_bicubic.F90 mpp.$(OBJEXT) fms.$(OBJEXT) constants.$(OBJEXT) horiz_interp_type.$(OBJEXT)
horiz_interp_bilinear.$(OBJEXT): horiz_interp_bilinear.F90 mpp.$(OBJEXT) fms.$(OBJEXT) constants.$(OBJEXT) horiz_interp_type.$(OBJEXT)
horiz_interp.$(OBJEXT): horiz_interp.F90 mpp.$(OBJEXT) fms.$(OBJEXT) constants.$(OBJEXT) horiz_interp_type.$(OBJEXT) horiz_interp_conserve.$(OBJEXT) horiz_interp_bilinear.$(OBJEXT) horiz_interp_bicubic.$(OBJEXT) horiz_interp_spherical.$(OBJEXT)
cr_lake_files.$(OBJEXT): cr_lake_files.f90 horiz_interp.$(OBJEXT)

CLEANFILES = *.mod


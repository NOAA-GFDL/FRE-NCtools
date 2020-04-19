bin_PROGRAMS = cr_lake_files 

AM_CFLAGS = -I$(top_srcdir)/tools/libfrencutils \
            $(NETCDF_CFLAGS)
LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS) -lnetcdff -lnetcdf -lhdf5_hl  -lhdf5 -lz -ldl

AM_FCFLAGS =  -I$(top_srcdir)/tools/libfrencutils \
              $(FC_DEFAULT_REAL_KIND8_FLAG) $(NETCDF_FCFLAGS)

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

cr_lake_files_SOURCES =  platform.F90 mpp_parameter.F90 mpp_data.F90 mpp.F90 mpp_memutils.F90 mpp_efp.F90 mpp_pset.F90 mpp_domains.F90 mpp_io.F90 fms_io.F90 memutils.F90 constants.F90 fms.F90 horiz_interp_type.F90 horiz_interp_conserve.F90 horiz_interp_spherical.F90 horiz_interp_bicubic.F90 horiz_interp_bilinear.F90 horiz_interp.F90 gradient.F90 mosaic.F90 mpp_utilities.F90 grid.F90 cr_lake_files.f90

cr_lake_files_LDADD = $(top_builddir)/tools/libfrencutils/libfrencutils.a $(LDADD) libcrlake.a   

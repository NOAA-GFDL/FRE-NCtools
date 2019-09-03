noinst_LIBRARIES = libfrencutils.a
if WITH_MPI
  noinst_LIBRARIES += libfrencutils_mpi.a
endif

AM_CFLAGS = $(NETCDF_CFLAGS)

libfrencutils_a_SOURCES = affinity.c \
                          constant.h \
                          create_xgrid.c \
                          create_xgrid.h \
                          gradient_c2l.c \
                          gradient_c2l.h \
                          interp.c \
                          interp.h \
                          mosaic_util.c \
                          mosaic_util.h \
                          mpp_domain.c \
                          mpp_domain.h \
                          mpp_io.c \
                          mpp_io.h \
                          mpp.c \
                          mpp.h \
                          read_mosaic.c \
                          read_mosaic.h \
                          tool_util.c \
                          tool_util.h

libfrencutils_mpi_a_SOURCES = $(libfrencutils_a_SOURCES)
libfrencutils_mpi_a_CFLAGS = -Duse_libMPI $(MPI_CFLAGS) $(AM_CFLAGS)
libfrencutils_mpi_a_CC = $(MPICC)

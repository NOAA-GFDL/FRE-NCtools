CLIBS_SITE  := -L/app/mpich2-1.2.1p1/lib -lmpich -lopa  -lpthread -lrt
CFLAGS_SITE := -I/app/mpich2-1.2.1p1/lib -DMPICH_IGNORE_CXX_SEEK

# GFDL uses the mpicc wrapper and Intel icc
MPICC    := icc
CC       := icc

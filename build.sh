#!/bin/bash
#Builds FRE-NCtools and runs regression tests
autoreconf -i configure.ac
./configure $MPI $QUAD_P
make -j check LOG_DRIVER_FLAGS=--comments

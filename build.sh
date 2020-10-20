#!/bin/bash
#Builds FRE-NCtools and runs regression tests
autoreconf -i configure.ac
./configure $MPI
make -j check LOG_DRIVER_FLAGS=--comments

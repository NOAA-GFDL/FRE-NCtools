#!/bin/bash
cd FRE-NCtools
autoreconf -i configure.ac
./configure $MPI
make -j check LOG_DRIVER_FLAGS=--comments

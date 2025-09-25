#!/bin/bash

set -ex

which nc-config
which nf-config

echo 'building FRE-NCtools conda package...'
echo "SRC_DIR / Build directory is:"
pwd
echo "Contents of SRC_DIR / Build directory are:"
ls

echo "PRE CONFIGURATION::"
echo ""
echo "PATH is:"
echo $PATH
echo "RPATH is:"
echo $RPATH
echo "LD_LIBRARY_PATH is:"
echo $LD_LIBRARY_PATH
echo "DYLD_LIBRARY_PATH is:"
echo $DYLD_LIBRARY_PATH
echo ""
echo ""
echo "CC is:"
echo $CC
echo "CF is:"
echo $CF
echo "CFLAGS is:"
echo $CFLAGS
echo "FCFLAGS is:"
echo $FCFLAGS
echo "LDFLAGS is:"
echo $LDFLAGS

#CC=
#FC=
#CFLAGS=
#FCFLAGS=
#LDFLAGS=#
#LD_LIBRARY_PATH=${PREFIX}/lib

## this is sufficient
autoreconf -iv
./configure --prefix=$PREFIX || cat config.log
#./configure --prefix=$PREFIX --with-mpi || cat config.log
#./configure --prefix=$PREFIX --enable-quad-precision || cat config.log
#./configure --prefix=$PREFIX --enable-quad-precision --with-mpi || cat config.log

echo "POST CONFIGURATION::"
echo ""
echo ""
echo "PATH is:"
echo $PATH
echo "LD_LIBRARY_PATH is:"
echo $LD_LIBRARY_PATH
echo ""

echo "compiling/building"
make

echo "installing into $PREFIX"
make install

### to test, build-dir option, ala README
#autoreconf -iv
#mkdir build && cd build
#../configure --prefix=$PREFIX || cat config.log
##../configure --prefix=$PREFIX --with-mpi || cat config.log
##../configure --prefix=$PREFIX --enable-quad-precision || cat config.log
##../configure --prefix=$PREFIX --enable-quad-precision --with-mpi || cat config.logecho "compiling/building"
#make
#echo "installing into $PREFIX"
#make install

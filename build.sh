#!/bin/bash

set -ex

echo ""
which nc-config
echo "we have nc-config"
echo ""

echo ""
which nf-config
echo "we have nf-config"
echo ""

#echo ""
#PATH=
#echo "PATH is:"
#echo $PATH
#echo ""
#RPATH=
#echo "RPATH is:"
#echo $RPATH
#echo ""
#LD_LIBRARY_PATH=${PREFIX}/lib
#echo "LD_LIBRARY_PATH is:"
#echo $LD_LIBRARY_PATH
#echo ""
#DYLD_LIBRARY_PATH=
#echo "DYLD_LIBRARY_PATH is:"
#echo $DYLD_LIBRARY_PATH
#echo ""

echo ""
#CC=mpicc
CC=`nc-config --cc`
echo "CC is:"
echo $CC
echo ""

echo ""
#FC=mpifc
FC=`nf-config --fc`
echo "FC is:"
echo $FC
echo ""

echo ""
CFLAGS=`nc-config --cflags`
echo "CFLAGS is:"
echo $CFLAGS
echo ""

echo ""
FCFLAGS=`nf-config --fflags`
echo "FCFLAGS is:"
echo $FCFLAGS
echo ""

echo ""
LDFLAGS=`nc-config --libs`
echo "LDFLAGS is:"
echo $LDFLAGS
echo ""

echo ""
echo "SRC_DIR / Build directory is:"
pwd
echo "Contents of SRC_DIR / Build directory are:"
ls
echo ""

echo ""
echo 'building FRE-NCtools conda package...'

## this is sufficient
autoreconf -iv
./configure --with-mpi || cat config.log
#./configure --prefix=$PREFIX --with-mpi || cat config.log
#./configure --prefix=$PREFIX --with-mpi || cat config.log
#./configure --prefix=$PREFIX --enable-quad-precision --with-mpi || cat config.log

echo ""
echo "compiling/building"
make

echo ""
#echo "installing into $PREFIX"
echo "installing no PREFIX"
make install

echo ""
echo "trying a make check in the build process, not advisable but i want info"
make check RUN_EXPENSIVE_TESTS=no \
	|| echo "make check failed- see test-suite.log, guarding against the failure to not clobber helpful output"

cp tests/test-suite.log /app/fre-nctools/tarball/test-suite.log || echo "copying test-suite log failed"


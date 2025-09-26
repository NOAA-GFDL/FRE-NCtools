#!/bin/bash

#set -e

## testing locally on ppan/ workstation
#module load gcc
#module load hdf5
#module load mpich
#module load hdf5
#module load netcdf-c
#module load nco
#module load netcdf-fortran/
#module list


echo ""
which nc-config
echo "we have nc-config"
echo ""

echo ""
which nf-config
echo "we have nf-config"
echo ""

echo ""
echo "BUILD_PREFIX/bin"
ls $BUILD_PREFIX/bin
echo ""
echo "BUILD_PREFIX/include"
ls $BUILD_PREFIX/include
echo ""


echo ""
echo "PREFIX/bin"
ls $PREFIX/bin
echo ""
echo "PREFIX/include"
ls $PREFIX/include
echo ""

#echo ""
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
#PREFIX=/home/inl/Working/fre-nctools/FRENCTOOLS
#PREFIX=/home/inl/FOO_BUILD/FRENCTOOLS
autoreconf -iv --include $BUILD_PREFIX/include
./configure --includedir $BUILD_PREFIX/include --prefix=$PREFIX --enable-quad-precision --with-mpi || cat config.log
#./configure --with-mpi || cat config.log
#./configure --prefix=$PREFIX --with-mpi || cat config.log
#./configure --prefix=$PREFIX --with-mpi || cat config.log

echo ""
echo "compiling/building"
make

echo ""
#echo "installing into $PREFIX"
echo "installing no PREFIX"
make install

cp -r lib/ $PREFIX
cp -r src/ $PREFIX
cp -r tests/ $PREFIX
#cp -r tools/ $PREFIX
cp -r man/ $PREFIX
cp -r m4/ $PREFIX
cp -r docs/ $PREFIX

#echo ""
#echo "trying a make check in the build process, not advisable but i want info"
#make check RUN_EXPENSIVE_TESTS=no \
#	|| echo "make check failed- see test-suite.log, guarding against the failure to not clobber helpful output"

#ls $SRC_DIR/tests/test-suite.log || echo "test-suite.log not found at $SRC_DIR/tests/test-suite.log"
#cp $SRC_DIR/tests/test-suite.log /app/fre-nctools-tarball || echo "copying test-suite log failed"



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
export CC=`nc-config --cc`
echo "CC is:"
echo $CC
echo ""

echo ""
#FC=mpifc
export FC=`nf-config --fc`
echo "FC is:"
echo $FC
echo ""

echo ""
export CFLAGS=`nc-config --cflags`
echo "CFLAGS is:"
echo $CFLAGS
echo ""

echo ""
export FCFLAGS=`nf-config --fflags`
echo "FCFLAGS is:"
echo $FCFLAGS
echo ""

echo ""
export LDFLAGS=`nc-config --libs`
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

cp aclocal.m4 $PREFIX          || echo "oops couldnt do it"
cp ar-lib $PREFIX               || echo "oops couldnt do it"
cp build.sh $PREFIX               || echo "oops couldnt do it"
cp CODE_OF_CONDUCT.md $PREFIX  || echo "oops couldnt do it"
cp compile $PREFIX               || echo "oops couldnt do it"
cp config.h $PREFIX               || echo "oops couldnt do it"
cp config.h.in $PREFIX           || echo "oops couldnt do it"
cp config.log $PREFIX           || echo "oops couldnt do it"
cp config.status $PREFIX       || echo "oops couldnt do it"
cp configure $PREFIX           || echo "oops couldnt do it"
cp configure.ac $PREFIX           || echo "oops couldnt do it"
cp CONTRIBUTING.md $PREFIX       || echo "oops couldnt do it"
cp depcomp $PREFIX               || echo "oops couldnt do it"
cp environment.yml $PREFIX       || echo "oops couldnt do it"
cp install-sh $PREFIX           || echo "oops couldnt do it"
cp LICENSE.md $PREFIX           || echo "oops couldnt do it"
cp Makefile $PREFIX               || echo "oops couldnt do it"
cp Makefile.am $PREFIX           || echo "oops couldnt do it"
cp Makefile.in $PREFIX           || echo "oops couldnt do it"
cp meta.yaml $PREFIX           || echo "oops couldnt do it"
cp missing $PREFIX               || echo "oops couldnt do it"
cp README.md $PREFIX           || echo "oops couldnt do it"
cp stamp-h1 $PREFIX               || echo "oops couldnt do it"
cp tap-driver.sh $PREFIX       || echo "oops couldnt do it"
cp test-driver $PREFIX           || echo "oops couldnt do it"
cp autom4te.cache/** $PREFIX   || echo "oops couldnt do it"
cp -r docs/* $PREFIX           || echo "oops couldnt do it"
cp -r lib/* $PREFIX               || echo "oops couldnt do it"
cp -r m4/*  $PREFIX               || echo "oops couldnt do it"
cp -r man/* $PREFIX               || echo "oops couldnt do it"
cp -r site-configs/* $PREFIX   || echo "oops couldnt do it"
cp -r src/* $PREFIX               || echo "oops couldnt do it"
cp -r tests/* $PREFIX           || echo "oops couldnt do it"
#cp -r tools/* $PREFIX          || echo "oops couldnt do it"

#echo ""
#echo "trying a make check in the build process, not advisable but i want info"
#make check RUN_EXPENSIVE_TESTS=no \
#    || echo "make check failed- see test-suite.log, guarding against the failure to not clobber helpful output"

#ls $SRC_DIR/tests/test-suite.log || echo "test-suite.log not found at $SRC_DIR/tests/test-suite.log"
#cp $SRC_DIR/tests/test-suite.log /app/fre-nctools-tarball || echo "copying test-suite log failed"

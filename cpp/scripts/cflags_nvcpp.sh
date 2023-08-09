#!/bin/bash

if [ "x$1" == "xgpu" ]
then
    olev=2
    cpphware=gpu
elif [ "x$1" == "xrelease" ]
then
    olev=2
    cpphware=multicore
else
    olev=0
    cpphware=multicore
fi

export CC=nvc
export CXX=nvc++
export CXXFLAGS="-std=c++20 -stdpar=$cpphware -O$olev -g"
export CFLAGS="-O$olev -g"

echo "CC:" $CC
echo "CFLAGS:" $CFLAGS
echo "CXX:" $CXX
echo "CXXFLAGS:" $CXXFLAGS

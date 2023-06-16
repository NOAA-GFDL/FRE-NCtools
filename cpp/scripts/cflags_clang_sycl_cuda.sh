#!/bin/bash

#settingd for compiling with the Intel/CodePlay Clang based SYCL compiler.
#In the Edison build, compiler was at
#/home/mzuniga/llvm/build/install/bin/clang++
#Todo: netcdf libs, etc?

if [ "x$1" == "xrelease" ]
then
    olev=2
else
    olev=0
fi

export CC=clang
export CXX=clang++
export CXXFLAGS="-std=c++20 -fsycl -fsycl-targets=nvptx64-nvidia-cuda"
export CFLAGS="-O$olev -g"

echo "CC:" $CC
echo "CFLAGS:" $CFLAGS
echo "CXX:" $CXX
echo "CXXFLAGS:" $CXXFLAGS

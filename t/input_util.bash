#!/bin/sh
# Ryan Mulhall 2021
# Utility functions for FRE-NCtools bats test script
# SETUP_FNCT should be set to a function that creates the tests input files
# such as generate_from_ncl, otherwise no input will be created

testDir="$( basename $BATS_TEST_FILENAME .sh )"
testDir="`expr substr $testDir 1 6`"

setup(){
  BASE_TEST_DIR=$( pwd )
  mkdir $testDir
  cd $testDir
  # run the tests setup function, if set
  [[ ! -z $SETUP_FNCT ]] && $SETUP_FNCT
}

teardown(){
  cd $BASE_TEST_DIR
  rm -rf $testDir
}

generate_from_ncl(){
  filename=$1
#mppnccombine.
  #Generate the netcdf files from the .ncl
  for f in $top_srcdir/t/Test??-input/*.ncl.????
  do
    [[ -z $filename ]] && filename=$(basename $f .ncl)
    echo $filename
    ncgen -o $filename.nc.${f##*.} $f
  done
}

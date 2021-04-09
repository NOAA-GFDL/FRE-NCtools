#!/bin/bash
# Ryan Mulhall 2021
# Utility functions for FRE-NCtools bats test scripts
# setup and teardown run implicitly at beg and end of tests
# SETUP_FNCT must set to a function(or cmd) that creates the tests input files
# such as generate_from_ncl, otherwise no input will be created

testDir="$( basename $BATS_TEST_FILENAME .sh )"
testDir="`expr substr $testDir 1 6`"

setup(){
  echo $testDir
  BASE_TEST_DIR=$( pwd )
  mkdir $testDir
  cd $testDir
  # run the tests setup function, if set
  [[ ! -z $SETUP_FNCT ]] && $SETUP_FNCT || echo "SETUP_FNCT not set, skipping input generation"
}

teardown(){
  cd $BASE_TEST_DIR
  rm -rf $testDir
}

generate_all_from_ncl_num(){
  test_no="`expr substr $(basename $BATS_TEST_FILENAME .sh) 5 2`"
  echo "Test no: $test_no"
  filename=$1
  #Generate the netcdf files from the .ncl
  for f in $top_srcdir/t/Test${test_no}-input/*.ncl.????
  do
    [[ -z $filename ]] && filename=$(basename $f .ncl)
    ncgen -o $filename.nc.${f##*.} $f
  done
}

generate_all_from_ncl(){
  if [[ -z $1 ]]; then 
    test_no="`expr substr $(basename $BATS_TEST_FILENAME .sh) 5 2`"
  else
    test_no="$1"
  fi
  #Generate the netcdf files from the .ncl
  for f in $top_srcdir/t/Test${test_no}-input/*.ncl
  do
    filename=$(basename $f .ncl)
    ncgen -o $filename.nc $f
  done
}

#!/usr/bin/env bats

#***********************************************************************
#                   GNU Lesser General Public License
#
# This file is part of the GFDL FRE NetCDF tools package (FRE-NCTools).
#
# FRE-NCTools is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# FRE-NCTools is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with FRE-NCTools.  If not, see
# <http://www.gnu.org/licenses/>.
#***********************************************************************
# Ryan Mulhall 2021
# Utility functions for FRE-NCtools bats test scripts
# setup and teardown run implicitly at beg and end of tests

# SETUP_FNCT is set to a function name or command in test files to specify how to create input

# SKIP_TESTS is to be set to a space separated list of one or two digit test numbers to be skipped

testDir="$( basename $BATS_TEST_FILENAME .sh )"
testDir="`expr substr $testDir 1 6`"

setup(){
  BASE_TEST_DIR=$( pwd )
  TEST_NUM="`expr substr $(basename $BATS_TEST_FILENAME .sh) 5 2`"

  # check if test is to be skipped first
  [[ ! -z $SKIP_TESTS ]] && skip_test $SKIP_TESTS

  if [ ! -d $testDir ]; then
    mkdir $testDir
  fi
  cd $testDir

  # run the tests setup function, if set
  [[ ! -z $SETUP_FNCT ]] && $SETUP_FNCT || echo "SETUP_FNCT not set, no input will be created"
}

teardown(){
  cd $BASE_TEST_DIR
  rm -rf $testDir
}

skip_test(){
  if [[ "$1" == "$TEST_NUM" || "0$1" == $TEST_NUM ]]; then
    skip "Set to skip in SKIP_TESTS"
  else
    if [[ $# -gt 1 ]]; then
      shift
      skip_test $@
    fi
  fi
}

#generates numbered nc files with a given name
generate_all_from_ncl_num(){
  filename=$1
  for f in $top_srcdir/t/Test${TEST_NUM}-input/*.ncl.????
  do
    [[ -z $filename ]] && filename=$(basename $f .ncl)
    ncgen -o $filename.nc.${f##*.} $f
  done
}
#generates nc files from any ncl files in input directory
#takes a two digit test number, defaults to current test
generate_all_from_ncl(){
  if [[ -z $1 ]]; then 
    test_no=$TEST_NUM
  else
    test_no="$1"
  fi
  for f in $top_srcdir/t/Test${test_no}-input/*.ncl
  do
    filename=$(basename $f .ncl)
    ncgen -o $filename.nc $f
  done
}

#generates from a list of names, TODO might delete
generate_from_ncl(){
  if [[ $# == 0 || $1 != *.ncl ]]; then
    echo "Error: no/invalid filenames provided"
    exit
  else
    filename=$(basename $1 .ncl)
    ncgen -o $filename.nc $top_srcdir/t/$1
    [[ $# -eq 1 ]] && return 0
    shift
    generate_from_ncl $@
  fi
}

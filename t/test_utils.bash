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

# SKIP_TESTS can be set to a space separated list of one or two digit test numbers
# to be skipped during make check, e.g. SKIP_TESTS="1 2 10", etc.

testDir="$( basename $BATS_TEST_FILENAME .sh )"
testDir=${testDir%-*}

setup(){
  BASE_TEST_DIR=$( pwd )
  TEST_NUM="${testDir:4:2}"

  # check if test is to be skipped first
  [[ ! -z $SKIP_TESTS ]] && skip_test $SKIP_TESTS

  if [ ! -d $testDir ]; then
    mkdir $testDir
  fi
  cd $testDir
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
# $1 = filename to generate, $2 = input dir to generate from
generate_all_from_ncl_num(){
  filename=$1
  for f in $top_srcdir/t/$2/*.ncl.????
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

# run the command (the argument list; 1st in list is the command) using the bats "run command",
# theck the status is 0, and echo the output if not.
function run_and_check()
{
    local cmd="$*" #Expands the list into a single string; spearating parms with space.
    run command $cmd
    if [ $status -ne 0 ];
    then
	    echo failed:  $cmd
	    echo $output && exit 1
    fi
}

# make a csv (comma separated value) list of
# all the file names in dir argument 1 (i.e. dir arg) that
# have the pattern of argument 2 (i.e. pattern arg)
function get_csv_filename_list()
{
    local dir=$1
    local pattern=$2
    local filelist=""
    for file in $dir"/"$pattern; do
	if [[ $filelist == "" ]]
	then
	    filelist="$file" #if this is the first file
	else
	    filelist=$filelist",""$file"
	fi
    done
    echo $filelist  #I.e. this return must be capture by a $(...)
}

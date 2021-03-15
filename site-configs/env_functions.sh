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
# License along with FRE-NCTools (LICENSE.md).  If not, see
# <http://www.gnu.org/licenses/>.
#***********************************************************************
#
# Copyright (c) 2021 - Seth Underwood (@underwoo)

# Functions to help configure the required environment
# These functions assume the system us using the Environment modules
# system (https://github.com/cea-hpc/modules) version 4 or later.

# echo_stderr `echo`s a string to stderr
echo_stderr ()
{
  echo "$*" 1>&2
}

# echo_error/echo_info will write a message to stderr
echo_error ()
{
  echo_stderr "ERROR: $*"
}

echo_info ()
{
  echo_stderr "INFO: $*"
}

# get_sh returns the parent shell (the shell used to call the script).
get_sh ()
{
  # Get the parent process/shell
  ppid=$( ps -o ppid= -p $$ )
  echo $( ps -o command= -p $ppid | tr -d '-' | xargs basename )
}

# module is a standin for the `module` command, and ensures the module command
# is echo'd to the screen with a trailing ';'.  The ';' is needed for `eval` to
# correctly interpret the commands.
#
# This `module` standin also contains the newer module subcommands `prepend-path`
# `remove_path` and `append_path` that are found in the 4.X version.
module ()
{
  case "$1" in
    prepend-path)
      shift
      prepend_path $*
      ;;
    append-path)
      shift
      append_path $*
      ;;
    remove-path)
    shift
      remove_path $*
      ;;
    load|remove|rm|swap|switch)
      echo "module $*;"
      ;;
    *)
      echo_error "unknown subcommand"
      exit 1
      ;;
  esac
}

# prepend_path prepends a string to an environment variable using `:` as a
# delimiter
prepend_path ()
{
  if [ "$( eval echo \${$1:-unset}; )" = 'unset' ]
  then
    setenv $*
    eval $1=1
  else
    case $( get_sh ) in
      bash|ash|zsh|sh|ksh)
        echo "$1=$2:\${$1}; export $1;"
        ;;
      tcsh|csh)
        echo "setenv $1 $2:\${$1};"
        ;;
      *)
        echo_error "unknown shell type" 1>&2
        exit 1
        ;;
    esac
  fi
}

# append_path appends a string to an environment variable using `:` as a
# delimiter
append_path ()
{
  if [ ${1:-unset} = 'unset' ]
  then
    setenv $*
    eval $1=1
  else
    case $( get_sh ) in
      bash|ash|zsh|sh|ksh)
        echo "$1=\${$1}:$2; export $1;"
        ;;
      tcsh|csh)
        echo "setenv $1 \${$1}\:$2;"
        ;;
      *)
        echo_error "unknown shell type" 1>&2
        exit 1
        ;;
    esac
  fi
}

# remove_path removes a string from an environment variable that uses `:` as
# a delimiter.
remove_path ()
{
  if [ ${1:+set} = 'set' ]
  then
    case $( get_sh ) in
      bash|ash|zsh|sh|ksh)
      echo "$1=\`echo :\${$1} | sed -e 's@:$2@@' -e 's/^://'\`; export $1;"
      ;;
    tcsh|csh)
      echo "setenv $1 \`echo :\${$1} | sed -e 's@:$2@@' -e 's/^://'\`;"
      ;;
    *)
      echo_error "unknown shell type" 1>&2
      exit 1
      ;;
    esac
  fi
}

# setenv is used similar to how the modulefiles will use setenv to set environment
# variable.  setenv will output the correct commands for the shell type (csh or sh).
setenv ()
{
  case $( get_sh ) in
    bash|ash|zsh|sh|ksh)
      echo "$1=$2; export $1;"
      ;;
    tcsh|csh)
      echo "setenv $1 $2;"
      ;;
    *)
      echo_error "unknown shell type" 1>&2
      exit 1
      ;;
  esac
}

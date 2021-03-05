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

# modpath_prepend uses the `module remote-path` and `module prepend-path`
# to ensure the module directory passed in is the first listed in MODULEPATH
modpath_prepend ()
{
  module remove-path MODULEPATH $1
  module prepend-path MODULEPATH $1
}

# module is a standin for the `module` command, and ensures the module command
# is echo'd to the screen with a trailing ';'.  The ';' is needed for `eval` to
# correctly interpret the commands.
module ()
{
  echo "module $*;"
}

# setenv is used similar to how the modulefiles will use setenv to set environment
# variable.  setenv will output the correct commands for the shell type (csh or sh).
setenv ()
{
  # Get the parent process/shell
  ppid=$( ps -o ppid= -p $$ )
  pcommand=$( ps -o command= -p $ppid | tr -d '-' | xargs basename )

  case $pcommand in
    bash|ash|zsh|sh|ksh)
      shtype=sh
      setenv=
      eq='='
      ;;
    tcsh|csh)
      shtype=csh
      setenv="setenv"
      eq=' '
      ;;
    *)
      echo "Unknown shell type" 1>&2
      exit 1
      ;;
  esac

  echo -n $setenv "$( echo $1 | tr '=' "$eq" );"
  if [ "$shtype" = "sh" ]
  then
    echo -n " export $( echo $1 | cut -d '=' -f 1 );"
  fi
  echo 
}

#!/bin/tcsh -f
#
# $Id: varlist.csh,v 1.1.2.1.2.2 2010/06/01 23:01:57 afy Exp $
# ------------------------------------------------------------------------------
# FMS/FRE Project: Script to Print Info About netCDF Variables
# ------------------------------------------------------------------------------
# afy    Ver   1.00  Copied from ~fms/local/ia64/netcdf4.fix        June 10
# afy    Ver   2.00  Don't source the 'init.csh' script             June 10
# afy    Ver   2.01  Use 'which' to locate the 'list_ncvars.csh'    June 10
# ------------------------------------------------------------------------------
# Copyright (C) NOAA Geophysical Fluid Dynamics Laboratory, 2000-2010
# Designed and written by V. Balaji, Amy Langenhorst and Aleksey Yakovlev
#

set file = $1

##################################################################

alias list_ncvars.csh `which list_ncvars.csh`

##################################################################

if (-e .tmpdump) then
   echo ERROR: File .tmpdump exists.
   exit 1
endif

ncdump -h $file > .tmpdump
set vars = `list_ncvars.csh -t0123 $file`

foreach var ($vars)
#   --- get the long name and units ---
    set string = `grep "	${var}:long_name" .dump`
    set name = `echo $string[3-] | sed -e 's/"/ /g' | sed -e 's/;/ /g'`
    set string = `grep "	${var}:units" .dump`
    set units = `echo $string[3-] | sed -e 's/"/ /g' | sed -e 's/;/ /g'`
    echo $var.nc $name "("$units")" | awk '{printf "%20s = ", $1}{for (i = 2; i < NF; ++i) printf "%s ", $i}{printf "%s\n", $NF}'
end

if (-e .dump) rm -f .tmpdump

#!/bin/csh -f
#
# $Id: split_ncvars.csh,v 1.1.2.1.2.2.2.1.2.2.4.4 2013/07/23 13:09:56 keo Exp $
# ------------------------------------------------------------------------------
# FMS/FRE Project: Script to Split netCDF Files
# ------------------------------------------------------------------------------
# afy    Ver   1.00  Copied from ~fms/local/ia64/netcdf4.fix        June 10
# afy    Ver   2.00  Don't source the 'init.csh' script             June 10
# afy    Ver   2.01  Add aliases ncrcat/ncks (from the 'init.csh')  June 10
# afy    Ver   2.02  Use 'which' to locate the 'list_ncvars.csh'    June 10
# afy    Ver   2.03  Use 'which' to locate the 'varlist.csh'        June 10
# ------------------------------------------------------------------------------
# Copyright (C) NOAA Geophysical Fluid Dynamics Laboratory, 2000-2012
# Designed and written by V. Balaji, Amy Langenhorst and Aleksey Yakovlev
#

set odir = .
set ifiles
set splitvarsstatus = 0

#  ----- parse input argument list ------

set argv = (`getopt sdlf:i:o:v: $*`)

while ("$argv[1]" != "--")
    switch ($argv[1])
        case -d:
            set debug; breaksw
        case -l:
            set log; breaksw
        case -f:
            set onefile = $argv[2]; shift argv; breaksw
        case -i:
            set idir = $argv[2]; shift argv; breaksw
        case -o:
            set odir = $argv[2]; shift argv; breaksw
        case -v:
            set vars = $argv[2]; shift argv; breaksw
        case -s:
            set static; breaksw
    endsw
    shift argv
end
shift argv
if ($?debug) set echo
set ifiles = ( $argv )
set first

##################################################################
#  ----- help message -----

if ($#ifiles == 0) then
set name = `basename $0`
cat << EOF

Split NetCDF Variables

Usage:  $name [-d] [-l] [-i idir] [-o odir] [-f file] [-v vars]  files.....

        -l       = include readme log in output directory
        -d       = debug option: turns on command echo
        -i idir  = input (archive) directory path
        -o odir  = output directory path
        -f file  = one file output option
                   file is the name of the output file
                   this option must be used with -v option
        -v vars  = comma separated list of variables
        files... = list netcdf files in directory idir
                   the files must be in chronological order

EOF
exit 1
endif

##################################################################

alias ncrcat `which ncrcat` -t 2 --header_pad 16384
alias ncks `which ncks` --header_pad 16384
alias list_ncvars.csh `which list_ncvars.csh`
alias varlist.csh `which varlist.csh`

##################################################################

# user supplied variable list
  if ($?vars) then
     set varlist = `echo $vars | sed -e "s/,/ /g"`
  else
     # vars list must be set to one file option
     if ($?onefile) then
        unset onefile
     endif
  endif

# -f option disallow >1 input files
if ($?onefile && $#ifiles > 1) then
    echo "ERROR: only one input file allowed when using -f option"
    exit 1
endif

##################################################################

#  need to make output directory

 if (! -e $odir) mkdir -p $odir

#  process each input file separately

 foreach file ($ifiles)
     if ($?idir) then
        dmcopy $idir/$file $file
     endif
     ncdump -h $file > .dump

   # generate variable list (if not user supplied)
     if (! $?varlist) then
        if (! $?static) then
           set varlist = `list_ncvars.csh -t0123 $file`
        else
           set varlist = `list_ncvars.csh -s0123 $file`
        endif
     endif

   # get time axis name
     set string = `grep UNLIMITED .dump`
     if ($#string > 0) then
         set timename = `echo $string[1]`
     else
         set timename = ""
     endif
     if ($?static) set timename = ""
#Balaji: add all variables that appear in a coordinate attribute
     set coords = ( `ncdump -h $file |grep :coordinates |cut -f2 -d\" |sed "s/ /\n/g" |sort -u` \
                    `ncdump -h $file |grep "float geol[oa][nt]" | awk '{print $2}' |cut -f1 -d\(` )
     set coords = `echo $coords |sed "s/ /\n/g" |sort -u`
#     exec echo coords = $coords
#### loop through variables ####

     set vlist = ( $coords )
     foreach var ($varlist)
        # create comma-separated list of variables to extract
         set vlist = ($vlist $var)
         if ($?echo) unset echo
     #--- first time: need to extract additional static fields ---
         if (! -e $odir/$var.nc) then
            # get axis names then edge names (may want to get bounds for CF compliance)
            set dimlist = `grep " $var(" .dump | sed -e "s/(/ /" -e "s/,/ /g" -e "s/)/ /" -e "s/;/ /"`
            if ($#dimlist > 2) then
                foreach dim ($dimlist[3-])
                    if ($dim == $timename) continue # skip time axis (non-static)
                    set string = `grep "	${dim}:edges" .dump`
                    if ($#string > 2) set vlist = ($vlist `echo $string[3] | sed -e 's/"/ /g'`)
                    set string = `grep "	${dim}:bounds" .dump`
                    if ($#string > 2) set vlist = ($vlist `echo $string[3] | sed -e 's/"/ /g'`)
                end
            endif
         endif
     #--- every time: need to extract time averaging strings ---
        # fms-style time avg variable
         if ($timename != "") then
            set string = `grep "	${var}:time_avg_info" .dump`
            if ($#string > 2) set vlist = ($vlist `echo $string[3] | sed -e 's/"/ /g'`)
           # CF compliant time avg variables (bounds or climatology)
            foreach tbound ( bounds climatology )
                set string = `grep "	${timename}:$tbound" .dump`
                if ($#string > 2) then
                    set name = `echo $string[3] | sed -e 's/"/ /g'`
                    if ($name == $var) goto END_OF_VAR_LOOP
                    set vlist = ($vlist $name)
                endif
            end
         endif

         if ($?debug) set echo
         if ($?onefile) continue

     #--- extract fields ---
         set vlist = `echo $vlist | sed -e  "s/ /,/g"`
         set appendopt = ''
         if ( -f ./.var.nc ) set appendopt = "-A"
         echo ncks -h $appendopt -v $vlist $file $odir/$var.nc
         ncks -h $appendopt -v $vlist $file .var.nc
         set thisstatus = $status
         if ( $thisstatus != 0 ) then
            echo ERROR: ncks returned status $thisstatus
            set splitvarsstatus = $thisstatus
            sleep 30
            ncks -h $appendopt -v $vlist $file .var.nc
            set thisstatus = $status
            if ( $thisstatus != 0 ) then
               echo ERROR ON RETRY: ncks returned status $thisstatus
               set splitvarsstatus = $thisstatus
            else
               echo RETRY SUCCESSFUL.
               set splitvarsstatus = $thisstatus
            endif
         endif

     #--- write to the output file ---
         if (-e $odir/$var.nc) then
             if ($timename == "") then
                echo "ERROR: try to append to unlimited dimension when no unlimited dimension found"
                exit 1
             endif
            #--- append to existing file ---
            #dmget $odir/$var.nc
             mv $odir/$var.nc $odir/.tmp.nc
             ncrcat -h $odir/.tmp.nc .var.nc $odir/$var.nc
             set thisstatus = $status
             if ( $thisstatus != 0 ) then
               echo ERROR: ncrcat returned status $thisstatus
               set splitvarsstatus = $thisstatus
               sleep 30
               ncrcat -h $odir/.tmp.nc .var.nc $odir/$var.nc
               set thisstatus = $status
               if ( $thisstatus != 0 ) then
                  echo ERROR ON RETRY: ncrcat returned status $thisstatus
                  set splitvarsstatus = $thisstatus
               else
                  echo RETRY SUCCESSFUL.
                  set splitvarsstatus = $thisstatus
               endif
             endif
             rm -f $odir/.tmp.nc .var.nc
         else
            #--- create new file ---
            #--- modify filename attribute ---
            #--- move to output directory ---
             ncatted -h -O -a filename,global,m,c,"$var.nc" .var.nc
             mv .var.nc $odir/$var.nc
         endif

         END_OF_VAR_LOOP:
         set vlist = ( $coords )
     end ### end of variable loop ###

    #--- extract fields for single file option ---
     if ($?onefile) then
         set vlist = `echo $vlist | sed -e  "s/ /,/g"`
         echo ncks -h -O -v $vlist $file $onefile
         ncks -h -O -v $vlist $file $onefile
         set thisstatus = $status
         if ( $thisstatus != 0 ) then
            echo ERROR: ncks returned status $thisstatus
            set splitvarsstatus = $thisstatus
            sleep 30
            ncks -h -O -v $vlist $file $onefile
            set thisstatus = $status
            if ( $thisstatus != 0 ) then
               echo ERROR ON RETRY: ncks returned status $thisstatus
               set splitvarsstatus = $thisstatus
            else
               echo RETRY SUCCESSFUL.
               set splitvarsstatus = $thisstatus
            endif
         endif
         ncatted -h -O -a filename,global,m,c,"$onefile:t" $onefile
     endif

    #--- create readme ---
     if ($?first && $?log) then
        varlist.csh $file > $odir/README
        unset first
     endif
    # remove file if copied from archive
     if ($?idir) rm -f $file

    # clean up
     if (-e .dump)    rm -f .dump
     if (-e ncks.out) rm -f ncks.out

 end ### end of file loop ###

exit $splitvarsstatus

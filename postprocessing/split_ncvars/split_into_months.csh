#!/bin/csh -f
# $Id: split_into_months.csh,v 1.3 2014/09/16 18:58:59 fms Exp $
# split FMS history file into months
# expected input filename pattern is yyyymmdd.*.nc
# returns files of the form yyyymm.*.nc
# also, month 12 of a year is synonymous with month 00 of next year
# beyond month 12 goes to next year
#  (to enable seasonal average for DJF)
#!!!First argument should be start year
#!!!(FMS convention is to treat Dec of first year as Dec of year before first year for DJF purposes)
set days_year = ( 31 28 31 30 31 30 31 31 30 31 30 31 )
set days_leap = ( 31 29 31 30 31 30 31 31 30 31 30 31 )
#set start = $1; shift

foreach in ( $* )
   set file = `basename $in`
   echo Processing $in ...

   set year = `echo $file |cut -c1-4`
   set next = `expr $year + 1`
   set next = `printf "%4.4d" $next`
   set rest = `echo $file |cut -c9-`
   set startmonth = `echo $file |cut -c5-6`
#find time axis name
   set tline = `ncdump -h $in |grep -F UNLIMITED`
#echo tline = $tline
   set taxis = `echo $tline |awk '{print $1}'`
   set ntime = `echo $tline |cut -f2 -d\( |awk '{print $1}'`
   set calendar = `ncdump -h $in |grep -F ${taxis}:calendar_type |cut -f2 -d\"`
   test $calendar = 'JULIAN' || test $calendar = 'NOLEAP' || exec echo ERROR: $0 only knows about JULIAN and NOLEAP calendars.
   @ m=$startmonth
   @ n=`expr $m - 12`           # month in next year; could be generalized for multiple years
   @ i=1
#how many time records in a month?
#if month file, 1
#if day file use days_year/days_leap
#if 8xdaily, multiply that by 8
   echo $file | grep -F -q '_month'   && set type = 'month'
   echo $file | grep -F -q '_daily'   && set type = 'daily'
   echo $file | grep -F -q '_8xdaily' && set type = '8xdaily'
   echo $file | grep -F -q '_4xdaily' && set type = '4xdaily'
   echo $file | grep -F -q '_tracer'  && set type = 'tracer'
   echo $file | grep -F -q '_scalar'  && set type = 'scalar' 

   set created     
   while( $i <= $ntime )
       @ k=`expr \( $m - 1 \) % 12 + 1`
       switch ( $type )
       case 'month':
         @ j = $i
         breaksw
       case 'daily':
         @ j = `expr $i + $days_year[$k] - 1`
         test $calendar = 'JULIAN' && test `expr $year % 4` -eq 0 && @ j = `expr $i + $days_leap[$k] - 1`
         breaksw
       case '8xdaily':
         @ j = `expr $i + $days_year[$k] \* 8 - 1`
         test $calendar = 'JULIAN' && test `expr $year % 4` -eq 0 && @ j = `expr $i + $days_leap[$k] \* 8 - 1`
         breaksw
       case '4xdaily':
         @ j = `expr $i + $days_year[$k] \* 4 - 1`
         test $calendar = 'JULIAN' && test `expr $year % 4` -eq 0 && @ j = `expr $i + $days_leap[$k] \* 4 - 1`
         breaksw
       case 'tracer':
         @ j = $i
	 breaksw
       case 'scalar':
         @ j = $i
	 breaksw
       default:
         exec echo ERROR filename $file does not resolve into monthly, daily or 8xdaily.
       endsw
       set mm = `printf "%2.2d" $m`
       set nn = `printf "%2.2d" $n`
	   echo "ncks -O -a -h -F -d $taxis,$i,$j $in $year$mm$rest"
       ncks -O -a -h -F -d $taxis,$i,$j $in $year$mm$rest || exec echo ERROR Unable to create $year$mm$rest.
       set created = "$created $year$mm$rest"
       test $m =   12 && ln -f $year$mm$rest $next$nn$rest && set created = "$created $next$nn$rest"
#       test $m =   12 && test $year = $start && ln $year$mm$rest $year$nn$rest && set created = "$created $year$nn$rest"
       test $m -gt 12 && mv -f $year$mm$rest $next$nn$rest && set created = "$created $next$nn$rest"
       @ i = `expr $j + 1`
       @ m++
       @ n++
   end
#special treatment for atmos... dubious practice:-(
#    set csave = $created
#    foreach foo ( $csave )
#      echo $foo | grep -q '^......\.atmos' && ln $foo $foo:r.modellevels.nc
#      set created = "$created $foo:r.modellevels.nc"
#    end
   echo `basename $0`: CREATED $created ...
#    test $i -ne `expr $ntime + 1` && exec echo \
#     Number of time records does not match expected! \
#     See if files created: $created ... are each one month long.
end

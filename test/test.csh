#!/bin/csh -f
#
# Script to copy and regrid cubed-sphere history data and test split_ncvars.pl
#
# USAGE
#   module load fre
#   test.csh subDir expDir yr1 yr2 yr3 ...

if ($#argv < 3) then
  echo "USAGE: $0:t subWorkDir experDir yr1 [yr2 yr3 ...]"
  exit 1
endif

set subDir = $1
set expDir = $2
set years = ( $argv[3-] )
if (! -e $expDir) then
  echo "ERROR: experiment directory does not exist"
  exit 1
endif

set FREGRID = `which fregrid` # /home/z1l/bin/20160428/fregrid # bronx-12
set SPLIT_NCVARS = /home/bw/sources/FMS/split_ncvars/split_ncvars.pl
set LIST_NCVARS = `which list_ncvars.csh`
if ($?) then
  echo "ERROR: fre not loaded"
  exit 1
endif

if (! $?TMPDIR) then
  echo "ERROR: TMPDIR not defined, probably not on supported platform"
  exit 1
endif

###############
#  find data
###############

set dmgetList = ()
foreach year ($years)
  if (-e $expDir/history/${year}0101.nc.tar) then
    set dmgetList = ($dmgetList $expDir/history/${year}0101.nc.tar)
  else
    echo "ERROR: No file: $expDir/history/${year}0101.tar"
    exit 1
  endif
end

set dmgetList2 = ()
if (-e $expDir/history_refineDiag) then
  foreach year ($years)
    if (-e $expDir/history_refineDiag/${year}0101.nc.tar) then
      set dmgetList2 = ($dmgetList2 $expDir/history_refineDiag/${year}0101.nc.tar)
    else
      echo "ERROR: No file: $expDir/history/${year}0101.tar"
      exit 1
    endif
  end
endif

###############
#  copy data
###############

if (! -e $TMPDIR/$subDir) mkdir $TMPDIR/$subDir
cd $TMPDIR/$subDir

if (! -e $years[1]0101.atmos_month.tile1.nc) then

@ nfiles = $#dmgetList + $#dmgetList2
echo "dmget $nfiles files"
dmget $dmgetList $dmgetList2

echo "copying $#dmgetList files"
gcp $dmgetList .
foreach year ($years)
  tar -xvf ${year}0101.nc.tar \*.atmos_\* \*.grid_spec.\*
 #tar -xvf ${year}0101.nc.tar \*.atmos_month.\* \*.atmos_level.\* \*.atmos_daily.\* \*.grid_spec.\*
  rm -f ${year}0101.nc.tar
end

if ($#dmgetList2 > 0) then
  echo "copying $#dmgetList2 refined files"
  gcp $dmgetList2 .
  foreach year ($years)
    tar -xvf ${year}0101.nc.tar
    rm -f ${year}0101.nc.tar
  end
endif

endif


############
#  regrid
############

# make the atmos exchange grid from scratch
set ncube = `ncdump -h $years[1]0101.grid_spec.tile1.nc | grep "grid_xt = " | awk '{print $3}'`
if (! -e grid/C${ncube}_mosaic.nc) then
  mkdir grid
  @ ngrid = 2 * $ncube
  make_hgrid --grid_type gnomonic_ed --nlon $ngrid --nlat $ngrid --grid_name grid/C${ncube}_grid
  set tilefiles = `echo C${ncube}_grid.tile{1,2,3,4,5,6}.nc | sed -e "s/ /,/g"`
  make_solo_mosaic --num_tiles 6 --dir grid --mosaic grid/C${ncube}_mosaic --tile_file $tilefiles
endif

set options = "--standard_dimension --nlon 288 --nlat 180 --input_mosaic grid/C${ncube}_mosaic.nc --interp_method conserve_order2 --remap_file .fregrid_remap_file.nc"
set filetypes = `/bin/ls *.*.tile1.nc`;
foreach file1 ($filetypes) 
  set base = $file1:r:r
  set comp = $base:e
  if ($comp != "grid_spec") then
    if (! -e $base.nc) then
      set vlist = `$LIST_NCVARS -st01234 $file1`
      set scalar_field = `echo $vlist | sed -e "s/ /,/g"`
      echo "$FREGRID $options --input_file $base --scalar_field $scalar_field"
      $FREGRID $options --input_file $base --scalar_field $scalar_field
      if ($?) exit 1
    endif
  endif
end

############################################
#  split files into one variable per file
############################################

set cmip = ""
if ($subDir == "cmip") set cmip = "-c"

foreach filetype ($filetypes)
  set comp = $filetype:r:r:e
  if (! -e $comp && $comp != "grid_spec") then
    echo "Spliting $comp files"
    mkdir $comp
    $SPLIT_NCVARS $cmip -l -VVV -o $comp `/bin/ls ????0101.$comp.nc` > $comp.out
  endif
end


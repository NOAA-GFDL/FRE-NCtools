#!/bin/csh -f

# =========================================================================
#   script creates a simple hydrography, with no manual intervention,
#     simplified lake fields, no elevation
#    -- run river_regrid tool
#    -- run post-processor on river_regrid output
#    -- remove parallel rivers
#    -- run post-processor again
#    -- add GLCC waterbod fractions
# =========================================================================

#set echo on

source $MODULESHOME/init/csh

set argv = (`getopt  hc:o:f:t:m: $*`)

while ("$argv[1]" != "--")
    switch ($argv[1])
        case -h:
            set help; breaksw
        case -o:
            set outdir = $argv[2]; shift argv; breaksw
        case -c:
            set compile = $argv[2]; shift argv; breaksw
        case -f:
            set min_frac = $argv[2]; shift argv; breaksw
        case -t:
            set land_thresh = $argv[2]; shift argv; breaksw
        case -m:
            set mosaic_file = $argv[2]; shift argv; breaksw
    endsw
    shift argv
end
shift argv

# argument error checking

if (! $?min_frac) then
   echo "ERROR: no argument given for min_frac (river_regrid)."
   set help
endif

if (! $?land_thresh) then
   echo "ERROR: no argument given for land_thresh (river_regrid)."
   set help
endif

if (! $?mosaic_file) then
   echo "ERROR: no argument given for mosaic_file."
   set help
endif

# usage message
if ($?help) then
   set cmdname = `basename $0`
   echo
   echo "USAGE:  $cmdname [-f min_frac] [-t land_thresh] [-m mosaic] [-c compile] [-o outdir]"
   echo
   echo "         -f min_frac:    required, river-regrid, min land fraction; set to 0 to retain small "
   echo "                         land fractions from gridspec"
   echo "         -t land_thresh: required, river-regrid, land fraction exceeding (1-land_thresh) is set to 1;"
   echo "                         recommend 1.e-5"
   echo "         -m mosaic:      required, mosaic file location"
   echo "         -c compile:     use '-c compile' to compile code (optional)"
   echo "         -o outdir:      output directory (optional)"
   echo
   exit 1
endif

set riv_regrd = river_regrid

# set file for GLCC data (waterbod)
set glcc_file = gigbp2_0ll.nc

# set file for the disaggregated, extended river network

set river_network_ext_file = river_network_fill_coast.nc

if (! -e OUTPUT) then
      mkdir -p OUTPUT/{river_regrid,post_regrid,rmv_parallel_rivers,post_rmvp}
endif

# ------------------------------------------------
#  RUN RIVER_REGRID
# ------------------------------------------------
echo ""
echo RUN RIVER-REGRID
$riv_regrd --mosaic $mosaic_file --river_src $river_network_ext_file --min_frac $min_frac --land_thresh $land_thresh
mv river_output*nc OUTPUT/river_regrid/

# ------------------------------------------------
#  POST-PROCESS OUTPUT FROM RIVER_REGRID
# ------------------------------------------------
set river_input_files = OUTPUT/river_regrid/river_output*nc 
echo $#river_input_files > fort.5
foreach file ($river_input_files)
   echo OUTPUT/river_regrid/$file:t >> fort.5
end
cho ""
echo RUN POST-PROCESSOR
cp_river_vars < fort.5
if ($status != 0) then
    echo ERROR in post-processing river-regrid output, exiting...
    exit
else
    mv river_network*nc OUTPUT/post_regrid/
    mv out.cp_river_vars OUTPUT/post_regrid/
endif

# ------------------------------------------------
#  REMOVE PARALLEL RIVERS
# ------------------------------------------------
set add_ocean_const = F
set river_input_files = OUTPUT/post_regrid/river_network*nc
echo $#river_input_files > fort.5
foreach file ($river_input_files)
   echo OUTPUT/post_regrid/$file:t >> fort.5
end
echo $add_ocean_const >> fort.5
echo ""
echo REMOVE PARALLEL RIVERS
rmv_parallel_rivers < fort.5
if ($status != 0) then
    echo ERROR in removal of parallel rivers, exiting...
    exit
else
    mv river_network*nc OUTPUT/rmv_parallel_rivers/
    mv out.rmv_parallel_rivers OUTPUT/rmv_parallel_rivers/
endif

# ------------------------------------------------
#  POST-PROCESS OUTPUT FROM REMOVE-PARALLEL-RIVERS
# ------------------------------------------------
set river_input_files = OUTPUT/rmv_parallel_rivers/river_network*nc
echo $#river_input_files > fort.5
foreach file ($river_input_files)
   echo OUTPUT/rmv_parallel_rivers/$file:t >> fort.5
end
echo ""
echo RUN POST-PROCESSOR
cp_river_vars < fort.5
if ($status != 0) then
    echo ERROR in post-processing parallel-river output, exiting...
    exit
else
    mv river_network*nc OUTPUT/post_rmvp/
    mv out.cp_river_vars OUTPUT/post_rmvp/
endif

# --------------------------------------
#  ADD GLCC WATERBOD DATA
# --------------------------------------
set travel_thresh = 2.
set river_input_files = OUTPUT/post_rmvp/river_network*nc
echo $#river_input_files > fort.5
foreach file ($river_input_files)
   echo OUTPUT/post_rmvp/$file:t >> fort.5
end
echo $glcc_file:t >> fort.5
echo $travel_thresh >> fort.5
touch input.nml
echo ""
echo ADD LAKES
cr_lake_files < fort.5
if ($status != 0) then
    echo ERROR in addition of lake data, exiting...
    exit
endif

@ k = 0
while ($k < $#river_input_files)
  @ k = $k + 1
  mv OUTPUT/post_rmvp/river_network.tile"$k".nc OUTPUT/hydrography.tile"$k".nc
end
set hydro_files = OUTPUT/hydrography*.nc
foreach file ($hydro_files)
   set tn = `echo "$file" | awk -Ftile '{print $2}'`
   ncks -A -a -v lake_frac,lake_depth_sill,lake_tau,WaterBod,PWetland,connected_to_next,whole_lake_area,max_slope_to_next \
     lake_frac.tile"$tn" "$file:t"
end



#!/usr/bin/env perl
#
# ------------------------------------------------------------------------------
# FMS/FRE Project: Script to Split netCDF Files
# ------------------------------------------------------------------------------
# Copyright (C) NOAA Geophysical Fluid Dynamics Laboratory, 2000-2019
# Designed and written by V. Balaji, Amy Langenhorst and Aleksey Yakovlev
#
use strict;
use Cwd;
use Getopt::Long;
use List::MoreUtils qw{uniq};
Getopt::Long::Configure("bundling");
use File::Path;
use File::Spec;
use File::Basename;

my $cwd = getcwd;
my $ncstatus = 0;
my $TEST = 0;
my $tmp_var_filename = ".split_ncvars.$$.var.nc";
my $VERSION = '20190813';
my %Opt = ( HELP=>0, VERBOSE=>0, QUIET=>0, LOG=>0, STATIC=>0, CMIP=>0, PS=>0, AUTO=>0, odir=>$cwd );

#  ----- parse input argument list ------

my $status = GetOptions ('h|help!'     => \$Opt{HELP},
                         'V|verbose+'  => \$Opt{VERBOSE},
                         'q|quiet!'    => \$Opt{QUIET},
                         'l|log!'      => \$Opt{LOG},
                         's|static!'   => \$Opt{STATIC},
                         'c|cmip!'     => \$Opt{CMIP},
			 'p|ps'        => \$Opt{PS},
                         'a|auto!'     => \$Opt{AUTO},
                         'f|onefile=s' => \$Opt{onefile},
                         'i|idir=s'    => \$Opt{idir},
                         'o|odir=s'    => \$Opt{odir},
                         'v|vars=s'    => \$Opt{vars});

usage() if ($Opt{HELP} || @ARGV == 0);
my @ifiles = @ARGV;
my $first = 1;
my $odir = $Opt{odir};

##################################################################
# external scripts 

my $ncrcat = `which ncrcat 2>&1`;
die "Unable to locate ncrcat, cannot continue.\n" if $ncrcat =~ /(: no)|(not found)/;
chomp $ncrcat; $ncrcat .= " --64 -t 2 --header_pad 16384";

my $ncks = `which ncks 2>&1`; 
die "Unable to locate ncks, cannot continue.\n" if $ncks =~ /(: no)|(not found)/;
chomp $ncks; 
# NCO changed options between 4.5.4 and 4.5.5 and 4.6. Hopefully they stopped those schenanigans.
my $tmp = qx{$ncks --version 2>&1};
$ncks .=  ($tmp =~ /4\.5\.4/ ? ' --64bit'
	   : $tmp =~ /4\.5\.5/ ? ' --64bit_offset'
	   : ' --fl_fmt=64bit_offset')
       . ' --header_pad 16384';

my $list_ncvars = `which list_ncvars.csh 2>&1`;
die "Unable to locate list_ncvars.csh, cannot continue.\n" if $list_ncvars =~ /(: no)|(not found)/;
chomp $list_ncvars;

my $ncdump = `which ncdump 2>&1`;
die "Unable to locate ncdump, cannot continue.\n" if $ncdump =~ /(: no)|(not found)/;
chomp $ncdump;

my $dmget_available = `which dmget 2>&1`;
$dmget_available = 0 if $dmget_available =~ /(: no)|(not found)/;

##################################################################
# user supplied variable list

my @varlist;
if ($Opt{vars}) {
    @varlist = split /,/,$Opt{vars};
}
else {
    # vars list must be set to one file option
    if ($Opt{onefile}) {
        $Opt{onefile} = "";
    }
}

##################################################################
#  need to make output directory

File::Path::make_path($odir) unless -d $odir;

print "split_ncvars.pl version $VERSION\n" if $Opt{VERBOSE};
#  process each input file separately

my %ps_includes;
my $file_copied = 0;
foreach my $file (@ifiles) {
    print STDERR "Cowardly refusing to process files without a .nc filename" and next
	unless $file =~ /\.nc/;
    my ($ext) = (basename $file) =~ /\.nc(.*)/;
    my ($tile) = (basename $file) =~ /(\.tile\d)/;
    if ($Opt{idir} ) {
	if ($dmget_available
	    and File::Spec->rel2abs($Opt{idir}) =~ m@^/arch@) {
	    print  "dmget ".$Opt{idir}."/$file\n" if $Opt{VERBOSE} > 2;
	    system("dmget ".$Opt{idir}."/$file");

	    if ( $Opt{idir} ne '.' ) {
		print  "gcp ".$Opt{idir}."/$file $file\n" if $Opt{VERBOSE} > 2;
		system("gcp ".$Opt{idir}."/$file $file");
		$file_copied++;
	    }
	}
	$file = "$Opt{idir}/$file";
    }

    print STDERR "Unable to read $file, skipping\n" and next unless -r $file;
    my $dump = qx{$ncdump -h $file};

    # generate variable list (if not user supplied)
    if (!@varlist) {
        if (!$Opt{STATIC}) {
	    print tailname($list_ncvars)." -t0123 $file\n" if $Opt{VERBOSE} > 1;
	    @varlist = split /\n/, `$list_ncvars -t0123 $file`;
        } else {
	    print tailname($list_ncvars)." -s0123 $file\n" if $Opt{VERBOSE} > 1;
	    @varlist = split /\n/, `$list_ncvars -s0123 $file`;
        }
    }

    # get time axis name
    my $timename = get_time_dimension($dump) if !$Opt{STATIC};

    #Balaji: add all variables that appear in a coordinate attribute
    #    set coords = ( `ncdump -h $file |grep :coordinates |cut -f2 -d\" |sed "s/ /\n/g" |sort -u` \
    #                   `ncdump -h $file |grep "float geol[oa][nt]" | awk '{print $2}' |cut -f1 -d\(` )
    #    set coords = `echo $coords |sed "s/ /\n/g" |sort -u`
    #     exec echo coords = $coords
    my @coords;
    foreach my $coord (qw/geolat geolon GEOLAT GEOLON/) {
	push @coords, $coord if ($dump =~ /\t\w+ $coord\(.+\)/);
    }
#-------------------------------------
#-----  loop through variables  ------
#-------------------------------------

     print "processing ".scalar(@varlist)." variables\n\n" if $Opt{VERBOSE} > 2;

     my @vlist = @coords;
     my @xlist = ();

     foreach my $var (@varlist) {
         $var =~ s/^\s+//; $var =~ s/\s+$//;

         # skip this variable if it does not exist
         if ($dump !~ /\t\w+ $var\(.+\)/) {
           print "WARNING: variable $var does not exist ... skipping.\n" if !$Opt{QUIET};
           next;
         }
         print "Processing variable: $var\n" if $Opt{VERBOSE} > 1;
         print "Size of vlist = ".scalar(@vlist)."  (INITIAL)\n" if $Opt{VERBOSE} > 3;

         # add variable to list of variables to extract
         push @vlist, $var;

         # Add the cell_measures to the list of variables to include in the external_variables global attribute
         my @cell_measures_variables = get_cell_measures($dump,$var);
         if (@cell_measures_variables) {
            push @xlist, @cell_measures_variables;
            print "   $var: external_variables = ".join(", ",@cell_measures_variables)."\n" if $Opt{VERBOSE} > 1;
         }

         # automatically detect if this a a cmip variable - if attrib standard_name exists
         my $CMIP = $Opt{CMIP};
         if ($Opt{AUTO} && !$Opt{CMIP}) {
           $CMIP = 1 if (get_variable_att($dump,$var,"standard_name"));
         }

	 # get the number of dimensions if ps is provided on the commandline
	 if ( $Opt{PS} ){
	     my $dimensions = get_variable_dimensions($dump,$var);
	     if ( grep ( /(pfull|phalf)/, @{$dimensions} ) and scalar @{$dimensions} >= 4 ) {
		 print "Setting ps_incluides for $var\n" if $Opt{VERBOSE} > 1;
		 $ps_includes{$var} = 1;
	     } else {
		 print "Not setting ps_includes of $var\n" if $Opt{VERBOSE} > 1;
		 $ps_includes{$var} = 0;
	     }
	 }

	 my $output = "$odir/$var$tile.nc";
	 my $output_barefile = "$var$tile.nc";
        #------------------------------------------------------------
        # FIRST TIME ONLY (when output file does not exist)
        # need to extract additional static fields
        # 1) dimension bounds, 2) coordinate variables, 3) formula terms
        #------------------------------------------------------------
         if (!-e $output) {

            # get variable names of all dimension bounds/edges (except for time)
            my @bounds = get_variable_bounds   ($dump,$var,$CMIP);
            if (@bounds) {
              print "   $var: bounds = ".join(", ",@bounds)."\n" if $Opt{VERBOSE} > 1;
              push @vlist, @bounds;
            }

            # get variables specified by the coordinates attribute
            my @coordinates = get_variables_from_att($dump,$var,"coordinates");
            if (@coordinates) {
              print "   $var: coordinates = ".join(", ",@bounds)."\n" if $Opt{VERBOSE} > 1;
              push @vlist, @coordinates;
              # add bounds for coordinate variables
              foreach my $coord (@coordinates) {
                push @vlist, get_dimension_bounds ($dump,$coord);
              }
            }

	    

            # formula terms
            # only add static terms to the vlist
            # time-varying terms are "external_variables"
            my ($vterms,$xterms) = get_formula_terms($dump,$var);
            if (@$vterms) {
              print "   $var: formula_terms(static) = ".join(", ",@$vterms)."\n" if $Opt{VERBOSE} > 1;
              push @vlist, @$vterms;
            }
            if (@$xterms) {
              print "   $var: external_variables = ".join(", ",@$xterms)."\n" if $Opt{VERBOSE} > 1;
              push @xlist, @$xterms;
            }
         }
	 if ( $ps_includes{$var} ){
	   print "Appending ps,pk,bk to variables to include.\n" if $Opt{VERBOSE} > 1;
	   push @vlist, qw{ps pk bk};
         }
        #-------------------------------------------------
        # EVERY TIME (when there is a time dimension)
        # need to extract time bounds/average_XX variables
        #-------------------------------------------------
         if ($timename) {

            # fms-style time avg variables
            if (!$CMIP) {
              my @time_avg_info = get_time_avg_info($dump,$var);
              if (@time_avg_info) {
                print "   $var: time_avg_info = ".join(", ",@time_avg_info)."\n" if $Opt{VERBOSE} > 1;
                push @vlist, @time_avg_info;
              }
            }

           # CF compliant time avg variables (bounds or climatology)
            my $SKIPVAR = 0;
            foreach my $bounds_type ( qw/bounds climatology/ ) {
               my $bounds_name = get_variable_att($dump,$timename,$bounds_type);
               if ($bounds_name) {
                  if ($bounds_name eq $var) {
                    # skip this variable because it is a time bounds variable
                    print "   $var: SKIPPING because it is a time bounds variable\n" if $Opt{VERBOSE} > 1;
                    $SKIPVAR = 1;
                    last;
                  }
                  print "   $var: $timename: $bounds_type = $bounds_name\n" if $Opt{VERBOSE} > 1;
                  push @vlist, $bounds_name;
               }
            }

            # skip this variable, reset the variable list
            if ($SKIPVAR) {
              @vlist = @coords;
              @xlist = ();
              next;
            }
         }

         # if a single file is desired then do not extract into single-field files
         next if $Opt{onefile};

     #-----------------------
     #--- extract fields ----
     #-----------------------
         print "Size of vlist = ".scalar(@vlist)."  (NCKS)\n" if $Opt{VERBOSE} > 3;
         my $vlist = join ",",@vlist;
         print "   var=$var; timename=$timename; vlist=$vlist\n" if $Opt{VERBOSE} > 2;
         my $appendopt = "";
         $appendopt = "-A" if (-f "$cwd/$tmp_var_filename");
         print   "$ncks -h $appendopt -v $vlist $file $tmp_var_filename\n" if $Opt{VERBOSE};
         next if $TEST;
         system ("$ncks -h $appendopt -v $vlist $file $tmp_var_filename");
         $ncstatus += $?;

         # remove dimensions called "scalar_axis" (i.e., length = 1)
         # these are typically in the near-surface field files
         my @commands = remove_degenerate_dimension("$tmp_var_filename","scalar_axis");
         foreach my $cmd (@commands) {
           print "$cmd\n" if $Opt{VERBOSE} > 0;
           system($cmd);
         }

     #--- write to the output file ---
         if (-e $output) {
             die "ERROR: try to append to unlimited dimension when no unlimited dimension found" if !$timename;

            #--- append to existing file ---
            #push @commands, "dmget ".$Opt{odir}."/$var.nc";
             print  "mv $output $odir/.split_ncvars.$$.tmp.nc\n" if $Opt{VERBOSE} > 2;
             system("mv $output $odir/.split_ncvars.$$.tmp.nc");
             print "ncrcat $odir/.split_ncvars.$$.tmp.nc $tmp_var_filename $output\n" if $Opt{VERBOSE} > 0;
             system ("$ncrcat -h $odir/.split_ncvars.$$.tmp.nc $tmp_var_filename $output");
             $ncstatus += $?;
             print  "rm -f $odir/.split_ncvars.$$.tmp.nc $tmp_var_filename\n" if $Opt{VERBOSE} > 2;
             system("rm -f $odir/.split_ncvars.$$.tmp.nc $tmp_var_filename");
         } else {
            #--- create new file ---
            #--- modify filename attribute ---
            #--- move to output directory ---

             my @ncatted_opts = set_ncatted_opts("$tmp_var_filename","$output_barefile",$var,$CMIP);
             # external_variables attribute (cell_measures & time-varying formula terms)
             if (@xlist) {
               push @ncatted_opts, "-a external_variables,global,c,c,\"".join(" ",@xlist)."\"";
             }

             # rename coordinate variables
             if (!$Opt{onefile}) {
               my $vdump = `ncdump -h $tmp_var_filename`;
               my @coordinates = get_variables_from_att($vdump,$var,"coordinates");

               # rename "height#" (coordinate) variables to just "height"
               my @height = grep{/height\d+/} @coordinates;
               if (@height == 1) {
                 print  "ncrename -h -v $height[0],height $tmp_var_filename\n";
                 system("ncrename -h -v $height[0],height $tmp_var_filename");
                 for (@coordinates) {
                   s/$height[0]/height/;
                 }
                 push @ncatted_opts, "-a coordinates,$var,m,c,\"@coordinates\"";
               } elsif (@height > 1) {
                 die "ERROR: can not have more than one height coordinate attribute/variable";
               }

               # Do some extra CMORizing if --cmip
               if ($CMIP) {
                 # rename "p###" (coordinate) variables to just "plev"
                 my @plev = grep{/p\d+/} @coordinates;
                 push @plev, grep{/pl700/} @coordinates;
                 if (@plev == 1) {
                   if (get_variable_att($vdump,$plev[0],"units") eq "Pa") {
                     print  "ncrename -h -v $plev[0],plev $tmp_var_filename\n";
                     system("ncrename -h -v $plev[0],plev $tmp_var_filename");
                     for (@coordinates) {
                       s/$plev[0]/plev/;
                     }
                     push @ncatted_opts, "-a coordinates,$var,m,c,\"@coordinates\"";
                   }
                   # also rename the bounds variable/attribute (if there is one)
                   my $plev_bnds = get_variable_att($vdump,$plev[0],"bounds");
                   if ($plev_bnds) {
                     if ($vdump =~ /\t\w+ $plev_bnds\(.+\)/) {
                       print  "ncrename -h -v $plev_bnds,plev_bnds $tmp_var_filename\n";
                       system("ncrename -h -v $plev_bnds,plev_bnds $tmp_var_filename");
                       push @ncatted_opts, "-a bounds,plev,m,c,\"plev_bnds\"";
                     }
                   }
                 } elsif (@plev > 1) {
                   die "ERROR: can not have more than one plev coordinate attribute/variable";
                 }

                 # rename vertical dimensions (plev## -> plev, levhalf -> lev)
                 foreach my $dim (get_vertical_dimensions($vdump,$var)) {
                   if ($dim =~ /^plev\d+/) {
                     print  "ncrename -h -d $dim,plev -v $dim,plev $tmp_var_filename\n";
                     system("ncrename -h -d $dim,plev -v $dim,plev $tmp_var_filename");
                   } elsif ($dim eq "levhalf") {
                     print  "ncrename -h -d $dim,lev -v $dim,lev $tmp_var_filename\n";
                     system("ncrename -h -d $dim,lev -v $dim,lev $tmp_var_filename");
                     # also rename the formula terms
                     if (get_variable_att($vdump,$dim,"formula_terms") eq "ap: ap_half b: b_half ps: ps") {
                       print  "ncrename -h -v ap_half,ap -v b_half,b $tmp_var_filename\n";
                       system("ncrename -h -v ap_half,ap -v b_half,b $tmp_var_filename");
                       push @ncatted_opts, "-a formula_terms,lev,m,c,\"ap: ap b: b ps: ps\"";
                     }
                   }
                 }
               }
             }

#	     my $destination = "$odir/$var.nc$ext"; #  . ($Opt{UNCOMB} ? $ext : '');
	     my $destination = "$output$ext";
	     print  "mv $tmp_var_filename $destination\n" if $Opt{VERBOSE} > 2;
	     system("mv $tmp_var_filename $destination");
	     
	     if (@ncatted_opts) {
		 print  "ncatted -h -O @ncatted_opts $destination\n" if $Opt{VERBOSE} > 1;
		 system("ncatted -h -O @ncatted_opts $destination");
	     }
         }

         # reset the variable lists
         @vlist = @coords;
         @xlist = ();
     }   ### end of variable loop ###

    #--- extract fields for single file option ---
     if ($Opt{onefile}) {
         my $vlist = join ",",@vlist;
         print  "$ncks -h -A -v $vlist $file $Opt{onefile}\n" if $Opt{VERBOSE} > 0;
         system ("$ncks -h -A -v $vlist $file $Opt{onefile}");
         $ncstatus += $?;

         my @ncatted_opts = set_ncatted_opts($Opt{onefile},tailname($Opt{onefile}));
         # external_variables attribute (missin & time-varying formula terms)
         if (@xlist) {
            push @ncatted_opts, "-a external_variables,global,c,c,\"".join(" ",@xlist)."\"";
         }
         if (@ncatted_opts) {
           print  "ncatted -h -O @ncatted_opts ".$Opt{onefile}."\n" if $Opt{VERBOSE} > 1;
           system("ncatted -h -O @ncatted_opts ".$Opt{onefile})
         }
     }

    #--- create readme ---
     if ($first && $Opt{LOG}) {
        variable_log($file,"$odir/README") if !$TEST;
        $first = 0;
     }
    # remove file if copied from archive
     unlink $file if $file_copied;

    # clean up ncks dropping
    unlink "ncks.out" if (-e "ncks.out");

 }   ### end of file loop ###

exit 1 if $ncstatus;

##################################################################
##################################################################
##################################################################
#  ----- help message -----

sub usage {
  my $name = substr($0,rindex($0,"/")+1);
  print "
Split NetCDF Variables, version $VERSION

Usage:  $name [-d] [-l] [-i idir] [-o odir] [-f file] [-v vars]  files.....

        -l       = include readme log in output directory
        -V       = verbose option: increases standard output, multiple -V allowed
        -c       = cmip option, time average info variable not written
        -p       = includes ps variable in output files that will be zInterpolated
        -i idir  = input (archive) directory path
        -o odir  = output directory path
        -u file  = option to operate on uncombined files
        -f file  = one file output option
                   file is the name of the output file
                   this option must be used with -v option
        -v vars  = comma separated list of variables
        files... = list netcdf files in directory idir
                   the files must be in chronological order

";
  exit 1;
}

#------------------------------------------------

sub get_time_dimension {
   my $dump = shift;
   my $timeName;
   my $timeSize;
   if ($dump =~ /\t(.+) = UNLIMITED ; \/\/ \((\d+) currently\)/) {
      $timeName = $1; 
      $timeSize = $1; 
   }   
   return $timeName;
}

#-------------------------------------------

sub get_variable_dimensions {
  my $dump = shift;
  my $var = shift;
  my %cartesian_coords;
  my @coords;
  # Get a list of the cartesian coordinates that are in a file
  while ( $dump =~ /\t\t(.*):(cartesian_axis|axis) = "(.*)"/g ){
      $cartesian_coords{$1} = $2;
  }
  # Compare the cartesian coords to those of the variables
  if ( $dump =~ /\t.*$var\((.+)\)/){
    @coords = uniq map { grep $_, split /,\s/, $1 } keys %cartesian_coords;
  }
  return \@coords;
}

#-------------------------------------------
# return variable dimensions with the attribute
#   axis or cartesian_axis ="Z"

sub get_vertical_dimensions {
  my $dump = shift;
  my $var = shift;
  my @vdims;
  if ( $dump =~ /\t.*$var\((.+)\)/){
    foreach my $dim (split /,\s/, $1) {
      if ($dump =~ /\t\w+ $dim\(.+\)/) {
        push @vdims, $dim if ($dump =~ /\t\t$dim:(axis|cartesian_axis) = "Z" ;/);
      }
    }
  }
  return @vdims;
}

#-------------------------------------------

sub get_variable_att {
  my $dump = shift;
  my $var = shift;
  my $att = shift;
  my $value;
  if ($dump =~ /\t\t$var:$att = "(.+)" ;/) {
    $value = $1;
  }   
  return $value;
}

#-------------------------------------------

sub get_time_avg_info {
  my $dump = shift;
  my $var = shift;
  my @info = split /,/, get_variable_att($dump,$var,"time_avg_info");
  return @info;
}

#-------------------------------------------

sub get_variables_from_att {
  my $dump = shift;
  my $var = shift;
  my $att = shift;
  my @attvars;
  foreach (split /\s+/, get_variable_att($dump,$var,$att)) {
    push @attvars, $_ if ($dump =~ /\t\w+ $_\(.+\)/ || $dump =~ /\t\w+ $_ ;/);
  }
  return @attvars;
}

#-------------------------------------------

sub get_cell_measures {
  my $dump = shift;
  my $var  = shift;
  my @cell_measures;
  if ( $dump =~ /\t\w+ $var\((.+)\)/ ){
    my $cell_measures_att = get_variable_att($dump,$var,"cell_measures");
    # e.g. volume: volcello area: areacello
    @cell_measures = grep { !/:$/ } split ' ', $cell_measures_att;
  }
  return @cell_measures;
}

#-------------------------------------------

sub get_variable_bounds {
  my $dump = shift;
  my $var = shift;
  my $CMIPVAR = shift; # cmip flag
  my $timename = get_time_dimension($dump);
  my @bounds;
  # find dimensions for this variable
  if ($dump =~ /\t\w+ $var\((.+)\)/ ) {
    my @axes = split /, /, $1;
    foreach my $dim (@axes) {
      next if ($dim eq $timename);
      push @bounds, get_dimension_bounds($dump,$dim,$CMIPVAR);
    }     
  }     
  return @bounds;
}

#-------------------------------------------

sub get_dimension_bounds ($$;$) {
  my ($dump,$dim,$CMIPVAR) = @_;
  my @bounds;
  # get the field pointed to by bounds attribute
  my $bnds = get_variable_att($dump,$dim,"bounds");
  push @bounds, $bnds if ($dump =~ /\t\w+ $bnds\(.+\)/);
  # get the field pointed to by edges attribute
  my $bnds = get_variable_att($dump,$dim,"edges");
  push @bounds, $bnds if ($dump =~ /\t\w+ $bnds\(.+\)/);
  return @bounds;
}

#-------------------------------------------
# find the formula_terms for this variable
# return them separately as formula_terms (fterms) and external_terms (xterms)
# external_terms are variables not found in the file and (for now) non-static variables

sub get_formula_terms {
  my $dump = shift;
  my $var = shift;

  #--- first find all formula terms ---
  my @allterms;
  # find dimensions for this variable
  if ($dump =~ /\t\w+ $var\((.+)\)/ ) {
    my @axes = split /, /, $1;
    # check each axis for formula terms
    foreach my $dim (@axes) {
      my $form = get_variable_att($dump,$dim,"formula_terms");
     #print "   $dim:formula_terms = \"$form\"\n" if $form;
      foreach my $term ($form =~ /\w+:\s+(\w+)/g) {
        push @allterms, $term if !grep{$_ eq $term} @allterms;
      }
      # also check bounds for formula terms
      my $bnds = get_variable_att($dump,$dim,"bounds");
      my $form = get_variable_att($dump,$bnds,"formula_terms");
     #print "   $bnds:formula_terms = \"$form\"\n" if $form;
      foreach my $term ($form =~ /\w+:\s+(\w+)/g) {
        push @allterms, $term if !grep{$_ eq $term} @allterms;
      }
    }
  }

  #--- now split into formula terms and external terms ---
  my (@fterms,@xterms);
  foreach my $term (@allterms) {
    if (is_static($dump,$term) && $dump =~ /\t\w+ $term\((.+)\)/) {  # is static and variable
      push @fterms, $term;
    } else {
      push @xterms, $term;
    }
  }

  return (\@fterms,\@xterms);
}

#-------------------------------------------

sub is_static {
  my ($dump,$var) = @_;
  my $timename = get_time_dimension($dump);
  my $static = 0;
  # find dimensions for this variable
  if ($dump =~ /\t\w+ $var\((.+)\)/ ) {
    my @axes = split /, /, $1;
    $static = !grep{$_ eq $timename} @axes;
  }
  return $static;
}
  
#-------------------------------------------

sub tailname {
   my $tail = shift;
   while ($tail =~ s/^.*\///) {}
   return $tail;
}

#-------------------------------------------
# cleanup attributes in the output file

sub set_ncatted_opts ($$;$$) {
  my ($file,$filename,$var,$CMIPVAR) = @_;
  my @opts;
  my $dump = `ncdump -h $file`;
  push @opts, "-a filename,global,m,c,\"$filename\"" if ($dump =~ /\t\t:filename = ".+" ;/);

  my $svar = "\w+"; # globally edit attributes
  $svar = $var if $var;

  if ($CMIPVAR) {
    push @opts, "-a time_avg_info,$var,d,," if ($dump =~ /\t\t$svar:time_avg_info = ".+" ;/);

    # loop over dimensions for this variable
    if ($dump =~ /\t\w+ $var\((.+)\)/ ) {
      my @axes = split /, /, $1;
      foreach my $dim (@axes) {
        push @opts, cleanup_dim_atts($dump,$dim);
        # also check bounds dimensions
        my $bnds = get_variable_att($dump,$dim,"bounds");
        if ($bnds) {
          push @opts, cleanup_dim_atts($dump,$bnds);
          if ($dump =~ /\t\w+ $bnds\((.+)\)/) {
            foreach my $dim2 (split /, /, $1) {
              push @opts, cleanup_dim_atts($dump,$dim2) if (!grep{$_ eq $dim2} @axes);
            }
          }
        }
      }
    }
  }
  # remove redundant values
  my %seen;
  my @unique = grep{!$seen{$_}++} @opts;
  return @unique;
}

#-------------------------------------------

sub cleanup_dim_atts {
  my $dump = shift;
  my $dim = shift;
  my @opts;

  # correct (dimension) attribute "cartesian_axis -> axis"
  my $cart_axis = get_variable_att($dump,$dim,"cartesian_axis");
  if ($cart_axis) {
    if ($cart_axis ne "N") {
      push @opts, "-a axis,$dim,c,c,\"$cart_axis\"" if !get_variable_att($dump,$dim,"axis");
    }
    push @opts, "-a cartesian_axis,$dim,d,,";
  }

  # correct (dimension) attribute "calendar_type -> calendar"
  my $cal = get_variable_att($dump,$dim,"calendar");
  my $calt = get_variable_att($dump,$dim,"calendar_type");
  if ($cal) {
     my $calc = lc($cal);
     push @opts, "-a calendar,$dim,m,c,\"$calc\"" if ($cal ne $calc);
  } elsif ($calt) {
     my $calc = lc($cal);
     push @opts, "-a calendar,$dim,c,c,\"$calc\"";
  }
  push @opts, "-a calendar_type,$dim,d,," if $calt;

  # correct units from "none" -> "1"
  my $units = get_variable_att($dump,$dim,"units");
  if ($units) {
    push @opts, "-a units,$dim,m,c,\"1\"" if (lc($units) eq "none");
  }

  # add standard name for lat/lon dimensions
  foreach my $name (qw/latitude longitude/) {
    if ($dim eq $name || $dim eq substr($name,0,3)) {
      my $stdname = get_variable_att($dump,$dim,"standard_name");
      push @opts, "-a standard_name,$dim,c,c,\"$name\"" if (!$stdname);
    }
  }

  return @opts;
}

#-------------------------------------------

sub remove_degenerate_dimension {
  my ($file,$dim) = @_;
  my $dump =`ncdump -h $file`;
  my @cmds;
  # is this a singleton dimension?
  if ($dump =~ /\t$dim = 1 ;/) {
    # remove dim by averaging, then (un)extract the unused dimension
    push @cmds, "ncwa -h -a $dim $file $file.ztmp";
    push @cmds, "ncks -h -O -x -v scalar_axis $file.ztmp $file";
    push @cmds, "rm -f $file.ztmp";
  }
  return @cmds;
}

#-------------------------------------------

sub variable_log {
  my $ncfile = shift;
  my $prtfile = shift;
  my $dump = `ncdump -h $ncfile`;

  open (OUT,"> $prtfile") || die "Cannot open $prtfile for output";

  foreach my $var (split /\n/, `$list_ncvars -t0123 $ncfile`) {
    my @out;
    $var =~ s/^\s+//; $var =~ s/\s+$//;
    push @out, $var;
    push @out, get_variable_att($dump,$var,"long_name");
    push @out, get_variable_att($dump,$var,"units");
    printf OUT "%20s = %s (%s)\n", @out;
  }

  close(OUT);
  return 0;
}


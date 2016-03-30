#!/usr/bin/perl
#
# $Id$
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
use strict;
use Cwd;
use Getopt::Long;
Getopt::Long::Configure("bundling");

my $cwd = getcwd;
my $ncstatus = 0;
my $TEST = 0;

my %Opt = ( HELP=>0, VERBOSE=>0, QUIET=>0, LOG=>0, STATIC=>0, CMIP=>0, odir=>$cwd );

#  ----- parse input argument list ------

my $status = GetOptions ('h|help!'     => \$Opt{HELP},
                         'V|verbose+'  => \$Opt{VERBOSE},
                         'q|quiet!'    => \$Opt{QUIET},
                         'l|log!'      => \$Opt{LOG},
                         's|static!'   => \$Opt{STATIC},
                         'c|cmip!'     => \$Opt{CMIP},
                         'f|onefile=s' => \$Opt{onefile},
                         'i|idir=s'    => \$Opt{idir},
                         'o|odir=s'    => \$Opt{odir},
                         'v|vars=s'    => \$Opt{vars});

#if ($Opt{DEBUG}) set echo
usage () if ($Opt{HELP} || @ARGV == 0);
my @ifiles = @ARGV;
my $first = 1;
my $odir = $Opt{odir};

##################################################################

my $ncrcat = `which ncrcat`; chomp $ncrcat; $ncrcat .= " --64bit -t 2 --header_pad 16384";
my $ncks   = `which ncks`; chomp $ncks; $ncks .= " --64bit --header_pad 16384";
my $list_ncvars = `which list_ncvars.csh`; chomp $list_ncvars;

##################################################################
# user supplied variable list

  my @varlist;
  if ($Opt{vars}) {
     @varlist = split /,/,$Opt{vars};
  } else {
     # vars list must be set to one file option
     if ($Opt{onefile}) {
        $Opt{onefile} = "";
     }
  }

##################################################################

#  need to make output directory

 system("mkdir -p $odir") if (!-e $odir);

#  process each input file separately

 foreach my $file (@ifiles) {
     if ($Opt{idir}) {
        my @commands;
        print  "dmget ".$Opt{idir}."/$file\n" if $Opt{VERBOSE} > 2;
        system("dmget ".$Opt{idir}."/$file");
        print  "gcp ".$Opt{idir}."/$file $file\n" if $Opt{VERBOSE} > 2;
        system("gcp ".$Opt{idir}."/$file $file");
     }
     my $dump = `ncdump -h $file`;

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
     foreach my $coord (qw/geolat geolon/) {
       push @coords, $coord if ($dump =~ /\t\w+ $coord\(.+\)/);
     }

#### loop through variables ####
     print "processing ".scalar(@varlist)." variables\n\n" if $Opt{VERBOSE} > 2;
     my @vlist;

     foreach my $var (@varlist) {
         @vlist = @coords;
         my @xlist;

         $var =~ s/^\s+//; $var =~ s/\s+$//;
         # skip this variable if it does not exist
         if ($dump !~ /\t\w+ $var\(.+\)/) {
           print "WARNING: variable $var does not exist ... skipping.\n" if !$Opt{QUIET};
           next;
         }
         print "Processing variable: $var\n" if $Opt{VERBOSE} > 1;

         # add variable to list of variables to extract
         push @vlist, $var;

     #--- first time: need to extract additional static fields ---
         if (!-e "$odir/$var.nc") {
            # get axis names then edge names (may want to get bounds for CF compliance)
            my @bounds = get_variable_bounds   ($dump,$var);
            if (@bounds) {
              print "   $var: bounds = ".join(", ",@bounds)."\n" if $Opt{VERBOSE} > 1;
              push @vlist, @bounds;
            }
            # get variables specified by the coordinates attribute
            my @coordinates = get_variables_from_att($dump,$var,"coordinates");
            if (@coordinates) {
              print "   $var: coordinates = ".join(", ",@bounds)."\n" if $Opt{VERBOSE} > 1;
              push @vlist, @coordinates;
            }

            # formula terms
            # only add static terms to the vlist
            # time-varying terms are "external_variables"
            my @all_terms = get_formula_terms ($dump,$var);
            my (@vterms,@xterms);
            foreach (@all_terms) {
              if (is_static($dump,$_)) {
                push @vterms, $_;
              } else {
                push @xterms, $_;
              }
            }
            if (@vterms) {
              print "   $var: formula_terms(static) = ".join(", ",@vterms)."\n" if $Opt{VERBOSE} > 1;
              push @vlist, @vterms;
            }
            if (@xterms) {
              print "   $var: external_variables = ".join(", ",@xterms)."\n" if $Opt{VERBOSE} > 1;
              push @xlist, @xterms;
            }
         }

     #--- every time: need to extract time average/bounds variables ---
         if ($timename) {

            # fms-style time avg variables
            if (!$Opt{CMIP}) {
              my @time_avg_info = get_time_avg_info($dump,$var);
              if (@time_avg_info) {
                print "   $var: time_avg_info = ".join(", ",@time_avg_info)."\n" if $Opt{VERBOSE} > 1;
                push @vlist, @time_avg_info;
              }
            }

           # CF compliant time avg variables (bounds or climatology)
            my $SKIPVAR = 0;
            foreach my $tbound ( qw/bounds climatology/ ) {
               my $name = get_variable_att($dump,$timename,$tbound);
               if ($name) {
                  if ($name eq $var) {
                    # skip this variable because it is a time bounds variable
                    print "   $var: SKIPPING because it is a time bounds variable\n" if $Opt{VERBOSE} > 1;
                    $SKIPVAR = 1;
                    last;
                  }
                  print "   $var: $timename: $tbound = $name\n" if $Opt{VERBOSE} > 1;
                  push @vlist, $name;
               }
            }
            # skip this variable, reset the variable list
            if ($SKIPVAR) {
              #@vlist = @coords;
              #my @xlist;
              next;
            }
         }

         # if a single file is desired then do not extract into single-field files
         next if $Opt{onefile};

     #-----------------------
     #--- extract fields ----
     #-----------------------
         my $vlist = join ",",@vlist;
         print "   var=$var; timename=$timename; vlist=$vlist\n" if $Opt{VERBOSE} > 2;
         my $appendopt = "";
         $appendopt = "-A" if (-f "$cwd/.var.nc");
         print "ncks -h $appendopt -v $vlist $file .var.nc\n" if !$Opt{QUIET};
         next if $TEST;
         $ncstatus += command_retry("ncks","-h $appendopt -v $vlist $file .var.nc");

         # remove dimensions called "scalar_axis" (i.e., length = 1)
         # these are typically in the near-surface field files
         my @commands = remove_degenerate_dimension(".var.nc","scalar_axis");
         foreach my $cmd (@commands) {
           print "$cmd\n" if $Opt{VERBOSE} > 0;
           system($cmd);
         }

     #--- write to the output file ---
         if (-e "$odir/$var.nc") {
             die "ERROR: try to append to unlimited dimension when no unlimited dimension found" if !$timename;

            #--- append to existing file ---
            #push @commands, "dmget ".$Opt{odir}."/$var.nc";
             print  "mv $odir/$var.nc $odir/.tmp.nc\n" if $Opt{VERBOSE} > 2;
             system("mv $odir/$var.nc $odir/.tmp.nc");
             print "ncrcat $odir/.tmp.nc .var.nc $odir/$var.nc\n" if $Opt{VERBOSE} > 0;
             $ncstatus += command_retry("ncrcat","-h $odir/.tmp.nc .var.nc $odir/$var.nc");
             print  "rm -f $odir/.tmp.nc .var.nc" if $Opt{VERBOSE} > 2;
             system("rm -f $odir/.tmp.nc .var.nc");
         } else {
            #--- create new file ---
            #--- modify filename attribute ---
            #--- move to output directory ---

             my @ncatted_opts = set_ncatted_opts(".var.nc","$var.nc",$var);
             # external_variables attribute (time-varying formula terms)
             if (@xlist) {
               push @ncatted_opts, "-a external_variables,global,c,c,\"".join(" ",@xlist)."\"";
             }

             print  "mv .var.nc $odir/$var.nc\n" if $Opt{VERBOSE} > 2;
             system("mv .var.nc $odir/$var.nc");

             if (@ncatted_opts) {
               print  "ncatted -h -O @ncatted_opts $odir/$var.nc\n" if $Opt{VERBOSE} > 1;
               system("ncatted -h -O @ncatted_opts $odir/$var.nc");
             }
         }

         #@vlist = @coords;
         #my @xlist;
     }   ### end of variable loop ###

    #--- extract fields for single file option ---
     if ($Opt{onefile}) {
         my $vlist = join ",",@vlist;
         print "ncks -h -A -v $vlist $file ".$Opt{onefile}."\n" if $Opt{VERBOSE} > 0;
         $ncstatus += command_retry("ncks","-h -A -v $vlist $file ".$Opt{onefile});
         my @ncatted_opts = set_ncatted_opts($Opt{onefile},tailname($Opt{onefile}),"");
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
     unlink $file if $Opt{idir};

    # clean up
     unlink "ncks.out" if (-e "ncks.out");

 }   ### end of file loop ###

exit $ncstatus;

##################################################################
##################################################################
##################################################################
#  ----- help message -----

sub usage {
  my $name = substr($0,rindex($0,"/")+1);
  print "
Split NetCDF Variables

Usage:  $name [-d] [-l] [-i idir] [-o odir] [-f file] [-v vars]  files.....

        -l       = include readme log in output directory
        -V       = verbose option: increases standard output, multiple -V allowed
        -c       = cmip option, time average info variable not written
        -i idir  = input (archive) directory path
        -o odir  = output directory path
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
    push @attvars, $_ if ($dump =~ /\t\w+ $_\(.+\)/);
  }
  return @attvars;
}

#-------------------------------------------

sub get_variable_bounds {
  my $dump = shift;
  my $var = shift;
  my $timename = get_time_dimension($dump);
  my @bounds;
  # find dimensions for this variable
  if ($dump =~ /\t\w+ $var\((.+)\)/ ) {
    my @axes = split /, /, $1;
    foreach my $dim (@axes) {
      next if ($dim eq $timename);
      my $bnds = get_variable_att($dump,$dim,"bounds");
      push @bounds, $bnds if ($dump =~ /\t\w+ $bnds\(.+\)/);
      if (!$Opt{CMIP}) {
        my $bnds = get_variable_att($dump,$dim,"edges");
        push @bounds, $bnds if ($dump =~ /\t\w+ $bnds\(.+\)/);
      }
    }     
  }     
  return @bounds;
}

#-------------------------------------------

sub get_formula_terms {
  my $dump = shift;
  my $var = shift;
  my @terms;
  # find dimensions for this variable
  if ($dump =~ /\t\w+ $var\((.+)\)/ ) {
    my @axes = split /, /, $1;
    # check each axis for formula terms
    foreach my $dim (@axes) {
      my $form = get_variable_att($dump,$dim,"formula_terms");
     #print "   $dim:formula_terms = \"$form\"\n" if $form;
      foreach my $term ($form =~ /\w+:\s+(\w+)/g) {
        push @terms, $term if !grep{$_ eq $term} @terms;
      }
      # also check bounds for formula terms
      my $bnds = get_variable_att($dump,$dim,"bounds");
      my $form = get_variable_att($dump,$bnds,"formula_terms");
     #print "   $bnds:formula_terms = \"$form\"\n" if $form;
      foreach my $term ($form =~ /\w+:\s+(\w+)/g) {
        push @terms, $term if !grep{$_ eq $term} @terms;
      }
    }
  }
  return @terms;
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

sub command_retry {
  my $cmd = shift;
  my $opts = shift;
  my $fullcmd = $cmd." ".$opts;
  system($fullcmd);
  my $status = $?;
  if ($status) {
    print STDERR "ERROR: $cmd returned status $status\n";
    sleep 30;
    system ($fullcmd);
    $status = $?;
    if ($status) {
      print STDERR "ERROR ON RETRY: $cmd returned status $status\n";
    } else {
      print STDERR "RETRY SUCCESSFUL\n";
    }
  }
  return $status;
}

#-------------------------------------------
# cleanup attributes in the output file

sub set_ncatted_opts {
  my ($file,$filename,$var) = @_;
  my @opts;
  my $dump = `ncdump -h $file`;
  push @opts, "-a filename,global,m,c,\"$filename\"" if ($dump =~ /\t\t:filename = ".+" ;/);
  my $svar = $var;
  $svar = "\w+" if !$svar; # globally edit attributes
  if ($Opt{CMIP}) {
    push @opts, "-a time_avg_info,$var,d,," if ($dump =~ /\t\t$svar:time_avg_info = ".+" ;/);
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


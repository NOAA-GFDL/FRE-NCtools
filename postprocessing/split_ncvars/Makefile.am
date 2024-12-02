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
# Perform several substitutions
do_subst = $(SED) -e 's,[@]prefix[@],$(prefix),g' \
               -e 's,[@]pkglibexecdir[@],$(pkglibexecdir),g' \
               -e 's,[@]NCDUMP[@],$(NCDUMP),g' \
               -e 's,[@]NCRCAT[@],$(NCRCAT),g' \
               -e 's,[@]NCATTED[@],$(NCATTED),g' \
               -e 's,[@]NCKS[@],$(NCKS),g' \
               -e 's,[@]SED[@],$(SED),g' \
               -e 's,[@]EGREP[@],$(EGREP),g'


bin_SCRIPTS = split_ncvars.pl

split_ncvars.pl: split_ncvars.pl.in Makefile
	$(do_subst) < $(srcdir)/split_ncvars.pl.in > split_ncvars.pl
	chmod +x split_ncvars.pl

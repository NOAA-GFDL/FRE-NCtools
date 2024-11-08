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
bin_SCRIPTS = iceberg_comb.sh
EXTRA_DIST = iceberg_comb.sh.in \
    iceberg_comb.sh.txt \
    asciidoc.conf

# Documentation settings
A2X = a2x
ASCIIDOC_CFG = $(srcdir)/asciidoc.conf

do_subst = $(SED) -e 's,[@]prefix[@],$(prefix),g' \
               -e 's,[@]pkglibexecdir[@],$(pkglibexecdir),g' \
               -e 's,[@]NCDUMP[@],$(NCDUMP),g' \
               -e 's,[@]NCRCAT[@],$(NCRCAT),g' \
               -e 's,[@]NCATTED[@],$(NCATTED),g'

if WITH_ASCIIDOC
man_MANS = iceberg_comb.sh.1
else
man_MANS =
endif

iceberg_comb.sh: iceberg_comb.sh.in Makefile
	$(do_subst) < $(srcdir)/iceberg_comb.sh.in > iceberg_comb.sh
	chmod +x iceberg_comb.sh

.txt.1:
	$(A2X) -v -L\
          -D $(builddir) \
          -d manpage \
          -f manpage \
          -afre_version=$(PACKAGE_VERSION) \
          --asciidoc-opts="-f $(ASCIIDOC_CFG)" \
          $<

clean-local:
	$(RM) *.1

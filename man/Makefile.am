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
EXTRA_DIST = \
  combine-ncc.txt \
  combine_restarts.txt \
  decompress-ncc.txt \
  iceberg_comb.sh.txt \
  is-compressed.txt \
  list_ncvars.txt \
  list_ncvars.csh.txt \
  mppnccombine.txt \
  scatter-ncc.txt \
  asciidoc.conf

# Documentation settings
A2X = a2x
ASCIIDOC_CFG = $(srcdir)/asciidoc.conf

man_MANS =

if WITH_ASCIIDOC
man_MANS += \
  combine-ncc.1 \
  combine_restarts.1 \
  decompress-ncc.1 \
  iceberg_comb.sh.1 \
  is-compressed.1 \
  list_ncvars.1 \
  list_ncvars.csh.1 \
  mppnccombine.1 \
  scatter-ncc.1
endif

.txt.1:
	$(A2X) -v -L\
          -D $(builddir) \
          -d manpage \
          -f manpage \
          -afre_version=$(PACKAGE_VERSION) \
          --asciidoc-opts="-f $(ASCIIDOC_CFG)" \
          $<

clean-local:
	$(RM) *.1 *.xml

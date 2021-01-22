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
# License along with FRE-NCTools (LICENSE.md).  If not, see
# <http://www.gnu.org/licenses/>.
#***********************************************************************
dist_bin_SCRIPTS = combine_restarts
EXTRA_DIST = combine_restarts.txt asciidoc.conf

# Documentation settings
A2X = a2x
ASCIIDOC_CFG = $(srcdir)/asciidoc.conf

if WITH_ASCIIDOC
man_MANS = combine_restarts.1
else
man_MANS =
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
	$(RM) *.1

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
if WITH_MPI_TESTS
skip_MPI=""
else
skip_MPI="skip"
endif

if SKIP_XARRAY_TESTS
skipflag="skip"
else
skipflag=""
endif

TESTS_ENVIRONMENT = top_srcdir=$(abs_top_srcdir); export top_srcdir; \
                    export skip_mpi=$(skip_MPI); \
                    export skipflag=$(skipflag); \
                    PATH="$(abs_top_builddir)/src"'$(PATH_SEPARATOR)'"$$PATH";

TESTS = combine-ncc/combine \
        decompress-ncc/decompress \
        is-compressed/compressed \
        is-compressed/not-compressed \
        mppnccombine/six-digit-extension \
        mppnccombine/xt-decomp \
        mppnccombine/xyt-decomp \
        mppnccombine/yt-decomp \
        mppnccombine/list-output-files \
        scatter-ncc/scatter_1x1 \
        scatter-ncc/scatter_1x2 \
        scatter-ncc/scatter_2x1

EXTRA_DIST = \
  init.sh \
  $(TESTS)

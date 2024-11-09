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
                    export PKGLIBEXECDIR="$(abs_top_builddir)/src" \
                    PATH="$(abs_top_builddir)/src"'$(PATH_SEPARATOR)'"$$PATH";

TESTS = combine-ncc/combine \
        combine_restarts/combine \
        combine_restarts/dryrun \
        combine_restarts/extra_options \
        decompress-ncc/decompress \
        iceberg_comb/check-valid-file \
        iceberg_comb/combine \
        iceberg_comb/combine-no-icebergs \
        iceberg_comb/missing-file \
        iceberg_comb/outfile-exists \
        iceberg_comb/setnum-mismatch \
        is-compressed/compressed \
        is-compressed/not-compressed \
        list_ncvars/list_ncvars_csh \
        list_ncvars/list_ncvars_ftn \
        mppnccombine/dry-run \
        mppnccombine/list-output-files \
        mppnccombine/missing-file-fail \
        mppnccombine/missing-file-force \
        mppnccombine/non-zero-start-number \
        mppnccombine/output-kinds \
        mppnccombine/print-memory-usage \
        mppnccombine/remove-decomp-files \
        mppnccombine/six-digit-extension \
        mppnccombine/with-end-number \
        mppnccombine/xt-decomp \
        mppnccombine/xyt-decomp \
        mppnccombine/yt-decomp \
        scatter-ncc/scatter_1x1 \
        scatter-ncc/scatter_1x2 \
        scatter-ncc/scatter_2x1

EXTRA_DIST = \
  init.sh \
  $(TESTS)

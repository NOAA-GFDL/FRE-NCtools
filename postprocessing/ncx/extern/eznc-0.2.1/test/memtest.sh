#!/bin/bash
# Specify specific binary to test, or none will test all known.

if [[ $# != 0 ]]; then
TESTS=$1
else
TESTS="test_nc3 test_ezNc test_ezNcCompare test_ezNcCF test_ezNcBuffers test_ezNcCopySlabTask test_ezNcDateTime test_ezNcFiles test_ezNcSlab test_ezNcSlabIterator test_ezNcSlicer test_ezNcUtil test_ezNcVariant test_ezPointerVector"
fi

TESTFLAGS="--tmpdir /tmp/ --verbose --color"

for test in $TESTS; do
  echo $test
  valgrind --tool=memcheck --leak-check=full --show-reachable=yes --log-file=tmp.valgrind.log ./$test $TESTFLAGS
  grep "leak" tmp.valgrind.log || echo "ERROR: Leaks are possible. See tmp.$test.log" && mv tmp.valgrind.log tmp.$test.log
done
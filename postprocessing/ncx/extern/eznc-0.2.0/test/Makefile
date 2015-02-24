include ../Makefile.inc

TESTS = test_nc3 test_ezNc test_ezNcBuffers test_ezNcCF test_ezNcCompare test_ezNcCopySlabTask test_ezNcDateTime test_ezNcSlab test_ezNcSlabIterator test_ezNcSlicer test_ezNcUtil test_ezNcVariant test_ezPointerVector

TESTFLAGS = --tmpdir /tmp/ --verbose --color

all: $(TESTS)

test: $(TESTS)
	python test.py $(TESTFLAGS)

memtest: $(TESTS)
	./memtest.sh

test_ezNc: test_ezNc.cpp ../src/ezNc.hpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(EZTEST) $(EZSTR) $< -o $@

test_ezNcBuffers: test_ezNcBuffers.cpp ../src/ezNcBuffers.hpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(EZTEST) $< -o $@

test_ezNcCF: test_ezNcCF.cpp ../src/ezNcCF.hpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(EZTEST) $(EZSTR) $(EZSLICE)  $(EZODO) $< -o $@

test_ezNcCompare: test_ezNcCompare.cpp ../src/ezNcCompare.hpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(EZTEST) $(EZSTR) $(EZSLICE)  $(EZODO) $< -o $@

test_ezNcCopySlabTask: test_ezNcCopySlabTask.cpp ../src/ezNcCopySlabTask.hpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(EZTEST) $(EZSLICE) $(EZODO) $(EZSTR) $< -o $@

test_ezNcDateTime: test_ezNcDateTime.cpp ../src/ezNcDateTime.hpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(EZTEST) $(EZSLICE) $(EZSTR) $(EZDATE) $< -o $@
  
test_ezNcSlab: test_ezNcSlab.cpp ../src/ezNcSlab.hpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(EZTEST) $(EZSLICE) $(EZSTR) $< -o $@

test_ezNcSlabIterator: test_ezNcSlabIterator.cpp ../src/ezNcSlabIterator.hpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(EZTEST) $(EZSLICE) $(EZODO) $(EZSTR) $< -o $@

test_ezNcSlicer: test_ezNcSlicer.cpp ../src/ezNcSlicer.hpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(EZTEST) $(EZSLICE) $(EZSTR) $< -o $@

test_ezNcUtil: test_ezNcUtil.cpp ../src/ezNcUtil.hpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(EZTEST) $< -o $@

test_ezNcVariant: test_ezNcVariant.cpp ../src/ezNcVariant.hpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(EZTEST) $< -o $@

test_ezPointerVector: test_ezPointerVector.cpp ../src/ezPointerVector.hpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(EZTEST) $< -o $@

test_nc3: test_nc3.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(EZTEST) $< -o $@

clean:
	rm -f $(TESTS) tmp* *.o *.exe

include ../Makefile.inc

all: bench

bench: ncbench
	./ncbench 10 200 100 50 32 0 0 0 4 0 4 0 0 0 0 /tmp/tmp.nc

ncbench: ncbench.cpp
	$(CXX) $(CXXFLAGS) $(EZSTR) $(LDFLAGS) ncbench.cpp -o ncbench

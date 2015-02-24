EZDATE=-I/home/rsz/projects/ezDateTime/ezDateTime-0.0.0
EZODO=-I/home/rsz/projects/ezOdometer/ezOdometer-0.1.1
EZSLICE=-I/home/rsz/projects/ezSlice/ezSlice-0.0.0
EZSTR=-I/home/rsz/projects/ezStringUtil/ezStringUtil-0.1.0
EZTEST=-I/home/rsz/projects/eztest/ezTest-0.0.0

NC=/usr/local/netcdf-4.1.3/
NCCLASSIC=$(NC)/classic

PYTHONPATH:=$(PYTHONPATH):/home/rsz/projects/eztest/ezTest-0.0.0

CXXFLAGS=-O0 -I$(NC)/include -I../src
LDFLAGS=-L$(NC)/lib -lnetcdf -Wl,-rpath,$(NC)/lib

CLASSIC_CXXFLAGS=-I$(CLASSIC)/include -I../src
CLASSIC_LDFLAGS=-L$(CLASSIC)/lib -lnetcdf -Wl,-rpath,$(CLASSIC)/lib

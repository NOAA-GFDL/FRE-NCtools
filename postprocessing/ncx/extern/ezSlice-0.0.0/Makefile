.PHONY = all clean dist memtest test

EZTEST ?= /home/rsz/projects/eztest/eztest-0.0.0
TESTS = test_ezSlice

all: test

clean:
	rm -f *.o $(TESTS) tmp* *.exe

dist:
  ifndef VER
		@echo "ERROR: VER is not defined. Try: make dist VER=0.0.0"
  else
		cd ..; tar zcvf ezSlice-$(VER).tgz ezSlice-$(VER)
  endif

memtest: $(TESTS)
	@valgrind --tool=memcheck --leak-check=full --show-reachable=yes --log-file=tmp.valgrind.log ./test_ezSlice > /dev/null
	@grep "leaks" tmp.valgrind.log

test: $(TESTS)
	./test_ezSlice

test_ezSlice: test_ezSlice.cpp ezSlice.hpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I$(EZTEST) $< -o $@

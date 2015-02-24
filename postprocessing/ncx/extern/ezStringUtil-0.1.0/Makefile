.PHONY = all clean dist memtest test

EZTEST ?= /home/rsz/projects/eztest/eztest-0.0.0
TESTS = test_ezStringUtil

all: test

clean:
	rm -f *.o $(TESTS) tmp* *.exe

dist:
  ifndef VER
		@echo "ERROR: VER is not defined. Try: make dist VER=0.0.0"
  else
		cd ..; tar zcvf ezStringUtil-$(VER).tgz ezStringUtil-$(VER)
  endif

memtest: $(TESTS)
	@valgrind --tool=memcheck --leak-check=full --show-reachable=yes --log-file=tmp.valgrind.log ./test_ezStringUtil > /dev/null
	@grep "leak" tmp.valgrind.log

test: $(TESTS)
	./test_ezStringUtil --verbose --color

test_ezStringUtil: test_ezStringUtil.cpp ezStringUtil.hpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I$(EZTEST) $< -o $@
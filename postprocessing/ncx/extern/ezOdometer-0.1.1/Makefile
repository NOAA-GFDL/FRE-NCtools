.PHONY: clean test memtest

all: test

test: ezOdometer.hpp test.cpp
	$(CXX) $(CXXFLAGS) test.cpp -o test
	./test
	
memtest: test
	@valgrind --tool=memcheck --leak-check=full --show-reachable=yes --log-file=tmp.valgrind.log ./test > /dev/null
	@grep "leak" tmp.valgrind.log

clean:
	rm -f *.o test tmp* *.exe

dist:
  ifndef VER
		@echo "ERROR: VER is not defined. Try: make dist VER=0.0.0"
  else
		cd ..; tar zcvf ezOdometer-$(VER).tgz ezOdometer-$(VER)
  endif

FULLDEMO = python hybrid.py --clean 1 --tmpdir /tmp --verbose --color 1 2 3 4

all: ezTestDemo

ezTestDemo: ezTestDemo.cpp ../ezTest.hpp
	$(CXX) -I.. $< -o $@

test:
	$(FULLDEMO)

memtest:
	@valgrind --tool=memcheck --leak-check=full --show-reachable=yes --log-file=tmp.valgrind.log ./ezTestDemo --clean 1 --tmpdir /tmp --verbose --color
	@grep "leak" tmp.valgrind.log

clean:
	rm -f ezTestDemo tmp.* *.pyc


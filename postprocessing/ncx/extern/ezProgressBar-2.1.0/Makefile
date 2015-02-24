# Windows.
#CXX=cl /D WIN32

.PHONY: html test clean memtest

all: examples test memtest

clean:
	rm -f example etaExample rateExample test *.exe* *.tmp *.obj

html:
	pygmentize -O full,linenos=1,style=manni -o ezProgressBar.html ezProgressBar.hpp
	pygmentize -O full,linenos=1,style=manni -o ezETAProgressBar.html ezETAProgressBar.hpp
	pygmentize -O full,linenos=1,style=manni -o ezRateProgressBar.html ezRateProgressBar.hpp
	pygmentize -O full,linenos=1,style=manni -o example.html example.cpp
	pygmentize -O full,linenos=1,style=manni -o etaExample.html etaExample.cpp
	pygmentize -O full,linenos=1,style=manni -o rateExample.html rateExample.cpp

examples: example.cpp etaExample.cpp rateExample.cpp ezProgressBar.hpp ezETAProgressBar.hpp ezRateProgressBar.hpp
	g++ example.cpp -o example 
	g++ etaExample.cpp -o etaExample
	g++ rateExample.cpp -o rateExample
	./example
	./etaExample
	./rateExample
	
test: test.cpp
	$(CXX) test.cpp -o test;
	./test > test.tmp
	diff test.truth test.tmp
	
memtest:
	valgrind --tool=memcheck --leak-check=full --show-reachable=yes ./example
	valgrind --tool=memcheck --leak-check=full --show-reachable=yes ./etaExample
	valgrind --tool=memcheck --leak-check=full --show-reachable=yes ./rateExample

dist:
	cd ..; tar zcvf ezProgressBar-$(VER).tgz ezProgressBar-$(VER)

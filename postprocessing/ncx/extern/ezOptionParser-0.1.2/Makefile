.PHONY: all examples html memtest test clean

all: examples

examples:
	cd examples && $(MAKE) all

clean:
	cd examples && $(MAKE) clean

html:
	cd examples && $(MAKE) html
	pygmentize -O full,linenos=1,style=manni -o html/ezOptionParser.html ezOptionParser.hpp

memtest:
	cd examples && $(MAKE) memtest

test:
	cd examples && $(MAKE) test

dist:
  ifndef VER
		@echo "ERROR: VER is not defined. Try: make dist VER=0.0.0"
  else
		cd ..; tar zcvf ezOptionParser-$(VER).tgz ezOptionParser-$(VER)
  endif

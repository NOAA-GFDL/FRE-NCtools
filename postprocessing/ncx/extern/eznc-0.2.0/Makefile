.PHONY: all bench clean dist test memtest

include Makefile.inc

all: test

bench:
	cd bench && $(MAKE) bench

clean:
	cd test && $(MAKE) clean

dist:
  ifndef VER
		@echo "ERROR: VER is not defined. Try: make dist VER=0.0.0"
  else
		cd ..; tar zcvf eznc-$(VER).tgz eznc-$(VER)
  endif

test:
	cd test && $(MAKE) test

memtest:
	cd test && $(MAKE) memtest

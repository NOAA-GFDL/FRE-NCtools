
all:
	cd examples && $(MAKE) all

test:
	cd examples && $(MAKE) test

memtest:
	cd examples && $(MAKE) memtest

clean:
	rm -f ezTest tmp.* *.pyc
	cd examples && $(MAKE) clean

dist:
  ifndef VER
		@echo "ERROR: VER is not defined. Try: make dist VER=0.0.0"
  else
		cd ..; tar zcvf ezTest-$(VER).tgz ezTest-$(VER)
  endif

export NORNIR_PATH           = /usr/local
export NORNIR_PATH_LIB       = $(NORNIR_PATH)/lib
export NORNIR_PATH_INCLUDE   = $(NORNIR_PATH)/include/nornir

export CC                    = gcc
export CXX                   = g++ 
export OPTIMIZE_FLAGS        = -finline-functions -O3 -g #-DPOOL
export DEBUG_FLAGS           = #-DDEBUG_MANAGER -DDEBUG_PREDICTORS -DDEBUG_SELECTORS -DDEBUG_DF_STREAM -DDEBUG_NODE -DDEBUG_KNOB 
export CXXFLAGS              = -Wall -pedantic --std=c++11 $(OPTIMIZE_FLAGS) $(DEBUG_FLAGS)
export LDLIBS                = -lnornir -pthread -lrt -lm -lmlpack -llapack -lblas -lgsl -lgslcblas -larmadillo -lanl
export INCS                  = -I$(realpath ./src/external/fastflow) -I$(realpath ./src/external/tclap-1.2.1/include) -I/usr/include/libxml2
export LDFLAGS               = -L$(realpath .)/src

.PHONY: all demo clean cleanall install uninstall microbench bin test testquick gcov develcheck

all:
	$(MAKE) -C src
	$(MAKE) -C microbench
	$(MAKE) -C microbench checksupported
clean: 
	$(MAKE) -C src clean
	$(MAKE) -C demo clean
	$(MAKE) -C examples clean
	$(MAKE) -C microbench clean
	$(MAKE) -C bin clean
	$(MAKE) -C test clean
demo:
	$(MAKE) -C demo 
	$(MAKE) -C examples
bin:
	$(MAKE) -C bin
cppcheck:
	cppcheck --xml --xml-version=2 --enable=warning,performance,information,style --inline-suppr --error-exitcode=1 --force  . --suppressions-list=./test/cppcheck/suppressions-list.txt -isrc/external -itest -iexamples 2> cppcheck-report.xml || (cat cppcheck-report.xml; exit 2) 
# Compiles and runs all the tests.
test:
	$(MAKE) cleanall
	$(MAKE) "NORNIR_COVERAGE_FLAGS=-fprofile-arcs -ftest-coverage"
	$(MAKE) testquick
testquick:
	cd test && ./installdep.sh 
	cd ..
	$(MAKE) "NORNIR_COVERAGE_LIBS=-lgcov" -C test && cd test && ./runtests.sh
	cd ..
gcov:
	./test/gcov/gcov.sh
# Performs all the checks
develcheck:
	$(MAKE) cppcheck && $(MAKE) test && $(MAKE) gcov && $(MAKE) cleanall && $(MAKE) && $(MAKE) bin
cleanall:
	$(MAKE) -C src cleanall
	$(MAKE) -C demo cleanall
	$(MAKE) -C examples cleanall
	$(MAKE) -C microbench cleanall
	$(MAKE) -C bin cleanall
	$(MAKE) -C test cleanall
install:
	$(MAKE) -C src install
uninstall:
	$(MAKE) -C src uninstall
microbench:
	$(MAKE) -C microbench microbench

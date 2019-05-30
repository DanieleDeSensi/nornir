.PHONY = all

all: leo.o

leo.o: leo.cpp
	$(CXX) $(CXXFLAGS) --std=c++11 -DARMA_DONT_USE_WRAPPER -DARMA_USE_BLAS -DARMA_USE_LAPACK -DARMA_USE_HDF5 -c leo.cpp
clean: 
	-rm -fr *.o *~
cleanall:
	make clean

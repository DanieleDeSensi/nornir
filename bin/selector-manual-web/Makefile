TARGET = server

all: $(TARGET)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(INCS) $(LDFLAGS) $(LDLIBS) -L$(realpath ../../src/external/riff/src/external/nanomsg/build/) -lnanomsg -lanl -pthread

clean: 
	-rm -fr *.o *~
cleanall:
	-rm -fr *.o *~ 
	-rm -fr server
TARGET = interpreter.o stream.o mdfg.o mdfi.o

.PHONY: all clean cleanall

all: $(TARGET)

%.o: %.cpp *.hpp
	$(CXX) $(CXXFLAGS) $(NORNIR_COVERAGE_FLAGS) $(INCS) $(OPTIMIZE_FLAGS) -c -o $@ $<
clean: 
	-rm -fr *.o *~
cleanall:
	-rm -fr *.o *~ *.a


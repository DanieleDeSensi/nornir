TARGET = interface.o instrumenter.o manager-multi.o manager.o knob.o trigger.o node.o parameters.o selectors.o predictors.o explorers.o configuration.o stats.o 
TARGET_OMP = ompt.o instrumenter.o parameters.o
CXXFLAGS += -fPIC

.PHONY: all clean cleanall install uninstall lib

all:
	$(MAKE) -C external
	$(MAKE) -C dataflow
	$(MAKE) lib
ompt.o: ompt.cpp 
	$(CXX) $(CXXFLAGS) $(NORNIR_COVERAGE_FLAGS) $(INCS) $(OPTIMIZE_FLAGS) -c -o ompt.o ompt.cpp
%.o: %.cpp *.hpp
	$(CXX) $(CXXFLAGS) $(NORNIR_COVERAGE_FLAGS) $(INCS) $(OPTIMIZE_FLAGS) -c -o $@ $<
clean: 
	-rm -fr *.o *~
	$(MAKE) -C dataflow clean
	$(MAKE) -C external clean
cleanall:
	-rm -fr *.o *~ *.a
	$(MAKE) -C dataflow cleanall
#	$(MAKE) -C external cleanall
lib: $(TARGET) $(TARGET_OMP)
	rm -rf libnornir.a
	ar rc libnornir.a $(TARGET) dataflow/*.o
	ranlib libnornir.a
	ar -M < buildlib.ar
	$(CXX) -shared -fPIC $(TARGET) dataflow/*.o external/leo/leo.o external/queues/hzdptr.o external/queues/xxhash.o external/graphite-c-client/graphite-client.o -Wl,--whole-archive external/mammut/mammut/libmammut.a external/riff/src/libriff.a -Wl,--no-whole-archive -o libnornir.so
	$(CXX) -shared -fPIC $(TARGET_OMP) -Wl,--whole-archive external/mammut/mammut/libmammut.a external/riff/src/libriff.a -Wl,--no-whole-archive -o libnornir_ompt.so
install:
	mkdir -p $(NORNIR_PATH_INCLUDE)
	mkdir -p $(NORNIR_PATH_LIB)
	find . -name '*.h' -exec cp --parents \{\} $(NORNIR_PATH_INCLUDE) \;
	find . -name '*.hpp' -exec cp --parents \{\} $(NORNIR_PATH_INCLUDE) \;
	find . -name '*.tpp' -exec cp --parents \{\} $(NORNIR_PATH_INCLUDE) \;
	cp ./libnornir.a $(NORNIR_PATH_LIB)/libnornir.a
	cp ./libnornir.so $(NORNIR_PATH_LIB)/libnornir.so
	cp ./libnornir_ompt.so $(NORNIR_PATH_LIB)/libnornir_ompt.so
uninstall:
	rm -r $(NORNIR_PATH_INCLUDE)/*
	rm $(NORNIR_PATH_LIB)/libnornir.a
	rm $(NORNIR_PATH_LIB)/libnornir.so
	rm $(NORNIR_PATH_LIB)/libnornir_ompt.so


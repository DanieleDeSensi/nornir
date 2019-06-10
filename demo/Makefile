TARGET               = pipe_ffcompat application_instrumented farm farm_ffcompat farm_ffcompat_new_api farm_accelerator_1 farm_accelerator_2 farm_accelerator_3 farm_accelerator_4 parallel_for
INCS += -I$(realpath ../src)

.PHONY: all clean cleanall

all: $(TARGET)
	$(MAKE) -C dataflow

%: %.cpp ../src/libnornir.a
	$(CXX) $(CXXFLAGS) -o $@ $< $(INCS) $(LDFLAGS) $(LDLIBS)
clean: 
	-rm -fr *.o *~
	$(MAKE) -C dataflow clean
cleanall:
	-rm -fr *.o *~ 
	-rm -fr $(TARGET)
	$(MAKE) -C dataflow cleanall


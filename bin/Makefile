TARGET               = manager-single manager-external manager-blackbox selector-manual manager-single-bb

.PHONY: all clean cleanall

all: $(TARGET)
	$(MAKE) -C selector-manual-web

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(INCS) $(LDFLAGS) $(LDLIBS) 
clean: 
	-rm -fr *.o *~
cleanall:
	-rm -fr *.o *~ 
	-rm -fr $(TARGET)

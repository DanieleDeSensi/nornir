INCS                 += -I$(realpath ../src/external/fastflow)
TARGET               = ticksPerNs check

.PHONY: all clean cleanall microbench

all: $(TARGET)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(INCS) $(LDFLAGS) $(LDLIBS) 
clean: 
	-rm -fr *.o *~
cleanall:
	-rm -fr *.o *~ 
	-rm -fr $(TARGET)
checksupported:
	@echo "=============================== ATTENTION ==============================="
	@echo "| Sudo rights are needed to check if all the mechanisms present on this |"
	@echo "| machine are supported by nornir.                                      |"
	@echo "========================================================================="
	sudo ./check
microbench:
	cp ../src/external/mammut/demo/energy/voltageTable.cpp .
	$(CXX) $(CXXFLAGS) -o voltageTable voltageTable.cpp $(INCS) -I../src/external/mammut/  -L../src/external/mammut/mammut $(LDFLAGS) $(LDLIBS)
	$(CXX) $(CXXFLAGS) -o idlePower idlePower.cpp $(INCS) -I../src/external/mammut/  -L../src/external/mammut/mammut $(LDFLAGS) $(LDLIBS)
	@echo ""
	@echo ""
	@echo "=============================== ATTENTION ==============================="
	@echo "| Sudo rights are needed to compute the voltage table. If voltage table |" 
	@echo "| is not computed, it will not be possible to predict power consumption |"
	@echo "| of the applications.                                                  |" 
	@echo "========================================================================="
	sudo sh runmicrobenchs.sh


/*
 * selector-external.cpp
 *
 * This tool is needed to manually control the reconfiguration when
 * STRATEGY_SELECTION_MANUAL_CLI is chosen as selector.
 *
 * Created on: 30/08/2017
 *
 * =========================================================================
 *  Copyright (C) 2015-, Daniele De Sensi (d.desensi.software@gmail.com)
 *
 *  This file is part of nornir.
 *
 *  nornir is free software: you can redistribute it and/or
 *  modify it under the terms of the Lesser GNU General Public
 *  License as published by the Free Software Foundation, either
 *  version 3 of the License, or (at your option) any later version.

 *  nornir is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  Lesser GNU General Public License for more details.
 *
 *  You should have received a copy of the Lesser GNU General Public
 *  License along with nornir.
 *  If not, see <http://www.gnu.org/licenses/>.
 *
 * =========================================================================
 */

#include "../src/nornir.hpp"
#include <tclap/CmdLine.h>

using namespace nornir;

int main(int argc, char** argv){
    TCLAP::CmdLine cmd("Manually controls Nornir execution by setting individually the knobs.", ' ', "1.0");
    // Flag Name Description Required DefaultValue TypeDesc 
    TCLAP::ValueArg<double> virtualCoresArg("v", "virtualcores", "Number of VirtualCores (Relative value between 0 and 100)", false, -1, "double", cmd);
    TCLAP::ValueArg<double> frequencyArg("f", "frequency", "Frequency (Relative value between 0 and 100)", false, -1, "double", cmd);
    TCLAP::ValueArg<double> clkmodArg("e", "clkmodemu", "Clock Modulation Simulation (Relative value between 0 and 100)", false, -1, "double", cmd);
    //TODO Add option to specify real values instead of relative ones.
    cmd.parse(argc, argv);
    double virtualCores = virtualCoresArg.getValue();
    double frequency = frequencyArg.getValue();
    double clkmod = clkmodArg.getValue();

    bool needChange = false;
    // Set knobs
    KnobsValues kv(KNOB_VALUE_RELATIVE);
 	// First load the old ones (we need to do it so to do not overwrite unchanged values)
	std::ifstream instream(getSelectorManualCliControlFile());
	instream >> kv;
	instream.close();

    if(virtualCores != -1){
    	if(virtualCores < 0 || virtualCores > 100){
    		std::cerr << "FATAL ERROR: virtualcores must be between 0 and 100" << std::endl;
    		return -1;
    	}
        kv[KNOB_VIRTUAL_CORES] = virtualCores;
        needChange = true;
    }
    if(frequency != -1){
		if(frequency < 0 || frequency > 100){
    		std::cerr << "FATAL ERROR: frequency must be between 0 and 100" << std::endl;
    		return -1;
    	}
    	kv[KNOB_FREQUENCY] = frequency;	
        needChange = true;
    }
    if(clkmod != -1){
        if(clkmod < 0 || clkmod > 100){
            std::cerr << "FATAL ERROR: clkmodemu must be between 0 and 100" << std::endl;
            return -1;
        }
        kv[KNOB_CLKMOD] = clkmod; 
        needChange = true;
    }
    
    if(needChange){   	
	    // Write knobs to file.
	    std::ofstream outstream;
	    outstream.open(getSelectorManualCliControlFile());
	    if(!outstream.is_open()){
	    	std::cerr << "FATAL ERROR: impossible to open " << getSelectorManualCliControlFile() << std::endl;
	    	return -1;
	    }
	    outstream << kv;
	    outstream.close();
	}
	return 0;
}
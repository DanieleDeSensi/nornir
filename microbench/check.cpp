/*
 * check.cpp
 *
 * Created on: 26/03/2016
 *
 * Checks if the architecture is supported.
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

#include <iostream>
#include "../src/external/mammut/mammut/mammut.hpp"

using namespace std;
using namespace mammut;
using namespace mammut::topology;
using namespace mammut::energy;

#ifndef __linux__
#error "For the moment, Nornir only supports linux machines."
#endif

int main(int argc, char** argv){
    Mammut m;
    Topology* t = m.getInstanceTopology();
    vector<VirtualCore*> virtualCores = t->getVirtualCores();
    for(size_t i = 0; i < virtualCores.size(); i++){
        if(!virtualCores.at(i)->hasFlag("constant_tsc")){
            cerr << "=======================ATTENTION=======================\n"
                    "| Processor must have the constant_tsc flag. This is  |\n"
                    "| needed since nornir gets timestamps considering     |\n "
                    "| cores clock ticks. This constraint is not imposed   |\n"
                    "| by the algorithm but it is needed since for the     |\n"
                    "| moment no other timestamping mechanisms are         |\n"
                    "|supported.                                           |\n"
                    "=======================================================\n";
            return -1;
        }
    }


    Energy* energy = m.getInstanceEnergy();
    Counter* counter = energy->getCounter();
    if(!counter){
        cerr << "========================WARNING========================\n"
                "| Power counters not available on this machine (or if |\n"
                "| available, they are still not supported by Nornir). |\n"
                "| Accordingly, no guarantees on power consumption can |\n"
                "| be provided.                                        |\n"
                "=======================================================\n";
    }
}

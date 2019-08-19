/*
 * manager-single.cpp
 *
 * Monitors and adapt an external application by using the riff library.
 * Needs the explicit channel name. This is deprecated. Use manager-external.
 *
 * Created on: 21/06/2016
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

#include <vector>
#include <iostream>
#include <fstream>
#include <nornir/nornir.hpp>

using namespace nornir;

int main(int argc, char * argv[]) {
    char* channelName = NULL;
    if(argc != 2) {
        std::cerr << "use: " 
                  << argv[0] 
                  << " channelName\n";
        return -1;
    }   
    channelName = argv[1];

    Parameters p("parameters.xml");
    ManagerInstrumented m(channelName, p);
    m.start();
    m.join();

    return 0;
}


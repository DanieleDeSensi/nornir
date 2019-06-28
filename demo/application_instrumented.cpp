/*
 * application_instrumented.cpp
 *
 * Created on: 05/07/2016
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

#include <nornir/nornir.hpp>

#include <unistd.h>

/**
 * Basic example on how to monitor an instrumented iterative application.
 * Original code was:
 * --------------------------------
 * int main(int argc, char** argv){
 *   size_t i = 0;
 *   while(i++ < 10){
 *     sleep(2);
 *   }
 *   return 0;
 * }
 * --------------------------------
 * Instrumented code is:
 */
int main(int argc, char** argv){
    nornir::Instrumenter instr("parameters.xml", 1, NULL, true);
    size_t i = 0;
    while(i++ < 10){
        instr.begin();
        sleep(2);
        instr.end();
    }
    instr.terminate();
    return 0;
}


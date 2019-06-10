/*
 * ticksPerNs.cpp
 *
 * Created on: 28/07/2015
 *
 * Tries to determine how many ticks are in a nanosecond.
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
#include <sstream>
#include <cmath>
#include <sched.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <nornir/external/fastflow/ff/utils.hpp>

unsigned long getNanoSeconds(){
    struct timespec spec;
    clock_gettime(CLOCK_MONOTONIC, &spec);
    return spec.tv_sec * 1000000000 + spec.tv_nsec;
}

double getTicksPerNanosec(){
    double x = 0.691812048120;

    unsigned long start = getNanoSeconds();
    ticks t1 = getticks();
    while(getticks() - t1 < 1000000){
        x = std::sin(x);
    }
    ticks t2 = getticks();
    unsigned long end = getNanoSeconds();

    std::ostringstream sstream;
    sstream << x;
    FILE *p = popen(std::string("echo " + sstream.str() +
                                " >/dev/null").c_str(), "r");
    pclose(p);

    return ((double) (t2 - t1) / (double) (end - start));
}

int main(int argc, char** argv){
#ifndef __linux__
#error "This only works on Linux"
#endif
    // Pin to a core because TSC may not be coherent between different cores.
    cpu_set_t mask;
    CPU_ZERO(&mask);
    CPU_SET(0, &mask);
    assert(!sched_setaffinity(0, sizeof(mask), &mask));
    setpriority(PRIO_PROCESS, 0, -20);

    uint numIterations = 1000;
    double x = 0;
    for(size_t i = 0; i < numIterations; i++){
        x += getTicksPerNanosec();
    }
    std::cout << "<ticksPerNs>"  << x / (double) numIterations
              << "</ticksPerNs>" << std::endl;
}


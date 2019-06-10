/*
 * farm_accelerator.hpp
 *
 * Created on: 27/02/2016
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

/**
 * Basic test for the nornir accelerator.
 */

#include <vector>
#include <iostream>
#include <fstream>
#include "../src/nornir.hpp"

using namespace ff;

#define MICROSECSSLEEP 1000000

static int maxTasks;

class Scheduler: public nornir::Scheduler<double, int>{
public:
    int* schedule(double* task){
        std::cout << "Scheduler received " << *task
                  << " and send: " << (int)(*task) << std::endl;
        int* newTask = new int((int)*task);
        delete task;
        return newTask;
    }
};

/**
 * Worker.
 */
class Worker: public nornir::Worker<int, int>{
public:
    int* compute(int * task) {
        usleep(MICROSECSSLEEP);
        std::cout << "Worker " << getId()
                  << " received task " << *task << std::endl;
        return task;
    }
};

class Gatherer: public nornir::Gatherer<int>{
public:
    void gather(int* task){
        std::cout << "Gatherer received task " << *task << std::endl;
        delete task;
    }
};

int main(int argc, char * argv[]) {
    int nworkers = 1;
    int streamlen = 10;

    if (argc>1) {
        if (argc!=3 || atoi(argv[1]) <= 0 || atoi(argv[2]) <= 0) {
            std::cerr << "use: "
                      << argv[0]
                      << " nworkers(>0) streamlen(>0)\n";
            return -1;
        }
        nworkers = atoi(argv[1]);
        streamlen = atoi(argv[2]);
    }

    if (!nworkers || !streamlen) {
        std::cerr << "Wrong parameters values\n";
        return -1;
    }

    maxTasks = streamlen;

    nornir::Parameters p("parameters.xml");
    nornir::FarmAccelerator<double, int, int> farm(&p);
    farm.addScheduler(new Scheduler);
    for(size_t i = 0; i < (size_t) nworkers; i++){
        farm.addWorker(new Worker);
    }
    farm.addGatherer(new Gatherer);
    farm.start();

    for(size_t i = 0; i < (size_t) streamlen; i++){
        double *task = new double(i + i * 0.1);
        farm.offload(task);
        std::cout << "Sent " << *task << std::endl;
    }

    farm.shutdown();
    farm.wait();
    std::cout << "Accelerator terminated." << std::endl;
    return 0;
}


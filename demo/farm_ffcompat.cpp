/*
 * farm_ffcompat.hpp
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
 * This is a demo on how to enable an existing FastFlow farm to operate
 * under the nornir management.
 *
 * The steps to be followed are:
 *  1. Emitter, Workers and Collector of the farm must extend
 *     nornir::AdaptiveNode instead of ff::ff_node.
 *     ATTENTION: svc_init and svc_end are now called after each rethreading.
 *     Accordingly, if those operations need to be performed only once, you
 *     should ensure that.
 *  2. In the emitter, when there are no more tasks to be produced, please
 *     just call the macro TERMINATE_APPLICATION instead of returing NULL.
 *  3. If the application wants to be aware of the changes in the number
 *     of workers, the nodes can implement the notifyRethreading virtual
 *     method.
 *  4. Pass the existing farm to the manager.
 */

#include <vector>
#include <iostream>
#include <fstream>
#include <nornir/nornir.hpp>

using namespace ff;

#define MICROSECSSLEEP 1000000

/**
 * The emitter of the farm.
 */
class Emitter: public nornir::AdaptiveNode {
public:
    explicit Emitter(int max_task):ntask(max_task){;}

    void * svc(void *) {
        usleep(MICROSECSSLEEP);
        int * task = new int(ntask);
        --ntask;
        if (ntask<0){
            std::cout << "Emitter finished" << std::endl;
            //cppcheck-suppress memleak
            TERMINATE_APPLICATION;
        }
        //cppcheck-suppress memleak
        return task;
    }
private:
    int ntask;
};

/**
 * Farm worker.
 */
class Worker: public nornir::AdaptiveNode{
public:
    int svc_init(){
        std::cout << "Worker svc_init called" << std::endl;
        return 0;
    }

    void notifyRethreading(size_t oldNumWorkers, size_t newNumWorkers){
        /**
         * Do operation to consolidate state after the number of workers
         * changed from oldNumWorkers to newNumWorkers.
         * This operation is executed while the farm is frozen so it's safe
         * to change the internal state.
         */
    }
    void * svc(void * task) {
        int * t = (int *)task;
        usleep(MICROSECSSLEEP);
        std::cout << "Worker " << ff_node::get_my_id()
                  << " received task " << *t << "\n";
        return task;
    }
};

/**
 * The collector of the farm.
 */
class Collector: public nornir::AdaptiveNode {
public:
    void * svc(void * task) {
        int * t = (int *)task;
        if (*t == -1) return NULL;
        return task;
    }
};

int main(int argc, char * argv[]) {
    int nworkers = 1;
    int streamlen = 10;

    if (argc>1) {
        if (argc!=3) {
            std::cerr << "use: "
                      << argv[0]
                      << " nworkers streamlen\n";
            return -1;
        }
        nworkers=atoi(argv[1]);
        streamlen=atoi(argv[2]);
    }

    if (!nworkers || !streamlen) {
        std::cerr << "Wrong parameters values\n";
        return -1;
    }

    ff::ff_farm<> farm; // farm object

    Emitter E(streamlen);
    farm.add_emitter(&E);

    std::vector<ff_node *> w;
    for(int i=0;i<nworkers;++i) w.push_back(new Worker);
    farm.add_workers(w); // add all workers to the farm

    Collector C;
    farm.add_collector(&C);

    /***************************************************************/
    /*  START - New code needed with respect to the existing code. */
    /***************************************************************/
    nornir::Parameters ap("parameters.xml"); // Load parameters.
    nornir::ManagerFastFlow amf(&farm, ap); // Create nornir manager.
    amf.start(); // Start farm.
    amf.join(); // Wait for farm end.
    /***************************************************************/
    /*  END - New code needed with respect to the existing code. */
    /***************************************************************/

    std::cout << "Farm end" << std::endl;
    std::cerr << "DONE, time= " << farm.ffTime() << " (ms)\n";
    farm.ffStats(std::cerr);

    return 0;
}


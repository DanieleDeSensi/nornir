/*
 * pipe_ffcompat.hpp
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
 * This is a demo on how to enable an existing FastFlow pipeline to operate
 * under the nornir management.
 * Currently we only support pipelines where stages are either sequential nodes
 *
 * The steps to be followed are:
 *  1. Emitter, Workers and Collector of the farms must extend
 *     nornir::AdaptiveNode instead of ff::ff_node.
 *     ATTENTION: svc_init and svc_end are now called after each rethreading.
 *     Accordingly, if those operations need to be performed only once, you
 *     should ensure that.
 *  2. If the application wants to be aware of the changes in the number
 *     of workers, the nodes can implement the notifyRethreading virtual
 *     method.
 *  3. Pass the existing pipeline to the manager.
 */

#include <vector>
#include <iostream>
#include <fstream>
#include "../src/nornir.hpp"

using namespace ff;

int dummyTask;

class Generator: public nornir::AdaptiveNode {
public:
    void * svc(void *) {
        return (void*) &dummyTask;
    }
};

/**
 * Farm worker.
 */
class Worker: public nornir::AdaptiveNode{
private:
    unsigned long _sleepTimeMs;
    uint _farmId;
    uint _workerId;
public:
    Worker(unsigned long sleepTimeMs, uint farmId, uint workerId):
        _sleepTimeMs(sleepTimeMs), _farmId(farmId), _workerId(workerId){;}

    int svc_init(){
        //std::cout << "Worker svc_init called" << std::endl;
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
        //std::cout << "Executed <" << _farmId << ", " << _workerId << ">" << std::endl;
        usleep(_sleepTimeMs * MAMMUT_MICROSECS_IN_MILLISEC);
        return task;
    }
};

class DummyNode: public nornir::AdaptiveNode {
public:
    void* svc(void* task){
        return task;
    }
};

int main(int argc, char * argv[]) {
    int nworkers = 4;
    ff::ff_farm<> farm1, farm2, farm3; // farm objects
    ff::ff_pipeline pipe;

    std::vector<ff_node *> w1;
    for(int i = 0; i < nworkers; i++) 
        w1.push_back(new Worker(100, 1, i));
    farm1.add_workers(w1); // add all workers to the farm
    farm1.add_emitter(new DummyNode());
    farm1.add_collector(new DummyNode());

    std::vector<ff_node *> w2;
    for(int i = 0; i < nworkers; i++) 
        w2.push_back(new Worker(200, 2, i));
    farm2.add_workers(w2); // add all workers to the farm
    farm2.add_emitter(new DummyNode());
    farm2.add_collector(new DummyNode());

    std::vector<ff_node *> w3;
    for(int i = 0; i < nworkers; i++) 
        w3.push_back(new Worker(400, 3, i));
    farm3.add_workers(w3); // add all workers to the farm
    farm3.add_emitter(new DummyNode());
    farm3.add_collector(new DummyNode());

    pipe.add_stage(new Generator());
    pipe.add_stage(&farm1);
    pipe.add_stage(&farm2);
    pipe.add_stage(&farm3);

    /***************************************************************/
    /*  START - New code needed with respect to the existing code. */
    /***************************************************************/
    nornir::Parameters ap("parameters.xml"); // Load parameters.
    std::vector<bool> farmFlags(4, true);
    farmFlags[0] = false;
    nornir::ManagerFastFlowPipeline amf(&pipe, farmFlags, ap); // Create nornir manager.
    amf.start(); // Start farm.
    amf.join(); // Wait for farm end.
    /***************************************************************/
    /*  END - New code needed with respect to the existing code. */
    /***************************************************************/
    pipe.ffStats(std::cerr);
    return 0;
}


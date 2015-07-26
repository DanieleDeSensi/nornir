/* -*- Mode: C++; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ***************************************************************************
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License version 2 as 
 *  published by the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *  As a special exception, you may use this file as part of a free software
 *  library without restriction.  Specifically, if other files instantiate
 *  templates or use macros or inline functions from this file, or you compile
 *  this file and link it with other files to produce an executable, this
 *  file does not by itself cause the resulting executable to be covered by
 *  the GNU General Public License.  This exception does not however
 *  invalidate any other reasons why the executable file might be covered by
 *  the GNU General Public License.
 *
 ****************************************************************************
 */

/**

 Very basic test for the FastFlow farm.
 
*/
#include <vector>
#include <iostream>
#include <fstream>
#include "../farm.hpp"
#include "../predictors_impl.hpp"

using namespace ff;

// generic worker
class Worker: public adpff::adp_ff_node{
public:
    int adp_svc_init(){
        std::cout << "Worker svc_init called" << std::endl;
        return 0;
    }

    void * adp_svc(void * task) {
        int * t = (int *)task;
        sleep(1);
        std::cout << "Worker " << ff_node::get_my_id()
                  << " received task " << *t << "\n";
        return task;
    }
    // I don't need the following functions for this test
    //int   svc_init() { return 0; }
    //void  svc_end() {}

};

// the gatherer filter
class Collector: public adpff::adp_ff_node {
public:
    void * adp_svc(void * task) {
        int * t = (int *)task;
        if (*t == -1) return NULL;
        return task;
    }
};

// the load-balancer filter
class Emitter: public adpff::adp_ff_node {
public:
    Emitter(int max_task):ntask(max_task) {};

    void * adp_svc(void *) {
        sleep(1);
        int * task = new int(ntask);
        --ntask;
        if (ntask<0){
            std::cout << "Emitter finished" << std::endl;
            return NULL;
        }
        return task;
    }
private:
    int ntask;
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
    
    adpff::Observer obs;
    adpff::AdaptivityParameters ap("demo-farm.xml");
    ap.observer = &obs;
    adpff::adp_ff_farm<> farm(ap); // farm object
    
    Emitter E(streamlen);
    farm.add_emitter(&E);

    std::vector<ff_node *> w;
    for(int i=0;i<nworkers;++i) w.push_back(new Worker);
    farm.add_workers(w); // add all workers to the farm

    Collector C;
    farm.add_collector(&C);
    
    if (farm.run_then_freeze()<0) {
    //if (farm.run_and_wait_end()<0) {
        error("running farm\n");
        return -1;
    }
    farm.wait();
    std::cout << "Farm end" << std::endl;
    std::cerr << "DONE, time= " << farm.ffTime() << " (ms)\n";
    farm.ffStats(std::cerr);

    return 0;
}


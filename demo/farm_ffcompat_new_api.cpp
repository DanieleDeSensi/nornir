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
 * under the nornir management. This demo is specific for the new typed
 * FastFlow interface.
 *
 * The steps to be followed are:
 *  1. If you defined any node extending ff_node_t, you need to extend nrnr_node_t instead
 *  2. Pass the existing farm to the manager.
 */

#include <vector>
#include <iostream>
#include <fstream>
#include "../src/nornir.hpp"

using namespace ff;

int main(int argc, char** argv){
   const int K = 10;
   int k = 0;
   auto lambda = [K, &k]() -> long* {
       if (k++ == K){
           printf("Going to terminate");
           return NULL;
       }
       sleep(1);
       return new long(k);
   };
   auto myF = [](long *t, ff_node* const node)->long* {
       printf("hello I've got one task my id is=%ld\n", node->get_my_id());
       delete t;
       return t;
   };
   struct Emitter: public nornir::nrnr_node_t<long> {
       std::function<long*()> F;
       Emitter(std::function<long*()> F):F(F) {}
       long *svc(long*) { return F(); }
   };
   ff_Farm<long> farm(myF, 4);

   farm.remove_collector();
   Emitter E(lambda);
   farm.add_emitter(E);
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
}

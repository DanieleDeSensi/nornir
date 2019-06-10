/*
 * mdfi.cpp
 *
 * Created on: 26/03/2016
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

#include <nornir/dataflow/mdfi.hpp>

namespace nornir{
namespace dataflow{

void Mdfi::compute(){
    Data d;
    auto sm = sourcesMap.begin();
    while(sm != sourcesMap.end()){
        d.setSource(tInput[sm->second].task, sm->first);
        sm++;
    }

    auto dm = destinationsMap.begin();
    while(dm != destinationsMap.end()){
        d.setDestination(&(tOutput[dm->second].result), dm->first);
        dm++;
    }

    comp->compute(&d);
}

void Mdfi::reset(ulong newId){
    setGid(newId);
    for(uint i = 0; i < dInput; i++){
        tInput[i].clear();
    }

    for(uint i = 0; i < dOutput; i++){
        TokenId ti = tOutput[i].getDest();
        ti.setGraphId(newId);
        tOutput[i] = OutputToken(NULL, ti);
    }
}

}
}



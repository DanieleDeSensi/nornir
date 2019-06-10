/*
 * mdfg.cpp
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

#include <nornir/dataflow/mdfg.hpp>

namespace nornir{
namespace dataflow{

Mdfg::Mdfg():_nextId(0), _id(0), _firstId(std::numeric_limits<ulong>::max()),
             _lastId(std::numeric_limits<ulong>::max()), _init(false){
    ;
}

Mdfg::Mdfg(Computable* c):_nextId(1), _id(0), _firstId(std::numeric_limits<ulong>::max()),
                          _lastId(std::numeric_limits<ulong>::max()), _init(false){
    _instructions.emplace_back(c, 0);
}

Mdfg::Mdfg(const Mdfg& g, ulong gid):_nextId(g._nextId),_id(gid),
        _firstId(g._firstId),
        _lastId(g._lastId), _init(false){
    _instructions.reserve(g._instructions.size());
    for(size_t i = 0; i < g._instructions.size(); i++){
        _instructions.emplace_back(g._instructions[i]);
        _instructions.back().reset(gid);
    }
}

Mdfg::~Mdfg(){
    ;
}

void Mdfg::reset(ulong newId){
    _id = newId;
    for(size_t i = 0; i < _instructions.size(); i++){
        _instructions[i].reset(newId);
    }
}

}
}


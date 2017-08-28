/*
 * instrumenter.hpp
 *
 * Created on: 03/07/2016
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

#ifndef NORNIR_INSTRUMENTER_HPP_
#define NORNIR_INSTRUMENTER_HPP_

#include "external/knarr/src/knarr.hpp"
#include "external/knarr/src/external/cppnanomsg/nn.hpp"
#include "external/knarr/src/external/nanomsg/src/pair.h"
#include "external/mammut/mammut/mammut.hpp"

#define INSTRUMENTATION_CONNECTION_CHANNEL "ipc:///tmp/nornir.ipc"

namespace nornir{

class InstrumenterHelper: public knarr::Application, mammut::utils::NonCopyable{
public:
    InstrumenterHelper(std::pair<nn::socket*, uint> p,
                       size_t numThreads = 1,
                       knarr::Aggregator* aggregator = NULL):
        knarr::Application(*p.first, p.second, numThreads, aggregator){;}
};

/**
 * Extends knarr::Application by providing a way to automatically connect
 * to the nornir manager server.
 */
class Instrumenter: public InstrumenterHelper{
private:
    std::pair<nn::socket*, uint> getChannel(const std::string& parametersFile) const;
public:
    /**
     * Creates a client for interaction with a local server.
     * @param parameters The file containing the Nornir parameters.
     */
    explicit Instrumenter(const std::string& parametersFile,
                          size_t numThreads = 1,
                          knarr::Aggregator* aggregator = NULL);
};

}

#endif /* NORNIR_INSTRUMENTER_HPP_ */

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

#include "utils.hpp"
#include <riff/riff.hpp>
#include <riff/external/cppnanomsg/nn.hpp>
#include <riff/external/nanomsg/src/pair.h>
#include <mammut/mammut.hpp>

namespace nornir{

// We cannot use XDG runtime directories since they are
// user specific. Manager and application could be run
// by different users (e.g. manager by sudo user
// and application by normal user).
inline std::string getInstrumentationChannelsPath(){
    return "/tmp/nornir_external_manager/";
}

inline std::string getInstrumentationConnectionChannelPath(){
    return getInstrumentationChannelsPath() + "/connection.ipc";
}

inline std::string getInstrumentationPidChannelPath(uint pid){
    return getInstrumentationChannelsPath() + "/" + mammut::utils::intToString(pid) + ".ipc";
}

inline std::string getInstrumentationConnectionChannel(){
    return std::string("ipc://") + getInstrumentationConnectionChannelPath();
}

inline std::string getInstrumentationPidChannel(uint pid){
    return std::string("ipc://") + getInstrumentationPidChannelPath(pid);
}

class InstrumenterHelper: public riff::Application, mammut::utils::NonCopyable{
public:
    InstrumenterHelper(std::pair<nn::socket*, uint> p,
                       size_t numThreads = 1,
                       riff::Aggregator* aggregator = NULL):
        riff::Application(*p.first, p.second, numThreads, aggregator){;}
};

/**
 * Extends riff::Application by providing a way to automatically connect
 * to the nornir manager server.
 */
class Instrumenter: public InstrumenterHelper{
private:
    std::pair<nn::socket*, uint> getChannel(const std::string& parametersFile) const;
    std::pair<nn::socket*, uint> connectPidChannel(const std::string& parametersFile, uint pid) const;
public:
    /**
     * Creates a client for interaction with a local server. The suggestion
     * is to create as soon as possible (i.e. before the beginning of critical
     * parts of application) such that in the meanwhile a connection
     * with the manager can be established.
     * @param parametersFile The file containing the Nornir parameters.
     * @param numThreads The number of threads calling the instrumentation calls.
     * @param aggregator The riff aggregator.
     */
    explicit Instrumenter(const std::string& parametersFile,
                          size_t numThreads = 1,
                          riff::Aggregator* aggregator = NULL);
};

/**
 * A server which can interact with the Instrumenter.
 **/
class InstrumenterServer: public mammut::utils::Thread{
private:
    bool _singleClient;
public:
    InstrumenterServer(bool singleClient = true):_singleClient(singleClient){;}
    void run();
};

}

#endif /* NORNIR_INSTRUMENTER_HPP_ */

/*
 * monitor.cpp
 *
 * Created on: 03/07/2016
 *
 * This file contains the client code for monitoring a local application
 * through instrumentation, and to send the monitoring data to a Nornir
 * manager running in a different process and acting like a server.
 * First of all, we try to open the InstrumentationConnectionChannel channel,
 * whose name  is fixed and known by all the clients and by the server. On
 * this channel, the client sends its PID. Then, all the communications between
 * this specific client (i.e. application) and the server will be done on a
 * separate channel (identified by the PID).
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

#include <sys/types.h>
#include <unistd.h>

#include "parameters.hpp"
#include "instrumenter.h"
#include "instrumenter.hpp"
#include "external/mammut/mammut/mammut.hpp"

#undef DEBUG
#undef DEBUGB

#ifdef DEBUG_INSTRUMENTER
#define DEBUG(x) do { std::cerr << "[Instrumenter] " << x << std::endl; } while (0)
#define DEBUGB(x) do {x} while (0)
#else
#define DEBUG(x)
#define DEBUGB(x)
#endif

namespace nornir{

using namespace std;
using namespace mammut;
using namespace mammut::utils;

std::pair<nn::socket*, uint> Instrumenter::connectPidChannel(const std::string& parametersFile, uint pid) const{
    nn::socket* channel = new nn::socket(AF_SP, NN_PAIR);
    int chid;

    /** Send content length. **/
    chid = channel->connect(getInstrumentationPidChannel(pid).c_str());
    DEBUG("Connected to application channel.");
    assert(chid >= 0);
    vector<string> lines = readFile(parametersFile);
    size_t length = 0;
    for(size_t i = 0; i < lines.size(); i++){
        length += lines.at(i).length();
    }
    length += 1;
    int ret = channel->send(&length, sizeof(length), 0);
    DEBUG("Parameters sent.");
    assert(ret == sizeof(length));

    /** Send content. **/
    string fileContent;
    fileContent.reserve(length);
    for(size_t i = 0; i < lines.size(); i++){
        fileContent.append(lines.at(i));
    }
    fileContent.push_back('\0');
    const char* fileContentC = fileContent.c_str();
    ret = channel->send(fileContentC, length, 0);
    assert(ret == (int) length);

    /** Receive validation result. **/
    ParametersValidation pv;
    ret = channel->recv(&pv, sizeof(pv), 0);
    assert(ret == sizeof(pv));
    DEBUG("Validation results received.");
    if(pv != VALIDATION_OK){
        delete channel;
        throw runtime_error("Invalid adaptivity parameters: " + std::to_string(pv));
    }

    return std::pair<nn::socket*, uint>(channel, chid);
}

std::pair<nn::socket*, uint> Instrumenter::getChannel(const std::string& parametersFile) const{
    /** Send pid, then switch to the pid channel. */
    nn::socket mainChannel(AF_SP, NN_PAIR);
    int mainChid;
    mainChid = mainChannel.connect(getInstrumentationConnectionChannel().c_str());
    if(mainChid < 0){
        throw std::runtime_error("Impossible to connect to Nornir.");
    }
    DEBUG("Connected to main channel.");
    pid_t pid = getpid();
    int ret = mainChannel.send(&pid, sizeof(pid), 0);
    assert(ret == sizeof(pid));
    DEBUG("PID sent.");
    // Wait ack. This is needed because we have to be sure that the pid channel
    // has been created before trying to connect.
    char ack;
    ret = mainChannel.recv(&ack, sizeof(ack), 0);
    assert(ret == sizeof(ack));
    DEBUG("Ack received.");
    mainChannel.shutdown(mainChid);
    DEBUG("Main channel closed.");

    return connectPidChannel(parametersFile, pid);
}

Instrumenter::Instrumenter(const std::string& parametersFile,
                           size_t numThreads,
                           riff::Aggregator *aggregator):
        InstrumenterHelper(getChannel(parametersFile), numThreads, aggregator){
    ;
}

}

extern "C"{
    NornirInstrumenter* nornir_instrumenter_create(const char* parametersFile){
        return reinterpret_cast<NornirInstrumenter*>(new nornir::Instrumenter(parametersFile));
    }

    NornirInstrumenter* nornir_instrumenter_create_with_threads(const char* parametersFile, size_t numThreads){
        return reinterpret_cast<NornirInstrumenter*>(new nornir::Instrumenter(parametersFile, numThreads));
    }

    void nornir_instrumenter_destroy(NornirInstrumenter* instrumenter){
        delete reinterpret_cast<nornir::Instrumenter*>(instrumenter);
    }

    void nornir_instrumenter_begin(NornirInstrumenter* instrumenter){
        reinterpret_cast<nornir::Instrumenter*>(instrumenter)->begin();
    }

    void nornir_instrumenter_begin_with_threads(NornirInstrumenter* instrumenter, size_t threadId){
        reinterpret_cast<nornir::Instrumenter*>(instrumenter)->begin(threadId);
    }

    void nornir_instrumenter_end(NornirInstrumenter* instrumenter){
        reinterpret_cast<nornir::Instrumenter*>(instrumenter)->end();
    }

    void nornir_instrumenter_end_with_threads(NornirInstrumenter* instrumenter, size_t threadId){
        reinterpret_cast<nornir::Instrumenter*>(instrumenter)->end(threadId);
    }

    void nornir_instrumenter_terminate(NornirInstrumenter* instrumenter){
        reinterpret_cast<nornir::Instrumenter*>(instrumenter)->terminate();
    }

    unsigned long nornir_instrumenter_get_execution_time(NornirInstrumenter* instrumenter){
        return reinterpret_cast<nornir::Instrumenter*>(instrumenter)->getExecutionTime();
    }

    unsigned long long nornir_instrumenter_get_total_tasks(NornirInstrumenter* instrumenter){
        return reinterpret_cast<nornir::Instrumenter*>(instrumenter)->getTotalTasks();
    }

    void nornir_instrumenter_set_total_threads(NornirInstrumenter* instrumenter, unsigned int totalThreads){
        return reinterpret_cast<nornir::Instrumenter*>(instrumenter)->setTotalThreads(totalThreads);
    }

    void nornir_instrumenter_set_phase_id(NornirInstrumenter* instrumenter, unsigned int phaseId){
        return reinterpret_cast<nornir::Instrumenter*>(instrumenter)->setPhaseId(phaseId);
    }

    void nornir_instrumenter_mark_inconsistent_samples(NornirInstrumenter* instrumenter){
        return reinterpret_cast<nornir::Instrumenter*>(instrumenter)->markInconsistentSamples();
    }
}

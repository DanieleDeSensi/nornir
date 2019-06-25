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
#include <nornir/nornir.hpp>
#include <nornir/parameters.hpp>
#include <nornir/instrumenter.h>
#include <nornir/instrumenter.hpp>
#include <mammut/mammut.hpp>

#include <sys/types.h>
#include <unistd.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <type_traits>

static_assert(nornir::REQUIREMENT_NUM <= RIFF_MAX_CUSTOM_FIELDS, "Please increase RIFF_MAX_CUSTOM_FIELDS");

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

void Instrumenter::changeRequirement(RequirementType type, double value){
  storeCustomValue(static_cast<size_t>(type), value, 0);
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

using namespace nornir;
using namespace nn;
using namespace std;

#undef DEBUG
#undef DEBUGB

#ifdef DEBUG_MEXT
#define DEBUG(x) do { cerr << "[External Server] " << x << endl << flush; } while (0)
#define DEBUGB(x) do {x;} while(0)
#else
#define DEBUG(x)
#define DEBUGB(x)
#endif

class ApplicationInstance{
public:
    nn::socket channel;
    int chid;
    ManagerInstrumented *manager;

    ApplicationInstance():channel(AF_SP, NN_PAIR), chid(0), manager(NULL){;}
};

static void managerCleanup(Manager* m, std::list<ApplicationInstance*>& instances){
    if(m){
        for(auto it = instances.begin(); it != instances.end(); it++){
            if((*it)->manager == m){
                auto newit = std::next(it);
                DEBUG("Application manager terminated, cleaning.");
                delete ((*it)->manager);
                (*it)->channel.shutdown((*it)->chid);
                ApplicationInstance* ai = (*it);
                instances.erase(it);
                it = newit;
                delete ai;
            }
        }
    }
}

void InstrumenterServer::run(){
    // TODO: at the moment we do not support concurrent applications.
    // In the future the value of this flag should be given by the user
    // when starts this executable.
    //bool multipleApplications = false;

    // Create directory where the channels will be placed and 
    // set permissions so that everyone can access it.
    if(!mammut::utils::existsDirectory(getInstrumentationChannelsPath())){
        if(mkdir(getInstrumentationChannelsPath().c_str(), ACCESSPERMS)){
            throw std::runtime_error("Impossible to create nornir instrumentation channels dir.");
        }
        if(system((std::string("chmod ugo+rwx ") + getInstrumentationChannelsPath()).c_str()) == -1){
            throw std::runtime_error("Impossible to set permission to nornir channel dir.");
        }
    }

    nn::socket mainChannel(AF_SP, NN_PAIR);
    int mainChid;
    mainChid = mainChannel.bind(getInstrumentationConnectionChannel().c_str());
    std::list<ApplicationInstance*> instances;
    // We need to change the rights of the channel because most likely this manager will be
    // executed with sudoers rights (since we need to change frequency etc..).
    // By changing the rights we allow non sudoers users to interact with 
    // the manager on this channel.
    if(system((std::string("chmod ugo+rwx ") + getInstrumentationConnectionChannelPath()).c_str()) == -1){
        throw std::runtime_error("Impossible to set permission to nornir channel.");
    }

    /*
    ManagerMulti mm;
    if(multipleApplications){
        mm.start();
        DEBUG("ManagerMulti started.");
    }
    */

    while(true){
        pid_t pid;
        size_t r = mainChannel.recv(&pid, sizeof(pid), 0);
        assert(r == sizeof(pid));

        DEBUG("Received a request from process " << pid);
        ApplicationInstance* ai = new ApplicationInstance;
        ai->chid = ai->channel.bind(getInstrumentationPidChannel(pid).c_str());
        // We change the rights for same reasons as before.
        if(system((std::string("chmod ugo+rwx ") + getInstrumentationPidChannelPath(pid)).c_str()) == -1){
            throw std::runtime_error("Impossible to set permission to nornir channel.");
        }
        DEBUG("Created app channel.");

        DEBUG("Sending ack.");
        char ack = 0;
        r = mainChannel.send(&ack, sizeof(ack), 0);
        assert(r == sizeof(ack));

        size_t length = 0;
        r = ai->channel.recv(&length, sizeof(length), 0);
        assert(r == sizeof(length));
        char* parameters = new char[length];
        DEBUG("Receiving parameters.");
        r = ai->channel.recv(parameters, length*sizeof(char), 0);
        assert(r == (sizeof(char)*length));
        std::string parametersString(parameters);
        std::ofstream out("parameters.xml");
        out << parametersString;
        out.close();
        delete[] parameters;
        DEBUG("Validating parameters");
        //TODO: If we will let this work for remote machines too, we will need
        // to also send archfile.xml
        Parameters p("parameters.xml");
        ParametersValidation pv = p.validate();
        DEBUG("Sending validation result.");
        r = ai->channel.send(&pv, sizeof(pv), 0);
        assert(r == sizeof(pv));
        ai->manager = new ManagerInstrumented(ai->channel, ai->chid, p);
        ai->manager->start();
        DEBUG("Manager started.");

        // Add to the list of running managers and to the multimanager.
        /*
        if(multipleApplications){
            instances.push_back(ai);
            mm.addManager(ai->manager);
        }
        */

        /** Try to join and delete already terminated managers. **/
        Manager* m;
        // Corunning applications.
        if(true /*!multipleApplications*/){
            // Single application.
            m = ai->manager;
            DEBUG("Joining manager.");
            m->join();
            DEBUG("Manager joined.");
            delete m;
        }else{
            //m = mm.getTerminatedManager();
        }
        managerCleanup(m, instances);
        if(_singleClient){
        	break;
        }
    }
    mainChannel.shutdown(mainChid);
}

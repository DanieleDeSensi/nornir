/*
 * manager.hpp
 *
 * Created on: 23/03/2015
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

/*!
 * @file manager.hpp
 * @brief Implementation of an adaptive fastflow farm.
 */

#ifndef NORNIR_FARM_HPP_
#define NORNIR_FARM_HPP_

#include "configuration.hpp"
#include "ffincs.hpp"
#include "knob.hpp"
#include "parameters.hpp"
#include "selectors.hpp"
#include "node.hpp"
#include "utils.hpp"

//// Nornir's versions of the nodes
#include "external/fastflow/ff/config.hpp"
namespace nornir{
template<typename IN_t, typename OUT_t = IN_t>
struct nrnr_node_t: public AdaptiveNode {
    typedef IN_t  in_type;
    typedef OUT_t out_type;

    nrnr_node_t():
        GO_ON((OUT_t*)ff::FF_GO_ON),
        EOS((OUT_t*)ff::FF_EOS),
        EOSW((OUT_t*)ff::FF_EOSW),
        GO_OUT((OUT_t*)ff::FF_GO_OUT),
        EOS_NOFREEZE((OUT_t*) ff::FF_EOS_NOFREEZE),
        BLK((OUT_t*)ff::FF_BLK), NBLK((OUT_t*)ff::FF_NBLK) {
    }
    OUT_t * const GO_ON,  *const EOS, *const EOSW, *const GO_OUT, *const EOS_NOFREEZE, *const BLK, *const NBLK;
    virtual ~nrnr_node_t()  {}
    virtual OUT_t* svc(IN_t*)=0;
    inline  void *svc(void *task) {
        void* r = svc(reinterpret_cast<IN_t*>(task));
        if(!r){
            TERMINATE_APPLICATION;
        }else{
            return r;
        }
    };
};

template<typename TIN, typename TOUT, typename FUNC>
struct nrnr_node_F: public nrnr_node_t<TIN, TOUT> {
   nrnr_node_F(FUNC f):F(f) {}
   TOUT* svc(TIN* task) {
       TOUT* r = F(task, this);
       if(!r){
           TERMINATE_APPLICATION_TYPED(TOUT);
       }else{
           return r;
       }
   }
   FUNC F;
};
}

#include "external/fastflow/ff/pipeline.hpp"
#include "external/mammut/mammut/module.hpp"
#include "external/mammut/mammut/utils.hpp"
#include "external/mammut/mammut/mammut.hpp"
#include "external/riff/src/riff.hpp"

#include <cmath>
#include <iostream>
#include <limits>

namespace nornir{

class Parameters;
class ManagerMulti;


// How to react to the calibration of other applications
// in order to do not interfere too much with them.
typedef enum{
    // Do nothing
    CALIBRATION_SHRINK_NONE = 0, //TODO: New managers arrived to the multimanager may not have space to calibrate.
    // Move all the threads on the same core.
    CALIBRATION_SHRINK_AGGREGATE,
    // Move all the threads on the same core and pause the process.
    CALIBRATION_SHRINK_PAUSE,
}CalibrationShrink;

/*!
 * \class Manager
 * \brief This class manages the adaptivity parallel applications.
 *
 * This class manages the adaptivity in parallel applications.
 */
class Manager: public mammut::utils::Thread{
    friend class ManagerMulti;
public:
    explicit Manager(Parameters nornirParameters);

    virtual ~Manager();
    /**
     * Function executed by this thread.
     * ATTENTION: The user must not call this one but 'start()'.
     */
    void run();

    /**
     * Forces the manager to terminate.
     **/
    void terminate();


    /**
     * Sets the parameters to be used when simulating nornir. It must be
     * called soon after the object creation.
     * @param samplesFileName The name of the file containing the application
     * samples.
     * ATTENTION: Only used for testing purposes.
     */
    void setSimulationParameters(std::string samplesFileName);
protected:
    // Flag for checking farm termination.
    volatile bool _terminated;

    // The parameters used to take management decisions.
    Parameters _p;

    // The energy counter.
    mammut::energy::Counter* _counter;

    // The task module.
    mammut::task::TasksManager* _task;

    // The topology module.
    mammut::topology::Topology* _topology;

    // The cpufreq module.
    mammut::cpufreq::CpuFreq* _cpufreq;

    // Monitored samples;
    Smoother<MonitoredSample>* _samples;

    // Variations
    Smoother<double>* _variations;

    // The number of tasks processed since the last reconfiguration.
    double _totalTasks;

    // When contract is CONTRACT_COMPLETION_TIME, represent the number of tasks
    // that still needs to be processed by the application.
    uint64_t _remainingTasks;

    // When contract is CONTRACT_COMPLETION_TIME, represent the deadline of
    // the application.
    time_t _deadline;

    // Milliseconds timestamp of the last store of a sample.
    double _lastStoredSampleMs;

    // Inhibition flag.
    bool _inhibited;

    // The current configuration of the application.
    Configuration* _configuration;

    // The configuration selector.
    Selector* _selector;

    // Pid of the monitored process.
    pid_t _pid;

    // Flag indicating if the execution must be simulated.
    bool _toSimulate;

    // Samples to be used for simulation.
    std::list<MonitoredSample> _simulationSamples;

    // When debugging, we print all the monitored samples on this stream
    // ATTENTION: Do NOT protect with DEBUG ifdefs.
    std::ofstream samplesFile;

    mammut::topology::RollbackPoint _topologyRollbackPoint;
    mammut::cpufreq::RollbackPoint _cpufreqRollbackPoint;

    /**
     * Wait for the application to start and
     * sets the ProcessHandler to KnobMappingExternal if necessary.
     */
    virtual void waitForStart() = 0;

    /**
     * Returns a monitored sample.
     * @return A monitored sample sample.
     */
    virtual MonitoredSample getSample() = 0;

    /**
     * Clears the currently stored sample and returns it.
     * @return The current sample.
     **/
    virtual MonitoredSample clearStoredSample(){return getSample();}

    /**
     * Returns the execution time of the application (milliseconds).
     * - For ManagerFastflow, it will be the execution time of the pattern.
     * - For ManagerInstrumented, it will be the execution time between
     *   the construction of the 'Instrumentation' object and the call of the
     *   'terminate' function on it (application side).
     * - For ManagerBlackBox, it will be the execution time from the point where
     *   the manager is created to the end of the application.
     */
    virtual ulong getExecutionTime() = 0;

    /**
     * Manages a configuration change.
     */
    virtual void postConfigurationManagement();

    /**
     * Cleaning after termination.
     */
    virtual void terminationManagement();

    /**
     * Updates the required throughput.
     */
    void updateRequiredThroughput();

    /**
     * Set a specified domain to the highest frequency.
     * @param domain The domain.
     */
    void setDomainToHighestFrequency(const mammut::cpufreq::Domain* domain);

    /**
     * Returns true if the manager doesn't have still to check for a new
     * configuration.
     * @return True if the manager doesn't have still to check for a new
     * configuration.
     */
    bool persist() const;

    /**
     * Locks the knobs according to the selector/predictor.
     */
    void lockKnobs() const;

    /**
     * Initializes the samples.
     * return A samples smoother with no recorded samples.
     */
    Smoother<MonitoredSample>* initSamples() const;

    /**
     * Updates the tasks count.
     * @param sample The workers sample to be used for the update.
     */
    void updateTasksCount(MonitoredSample& sample);

    /**
     * Observes.
     **/
    void observe();

    /**
     * Decides the next configuration.
     */
    KnobsValues decide();

    /**
     * Move to the specified configuration.
     * @param kv The configuration
     */
    void act(KnobsValues kv, bool force = false);

    /**
     * Gets the consumed joules since the last reset and
     * resets the counter.
     * @return The joules consumed since the last reset.
     */
    mammut::energy::Joules getAndResetJoules();

    /**
     * Logs the last observation.
     **/
    void logObservation();

    /**
     * Creates the selector.
     */
    Selector* createSelector() const;
private:
    void inhibit();

    void disinhibit();

    void shrink(CalibrationShrink type);

    void stretch(CalibrationShrink type);

    virtual void shrinkPause() = 0;

    virtual void stretchPause() = 0;

    /**
     * Updates the prediction models by letting them now
     * that there is an interference caused by another application.
     */
    void updateModelsInterference();


    /**
     * Wait until the models have been updated and they are
     * ready to be used.
     */
    void waitModelsInterferenceUpdate();

    /**
     * Returns the vector of physical cores used by the manager.
     * If the manager is still calibrating, it waits until
     * the manager finishes the calibration.
     *
     * @return The vector of physical cores used by the manager.
     */
    std::vector<mammut::topology::PhysicalCoreId> getUsedCores();

    void allowCores(std::vector<mammut::topology::VirtualCoreId> ids);
};

/**
 * Manager for instrumented applications.
 **/
class ManagerInstrumented: public Manager{
private:
    riff::Monitor _monitor;
    void shrinkPause();
    void stretchPause();
public:
    /**
     * Creates an adaptivity manager for an instrumented application.
     * @param riffChannel The name of the riff channel.
     * @param nornirParameters The parameters to be used for
     * adaptivity decisions.
     */
    ManagerInstrumented(const std::string& riffChannel,
                        Parameters nornirParameters);

    /**
     * Creates an adaptivity manager for an instrumented application.
     * @param riffSocket The riff socket.
     * @param chid The channel id.
     * @param nornirParameters The parameters to be used for
     * adaptivity decisions.
     */
    ManagerInstrumented(nn::socket& riffSocket,
                        int chid,
                        Parameters nornirParameters);

    /**
     * Destroyes this adaptivity manager.
     */
    ~ManagerInstrumented();
protected:
    void waitForStart();
    MonitoredSample getSample();
    MonitoredSample clearStoredSample();
    ulong getExecutionTime();
    void terminationManagement();
};

/**
 * Manager for non-instrumented applications,
 * monitored thorugh hardware counters.
 **/
class ManagerBlackBox: public Manager{
private:
    mammut::task::ProcessHandler* _process;
    double _startTime;
    void shrinkPause();
    void stretchPause();
public:
    /**
     * Creates an adaptivity manager for an external NON-INSTRUMENTED
     * application (blackbox).
     * @param pid The identifier of an already running process.
     * @param nornirParameters The parameters to be used for
     * adaptivity decisions.
     */
    ManagerBlackBox(pid_t pid, Parameters nornirParameters);

    /**
     * Destroyes this adaptivity manager.
     */
    ~ManagerBlackBox();

    /**
     * Returns the pid of the monitored process.
     * @return The pid of the monitored process.
     */
    pid_t getPid() const;
protected:
    void waitForStart();
    MonitoredSample getSample();
    ulong getExecutionTime();
};

/**
 * Dummy manager for testing purposes.
 **/
class ManagerTest: public Manager{
private:
    void shrinkPause(){;}
    void stretchPause(){;}
public:
    explicit ManagerTest(Parameters nornirParameters, uint numthreads):Manager(nornirParameters){
        //TODO: Avoid this initialization phase which is common to all the managers.
        Manager::_configuration = new ConfigurationExternal(_p);
        dynamic_cast<KnobVirtualCores*>(_configuration->getKnob(KNOB_VIRTUAL_CORES))->changeMax(numthreads);
        // For instrumented application we do not care if synchronous of not
        // (we count iterations).
        _p.synchronousWorkers = false;
    }
    ~ManagerTest(){;}
protected:
    void waitForStart(){
        dynamic_cast<KnobMappingExternal*>(_configuration->getKnob(KNOB_MAPPING))->setPid(getpid());
    }
    MonitoredSample getSample(){return MonitoredSample();}
    ulong getExecutionTime(){return 0;}
};

/*!
 * \class ManagerFastFlow
 * \brief This class manages the adaptivity in applications written
 * with FastFlow programming framework.
 *
 * This class manages the adaptivity in applications written
 * with FastFlow programming framework.
 */
class ManagerFastFlow: public Manager{
    template <typename I, typename O> friend class FarmBase;
public:
    /**
     * Creates a farm adaptivity manager.
     * @param farm The farm to be managed.
     * @param nornirParameters The parameters to be used for
     * adaptivity decisions.
     */
    ManagerFastFlow(ff::ff_farm<>* farm, Parameters nornirParameters);

    /**
     * Destroyes this adaptivity manager.
     */
    ~ManagerFastFlow();
private:
    // The managed farm.
    ff::ff_farm<>* _farm;

    // The emitter (if present).
    AdaptiveNode* _emitter;

    // The collector (if present).
    AdaptiveNode* _collector;

    // The vector of active workers.
    std::vector<AdaptiveNode*> _activeWorkers;

    void waitForStart();
    MonitoredSample getSample();
    void postConfigurationManagement();
    void terminationManagement();
    ulong getExecutionTime();
    void shrinkPause();
    void stretchPause();
};

/*!
 * \class ManagerFastFlowPipeline
 * \brief This class manages the adaptivity in applications written
 * with pipeline pattern in the FastFlow programming framework.
 *
 * This class manages the adaptivity in applications written
 * with pipeline pattern in the FastFlow programming framework.
 */
class ManagerFastFlowPipeline: public Manager{
    template <typename I, typename O> friend class FarmBase;
public:
    /**
     * Creates a pipeline adaptivity manager.
     * @param pipe The pipeline to be managed.
     * @param nornirParameters The parameters to be used for
     * adaptivity decisions.
     */
    ManagerFastFlowPipeline(ff::ff_pipeline* pipe, std::vector<bool> farmsFlags, Parameters nornirParameters);

    /**
     * Destroyes this adaptivity manager.
     */
    ~ManagerFastFlowPipeline();
private:
    // The managed farm.
    ff::ff_pipeline* _pipe;
    std::vector<bool> _farmsFlags;
    std::vector<KnobVirtualCoresFarm*> _farmsKnobs;
    std::vector<AdaptiveNode*> _activeWorkers;
    std::vector<std::vector<AdaptiveNode*>> _allWorkers;
    std::vector<std::vector<double>> _allowedValues;

    
    void waitForStart();
    MonitoredSample getSample();
    void postConfigurationManagement();
    ulong getExecutionTime();
    void shrinkPause();
    void stretchPause();
};

}

#endif /* NORNIR_FARM_HPP_ */

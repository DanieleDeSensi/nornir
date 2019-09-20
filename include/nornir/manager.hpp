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
 * @brief Implementation of an adaptive manager.
 */

#ifndef NORNIR_MANAGER_HPP_
#define NORNIR_MANAGER_HPP_

#include <nornir/parameters.hpp>
#include <nornir/node.hpp>
#include <nornir/utils.hpp>

#include <mammut/module.hpp>
#include <mammut/utils.hpp>
#include <mammut/mammut.hpp>
#include <riff/riff.hpp>

#include <cmath>
#include <iostream>
#include <limits>
#include <list>

namespace nornir{
class Configuration;
class Selector;
class KnobVirtualCoresFarm;
class Parameters;
class ManagerMulti;

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

    // Flag indicating if the underlying hardware is an Heterogeneous Multiprocessor
    bool _isHMP;

    // Number of different HMP domains
    uint _numHMP;

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
     * Decides the next configuration and moves to the specified configuration.
     */
    void decideAndAct(bool force = false);

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
};

/**
 * Manager for instrumented applications.
 **/
class ManagerInstrumented: public Manager{
private:
    riff::Monitor _monitor;
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
public:
    explicit ManagerTest(Parameters nornirParameters, uint numthreads);
    ~ManagerTest(){;}
protected:
    void waitForStart();
    MonitoredSample getSample(){return MonitoredSample();}
    ulong getExecutionTime(){return 0;}
};

}

#endif /* NORNIR_MANAGER_HPP_ */

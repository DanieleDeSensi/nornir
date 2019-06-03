/*
 * manager-multi.hpp
 *
 * Created on: 22/07/2016
 *
 * This file contains the code of a Global Manager capable of controlling a
 * set of managers. It tries to synchronise their calibration phases and to
 * let them run at the same time on the same resources.
 * Since the resources are limited and not all the managers requirements may
 * be satisfied, this global manager tries to be fair and to minimise the
 * disadvantages for each single manager.
 *
 * The rationale is that individual applications can require PERF_* type
 * contracts (would not have sense for individual applications to require
 * POWER_* type contracts in multiprogrammed environment, since it has sense
 * to talk about power cap only at the level of the entire system. If the
 * application is the only one running into the system of course it makes
 * sense to require POWER_* type contracts).
 *
 * On top of the individual contracts by the application, the global manager
 * can add a POWER_* type contract (if power cap is +infinity only the single
 * applications contract will be considered).
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

#ifndef NORNIR_MANAGER_MULTI_HPP_
#define NORNIR_MANAGER_MULTI_HPP_

#include "manager.hpp"
#include "utils.hpp"

namespace nornir{

using AllocationIndexes = std::vector<size_t>;
using Allocation = std::pair<KnobsValues, double>;
using AllocationFlip = std::pair<double, KnobsValues>;
using ValidAllocations = std::vector<std::pair<double, const AllocationIndexes*> >;

typedef enum{
    QUALITY_ESTIMATION_PERFORMANCE = 0,
    QUALITY_ESTIMATION_PREFERENCE,
}QualityEstimation;

typedef struct ManagerMultiConfiguration{
    double powerCap;
    bool useVirtualCores;
    CalibrationShrink shrink;
    QualityEstimation qualityEstimation;

    ManagerMultiConfiguration():
        powerCap(0), useVirtualCores(false),
        shrink(CALIBRATION_SHRINK_AGGREGATE), qualityEstimation(QUALITY_ESTIMATION_PERFORMANCE)
    {;}
}ManagerMultiConfiguration;

// Contains data associated to each running manager.
typedef struct ManagerData{
    // Minimum performance requirement percentage [0, 100].
    double minPerfReqPerc;
    // Minimum performance (absolute value).
    double minPerf;
    // List of virtual cores identifiers used by the manager.
    std::vector<mammut::topology::VirtualCoreId> allocatedCores;
    // A vector of desirable allocations. This vector is sorted
    // from the most preferred to the least preferred.
    // To each KnobValue, we associate the corresponding predicted
    // performance.
    std::vector<Allocation> allocations;
}ManagerData;

class ManagerMulti: public mammut::utils::Thread{
private:
    // Mammut handlers.
    mammut::Mammut _m;
    mammut::topology::Topology* _topology;
    // Configuration parameters.
    ManagerMultiConfiguration _configuration;
    // Queue between the external world and the global manager.
    ff::SWSR_Ptr_Buffer _qIn;
    // Queue between the global manager and the external world.
    ff::SWSR_Ptr_Buffer _qOut;
    // A list of identifiers of all the cores available on the machine.
    std::vector<mammut::topology::VirtualCoreId> _allCores;
    // For each manager, we keep some data.
    std::map<Manager*, ManagerData> _managerData;
    // Contains all the combinations between the possible allocations.
    // It is not sorted according to any specific order. Each combination
    // will be evaluated in order to find the best one according to some metric.
    // For example, suppose to have 3 managers. This vector will contains a
    // certain number of entries. Each entry is composed by 3 numbers, for example:
    // ...
    // 2 9 4
    // ...
    // This specific allocation corresponds to the situation where:
    // - To the first manager in _managerData map, its 2° preferred allocation is assigned
    // - To the second manager in _managerData map, its 9° preferred allocation is assigned
    // - To the third manager in _managerData map, its 4° preferred allocation is assigned
    std::vector<AllocationIndexes> _allocationsCombinations;
    // Moving average on power consumption.
    Smoother<double>* _power;

    /**
     * Allows a specific manager to calibrate.
     * @param m The manager.
     */
    void allowCalibration(Manager* m);

    /**
     * Waits for a manager calibration.
     * @param m The manger.
     */
    void waitForCalibration(Manager* m);

    /**
     * Inhibits all the active managers except the one specified.
     * While inhibited, a manager ignores all the fluctuations on the
     * performance and/or power consumption.
     * @param except The manager that MUST NOT be inhibited. If NULL,
     * all the mangers will be inhibited.
     **/
    void inhibitAll(Manager* except = NULL);

    /**
     * Disinhibits all the active managers.
     **/
    void disinhibitAll();

    /**
     * Shrinks all the managers except the specified one.
     * @param except The manager that MUST NOT be shrunk.
     */
    void shrinkAll(Manager* except);

    /**
     * Returns a list of available virtual cores identifiers.
     * @return A list of available virtual cores identifiers.
     **/
    std::vector<mammut::topology::VirtualCoreId> getAvailableCores() const;

    /**
     * Updates the list of preferred allocations of a specific manager.
     * @param m The manager for which the list of preferred allocations should be updated.
     **/
    void updateAllocations(Manager* const m);

    /**
     * Updates the list of preferred allocations for all the managers.
     **/
    void updateAllocations();

    /**
     * Updates the performance models of all the mangers after that a manager
     * started/terminated its execution.
     */
    void updateModels();

    /**
     * Computes all the possible allocations.
     * @param array A vector containing all the possible values for each element.
     * @param i This is a recursive function. When called for the first time it must be 0.
     * @param accum This is a recursive functin. When called for the first time it must be empty.
     **/
    void combinations(std::vector<AllocationIndexes> array, size_t i, AllocationIndexes accum);

    /**
     * Computes a map of valid allocations, sorted from the worst to the best.
     * @param validAllocations A vector of valid allocations, it will be sorted from
     * the best (highest quality) to the worst (lowest quality).
     */
    void getValidAllocations(ValidAllocations& validAllocations) const;

    /**
     * Finds the best allocation.
     * @return A vector of indexes (one per manager). Each index is the position of
     * the chosen allocation inside the corresponding allocations vector.
     */
    AllocationIndexes findBestAllocation();

    /**
     * Applies the best found allocation.
     **/
    void applyNewAllocation();

    /**
     * Given an allocation vector. Returns its quality.
     * It must be the highest the better.
     * @param indexes The allocation vector.
     * @return The quality of the allocation.
     */
    double getQuality(const AllocationIndexes& indexes) const;

    /**
     * Estimates the power consumption of a given allocation.
     * @param indexes The allocation vector.
     * @return The power consumption estimation.
     */
    double estimatePower(const AllocationIndexes& indexes) const;

    /**
     * Calibrates the application associated to a specific manager
     * (starting it if required).
     * @param m The manager to be calibrated.
     * @param start If true, the manager will be started. If false,
     * the manager is supposed to be already running.
     */
    void calibrate(Manager* m, bool start = false);

    /**
     * Checks if a specified manager is supported.
     * @param m The manager to be checked.
     */
    void checkManagerSupported(Manager* m);

    /**
     * Allows a manager to run on some specific cores.
     * @param m The manager.
     * @param cores The vector containing the virtual cores identifiers.
     */
    void allowCores(Manager* m, const std::vector<mammut::topology::VirtualCoreId>& cores);

    /**
     * Checks if a specified allocation is valid.
     * @param indexes The allocation.
     * @return true if the allocation is valid, false otherwise.
     */
    bool isValidAllocation(const AllocationIndexes& indexes) const;
public:
    /**
     * Creates a global manager.
     * @param configuration Configuration parameters for the global manager.
     **/
    explicit ManagerMulti(ManagerMultiConfiguration configuration = ManagerMultiConfiguration());

    /**
     * Destroys the global manager.
     **/
    ~ManagerMulti();

    /**
     * Adds a manager to be coordinated.
     * @param m The new manager to be coordinated.
     * @param minPerformanceRequired A value in the range [0, 100]
     *        representing the minimum required performance, expressed
     *        as a percentage of the maximum achievable performance.
     *        For example, a value of 80 means that we would like
     *        to do not decrease the performance of the application
     *        below the 80% of its explicit performance requirement.
     *        0 means no specific requirement, i.e. can be decreased
     *        by any percentage.
     **/
    void addManager(Manager* m, double minPerformanceRequired = 0);

    /**
     * Returns a terminated manager if present, NULL otherwise.
     * @return A terminated manager if present, NULL otherwise.
     */
    Manager* getTerminatedManager();

    /**
     * This function contains the code of the global manager.
     * ATTENTION: Do not use this function directly. Start
     * the manager through the 'start' function.
     * TODO: Make private?
     **/
    void run();
};

}

#endif /* NORNIR_MANAGER_MULTI_HPP_ */

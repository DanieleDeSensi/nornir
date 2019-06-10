/*
 * configuration.hpp
 *
 * Created on: 05/12/2015
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

#ifndef NORNIR_CONFIGURATION_HPP_
#define NORNIR_CONFIGURATION_HPP_

#include "ffincs.hpp"
#include "knob.hpp"
#include "node.hpp"
#include "parameters.hpp"
#include "stats.hpp"
#include "trigger.hpp"

#include "external/mammut/mammut/utils.hpp"
#include "external/fastflow/ff/pipeline.hpp"

namespace nornir{

class ManagerTest;

class Configuration: public mammut::utils::NonCopyable {
    friend class ManagerTest;
protected:
    Knob* _knobs[KNOB_NUM];
    Trigger* _triggers[TRIGGER_TYPE_NUM];
    uint _numServiceNodes;
private:
    const Parameters& _p;
    bool _combinationsCreated;
    std::vector<KnobsValues> _combinations;
    bool _knobsChangeNeeded;
    ReconfigurationStats _reconfigurationStats;

    void combinations(std::vector<std::vector<double> > array, size_t i,
                      std::vector<double> accum);

    bool virtualCoresWillChange(const KnobsValues& values) const;

    ticks startReconfigurationStatsKnob() const;

    void stopReconfigurationStatsKnob(ticks start, KnobType type, bool vcChanged);

    ticks startReconfigurationStatsTotal() const;

    void stopReconfigurationStatsTotal(ticks start);
public:
    explicit Configuration(const Parameters& p);

    virtual ~Configuration() = 0;

    /**
     * Returns true if the values of this configuration are equal to those
     * passed as parameters, false otherwise.
     * @return true if the values of this configuration are equal to those
     * passed as parameters, false otherwise.
     */
    bool equal(const KnobsValues& values) const;

    /**
     * Given some knobs values, returns the corresponding real values.
     * @param values The knobs values.
     * @return The corresponding real values.
     */
    KnobsValues getRealValues(const KnobsValues& values) const;

    /**
     * Returns true if the knobs values need to be changed, false otherwise.
     * @return True if the knobs values need to be changed, false otherwise.
     */
    bool knobsChangeNeeded() const;

    /**
     * Creates all the possible knobs combinations.
     */
    void createAllRealCombinations();

    /**
     * Gets all the possible combinations of knobs values.
     * @return A std::vector containing all the possible combinations
     *         of knobs values.
     */
    const std::vector<KnobsValues>& getAllRealCombinations() const;

    /**
     * Sets the highest frequency to reduce the reconfiguration time.
     */
    void setFastReconfiguration();

    /**
     * Returns a specified knob.
     * @param t The type of the knob to return.
     * @return The specified knob.
     */
    Knob* getKnob(KnobType t) const;

    /**
     * Sets all the knobs to their maximum.
     */
    void maxAllKnobs();

    /**
     * Returns the real value of a specific knob.
     * @param t The type of the knob.
     * @return The real value of the specified knob.
     */
    double getRealValue(KnobType t) const;

    /**
     * Returns the real values for all the knobs.
     * @return The real values for all the knobs.
     */
    KnobsValues getRealValues() const;

    /**
     * Sets values for the knobs (may be relative or real).
     * @param values The values of the knobs.
     */
    void setValues(const KnobsValues& values);

    /**
     * Triggers the triggers.
     */
    void trigger();

    /**
     * Returns the reconfiguration statistics.
     * @return The reconfiguration statistics.
     */
    inline ReconfigurationStats getReconfigurationStats() const{
        return _reconfigurationStats;
    }

    inline uint getNumServiceNodes() const{return _numServiceNodes;}
};

class ConfigurationExternal: public Configuration{
public:
    explicit ConfigurationExternal(const Parameters& p);
};

class ConfigurationFarm: public Configuration{
public:
    ConfigurationFarm(const Parameters& p,
                      Smoother<MonitoredSample> const* samples,
                      AdaptiveNode* emitter,
                      std::vector<AdaptiveNode*> workers,
                      AdaptiveNode* collector,
                      ff::ff_gatherer* gt,
                      volatile bool* terminated);
};

class ConfigurationPipe: public Configuration{
public:
    ConfigurationPipe(const Parameters& p,
                      Smoother<MonitoredSample> const* samples,
                      std::vector<KnobVirtualCoresFarm*> farms,
                      std::vector<std::vector<double>> allowedValues);
};

std::vector<AdaptiveNode*> convertWorkers(ff::svector<ff::ff_node*> w);

}
#endif /* NORNIR_CONFIGURATION_HPP_ */

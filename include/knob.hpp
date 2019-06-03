/*
 * knob.hpp
 *
 * Created on: 02/11/2015
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

#ifndef NORNIR_KNOB_HPP_
#define NORNIR_KNOB_HPP_

#include "ffincs.hpp"
#include "node.hpp"

namespace nornir{

class Knob: public mammut::utils::NonCopyable{
public:
    Knob():_realValue(-1), _locked(false){;}

    /**
     * Computes the real value corresponding to a specific
     * relative value.
     * @param relative The relative value.
     * @param real The real value.
     * @return True if there is a real value, false otherwise (e.g. because this
     * knob has no possible values to be set).
     */
    bool getRealFromRelative(double relative, double& real) const;

    /**
     * Changes the value of this knob.
     * @param v Is a value in the range [0, 100], representing the value to be
     *          set for the knob.
     */
    void setRelativeValue(double v);

    /**
     * Changes the value of this knob.
     * @param v Is the real value of the knob.
     */
    void setRealValue(double v);

    /**
     * Sets this knob to its maximum.
     */
    void setToMax();

    /**
     * Locks this knob to a relative value in the range [0, 100].
     * @param value The value in the range [0, 100].
     */
    void lock(double v);

    /**
     * Locks this knob to the maximum value.
     */
    void lockToMax();

    /**
     * Locks this knob to the minimum value.
     */
    void lockToMin();

    /**
     * Checks if the knob is locked.
     * @return True if the knob is locked, false otherwise.
     */
    bool isLocked() const;

    /**
     * Returns the current real value of this knob.
     * @return The current real value of this knob.
     */
    double getRealValue() const;

    /**
     * Returns a vector of allowed values for this knob.
     * @return A vector of allowed values for this knob.
     */
    std::vector<double> getAllowedValues() const;

    virtual ~Knob(){;}
protected:
    /**
     * Changes the value of this knob.
     * @param v Is the real value that this knob will have.
     */
    virtual void changeValue(double v) = 0;

    double _realValue;
    bool _locked;
    std::vector<double> _knobValues;
};

class KnobVirtualCores: public Knob{
protected:
    Parameters _p;
public:
    explicit KnobVirtualCores(Parameters p);
    void changeValue(double v);
    /**
     * Changes the maximum allowed value for this knob.
     * If it is lower than the current value, the current value is changed as
     * well.
     * @param v The new maximum allowed value.
     */
    void changeMax(double v);
};

class KnobVirtualCoresFarm: public KnobVirtualCores{
    friend class ManagerFastFlow;
    template <typename I, typename O> friend class FarmBase;
public:
    KnobVirtualCoresFarm(Parameters p,
                  AdaptiveNode* emitter, AdaptiveNode* collector,
                  ff::ff_gatherer* gt,
                  const std::vector<AdaptiveNode*>& workers,
                  const volatile bool* terminated);

    void changeValue(double v);
    std::vector<double> getAllowedValues() const;
    std::vector<AdaptiveNode*> getActiveWorkers() const;
private:
    /**
     * Prepares the nodes to freeze.
     */
    void prepareToFreeze();

    /**
     * Freezes all the nodes.
     */
    void freeze();

    /**
     * Prepares the nodes to run.
     * @param numWorkers The new number of workers.
     */
    void prepareToRun(uint numWorkers);

    /**
     * Runs the farm with a new number of workers.
     * @param numWorkers The new number of workers.
     */
    void run(uint numWorkers);

    /**
     * Notifies a change in the number of workers to all the nodes.
     * @param numWorkers The new number of workers.
     */
    void notifyNewConfiguration(uint numWorkers);

    AdaptiveNode* _emitter;
    AdaptiveNode* _collector;
    ff::ff_gatherer* _gt;
    const std::vector<AdaptiveNode*> _allWorkers;
    std::vector<AdaptiveNode*> _activeWorkers;
    const volatile bool* _terminated;
};

class KnobVirtualCoresPipe: public KnobVirtualCores{
    friend class ManagerFastFlow;
    template <typename I, typename O> friend class FarmBase;
public:
    KnobVirtualCoresPipe(Parameters p,
                         std::vector<KnobVirtualCoresFarm*> farms,
                         std::vector<std::vector<double> > allowedValues);

    void changeValue(double v);
    std::vector<AdaptiveNode*> getActiveWorkers() const;
private:
    std::vector<KnobVirtualCoresFarm*> _farms;
    std::vector<std::vector<double>> _allowedValues;
};

class KnobHyperThreading: public Knob{
public:
    explicit KnobHyperThreading(Parameters p);
    void changeValue(double v);
};

// ATTENTION: Update enumString in knob.cpp
typedef enum{
    MAPPING_TYPE_LINEAR = 0,
    MAPPING_TYPE_INTERLEAVED, // One per CPU, round robin
    //MAPPING_TYPE_CACHE_OPTIMAL,
    MAPPING_TYPE_NUM // ATTENTION: This must be the last value.
}MappingType;

class KnobMapping: public Knob{
public:
    KnobMapping(const Parameters& p,
                const KnobVirtualCores& knobCores,
                const KnobHyperThreading& knobHyperThreading);
    void changeValue(double v);

    virtual void move(const std::vector<mammut::topology::VirtualCore*>& vcOrder) = 0;

    void setAllowedCores(std::vector<mammut::topology::VirtualCore*> vc);

    std::vector<mammut::topology::VirtualCore*> getAllowedCores() const;

    bool isAllowed(mammut::topology::VirtualCore*) const;

    const std::vector<mammut::topology::VirtualCore*>& getActiveVirtualCores() const;
    const std::vector<mammut::topology::VirtualCore*>& getUnusedVirtualCores() const;
protected:
    const Parameters& _p;
    const KnobVirtualCores& _knobCores;
    const KnobHyperThreading& _knobHyperThreading;

    virtual size_t getNumVirtualCores();
private:
    std::vector<mammut::topology::VirtualCore*> _activeVirtualCores;
    std::vector<mammut::topology::VirtualCore*> _unusedVirtualCores;
    mammut::topology::Topology* _topologyHandler;
    std::vector<mammut::topology::VirtualCore*> _allowedVirtualCores;

    std::vector<mammut::topology::VirtualCore*> computeVcOrderLinear();
    std::vector<mammut::topology::VirtualCore*> computeVcOrderInterleaved();
};

class KnobMappingExternal: public KnobMapping{
private:
    mammut::task::ProcessHandler* _processHandler;
public:
    KnobMappingExternal(const Parameters& p,
                const KnobVirtualCores& knobCores,
                const KnobHyperThreading& knobHyperThreading);

    void setPid(pid_t pid);
    void setProcessHandler(mammut::task::ProcessHandler* processHandler);
    void move(const std::vector<mammut::topology::VirtualCore*>& vcOrder);
};

class KnobMappingFarm: public KnobMapping{
private:
    AdaptiveNode* _emitter;
    AdaptiveNode* _collector;
protected:
    size_t getNumVirtualCores();
public:
    KnobMappingFarm(const Parameters& p,
                const KnobVirtualCoresFarm& knobCores,
                const KnobHyperThreading& knobHyperThreading,
                AdaptiveNode* emitter,
                AdaptiveNode* collector);

    void move(const std::vector<mammut::topology::VirtualCore*>& vcOrder);
};

class KnobFrequency: public Knob{
public:
    KnobFrequency(Parameters p, const KnobMapping& knobMapping);
    ~KnobFrequency();
    void changeValue(double v);
private:
    void applyUnusedVCStrategySame(const std::vector<mammut::topology::VirtualCore*>& unusedVc, mammut::cpufreq::Frequency v);
    void applyUnusedVCStrategyOff(const std::vector<mammut::topology::VirtualCore*>& unusedVc);
    void applyUnusedVCStrategyLowestFreq(const std::vector<mammut::topology::VirtualCore*>& unusedVc);
    void applyUnusedVCStrategy(mammut::cpufreq::Frequency v);

    Parameters _p;
    const KnobMapping& _knobMapping;
    mammut::cpufreq::CpuFreq* _frequencyHandler;
    mammut::topology::Topology* _topologyHandler;
};

class KnobClkMod: public Knob{
private:
    const KnobMapping& _knobMapping;
public:
    explicit KnobClkMod(Parameters p, const KnobMapping& knobMapping);
    void changeValue(double v);
};

class KnobClkModEmulated: public Knob{
private:
    Parameters _p;
    mammut::task::ProcessHandler* _processHandler;
public:
    explicit KnobClkModEmulated(Parameters p);
    void setPid(pid_t pid);
    void setProcessHandler(mammut::task::ProcessHandler* processHandler);
    void changeValue(double v);
};

/**
 * Knobs values can be:
 *  - Real: e.g. for workers it will get value between 1 and the maximum number of cores.
 *  - Relative: They will always assume values in the range [0.0, 100.0]
 */
typedef enum{
    KNOB_VALUE_UNDEF = 0,
    KNOB_VALUE_REAL,
    KNOB_VALUE_RELATIVE
}KnobValueType;

std::string knobTypeToString(KnobType kv);

class KnobsValues{
private:
    KnobValueType _type;
    double _values[KNOB_NUM];
public:
    inline void reset(){
        for(size_t i = 0; i < KNOB_NUM; i++){
            _values[i] = 0;
        }
    }

    explicit KnobsValues(KnobValueType type = KNOB_VALUE_UNDEF):_type(type){
        reset();
    }

    void swap(KnobsValues& x){
        using std::swap;

        swap(_type, x._type);
        swap(_values, x._values);
    }

    inline KnobsValues(const KnobsValues& other){
        _type = other._type;
        for(size_t i = 0; i < KNOB_NUM; i++){
            _values[i] = other._values[i];
        }
    }

    inline KnobsValues& operator=(KnobsValues other){
        swap(other);
        return *this;
    }

    inline bool areUndefined() const{return _type == KNOB_VALUE_UNDEF;}

    inline bool areRelative() const{return _type == KNOB_VALUE_RELATIVE;}

    inline bool areReal() const{return _type == KNOB_VALUE_REAL;}

    inline double& operator[](KnobType idx){
        return _values[idx];
    }

    inline double operator[](KnobType idx) const{
        return _values[idx];
    }
};


inline std::ostream& operator<<(std::ostream& os, const KnobsValues& obj){
    os << "[";
    for(size_t i = 0; i < KNOB_NUM; i++){
        os << obj[(KnobType) i] << ", ";
    }
    os << "]";
    return os;
}

inline std::istream& operator>>(std::istream& is, KnobsValues& sample){
    is.ignore(std::numeric_limits<std::streamsize>::max(), '[');
    for(size_t i = 0; i < KNOB_NUM; i++){
        is >> sample[(KnobType) i];
        is.ignore(std::numeric_limits<std::streamsize>::max(), ',');
    }
    is.ignore(std::numeric_limits<std::streamsize>::max(), ']');
    return is;
}

inline bool operator==(const KnobsValues& lhs,
                       const KnobsValues& rhs){
    for(size_t i = 0; i < KNOB_NUM; i++){
        if(lhs[(KnobType) i] !=
           rhs[(KnobType) i]){
            return false;
        }
    }
    return true;
}

inline bool operator!=(const KnobsValues& lhs,
                       const KnobsValues& rhs){
    return !operator==(lhs,rhs);
}

inline bool operator<(const KnobsValues& lhs,
                      const KnobsValues& rhs){
    for(size_t i = 0; i < KNOB_NUM; i++){
        if(lhs[(KnobType) i] <
           rhs[(KnobType) i]){
            return true;
        }else if(lhs[(KnobType) i] >
                 rhs[(KnobType) i]){
            return false;
        }
    }
    return false;
}

inline bool operator>(const KnobsValues& lhs,
                      const KnobsValues& rhs){
    return operator< (rhs,lhs);
}

inline bool operator<=(const KnobsValues& lhs,
                       const KnobsValues& rhs){
    return !operator> (lhs,rhs);
}

inline bool operator>=(const KnobsValues& lhs,
                       const KnobsValues& rhs){
    return !operator< (lhs,rhs);
}

}

#endif /* NORNIR_KNOB_HPP_ */

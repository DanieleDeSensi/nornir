/*
 * configuration.cpp
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

#include "configuration.hpp"
#include <iostream>

#ifdef DEBUG_CONFIGURATION
#define DEBUG(x) do { std::cerr << "[Configuration] " << x << std::endl; } while (0)
#define DEBUGB(x) do {x;} while(0)
#else
#define DEBUG(x)
#define DEBUGB(x)
#endif

namespace nornir{

Configuration::Configuration(const Parameters& p):
    _numServiceNodes(0),
    _p(p), _combinationsCreated(false),
    _knobsChangeNeeded(false){
    memset(_knobs, 0, sizeof(_knobs));
    memset(_triggers, 0, sizeof(_triggers));
}

Configuration::~Configuration(){
    for(size_t i = 0; i < KNOB_NUM; i++){
        if(_knobs[i]){
            delete _knobs[i];
        }
    }
    for(size_t i = 0; i < TRIGGER_TYPE_NUM; i++){
        if(_triggers[i]){
            delete _triggers[i];
        }
    }
}

//TODO: Works even if a std::vector is empty? (i.e. a knob has no values)
void Configuration::combinations(std::vector<std::vector<double> > array,
                                 size_t i, std::vector<double> accum){
    if(i == array.size()){
        KnobsValues kv(KNOB_VALUE_REAL);
        for(size_t i = 0; i < KNOB_NUM; i++){
            kv[(KnobType) i] = accum.at(i);
        }
        _combinations.push_back(kv);
    }else{
        std::vector<double> row = array.at(i);
        for(size_t j = 0; j < row.size(); ++j){
            std::vector<double> tmp(accum);
            tmp.push_back(row[j]);
            combinations(array, i+1, tmp);
        }
    }
}

bool Configuration::equal(const KnobsValues& values) const{
    KnobsValues real = getRealValues(values);
    for(size_t i = 0; i < KNOB_NUM; i++){
        KnobType kt = static_cast<KnobType>(i);
        if(real[kt] != getRealValue(kt)){
            return false;
        }
    }
    return true;
}

KnobsValues Configuration::getRealValues(const KnobsValues& values) const{
    if(values.areReal()){
        return values;
    }else{
        KnobsValues r(KNOB_VALUE_REAL);
        double real;
        for(size_t i = 0; i < KNOB_NUM; i++){
            if(_knobs[(KnobType) i]->getRealFromRelative(values[(KnobType) i], real)){
                r[(KnobType) i] = real;
            }
        }
        return r;
    }
}

bool Configuration::knobsChangeNeeded() const{
    return _knobsChangeNeeded;
}

void Configuration::createAllRealCombinations(){
    std::vector<std::vector<double>> values;
    std::vector<double> accum;
    _combinations.clear();
    for(size_t i = 0; i < KNOB_NUM; i++){
        values.push_back(_knobs[i]->getAllowedValues());
        if(!_knobs[i]->isLocked()){
            _knobsChangeNeeded = true;
        }
    }
    combinations(values, 0, accum);
    _combinationsCreated = true;
}

const std::vector<KnobsValues>& Configuration::getAllRealCombinations() const{
    if(!_combinationsCreated){
        throw std::runtime_error("[configuration.cpp] Combinations not created yet.");
    }
    return _combinations;
}

void Configuration::setFastReconfiguration(){
    if(_p.fastReconfiguration){
        dynamic_cast<KnobFrequency*>(_knobs[KNOB_FREQUENCY])->setRelativeValue(100.0);
    }
}

Knob* Configuration::getKnob(KnobType t) const{
    return _knobs[t];
}

void Configuration::maxAllKnobs(){
    DEBUG("Maxing all the knobs.");
    for(size_t i = 0; i < KNOB_NUM; i++){
        _knobs[(KnobType) i]->setToMax();
    }
}

double Configuration::getRealValue(KnobType t) const{
    return _knobs[t]->getRealValue();
}

KnobsValues Configuration::getRealValues() const{
    KnobsValues kv(KNOB_VALUE_REAL);
    for(size_t i = 0; i < KNOB_NUM; i++){
        kv[(KnobType) i] = getRealValue((KnobType) i);
    }
    return kv;
}

bool Configuration::virtualCoresWillChange(const KnobsValues& values) const{
    KnobsValues real = getRealValues(values);
    return real[KNOB_VIRTUAL_CORES] != _knobs[KNOB_VIRTUAL_CORES]->getRealValue();
}

ticks Configuration::startReconfigurationStatsKnob() const{
    if(_p.statsReconfiguration){
        return getticks();
    }
    return 0;
}

ticks Configuration::startReconfigurationStatsTotal() const{
    return startReconfigurationStatsKnob();
}

void Configuration::stopReconfigurationStatsKnob(ticks start, KnobType type,
                                                  bool vcChanged){
    if(_p.statsReconfiguration){
        if(type == KNOB_VIRTUAL_CORES && !vcChanged){
            // We do not add statistics about workers reconfiguration
            // since they did not changed.
            ;
        }else{
            double ms = ticksToMilliseconds(getticks() - start,
                                            _p.archData.ticksPerNs);
            _reconfigurationStats.addSample(type, ms);
        }
    }
}

void Configuration::stopReconfigurationStatsTotal(ticks start){
    if(_p.statsReconfiguration){
        double ms = ticksToMilliseconds(getticks() - start,
                                        _p.archData.ticksPerNs);
        _reconfigurationStats.addSampleTotal(ms);
    }
}

void Configuration::setValues(const KnobsValues& values){
    bool vcChanged = virtualCoresWillChange(values);

    ticks totalStart = startReconfigurationStatsTotal();
    setFastReconfiguration();
    for(size_t i = 0; i < KNOB_NUM; i++){
        ticks reconfigurationStart = startReconfigurationStatsKnob();

        // Start of the real reconfiguration
        if(values.areReal()){
            if(_knobs[i]->isLocked()){
                _knobs[i]->setRealValue(_knobs[i]->getRealValue());
            }else{
                _knobs[i]->setRealValue(values[(KnobType)i]);
            }
        }else if(values.areRelative()){
            if(_knobs[i]->isLocked()){
                // Any relative value would be ok since the knob only has one
                // possible value.
                _knobs[i]->setRelativeValue(0);
            }else{
                _knobs[i]->setRelativeValue(values[(KnobType)i]);
            }
        }else{
            throw std::runtime_error("KnobsValues with undefined type.");
        }
        // End of the real reconfiguration
        stopReconfigurationStatsKnob(reconfigurationStart, (KnobType) i, vcChanged);
    }

    stopReconfigurationStatsTotal(totalStart);

    DEBUG("Changed knobs values.");
}

void Configuration::trigger(){
    for(size_t i = 0; i < TRIGGER_TYPE_NUM; i++){
        if(_triggers[i]){
            _triggers[i]->trigger();
        }
    }
}


ConfigurationExternal::ConfigurationExternal(const Parameters& p):
        Configuration(p){
    _knobs[KNOB_VIRTUAL_CORES] = new KnobVirtualCores(p);
    _knobs[KNOB_HYPERTHREADING] = new KnobHyperThreading(p);
    _knobs[KNOB_MAPPING] = new KnobMappingExternal(p,
                                                *dynamic_cast<KnobVirtualCores*>(_knobs[KNOB_VIRTUAL_CORES]),
                                                *dynamic_cast<KnobHyperThreading*>(_knobs[KNOB_HYPERTHREADING]));
    _knobs[KNOB_FREQUENCY] = new KnobFrequency(p,
                                               *dynamic_cast<KnobMappingExternal*>(_knobs[KNOB_MAPPING]));
    if(p.clockModulationEmulated){
        _knobs[KNOB_CLKMOD] = new KnobClkModEmulated(p);
    }else{
        _knobs[KNOB_CLKMOD] = new KnobClkMod(p, *dynamic_cast<KnobMappingExternal*>(_knobs[KNOB_MAPPING]));
    }

    _triggers[TRIGGER_TYPE_Q_BLOCKING] = NULL;
}

ConfigurationFarm::ConfigurationFarm(const Parameters& p,
                                     Smoother<MonitoredSample> const* samples,
                                     AdaptiveNode* emitter,
                                     std::vector<AdaptiveNode*> workers,
                                     AdaptiveNode* collector,
                                     ff::ff_gatherer* gt,
                                     volatile bool* terminated):
        Configuration(p){
    _knobs[KNOB_VIRTUAL_CORES] = new KnobVirtualCoresFarm(p,
                                          emitter, collector, gt, workers,
                                          terminated);

    _knobs[KNOB_HYPERTHREADING] = new KnobHyperThreading(p);
    _knobs[KNOB_MAPPING] = new KnobMappingFarm(p,
                                               *dynamic_cast<KnobVirtualCoresFarm*>(_knobs[KNOB_VIRTUAL_CORES]),
                                               *dynamic_cast<KnobHyperThreading*>(_knobs[KNOB_HYPERTHREADING]),
                                               emitter, collector);
    _knobs[KNOB_FREQUENCY] = new KnobFrequency(p,
                                               *dynamic_cast<KnobMappingFarm*>(_knobs[KNOB_MAPPING]));
    if(p.clockModulationEmulated){
        _knobs[KNOB_CLKMOD] = new KnobClkModEmulated(p);
    }else{
        _knobs[KNOB_CLKMOD] = new KnobClkMod(p, *dynamic_cast<KnobMappingExternal*>(_knobs[KNOB_MAPPING]));
    }

    _triggers[TRIGGER_TYPE_Q_BLOCKING] = new TriggerQBlocking(p.triggerQBlocking,
                                                              p.thresholdQBlocking,
                                                              p.thresholdQBlockingBelt,
                                                              samples,
                                                              emitter);

    if(emitter){++_numServiceNodes;}
    if(collector){++_numServiceNodes;}
}

ConfigurationPipe::ConfigurationPipe(const Parameters& p,
                                    Smoother<MonitoredSample> const* samples,
                                    std::vector<KnobVirtualCoresFarm*> farms,
                                    std::vector<std::vector<double>> allowedValues):
        Configuration(p){
    _knobs[KNOB_VIRTUAL_CORES] = new KnobVirtualCoresPipe(p, farms, allowedValues);

    _knobs[KNOB_HYPERTHREADING] = new KnobHyperThreading(p);
    _knobs[KNOB_MAPPING] = new KnobMappingExternal(p,
                                                *dynamic_cast<KnobVirtualCores*>(_knobs[KNOB_VIRTUAL_CORES]),
                                                *dynamic_cast<KnobHyperThreading*>(_knobs[KNOB_HYPERTHREADING]));
    _knobs[KNOB_FREQUENCY] = new KnobFrequency(p,
                                               *dynamic_cast<KnobMappingExternal*>(_knobs[KNOB_MAPPING]));

    if(p.clockModulationEmulated){
        _knobs[KNOB_CLKMOD] = new KnobClkModEmulated(p);
    }else{
        _knobs[KNOB_CLKMOD] = new KnobClkMod(p, *dynamic_cast<KnobMappingExternal*>(_knobs[KNOB_MAPPING]));
    }

    _triggers[TRIGGER_TYPE_Q_BLOCKING] = NULL;
}


std::vector<AdaptiveNode*> convertWorkers(ff::svector<ff::ff_node*> w){
    std::vector<AdaptiveNode*> r;
    for(size_t i = 0; i < w.size(); i++){
        r.push_back(dynamic_cast<AdaptiveNode*>(w[i]));
    }
    return r;
}

}

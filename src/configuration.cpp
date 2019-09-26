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

#include <nornir/configuration.hpp>
#include <iostream>

#ifdef DEBUG_CONFIGURATION
#define DEBUG(x) do { std::cerr << "[Configuration] " << x << std::endl; } while (0)
#define DEBUGB(x) do {x;} while(0)
#else
#define DEBUG(x)
#define DEBUGB(x)
#endif

namespace nornir{

Configuration::Configuration(const Parameters& p, uint numHMPs):
    _numServiceNodes(0),
    _numHMPs(numHMPs),
    _p(p), _combinationsCreated(false){
    _knobs.resize(numHMPs);
    memset(_triggers, 0, sizeof(_triggers));
}

Configuration::~Configuration(){
    for(auto k : _knobs){
      for(size_t i = 0; i < KNOB_NUM; i++){
          if(k[i]){
              delete k[i];
          }
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

uint Configuration::getNumHMP() const{
  return _numHMPs;
}

bool Configuration::equal(const KnobsValues& values) const{
    KnobsValues real = getRealValues(values);
    for(size_t c = 0; c < _numHMPs; c++){
      for(size_t i = 0; i < KNOB_NUM; i++){
          KnobType kt = static_cast<KnobType>(i);
          if(real(c, kt) != getRealValue(c, kt)){
              return false;
          }
      }
    }
    return true;
}

KnobsValues Configuration::getRealValues(const KnobsValues& values) const{
    if(values.areReal()){
        return values;
    }else{
        KnobsValues r(KNOB_VALUE_REAL, _numHMPs);
        double real;
        for(uint c = 0; c < _knobs.size(); c++){
          for(size_t i = 0; i < KNOB_NUM; i++){
              double relative = values(c, (KnobType) i);
              if(_knobs[c][(KnobType) i]->getRealFromRelative(relative, real)){
                  r(c, (KnobType) i) = real;
              }
          }
        }
        return r;
    }
}

bool Configuration::knobsChangeNeeded() const{
  for(uint c = 0; c < _knobs.size(); c++){
    for(size_t i = 0; i < KNOB_NUM; i++){
        if(!_knobs[c][i]->isLocked()){
            return true;
        }
    }
  }
  return false;
}

void Configuration::createAllRealCombinations(){
    if(_numHMPs > 1){
      throw std::runtime_error("createAllRealCombinations() cannot be used on HMP systems.");
    }
    std::vector<std::vector<double>> values;
    std::vector<double> accum;
    _combinations.clear();
    for(size_t i = 0; i < KNOB_NUM; i++){
        values.push_back(_knobs[0][i]->getAllowedValues());
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
        for(auto k : _knobs){
          dynamic_cast<KnobFrequency*>(k[KNOB_FREQUENCY])->setRelativeValue(100.0);
        }
    }
}

Knob* Configuration::getKnob(KnobType t) const{
    if(_numHMPs > 1){
       throw std::runtime_error("getKnob: When running on HMP please specify the HMPid.");
    }
    return getKnob(0, t);
}

Knob* Configuration::getKnob(uint cpuId, KnobType t) const{
    return _knobs[cpuId][t];
}

void Configuration::maxAllKnobs(){
    DEBUG("Maxing all the knobs.");
    for(auto k: _knobs){
      for(size_t i = 0; i < KNOB_NUM; i++){
        k[(KnobType) i]->setToMax();
      }
    }
}

double Configuration::getRealValue(KnobType t) const{
  if(_numHMPs > 1){
    throw std::runtime_error("Please use getRealValue(cpuId, t) when running on HMP.");
  }
  return getRealValue(0, t);
}

double Configuration::getRealValue(uint cpuId, KnobType t) const{
    return _knobs[cpuId][t]->getRealValue();
}

KnobsValues Configuration::getRealValues() const{
    KnobsValues kv(KNOB_VALUE_REAL, _numHMPs);
    for(size_t c = 0; c < _numHMPs; c++){
      for(size_t i = 0; i < KNOB_NUM; i++){
        kv(c, (KnobType) i) = getRealValue(c, (KnobType) i);
      }
    }
    return kv;
}

bool Configuration::virtualCoresWillChange(const KnobsValues& values) const{
    KnobsValues real = getRealValues(values);
    for(size_t i = 0; i < _knobs.size(); i++){
      if(real(i, KNOB_VIRTUAL_CORES) != _knobs[i][KNOB_VIRTUAL_CORES]->getRealValue()){
        return true;
      }
    }
    return false;
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
    for(size_t k = 0; k < _knobs.size(); k++){
      for(size_t i = 0; i < KNOB_NUM; i++){
          ticks reconfigurationStart = startReconfigurationStatsKnob();

          // Start of the real reconfiguration
          if(values.areReal()){
              if(_knobs[k][i]->isLocked()){
                  _knobs[k][i]->setRealValue(_knobs[k][i]->getRealValue());
              }else{
                  _knobs[k][i]->setRealValue(values(k, (KnobType)i));
              }
          }else if(values.areRelative()){
              if(_knobs[k][i]->isLocked()){
                  // Any relative value would be ok since the knob only has one
                  // possible value.
                  _knobs[k][i]->setRelativeValue(0);
              }else{
                  _knobs[k][i]->setRelativeValue(values(k, (KnobType)i));
              }
          }else{
              throw std::runtime_error("KnobsValues with undefined type.");
          }
          // End of the real reconfiguration
          stopReconfigurationStatsKnob(reconfigurationStart, (KnobType) i, vcChanged);
      }
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


ConfigurationExternal::ConfigurationExternal(const Parameters& p, uint numHMPs):
        Configuration(p, numHMPs){
    for(size_t c = 0; c < _numHMPs; c++){
      _knobs[c][KNOB_VIRTUAL_CORES] = new KnobVirtualCores(p, _numHMPs, c);
      _knobs[c][KNOB_HYPERTHREADING] = new KnobHyperThreading(p, _numHMPs > 1, c);
      _knobs[c][KNOB_MAPPING] = new KnobMappingExternal(p,
                                                  *dynamic_cast<KnobVirtualCores*>(_knobs[c][KNOB_VIRTUAL_CORES]),
                                                  *dynamic_cast<KnobHyperThreading*>(_knobs[c][KNOB_HYPERTHREADING]),
                                                  _numHMPs > 1, c);
      _knobs[c][KNOB_FREQUENCY] = new KnobFrequency(p,
                                                    *dynamic_cast<KnobMappingExternal*>(_knobs[c][KNOB_MAPPING]),
                                                    _numHMPs > 1, c);
      if(p.clockModulationEmulated){
          _knobs[c][KNOB_CLKMOD] = new KnobClkModEmulated(p);
      }else{
          _knobs[c][KNOB_CLKMOD] = new KnobClkMod(p, *dynamic_cast<KnobMappingExternal*>(_knobs[c][KNOB_MAPPING]), _numHMPs > 1, c);
      }
    }

    _triggers[TRIGGER_TYPE_Q_BLOCKING] = NULL;
}

ConfigurationFarm::ConfigurationFarm(const Parameters& p,
                                     Smoother<MonitoredSample> const* samples,
                                     AdaptiveNode* emitter,
                                     std::vector<AdaptiveNode*> workers,
                                     AdaptiveNode* collector,
                                     ff::ff_gatherer* gt,
                                     volatile bool* terminated, uint numHMPs):
        Configuration(p, numHMPs){
    for(size_t c = 0; c < _knobs.size(); c++){
      _knobs[c][KNOB_VIRTUAL_CORES] = new KnobVirtualCoresFarm(p,
                                            emitter, collector, gt, workers,
                                            terminated);

      _knobs[c][KNOB_HYPERTHREADING] = new KnobHyperThreading(p);
      _knobs[c][KNOB_MAPPING] = new KnobMappingFarm(p,
                                                 *dynamic_cast<KnobVirtualCoresFarm*>(_knobs[c][KNOB_VIRTUAL_CORES]),
                                                 *dynamic_cast<KnobHyperThreading*>(_knobs[c][KNOB_HYPERTHREADING]),
                                                 emitter, collector);
      _knobs[c][KNOB_FREQUENCY] = new KnobFrequency(p,
                                                 *dynamic_cast<KnobMappingFarm*>(_knobs[c][KNOB_MAPPING]));
      if(p.clockModulationEmulated){
          _knobs[c][KNOB_CLKMOD] = new KnobClkModEmulated(p);
      }else{
          _knobs[c][KNOB_CLKMOD] = new KnobClkMod(p, *dynamic_cast<KnobMappingExternal*>(_knobs[c][KNOB_MAPPING]));
      }
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
                                    std::vector<std::vector<double>> allowedValues, uint numHMPs):
        Configuration(p, numHMPs){
    for(size_t c = 0; c < _knobs.size(); c++){
      _knobs[c][KNOB_VIRTUAL_CORES] = new KnobVirtualCoresPipe(p, farms, allowedValues);

      _knobs[c][KNOB_HYPERTHREADING] = new KnobHyperThreading(p);
      _knobs[c][KNOB_MAPPING] = new KnobMappingExternal(p,
                                                  *dynamic_cast<KnobVirtualCores*>(_knobs[c][KNOB_VIRTUAL_CORES]),
                                                  *dynamic_cast<KnobHyperThreading*>(_knobs[c][KNOB_HYPERTHREADING]));
      _knobs[c][KNOB_FREQUENCY] = new KnobFrequency(p,
                                                 *dynamic_cast<KnobMappingExternal*>(_knobs[c][KNOB_MAPPING]));

      if(p.clockModulationEmulated){
          _knobs[c][KNOB_CLKMOD] = new KnobClkModEmulated(p);
      }else{
          _knobs[c][KNOB_CLKMOD] = new KnobClkMod(p, *dynamic_cast<KnobMappingExternal*>(_knobs[c][KNOB_MAPPING]));
      }
    }

    _triggers[TRIGGER_TYPE_Q_BLOCKING] = NULL;
}

}

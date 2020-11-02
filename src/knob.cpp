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

#include <nornir/knob.hpp>
#include <nornir/parameters.hpp>

#include <mammut/cpufreq/cpufreq.hpp>
#include <mammut/mammut.hpp>
#include <mammut/utils.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

#undef DEBUG
#undef DEBUGB

#ifdef DEBUG_KNOB
#define DEBUG(x)                                                               \
  do {                                                                         \
    cerr << "[Knob] " << x << endl;                                            \
  } while (0)
#define DEBUGB(x)                                                              \
  do {                                                                         \
    x;                                                                         \
  } while (0)
#else
#define DEBUG(x)
#define DEBUGB(x)
#endif

namespace nornir {

using namespace mammut;
using namespace mammut::cpufreq;
using namespace mammut::utils;
using namespace mammut::topology;

using namespace ff;

using namespace std;

std::string knobTypeToString(KnobType kv) {
  switch (kv) {
  case KNOB_VIRTUAL_CORES: {
    return "Cores";
  } break;
  case KNOB_HYPERTHREADING: {
    return "HTLevel";
  } break;
  case KNOB_MAPPING: {
    return "Mapping";
  } break;
  case KNOB_FREQUENCY: {
    return "Frequency";
  } break;
  case KNOB_CLKMOD: {
    return "ClockModulation";
  } break;
  case KNOB_PFOR_CHUNK: {
    return "FarmGrain";
  } break;
  default: { return "Unknown"; } break;
  }
}

bool Knob::getRealFromRelative(double relative, double &real) const {
  // Maps from the range [0, 100] to the real range.
  vector<double> values = getAllowedValues();
  if (values.size()) {
    uint index = round((double) (values.size() - 1) * (relative / 100.0));
    real = values.at(index);
    return true;
  } else {
    return false;
  }
}

double Knob::getRelativeFromReal(double real) const {
  vector<double> values = getAllowedValues();
  for (size_t i = 0; i < values.size(); i++) {
    if (values[i] == real) {
      if(values.size() == 1){
	return 100.0;
      }else{
	return (((double) i) / (values.size() - 1)) * 100.0;
      }
    }
  }
  return -1;
}

double Knob::getPreviousRealValue(double realValue, uint step) const {
  vector<double> values = getAllowedValues();
  for (size_t i = 0; i < values.size(); i++) {
    if (values[i] == realValue) {
      if (i >= step) {
        return values[i - step];
      } else {
        return values[0];
      }
    }
  }
  std::ostringstream strs;
  strs << realValue;
  throw std::runtime_error(
      "getPreviousRealValue called with non-existing real value " + strs.str());
}

double Knob::getNextRealValue(double realValue, uint step) const {
  vector<double> values = getAllowedValues();
  for (size_t i = 0; i < values.size(); i++) {
    if (values[i] == realValue) {
      if (i + step < values.size()) {
        return values[i + step];
      } else {
        return values[values.size() - 1];
      }
    }
  }
  std::ostringstream strs;
  strs << realValue;
  throw std::runtime_error(
      "getNextRealValue called with non-existing real value " + strs.str());
}

void Knob::setRelativeValue(double v) {
  double real;
  if (getRealFromRelative(v, real)) {
    setRealValue(real);
  }
}

void Knob::setRealValue(double v) {
  if (getAllowedValues().size()) {
    changeValue(v);
    _realValue = v;
  }
}

void Knob::setToMax() {
  setRelativeValue(100.0);
}

void Knob::lock(double v) {
  setRelativeValue(v);
  double real;
  if (getRealFromRelative(v, real)) {
    _knobValues.clear();
    _knobValues.push_back(real);
  }
  _locked = true;
}

void Knob::lockToMax() {
  if (_knobValues.size()) {
    lock(100.0);
  } else {
    _locked = true;
  }
}

void Knob::lockToMin() {
  if (_knobValues.size()) {
    lock(0);
  } else {
    _locked = true;
  }
}

bool Knob::isLocked() const {
  return _locked;
}

double Knob::getRealValue() const {
  return _realValue;
}

std::vector<double> Knob::getAllowedValues() const {
  return _knobValues;
}

KnobVirtualCores::KnobVirtualCores(Parameters p, size_t numHMP, uint cpuId)
    : _cpuId(cpuId), _p(p) {
  bool isHMP = numHMP > 1;
  size_t numVirtualCores = 0;
  if (_p.activeThreads) {
    if (isHMP) {
      numVirtualCores = std::ceil(_p.activeThreads / (double) numHMP);
    } else {
      numVirtualCores = _p.activeThreads;
    }
  } else {
    if (_p.knobHyperthreadingEnabled) {
      if (isHMP) {
        numVirtualCores = _p.mammut.getInstanceTopology()
                              ->getCpu(cpuId)
                              ->getVirtualCores()
                              .size();
      } else {
        numVirtualCores =
            _p.mammut.getInstanceTopology()->getVirtualCores().size();
      }
    } else {
      if (isHMP) {
        numVirtualCores = _p.mammut.getInstanceTopology()
                              ->getCpu(cpuId)
                              ->getPhysicalCores()
                              .size();
      } else {
        numVirtualCores =
            _p.mammut.getInstanceTopology()->getPhysicalCores().size();
      }
    }
  }
  changeMax(numVirtualCores);
  _realValue = numVirtualCores;
}

void KnobVirtualCores::changeValue(double v) {
  ;
}

void KnobVirtualCores::changeMax(double v) {
  _knobValues.clear();
  for (size_t i = 0; i < v; i++) {
    if (!utils::contains(_p.disallowedNumCores, (uint) i + 1)) {
      _knobValues.push_back(i + 1);
    }
  }
  if (_realValue > v) {
    _realValue = v;
  }
}

#ifdef DEBUG_KNOB
std::ostream &operator<<(std::ostream &out,
                         const std::vector<AdaptiveNode *> &v) {
  out << "[";
  for (size_t i = 0; i < v.size(); i++) {
    out << (v.at(i))->getOSThreadId() << ", ";
  }
  out << "]";
  return out;
}

std::ostream &operator<<(std::ostream &out,
                         const std::vector<VirtualCore *> &v) {
  out << "[";
  for (size_t i = 0; i < v.size(); i++) {
    out << (v.at(i))->getVirtualCoreId() << ", ";
  }
  out << "]";
  return out;
}

std::ostream &operator<<(std::ostream &out, const std::vector<Domain *> &v) {
  out << "[";
  for (size_t i = 0; i < v.size(); i++) {
    out << (v.at(i))->getId() << ", ";
  }
  out << "]";
  return out;
}
#endif

KnobVirtualCoresFarm::KnobVirtualCoresFarm(
    Parameters p, AdaptiveNode *emitter, AdaptiveNode *collector,
    ff::ff_gatherer *gt, const std::vector<AdaptiveNode *> &workers,
    const volatile bool *terminated)
    : KnobVirtualCores(p), _emitter(emitter), _collector(collector), _gt(gt),
      _allWorkers(workers), _terminated(terminated) {
  _realValue = _allWorkers.size();
  _knobValues.clear();
  for (size_t i = 0; i < _allWorkers.size(); i++) {
    if (!utils::contains(p.disallowedNumCores, (uint) i + 1)) {
      _knobValues.push_back(i + 1);
    }
  }

  _activeWorkers = _allWorkers;
  DEBUG("Knob workers created.");
}

void KnobVirtualCoresFarm::changeValue(double v) {
#ifdef BLOCKING_MODE
#error                                                                         \
    "If you want to use concurrency throttling, you need to DO NOT define the FastFlow's BLOCKING_MODE macro."
#endif
  if (_p.useConcurrencyThrottling && v != _realValue) {
    DEBUG("[Workers] Changing real value to: " << v);
    freeze();

    if (!*_terminated) {
      _activeWorkers =
          vector<AdaptiveNode *>(_allWorkers.begin(), _allWorkers.begin() + v);

      notifyNewConfiguration(v);
      run(v);
      DEBUG("[Workers] Active Workers: " << _activeWorkers);
    }
  }
}

std::vector<double> KnobVirtualCoresFarm::getAllowedValues() const {
  return _knobValues;
}

std::vector<AdaptiveNode *> KnobVirtualCoresFarm::getActiveWorkers() const {
  return _activeWorkers;
}

void KnobVirtualCoresFarm::prepareToFreeze() {
  DEBUG("[Workers] Preparing the farm for freezing.");
  if (_emitter) {
    _emitter->prepareToFreeze();
  }
  for (size_t i = 0; i < _activeWorkers.size(); i++) {
    _activeWorkers.at(i)->prepareToFreeze();
  }
  if (_collector) {
    _collector->prepareToFreeze();
  }
}

void KnobVirtualCoresFarm::freeze() {
  prepareToFreeze();
  DEBUG("[Workers] Freezing the farm.");
  if (_emitter) {
    _emitter->freezeAll((void *) (_collector ? FF_EOSW : FF_GO_OUT));
  }

  for (size_t i = 0; i < _activeWorkers.size(); i++) {
    DEBUG("Waiting Worker " << i << " to be frozen:");
    _activeWorkers.at(i)->wait_freezing();
    assert(_activeWorkers.at(i)->isfrozen());
    DEBUG("Worker " << i << " frozen:");
  }
  if (_collector) {
    _gt->wait_freezing();
  }
}

void KnobVirtualCoresFarm::prepareToRun(uint numWorkers) {
  DEBUG("[Workers] Running the farm.");
  if (_emitter) {
    _emitter->prepareToRun();
  }
  for (size_t i = 0; i < numWorkers; i++) {
    _activeWorkers.at(i)->prepareToRun();
  }
  if (_collector) {
    _collector->prepareToRun();
  }
}

void KnobVirtualCoresFarm::run(uint numWorkers) {
  prepareToRun(numWorkers);
  DEBUG("[Workers] Running the emitter.");
  if (_emitter) {
    _emitter->thawAll(numWorkers);
  }
  DEBUG("[Workers] Running the collector.");
  if (_collector) {
    _gt->thaw(true, numWorkers);
  }
}

void KnobVirtualCoresFarm::notifyNewConfiguration(uint numWorkers) {
  if (_emitter) {
    _emitter->notifyRethreading(_realValue, numWorkers);
  }
  for (size_t i = 0; i < numWorkers; i++) {
    AdaptiveNode *w = _activeWorkers.at(i);
    w->notifyRethreading(_realValue, numWorkers);
  }
  if (_collector) {
    _collector->notifyRethreading(_realValue, numWorkers);
  }
}

KnobVirtualCoresPipe::KnobVirtualCoresPipe(
    Parameters p, std::vector<KnobVirtualCoresFarm *> farms,
    std::vector<std::vector<double>> allowedValues)
    : KnobVirtualCores(p), _farms(farms), _allowedValues(allowedValues) {
  _knobValues.clear();
  for (size_t i = 0; i < allowedValues.size(); i++) {
    double sum = 0;
    for (auto d : allowedValues[i]) {
      sum += d;
    }
    _knobValues.push_back(sum);
  }
  _realValue = _knobValues.front();
  DEBUG("KnobVirtualCoresPipe created.");
}

void KnobVirtualCoresPipe::changeValue(double v) {
  DEBUG("[Workers Pipe] Setting value to: " << v);
  // Search the value
  std::vector<double> allocation;
  for (auto av : _allowedValues) {
    double sum = 0;
    for (auto d : av) {
      sum += d;
    }
    if (sum == v) {
      allocation = av;
      break;
    }
  }
  assert(allocation.size() == _farms.size());
  for (size_t i = 0; i < allocation.size(); i++) {
    DEBUG("[Workers Pipe] Setting num workers to " << allocation[i]
                                                   << " for farm " << i);
    _farms[i]->setRealValue(allocation[i]);
    DEBUG("[Workers Pipe] Num workers for farm " << i << " changed.");
  }
}

std::vector<AdaptiveNode *> KnobVirtualCoresPipe::getActiveWorkers() const {
  std::vector<AdaptiveNode *> r;
  for (auto f : _farms) {
    auto tmp = f->getActiveWorkers();
    r.insert(r.end(), tmp.begin(), tmp.end());
  }
  return r;
}

KnobHyperThreading::KnobHyperThreading(Parameters p, size_t hmp, uint cpuId) {
  vector<PhysicalCore*> physical =
    p.mammut.getInstanceTopology()->getCpu(cpuId)->getPhysicalCores();
  size_t maxHtLevel = physical[0]->getVirtualCores().size();
  for (size_t i = 0; i < maxHtLevel; i++) {
    _knobValues.push_back(i + 1);
  }
  _realValue = 1;
}

void KnobHyperThreading::changeValue(double v) {
  ;
}

KnobMapping::KnobMapping(const Parameters &p, const KnobVirtualCores &knobCores,
                         const KnobHyperThreading &knobHyperThreading, size_t hmp,
                         uint cpuId)
    : _p(p), _knobCores(knobCores), _knobHyperThreading(knobHyperThreading),
      _hmp(hmp), _cpuId(cpuId),
      _topologyHandler(p.mammut.getInstanceTopology()) {
  for (size_t i = 0; i < MAPPING_TYPE_NUM; i++) {
    _knobValues.push_back((MappingType) i);
  }
  _realValue =
      MAPPING_TYPE_LINEAR; // This is just for initialization, is not locked.
}

template <>
char const *enumStrings<MappingType>::data[] = {"LINEAR", "INTERLEAVED", "NUM"};

void KnobMapping::changeValue(double v) {
  DEBUG("[Mapping] Changing real value to: "
        << enumToString<MappingType>((MappingType) v));

  // Its length will be equal to _knobCores.getRealValue()
  vector<VirtualCoreId> vcOrder;
  switch ((MappingType) v) {
  case MAPPING_TYPE_LINEAR: {
    vcOrder = computeVcOrderLinear();
  } break;
  case MAPPING_TYPE_INTERLEAVED: {
    vcOrder = computeVcOrderInterleaved();
  } break;
  default: {
    throw runtime_error("KnobMapping: Mapping type still not supported.");
  } break;
  }

  /** Performs mapping. **/
  _activeVirtualCores.clear();
  for (auto id : vcOrder) {
    _activeVirtualCores.push_back(_topologyHandler->getVirtualCore(id));
  }
  move(vcOrder);

  /** Updates unused virtual cores. **/
  _unusedVirtualCores.clear();
  vector<VirtualCore *> allVcs;
  if (_hmp == 1) {
    allVcs = _topologyHandler->getVirtualCores();
  } else {
    allVcs = _topologyHandler->getCpus()[_cpuId]->getVirtualCores();
  }
  for (size_t i = 0; i < allVcs.size(); i++) {
    VirtualCore *vc = allVcs.at(i);
    if (!contains(_activeVirtualCores, vc)) {
      _unusedVirtualCores.push_back(vc);
    }
  }

  DEBUG("[Mapping] Active VCs: " << _activeVirtualCores);
  DEBUG("[Mapping] Unused VCs: " << _unusedVirtualCores);
}

void KnobMapping::setAllowedCores(
    std::vector<mammut::topology::VirtualCore *> vc) {
  _allowedVirtualCores = vc;
}

std::vector<mammut::topology::VirtualCore *>
KnobMapping::getAllowedCores() const {
  return _allowedVirtualCores;
}

bool KnobMapping::isAllowed(mammut::topology::VirtualCore *v) const {
  if (!_allowedVirtualCores.size()) {
    return true;
  } else {
    return utils::contains(_allowedVirtualCores, v);
  }
}

const vector<VirtualCore *> &KnobMapping::getActiveVirtualCores() const {
  return _activeVirtualCores;
}

const vector<VirtualCore *> &KnobMapping::getUnusedVirtualCores() const {
  return _unusedVirtualCores;
}

size_t KnobMapping::getNumVirtualCores() {
  return _knobCores.getRealValue();
}

std::vector<mammut::topology::VirtualCoreId>
KnobMapping::computeVcOrderLinear() {
  /*
   * Generates a vector of virtual cores to be used for linear
   * mapping. It contains first one virtual core per physical
   * core (virtual cores on the same CPU are consecutive).
   * Then, the other groups of virtual cores follow.
   */
  vector<VirtualCoreId> vcOrder;
  size_t virtualPerPhysical = _knobHyperThreading.getRealValue();

  vector<Cpu *> cpus = _topologyHandler->getCpus();
  for (size_t k = 0; k < virtualPerPhysical; k++) {
    for (size_t i = 0; i < cpus.size(); i++) {
      if (_hmp > 1 && i != _cpuId) {
        continue;
      }
      vector<PhysicalCore *> phyCores = cpus.at(i)->getPhysicalCores();
      for (size_t j = 0; j < phyCores.size(); j++) {
        vector<VirtualCore *> virtCores = phyCores.at(j)->getVirtualCores();
        if ((!_p.isolateManager || virtCores.at(k)->getVirtualCoreId() !=
                                       NORNIR_MANAGER_VIRTUAL_CORE) &&
            isAllowed(virtCores.at(k))) {
          vcOrder.push_back(virtCores.at(k)->getVirtualCoreId());
          if (vcOrder.size() == getNumVirtualCores()) {
            return vcOrder;
          }
        }
      }
    }
  }
  return vcOrder;
}

std::vector<mammut::topology::VirtualCoreId>
KnobMapping::computeVcOrderInterleaved() {
  /*
   * Generates a vector of virtual cores to be used for interleaved
   * mapping.
   */
  vector<VirtualCoreId> vcOrder;
  size_t virtualPerPhysical = _knobHyperThreading.getRealValue();
  vector<Cpu *> cpus = _topologyHandler->getCpus();
  size_t physicalPerCpu = cpus[0]->getPhysicalCores().size();
  for (size_t k = 0; k < virtualPerPhysical; k++) {
    for (size_t j = 0; j < physicalPerCpu; j++) {
      for (size_t i = 0; i < cpus.size(); i++) {
        if (_hmp > 1 && i != _cpuId) {
          continue;
        }
        vector<PhysicalCore *> phyCores = cpus.at(i)->getPhysicalCores();
        vector<VirtualCore *> virtCores = phyCores.at(j)->getVirtualCores();
        if ((!_p.isolateManager || virtCores.at(k)->getVirtualCoreId() !=
                                       NORNIR_MANAGER_VIRTUAL_CORE) &&
            isAllowed(virtCores.at(k))) {
          vcOrder.push_back(virtCores.at(k)->getVirtualCoreId());
          if (vcOrder.size() == getNumVirtualCores()) {
            return vcOrder;
          }
        }
      }
    }
  }
  return vcOrder;
}

KnobMappingExternal::KnobMappingExternal(
    const Parameters &p, const KnobVirtualCores &knobCores,
    const KnobHyperThreading &knobHyperThreading, size_t hmp, uint cpuId)
    : KnobMapping(p, knobCores, knobHyperThreading, hmp, cpuId),
      _processHandler(NULL) {
  ;
}

void KnobMappingExternal::setPid(pid_t pid) {
  if (_processHandler) {
    _p.mammut.getInstanceTask()->releaseProcessHandler(_processHandler);
  }
  setProcessHandler(_p.mammut.getInstanceTask()->getProcessHandler(pid));
}

void KnobMappingExternal::setProcessHandler(
    task::ProcessHandler *processHandler) {
    _processHandler = processHandler;
    auto threads = _processHandler->getActiveThreadsIdentifiers();
    _lastThreads.clear();
    for(auto t : threads){
        _lastThreads.push_back(_processHandler->getThreadHandler(t));
    }
}

void KnobMappingExternal::move(
    const std::vector<mammut::topology::VirtualCoreId> &vcOrder) {
  if (_processHandler) {
    if(_p.fixedPinning){     
      std::vector<mammut::task::TaskId> threadsSubset;
      size_t goodId = 0;
      for(size_t i = 0; i < _lastThreads.size(); i++){
	    double usage = 0.0;
	    auto th = _lastThreads[i];
	    bool active = th->getCoreUsage(usage);
	    DEBUG("Thread " << th->getId() << " active " << active << " usage: " << usage);
	    if(active && usage > 5.0){
	      DEBUG("mod " << goodId % _hmp << " right " << (_hmp - _cpuId - 1));
	      if(goodId % _hmp == (_hmp - _cpuId - 1)){ // _hmp - _cpuId - 1 rather than _hmp, so that we start assignign from the last (fastest) HMP domain
	        threadsSubset.push_back(th->getId());
	      }
	      ++goodId;
	    }
      }

	  for(auto t : _lastThreads){
        _processHandler->releaseThreadHandler(t);
      }
      auto threads = _processHandler->getActiveThreadsIdentifiers();
	  _lastThreads.clear();
      for(auto t : threads){
        _lastThreads.push_back(_processHandler->getThreadHandler(t));
	  }

      size_t i = 0;
      for(auto t : threadsSubset){
        auto th = _processHandler->getThreadHandler(t);
        th->move(vcOrder[i]);
        DEBUG("CpuID: " << _cpuId << " Thread " << t << " pinned to core " << vcOrder[i]);
        i = (i + 1) % vcOrder.size();
        _processHandler->releaseThreadHandler(th);
      }
    }else{
      if (_cpuId == 0) {
        _processHandler->move(vcOrder);
      } else {
        if (_hmp == 1) {
          throw std::runtime_error("move called on cpuId != 0 but no HMP.");
        }
        std::vector<mammut::topology::VirtualCoreId> old;
        _processHandler->getVirtualCoreIds(old);
        old.insert(old.end(), vcOrder.begin(), vcOrder.end());
        _processHandler->move(old);
      }
    }
  } else {
    throw std::runtime_error("setPid or setProcessHandler must be called "
                             "before using KnobMappingExternal.");
  }
}

KnobMappingFarm::KnobMappingFarm(const Parameters &p,
                                 const KnobVirtualCoresFarm &knobCores,
                                 const KnobHyperThreading &knobHyperThreading,
                                 AdaptiveNode *emitter, AdaptiveNode *collector,
                                 size_t hmp, uint cpuId)
    : KnobMapping(p, knobCores, knobHyperThreading, hmp, cpuId),
      _emitter(emitter), _collector(collector) {
  if (hmp > 1) {
    throw std::runtime_error("HMP not supported for KnobMappingFarm.");
  }
}

size_t KnobMappingFarm::getNumVirtualCores() {
  size_t v = _knobCores.getRealValue();
  if (_emitter)
    ++v;
  if (_collector)
    ++v;
  return v;
}

void KnobMappingFarm::move(
    const std::vector<mammut::topology::VirtualCoreId> &vcOrder) {
  vector<AdaptiveNode *> workers =
      ((KnobVirtualCoresFarm *) &_knobCores)->getActiveWorkers();
  size_t numServiceNodes = 0;
  if (_emitter)
    ++numServiceNodes;
  if (_collector)
    ++numServiceNodes;
  if (workers.size() + numServiceNodes <= vcOrder.size()) {
    size_t nextIndex = 0;
    size_t emitterIndex = 0, collectorIndex = 0;
    if (_emitter) {
      _emitter->move(vcOrder[nextIndex]);
      emitterIndex = nextIndex;
      nextIndex = (nextIndex + 1) % vcOrder.size();
    }

    if (_collector) {
      _collector->move(vcOrder[nextIndex]);
      collectorIndex = nextIndex;
      nextIndex = (nextIndex + 1) % vcOrder.size();
    }

    for (size_t i = 0; i < workers.size(); i++) {
      // TODO Decide to do or not to do this check it with a parameter
      if (nextIndex == emitterIndex || nextIndex == collectorIndex) {
        nextIndex = (nextIndex + 1) % vcOrder.size();
      }
      workers[i]->move(vcOrder[nextIndex]);
      nextIndex = (nextIndex + 1) % vcOrder.size();
    }
  } else {
    _p.mammut.getInstanceTask()->getProcessHandler(getpid())->move(vcOrder);
  }
}

KnobFrequency::KnobFrequency(Parameters p, const KnobMapping &knobMapping,
                             size_t hmp, uint cpuId)
    : _p(p), _knobMapping(knobMapping),
      _frequencyHandler(_p.mammut.getInstanceCpuFreq()),
      _topologyHandler(_p.mammut.getInstanceTopology()) {
  _frequencyHandler->removeTurboFrequencies();
  std::vector<mammut::cpufreq::Frequency> availableFrequencies;
  availableFrequencies =
      _frequencyHandler->getDomains()[cpuId]->getAvailableFrequencies();
  for (Frequency f : availableFrequencies) {
    if((_p.minClockFrequency != -1 && f < _p.minClockFrequency) ||
       (_p.maxClockFrequency != -1 && f > _p.maxClockFrequency)){
      continue;
    }else{
      _knobValues.push_back(f);
    }
  }
  if (_p.knobFrequencyEnabled) {
    if (_knobValues.empty()) {
      throw std::runtime_error("Frequencies not available. Please set "
                               "knobFrequencyEnabled to false.");
    } else {
      std::vector<mammut::cpufreq::Domain *> scalableDomains;
      if (hmp > 1) {
        scalableDomains.push_back(_frequencyHandler->getDomains()[cpuId]);
      } else {
        scalableDomains = _frequencyHandler->getDomains();
      }
      _changeWithOnDemand = false;
      for (Domain *currentDomain : scalableDomains) {
        if (!currentDomain->setGovernor(GOVERNOR_USERSPACE)) {
          _changeWithOnDemand = true;
          // If userspace is not available, set ondemand with min/max equal to the selected frequency
          if (!currentDomain->setGovernor(GOVERNOR_ONDEMAND)) {
            throw runtime_error("KnobFrequency: Impossible to set the frequency manually.");
          }
        }
        // I set the minimum and maximum frequencies to the min and max
        // of this system.
        if (!currentDomain->setGovernorBounds(_knobValues.front(), _knobValues.back())) {
          throw runtime_error("KnobFrequency: Impossible "
                              "to set the governor bounds.");
        }
      }
    }
    _realValue = _knobValues.front();
  } else {
    _realValue = 1;
    _knobValues.push_back(_realValue);
  }
}

KnobFrequency::~KnobFrequency() {
  for (Domain *currentDomain : _frequencyHandler->getDomains()) {
    currentDomain->reinsertTurboFrequencies();
  }
}

void KnobFrequency::changeValue(double v) {
  if (!_p.knobFrequencyEnabled) {
    return;
  }
  DEBUG("[Frequency] Changing real value to: " << v);
  std::vector<mammut::cpufreq::Domain *> scalableDomains;
  scalableDomains =
      _frequencyHandler->getDomains(_knobMapping.getActiveVirtualCores());
  for (size_t i = 0; i < scalableDomains.size(); i++) {
    Domain *currentDomain = scalableDomains.at(i);
    setFrequency(currentDomain, v);
  }
  DEBUG("[Frequency] Frequency changed for domains: " << scalableDomains);
  applyUnusedVCStrategy(v);
  DEBUG("[Frequency] Active VC: " << _knobMapping.getActiveVirtualCores());
  DEBUG("[Frequency] Unused VC: " << _knobMapping.getUnusedVirtualCores());
}

void KnobFrequency::setFrequency(Domain *domain, Frequency frequency){
  if(_changeWithOnDemand){
    if (!domain->setGovernorBounds(frequency, frequency)) {
      throw runtime_error("KnobFrequency: Impossible to set the specified frequency.");
    }
  }else{
    if (!domain->setFrequencyUserspace((uint) frequency)) {      
      throw runtime_error("KnobFrequency: Impossible to set the specified frequency.");
    }
  }
}

void KnobFrequency::applyUnusedVCStrategySame(
    const vector<VirtualCore *> &unusedVc, Frequency v) {
  vector<Domain *> unusedDomains =
      _frequencyHandler->getDomainsComplete(unusedVc);
  DEBUG("[Frequency] " << unusedDomains.size() << " unused domains.");
  for (size_t i = 0; i < unusedDomains.size(); i++) {
    Domain *domain = unusedDomains.at(i);
    DEBUG("[Frequency] Setting unused domain " << domain->getId()
                                               << " to: " << v);
    setFrequency(domain, v);
  }
}

void KnobFrequency::applyUnusedVCStrategyOff(
    const vector<VirtualCore *> &unusedVc) {
  for (size_t i = 0; i < unusedVc.size(); i++) {
    VirtualCore *vc = unusedVc.at(i);
    if (vc->isHotPluggable() && vc->isHotPlugged()) {
      vc->hotUnplug();
    }
  }
}

void KnobFrequency::applyUnusedVCStrategyLowestFreq(
    const vector<VirtualCore *> &unusedVc) {
  vector<Domain *> unusedDomains =
      _frequencyHandler->getDomainsComplete(unusedVc);
  for (size_t i = 0; i < unusedDomains.size(); i++) {
    Domain *domain = unusedDomains.at(i);
    if(_changeWithOnDemand){
      if (!domain->setGovernorBounds(_knobValues.front(), _knobValues.front())) {
        throw runtime_error("KnobFrequency: Impossible to "
                            "set lowest frequency for unused "
                            "virtual cores.");
      }
    }else{
      if (!domain->setGovernor(GOVERNOR_USERSPACE) ||
          (_p.knobFrequencyEnabled && !domain->setLowestFrequencyUserspace())) {
        throw runtime_error("KnobFrequency: Impossible to "
                            "set lowest frequency for unused "
                            "virtual cores.");
      }
    }
  }
}

void KnobFrequency::applyUnusedVCStrategy(Frequency v) {
  /**
   * OFF 'includes' LOWEST_FREQUENCY. i.e. If we shutdown all the
   * virtual cores on a domain, we can also lower its frequency to
   * the minimum.
   */
  vector<VirtualCore *> unusedVc = _knobMapping.getUnusedVirtualCores();

  if (_p.strategyUnusedVirtualCores == STRATEGY_UNUSED_VC_SAME) {
    applyUnusedVCStrategySame(unusedVc, v);
  } else {
    vector<VirtualCore *> virtualCores;
    if (_p.strategyUnusedVirtualCores != STRATEGY_UNUSED_VC_NONE) {
      insertToEnd(unusedVc, virtualCores);
    }
    applyUnusedVCStrategyLowestFreq(virtualCores);

    virtualCores.clear();
    if (_p.strategyUnusedVirtualCores == STRATEGY_UNUSED_VC_OFF) {
      insertToEnd(unusedVc, virtualCores);
    }
    applyUnusedVCStrategyOff(virtualCores);
  }
}

KnobClkMod::KnobClkMod(Parameters p, const KnobMapping &knobMapping, size_t hmp,
                       uint cpuId)
    : _knobMapping(knobMapping) {
  for (double d : p.mammut.getInstanceTopology()
                      ->getCpus()[cpuId]
                      ->getClockModulationValues()) {
    _knobValues.push_back(d);
  }
  _realValue = 100.0;
}

void KnobClkMod::changeValue(double v) {
  for (VirtualCore *vc : _knobMapping.getActiveVirtualCores()) {
    vc->setClockModulation(v);
  }
}

KnobClkModEmulated::KnobClkModEmulated(Parameters p)
    : _p(p), _processHandler(NULL) {
  _realValue = 100.0;
  _knobValues.clear();
  for (double i = _p.clockModulationMin; i < 100;
       i += _p.clockModulationResolution) {
    _knobValues.push_back(i);
  }
  _knobValues.push_back(100.0);
}

void KnobClkModEmulated::setPid(pid_t pid) {
  setProcessHandler(_p.mammut.getInstanceTask()->getProcessHandler(pid));
}

void KnobClkModEmulated::setProcessHandler(
    task::ProcessHandler *processHandler) {
  _processHandler = processHandler;
  _p.mammut.getInstanceTask()->setThrottlingInterval(
      _p.clockModulationInterval);
}

void KnobClkModEmulated::changeValue(double v) {
  if (v != _realValue) {
    if (_processHandler) {
      _processHandler->throttle(v);
    } else {
      throw std::runtime_error("setPid or setProcessHandler must be called "
                               "before using KnobClkModEmulated.");
    }
  }
}

KnobPforChunk::KnobPforChunk(Parameters p):_p(p){
  ;
}

void KnobPforChunk::setChunkPointer(long int* chunkPointer){
  _chunkPointer = chunkPointer;
}

void KnobPforChunk::setPossibleValues(long long int start, long long end, long long step, size_t numThreads){
  _knobValues.clear();
  double numElements = std::ceil((end - start)/(double) step);
  double lastDenominator = numThreads;
  double v = 0;
  while((v = std::ceil(numElements / lastDenominator)) > 1){
    _knobValues.push_back(v);
    lastDenominator *= 2;
  }
  _knobValues.push_back(1);
  std::reverse(_knobValues.begin(), _knobValues.end());
#ifdef DEBUG_KNOB
  for(auto v : _knobValues){
    DEBUG("[PforChunk] Value " << v);
  }
#endif
  _realValue = _knobValues[0];
  *_chunkPointer = _realValue;
}

void KnobPforChunk::changeValue(double v){
  *_chunkPointer = v;
}


} // namespace nornir

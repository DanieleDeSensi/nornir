/*
 * manager-multi.cpp
 *
 * Created on: 22/07/2016
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

#include "manager-multi.hpp"
#include "external/mammut/mammut/mammut.hpp"

#include <map>
#include <signal.h>
#include <sys/types.h>

using namespace std;
using namespace mammut;
using namespace mammut::topology;
using namespace mammut::cpufreq;

/**
 * It works only if we have htlevel == 1.
 */

#undef DEBUG
#undef DEBUGB

#ifdef DEBUG_MANAGER_MULTI
#define DEBUG(x) do { cerr << "[ManagerMulti] " << x << endl; } while (0)
#define DEBUGB(x) do {x;} while(0)
#else
#define DEBUG(x)
#define DEBUGB(x)
#endif

namespace nornir{

class SubmittedManager{
public:
    Manager* manager;
    double minPerf;
    SubmittedManager(Manager* m, double minP):manager(m), minPerf(minP){;}
};

#define MAX_POWER_VIOLATION_SECONDS 10

static inline const std::vector<Allocation>& getAllocations(const std::pair<Manager* const, ManagerData>& it){
    return it.second.allocations;
}

static inline const std::vector<mammut::topology::VirtualCoreId>& getAllocatedCores(const std::pair<Manager* const, ManagerData>& it){
    return it.second.allocatedCores;
}

static inline const double& getMinPerf(const std::pair<Manager* const, ManagerData>& it){
    return it.second.minPerf;
}

static inline Manager* const getManager(const std::pair<Manager* const, ManagerData>& it){
    return it.first;
}

static inline const std::vector<Allocation>& getAllocations(const std::map<Manager*, ManagerData>::const_iterator& it){
    return getAllocations(*it);
}

static inline const std::vector<mammut::topology::VirtualCoreId>& getAllocatedCores(const std::map<Manager*, ManagerData>::const_iterator& it){
    return getAllocatedCores(*it);
}

static inline const double& getMinPerf(const std::map<Manager*, ManagerData>::const_iterator& it){
    return getMinPerf(*it);
}

static inline Manager* const getManager(const std::map<Manager*, ManagerData>::const_iterator& it){
    return getManager(*it);
}

static inline const KnobsValues& getKnobs(const Allocation& a){return a.first;}
static inline const double& getPrediction(const Allocation& a){return a.second;}

static inline const KnobsValues& getKnobs(const std::pair<Manager* const, ManagerData>& it, size_t pos){return getKnobs(getAllocations(it).at(pos));}
static inline const double& getPrediction(const std::pair<Manager* const, ManagerData>& it, size_t pos){return getPrediction(getAllocations(it).at(pos));}

ManagerMulti::ManagerMulti(ManagerMultiConfiguration configuration):
        _configuration(configuration), _qIn(10), _qOut(10),
        _power(new MovingAverageSimple<double>(MAX_POWER_VIOLATION_SECONDS)){
    if(_configuration.powerCap < 0){
        throw std::runtime_error("[ManagerMulti]: powerCap must be >= 0.");
    }
    _qIn.init();
    _qOut.init();
    _topology = _m.getInstanceTopology();

    if(configuration.useVirtualCores){
        vector<VirtualCore*> allCores = _topology->getVirtualCores();
        for(size_t i = 0; i < allCores.size(); i++){
            _allCores.push_back(allCores.at(i)->getVirtualCoreId());
        }
    }else{
        vector<PhysicalCore*> allCores = _topology->getPhysicalCores();
        for(size_t i = 0; i < allCores.size(); i++){
            _allCores.push_back(allCores.at(i)->getVirtualCore()->getVirtualCoreId());
        }
    }
}

ManagerMulti::~ManagerMulti(){
    delete _power;
}

void ManagerMulti::addManager(Manager* m, double minPerformanceRequired){
    // Send the new manager
    checkManagerSupported(m);
    SubmittedManager* sm = new SubmittedManager(m, minPerformanceRequired);
    while(!_qIn.push(sm)){;}
}

Manager* ManagerMulti::getTerminatedManager(){
    Manager *m;
    if(_qOut.pop((void**) &m)){
        return m;
    }else{
        return NULL;
    }
}

void ManagerMulti::allowCalibration(Manager* m){
    DEBUG("Manager (" << m << ") wants to calibrate.");
    _managerData[m].allocatedCores.clear();
    vector<VirtualCoreId> cores = getAvailableCores();
    allowCores(m, cores);
    DEBUG("Manager (" << m << ") can calibrate on: " << cores);
    m->_selector->acceptViolations();
    m->_selector->allowCalibration();
}

void ManagerMulti::waitForCalibration(Manager* m){
    // We need to check the total calibration time
    // to be sure that the calibration at least
    // started the first time.
    while(m->_selector->isCalibrating() ||
          !m->_selector->getTotalCalibrationTime()){
        ;
    }
    m->_selector->ignoreViolations();
    DEBUG("Manager (" << m << ") terminated its calibration.");
}

void ManagerMulti::inhibitAll(Manager* except){
    for(const auto& it : _managerData){
        Manager* const currentManager = getManager(it);
        if(currentManager != except){
            currentManager->inhibit();
        }
    }
    DEBUG("Inhibition done.");
}

void ManagerMulti::disinhibitAll(){
    for(const auto& it : _managerData){
        getManager(it)->disinhibit();
    }
    DEBUG("All managers disinhibited.");
}

void ManagerMulti::shrinkAll(Manager* except){
    for(auto& it : _managerData){
        Manager* const currentManager = getManager(it);
        std::vector<mammut::topology::VirtualCoreId>& aCores = it.second.allocatedCores;
        if(currentManager != except){
            currentManager->shrink(_configuration.shrink);
            switch(_configuration.shrink){
                case CALIBRATION_SHRINK_AGGREGATE:{
                    aCores.clear();
                    // TODO: At the moment I put everything on the same core.
                    // It would maybe be better
                    // to have each different application on a different core.
                    aCores.push_back(_topology->getVirtualCores().back()->getVirtualCoreId());
                }break;
                case CALIBRATION_SHRINK_PAUSE:{
                    aCores.clear();
                }break;
                default:{
                    ;
                }break;
            }
        }
    }
    DEBUG("Shrinking done.");
}

vector<VirtualCoreId> ManagerMulti::getAvailableCores() const{
    vector<VirtualCoreId> availableCores;
    for(mammut::topology::VirtualCoreId vid : _allCores){
        bool alreadyAllocated = false;
        for(const auto& it : _managerData){
            if(utils::contains(getAllocatedCores(it), vid)){
                alreadyAllocated = true;
                DEBUG("Core " << vid << " already used by manager " << getManager(it));
                break;
            }
        }
        if(!alreadyAllocated){
            availableCores.push_back(vid);
        }
    }
    return availableCores;
}

static inline double getMaxPerformance(const std::map<KnobsValues, double>& primaryValues){
    double maxPerformance = 0.0;
    for(const auto& it : primaryValues){
        const double& predictedPerformance = getPrediction(it);
        if(predictedPerformance > maxPerformance){
            maxPerformance = predictedPerformance;
        }
    }
    return maxPerformance;
}

void ManagerMulti::updateAllocations(Manager* const m){
    const std::map<KnobsValues, double>& primaryValues = dynamic_cast<SelectorPredictive*>(m->_selector)->getPrimaryPredictions();
    const std::map<KnobsValues, double>& secondaryValues = dynamic_cast<SelectorPredictive*>(m->_selector)->getSecondaryPredictions();
    double primaryBound = m->_p.requirements.throughput;
    ManagerData& md = _managerData[m];
    double referencePerformance;
    if(isMinMaxRequirement(m->_p.requirements.throughput)){
        referencePerformance = getMaxPerformance(primaryValues);
    }else{
        referencePerformance = m->_p.requirements.throughput;
    }
    md.minPerf = (md.minPerfReqPerc / 100.0) * referencePerformance;

    std::multimap<double, KnobsValues> unfeasible, sortedSecondary;
    for(auto it : primaryValues){
        const double& prediction = getPrediction(it);
        if(prediction >= primaryBound){
            // Insert the corresponding entry in secondaryValues
            // since is a map, they will be kept sorted from the lower power consuming
            // to the higher power consuming.
            sortedSecondary.insert(AllocationFlip(secondaryValues.at(it.first), it.first));
        }else{
            // Insert unfeasible solutions according to their relative performance in
            // percentage (from lowest to highest).
            double relativePerf = (prediction / referencePerformance) * 100;
            unfeasible.insert(AllocationFlip(relativePerf, it.first));
        }
    }

    md.allocations.clear();
    // First insert the solutions that satisfies the primary bound (sorted from
    // the best (lowest) to the worst (highest) secondary value).
    for(const auto& it : sortedSecondary){
        // Here all these solutions satisfy the constraint, and we only want to
        // order according to the secondary value. So we are not interested on
        // the relative performance. So we put 100%.
        double relativePerf = 100;
        md.allocations.push_back(Allocation(it.second, relativePerf));
    }
    // Then we insert the unfeasible solutions (i.e. that violate the primary bound)
    // sorted from the most performing to the least performing.
    // We scan on the reverse direction since they are ordered from the least
    // to the most performing.
    for(auto it = unfeasible.rbegin(); it != unfeasible.rend(); it++){
        double relativePerf = (it->first / referencePerformance) * 100;
        md.allocations.push_back(Allocation(it->second, relativePerf));
    }
}

void ManagerMulti::updateAllocations(){
    for(const auto& it : _managerData){
        updateAllocations(getManager(it));
    }
    std::vector<AllocationIndexes> values;
    AllocationIndexes accum;
    for(auto& it : _managerData){
        AllocationIndexes tmp;
        for(size_t j = 0; j < getAllocations(it).size(); j++){
            tmp.push_back(j);
        }
        values.push_back(tmp);
    }
    _allocationsCombinations.clear();
    combinations(values, 0, accum);
}

void ManagerMulti::combinations(std::vector<AllocationIndexes> array,
                                size_t i,
                                AllocationIndexes accum){
    if(i == array.size()){
        _allocationsCombinations.push_back(accum);
    }else{
        AllocationIndexes row = array.at(i);
        for(size_t j = 0; j < row.size(); ++j){
            AllocationIndexes tmp(accum);
            tmp.push_back(row[j]);
            combinations(array, i+1, tmp);
        }
    }
}

double ManagerMulti::getQuality(const AllocationIndexes& indexes) const{
    switch(_configuration.qualityEstimation){
        case QUALITY_ESTIMATION_PERFORMANCE:{
            // For each manager, we check how much performance degradation (in percentage)
            // we would get. Then, we return the average between these degradations.
            double avgRelPerf = 0;
            size_t pos = 0;
            for(const auto& it : _managerData){
                avgRelPerf += getPrediction(it, indexes.at(pos));
                ++pos;
            }
            avgRelPerf /= _managerData.size();
            return avgRelPerf;
        }break;
        case QUALITY_ESTIMATION_PREFERENCE:{
            // At the moment we are considering as quality metric the sum of
            // the positions of the chosen configurations in the corresponding priority lists
            // (the lower the better).
            // If we have 3 managers and we select the second choice for manager A, the
            // tenth choice for manager B and the fifth choice for manager C, then this
            // specific allocation will have a quality of 2+10+5 = 17. Accordingly, in
            // this specific case the best allocation would be the one with weight 3, i.e.
            // the one in which to each manager we assign its first (most preferred) choice.
            // We could also consider other metrics (e.g. average).
            double r = 0;
            for(size_t x : indexes){r += x;}
            // Since it is the lowest the better but we want the higher the better.
            return std::numeric_limits<double>::max() - r;
        }break;
        default:{
            throw std::runtime_error("Unsupported quality estimation algorithm.");
        }
    }
}

double ManagerMulti::estimatePower(const AllocationIndexes& indexes) const{
    return 0; //TODO
}

void ManagerMulti::updateModels(){
    for(const auto& it : _managerData){
        Manager* const m = getManager(it);
        m->updateModelsInterference();
        m->waitModelsInterferenceUpdate();
    }
}

void ManagerMulti::calibrate(Manager* m, bool start){
    inhibitAll(m);
    shrinkAll(m);
    allowCalibration(m);
    if(start){
        m->start(); DEBUG("Manager started.");
    }
    waitForCalibration(m);
    applyNewAllocation();
    m->inhibit(); DEBUG(m << " inhibited.");
    for(const auto& it : _managerData){
        Manager* const currentMan = getManager(it);
        if(currentMan != m){
            currentMan->stretch(_configuration.shrink);
        }
    }
    DEBUG("Everyone stretched.");
    sleep(1);
    DEBUG("Ready to update the models.");
    updateModels();
    DEBUG("Models updated.");
    disinhibitAll();
}

// We want the elements sorted from highest to lowest quality.
bool validAllocationComp(const std::pair<double const, const AllocationIndexes*>& i,
                         const std::pair<double const, const AllocationIndexes*>& j){
    return i.first > j.first;
}

void ManagerMulti::getValidAllocations(ValidAllocations& validAllocations) const{
    validAllocations.clear();
    validAllocations.reserve(_allocationsCombinations.size());
    for(size_t i = 0; i < _allocationsCombinations.size(); i++){
        const AllocationIndexes* indexes = &(_allocationsCombinations.at(i));
        double quality = -1;
        if(isValidAllocation(*indexes)){
            quality = getQuality(*indexes);
        }
        //DEBUG("Allocation " << *indexes << " has quality " << quality);
        validAllocations.emplace_back(quality, indexes);
    }
    std::sort(validAllocations.begin(), validAllocations.end(), validAllocationComp);
}

AllocationIndexes ManagerMulti::findBestAllocation(){
    DEBUG("Searching for best allocation...");
    ValidAllocations validAllocations;
    updateAllocations();
    getValidAllocations(validAllocations);
    // First try to find a solution that doesn't violate
    // any additional requirement. If does not exists,
    // just return the one with maximum quality.
    for(auto indexes : validAllocations){
        size_t pos = 0;
        bool feasible = true;
        for(auto& it : _managerData){
            // Check that we still satisfy the additional requirement
            // set on the global manager.
            if(getPrediction(it, indexes.second->at(pos)) < getMinPerf(it)){
                feasible = false;
                break;
            }
            ++pos;
        }
        // Since we are iterating from the best to the worst,
        // as soon as we find a feasible one we return.
        if(feasible){return *(indexes.second);}
    }

    DEBUG(validAllocations.size() << " allocations evaluated. All of them violates the additional performance requirements.");
    // If we are here, there are no solutions that satisfy the
    // additional performance constraints, so we return the
    // one with maximum quality.
    return *(validAllocations.begin()->second);
}

void ManagerMulti::applyNewAllocation(){
    AllocationIndexes alloc = findBestAllocation();
    size_t pos = 0, nextCoreId = 0;
    DEBUG("Best allocation found: " << alloc);
    for(const auto& it : _managerData){
        Manager* const m = getManager(it);
        const KnobsValues real = m->_configuration->getRealValues(getKnobs(it, alloc.at(pos)));
        if(!m->running()){++pos; continue;}
        DEBUG("Allocation: " << m << " " << real << " " << dynamic_cast<SelectorPredictive*>(m->_selector)->getPrimaryPrediction(real) << " " << dynamic_cast<SelectorPredictive*>(m->_selector)->getSecondaryPrediction(real));
        size_t numCores = real[KNOB_VIRTUAL_CORES] + m->_configuration->getNumServiceNodes();
        vector<VirtualCoreId> cores;
        cores.reserve(numCores);
        for(size_t i = 0; i < numCores; i++){
            cores.push_back(_allCores.at((i + nextCoreId) % _allCores.size()));
        }
        allowCores(m, cores);
        //man->_selector->forceConfiguration(real);
        m->act(real, true);
        dynamic_cast<SelectorPredictive*>(m->_selector)->updatePredictions(real);
        ++pos;
        nextCoreId += numCores;
    }
    DEBUG("Best allocation applied.");
}

void ManagerMulti::checkManagerSupported(Manager* m){
    // It is meaningless to control power for individual applications if
    // other applications are running on the system and want to be
    // controlled as well.
    assert(!isPrimaryRequirement(m->_p.requirements.powerConsumption));

    assert(m->_p.strategySelection == STRATEGY_SELECTION_LEARNING);
    assert(m->_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_AMDAHL ||
            m->_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_USL ||
            m->_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_USLP);
}

void ManagerMulti::allowCores(Manager* m, const std::vector<VirtualCoreId>& cores){
    DEBUG("Manager " << m << " uses cores: " << cores);
    _managerData[m].allocatedCores = cores;
    m->allowCores(cores);
}

bool ManagerMulti::isValidAllocation(const AllocationIndexes& indexes) const{
    if(indexes.empty()){
        throw std::runtime_error("FATAL ERROR: No indexes.");
    }
    size_t pos = 0, totalCores = 0;
    Frequency previousFreq = 0;
    bool validAllocation = true;
    for(const auto& it : _managerData){
        Manager* const currentManager = getManager(it);
        size_t allocationPosition = indexes.at(pos);
        const KnobsValues real = currentManager->_configuration->getRealValues(getKnobs(it, allocationPosition));
        size_t numCores = real[KNOB_VIRTUAL_CORES] + currentManager->_configuration->getNumServiceNodes();
        Frequency currentFreq = real[KNOB_FREQUENCY];
        // Only keep combinations on the same frequency.
        if(previousFreq && currentFreq != previousFreq){
            validAllocation = false;
            break;
        }
        previousFreq = currentFreq;
        totalCores += numCores;
        ++pos;
    }
    if(totalCores > _allCores.size() ||
       estimatePower(indexes) > _configuration.powerCap){
        validAllocation = false;
    }
    return validAllocation;
}

void ManagerMulti::run(){
    mammut::energy::Counter* joulesCounter = _m.getInstanceEnergy()->getCounter();
    double lastSampleTime = mammut::utils::getMillisecondsTime();
    double lastJoules = joulesCounter->getJoules();
    while(true){
        SubmittedManager* sm;
        if(_qIn.pop((void**) &sm)){
            ManagerData md;
            Manager* m = sm->manager;
            DEBUG("Manager (" << m << ") arrived.");
            md.minPerf = 0;
            md.minPerfReqPerc = sm->minPerf;
            md.allocatedCores = std::vector<VirtualCoreId>();
            md.allocations = std::vector<Allocation>();
            _managerData[m] = md;
            m->_selector->setCalibrationCoordination();
            calibrate(m, true);
            delete sm;
        }else{
            // Manage already present managers.
            for(auto it = _managerData.begin(); it != _managerData.end(); ){
                Manager* const m = getManager(it);
                if(!m->running()){
                    // Manager terminated
                    it = _managerData.erase(it);
                    inhibitAll();
                    updateModels();
                    if(!_managerData.empty()){
                        applyNewAllocation();
                        disinhibitAll();
                    }
                    while(!_qOut.push((void*) m)){;}
                }else{
                    // A manager requested permission to calibrate.
                    if(m->_selector->isCalibrating()){
                        calibrate(m);
                    }
                    ++it;
                }
            }
        }
        double currentSampleTime, currentJoules, currentWatts;
        currentJoules = joulesCounter->getJoules();
        currentSampleTime = mammut::utils::getMillisecondsTime();
        currentWatts = (currentJoules - lastJoules)/((currentSampleTime - lastSampleTime)/1000.0);
        lastSampleTime = currentSampleTime;
        lastJoules = currentJoules;
        _power->add(currentWatts);
        // TODO In realtà bisognerebbe controllare che non lo sforiamo
        // per un periodo consecutivo sostenuto
        if(_power->average() > _configuration.powerCap){
            DEBUG("Cap violated (" << _power->average() << ">" << _configuration.powerCap << ". "
                  "Falling back to RAP.");
            ; //TODO Fallback to RAPL
        }
        sleep(1);
    }
}

} // End namespace

/*
 * predictors.cpp
 *
 * Created on: 10/07/2015
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
 * \file predictors.cpp
 * \brief Predictors used by adaptive farm.
 */

#include "predictors.hpp"
#include "manager.hpp"

#include "external/mammut/mammut/cpufreq/cpufreq.hpp"

#undef DEBUG
#undef DEBUGB

#ifdef DEBUG_PREDICTORS
#define DEBUG(x) do { std::cerr << "[Predictors] " << x << std::endl; } while (0)
#define DEBUGB(x) do {x;} while (0)
#else
#define DEBUG(x)
#define DEBUGB(x)
#endif


namespace nornir{

using namespace mammut::cpufreq;
using namespace mammut::utils;
using namespace mammut::topology;

using namespace mlpack::regression;

using namespace arma;

using namespace std;

static double getVoltage(VoltageTable table, uint workers, Frequency frequency){
    VoltageTableKey key(workers, frequency);
    VoltageTableIterator it = table.find(key);
    if(it != table.end()){
        double v = it->second;
	if(!v){
	  v = 1; // Voltages not computed on this machine
	}
	return v;
    }else{
        throw runtime_error("Frequency and/or number of virtual cores "
                            "not found in voltage table.");
    }
}

RegressionData::RegressionData(const Parameters &p,
                               const Configuration &configuration,
                               const Smoother<MonitoredSample>* samples):
        _p(p), _configuration(configuration), _samples(samples){
    mammut::topology::Topology* t = _p.mammut.getInstanceTopology();
    _cpus = t->getCpus().size();
    _domains = _p.mammut.getInstanceCpuFreq()->getDomains().size();
    _phyCores = t->getPhysicalCores().size();
    uint virtCoresPerPhyCores = t->getPhysicalCore(0)->getVirtualCores().size();
    _phyCoresPerDomain = _p.mammut.getInstanceCpuFreq()->getDomains().at(0)->getVirtualCores().size() / virtCoresPerPhyCores;
    _phyCoresPerCpu = t->getCpus().at(0)->getPhysicalCores().size();
    if(!_phyCoresPerDomain){ // Fallback
      _phyCoresPerDomain = 1; //_phyCoresPerCpu;
    }
}

double RegressionData::getUsedPhysicalCores(double numVirtualCores){
    uint realVirtualCores = numVirtualCores + _configuration.getNumServiceNodes();
    return std::min(realVirtualCores, _phyCores);
}

void RegressionData::init(){init(_configuration.getRealValues());}

void RegressionDataServiceTime::init(const KnobsValues& values){
    _numPredictors = 0;
    uint numVirtualCores = values[KNOB_VIRTUAL_CORES];

    if(_p.knobFrequencyEnabled){
        double frequency = values[KNOB_FREQUENCY];

        _invScalFactorSpeed = (double) _minSpeed / frequency;
        ++_numPredictors;

        _invScalFactorSpeedAndCores = (double) _minSpeed /
                                      (numVirtualCores * frequency);
        ++_numPredictors;

        if(_p.knobClkModEnabled){
            double clkMod = values[KNOB_CLKMOD] / 100.0;
            _invScalFactorSpeed /= clkMod;
            _invScalFactorSpeedAndCores /= clkMod;
        }
    }else{
        throw std::runtime_error("Impossible to use Amdahl regression with no frequency.");
    }
}

RegressionDataServiceTime::RegressionDataServiceTime(const Parameters &p,
                                                     const Configuration &configuration,
                                                     const Smoother<MonitoredSample>* samples):
        RegressionData(p, configuration, samples),
        _invScalFactorSpeed(0),
        _invScalFactorSpeedAndCores(0),
        _numPredictors(0){
    _phyCores = _p.mammut.getInstanceTopology()->getPhysicalCores().size();
    std::vector<double> frequencies = _configuration.getKnob(KNOB_FREQUENCY)->getAllowedValues();
    if(!frequencies.empty()){
        _minSpeed = frequencies[0];
        if(_p.knobClkModEnabled){
            _minSpeed *= _configuration.getKnob(KNOB_CLKMOD)->getAllowedValues().front() / 100.0;
        }
    }else if(_p.knobFrequencyEnabled){
        throw std::runtime_error("Please set knobFrequencyEnabled to false.");
    }
    init(_configuration.getRealValues());
}

uint RegressionDataServiceTime::getNumPredictors() const{
    return _numPredictors;
}

void RegressionDataServiceTime::toArmaRow(size_t columnId, arma::mat& matrix) const{
    if(_p.knobFrequencyEnabled){
        size_t rowId = 0;
        matrix(rowId++, columnId) = _invScalFactorSpeed;
        matrix(rowId++, columnId) = _invScalFactorSpeedAndCores;
    }
}

static void getStaticDynamicPower(double usedPhysicalCores, Frequency frequency,
                                  MappingType mt, double clkMod,
                                  const Parameters &p, uint numCpus,
                                  uint numDomains, uint phyCoresPerCpu,
                                  uint phyCoresPerDomain,
                                  uint& unusedDomains, double& staticPower, double& dynamicPower){
    unusedDomains = 0;
    staticPower = 0;
    dynamicPower = 0;
    std::vector<uint> coresPerCpu;
    uint remainingCores = usedPhysicalCores;
    uint index = 0;

    for(size_t i = 0; i < numCpus; i++){
        coresPerCpu.push_back(0);
    }

    if(mt == MAPPING_TYPE_LINEAR){
        unusedDomains = numDomains - (usedPhysicalCores / (double) phyCoresPerDomain);
        // Simulate linear distribution
        while(remainingCores){
            ++coresPerCpu[index];
            --remainingCores;
            if(coresPerCpu[index] == phyCoresPerCpu){
                index = (index + 1) % numCpus;
            }
        }
    }else if(mt == MAPPING_TYPE_INTERLEAVED){
        // Simulate interleaved distribution
        while(remainingCores){
            ++coresPerCpu[index];
            --remainingCores;
            index = (index + 1) % numCpus;
        }
        uint usedDomains = 0;
        for(uint x : coresPerCpu){
            usedDomains += std::ceil(x / (double) phyCoresPerDomain);
        }
        unusedDomains = numDomains - usedDomains;
    }else{
        throw std::runtime_error("Unsupported mapping.");
    }

    for(uint x : coresPerCpu){
        uint fullDomains = std::floor(x / phyCoresPerDomain);
        uint coresOnSpuriousDomain = x % phyCoresPerDomain;

        // For full domains
        double voltage = 1.0;
        if(p.knobFrequencyEnabled){
            voltage = getVoltage(p.archData.voltageTable, phyCoresPerDomain, frequency);
        }
        staticPower += (voltage)*fullDomains;
        dynamicPower += (phyCoresPerDomain*frequency*voltage*voltage*clkMod)*fullDomains;

        // For spurious domain
        if(coresOnSpuriousDomain){
            if(p.knobFrequencyEnabled){
                voltage = getVoltage(p.archData.voltageTable, coresOnSpuriousDomain, frequency);
            }
            staticPower += voltage;
            dynamicPower += coresOnSpuriousDomain*frequency*voltage*voltage*clkMod;
        }
    }
}

void RegressionDataPower::init(const KnobsValues& values){
    assert(values.areReal());
    _numPredictors = 0;

    uint numVirtualCores = values[KNOB_VIRTUAL_CORES];
    double usedPhysicalCores = getUsedPhysicalCores(numVirtualCores);

    if(_p.knobFrequencyEnabled || _p.knobClkModEnabled){
        Frequency frequency = 1.0;
        double clkMod = 1.0;
        if(_p.knobFrequencyEnabled){
            frequency = values[KNOB_FREQUENCY];
        }
        if(_p.knobClkModEnabled){
            clkMod = values[KNOB_CLKMOD] / 100.0;
        }
        uint unusedDomains;
        double staticPowerUsedDomains = 0, dynamicPower = 0;

        getStaticDynamicPower(usedPhysicalCores, frequency,
                              (MappingType) values[KNOB_MAPPING],
                              clkMod, _p, _cpus,
                              _domains, _phyCoresPerCpu,
                              _phyCoresPerDomain,
                              unusedDomains, staticPowerUsedDomains, dynamicPower);

        //TODO: Potrebbe non esserci una frequenza se FREQUENCY_NO. In tal
        //caso non possiamo predirre nulla.

        //_voltagePerUsedDomains = staticPowerProp;
        //++_numPredictors;

        Frequency frequencyUnused;
        switch(_strategyUnused){
            case STRATEGY_UNUSED_VC_SAME:{
                frequencyUnused = frequency;
            }break;
            case STRATEGY_UNUSED_VC_LOWEST_FREQUENCY:{
                frequencyUnused = _lowestFrequency;
            }break;
            default:{
                throw std::runtime_error("RegressionDataPower: init: strategyUnused unsupported.");
            }
        }
        double voltage = 1.0;
        if(_p.knobFrequencyEnabled){
            voltage = getVoltage(_p.archData.voltageTable, 0, frequencyUnused);
        }
        // See TACO2016 paper, equations 6. and 7.
        _voltagePerUnusedDomains = (voltage * unusedDomains) + staticPowerUsedDomains;
        ++_numPredictors;

        _dynamicPowerModel = dynamicPower;
        ++_numPredictors;
    }else{
        /**
         * Since I do not control the frequency, this model can be very
         * inaccurate.
         */
        _dynamicPowerModel = usedPhysicalCores;
        ++_numPredictors;
    }
}

RegressionDataPower::RegressionDataPower(const Parameters &p,
                                         const Configuration &configuration,
                                         const Smoother<MonitoredSample>* samples):
        RegressionData(p, configuration, samples), _dynamicPowerModel(0),
        _voltagePerUsedDomains(0), _voltagePerUnusedDomains(0),
        _additionalContextes(0), _numPredictors(0){
    _strategyUnused = _p.strategyUnusedVirtualCores;
    std::vector<double> frequencies = _configuration.getKnob(KNOB_FREQUENCY)->getAllowedValues();
    if(!frequencies.empty()){
        _lowestFrequency = frequencies[0];
    }else if(_p.knobFrequencyEnabled){
        throw std::runtime_error("Please set knobFrequencyEnabled to false.");
    }

    init(_configuration.getRealValues());
}

uint RegressionDataPower::getNumPredictors() const{
    return _numPredictors;
}

void RegressionDataPower::toArmaRow(size_t columnId, arma::mat& matrix) const{
    size_t rowId = 0;
    matrix(rowId++, columnId) = _dynamicPowerModel;
    if(_p.knobFrequencyEnabled){
        //matrix(rowId++, columnId) = _voltagePerUsedDomains;
        matrix(rowId++, columnId) = _voltagePerUnusedDomains;
    }
}

bool RegressionDataPower::getInactivePowerPosition(size_t& pos) const{
    // ATTENTION: pos = 0 is the intercept (constant value).
    if(_p.knobFrequencyEnabled){
        pos = _numPredictors;
        return true;
    }else{
        return false;
    }
}

Predictor::Predictor(PredictorType type,
                     const Parameters &p,
                     const Configuration &configuration,
                     const Smoother<MonitoredSample>* samples):
        _type(type), _p(p), _configuration(configuration),
        _samples(samples), _modelError(0){
    ;
}

Predictor::~Predictor(){
    ;
}

double Predictor::getMaximumThroughput() const{
    if(_samples->average().loadPercentage >= MAX_RHO ||
       _samples->average().inconsistent){
        return _samples->average().throughput;
    }else{
        return _samples->average().getMaximumThroughput();
    }
}

double Predictor::getCurrentPower() const{
    return _samples->average().watts;
}

PredictorLinearRegression::PredictorLinearRegression(PredictorType type,
                                                     const Parameters& p,
                                                     const Configuration& configuration,
                                                     const Smoother<MonitoredSample>* samples):
        Predictor(type, p, configuration, samples), _preparationNeeded(true){
    switch(_type){
        case PREDICTION_THROUGHPUT:{
            _predictionInput = new RegressionDataServiceTime(p, configuration, samples);
        }break;
        case PREDICTION_POWER:{
            _predictionInput = new RegressionDataPower(p, configuration, samples);
        }break;
    }
    clear();
}

PredictorLinearRegression::~PredictorLinearRegression(){
    clear();
    delete _predictionInput;
}

void PredictorLinearRegression::clear(){
    for(obs_it iterator = _observations.begin();
               iterator != _observations.end();
               ++iterator){
        delete iterator->second.data;
    }
    _observations.clear();
    _agingVector.clear();
    _agingVector.reserve(_p.regressionAging);
    _currentAgingId = 0;
}

bool PredictorLinearRegression::readyForPredictions(){
    _predictionInput->init(_configuration.getRealValues());
    uint minPoints = std::max<uint>(_predictionInput->getNumPredictors(), 2);
    return _observations.size() >= minPoints;
}

double PredictorLinearRegression::getCurrentResponse() const{
    double r = 0.0;
    switch(_type){
        case PREDICTION_THROUGHPUT:{
            r = 1.0 / getMaximumThroughput();
        }break;
        case PREDICTION_POWER:{
            r = getCurrentPower();
        }break;
    }
    return r;
}

void PredictorLinearRegression::refine(){
    KnobsValues currentValues = _configuration.getRealValues();
    _preparationNeeded = true;
    if(_type == PREDICTION_POWER){
        // Add service nodes
        currentValues[KNOB_VIRTUAL_CORES] += _configuration.getNumServiceNodes();
    }
    obs_it lb = _observations.lower_bound(currentValues);

    if(_p.regressionAging && !contains(_agingVector, currentValues)){
        if(_agingVector.size() < _p.regressionAging){
            _agingVector.push_back(currentValues);
        }else{
            _agingVector.at(_currentAgingId) = currentValues;
        }
        _currentAgingId = (_currentAgingId + 1) % _p.regressionAging;
    }
    double response = getCurrentResponse();
    DEBUG("Refining with configuration " << currentValues << ": "
                                         << response);
    if(lb != _observations.end() &&
       !(_observations.key_comp()(currentValues, lb->first))){
        // Key already exists
        DEBUG("Replacing " << currentValues);
        lb->second.data->init();
        lb->second.response = response;
    }else{
        // The key does not exist in the map
        Observation o;
        switch(_type){
            case PREDICTION_THROUGHPUT:{
                o.data = new RegressionDataServiceTime(_p, _configuration, _samples);
            }break;
            case PREDICTION_POWER:{
                o.data = new RegressionDataPower(_p, _configuration, _samples);
            }break;
        }
        o.response = response;
        _observations.insert(lb, Observations::value_type(currentValues, o));
    }
}

void PredictorLinearRegression::prepareForPredictions(){
    if(_preparationNeeded){
        if(!readyForPredictions()){
            throw std::runtime_error("prepareForPredictions: Not enough "
                                     "points are present");
        }

        // One observation per column.
        arma::mat dataMl(_observations.begin()->second.data->getNumPredictors(),
                         _observations.size());
        arma::vec responsesMl(_observations.size());

        size_t i = 0;
        for(obs_it iterator = _observations.begin();
                   iterator != _observations.end();
                   ++iterator){
            const Observation& obs = iterator->second;

            if(!_p.regressionAging || contains(_agingVector, iterator->first)){
                obs.data->toArmaRow(i, dataMl);
                responsesMl(i) = obs.response;
                ++i;
            }
        }

        // Remove rows in dataMl if they contain all the same value.
        // In this way we ensure that the matrix can be inverted and that
        // the regression has a solution.
        size_t id = 0;
        _removedRows.clear();
        for(size_t i = 0; i < dataMl.n_rows; ){
            double firstElem = dataMl(i, 0);
            bool removeRow = true;
            for(size_t j = 1; j < dataMl.n_cols; j++){
                if(dataMl(i, j) != firstElem){
                    removeRow = false;
                    break;
                }
            }
            if(removeRow){
                dataMl.shed_row(i);
                _removedRows.push_back(id);
            }else{
                i++;
            }
            ++id;
        }
        //// End of row removal ////

        if(_p.regressionAging && _agingVector.size() != _observations.size()){
            dataMl.resize(_observations.begin()->second.data->getNumPredictors(),
                          _agingVector.size());
            responsesMl.resize(_agingVector.size());
        }

        std::ostringstream x1;
        std::ostringstream x2;
        arma::set_stream_err1(x1);
        arma::set_stream_err2(x2);
        LinearRegression newlr = LinearRegression(dataMl, responsesMl);
        if(x1.str().compare("") || x2.str().compare("")){
            ;
        }else{
            _lr = newlr;
            _modelError =  _lr.ComputeError(dataMl, responsesMl);
        }

        DEBUG("Error in model: " << _modelError);
        _preparationNeeded = false;
    }
}

double PredictorLinearRegression::predict(const KnobsValues& values){
    _predictionInput->init(values);

    // One observation per column.
    arma::mat predictionInputMl(_predictionInput->getNumPredictors(), 1);
    arma::vec result(1);
    double res = 0;
    _predictionInput->toArmaRow(0, predictionInputMl);

    // Remove from prediction input the corresponding positions
    //for(size_t i = 0; i < _removedRows.size(); i++){
    //    predictionInputMl.shed_row(_removedRows.at(i));
    //}
    arma::mat predictionInputMlShed(_predictionInput->getNumPredictors(), 1);
    size_t realSize = 0;
    for(size_t i = 0; i < _predictionInput->getNumPredictors(); i++){
        if(!contains(_removedRows, i)){
            predictionInputMlShed(realSize, 0) = predictionInputMl(i, 0);
            ++realSize;
        }
    }
    predictionInputMlShed.resize(realSize, 1);

    _lr.Predict(predictionInputMlShed, result);
    if(_type == PREDICTION_THROUGHPUT){
        res = 1.0 / result.at(0);
    }else{
        res = result.at(0);
    }

    return res;
}

double PredictorLinearRegression::getInactivePowerParameter() const{
    if(_type == PREDICTION_POWER){
        size_t pos;
        if(!dynamic_cast<RegressionDataPower*>(_predictionInput)->getInactivePowerPosition(pos)){
            throw std::runtime_error("Impossible to get inactive power parameter.");
        }
        return _lr.Parameters().at(pos);
    }else{
        throw std::runtime_error("getInactivePowerParameter can only be called on POWER predictors.");
    }
}

/**************** PredictorUSL ****************/
PredictorUsl::PredictorUsl(PredictorType type,
             const Parameters &p,
             const Configuration &configuration,
             const Smoother<MonitoredSample>* samples):
                    Predictor(type, p, configuration, samples),
                    _maxPolDegree(POLYNOMIAL_DEGREE_USL),
                    _ws(NULL), _x(NULL), _y(NULL), _chisq(0),
                    _preparationNeeded(true), _maxSpeedThr(0), _minSpeedThr(0),
                    _minSpeedCoresThr(0), _minSpeedCoresThrNew(0), _n1(0), _n2(0),
                    _minSpeedN1Thr(0), _minSpeedN2Thr(0){
    if(type != PREDICTION_THROUGHPUT){
        throw std::runtime_error("PredictorUsl can only be used for throughput predictions.");
    }
    if(_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_USLP){
        _maxPolDegree -= 1; //To remove x^0 value.
    }
    _c = gsl_vector_alloc(_maxPolDegree);
    _cov = gsl_matrix_alloc(_maxPolDegree, _maxPolDegree);

    _minSpeed = 1;
    _maxSpeed = 1;
    if(_p.knobFrequencyEnabled){
        std::vector<double> frequencies = _configuration.getKnob(KNOB_FREQUENCY)->getAllowedValues();
        _minSpeed = frequencies.front();
        _maxSpeed = frequencies.back();
    }
    if(_p.knobClkModEnabled){
        _minSpeed *= _configuration.getKnob(KNOB_CLKMOD)->getAllowedValues().front() / 100.0;
        _maxSpeed *= _configuration.getKnob(KNOB_CLKMOD)->getAllowedValues().back() / 100.0;
    }
}

PredictorUsl::~PredictorUsl(){
    gsl_matrix_free(_cov);
    gsl_vector_free(_c);
}

void PredictorUsl::clear(){
    _xs.clear();
    _ys.clear();
}

bool PredictorUsl::readyForPredictions(){
    return _xs.size() >= _maxPolDegree;
}

void PredictorUsl::refine(){
    double numCores = _configuration.getKnob(KNOB_VIRTUAL_CORES)->getRealValue();
    double throughput = getMaximumThroughput();
    double speed = 1;

    if(_p.knobFrequencyEnabled){
        speed = _configuration.getKnob(KNOB_FREQUENCY)->getRealValue();
    }

    if(_p.knobClkModEnabled){
        speed *= _configuration.getKnob(KNOB_CLKMOD)->getRealValue() / 100.0;
    }

    double maxCores = _configuration.getKnob(KNOB_VIRTUAL_CORES)->getAllowedValues().size();
    if(speed == _maxSpeed && numCores == maxCores){
        _maxSpeedThr = throughput;
        return;
    }else if(speed == _minSpeed){
        if(numCores == maxCores){
            _minSpeedThr = throughput;
        }else if(numCores == 1){
            _minSpeedCoresThr = throughput;
        }
    }else if(speed != _minSpeed){
        double minMaxScaling = _maxSpeedThr / _minSpeedThr;
        // We get the expected throughput at minimum frequency starting from the actual
        // throughput at generic frequency. To do so we apply the inverse of the function
        // used in predict() to get the throughput at a generic frequency starting
        // from that at minimum frequency.
        throughput = (throughput * (_maxSpeed - _minSpeed)) / (_maxSpeed - speed + minMaxScaling * (speed - _minSpeed));
        //speed = _minSpeed;
    }

    double x, y;
    if(_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_USLP){
        x = numCores - 1;
        y = ((_minSpeedCoresThr * numCores) / throughput) - 1;
    }else{
        x = numCores - 1;
        y = numCores / throughput;
    }

    // Checks if a y is already present for this x
    int pos = -1;
    for(size_t i = 0; i < _xs.size(); i++){
        if(_xs[i] == x){
            pos = i;
            break;
        }
    }
    if(pos == -1){
        _xs.push_back(x);
        _ys.push_back(y);
    }else{
        _ys.at(pos) = y;
    }
    _preparationNeeded = true;
}

void PredictorUsl::prepareForPredictions(){
    if(_preparationNeeded){
        if(!readyForPredictions()){
            throw std::runtime_error("PredictorUsl: Not yet ready for predictions.");
        }

        _x = gsl_matrix_alloc(_xs.size(), _maxPolDegree);
        _y = gsl_vector_alloc(_xs.size());
        size_t i = 0, j = 0;
        for(i = 0; i < _xs.size(); i++) {
            for(j = 0; j < _maxPolDegree; j++) {
                double degree = j;
                if(_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_USLP){
                    degree += 1; //To skip x^0 value.
                }
                gsl_matrix_set(_x, i, j, pow(_xs.at(i), degree));
            }
            gsl_vector_set(_y, i, _ys.at(i));
        }

        _ws = gsl_multifit_linear_alloc(_xs.size(), _maxPolDegree);
        gsl_multifit_linear(_x, _y, _c, _cov, &_chisq, _ws);

        _coefficients.clear();
        for(i = 0; i < _maxPolDegree; i++){
            _coefficients.push_back(gsl_vector_get(_c, i));
        }
        gsl_matrix_free(_x);
        gsl_vector_free(_y);
        gsl_multifit_linear_free(_ws);
        _preparationNeeded = false;
    }
}

double PredictorUsl::predict(const KnobsValues& knobsValues){
    const KnobsValues real = _configuration.getRealValues(knobsValues);
    double result = 0;
    double numCores = real[KNOB_VIRTUAL_CORES];
    double speed = 1;

    if(_p.knobFrequencyEnabled){
        speed = real[KNOB_FREQUENCY];
    }

    if(_p.knobClkModEnabled){
        speed *= real[KNOB_CLKMOD] / 100.0;
    }

    for(size_t i = 0; i < _coefficients.size(); i++){
        double exp = i;
        if(_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_USLP){
            // I do not have the constant factor b0*(x^0)
            exp += 1;
        }
        result += _coefficients.at(i)*std::pow(numCores - 1, exp);
    }

    double throughput;
    if(_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_USLP){
        throughput = (numCores * _minSpeedCoresThr)/(result + 1);
    }else{
        throughput = (numCores / result);
    }

    // If we are trying to predict the throughput with the
    // maximum number of cores, as value for throughput at minimum frequency
    // use the value we stored instead of the predicted one.
    if(numCores == _configuration.getKnob(KNOB_VIRTUAL_CORES)->getAllowedValues().size()){
        throughput = _minSpeedThr;
    }

    if(speed != _minSpeed){
        double minMaxScaling = _maxSpeedThr / _minSpeedThr;
        double maxFreqPred = throughput * minMaxScaling;
        return ((maxFreqPred - throughput)/(_maxSpeed - _minSpeed))*speed + ((throughput*_maxSpeed)-(maxFreqPred*_minSpeed))/(_maxSpeed - _minSpeed);
    }else{
        return throughput;
    }
}

double PredictorUsl::getMinSpeedCoresThr() const{
    if(_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_USLP){
        return _minSpeedCoresThr;
    }else{
        return 1.0 / _coefficients[0];
    }
}

double PredictorUsl::geta() const{
    if(_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_USLP){
        return _coefficients[1];
    }else{
        return _coefficients[2];
    }
}

double PredictorUsl::getb() const{
    if(_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_USLP){
        return _coefficients[0];
    }else{
        return _coefficients[1];
    }
}

void PredictorUsl::seta(double a){
    if(_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_USLP){
        _coefficients[1] = a;
    }else{
        _coefficients[2] = a;
    }
}

void PredictorUsl::setb(double b){
    if(_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_USLP){
        _coefficients[0] = b;
    }else{
        _coefficients[1] = b;
    }
}

double PredictorUsl::getTheta(RecordInterferenceArgument n){
    int numCores;
    int throughput;
    if(n == RECONFARG_N1){
        numCores = _n1;
        throughput = _minSpeedN1Thr;
    }else{
        numCores = _n2;
        throughput = _minSpeedN2Thr;
    }
    KnobsValues kv(KNOB_VALUE_REAL);
    kv[KNOB_FREQUENCY] = _minSpeed;
    kv[KNOB_VIRTUAL_CORES] = numCores;
    // TODO Does not consider clkMod
    if(_p.knobClkModEnabled){
        throw std::runtime_error("Interference cannot be updated when clock modulation is used.");
    }
    return (((getMinSpeedCoresThr()*numCores)/predict(kv)) - 1) - (((_minSpeedCoresThrNew*numCores)/throughput) - 1);
}

void PredictorUsl::updateInterference(){
    double numCores = _configuration.getRealValue(KNOB_VIRTUAL_CORES);
    double frequency = _configuration.getRealValue(KNOB_FREQUENCY);
    // TODO Does not consider clkMod
    if(_p.knobClkModEnabled){
        throw std::runtime_error("Interference cannot be updated when clock modulation is used.");
    }
    double throughput = _samples->average().throughput;
    // TODO: At the moment I'm assuming that interferences does not change the
    // slope when we fix number of cores and we change the frequency.
    if(frequency != _minSpeed){
        // We compute the 'old' throughput at minimum frequency
        KnobsValues kv(KNOB_VALUE_REAL);
        kv[KNOB_FREQUENCY] = _minSpeed;
        kv[KNOB_VIRTUAL_CORES] = numCores;
        double minFreqPred = predict(kv);
        kv[KNOB_FREQUENCY] = frequency;
        double currentFreqPred = predict(kv);
        double prop = throughput / currentFreqPred;
        throughput = prop*minFreqPred;
    }

    if(numCores == 1){
        _minSpeedCoresThrNew = throughput;
    }else if(!_n1){
        _n1 = numCores;
        _minSpeedN1Thr = throughput;
    }else{
        _n2 = numCores;
        _minSpeedN2Thr = throughput;
    }
}

void PredictorUsl::updateCoefficients(){
    double newB, newA, firstPartB, secondPartB;
    firstPartB = ((_n2 - 1)/(_n2 - _n1)*(_n1 - 1))*
                 (getLambda(RECONFARG_N1) - getTheta(RECONFARG_N1));
    secondPartB = ((_n1 - 1)/(_n2 - _n1)*(_n2 - 1))*
                  (getLambda(RECONFARG_N2) - getTheta(RECONFARG_N2));
    newB = firstPartB - secondPartB;
    newA = (getLambda(RECONFARG_N2) - getTheta(RECONFARG_N2) - newB*(_n2 - 1))/
            pow(_n2 - 1, 2);
    seta(newA);
    setb(newB);
    if(_p.strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_USL){
        _minSpeedCoresThr = _minSpeedCoresThrNew;
        _n1 = 0;
        _n2 = 0;
        _minSpeedCoresThrNew = 0;
    }
}

double PredictorUsl::getLambda(RecordInterferenceArgument n) const{
    int numCores;
    if(n == RECONFARG_N1){
        numCores = _n1;
    }else{
        numCores = _n2;
    }
    return getb()*(numCores - 1) + geta()*(numCores - 1)*(numCores - 1);
}

/**************** PredictorSimple ****************/

PredictorAnalytical::PredictorAnalytical(PredictorType type,
                                 const Parameters &p,
                                 const Configuration &configuration,
                                 const Smoother<MonitoredSample>* samples):
            Predictor(type, p, configuration, samples) {
    Topology* t = _p.mammut.getInstanceTopology();
    _phyCores = t->getPhysicalCores().size();
    _cpus = t->getCpus().size();
    _domains = _p.mammut.getInstanceCpuFreq()->getDomains().size();
    uint virtCoresPerPhyCores = t->getPhysicalCore(0)->getVirtualCores().size();
    _phyCoresPerDomain = _p.mammut.getInstanceCpuFreq()->getDomains().at(0)->getVirtualCores().size() / virtCoresPerPhyCores;
    _phyCoresPerCpu = t->getCpus().at(0)->getPhysicalCores().size();
    if(!_phyCoresPerDomain){ // Fallback
        _phyCoresPerDomain = _phyCoresPerCpu;
    }
}

double PredictorAnalytical::getScalingFactor(const KnobsValues& values){
    double usedVirtualCores = values[KNOB_VIRTUAL_CORES];
    return (double)(values[KNOB_FREQUENCY] * usedVirtualCores) /
           (double)(_configuration.getRealValue(KNOB_FREQUENCY) *
                    _configuration.getRealValue(KNOB_VIRTUAL_CORES));
}

double PredictorAnalytical::getPowerPrediction(const KnobsValues& values){
    assert(values.areReal());
    double usedPhysicalCores = values[KNOB_VIRTUAL_CORES];
    Frequency frequency = 1.0;
    double clkMod = 1.0;
    if(_p.knobFrequencyEnabled){
        frequency = values[KNOB_FREQUENCY];
    }
    if(_p.knobClkModEnabled){
        clkMod = values[KNOB_CLKMOD] / 100.0;
    }
    uint unusedDomains;
    double staticPower = 0, dynamicPower = 0;
    getStaticDynamicPower(usedPhysicalCores, frequency,
                          (MappingType) values[KNOB_MAPPING],
                          clkMod, _p, _cpus,
                          _domains, _phyCoresPerCpu,
                          _phyCoresPerDomain,
                          unusedDomains, staticPower, dynamicPower);
    return dynamicPower;
}

void PredictorAnalytical::prepareForPredictions(){
    ;
}

bool PredictorAnalytical::readyForPredictions(){
    return true;
}

void PredictorAnalytical::refine(){
    ;
}

void PredictorAnalytical::clear(){
    ;
}

double PredictorAnalytical::predict(const KnobsValues& values){
    switch(_type){
        case PREDICTION_THROUGHPUT:{
            double scalingFactor = getScalingFactor(values);
            return getMaximumThroughput() * scalingFactor;
        }break;
        case PREDICTION_POWER:{
            return getPowerPrediction(values);
        }break;
    }
    return 0.0;
}


PredictorLeo::PredictorLeo(PredictorType type,
              const Parameters &p,
              const Configuration &configuration,
              const Smoother<MonitoredSample>* samples):
                  Predictor(type, p, configuration, samples),
                  _preparationNeeded(true){
    const std::vector<KnobsValues>& combinations = _configuration.getAllRealCombinations();
    DEBUG("Found: " << combinations.size() << " combinations.");
    for(size_t i = 0; i < combinations.size(); i++){
        _confIndexes[combinations.at(i)] = i;
    }
    _values.resize(combinations.size());
    _values.zeros();
    std::vector<std::string> names = mammut::utils::readFile(p.leo.namesData);
    bool appFound = false;
    for(size_t i = 0; i < names.size(); i++){
        if(names.at(i).compare(p.leo.applicationName) == 0){
            _appId = i;
            appFound = true;
        }
    }
    if(!appFound){
        throw std::runtime_error("Impossible to find application " + p.leo.applicationName +
                                 " in " + p.leo.namesData);
    }
}

PredictorLeo::~PredictorLeo(){
    ;
}

bool PredictorLeo::readyForPredictions(){
    return true;
}

void PredictorLeo::clear(){
    _values.zeros();
    _predictions.zeros();
}

void PredictorLeo::refine(){
    _preparationNeeded = true;
    auto it = _confIndexes.find(_configuration.getRealValues());
    if(it == _confIndexes.end()){
        throw std::runtime_error("[Leo] Impossible to find index for configuration.");
    }
    size_t confId = it->second;
    if(confId >= _values.size()){
        throw std::runtime_error("[Leo] Invalid configuration index: " + confId);
    }
    switch(_type){
        case PREDICTION_THROUGHPUT:{
            // TODO In the offline profiling data we have throughputMax, not throughput
            // According, if we do not have PERF_UTILIZATION contracts and we have
            // a utilization < 1, results may be wrong.
#if 0
            if(_samples->average().utilisation < MAX_RHO){
                throw std::runtime_error("[Leo] throughput throughputMax problem.");
            }
#endif
            _values.at(confId) = getMaximumThroughput();
        }break;
        case PREDICTION_POWER:{
            _values.at(confId) = getCurrentPower();
        }break;
        default:{
            throw std::runtime_error("[Leo] Unknown predictor type.");
        }
    }
}

void PredictorLeo::prepareForPredictions(){
    if(_preparationNeeded){
        if(!readyForPredictions()){
            throw std::runtime_error("[Leo] prepareForPredictions: Not enough "
                                     "points are present");
        }

        string dataFile;
        bool perColumnNormalization;
        switch(_type){
            case PREDICTION_THROUGHPUT:{
                dataFile = _p.leo.throughputData;
                perColumnNormalization = true;
            }break;
            case PREDICTION_POWER:{
                dataFile = _p.leo.powerData;
                perColumnNormalization = false;
            }break;
            default:{
                throw std::runtime_error("Unknown predictor type.");
            }
        }

        leo::PredictionResults pr = leo::compute(_appId, dataFile,
                                                 &_values, perColumnNormalization);
        _predictions = arma::conv_to<std::vector<double> >::from(pr.predictions);
        _preparationNeeded = false;
    }
}

double PredictorLeo::predict(const KnobsValues& realValues){
    auto it = _confIndexes.find(realValues);
    if(it == _confIndexes.end()){
        throw std::runtime_error("[Leo] Impossible to find index for configuration.");
    }
    size_t confId = it->second;
    if(confId >= _predictions.size()){
        throw std::runtime_error("[Leo] Invalid configuration index: " + confId);
    }
    return _predictions.at(confId);
}

PredictorFullSearch::PredictorFullSearch(PredictorType type,
          const Parameters &p,
          const Configuration &configuration,
          const Smoother<MonitoredSample>* samples):
      Predictor(type, p, configuration, samples),
      _allConfigurations(_configuration.getAllRealCombinations()){
    ;
}

PredictorFullSearch::~PredictorFullSearch(){
    ;
}

bool PredictorFullSearch::readyForPredictions(){
    return _values.size() == _allConfigurations.size();
}

void PredictorFullSearch::clear(){
    _values.clear();
}

void PredictorFullSearch::refine(){
    double value = 0;
    switch(_type){
        case PREDICTION_THROUGHPUT:{
            value = getMaximumThroughput();
        }break;
        case PREDICTION_POWER:{
            value = getCurrentPower();
        }break;
        default:{
            throw std::runtime_error("Unknown predictor type.");
        }
    }
    _values[_configuration.getRealValues()] = value;
}

void PredictorFullSearch::prepareForPredictions(){
    ;
}

double PredictorFullSearch::predict(const KnobsValues& realValues){
    if(!readyForPredictions()){
        throw std::runtime_error("prepareForPredictions: Not enough "
                                 "points are present");
    }
    return _values.at(realValues);
}

	/**************** PredictorSMT****************/

	PredictorSMT::PredictorSMT(PredictorType type,
				const Parameters &p,
				const Configuration &configuration,
				const Smoother<MonitoredSample>* samples) :
				Predictor(type, p, configuration, samples)
    {

		//One observation per column
		_xs2.set_size(4, 0);  //Four regressors (3 for USL/Freq + 1 for HT)
		_xs1.set_size(3, 0);  //Three regressors (USL/Freq only)

		//getting min and max clock frequency
		_minFreq = 1;
		if (_p.knobFrequencyEnabled) {
			std::vector<double> frequencies = _configuration.getKnob(KNOB_FREQUENCY)->getAllowedValues();
			_minFreq = frequencies.front();
		}

		//getting number of CPU and cores per CPU
		mammut::topology::Topology* t = _p.mammut.getInstanceTopology();
		_phyCoresPerCpu = t->getCpus().at(0)->getPhysicalCores().size();
		_cpus = t->getCpus().size();
	}

        //Add (or update existing) y to _ys and x to _xs
    void PredictorSMT::addObservation(mat& _xs,
                       		 rowvec& _ys,
                       		 const vec& x,
                        	 const double y)
    {
        int col = -1;
        for (size_t i = 0; i < _xs.n_cols; i++) {
            bool found = true;

            for (size_t j = 0; j < _xs.n_rows; j++)
            {
                found &= (x(j) == _xs(j, i));
            }

            if (found) {
                col = i;
                break;
            }
        }
        if (col == -1) {
            _xs.insert_cols(_xs.n_cols, x);
            _ys.resize(_ys.size() + 1);
            _ys.at(_ys.size() - 1) = y;
        }
        else {
            _ys.at(col) = y;
        }
    }

    double PredictorSMT::getSigma(double numCores, double freq) const
    {
        return ((numCores - 1) * _minFreqCoresExtime * _minFreq / (numCores*freq));
    }

    double PredictorSMT::getKi(double numCores, double freq) const
    {
        return ((numCores - 1) * _minFreqCoresExtime * _minFreq / (freq));
    }

    double PredictorSMT::getGamma(double numCores, double freq) const
    {
        return (_minFreqCoresExtime*_minFreq / (numCores*freq));
    }

    double PredictorSMT::getHT(double numContexts) const
    {
        return (1 / numContexts) - 1;
    }

	void PredictorSMT::refine() {
		double numContexts = _configuration.getKnob(KNOB_HYPERTHREADING)->getRealValue();
		double numCores = (_configuration.getKnob(KNOB_VIRTUAL_CORES)->getRealValue()) / numContexts;
		double executionTime = 1.0 /(getMaximumThroughput());
		double power = getCurrentPower();
		double freq = 1.0;
		double nActiveCpu = ceil(numCores / _phyCoresPerCpu);

		if (_p.knobFrequencyEnabled) {
			freq = _configuration.getKnob(KNOB_FREQUENCY)->getRealValue();
		}

		switch (_type) {
		case PREDICTION_THROUGHPUT: {
		    if (numContexts == 1)
			{
			if (numCores == 1 && freq == _minFreq)
			    _minFreqCoresExtime = executionTime;

			vec x1(3);
			x1(0) = getSigma(numCores, freq);
			x1(1) = getKi(numCores, freq);
			x1(2) = getGamma(numCores, freq);

			addObservation(_xs1, _ys1, x1, executionTime);

			}

			vec x2(4);
			x2(0) = getSigma(numCores, freq);
			x2(1) = getKi(numCores, freq);
			x2(2) = getGamma(numCores, freq);
			x2(3) = getHT(numContexts);

			addObservation(_xs2, _ys2, x2, executionTime);

		}break;
		case PREDICTION_POWER: {

			double fullCPUvoltage = getVoltage(_p.archData.voltageTable, _phyCoresPerCpu, freq);

			vec x(4);
			x(0) = (_cpus - nActiveCpu) * getVoltage(_p.archData.voltageTable, 0, freq); //Static power of inutilized cpus
			x(1) = nActiveCpu * fullCPUvoltage; //Static power of fully utilized cpus
			x(2) = nActiveCpu * fullCPUvoltage * fullCPUvoltage* freq *_phyCoresPerCpu; //Dynamic power of fully utilized cpus
			x(3) = 1 - (1 / numContexts); //Hyper Threading power overhead

			addObservation(_xs2, _ys2, x, power);

		}break;
		default: {
			throw std::runtime_error("Unknown predictor type.");
		}
		}

		_preparationNeeded = true;
	}
	bool PredictorSMT::readyForPredictions() {
		
		uint minPoints = 6;
		return _xs2.n_cols >= minPoints;

	}
	void PredictorSMT::prepareForPredictions() {
		if (_preparationNeeded) {
			if (!readyForPredictions()) {
				throw std::runtime_error("PredictorSMT: Not yet ready for predictions.");
			}

			switch (_type) {
			case PREDICTION_THROUGHPUT: {

				//Regression for USL/Freq only
				_lr1 = new LinearRegression(_xs1, _ys1, 30.1, true);

				mat usl_freq = _xs2.submat(0, 0, _xs2.n_rows - 2, _xs2.n_cols - 1);
				mat ht = _xs2.submat(_xs2.n_rows - 1, 0, _xs2.n_rows - 1, _xs2.n_cols - 1);
				vec extime_noht(_xs2.n_cols);

				//Using trained model to obtain data needed for the second regression
				_lr1->Predict(usl_freq, extime_noht);

				rowvec extime_ht(_xs2.n_cols);
				for (size_t i = 0; i < _xs2.n_cols; i++) extime_ht(i) = _ys2(i) / extime_noht(i);

				//Regression for HT
				_lr2 = new LinearRegression(ht, extime_ht, 0.04, true);

			}break;
			case PREDICTION_POWER: {

				_lr2 = new LinearRegression(_xs2, _ys2, 0, true);

			}break;
			default: {
				throw std::runtime_error("Unknown predictor type.");
			}
			}

			_preparationNeeded = false;
		}
	}

	double PredictorSMT::predict(const KnobsValues& knobValues)
	{
		const KnobsValues realValues = _configuration.getRealValues(knobValues);

		double numContexts = realValues[KNOB_HYPERTHREADING];
		double numCores = (realValues[KNOB_VIRTUAL_CORES]) / numContexts;
		double freq = 1;
		double nActiveCpu = ceil(numCores / _phyCoresPerCpu);


		if (_p.knobFrequencyEnabled) {
			freq = realValues[KNOB_FREQUENCY];
		}

		switch (_type) {
		case PREDICTION_THROUGHPUT: {

			mat usl_freq(3, 1);
			usl_freq(0, 0) = getSigma(numCores, freq);
			usl_freq(1, 0) = getKi(numCores, freq);
			usl_freq(2, 0) = getGamma(numCores, freq);

			vec extime_noht(1);

			_lr1->Predict(usl_freq, extime_noht);

			mat ht(1, 1);
			ht(0, 0) = getHT(numContexts);
			vec extime_ht(1);

			_lr2->Predict(ht, extime_ht);

			extime_ht(0) *= extime_noht(0);

			return (1/extime_ht(0));


		}break;
		case PREDICTION_POWER: {

			double fullCPUvoltage = getVoltage(_p.archData.voltageTable, _phyCoresPerCpu, freq);

			mat p_regr(4, 1);
			p_regr(0, 0) = (_cpus - nActiveCpu) * getVoltage(_p.archData.voltageTable, 0, freq);
			p_regr(1, 0) = nActiveCpu * fullCPUvoltage;
			p_regr(2, 0) = nActiveCpu * fullCPUvoltage * fullCPUvoltage* freq *_phyCoresPerCpu;
			p_regr(3, 0) = 1 - (1 / numContexts);

			vec power_ht(1);

			_lr2->Predict(p_regr, power_ht);

			return power_ht(0);

		}break;
		default: {
			throw std::runtime_error("Unknown predictor type.");
		}
		}

	}

	void PredictorSMT::clear() {
		if (_type == PREDICTION_THROUGHPUT) delete _lr1;
		delete _lr2;

		_xs1.reset();	_xs2.reset();
		_ys1.reset();	_ys2.reset();

		_xs2.set_size(4, 0);
		_xs1.set_size(3, 0);

	}

	PredictorSMT::~PredictorSMT() {
		clear();
	}

}


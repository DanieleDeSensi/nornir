/*
 * selectors.hpp
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

#ifndef NORNIR_SELECTORS_HPP_
#define NORNIR_SELECTORS_HPP_

#include "predictors.hpp"
#include "explorers.hpp"
#include "configuration.hpp"
#include "knob.hpp"
#include <mammut/mammut.hpp>
#include <nornir/external/nelder-mead.h>

#include <memory>

namespace nornir{

class Selector: public mammut::utils::NonCopyable{
private:
    std::vector<CalibrationStats> _calibrationStats;
    uint _calibrationStartMs;
    uint64_t _calibrationStartTasks;
    mammut::Mammut _localMammut;
    mammut::energy::Counter* _joulesCounter;
    double _totalCalibrationTime;
    bool _calibrating;
    bool _ignoreViolations;
    uint _numPhyCores;  
protected:
    const Parameters& _p;
    const Configuration& _configuration;
    const Smoother<MonitoredSample>* _samples;
    size_t _numCalibrationPoints;
    KnobsValues _previousConfiguration;
    Smoother<double>* _bandwidthIn;
    bool _forced;
    bool _forcedReturned;
    KnobsValues _forcedConfiguration;
    bool _calibrationCoordination;
    bool _calibrationAllowed;
    u_int64_t _totalTasks;
    uint64_t _remainingTasks;

    /**
     * Checks if the throughput respects the required contract.
     * @param value The value to be checked.
     * @param conservative If true applies the conservativeValue.
     */
    bool isFeasibleThroughput(double value, bool conservative) const;

    /**
     * Checks if the latency respects the required contract.
     * @param value The value to be checked.
     * @param conservative If true applies the conservativeValue.
     */
    bool isFeasibleLatency(double value, bool conservative) const;

    /**
     * Checks if the power consumption respects the required contract.
     * @param value The value to be checked.
     * @param conservative If true applies the conservativeValue.
     */
    bool isFeasiblePower(double value, bool conservative) const;

    /**
     * Checks if the execution time respects the required contract.
     * @param value The value to be checked.
     * @param conservative If true applies the conservativeValue.
     */
    bool isFeasibleTime(double value, bool conservative) const;

    /**
     * Checks if the energy respects the required contract.
     * @param value The value to be checked.
     * @param conservative If true applies the conservativeValue.
     */
    bool isFeasibleEnergy(double value, bool conservative) const;

    /**
     * Checks if the utilization respects the required contract.
     * @param value The value to be checked.
     * @param conservative If true applies the conservativeValue.
     */
    bool isFeasibleUtilization(double value, bool conservative) const;

    /**
     * Initializes the best value.
     * @return The best value seed.
     */
    double initBestValue() const;

    /**
     * Initializes the suboptimal value.
     * @return The suboptimal value seed.
     */
    double initBestSuboptimalValue() const;

    /**
     * Checks if the current configuration is the most performing one.
     * @return true if the current configuration is the most performing one,
     *         false otherwise.
     */
    virtual bool isMaxPerformanceConfiguration() const = 0;

    /**
     * Checks if the contract is violated.
     * @return true if the contract has been violated, false otherwise.
     */
    bool isContractViolated() const;

    /**
     * Starts the recording of calibration stats.
     */
    void startCalibration();

    bool areKnobsValid(const KnobsValues& kv);  
public:
    Selector(const Parameters& p,
             const Configuration& configuration,
             const Smoother<MonitoredSample>* samples);

    virtual ~Selector(){;}


    /**
     * Forces the selector on a specific configuration.
     * @param kv The configuration.
     */
    void forceConfiguration(KnobsValues& kv);

    /**
     * Updates the total number of tasks.
     * MUST be called before calling getNextKnobsValues().
     * @param totalTasks The total processed tasks.
     */
    void updateTotalTasks(u_int64_t totalTasks);

    /**
     * Updates the input bandwidth history with the current value.
     * MUST be called before calling getNextKnobsValues().
     */
    virtual void updateBandwidthIn();

    /**
     * Returns the next values to be set for the knobs.
     *
     * @return The next values to be set for the knobs.
     */
    virtual KnobsValues getNextKnobsValues() = 0;

    /**
     * Returns the calibration statistics.
     * @return A vector containing the calibration statistics.
     */
    std::vector<CalibrationStats> getCalibrationsStats() const;

    /**
     * Returns true if the calibrator is in the calibration phase,
     * false otherwise.
     * @return true if the calibrator is in the calibration phase,
     * false otherwise.
     */
    bool isCalibrating() const;

    /**
     * Stops the recording of calibration stats.
     */
    void stopCalibration();

    /**
     * Returns the total calibration time.
     * @return The total calibration time.
     */
    double getTotalCalibrationTime() const;

    /**
     * Resets the total calibration time.
     */
    void resetTotalCalibrationTime();

    /**
     * If this function is called, the selector needs to coordinate with a
     * centralised manager before performing calibrations.
     */
    void setCalibrationCoordination();

    /**
     * Allows the selector to start calibration.
     */
    void allowCalibration();

    /**
     * Ignore the contract violations.
     */
    void ignoreViolations();

    /**
     * Considers the contract violations.
     */
    void acceptViolations();
};

/**
 * Returns the file from which the configuration is read.
 * The contained knobs values must be in relative form.
 **/
inline std::string getSelectorManualCliControlFile(){
    return getRuntimeDir(false) + std::string("/selector_manual.dat");
}

/**
 * Configuration is manually chosen by an
 * external entity (selector-manual) command line interface.
 */
class SelectorManualCli: public Selector{
public:
    SelectorManualCli(const Parameters& p,
                   const Configuration& configuration,
                   const Smoother<MonitoredSample>* samples);

    ~SelectorManualCli();
    KnobsValues getNextKnobsValues();
    bool isMaxPerformanceConfiguration() const{return false;} // Never used by this selector
};

/**
 * Configuration is manually chosen bythe
 * selector-manual web interface.
 */
class SelectorManualWeb: public Selector{
private:
    int _socket;
    bool _connected;
public:
    SelectorManualWeb(const Parameters& p,
                   const Configuration& configuration,
                   const Smoother<MonitoredSample>* samples);

    ~SelectorManualWeb();
    KnobsValues getNextKnobsValues();
    bool isMaxPerformanceConfiguration() const{return false;} // Never used by this selector
};

/**
 * Always returns the current configuration.
 */
class SelectorFixed: public Selector{
protected:
    bool isMaxPerformanceConfiguration() const{return false;} // Never used by this selector
public:
    SelectorFixed(const Parameters& p,
                  const Configuration& configuration,
                  const Smoother<MonitoredSample>* samples);

    ~SelectorFixed();

    KnobsValues getNextKnobsValues();
};

class ManagerMulti;

/**
 * A generic selector that uses the predictions of ALL the configurations
 * to find the best one.
 */
class SelectorPredictive: public Selector{
    friend class ManagerMulti;
private:
    std::unique_ptr<Predictor> _throughputPredictor;
    std::unique_ptr<Predictor> _powerPredictor;
    bool _feasible;
    KnobsValues _maxPerformanceConfiguration;
    double _maxPerformance;
    // Association between REAL values and observed data.
    std::map<KnobsValues, MonitoredSample> _observedValues;
    std::map<KnobsValues, double> _performancePredictions;
    std::map<KnobsValues, double> _powerPredictions;

    /**
     * Checks if the specified value to maximize/minimize
     * is better than the best found
     * up to now. If so, the new best value is stored.
     * @param throughput The throughput value.
     * @param latency The latency value.
     * @param utilization The utilization value.
     * @param power The power consumption value.
     * @param best The best found up to now.
     * @return true if it was a better value (best is modified).
     */
    bool isBestMinMax(double throughput, double latency, double utilization,
                      double power, double time, double energy, double& best);

    /**
     * Checks if the specified value to control
     * is better than the best suboptimal value found
     * up to now. If so, the new best value is stored.
     * @param throughput The throughput value.
     * @param latency The latency value.
     * @param utilization The utilization value.
     * @param power The power consumption value.
     * @param best The best found up to now.
     * @return true if it was a better value (best is modified).
     */
    bool isBestSuboptimal(double throughput, double latency, double utilization,
                          double power, double time, double energy, double& best);
protected:
    double _throughputPrediction;
    double _powerPrediction;

    /**
     * Updates the most performing configuration considering the predictions
     * made for a new configuration.
     * @param values The values of the new configuration.
     * @param primaryValue The primary predicted value.
     * @param secondaryValue The secondary predicted value.
     */
    void updateMaxPerformanceConfiguration(KnobsValues values, double performance);

    /**
     * Computes the best relative knobs values for the farm.
     * @return The best relative knobs values.
     */
    KnobsValues getBestKnobsValues();

    /**
     * Refines the models with current data.
     */
    void refine();

    /**
     * Updates the predictions for the next configuration.
     * @param next The next configuration.
     */
    void updatePredictions(const KnobsValues& next);

    /**
     * Returns true if the predictors are ready to perform predictions,
     * false otherwise.
     * @return True if the predictors are ready to perform predictions,
     * false otherwise.
     */
    bool predictorsReady() const;

    /**
     * Returns true if predictions have been performed, false otherwise.
     * @return True if predictions have been performed, false otherwise.
     */
    bool predictionsDone() const;

    /**
     * Returns true if the best solution found is a feasible solution,
     * false otherwise.
     */
    bool isBestSolutionFeasible() const;

    /**
     * Clears the predictors.
     */
    void clearPredictors();

    bool isMaxPerformanceConfiguration() const;
public:
    SelectorPredictive(const Parameters& p,
                       const Configuration& configuration,
                       const Smoother<MonitoredSample>* samples,
                       std::unique_ptr<Predictor> throughputPredictor,
                       std::unique_ptr<Predictor> powerPredictor);

    virtual ~SelectorPredictive();

    /**
     * Returns the next values to be set for the knobs.
     *
     * @return The next values to be set for the knobs.
     */
    virtual KnobsValues getNextKnobsValues() = 0;

    /**
     * Given a prediction, returns the real predicted throughput, i.e.
     * the minimum between the prediction and the input bandwidth.
     * @param prediction The predicted throughput.
     * @return The real predicted throughput, i.e.
     * the minimum between the prediction and the input bandwidth.
     */
    double getRealThroughput(double predicted) const;

    /**
     * Return the throughput prediction for a given configuration.
     * @param values The knobs values.
     * @return The throughput prediction for a given configuration.
     */
    double getThroughputPrediction(const KnobsValues& values);

    /**
     * Return the power consumption prediction for a given configuration.
     * @param values The knobs values.
     * @return The power consumption prediction for a given configuration.
     */
    double getPowerPrediction(const KnobsValues& values);

    /**
     * Returns a map with all the primary predictions.
     * @return A map with all the primary predictions.
     */
    const std::map<KnobsValues, double>& getPrimaryPredictions() const;

    /**
     * Returns a map with all the secondary predictions.
     * @return A map with all the secondary predictions.
     */
    const std::map<KnobsValues, double>& getSecondaryPredictions() const;

    Predictor* getPrimaryPredictor() const{return _throughputPredictor.get();}

    Predictor* getSecondaryPredictor() const{return _powerPredictor.get();}
};

/**
 * Selector described in PPL2016 paper (a simpler version was present in PDP2015 paper).
 */
class SelectorAnalytical: public SelectorPredictive{
private:
    uint _violations;
    bool _firstPointGenerated;
public:
    SelectorAnalytical(const Parameters& p,
                   const Configuration& configuration,
                   const Smoother<MonitoredSample>* samples);

    KnobsValues getNextKnobsValues();
    virtual void updateBandwidthIn();
};

/**
 * Selector described in PPAM2019 paper.
 */
class SelectorAnalyticalFull: public SelectorPredictive{
public:
    SelectorAnalyticalFull(const Parameters& p,
                       const Configuration& configuration,
                       const Smoother<MonitoredSample>* samples);

    KnobsValues getNextKnobsValues();
    virtual void updateBandwidthIn();
};

/**
 * A generic online learner selector.
 */
class SelectorLearner: public SelectorPredictive{
    friend class Manager;
private:
    Explorer* _explorer;
    bool _firstPointGenerated;
    uint _contractViolations;
    uint _accuracyViolations;
    uint _totalCalPoints;

    // Stuff used for interference adjustment.
    KnobsValues _beforeInterferenceConf;
    bool _updatingInterference;
    std::vector<KnobsValues> _interferenceUpdatePoints;

    KnobsValues getNextMeaningfulKnobsValues(Explorer* explorer);
  
    std::unique_ptr<Predictor> getPredictor(PredictorType type,
                                            const Parameters& p,
                                            const Configuration& configuration,
                                            const Smoother<MonitoredSample>* samples) const;
    
    /**
     * Starts producing data in order to update the models by
     * considering that another application is interfering with
     * this one.
     */
    void updateModelsInterference();

    /**
     * Returns true if the models have been updated and they are
     * ready to be used.
     */
    bool areModelsUpdated() const;
public:
    SelectorLearner(const Parameters& p,
                       const Configuration& configuration,
                       const Smoother<MonitoredSample>* samples);

    ~SelectorLearner();

    KnobsValues getNextKnobsValues();

    /**
     * Checks if the application phase changed.
     * @return true if the phase changed, false otherwise.
     */
    bool phaseChanged() const;

    /**
     * Checks the accuracy of the predictions.
     * @return True if the predictions were accurate, false otherwise.
     */
    bool isAccurate();
};

/**
 * Explores a fixed number of configurations before starting
 * making predictions.
 * Configuration is never changed once it has been found.
 */
class SelectorFixedExploration: public SelectorPredictive{
private:
    std::vector<KnobsValues> _confToExplore;
public:
    SelectorFixedExploration(const Parameters& p,
                   const Configuration& configuration,
                   const Smoother<MonitoredSample>* samples,
                   std::unique_ptr<Predictor> throughputPredictor,
                   std::unique_ptr<Predictor> powerPredictor,
                   size_t numSamples);

    ~SelectorFixedExploration();


    KnobsValues getNextKnobsValues();
};

#ifdef ENABLE_ARMADILLO
/**
 * Applies the algorithm described in:
 * "A Probabilistic Graphical Model-based Approach for Minimizing
 * Energy Under Performance Constraints" - Mishra, Nikita and Zhang, Huazhe
 * and Lafferty, John D. and Hoffmann, Henry
 */
class SelectorLeo: public SelectorFixedExploration{
public:
    SelectorLeo(const Parameters& p,
                   const Configuration& configuration,
                   const Smoother<MonitoredSample>* samples);

    ~SelectorLeo();
};
#endif

/**
 * This selector does a full exploration of the search space.
 */
class SelectorFullSearch: public SelectorFixedExploration{
public:
    SelectorFullSearch(const Parameters& p,
                   const Configuration& configuration,
                   const Smoother<MonitoredSample>* samples);

    ~SelectorFullSearch();
};

/**
 * A selector that implements the algorithm described in:
 * "Dynamic Power-Performance Adaptation of Parallel Computation
 * on Chip Multiprocessors" - Jian Li and Jose F. Martınez
 */
class SelectorLiMartinez: public Selector{
private:
    bool _firstPointGenerated;
    uint _low1, _mid1, _high1;
    uint _low2, _mid2, _high2;
    uint _midId;
    std::vector<mammut::cpufreq::Frequency> _availableFrequencies;
    double _currentWatts;
    double _optimalWatts;
    mammut::cpufreq::Frequency _optimalFrequency;
    uint _optimalWorkers;
    double _currentBw, _leftBw, _rightBw;
    KnobsValues _optimalKv;
    bool _improved;
    std::vector<double> _allowedCores;
    bool _optimalFound;

    mammut::cpufreq::Frequency findNearestFrequency(mammut::cpufreq::Frequency f) const;
    void goRight();
    void goLeft();
protected:
    bool isMaxPerformanceConfiguration() const;
public:
    SelectorLiMartinez(const Parameters& p,
                         const Configuration& configuration,
                         const Smoother<MonitoredSample>* samples);
    ~SelectorLiMartinez();

    KnobsValues getNextKnobsValues();
};

class SelectorHMPNelderMead: public Selector{
private:
    neme::NelderMeadOptimizer _opt;
    KnobsValues nmVectorToKv(neme::Vector v) const;
    neme::Vector kvToNmVector(KnobsValues kv) const;
    KnobsValues getFirstConfiguration() const;
    double nmScore() const;
    bool _firstGenerated;
    KnobsValues _lastRelative;
protected:
    bool isMaxPerformanceConfiguration() const{return false;} // Never used by this selector
public:
  SelectorHMPNelderMead(const Parameters& p,
                         const Configuration& configuration,
                         const Smoother<MonitoredSample>* samples);
  ~SelectorHMPNelderMead();

  KnobsValues getNextKnobsValues();
};

}



#endif /* NORNIR_SELECTORS_HPP_ */

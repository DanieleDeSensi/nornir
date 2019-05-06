/*
 * parameters.hpp
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
 *
 * Strategies used/developed for different papers:
 * - PDP2015:
 *         strategySelection = STRATEGY_SELECTION_ANALYTICAL
 *         knobMappingEnabled = false
 *         (Actually was a simpler selection strategy than the current one).
 * - PPL2016:
 *         strategySelection = STRATEGY_SELECTION_ANALYTICAL
 *         knobMappingEnabled = false
 * - PDP2016:
 *         strategySelection =  STRATEGY_SELECTION_LEARNING
 *         strategyPredictionPerformance = STRATEGY_PREDICTION_PERFORMANCE_AMDAHL
 *         strategyPredictionPower = STRATEGY_PREDICTION_POWER_LINEAR
 *         strategyExploration = (Many of them)
 *         knobMappingEnabled = false
 * - TACO2016:
 *         strategySelection = STRATEGY_SELECTION_LEARNING
 *         strategyPredictionPerformance = STRATEGY_PREDICTION_PERFORMANCE_USL
 *         strategyPredictionPower = STRATEGY_PREDICTION_POWER_LINEAR
 *         strategyUnusedVirtualCores = STRATEGY_UNUSED_VC_LOWEST_FREQUENCY
 *         strategyExploration = STRATEGY_EXPLORATION_HALTON
 *         knobMappingEnabled = true
 *
 *         In addition to that, we also developed and compared with:
 *             strategySelection = STRATEGY_SELECTION_LEO
 *             strategySelection = STRATEGY_SELECTION_FULLSEARCH
 *             strategySelection = STRATEGY_SELECTION_LIMARTINEZ
 *         and also with
 *             knobMappingEnabled = false
 *
 */


#ifndef NORNIR_PARAMETERS_HPP_
#define NORNIR_PARAMETERS_HPP_

#include "external/rapidXml/rapidxml.hpp"

#include <cmath>
#include <limits>
#include <fstream>
#include "external/mammut/mammut/utils.hpp"
#include "external/mammut/mammut/mammut.hpp"

#define NORNIR_REQUIREMENT_UNDEF 0
#define NORNIR_REQUIREMENT_MIN std::numeric_limits<double>::min()
#define NORNIR_REQUIREMENT_MAX std::numeric_limits<double>::max()

namespace nornir{

#define NORNIR_MANAGER_VIRTUAL_CORE (VirtualCoreId) 0

class Logger;

typedef enum{
    LOGGER_FILE = 0, // Log on file
    LOGGER_GRAPHITE // Log on graphite
}LoggerType;

// Possible knobs
typedef enum{
    KNOB_VIRTUAL_CORES = 0, // Number of contexts to be used.
    KNOB_HYPERTHREADING, // Number of contexts to be used on each physical core.
    KNOB_MAPPING, // Mapping of threads on physical cores.
    KNOB_FREQUENCY, // Clock frequency of the cores.
    KNOB_CLKMOD, // Clock modulation.
    KNOB_NUM  // <---- This must always be the last value
}KnobType;

/// Communication queues blocking/nonblocking.
typedef enum{
    // Non blocking queue.
    TRIGGER_Q_BLOCKING_NO = 0,

    // Blocking queue.
    TRIGGER_Q_BLOCKING_YES,

    // Automatically choose between blocking and nonblocking.
    TRIGGER_Q_BLOCKING_AUTO
}TriggerConfQBlocking;

/// Possible strategies to apply for unused virtual cores.
typedef enum{
    // Automatically choose one of the other strategies.
    STRATEGY_UNUSED_VC_AUTO = 0,

    // The unused virtual cores will run at the same
    // frequency of the used ones.
    STRATEGY_UNUSED_VC_SAME,

    // Nothing is done on unused virtual cores.
    STRATEGY_UNUSED_VC_NONE,

    // Set the virtual cores to the lowest frequency (only
    // possible if all the other virtual cores on the same
    // domain are unused).
    STRATEGY_UNUSED_VC_LOWEST_FREQUENCY,

    // Turn off the virtual cores. They will not be anymore seen by the
    // operating system and it will not schedule anything on them.
    STRATEGY_UNUSED_VC_OFF
}StrategyUnusedVirtualCores;

// Possible strategies for the selection of the best configuration.
typedef enum{
	// Best configuration is manually chosen by using the selector-manual command line interface.
	STRATEGY_SELECTION_MANUAL_CLI = 0,

    // Best configuration is manually chosen by using the selector-manual web interface.
    STRATEGY_SELECTION_MANUAL_WEB,

    // Applies an online learning algorithm
    STRATEGY_SELECTION_LEARNING,

    // Applies a simple analytical model.
    STRATEGY_SELECTION_ANALYTICAL,

    // Tries all the configurations in order to find the best one.
    STRATEGY_SELECTION_FULLSEARCH,

    // Applies the algorithm described in:
    // "Dynamic Power-Performance Adaptation of Parallel Computation
    // on Chip Multiprocessors" - Jian Li and Jose F. Martınez
    STRATEGY_SELECTION_LIMARTINEZ,

    // Applies the algorithm described in:
    // "A Probabilistic Graphical Model-based Approach for Minimizing
    // Energy Under Performance Constraints" - Mishra, Nikita and Zhang, Huazhe
    // and Lafferty, John D. and Hoffmann, Henry
    STRATEGY_SELECTION_LEO,

    STRATEGY_SELECTION_NUM // <- Must always be the last.
}StrategySelection;

// Possible prediction strategies for performance. Can only be specified if the
// selection strategy is "LEARNING".
typedef enum{
    STRATEGY_PREDICTION_PERFORMANCE_AMDAHL = 0,
    STRATEGY_PREDICTION_PERFORMANCE_USL,
    STRATEGY_PREDICTION_PERFORMANCE_USLP, // <- More precise than USL but needs one additional calibration point.
    STRATEGY_PREDICTION_PERFORMANCE_LEO,
    STRATEGY_PREDICTION_PERFORMANCE_SMT, // <-  Support for Simultaneous Multi-Threading
    STRATEGY_PREDICTION_PERFORMANCE_NUM // <- This must always be the last.
}StrategyPredictionPerformance;

// Possible prediction strategies for power consumption. Can only be specified
// if the selection strategy is "LEARNING".
typedef enum{
    STRATEGY_PREDICTION_POWER_LINEAR = 0,
    STRATEGY_PREDICTION_POWER_LEO,
    STRATEGY_PREDICTION_POWER_SMT, // <-  Support for Simultaneous Multi-Threading
    STRATEGY_PREDICTION_POWER_NUM // <- This must always be the last.
}StrategyPredictionPower;

/// Possible ways to select the calibration points. Can only be specified if
/// the selection strategy is "LEARNING".
typedef enum{
    // Random choice of the points.
    STRATEGY_EXPLORATION_RANDOM = 0,
    // Bratley, Fox, Niederreiter, ACM Trans. Model. Comp. Sim. 2, 195 (1992)
    STRATEGY_EXPLORATION_NIEDERREITER,
    // Antonov, Saleev, USSR Comput. Maths. Math. Phys. 19, 252 (1980)
    STRATEGY_EXPLORATION_SOBOL,
    // J.H. Halton, Numerische Mathematik 2, 84-90 (1960) and B. Vandewoestyne
    // R. Cools Computational and Applied Mathematics 189, 1&2, 341-361 (2006)
    STRATEGY_EXPLORATION_HALTON,
    STRATEGY_EXPLORATION_HALTON_REVERSE,
}StrategyExploration;

/// Possible ways to smooth the values.
typedef enum{
    // Simple moving average
    STRATEGY_SMOOTHING_MOVING_AVERAGE = 0,

    // Exponential moving average
    STRATEGY_SMOOTHING_EXPONENTIAL
}StrategySmoothing;

/// Possible ways to change the smoothing factor at runtime.
typedef enum{
    // Constant smoothing factor
    STRATEGY_SMOOTHING_FACTOR_CONST = 0,

    // Changes the smoothing factor according to the variation
    // of the data.
    STRATEGY_SMOOTHING_FACTOR_DYNAMIC
}StrategySmoothingFactor;

/// Possible ways to act when the manager finds no response on a queue
/// from a farm node.
typedef enum{
    // Spins.
    STRATEGY_POLLING_SPINNING = 0,

    // Suspend the execution on the processor for a short period of time.
    // However, the OS will see the thread as doing useful work and probably
    // will not get descheduled.
    STRATEGY_POLLING_PAUSE,

    // Sleeps for 0 nanoseconds. The OS may deschedule the thread.
    STRATEGY_POLLING_SLEEP_SMALL,

    // Sleeps for a period of time equal to the node average latency. The
    // OS may deschedule the thread.
    STRATEGY_POLLING_SLEEP_LATENCY
}StrategyPolling;

// How to decide if we have to consider new configurations
typedef enum{
    // Number of samples.
    STRATEGY_PERSISTENCE_SAMPLES = 0,

    // Coefficient of variation.
    STRATEGY_PERSISTENCE_VARIATION
}StrategyPersistence;

// Phase detection algorithm.
typedef enum{
    STRATEGY_PHASE_DETECTION_NONE = 0,
    STRATEGY_PHASE_DETECTION_TRIVIAL
}StrategyPhaseDetection;

/// Possible parameters validation results.
typedef enum{
    // Parameters are ok.
    VALIDATION_OK = 0,

    // Generic error
    VALIDATION_NO,

    // Wrong requirement or combination of requirements.
    VALIDATION_WRONG_REQUIREMENT,

    // Some of the knobs enabled by the user is not supported
    // by the selection algorithm he chosen.
    VALIDATION_UNSUPPORTED_KNOBS,

    // On this architecture is not possible to manually change
    // the frequency.
    VALIDATION_NO_MANUAL_DVFS,

    // Frequency can be changed by the operating system and the flag
    // "constant_tsc" is not present on the CPU. Accordingly, since we
    // realy on getticks() to perform measurements, the amount of ticks
    // per second may change with the frequency, probably causing inaccuracy
    // in the measurements.
    VALIDATION_NO_CONSTANT_TSC,

    // strategyFrequencies can be different from STRATEGY_FREQUENCY_NO
    // only if strategyMapping is different from STRATEGY_MAPPING_NO.
    VALIDATION_STRATEGY_FREQUENCY_REQUIRES_MAPPING,

    // The specified frequency strategy is not supported on this machine.
    VALIDATION_STRATEGY_FREQUENCY_UNSUPPORTED,

    // Specified governor not supported on this machine.
    VALIDATION_GOVERNOR_UNSUPPORTED,

    // sensitiveEmitter or sensitiveCollector specified but frequency
    // strategy is STRATEGY_FREQUENCY_NO.
    VALIDATION_EC_SENSITIVE_WRONG_F_STRATEGY,

    // sensitiveEmitter or sensitiveCollector specified but highest
    // frequency can't be set.
    VALIDATION_EC_SENSITIVE_MISSING_GOVERNORS,

    // The bounds are invalid or the frequency strategy is not
    // STRATEGY_FREQUENCY_OS.
    VALIDATION_INVALID_FREQUENCY_BOUNDS,

    // Strategy for unused virtual cores requires turning off the virtual
    // cores but they can't be turned off.
    VALIDATION_UNUSED_VC_NO_OFF,

    // Strategy for unused virtual cores requires lowering the frequency but
    // frequency scaling not available.
    VALIDATION_UNUSED_VC_NO_FREQUENCIES,

    // Blocking threshold needs to be specified.
    VALIDATION_BLOCKING_PARAMETERS,

    // Parameters for Leo predictors not specified.
    VALIDATION_NO_LEO_PARAMETERS,

    // Clock modulation not available.
    VALIDATION_NO_CLKMOD,
}ParametersValidation;

/**
 * This class represents an XML tree.
 */
class XmlTree{
private:
    char* _fileContentChars;
    rapidxml::xml_node<>* _root;

    void init(const std::string& content, const std::string& rootName);
public:
    /**
     * Reads the XML tree contained in 'fileName' and starting
     * with a root named 'rootName'.
     * @param fileName The name of the XML file.
     * @param rootName The name of the root.
     */
    XmlTree(const std::string& fileName, const std::string& rootName);

    ~XmlTree();


    /**
     * Returns a node with name valueName.
     * @param valueName The name of the node. It can be of the form
     *        field.subfield.subsubfield. ... (similar to a C struct).
     * @return The node with name valueName.
     */
    rapidxml::xml_node<>* getNode(const char* valueName) const;

    /**
     * If an element named 'valueName' is present, its value is
     * copied into 'value' and interpreted as a bool.
     * @param valueName The name of the XML element.
     * @param value The value that will contain the element
     *              valueName (if present).
     */
    void getBool(const char* valueName, bool& value);

    /**
     * If an element named 'valueName' is present, its value is
     * copied into 'value' and interpreted as an int.
     * @param valueName The name of the XML element.
     * @param value The value that will contain the element
     *              valueName (if present).
     */
    void getInt(const char* valueName, int& value);

    /**
     * If an element named 'valueName' is present, its value is
     * copied into 'value' and interpreted as an uint.
     * @param valueName The name of the XML element.
     * @param value The value that will contain the element
     *              valueName (if present).
     */
    void getUint(const char* valueName, uint& value);

    /**
     * If an element named 'valueName' is present, its value is
     * copied into 'value' and interpreted as an ulong.
     * @param valueName The name of the XML element.
     * @param value The value that will contain the element
     *              valueName (if present).
     */
    void getUlong(const char* valueName, ulong& value);

    /**
     * If an element named 'valueName' is present, its value is
     * copied into 'value' and interpreted as a double.
     * @param valueName The name of the XML element.
     * @param value The value that will contain the element
     *              valueName (if present).
     */
    void getDouble(const char* valueName, double& value);

    /**
     * If an element named 'valueName' is present, its value is
     * copied into 'value' and interpreted as a double. If the value
     * of the element was the string "MIN", the minimum possible
     * double is copied into 'value'.
     * @param valueName The name of the XML element.
     * @param value The value that will contain the element
     *              valueName (if present).
     */
    void getDoubleOrMin(const char* valueName, double& value);

    /**
     * If an element named 'valueName' is present, its value is
     * copied into 'value' and interpreted as a double. If the value
     * of the element was the string "MAX", the maximum possible
     * double is copied into 'value'.
     * @param valueName The name of the XML element.
     * @param value The value that will contain the element
     *              valueName (if present).
     */
    void getDoubleOrMax(const char* valueName, double& value);

    /**
     * If an element named 'valueName' is present, its value is
     * copied into 'value' and interpreted as a string.
     * @param valueName The name of the XML element.
     * @param value The value that will contain the element
     *              valueName (if present).
     */
    void getString(const char* valueName, std::string& value);

    /**
     * If an element named 'valueName' is present, its value is
     * copied into 'value' and interpreted as an array of uint
     * (separated by ':').
     * @param valueName The name of the XML element.
     * @param value The value that will contain the element
     *              valueName (if present).
     */
    void getArrayUint(const char* valueName, std::vector<uint>& value);

    /**
     * If an element named 'valueName' is present, its value is
     * copied into 'value' and interpreted as an array of enums
     * (separated by ':').
     * @param valueName The name of the XML element.
     * @param value The value that will contain the element
     *              valueName (if present).
     */
    template<typename T> void getArrayEnums(const char* valueName,
                                            std::vector<T>& value);

    /**
     * If an element named 'valueName' is present, its value is
     * copied into 'value' and interpreted as an enumeration.
     * @param T The type of the enumeration.
     * @param valueName The name of the XML element.
     * @param value The value that will contain the element
     *              valueName (if present).
     */
    template<typename T> void getEnum(const char* valueName, T& value);
};

#define SETVALUE(XML, TYPE, NAME) XML.get##TYPE(#NAME, NAME)

/*!
 * \class ArchData
 * \brief This class contains data specific to the architecture,
 *        obtained by microbenchmarks.
 *
 * This class contains data specific to the architecture,
 * obtained by microbenchmarks.
 */
typedef struct ArchData{
    // Number of ticks for a nanosecond.
    double ticksPerNs;

    // Number of ticks spent by a worker to reply to a monitoring request.
    double monitoringCost;

    mammut::cpufreq::VoltageTable voltageTable;

    ArchData():ticksPerNs(0),
               monitoringCost(0){;}

    void loadXml(const std::string& archFileName);
}ArchData;

/*!
 * \class Requirements
 * \brief This class contains the requirements specified by the user.
 *
 * This class contains the requirements specified by the user.
 */
typedef struct Requirements{
    // The throughput required for the application (expressed as tasks/sec).
    // It must be greater or equal than 0.
    // [default = unused].
    double throughput;

    // The maximum cores power to be used.
    // It must be greater or equal than 0.
    // [default = unused].
    double powerConsumption;

    // The underload threshold for the entire farm.
    // It must be in [0, 100].
    // [default = unused].
    double minUtilization;

    // The overload threshold for the entire farm.
    // It must be in [0, 100].
    // [default = unused].
    double maxUtilization;

    // The required completion time for the application (in seconds).
    // It must be greater or equal than 0. When used, 'expectedTasksNumber'
    // must be specified.
    // [default = unused].
    double executionTime;

    // The maximum latency required for each input element processed by the
    // application (in milliseconds).
    // It must be greater or equal than 0.
    // NOT AVAILABLE AT THE MOMENT.
    double latency;

    // The number of task expected for this computation.
    // [default = unused].
    double expectedTasksNumber;

    Requirements();

    /**
     * Checks if any requirement has been specified.
     * @return true if some requirement have been specified.
     */
    bool anySpecified() const;
}Requirements;

typedef struct{
    /**
     * The name of this application (it must matches with one of the names
     * in applicationNames file).
     */
    std::string applicationName;
    /**
     * File containing the names of the applications, in the same
     * order of the columns in throughput and power data.
     */
    std::string namesData;
    /**
     * File containing the data about the throughput. One column per application,
     * one row per configuration. Data doesn't need to be normalized.
     */
    std::string throughputData;
    /**
     * File containing the data about the power. One column per application,
     * one row per configuration. Data doesn't need to be normalized.
     */
    std::string powerData;
    /**
     * Number of samples to be used [default = 20].
     */
    uint numSamples;
}LeoParameters;

typedef struct{
    /**
     * If true, the interpreter will ensure that instructions belonging
     * to different stream elements will be processed in the same order
     * they are received [default = false].
     */
    bool orderedProcessing;

    /**
     * If true, the interpreter will ensure that the output will be produced
     * in the same order of the received data [default = false]. If
     * orderedProcessing is true, this property is automatically guaranteed.
     */
    bool orderedOutput;

    /**
     * Maximum number of graphs to keep in the system [default = 1000].
     */
    uint maxGraphs;

    /**
     * Maximum number of interpreters to be used [default = #Physical cores
     * in the system - 2].
     */
    uint maxInterpreters;
}DataflowParameters;

/*!
 * \class Parameters
 * \brief This class contains nornir parameters.
 *
 * This class contains nornir parameters.
 */
class Parameters{
    friend class Manager;
private:
    // The type of loggers to be activated.
    // It MUST be private because only the manager
    // needs to access it. The loggers created by
    // the manager will be deleted by the manager.
    // If the user needs a new logger it must add it directly
    // to the loggers vector, not to the loggersType.
    std::vector<LoggerType> loggersTypes;

    bool _knobEnabled[KNOB_NUM];

    /**
     * Sets default parameters
     */
    void setDefault();

    /**
     * Computes the sampling interval such to have a low
     * performance overhead.
     * @return The sampling interval.
     */
    uint getLowOverheadSamplingInterval() const;

    /**
     * Sets the default values for parameters that depends
     * from others.
     */
    void setDefaultPost();

    /**
     * Checks if a governor is available on the current machine.
     * @param g The governor.
     * @return true if the governor is available, false otherwise.
     */
    bool isGovernorAvailable(mammut::cpufreq::Governor g);

    /**
     * Gets the frequencies available on this machine.
     * @return A vector containing the frequencies available on this machine.
     */
    std::vector<mammut::cpufreq::Frequency> getAvailableFrequencies();

    /**
     * Checks if it is possible to turn off the cores.
     * @return true if it is possible to turn off the cores,
     *         false otherwise.
     */
    bool isUnusedVcOffAvailable();

    /**
     * Checks if it is possible to set the frequency on this machine.
     * @return true if it is possible to set the frequency on this machine,
     *         false otherwise.
     */
    bool isFrequencySettable();

    /**
     * Checks if it is possible to set the lowest frequency on this machine.
     * @return true if it is possible to set the lowest frequency on this
     *         machine, false otherwise.
     */
    bool isLowestFrequencySettable();

    /**
     * Checks if it is possible to set the highest frequency on this machine.
     * @return true if it is possible to set the highest frequency on this
     *         machine, false otherwise.
     */
    bool isHighestFrequencySettable();

    /**
     * Validates the strategy required for unused virtual cores.
     * @param s The strategy required for unused virtual cores.
     * @return The result of the validation.
     */
    ParametersValidation validateUnusedVc(StrategyUnusedVirtualCores& s);

    /**
     * Validates the frequency knob.
     * @return The result of the validation.
     */
    ParametersValidation validateKnobFrequencies();

    /**
     * Validates the clock modulation knob.
     * @return The result of the validation.
     **/
    ParametersValidation validateKnobClkMod();

    /**
     * Validates the triggers.
     * @return The result of the validation.
     */
    ParametersValidation validateTriggers();

    /**
     * Validates the required contract.
     * @return The result of the validation.
     */
    ParametersValidation validateRequirements();

    /**
     * Validates the selector.
     * @return The result of the validation.
     */
    ParametersValidation validateSelector();

    /**
     * Loads the content of the parameters with the content
     * of an XML file.
     * @param fileName The name of the XML file.
     */
    void loadXml(const std::string& fileName);

public:
    // The mammut modules handler.
    mammut::Mammut mammut;

    // Architecture's specific data.
    ArchData archData;

    // Requirements specified by the user.
    Requirements requirements;

    // The Q blocking knob [default = KNOB_Q_BLOCKING_NO].
    TriggerConfQBlocking triggerQBlocking;

    // Strategy for unused virtual cores [default = STRATEGY_UNUSED_VC_SAME].
    StrategyUnusedVirtualCores strategyUnusedVirtualCores;

    // Strategy to be used to select the best configuration according
    // to the requirements [default = STRATEGY_SELECTION_LEARNING].
    StrategySelection strategySelection;

    // Strategy to be used to predict performance values.
    // Only valid when strategySelection is LEARNING
    // [default = STRATEGY_PREDICTION_PERFORMANCE_USL].
    StrategyPredictionPerformance strategyPredictionPerformance;

    // Strategy to be used to predict power values.
    // Only valid when strategySelection is LEARNING
    // [default = STRATEGY_PREDICTION_POWER_LINEAR].
    StrategyPredictionPower strategyPredictionPower;

    // Strategy to be used to select the points to be explored during
    // calibration. Only valid when strategySelection is LEARNING
    // [default = STRATEGY_EXPLORATION_HALTON].
    StrategyExploration strategyExploration;

    // Smoothing strategy [default = STRATEGY_SMOOTHING_EXPONENTIAL].
    StrategySmoothing strategySmoothing;

    // Polling strategy [default = STRATEGY_POLLING_SLEEP_SMALL].
    StrategyPolling strategyPolling;

    // Persistence strategy [default = STRATEGY_PERSISTENCE_SAMPLES].
    StrategyPersistence strategyPersistence;

    // Phase detection strategy [default = STRATEGY_PHASE_DETECTION_NONE].
    StrategyPhaseDetection strategyPhaseDetection;

    // Flag to enable/disable cores knobs autotuning [default = true].
    bool knobCoresEnabled;

    // Flag to enable/disable mapping knob autotuning [default = true].
    bool knobMappingEnabled;

    // Flag to enable/disable frequency knob autotuning [default = true].
    bool knobFrequencyEnabled;

    // Flag to enable/disable clock modulation autotuning [default = false].
    bool knobClkModEnabled;

    // Flag to enable/disable hyperthreading knob autotuning [default = false].
    bool knobHyperthreadingEnabled;

    // Specifies a fixed value to apply to a knob (in the range [0, 100])
    // It will only be considered if knobHyperthreadingEnabled = false.
    // Otherwise, value of the knob will be autotuned. [default = 0]
    double knobHyperthreadingFixedValue; // TODO Do also for other knobs

    // Number of active threads in the application. Useful for
    // external managers. If 0, number of active threads is unknown [default = 0].
    uint32_t activeThreads;

    // If true, concurrency throttling is used on fastflow/nornir applications [default = true].
    bool useConcurrencyThrottling;

    // Parameters for LEO predictor.
    LeoParameters leo;

    // If true, before changing the number of workers the frequency will be
    // set to maximum to reduce the latency of the reconfiguration. The
    // frequency will be be set again to the correct value after the farm
    // is restarted [default = true].
    bool fastReconfiguration;

    // If true, when a reconfiguration occur, the collector is migrated to a
    // different virtual core (if needed) [default = false].
    bool migrateCollector;

    // The smoothing factor. It's meaning changes according to the smoothing
    // strategy adopted. If 0, default value will be set.
    // [default = 10 for moving average, 0.5 for exponential].
    double smoothingFactor;

    // The persistence value. It's meaning changes according to the persistence
    // strategy adopted. If 0, default value will be set.
    // [default = 1 for samples, 5 for variation].
    double persistenceValue;

    // The number of milliseconds required for the system to cooldown
    // after a reconfiguration. Since this is a transient state, the statistics
    // collected in the cooldownPeriod will be discarded [default = 200].
    double cooldownPeriod;

    // The length of the sampling interval (in milliseconds) for the data
    // reading during calibration phase. If 0, it will be automatically computed
    // such to have a low performance overhead [default = 100].
    uint32_t samplingIntervalCalibration;

    // The length of the sampling interval (in milliseconds) for the data
    // reading during steady phase. If 0, it will be automatically computed
    // such to have a low performance overhead [default = 1000].
    uint32_t samplingIntervalSteady;

    // If the application stays in the steady phase for at least
    // steadyThreshold samples, then we consider the application to be stedy
    // [default = 4].
    uint32_t steadyThreshold;

    // The minimum number of tasks in a worker sample. If 0, no minimum.
    // [default = 0].
    uint minTasksPerSample;

    // If true, this is an application when the workers works synchronously,
    // i.e. the emitter works as a barrier between two successive iterations
    // of the workers. In this case, the emitter always broadcasts tasks
    // to all the workers [default = false].
    bool synchronousWorkers;

    // Maximum calibration time (milliseconds). 0 is no limit.
    // We will keep calibrating until the error is higher than the
    // max*PredictionError AND calibration time is lower than the
    // maxCalibrationTime  AND the number of visited configurations
    // is lower than maxCalibrationConfigurations [default = 0.0].
    double maxCalibrationTime;

    // Maximum number of configurations to explore during calibration
    // phase. 0 is no limit. We will keep calibrating until the error is
    // higher than the max*PredictionError AND calibration time is lower
    // than the maxCalibrationTime AND the number of visited configurations
    // is lower than maxCalibrationSteps [default = 0.0].
    uint maxCalibrationSteps;

    // Maximum error percentage allowed for performance prediction.
    // [default = 10.0].
    double maxPerformancePredictionError;

    // Maximum error percentage allowed for power consumption prediction.
    // [default = 5.0].
    double maxPowerPredictionError;

    // The aging value for the linear regression predictor. If it has a value
    // of n, then we will only consider the last n configurations we collected
    // in order to perform predictions. If 0, no aging will be applied and all
    // the previous samples will be considered [default = 0].
    uint regressionAging;

    // The maximum percentage of monitoring overhead, in the range (0, 100).
    // [default = 1.0].
    double maxMonitoringOverhead;

    // If true, clock modulation will be emulated, instead of using hardware
    // clock modulation [default = true].
    bool clockModulationEmulated;

    // The minimum value that can be used for clock modulation, in the range (0, 100).
    // Only valid if clockModulationEmulated = true.
    // [default = 1].
    double clockModulationMin;

    // Resolution of the clock modulation, in the range (0, 100).
    // Only valid if clockModulationEmulated = true.
    // [default = 2].
    double clockModulationResolution;

    // Interval of the clock modulation, in microseconds.
    // Only valid if clockModulationEmulated = true.
    // [default = 100000].
    double clockModulationInterval;

    // Idle threshold (in microseconds) to switch from non blocking to blocking
    // runtime support. If the current idle time goes above this value,
    // and if the runtime has been configured to do so, it will switch
    // from non blocking to blocking support (and viceversa). If -1.0,
    // the runtime will never switch [default = -1.0].
    double thresholdQBlocking;

    // Safe range size around the threshold. The switch will occur
    // when idle time (or utilisation) is lesser than
    // thresholdQBlocking - thresholdQBlocking*thresholdQBlockingBelt or
    // greater than
    // thresholdQBlocking + thresholdQBlocking*thresholdQBlockingBelt.
    // This is used to prevent switching too often when the metric stays
    // for a while close to the threshold. It must be a value between
    // 0 and 1 [default = 0.05].
    double thresholdQBlockingBelt;

    // Max number of samples that can be violated (either accuracy or contract
    // violations) [default = 0].
    uint tolerableSamples;

    // Maximum size for the internal queues of the farm.
    // 0 corresponds to infinite size [default = 1].
    ulong qSize;

    // If different from zero:
    //      - If we have a PERF_* contract, we will search for a solution which
    //        is the conservativeValue% more performing than the requirement.
    //      - If we have a POWER_* contract, we will search for a solution
    //        which consumes conservativeValue% less power than the requirement.
    // This is done to amortize fluctuations. [default = 0.0]
    double conservativeValue;

    // A vector containing the number of cores not allowed to be used.
    // E.g. if it contains the number 3, then the runtime will
    // never use  3 cores. It can only be specified when
    // knobCoresEnabled is true [default = empty].
    std::vector<uint> disallowedNumCores;

    // If true, the manager will run on a physical core by itself.
    // If false, it will run on the same physical core used by the emitter
    // of the farm [default = false].
    bool isolateManager;

    // If true, computes the statistics about the cost
    // of the reconfigurations [default = false].
    bool statsReconfiguration;

    // Parameters for dataflow applications.
    DataflowParameters dataflow;

    // This is the path of the file used to signal begin/end of a
    // (region of interest). This file should be created when the
    // application enters the ROI and deleted when the application
    // exits the ROI. It can be used for example on PARSEC applications
    // in order to control only the region of interest of the application
    // excluding initialization and cleanup phases. In that case, the
    // file would be created in the roi_start hook and deleted in the
    // roi_end hook. If not specified or if empty string, the application
    // will be monitored and controlled throughout its entire execution
    // [default = ""].
    std::string roiFile;

    // The loggers objects. They will be called every samplingInterval
    // milliseconds to monitor the application [default = NULL].
    std::vector<Logger*> loggers;

    /**
     * Creates the nornir paramters.
     * @param communicator The communicator used to instantiate the other
     *        modules. If NULL, the modules will be created as local modules.
     */
    explicit Parameters(mammut::Communicator* const communicator = NULL);

    /**
     * Creates the nornir parameters.
     * @param paramFileName The name of the XML file containing the nornir
     *        parameters.
     * @param communicator The communicator used to instantiate the other
     *        modules. If NULL, the modules will be created as local modules.
     */
    Parameters(const std::string& paramFileName,
               mammut::Communicator* const communicator = NULL);

    /**
     * Destroyes these parameters.
     */
    ~Parameters();

    /**
     * Loads the parameters from a file.
     * @param paramFileName The name of the XML file containing the nornir
     *        parameters.
     */
    void load(const std::string& paramFileName);

    /**
     * Validates these parameters.
     * @return The validation result.
     */
    ParametersValidation validate();

    /**
     * Check if a specific knob is enabled.
     * ATTENTION: Can only be called after validate().
     * @param k The knob to check.
     * @return true if the specified knob is enabled, false otherwise.
     */
    bool isKnobEnabled(KnobType k) const;
};

/**
 * Checks if a requirement is a min/max one.
 * @param r The requirement.
 * @return true if the requirement is a min/max one.
 */
bool isMinMaxRequirement(double r);

/**
 * Checks if a requirement is a primary one.
 * @param r The requirement.
 * @return true if the requirement is a primary one.
 */
bool isPrimaryRequirement(double r);

}

#endif /* NORNIR_PARAMETERS_HPP_ */

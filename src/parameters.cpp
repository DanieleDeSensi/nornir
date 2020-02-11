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
 */

#include <nornir/parameters.hpp>
#include <nornir/utils.hpp>

#include <cstring>
#include <limits>
#include <mammut/mammut.hpp>

namespace nornir {

using namespace std;
using namespace mammut::cpufreq;
using namespace mammut::topology;
using namespace mammut::utils;
using mammut::Communicator;
using mammut::Mammut;
using mammut::utils::enumStrings;

#define CONFIGURATION_VERSION "1.1.0"
#define CONFPATH_LEN_MAX 512
#define CONFFILE_VERSION "/version.csv"
#define CONFFILE_ARCH "/archdata.xml"
#define CONFFILE_VOLTAGE "/voltage.csv"

void XmlTree::init(const std::string &content, const std::string &rootName) {
  rapidxml::xml_document<> xmlContent;
  _fileContentChars = new char[content.size() + 1];
  copy(content.begin(), content.end(), _fileContentChars);
  _fileContentChars[content.size()] = '\0';
  xmlContent.parse<0>(_fileContentChars);
  _root = xmlContent.first_node(rootName.c_str());
}

XmlTree::XmlTree(const string &fileName, const string &rootName) {
  ifstream file(fileName.c_str());
  if (!file) {
    throw runtime_error("Impossible to read xml file " + fileName);
  }

  string fileContent;
  file.seekg(0, ios::end);
  fileContent.reserve(file.tellg());
  file.seekg(0, ios::beg);
  fileContent.assign((istreambuf_iterator<char>(file)),
                     istreambuf_iterator<char>());
  init(fileContent, rootName);
}

/*XmlTree::XmlTree(const std::string& content, const std::string& rootName){
    init(content, rootName);
}*/

XmlTree::~XmlTree() {
  delete[] _fileContentChars;
}

rapidxml::xml_node<> *XmlTree::getNode(const char *valueName) const {
  rapidxml::xml_node<> *node = _root;
  string vName(valueName);
  size_t i = 0;
  vector<string> tokens = mammut::utils::split(vName, '.');
  if (!tokens.size()) {
    return NULL;
  }
  while (i < tokens.size() && node != NULL) {
    node = node->first_node(tokens.at(i).c_str());
    i++;
  }
  return node;
}

void XmlTree::getBool(const char *valueName, bool &value) {
  rapidxml::xml_node<> *node = getNode(valueName);
  if (node) {
    value = string(node->value()).compare("true") ? false : true;
  }
}

void XmlTree::getInt(const char *valueName, int &value) {
  rapidxml::xml_node<> *node = getNode(valueName);
  if (node) {
    value = stringToInt(node->value());
  }
}

void XmlTree::getUint(const char *valueName, uint &value) {
  rapidxml::xml_node<> *node = getNode(valueName);
  if (node) {
    value = stringToUint(node->value());
  }
}

void XmlTree::getUlong(const char *valueName, ulong &value) {
  rapidxml::xml_node<> *node = getNode(valueName);
  if (node) {
    value = stringToUlong(node->value());
  }
}

void XmlTree::getDouble(const char *valueName, double &value) {
  rapidxml::xml_node<> *node = getNode(valueName);
  if (node) {
    value = stringToDouble(node->value());
  }
}

void XmlTree::getDoubleOrMin(const char *valueName, double &value) {
  rapidxml::xml_node<> *node = getNode(valueName);
  if (node) {
    if (std::string(node->value()).compare("MIN") == 0) {
      value = NORNIR_REQUIREMENT_MIN;
    } else {
      value = stringToDouble(node->value());
    }
  }
}

void XmlTree::getDoubleOrMax(const char *valueName, double &value) {
  rapidxml::xml_node<> *node = getNode(valueName);
  if (node) {
    if (std::string(node->value()).compare("MAX") == 0) {
      value = NORNIR_REQUIREMENT_MAX;
    } else {
      value = stringToDouble(node->value());
    }
  }
}

void XmlTree::getString(const char *valueName, string &value) {
  rapidxml::xml_node<> *node = getNode(valueName);
  if (node) {
    value = node->value();
  }
}

void XmlTree::getArrayUint(const char *valueName, std::vector<uint> &value) {
  rapidxml::xml_node<> *node = getNode(valueName);
  value.clear();
  if (node) {
    vector<string> strValues = mammut::utils::split(node->value(), ':');
    for (size_t i = 0; i < strValues.size(); i++) {
      value.push_back(stringToUint(strValues.at(i)));
    }
  }
}

template <typename T>
void XmlTree::getArrayEnums(const char *valueName, std::vector<T> &values) {
  rapidxml::xml_node<> *node = getNode(valueName);
  values.clear();
  if (node) {
    vector<string> strValues = mammut::utils::split(node->value(), ':');
    for (size_t i = 0; i < strValues.size(); i++) {
      T value;
      stringstream line(strValues.at(i));
      line >> enumFromStringInternal(value);
      values.push_back(value);
    }
  }
}

template <typename T> void XmlTree::getEnum(const char *valueName, T &value) {
  rapidxml::xml_node<> *node = getNode(valueName);
  if (node) {
    stringstream line(node->value());
    line >> enumFromStringInternal(value);
  }
}

void ArchData::loadXml(const string &archFileName) {
  XmlTree xc(archFileName, "archData");
  SETVALUE(xc, Double, ticksPerNs);
  SETVALUE(xc, Double, idlePower);
  SETVALUE(xc, Double, monitoringCost);
}

Requirements::Requirements() {
  throughput = NORNIR_REQUIREMENT_UNDEF;
  powerConsumption = NORNIR_REQUIREMENT_UNDEF;
  minUtilization = NORNIR_REQUIREMENT_UNDEF;
  maxUtilization = NORNIR_REQUIREMENT_UNDEF;
  executionTime = NORNIR_REQUIREMENT_UNDEF;
  energy = NORNIR_REQUIREMENT_UNDEF;
  expectedTasksNumber = NORNIR_REQUIREMENT_UNDEF;
  latency = NORNIR_REQUIREMENT_UNDEF;
}

bool Requirements::anySpecified() const {
  return throughput != NORNIR_REQUIREMENT_UNDEF ||
         powerConsumption != NORNIR_REQUIREMENT_UNDEF ||
         minUtilization != NORNIR_REQUIREMENT_UNDEF ||
         maxUtilization != NORNIR_REQUIREMENT_UNDEF ||
         executionTime != NORNIR_REQUIREMENT_UNDEF ||
         energy != NORNIR_REQUIREMENT_UNDEF ||
         latency != NORNIR_REQUIREMENT_UNDEF;
}

void Parameters::setDefault() {
  triggerQBlocking = TRIGGER_Q_BLOCKING_NO;
  strategyUnusedVirtualCores = STRATEGY_UNUSED_VC_SAME;
  strategySelection = STRATEGY_SELECTION_LEARNING;
  strategyPredictionPerformance = STRATEGY_PREDICTION_PERFORMANCE_USL;
  strategyPredictionPower = STRATEGY_PREDICTION_POWER_LINEAR;
  strategyExploration = STRATEGY_EXPLORATION_HALTON;
  strategySmoothing = STRATEGY_SMOOTHING_EXPONENTIAL;
  strategyPolling = STRATEGY_POLLING_SLEEP_SMALL;
  strategyPersistence = STRATEGY_PERSISTENCE_SAMPLES;
  strategyPhaseDetection = STRATEGY_PHASE_DETECTION_NONE;
  knobCoresEnabled = true;
  knobMappingEnabled = true;
  knobFrequencyEnabled = true;
  knobClkModEnabled = false;
  knobHyperthreadingEnabled = false;
  knobHyperthreadingFixedValue = 0;
  knobPforChunkEnabled = false;
  activeThreads = 0;
  useConcurrencyThrottling = true;
  fastReconfiguration = true;
  migrateCollector = false;
  smoothingFactor = 0;
  persistenceValue = 0;
  cooldownPeriod = 200;
  samplingIntervalCalibration = 100;
  samplingIntervalSteady = 1000;
  steadyThreshold = 4;
  minTasksPerSample = 0;
  synchronousWorkers = false;
  maxCalibrationTime = 0;
  maxCalibrationSteps = 0;
  maxPerformancePredictionError = 10.0;
  maxPowerPredictionError = 5.0;
  regressionAging = 0;
  maxMonitoringOverhead = 1.0;
  clockModulationEmulated = true;
  clockModulationMin = 1.0;
  clockModulationResolution = 2.0;
  clockModulationInterval = 1000;
  thresholdQBlocking = -1;
  thresholdQBlockingBelt = 0.05;
  tolerableSamples = 0;
  qSize = 1;
  conservativeValue = 0;
  isolateManager = false;
  statsReconfiguration = false;
  nelderMeadRange = 2;
  fixedPinning = false;
  powerDomain = mammut::energy::COUNTER_CPUS;
  perPidLog = false;
  roiFile = "";

  leo.applicationName = "";
  leo.namesData = "";
  leo.throughputData = "";
  leo.powerData = "";
  leo.numSamples = 20;

  dataflow.orderedProcessing = false;
  dataflow.orderedOutput = false;
  dataflow.maxGraphs = 1000;
  dataflow.maxInterpreters = 0;

  /** Retrieving global configuration files. **/
  bool found = false;
  for (string confHome : getXdgConfigDirs()) {
    string confFileArch = confHome + string(CONFFILE_ARCH);
    string confFileVoltage = confHome + string(CONFFILE_VOLTAGE);
    string confFileVersion = confHome + string(CONFFILE_VERSION);

    if (existsFile(confFileArch) && existsFile(confFileVoltage) &&
        existsFile(confFileVersion) &&
        readFirstLineFromFile(confFileVersion).compare(CONFIGURATION_VERSION) ==
            0) {
      archData.loadXml(confFileArch);
      archData.voltageTable = new VoltageTable();
      loadVoltageTable(*static_cast<VoltageTable *>(archData.voltageTable),
                       confFileVoltage);
      found = true;
      break;
    }
  }

  if (!found) {
    throw runtime_error("Impossible to find configuration files. Please run "
                        "'make microbench' from the nornir build folder.");
  }
}

uint Parameters::getLowOverheadSamplingInterval() const {
  // TODO Se questo sampling interval è molto minore della latenza
  // media di un task potrei settare il sampling interval alla
  // latenza media di un task.

  // We are setting the sampling interval such that the latency
  // due to the communication between manager and node is the
  // maxMonitoringOverhead% of the interval.
  double msMonitoringCost =
      (archData.monitoringCost / archData.ticksPerNs * 0.000001);
  return ceil(msMonitoringCost * (100.0 - maxMonitoringOverhead));
}

/**
 * Sets the default values for parameters that depends
 * from others.
 */
void Parameters::setDefaultPost() {
  if (!samplingIntervalCalibration) {
    samplingIntervalCalibration = getLowOverheadSamplingInterval();
  }

  if (!samplingIntervalSteady) {
    samplingIntervalSteady = getLowOverheadSamplingInterval();
  }

  if (!smoothingFactor) {
    switch (strategySmoothing) {
    case STRATEGY_SMOOTHING_MOVING_AVERAGE: {
      smoothingFactor = 10;
    } break;
    case STRATEGY_SMOOTHING_EXPONENTIAL: {
      smoothingFactor = 0.5;
    } break;
    }
  }

  if (!persistenceValue) {
    switch (strategyPersistence) {
    case STRATEGY_PERSISTENCE_SAMPLES: {
      persistenceValue = 1;
    } break;
    case STRATEGY_PERSISTENCE_VARIATION: {
      persistenceValue = 5;
    } break;
    }
  }
}

static vector<Frequency> getAvailableFrequencies(const Mammut &mammut) {
  vector<Frequency> frequencies;
  CpuFreq *cpuFreq = mammut.getInstanceCpuFreq();
  if (cpuFreq) {
    cpuFreq->removeTurboFrequencies();
    vector<Domain *> fDomains = cpuFreq->getDomains();
    if (fDomains.size()) {
      frequencies = fDomains.front()->getAvailableFrequencies();
    }
  }
  return frequencies;
}

bool Parameters::isUnusedVcOffAvailable() {
  vector<VirtualCore *> vc = mammut.getInstanceTopology()->getVirtualCores();
  for (size_t i = 0; i < vc.size(); i++) {
    if (vc.at(i)->isHotPluggable()) {
      return true;
    }
  }
  return false;
}

bool Parameters::isFrequencySettable() {
  vector<Frequency> frequencies = getAvailableFrequencies(mammut);
  return mammut.getInstanceCpuFreq()->isGovernorAvailable(GOVERNOR_USERSPACE) &&
         frequencies.size();
}

bool Parameters::isLowestFrequencySettable() {
  return mammut.getInstanceCpuFreq()->isGovernorAvailable(GOVERNOR_POWERSAVE) ||
         isFrequencySettable();
}

bool Parameters::isHighestFrequencySettable() {
  return mammut.getInstanceCpuFreq()->isGovernorAvailable(
             GOVERNOR_PERFORMANCE) ||
         isFrequencySettable();
}

ParametersValidation
Parameters::validateUnusedVc(StrategyUnusedVirtualCores &s) {
  switch (s) {
  case STRATEGY_UNUSED_VC_OFF: {
    if (!isUnusedVcOffAvailable()) {
      return VALIDATION_UNUSED_VC_NO_OFF;
    }
  } break;
  case STRATEGY_UNUSED_VC_LOWEST_FREQUENCY: {
    if (!isLowestFrequencySettable()) {
      return VALIDATION_UNUSED_VC_NO_FREQUENCIES;
    }
  } break;
  case STRATEGY_UNUSED_VC_AUTO: {
    if (isUnusedVcOffAvailable()) {
      s = STRATEGY_UNUSED_VC_OFF;
    } else if (isLowestFrequencySettable()) {
      s = STRATEGY_UNUSED_VC_LOWEST_FREQUENCY;
    } else {
      s = STRATEGY_UNUSED_VC_SAME;
    }
  } break;
  default:
    break;
  }
  return VALIDATION_OK;
}

ParametersValidation Parameters::validateKnobFrequencies() {
  vector<Frequency> availableFrequencies = getAvailableFrequencies(mammut);

  if (knobFrequencyEnabled &&
      (!mammut.getInstanceCpuFreq()->isGovernorAvailable(GOVERNOR_USERSPACE) ||
       availableFrequencies.empty())) {
    return VALIDATION_NO_MANUAL_DVFS;
  }

#if defined(__x86_64__)
  vector<VirtualCore *> virtualCores;
  virtualCores = mammut.getInstanceTopology()->getVirtualCores();
  for (size_t i = 0; i < virtualCores.size(); i++) {
    if (!virtualCores.at(i)->hasFlag("constant_tsc")) {
      return VALIDATION_NO_CONSTANT_TSC;
    }
  }
#endif

  if (fastReconfiguration &&
      (!isHighestFrequencySettable() ||
       strategyUnusedVirtualCores == STRATEGY_UNUSED_VC_NONE)) {
    fastReconfiguration = false;
  }

  return VALIDATION_OK;
}

ParametersValidation Parameters::validateKnobClkMod() {
  if (knobClkModEnabled && !clockModulationEmulated) {
    vector<Cpu *> cpus = mammut.getInstanceTopology()->getCpus();
    for (Cpu *c : cpus) {
      if (!c->hasClockModulation()) {
        return VALIDATION_NO_CLKMOD;
      }
    }
  }
  return VALIDATION_OK;
}

ParametersValidation Parameters::validateTriggers() {
  if (triggerQBlocking == TRIGGER_Q_BLOCKING_AUTO &&
      (thresholdQBlocking == -1 || thresholdQBlockingBelt < 0 ||
       thresholdQBlockingBelt > 1)) {
    return VALIDATION_BLOCKING_PARAMETERS;
  }
  return VALIDATION_OK;
}

ParametersValidation Parameters::validateRequirements() {
  if (requirements.minUtilization != NORNIR_REQUIREMENT_UNDEF ||
      requirements.maxUtilization != NORNIR_REQUIREMENT_UNDEF) {
    if (requirements.minUtilization == NORNIR_REQUIREMENT_UNDEF) {
      requirements.minUtilization = 0.1;
    }
    if (requirements.maxUtilization == NORNIR_REQUIREMENT_UNDEF) {
      requirements.maxUtilization = 100;
    }
    if (requirements.minUtilization > requirements.maxUtilization ||
        requirements.minUtilization < 0 || requirements.maxUtilization > 100) {
      return VALIDATION_WRONG_REQUIREMENT;
    }
  }
  if (requirements.throughput != NORNIR_REQUIREMENT_UNDEF &&
      requirements.throughput < 0) {
    return VALIDATION_WRONG_REQUIREMENT;
  }
  if (requirements.executionTime != NORNIR_REQUIREMENT_UNDEF &&
      requirements.executionTime != NORNIR_REQUIREMENT_MIN &&
      requirements.executionTime < 0) {
    return VALIDATION_WRONG_REQUIREMENT;
  }
  if (requirements.energy != NORNIR_REQUIREMENT_UNDEF &&
      requirements.energy != NORNIR_REQUIREMENT_MIN &&
      requirements.energy < 0) {
    return VALIDATION_WRONG_REQUIREMENT;
  }
  if (requirements.powerConsumption != NORNIR_REQUIREMENT_UNDEF &&
      requirements.powerConsumption != NORNIR_REQUIREMENT_MIN &&
      requirements.powerConsumption < 0) {
    return VALIDATION_WRONG_REQUIREMENT;
  }
  if (requirements.latency != NORNIR_REQUIREMENT_UNDEF &&
      requirements.latency != NORNIR_REQUIREMENT_MIN &&
      requirements.latency < 0) {
    return VALIDATION_WRONG_REQUIREMENT;
  }
  if ((requirements.executionTime != NORNIR_REQUIREMENT_UNDEF ||
       requirements.energy != NORNIR_REQUIREMENT_UNDEF) &&
      requirements.expectedTasksNumber == NORNIR_REQUIREMENT_UNDEF) {
    return VALIDATION_WRONG_REQUIREMENT;
  }
  if (isPrimaryRequirement(requirements.powerConsumption)) {
    if (strategySelection == STRATEGY_SELECTION_ANALYTICAL ||
        strategySelection == STRATEGY_SELECTION_LIMARTINEZ) {
      return VALIDATION_WRONG_REQUIREMENT;
    }
  }
  if(strategySelection == STRATEGY_SELECTION_RAPL &&
     !isPrimaryRequirement(requirements.powerConsumption)) {
    return VALIDATION_WRONG_REQUIREMENT;
  }
  
  if (requirements.energy != NORNIR_REQUIREMENT_UNDEF &&
      strategySelection != STRATEGY_SELECTION_ANALYTICAL_FULL &&
      strategySelection != STRATEGY_SELECTION_LEARNING &&
      strategySelection != STRATEGY_SELECTION_HMP_NELDERMEAD) {
    return VALIDATION_WRONG_REQUIREMENT;
  }
  if (maxCalibrationTime == 0 && maxCalibrationSteps &&
      (maxPerformancePredictionError <= 0 ||
       maxPerformancePredictionError > 100.0 || maxPowerPredictionError <= 0 ||
       maxPowerPredictionError > 100.0)) {
    return VALIDATION_WRONG_REQUIREMENT;
  }

  uint maxMinRequirements = 0;
  if (requirements.throughput == NORNIR_REQUIREMENT_MAX) {
    ++maxMinRequirements;
  }
  if (requirements.executionTime == NORNIR_REQUIREMENT_MIN) {
    ++maxMinRequirements;
  }
  if (requirements.energy == NORNIR_REQUIREMENT_MIN) {
    ++maxMinRequirements;
  }
  if (requirements.latency == NORNIR_REQUIREMENT_MIN) {
    ++maxMinRequirements;
  }
  if (requirements.powerConsumption == NORNIR_REQUIREMENT_MIN) {
    ++maxMinRequirements;
  }

  // TODO: Remove Utilization requirements (can be obtained through
  // throughput + tolerance/conservative
  // At most 1 min/max requirement can be specified.
  // TODO: At most one or exactly one?
  if (maxMinRequirements > 1) {
    return VALIDATION_WRONG_REQUIREMENT;
  }

  return VALIDATION_OK;
}

// clang-format off
ParametersValidation Parameters::validateSelector() {
  bool knobsSupportSelector[STRATEGY_SELECTION_NUM][KNOB_NUM];
  if (strategySelection == STRATEGY_SELECTION_NUM) {
    return VALIDATION_NO;
  }

  // MANUAL CLI
  knobsSupportSelector[STRATEGY_SELECTION_MANUAL_CLI][KNOB_VIRTUAL_CORES] =
      true;
  knobsSupportSelector[STRATEGY_SELECTION_MANUAL_CLI][KNOB_FREQUENCY] = true;
  knobsSupportSelector[STRATEGY_SELECTION_MANUAL_CLI][KNOB_MAPPING] = true;
  knobsSupportSelector[STRATEGY_SELECTION_MANUAL_CLI][KNOB_HYPERTHREADING] =
      true;
  knobsSupportSelector[STRATEGY_SELECTION_MANUAL_CLI][KNOB_CLKMOD] = true;
  knobsSupportSelector[STRATEGY_SELECTION_MANUAL_CLI][KNOB_PFOR_CHUNK] = false;

  // MANUAL WEB
  knobsSupportSelector[STRATEGY_SELECTION_MANUAL_WEB][KNOB_VIRTUAL_CORES] =
      true;
  knobsSupportSelector[STRATEGY_SELECTION_MANUAL_WEB][KNOB_FREQUENCY] = true;
  knobsSupportSelector[STRATEGY_SELECTION_MANUAL_WEB][KNOB_MAPPING] = false;
  knobsSupportSelector[STRATEGY_SELECTION_MANUAL_WEB][KNOB_HYPERTHREADING] =
      false;
  knobsSupportSelector[STRATEGY_SELECTION_MANUAL_WEB][KNOB_CLKMOD] = false;
  knobsSupportSelector[STRATEGY_SELECTION_MANUAL_WEB][KNOB_PFOR_CHUNK] = false;

  // ANALYTICAL
  knobsSupportSelector[STRATEGY_SELECTION_ANALYTICAL][KNOB_VIRTUAL_CORES] =
      true;
  knobsSupportSelector[STRATEGY_SELECTION_ANALYTICAL][KNOB_FREQUENCY] = true;
  knobsSupportSelector[STRATEGY_SELECTION_ANALYTICAL][KNOB_MAPPING] = false;
  knobsSupportSelector[STRATEGY_SELECTION_ANALYTICAL][KNOB_HYPERTHREADING] =
      false;
  knobsSupportSelector[STRATEGY_SELECTION_ANALYTICAL][KNOB_CLKMOD] = true;
  knobsSupportSelector[STRATEGY_SELECTION_ANALYTICAL][KNOB_PFOR_CHUNK] = false;

  // ANALYTICAL_FULL
  knobsSupportSelector[STRATEGY_SELECTION_ANALYTICAL_FULL][KNOB_VIRTUAL_CORES] =
      true;
  knobsSupportSelector[STRATEGY_SELECTION_ANALYTICAL_FULL][KNOB_FREQUENCY] =
      true;
  knobsSupportSelector[STRATEGY_SELECTION_ANALYTICAL_FULL][KNOB_MAPPING] =
      false;
  knobsSupportSelector[STRATEGY_SELECTION_ANALYTICAL_FULL]
                      [KNOB_HYPERTHREADING] = false;
  knobsSupportSelector[STRATEGY_SELECTION_ANALYTICAL_FULL][KNOB_CLKMOD] = true;
  knobsSupportSelector[STRATEGY_SELECTION_ANALYTICAL_FULL][KNOB_PFOR_CHUNK] = false;

  // FULLSEARCH
  knobsSupportSelector[STRATEGY_SELECTION_FULLSEARCH][KNOB_VIRTUAL_CORES] =
      true;
  knobsSupportSelector[STRATEGY_SELECTION_FULLSEARCH][KNOB_FREQUENCY] = true;
  knobsSupportSelector[STRATEGY_SELECTION_FULLSEARCH][KNOB_MAPPING] = true;
  knobsSupportSelector[STRATEGY_SELECTION_FULLSEARCH][KNOB_HYPERTHREADING] =
      true;
  knobsSupportSelector[STRATEGY_SELECTION_FULLSEARCH][KNOB_CLKMOD] = true;
  knobsSupportSelector[STRATEGY_SELECTION_FULLSEARCH][KNOB_PFOR_CHUNK] = false;

  // For learning we do not check since it depends from the predictors choice.
  // (we will check in validatePredictors())

  // LIMARTINEZ
  knobsSupportSelector[STRATEGY_SELECTION_LIMARTINEZ][KNOB_VIRTUAL_CORES] =
      true;
  knobsSupportSelector[STRATEGY_SELECTION_LIMARTINEZ][KNOB_FREQUENCY] = true;
  knobsSupportSelector[STRATEGY_SELECTION_LIMARTINEZ][KNOB_MAPPING] = false;
  knobsSupportSelector[STRATEGY_SELECTION_LIMARTINEZ][KNOB_HYPERTHREADING] =
      false;
  knobsSupportSelector[STRATEGY_SELECTION_LIMARTINEZ][KNOB_CLKMOD] = false;
  knobsSupportSelector[STRATEGY_SELECTION_LIMARTINEZ][KNOB_PFOR_CHUNK] = false;

  // LEO
  knobsSupportSelector[STRATEGY_SELECTION_LEO][KNOB_VIRTUAL_CORES] = true;
  knobsSupportSelector[STRATEGY_SELECTION_LEO][KNOB_FREQUENCY] = true;
  knobsSupportSelector[STRATEGY_SELECTION_LEO][KNOB_MAPPING] = false;
  knobsSupportSelector[STRATEGY_SELECTION_LEO][KNOB_HYPERTHREADING] = false;
  knobsSupportSelector[STRATEGY_SELECTION_LEO][KNOB_CLKMOD] = false;
  knobsSupportSelector[STRATEGY_SELECTION_LEO][KNOB_PFOR_CHUNK] = false;

  if (strategySelection == STRATEGY_SELECTION_LEO &&
      (leo.throughputData.compare("") == 0 || leo.powerData.compare("") == 0 ||
       leo.applicationName.compare("") == 0 ||
       leo.namesData.compare("") == 0)) {
    return VALIDATION_NO_LEO_PARAMETERS;
  }

  // HMP - Local Search
  knobsSupportSelector[STRATEGY_SELECTION_HMP_NELDERMEAD][KNOB_VIRTUAL_CORES] =
      true;
  knobsSupportSelector[STRATEGY_SELECTION_HMP_NELDERMEAD][KNOB_FREQUENCY] =
      true;
  knobsSupportSelector[STRATEGY_SELECTION_HMP_NELDERMEAD][KNOB_MAPPING] = false;
  knobsSupportSelector[STRATEGY_SELECTION_HMP_NELDERMEAD][KNOB_HYPERTHREADING] =
      false;
  knobsSupportSelector[STRATEGY_SELECTION_HMP_NELDERMEAD][KNOB_CLKMOD] = false;
  knobsSupportSelector[STRATEGY_SELECTION_HMP_NELDERMEAD][KNOB_PFOR_CHUNK] = false;

  // RAPL
  knobsSupportSelector[STRATEGY_SELECTION_RAPL][KNOB_VIRTUAL_CORES] = false;
  knobsSupportSelector[STRATEGY_SELECTION_RAPL][KNOB_FREQUENCY] = false;
  knobsSupportSelector[STRATEGY_SELECTION_RAPL][KNOB_MAPPING] = false;
  knobsSupportSelector[STRATEGY_SELECTION_RAPL][KNOB_HYPERTHREADING] = false;
  knobsSupportSelector[STRATEGY_SELECTION_RAPL][KNOB_CLKMOD] = false;
  knobsSupportSelector[STRATEGY_SELECTION_RAPL][KNOB_PFOR_CHUNK] = false;

  // PFOR_CHUNK
  knobsSupportSelector[STRATEGY_SELECTION_PFOR_CHUNK][KNOB_VIRTUAL_CORES] = false;
  knobsSupportSelector[STRATEGY_SELECTION_PFOR_CHUNK][KNOB_FREQUENCY] = false;
  knobsSupportSelector[STRATEGY_SELECTION_PFOR_CHUNK][KNOB_MAPPING] = false;
  knobsSupportSelector[STRATEGY_SELECTION_PFOR_CHUNK][KNOB_HYPERTHREADING] = false;
  knobsSupportSelector[STRATEGY_SELECTION_PFOR_CHUNK][KNOB_CLKMOD] = false;
  knobsSupportSelector[STRATEGY_SELECTION_PFOR_CHUNK][KNOB_PFOR_CHUNK] = true;

  if (strategySelection == STRATEGY_SELECTION_HMP_NELDERMEAD &&
      (firstConfiguration.virtualCores.empty() ||
       firstConfiguration.frequency.empty())) {
    return VALIDATION_NO_FIRSTCONF_PARAMETERS;
  }

  if (!(firstConfiguration.frequency.size() ==
        firstConfiguration.virtualCores.size())) {
    return VALIDATION_NO;
  }

  if (strategySelection != STRATEGY_SELECTION_LEARNING) {
    // Check if the knob enabled can be managed by the selector specified.
    for (size_t i = 0; i < KNOB_NUM; i++) {
      if (_knobEnabled[i] && !knobsSupportSelector[strategySelection][i]) {
        return VALIDATION_UNSUPPORTED_KNOBS;
      }
    }
  } else {
    /*********************************************/
    /*            Validate predictors.           */
    /*********************************************/
    bool knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_NUM][KNOB_NUM];
    bool knobsSupportPower[STRATEGY_PREDICTION_POWER_NUM][KNOB_NUM];

    if (strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_NUM ||
        strategyPredictionPower == STRATEGY_PREDICTION_POWER_NUM) {
      return VALIDATION_NO;
    }
    /******************************************/
    /*          Performance models.           */
    /******************************************/
    // AMDAHL
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_AMDAHL][KNOB_VIRTUAL_CORES] = true;
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_AMDAHL][KNOB_FREQUENCY] = true;
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_AMDAHL][KNOB_MAPPING] = true;
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_AMDAHL][KNOB_HYPERTHREADING] = false;
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_AMDAHL][KNOB_CLKMOD] = true;
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_AMDAHL][KNOB_PFOR_CHUNK] = false;
    // LEO
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_LEO][KNOB_VIRTUAL_CORES] = true;
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_LEO][KNOB_FREQUENCY] = true;
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_LEO][KNOB_MAPPING] = false;
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_LEO][KNOB_HYPERTHREADING] = false;
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_LEO][KNOB_CLKMOD] = false;
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_LEO][KNOB_PFOR_CHUNK] = false;
    // USL
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_USL][KNOB_VIRTUAL_CORES] = true;
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_USL][KNOB_FREQUENCY] = true;
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_USL][KNOB_MAPPING] = true;
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_USL][KNOB_HYPERTHREADING] = false;
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_USL][KNOB_CLKMOD] = true;
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_USL][KNOB_PFOR_CHUNK] = false;
    // USLP
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_USLP][KNOB_VIRTUAL_CORES] = true;
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_USLP][KNOB_FREQUENCY] = true;
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_USLP][KNOB_MAPPING] = true;
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_USLP][KNOB_HYPERTHREADING] = false;
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_USLP][KNOB_CLKMOD] = true;
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_USLP][KNOB_PFOR_CHUNK] = false;
    // SMT
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_SMT][KNOB_VIRTUAL_CORES] = true;
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_SMT][KNOB_FREQUENCY] = true;
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_SMT][KNOB_MAPPING] = false;
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_SMT][KNOB_HYPERTHREADING] = true;
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_SMT][KNOB_CLKMOD] = false;
    knobsSupportPerformance[STRATEGY_PREDICTION_PERFORMANCE_SMT][KNOB_PFOR_CHUNK] = false;

    /******************************************/
    /*              Power models.             */
    /******************************************/
    // LINEAR
    knobsSupportPower[STRATEGY_PREDICTION_POWER_LINEAR][KNOB_VIRTUAL_CORES] = true;
    knobsSupportPower[STRATEGY_PREDICTION_POWER_LINEAR][KNOB_FREQUENCY] = true;
    knobsSupportPower[STRATEGY_PREDICTION_POWER_LINEAR][KNOB_MAPPING] = true;
    knobsSupportPower[STRATEGY_PREDICTION_POWER_LINEAR][KNOB_HYPERTHREADING] = false;
    knobsSupportPower[STRATEGY_PREDICTION_POWER_LINEAR][KNOB_CLKMOD] = true;
    knobsSupportPower[STRATEGY_PREDICTION_POWER_LINEAR][KNOB_PFOR_CHUNK] = false;
    // LEO
    knobsSupportPower[STRATEGY_PREDICTION_POWER_LEO][KNOB_VIRTUAL_CORES] = true;
    knobsSupportPower[STRATEGY_PREDICTION_POWER_LEO][KNOB_FREQUENCY] = true;
    knobsSupportPower[STRATEGY_PREDICTION_POWER_LEO][KNOB_MAPPING] = false;
    knobsSupportPower[STRATEGY_PREDICTION_POWER_LEO][KNOB_HYPERTHREADING] = false;
    knobsSupportPower[STRATEGY_PREDICTION_POWER_LEO][KNOB_CLKMOD] = false;
    knobsSupportPower[STRATEGY_PREDICTION_POWER_LEO][KNOB_PFOR_CHUNK] = false;
    // SMT
    knobsSupportPower[STRATEGY_PREDICTION_POWER_SMT][KNOB_VIRTUAL_CORES] = true;
    knobsSupportPower[STRATEGY_PREDICTION_POWER_SMT][KNOB_FREQUENCY] = true;
    knobsSupportPower[STRATEGY_PREDICTION_POWER_SMT][KNOB_MAPPING] = false;
    knobsSupportPower[STRATEGY_PREDICTION_POWER_SMT][KNOB_HYPERTHREADING] = true;
    knobsSupportPower[STRATEGY_PREDICTION_POWER_SMT][KNOB_CLKMOD] = false;
    knobsSupportPower[STRATEGY_PREDICTION_POWER_SMT][KNOB_PFOR_CHUNK] = false;

    // Check if the knob enabled can be managed by the predictors specified.
    for (size_t i = 0; i < KNOB_NUM; i++) {
      if (_knobEnabled[i] &&
          (!knobsSupportPerformance[strategyPredictionPerformance][i] ||
           !knobsSupportPower[strategyPredictionPower][i])) {
        return VALIDATION_UNSUPPORTED_KNOBS;
      }
    }

    if (strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_LEO &&
        (leo.throughputData.compare("") == 0 ||
         leo.powerData.compare("") == 0 ||
         leo.applicationName.compare("") == 0 ||
         leo.namesData.compare("") == 0)) {
      return VALIDATION_NO_LEO_PARAMETERS;
    }

    // Currently, USL predictors only works with low discrepancy explorators.
    // TODO: This is because the additional exploration points at the moment
    // can only be added to the low discrepancy generators.
    if ((strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_USL ||
         strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_USLP ||
         strategyPredictionPerformance == STRATEGY_PREDICTION_PERFORMANCE_SMT ||
         strategyPredictionPower == STRATEGY_PREDICTION_POWER_SMT) &&
        (strategyExploration != STRATEGY_EXPLORATION_HALTON &&
         strategyExploration != STRATEGY_EXPLORATION_HALTON_REVERSE &&
         strategyExploration != STRATEGY_EXPLORATION_RANDOM &&
         strategyExploration != STRATEGY_EXPLORATION_SOBOL)) {
      return VALIDATION_NO;
    }
  }

  return VALIDATION_OK;
}

template <> char const *enumStrings<LoggerType>::data[] = {
  "FILE",
  "GRAPHITE"
};

template <> char const *enumStrings<TriggerConfQBlocking>::data[] = {
  "NO",
  "YES",
  "AUTO"
};

template <> char const *enumStrings<StrategyUnusedVirtualCores>::data[] = {
  "AUTO",
  "SAME",
  "NONE",
  "LOWEST_FREQUENCY",
  "OFF"
};

template <> char const *enumStrings<StrategySelection>::data[] = {
  "MANUAL_CLI",
  "MANUAL_WEB",
  "LEARNING",
  "ANALYTICAL",
  "ANALYTICAL_FULL",
  "FULLSEARCH",
  "LIMARTINEZ",
  "LEO",
  "HMP_NELDERMEAD",
  "RAPL",
  "PFOR_CHUNK",
  "NUM" // <- Must always be the last
};

template <> char const *enumStrings<StrategyPredictionPerformance>::data[] = {
  "AMDAHL",
  "USL",
  "USLP",
  "LEO",
  "SMT",
  "NUM" // <- Must always be the last
};

template <> char const *enumStrings<StrategyPredictionPower>::data[] = {
  "LINEAR",
  "LEO",
  "SMT",
  "NUM" // <- Must always be the last
};

template <> char const *enumStrings<StrategyExploration>::data[] = {
  "RANDOM",
  "NIEDERREITER",
  "SOBOL",
  "HALTON",
  "HALTON_REVERSE"
};

template <> char const *enumStrings<StrategySmoothing>::data[] = {
  "MOVING_AVERAGE",
  "EXPONENTIAL"
};

template <> char const *enumStrings<StrategySmoothingFactor>::data[] = {
  "CONST",
  "DYNAMIC"
};

template <> char const *enumStrings<StrategyPolling>::data[] = {
  "SPINNING",
  "PAUSE",
  "SLEEP_SMALL",
  "SLEEP_LATENCY"
};

template <> char const *enumStrings<StrategyPersistence>::data[] = {
  "SAMPLES",
  "VARIATION"
};

template <> char const *enumStrings<StrategyPhaseDetection>::data[] = {
  "NONE",
  "TRIVIAL"
};

template <> char const *enumStrings<mammut::energy::CounterType>::data[] = {
  "CPUS",
  "MEMORY",
  "PLUG"
};
// clang-format on

void Parameters::loadXml(const string &paramFileName) {
  XmlTree xt(paramFileName, "nornirParameters");

  SETVALUE(xt, DoubleOrMax, requirements.throughput);
  SETVALUE(xt, DoubleOrMin, requirements.powerConsumption);
  SETVALUE(xt, Double, requirements.minUtilization);
  SETVALUE(xt, Double, requirements.maxUtilization);
  SETVALUE(xt, Double, requirements.expectedTasksNumber);
  SETVALUE(xt, DoubleOrMin, requirements.executionTime);
  SETVALUE(xt, DoubleOrMin, requirements.energy);
  SETVALUE(xt, DoubleOrMin, requirements.latency);

  SETVALUE(xt, Enum, strategyUnusedVirtualCores);
  SETVALUE(xt, Enum, strategySelection);
  SETVALUE(xt, Enum, strategyPredictionPerformance);
  SETVALUE(xt, Enum, strategyPredictionPower);
  SETVALUE(xt, Enum, strategyExploration);
  SETVALUE(xt, Enum, strategySmoothing);
  SETVALUE(xt, Enum, strategyPolling);
  SETVALUE(xt, Enum, strategyPersistence);
  SETVALUE(xt, Enum, strategyPhaseDetection);
  SETVALUE(xt, Enum, triggerQBlocking);
  SETVALUE(xt, Bool, knobCoresEnabled);
  SETVALUE(xt, Bool, knobMappingEnabled);
  SETVALUE(xt, Bool, knobFrequencyEnabled);
  SETVALUE(xt, Bool, knobClkModEnabled);
  SETVALUE(xt, Bool, knobHyperthreadingEnabled);
  SETVALUE(xt, Double, knobHyperthreadingFixedValue);
  SETVALUE(xt, Bool, knobPforChunkEnabled);

  SETVALUE(xt, Uint, activeThreads);
  SETVALUE(xt, Bool, useConcurrencyThrottling);
  SETVALUE(xt, Bool, fastReconfiguration);
  SETVALUE(xt, Double, smoothingFactor);
  SETVALUE(xt, Double, persistenceValue);
  SETVALUE(xt, Double, cooldownPeriod);
  SETVALUE(xt, Uint, samplingIntervalCalibration);
  SETVALUE(xt, Uint, samplingIntervalSteady);
  SETVALUE(xt, Uint, steadyThreshold);
  SETVALUE(xt, Uint, minTasksPerSample);
  SETVALUE(xt, Bool, migrateCollector);
  SETVALUE(xt, Bool, synchronousWorkers);
  SETVALUE(xt, Double, maxCalibrationTime);
  SETVALUE(xt, Uint, maxCalibrationSteps);
  SETVALUE(xt, Double, maxPerformancePredictionError);
  SETVALUE(xt, Double, maxPowerPredictionError);
  SETVALUE(xt, Uint, regressionAging);
  SETVALUE(xt, Double, maxMonitoringOverhead);
  SETVALUE(xt, Bool, clockModulationEmulated);
  SETVALUE(xt, Double, clockModulationMin);
  SETVALUE(xt, Double, clockModulationResolution);
  SETVALUE(xt, Double, clockModulationInterval);
  SETVALUE(xt, Double, thresholdQBlocking);
  SETVALUE(xt, Double, thresholdQBlockingBelt);
  SETVALUE(xt, Uint, tolerableSamples);
  SETVALUE(xt, Ulong, qSize);
  SETVALUE(xt, Double, conservativeValue);
  SETVALUE(xt, ArrayUint, disallowedNumCores);
  SETVALUE(xt, Bool, isolateManager);
  SETVALUE(xt, Bool, statsReconfiguration);
  SETVALUE(xt, Enum, powerDomain);
  SETVALUE(xt, Uint, nelderMeadRange);
  SETVALUE(xt, Bool, fixedPinning);
  SETVALUE(xt, Bool, perPidLog);
  SETVALUE(xt, String, roiFile);
  SETVALUE(xt, ArrayEnums, loggersTypes);
  // xt.getArrayEnums<LoggerType>("loggersTypes", loggersTypes);

  SETVALUE(xt, String, leo.applicationName);
  SETVALUE(xt, String, leo.namesData);
  SETVALUE(xt, String, leo.throughputData);
  SETVALUE(xt, String, leo.powerData);
  SETVALUE(xt, Uint, leo.numSamples);

  SETVALUE(xt, ArrayUint, firstConfiguration.frequency);
  SETVALUE(xt, ArrayUint, firstConfiguration.virtualCores);

  SETVALUE(xt, Bool, dataflow.orderedProcessing);
  SETVALUE(xt, Bool, dataflow.orderedOutput);
  SETVALUE(xt, Uint, dataflow.maxGraphs);
  SETVALUE(xt, Uint, dataflow.maxInterpreters);
}

Parameters::Parameters(Communicator *const communicator) {
  setDefault();
}

Parameters::Parameters(const string &paramFileName,
                       Communicator *const communicator) {
  setDefault();
  /** Loading parameters. **/
  loadXml(paramFileName);
}

Parameters::~Parameters() {
  ;
}

void Parameters::load(const string &paramFileName) {
  setDefault();
  /** Loading parameters. **/
  loadXml(paramFileName);
}

ParametersValidation Parameters::validate() {
  setDefaultPost();

  _knobEnabled[KNOB_FREQUENCY] = knobFrequencyEnabled;
  _knobEnabled[KNOB_VIRTUAL_CORES] = knobCoresEnabled;
  _knobEnabled[KNOB_MAPPING] = knobMappingEnabled;
  _knobEnabled[KNOB_HYPERTHREADING] = knobHyperthreadingEnabled;
  _knobEnabled[KNOB_CLKMOD] = knobClkModEnabled;
  _knobEnabled[KNOB_PFOR_CHUNK] = knobPforChunkEnabled;

  /** Validate frequency knob. **/
  ParametersValidation r = validateKnobFrequencies();
  if (r != VALIDATION_OK) {
    return r;
  }

  r = validateKnobClkMod();
  if (r != VALIDATION_OK) {
    return r;
  }

  /** Validate triggers. **/
  r = validateTriggers();
  if (r != VALIDATION_OK) {
    return r;
  }

  /** Validate unused cores strategy. **/
  r = validateUnusedVc(strategyUnusedVirtualCores);
  if (r != VALIDATION_OK) {
    return r;
  }

  /** Validate requirements. **/
  r = validateRequirements();
  if (r != VALIDATION_OK) {
    return r;
  }

  /** Validate selectors. **/
  r = validateSelector();
  if (r != VALIDATION_OK) {
    return r;
  }

  return VALIDATION_OK;
}

bool Parameters::isKnobEnabled(KnobType k) const {
  return _knobEnabled[k];
}

bool isMinMaxRequirement(double r) {
  return r == NORNIR_REQUIREMENT_MAX || r == NORNIR_REQUIREMENT_MIN;
}

bool isPrimaryRequirement(double r) {
  return r != NORNIR_REQUIREMENT_UNDEF && !isMinMaxRequirement(r);
}

} // namespace nornir

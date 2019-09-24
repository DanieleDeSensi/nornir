/*
 * stats.cpp
 *
 * Created on: 09/07/2015
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
 * \file stats.cpp
 * \brief Statistics collection utilities.
 **/

#include <nornir/stats.hpp>
#include <nornir/configuration.hpp>
#include <nornir/selectors.hpp>
extern "C"{
#include <nornir/graphite-client.h>
}

#include <iomanip>

namespace nornir{

using namespace std;
using namespace mammut;
using namespace mammut::cpufreq;
using namespace mammut::energy;
using namespace mammut::topology;

Logger::Logger(unsigned int timeOffset):
    _timeOffset(timeOffset), _startMonitoring(0){;}

double Logger::getAbsoluteTimestamp(){
    return mammut::utils::getMillisecondsTime();
}

double Logger::getRelativeTimestamp(){
    double absolute = getAbsoluteTimestamp();
    if(!_startMonitoring){
        _startMonitoring = absolute;
    }

    return absolute - _startMonitoring + _timeOffset;
}

void Logger::setStartTimestamp(){
    _startMonitoring = getAbsoluteTimestamp();
}

LoggerStream::LoggerStream(std::ostream *statsStream,
                           std::ostream *calibrationStream,
                           std::ostream *summaryStream,
                           unsigned int timeOffset):
        Logger(timeOffset),
        _statsStream(statsStream), _calibrationStream(calibrationStream),
        _summaryStream(summaryStream), _timeOffset(timeOffset),
        _steadySamples(0), _steadyThroughput(0), _steadyWatts(0){
    if(!*_statsStream ||
       !*_calibrationStream ||
       !*_summaryStream){
        throw runtime_error("LoggerOutStream: Impossible to use stream.");
    }
    *_statsStream << "TimestampMillisecs" << "\t";
    *_statsStream << "[VirtualCores]" << "\t";
    *_statsStream << "Workers" << "\t";
    *_statsStream << "HT" << "\t";
    *_statsStream << "Frequency" << "\t";
    *_statsStream << "ClockModulation" << "\t";
    *_statsStream << "CurrentThroughput" << "\t";
    *_statsStream << "SmoothedThroughput" << "\t";
    *_statsStream << "CoeffVarThroughput" << "\t";
    *_statsStream << "SmoothedLatency" << "\t";
    *_statsStream << "SmoothedUtilization" << "\t";
    *_statsStream << "CurrentWatts" << "\t";
    *_statsStream << "SmoothedWatts" << "\t";
    *_statsStream << endl;

    *_calibrationStream << "NumSteps" << "\t";
    *_calibrationStream << "TimeMs" << "\t";
    *_calibrationStream << "Time%" << "\t";
    *_calibrationStream << "TasksNum" << "\t";
    *_calibrationStream << "Tasks%" << "\t";
    *_calibrationStream << "Watts" << "\t";
    *_calibrationStream << endl;

    *_summaryStream << "Watts" << "\t";
    *_summaryStream << "Throughput" << "\t";
    *_summaryStream << "CompletionTimeSec" << "\t";
    *_summaryStream << "CalibrationSteps" << "\t";
    *_summaryStream << "CalibrationTimeMs" << "\t";
    *_summaryStream << "CalibrationTime%" << "\t";
    *_summaryStream << "CalibrationTasksNum" << "\t";
    *_summaryStream << "CalibrationTasks%" << "\t";
    *_summaryStream << "CalibrationWatts" << "\t";
    for(size_t i = 0; i < KNOB_NUM; i++){
        *_summaryStream << "Reconfigurations" << knobTypeToString((KnobType) i) << "Average" << "\t";
        *_summaryStream << "Reconfigurations" << knobTypeToString((KnobType) i) << "Stddev" << "\t";
    }
    *_summaryStream << "ReconfigurationsTotalAverage" << "\t";
    *_summaryStream << "ReconfigurationsTotalStddev" << "\t";
    *_summaryStream << endl;
}

void LoggerStream::log(bool isCalibrationPhase,
                       const Configuration& configuration,
                       const Smoother<MonitoredSample>& samples,
                       const Requirements& requirements){
    MonitoredSample ms = samples.average();
    *_statsStream << getRelativeTimestamp() << "\t";
    *_statsStream << "[";
    for(size_t c = 0; c < configuration.getNumHMP(); c++){
      const vector<VirtualCore*>& virtualCores = dynamic_cast<const KnobMapping*>(configuration.getKnob(c, KNOB_MAPPING))->getActiveVirtualCores();
      for(size_t i = 0; i < virtualCores.size(); i++){
          *_statsStream << virtualCores.at(i)->getVirtualCoreId() << ",";
      }
      if(configuration.getNumHMP() > 1){
        *_statsStream << "|";
      }
    }
    *_statsStream << "]" << "\t";

    for(size_t c = 0; c < configuration.getNumHMP(); c++){
      *_statsStream << configuration.getRealValue(c, KNOB_VIRTUAL_CORES);
      if(configuration.getNumHMP() > 1){
        *_statsStream << "|";
      }
    }
    *_statsStream << "\t";

    for(size_t c = 0; c < configuration.getNumHMP(); c++){
      *_statsStream << configuration.getRealValue(c, KNOB_HYPERTHREADING);
      if(configuration.getNumHMP() > 1){
        *_statsStream << "|";
      }
    }
    *_statsStream << "\t";
    
    for(size_t c = 0; c < configuration.getNumHMP(); c++){
      Frequency frequency = configuration.getRealValue(c, KNOB_FREQUENCY);
      // Print frequency as string to avoid conversion to exp notation.
      std::ostringstream strs;
      strs << std::fixed << std::setprecision(0) << frequency;
      *_statsStream << strs.str();
      if(configuration.getNumHMP() > 1){
        *_statsStream << "|";
      }
    }
    *_statsStream << "\t";

    for(size_t c = 0; c < configuration.getNumHMP(); c++){
      *_statsStream << configuration.getRealValue(c, KNOB_CLKMOD);
      if(configuration.getNumHMP() > 1){
        *_statsStream << "|";
      }
    }
    *_statsStream << "\t";


    *_statsStream << samples.getLastSample().throughput << "\t";
    *_statsStream << ms.throughput << "\t";
    *_statsStream << samples.coefficientVariation().throughput << "\t";
    *_statsStream << ms.latency << "\t";
    *_statsStream << ms.loadPercentage << "\t";

    *_statsStream << samples.getLastSample().watts << "\t";
    *_statsStream << ms.watts << "\t";

    *_statsStream << endl;

    if(!isCalibrationPhase){
        ++_steadySamples;
        _steadyThroughput += samples.getLastSample().throughput;
        _steadyWatts += samples.getLastSample().watts;
    }
}

void LoggerStream::logSummary(const Configuration& configuration, Selector* selector, ulong durationMs, double totalTasks){
    vector<CalibrationStats> calibrationStats;
    ReconfigurationStats reconfigurationStats = configuration.getReconfigurationStats();
    if(selector){
        selector->updateTotalTasks(totalTasks);
        selector->stopCalibration();
        calibrationStats = selector->getCalibrationsStats();
    }

    for(size_t i = 0; i < calibrationStats.size(); i++){
        const CalibrationStats& cs = calibrationStats.at(i);
        *_calibrationStream << cs.numSteps << "\t";
        *_calibrationStream << cs.duration << "\t";
        *_calibrationStream << ((double) cs.duration / (double)durationMs) * 100.0 << "\t";
        *_calibrationStream << cs.numTasks << "\t";
        *_calibrationStream << ((double) cs.numTasks / totalTasks) * 100.0 << "\t";
        *_calibrationStream << cs.joules / (cs.duration / 1000.0)<< "\t";
        *_calibrationStream << endl;
    }

    CalibrationStats totalCalibration;
    for(size_t i = 0; i < calibrationStats.size(); i++){
        totalCalibration += calibrationStats.at(i);
    }

    if(durationMs == totalCalibration.duration){
        durationMs += 0.001; // Just to avoid division by 0
    }

    *_summaryStream << _steadyWatts / _steadySamples << "\t";
    *_summaryStream << _steadyThroughput / _steadySamples << "\t";
    *_summaryStream << (double) durationMs / 1000.0 << "\t";
    *_summaryStream << totalCalibration.numSteps << "\t";
    *_summaryStream << totalCalibration.duration << "\t";
    *_summaryStream << ((double) totalCalibration.duration / (double)durationMs) * 100.0 << "\t";
    *_summaryStream << totalCalibration.numTasks << "\t";
    *_summaryStream << ((double) totalCalibration.numTasks / totalTasks) * 100.0 << "\t";
    *_summaryStream << totalCalibration.joules / (totalCalibration.duration / 1000.0) << "\t";
    for(size_t i = 0; i < KNOB_NUM; i++){
        if(reconfigurationStats.storedKnob((KnobType) i)){
            *_summaryStream << reconfigurationStats.getAverageKnob((KnobType) i) << "\t";
            *_summaryStream << reconfigurationStats.getStdDevKnob((KnobType) i) << "\t";
        }else{
            *_summaryStream << "N.D." << "\t";
            *_summaryStream << "N.D." << "\t";
        }
    }

    if(reconfigurationStats.storedTotal()){
        *_summaryStream << reconfigurationStats.getAverageTotal() << "\t";
        *_summaryStream << reconfigurationStats.getStdDevTotal() << "\t";
    }else{
        *_summaryStream << "N.D." << "\t";
        *_summaryStream << "N.D." << "\t";
    }

    *_summaryStream << endl;
}

LoggerGraphite::LoggerGraphite(const std::string& host, unsigned int port){
    if(graphite_init(host.c_str(), port)){
        throw std::runtime_error("Impossible to connect to Graphite server.");
    }
}

LoggerGraphite::~LoggerGraphite(){
    graphite_finalize();
}

void LoggerGraphite::log(bool isCalibrationPhase,
                         const Configuration& configuration,
                         const Smoother<MonitoredSample>& samples,
                         const Requirements& requirements){
    unsigned int timestamp = time(NULL);


    /*************************************************/
    /*               Resources info                  */
    /*************************************************/
    // TODO: Fix for HMP (getKnob and getRealValue)
    for(auto vc : dynamic_cast<const KnobMapping*>(configuration.getKnob(0, KNOB_MAPPING))->getActiveVirtualCores()){
        graphite_send_plain((std::string("nornir.resources.cores.") + utils::intToString(vc->getVirtualCoreId())).c_str(), 1, timestamp);
    }
    // Set all other cores to zero.
    for(auto vc : dynamic_cast<const KnobMapping*>(configuration.getKnob(0, KNOB_MAPPING))->getUnusedVirtualCores()){
        graphite_send_plain((std::string("nornir.resources.cores.") + utils::intToString(vc->getVirtualCoreId())).c_str(), 0, timestamp);
    }
    graphite_send_plain("nornir.resources.cores.num", configuration.getRealValue(KNOB_VIRTUAL_CORES), timestamp);
    graphite_send_plain("nornir.resources.frequency", configuration.getRealValue(KNOB_FREQUENCY), timestamp);

    /*************************************************/
    /*                Monitor info                   */
    /*************************************************/
    graphite_send_plain("nornir.monitor.throughput.current", samples.getLastSample().throughput, timestamp);
    graphite_send_plain("nornir.monitor.throughput.average", samples.average().throughput, timestamp);
    graphite_send_plain("nornir.monitor.power.current", samples.getLastSample().watts, timestamp);
    graphite_send_plain("nornir.monitor.power.average", samples.average().watts, timestamp);
    graphite_send_plain("nornir.monitor.latency.current", samples.getLastSample().latency, timestamp);
    graphite_send_plain("nornir.monitor.latency.average", samples.average().latency, timestamp);
    graphite_send_plain("nornir.monitor.utilization.current", samples.getLastSample().loadPercentage, timestamp);
    graphite_send_plain("nornir.monitor.utilization.average", samples.average().loadPercentage, timestamp);

    /*************************************************/
    /*              Requirements info                */
    /*************************************************/
    if(requirements.throughput != NORNIR_REQUIREMENT_MAX &&
       requirements.throughput != NORNIR_REQUIREMENT_UNDEF){
        graphite_send_plain("nornir.requirements.throughput", requirements.throughput, timestamp);
    }

    if(requirements.powerConsumption != NORNIR_REQUIREMENT_MIN &&
       requirements.powerConsumption != NORNIR_REQUIREMENT_UNDEF){
        graphite_send_plain("nornir.requirements.power", requirements.powerConsumption, timestamp);
    }
}

}

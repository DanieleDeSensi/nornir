/*
 * stats.hpp
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
 * \file stats.hpp
 * \brief Statistics collection utilities.
 **/

#ifndef NORNIR_STATS_HPP_
#define NORNIR_STATS_HPP_

#include "knob.hpp"
#include "utils.hpp"

namespace nornir{

class ReconfigurationStats{
private:
    std::vector<double> _knobs[KNOB_NUM];
    std::vector<double> _total;
    bool _storedKnob[KNOB_NUM];
    bool _storedTotal;
public:
    ReconfigurationStats(){
        for(size_t i = 0; i < KNOB_NUM; i++){
            _storedKnob[i] = false;
        }
        _storedTotal = false;
    }

    void swap(ReconfigurationStats& x){
        using std::swap;
        swap(_knobs, x._knobs);
        swap(_total, x._total);
        swap(_storedKnob, x._storedKnob);
        swap(_storedTotal, x._storedTotal);
    }

    inline ReconfigurationStats(const ReconfigurationStats& other){
        for(size_t i = 0; i < KNOB_NUM; i++){
            _knobs[i] = other._knobs[i];
            _storedKnob[i] = other._storedKnob[i];
        }
        _total = other._total;
        _storedTotal = other._storedTotal;
    }

    inline ReconfigurationStats& operator=(ReconfigurationStats other){
        swap(other);
        return *this;
    }

    inline void addSample(KnobType idx, double sample){
        _storedKnob[idx] = true;
        _knobs[idx].push_back(sample);
    }

    inline void addSampleTotal(double total){
        _storedTotal = true;
        _total.push_back(total);
    }

    inline double getAverageKnob(KnobType idx){
        return average(_knobs[idx]);
    }

    inline double getStdDevKnob(KnobType idx){
        return stddev(_knobs[idx]);
    }

    inline double getAverageTotal(){
        return average(_total);
    }

    inline double getStdDevTotal(){
        return stddev(_total);
    }

    inline bool storedKnob(KnobType idx){
        return _storedKnob[idx];
    }

    inline bool storedTotal(){
        return _storedTotal;
    }
};

typedef struct CalibrationStats{
    uint numSteps;
    uint duration;
    uint64_t numTasks;
    mammut::energy::Joules joules;

    CalibrationStats():numSteps(0), duration(0),numTasks(0),joules(0){
        ;
    }

    void swap(CalibrationStats& x){
        using std::swap;

        swap(numSteps, x.numSteps);
        swap(duration, x.duration);
        swap(numTasks, x.numTasks);
        swap(joules, x.joules);
    }

    CalibrationStats& operator=(CalibrationStats rhs){
        swap(rhs);
        return *this;
    }

    CalibrationStats& operator+=(const CalibrationStats& rhs){
        numSteps += rhs.numSteps;
        duration += rhs.duration;
        numTasks += rhs.numTasks;
        joules += rhs.joules;
        return *this;
    }
}CalibrationStats;

class Configuration;
class Selector;

/*!
 * This class can be used to obtain statistics about reconfigurations
 * performed by the manager.
 * It can be extended by a user defined class to customize action to take
 * for each observed statistic.
 */
class Logger{
private:
    unsigned int _timeOffset;
protected:
    double _startMonitoring;

    double getRelativeTimestamp();
    double getAbsoluteTimestamp();
public:
    explicit Logger(unsigned int timeOffset = 0);
    virtual ~Logger(){;}
    void setStartTimestamp();
    virtual void log(bool isCalibrationPhase,
                     const Configuration& configuration,
                     const Smoother<MonitoredSample>& samples,
                     const Requirements& requirements) = 0;
    virtual void logSummary(const Configuration& configuration,
                            Selector* selector, ulong duration, double totalTasks) = 0;
};

/**
 * A logger to store data on C++ streams (e.g. ofstream, etc...).
 */
class LoggerStream: public Logger{
protected:
    std::ostream* _statsStream;
    std::ostream* _calibrationStream;
    std::ostream* _summaryStream;
    unsigned int _timeOffset;
    unsigned long long _steadySamples;
    double _steadyThroughput;
    double _steadyWatts;
public:

    LoggerStream(std::ostream* statsStream,
                 std::ostream* calibrationStream,
                 std::ostream* summaryStream,
                 unsigned int timeOffset = 0);

    void log(bool isCalibrationPhase,
             const Configuration& configuration,
             const Smoother<MonitoredSample>& samples,
             const Requirements& requirements);
    void logSummary(const Configuration& configuration,
                    Selector* selector, ulong durationMs, double totalTasks);
};

/**
 * This logger logs data on files.
 */
class LoggerFile: public LoggerStream{
public:
    LoggerFile(std::string prefix = "",
               std::string folder = ".",
               unsigned int timeOffset = 0):
        LoggerStream(new std::ofstream(folder + std::string("/") + prefix + "stats.csv"),
                     new std::ofstream(folder + std::string("/") + prefix + "calibration.csv"),
                     new std::ofstream(folder + std::string("/") + prefix + "summary.csv"),
                     timeOffset){;}

    /*
    LoggerFile(std::string statsFile = "stats.csv",
               std::string calibrationFile = "calibration.csv",
               std::string summaryFile = "summary.csv",
               unsigned int timeOffset = 0):
        LoggerStream(new std::ofstream(statsFile),
                     new std::ofstream(calibrationFile),
                     new std::ofstream(summaryFile),
                     timeOffset){;}
    */

    ~LoggerFile(){
        dynamic_cast<std::ofstream*>(_statsStream)->close();
        dynamic_cast<std::ofstream*>(_calibrationStream)->close();
        dynamic_cast<std::ofstream*>(_summaryStream)->close();
        delete _statsStream;
        delete _calibrationStream;
        delete _summaryStream;
    }
};

/**
 * This logger can be used to send data to a Graphite (https://graphiteapp.org/)
 * monitoring system (and maybe to show monitored data through a Grafana
 * (https://grafana.com/) front end. To install Graphite and Grafana, you
 * can use an already configured virtual machine, like the one present
 * at (https://github.com/pellepelster/graphite-grafana-vagrant-box).
 *
 * By default, graphite will store data at 60 seconds granularity.
 * If your sampling interval is shorter than 60 seconds, graphite will
 * average the data received by Nornir.
 * All the metrics exported by Nornir are prefixed by "nornir." string.
 *
 * To increase the graphite granularity, please refer to the graphite
 * documentation. However, in most cases it should be sufficient to execute
 * the following steps:
 *  1. Modify the retention by adding an appropriate rule to the
 *     /etc/carbon/storage-schemas.conf (path may be different).
 *     For example, to store one sample for each second and to keep stored
 *     only the samples received during the last day the rule is the following:
 *     (ATTENTION: Put at the top of the file to overwrite other rules)
 *
 *         [default_nornir]
 *         pattern = nornir.*
 *         retentions = 1s:1d
 * 2. Change the scale of the data already stored by Nornir in graphite. You
 *    should execute the following command for all the .wsp files present
 *    in /var/lib/graphite/whisper/nornir folder (path may be different):
 *         $ whisper-resize.py filename.wsp
 *    If you do not need the old data, you can simply remove them by doing:
 *         $ rm -r /var/lib/graphite/whisper/nornir/
 * 3. Change graphite's web app caching, by modifying the value of the
 *    DEFAULT_CACHE_DURATION variable from 60 (default) to 1 (if you want a
 *    1 second granularity) in [graphite/settings.py] file (specific path
 *    depends on your installation).
 * After 60 seconds circa the modification should be loaded automatically by
 * graphite and your data can now be displayed at an higher resolution.
 */
class LoggerGraphite: public Logger{
public:
    LoggerGraphite(const std::string& host, unsigned int port);
    ~LoggerGraphite();

    void log(bool isCalibrationPhase,
             const Configuration& configuration,
             const Smoother<MonitoredSample>& samples,
             const Requirements& requirements);

    // Doesn't log anything.
    void logSummary(const Configuration& configuration, Selector* selector,
                    ulong durationMs, double totalTasks){;}
};

}

#endif

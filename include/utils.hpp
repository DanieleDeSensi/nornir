/*
 * utils.hpp
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
 * \file utils.hpp
 * \brief Implementation of various utilities.
 **/

#ifndef NORNIR_UTILS_HPP_
#define NORNIR_UTILS_HPP_

#include "external/mammut/mammut/mammut.hpp"
#include "external/riff/src/riff.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <sys/stat.h>
#include <sys/types.h>

#define MSECS_IN_SECS 1000.0 // Milliseconds in 1 second
#define NSECS_IN_SECS 1000000000.0 // Nanoseconds in 1 second
#define MAX_RHO 2 //0.0001 (usato per i testcase) //96 //TODO Fix. (occhio ai test)

#define XDG_CONFIG_DIR_FALLBACK "/etc/xdg"

namespace nornir{

inline std::vector<std::string> getXdgConfigDirs(){
    char* confHome_c = getenv("XDG_CONFIG_DIRS");
    std::vector<std::string> confHomes;
    if(!confHome_c || strcmp(confHome_c, "") == 0){
        confHomes.push_back(std::string(XDG_CONFIG_DIR_FALLBACK));
    }else{
        confHomes = mammut::utils::split(std::string(confHome_c), ':');
    }
    for(std::string& s : confHomes){
        s += "/nornir/";
    }
    return confHomes;
}

inline std::string getRuntimeDir(bool userSpecific = true){
    char* runtimeDir_c = getenv("XDG_RUNTIME_DIR");
    if(!runtimeDir_c || strcmp(runtimeDir_c, "") == 0 || !userSpecific){
        runtimeDir_c = (char*) "/tmp/";
    }
    std::string runtimeDir = std::string(runtimeDir_c) + std::string("/nornir/");

    // Create dir if it does not exist
    if(!mammut::utils::existsDirectory(runtimeDir)){
        if(system((std::string("mkdir -p ") + runtimeDir).c_str())){
            throw std::runtime_error("Impossible to create nornir runtime dir.");
        }
        if(!userSpecific){
            if(system((std::string("chmod ugo+rwx ") + runtimeDir).c_str())){
                throw std::runtime_error("Impossible to set permissions on nornir runtime dir.");
            }
        }
    }
    return runtimeDir;
}

inline double ticksToSeconds(double ticks, double ticksPerNs){
    return (ticks/ticksPerNs)/NSECS_IN_SECS;
}

inline double ticksToMilliseconds(double ticks, double ticksPerNs){
    return ticksToSeconds(ticks, ticksPerNs)*1000;
}

inline double average(const std::vector<double>& v){
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    return sum / v.size();
}

inline double stddev(const std::vector<double>& v, double average){
    std::vector<double> diff(v.size());
    std::transform(v.begin(), v.end(), diff.begin(),
                   bind2nd(std::minus<double>(), average));
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    return sqrt(sq_sum / v.size());
}

inline double stddev(const std::vector<double>& v){
    return stddev(v, average(v));
}

typedef struct MonitoredSample: public riff::ApplicationSample{
    mammut::energy::Joules watts; ///< Consumed watts.

    MonitoredSample():riff::ApplicationSample(), watts(0){;}

    MonitoredSample(MonitoredSample const& sample):
        riff::ApplicationSample(sample), watts(sample.watts){;}

    double getMaximumThroughput(){
        if(loadPercentage < MAX_RHO &&
           !inconsistent){
            return throughput / (loadPercentage / 100.0);
        }else{
            return throughput;
        }
    }

    void swap(MonitoredSample& x){
        using std::swap;

        riff::ApplicationSample::swap(x);
        swap(watts, x.watts);
    }

    MonitoredSample& operator=(MonitoredSample rhs){
        swap(rhs);
        return *this;
    }

    MonitoredSample& operator+=(const MonitoredSample& rhs){
        riff::ApplicationSample::operator+=(rhs);
        watts += rhs.watts;
        return *this;
    }

    MonitoredSample& operator-=(const MonitoredSample& rhs){
        riff::ApplicationSample::operator-=(rhs);
        watts -= rhs.watts;
        return *this;
    }

    MonitoredSample& operator*=(const MonitoredSample& rhs){
        riff::ApplicationSample::operator*=(rhs);
        watts *= rhs.watts;
        return *this;
    }

    MonitoredSample& operator/=(const MonitoredSample& rhs){
        riff::ApplicationSample::operator/=(rhs);
        watts /= rhs.watts;
        return *this;
    }

    MonitoredSample operator/=(double x){
        riff::ApplicationSample::operator/=(x);
        watts /= x;
        return *this;
    }

    MonitoredSample operator*=(double x){
        riff::ApplicationSample::operator*=(x);
        watts *= x;
        return *this;
    }
}MonitoredSample;

inline MonitoredSample operator+(const MonitoredSample& lhs,
                                 const MonitoredSample& rhs){
    MonitoredSample r = lhs;
    r += rhs;
    return r;
}

inline MonitoredSample operator-(const MonitoredSample& lhs,
                                 const MonitoredSample& rhs){
    MonitoredSample r = lhs;
    r -= rhs;
    return r;
}

inline MonitoredSample operator*(const MonitoredSample& lhs,
                                 const MonitoredSample& rhs){
    MonitoredSample r = lhs;
    r *= rhs;
    return r;
}

inline MonitoredSample operator*(const MonitoredSample& lhs, double x){
    MonitoredSample r = lhs;
    r *= x;
    return r;
}

inline MonitoredSample operator/(const MonitoredSample& lhs,
                                 const MonitoredSample& rhs){
    MonitoredSample r = lhs;
    r /= rhs;
    return r;
}

inline MonitoredSample operator/(const MonitoredSample& lhs, double x){
    MonitoredSample r = lhs;
    r /= x;
    return r;
}

inline std::ostream& operator<<(std::ostream& os, const MonitoredSample& sample){
    os << "[";
    os << "Watts: " << sample.watts << " ";
    os << "Knarr Sample: " << static_cast<const riff::ApplicationSample&>(sample) << " ";
    os << "]";
    return os;
}

inline std::istream& operator>>(std::istream& is, MonitoredSample& sample){
    is.ignore(std::numeric_limits<std::streamsize>::max(), '[');
    is.ignore(std::numeric_limits<std::streamsize>::max(), ':');
    is >> sample.watts;
    is.ignore(std::numeric_limits<std::streamsize>::max(), ':');
    riff::operator >> (is, sample);
    is.ignore(std::numeric_limits<std::streamsize>::max(), ']');
    return is;
}

inline MonitoredSample squareRoot(const MonitoredSample& x){
    MonitoredSample r;
    if(x.inconsistent){
        r.inconsistent = true;
    }
    r.loadPercentage = sqrt(x.loadPercentage);
    r.throughput = sqrt(x.throughput);
    r.latency = sqrt(x.latency);
    r.numTasks = sqrt(x.numTasks);
    for(size_t i = 0; i < RIFF_MAX_CUSTOM_FIELDS; i++){
        r.customFields[i] = sqrt(x.customFields[i]);
    }
    r.watts = sqrt(x.watts);
    return r;
}

inline void zero(MonitoredSample& x){
    x.loadPercentage = 0;
    x.throughput = 0;
    x.latency = 0;
    x.numTasks = 0;
    for(size_t i = 0; i < RIFF_MAX_CUSTOM_FIELDS; i++){
        x.customFields[i] = 0;
    }
    x.watts = 0;
}

inline void regularize(MonitoredSample& x){
    if(x.loadPercentage < 0){
        x.loadPercentage = 0;
    }
    if(x.throughput < 0){
        x.throughput = 0;
    }
    if(x.latency < 0){
        x.latency = 0;
    }
    if(x.numTasks < 0){
        x.numTasks = 0;
    }
    for(size_t i = 0; i < RIFF_MAX_CUSTOM_FIELDS; i++){
        if(x.customFields[i] < 0){
            x.customFields[i] = 0;
        }
    }
    if(x.watts < 0){
        x.watts = 0;
    }
}

inline MonitoredSample minimum(const MonitoredSample& a,
                               const MonitoredSample& b){
    MonitoredSample ms;
    if(a.inconsistent || b.inconsistent){
        ms.inconsistent = true;
    }

    ms.loadPercentage = std::min(a.loadPercentage, b.loadPercentage);
    ms.throughput = std::min(a.throughput, b.throughput);
    ms.latency = std::min(a.latency, b.latency);
    ms.numTasks = std::min(a.numTasks, b.numTasks);

    for(size_t i = 0; i < RIFF_MAX_CUSTOM_FIELDS; i++){
        ms.customFields[i] = std::min(a.customFields[i], b.customFields[i]);
    }
    ms.watts = std::min(a.watts, b.watts);
    return ms;
}

inline MonitoredSample maximum(const MonitoredSample& a,
                               const MonitoredSample& b){
    MonitoredSample ms;
    if(a.inconsistent || b.inconsistent){
        ms.inconsistent = true;
    }
    ms.loadPercentage = std::max(a.loadPercentage, b.loadPercentage);
    ms.throughput = std::max(a.throughput, b.throughput);
    ms.latency = std::max(a.latency, b.latency);
    ms.numTasks = std::max(a.numTasks, b.numTasks);

    for(size_t i = 0; i < RIFF_MAX_CUSTOM_FIELDS; i++){
        ms.customFields[i] = std::max(a.customFields[i], b.customFields[i]);
    }
    ms.watts = std::max(a.watts, b.watts);
    return ms;
}

inline double squareRoot(const double& x){
    return x?sqrt(x):0;
}

inline void zero(double& x){
    x = 0;
}

inline void regularize(double& x){
    if(x < 0){x = 0;}
}

inline double minimum(const double& a, const double& b){
    return a<b?a:b;
}

inline double maximum(const double& a, const double& b){
    return a>b?a:b;
}

/**
 * Represents a smoothing technique.
 * Requirement: There must exists the following functions:
 *   - 'T squareRoot(const T&)' to compute the square root.
 *   - 'void regularize(T&)' to set the values < 0 to zero.
 *   - 'void zero(T&)' to set to zero.
 *   - 'T minimum(const T& a, const T& b)' returns the minimum between a and b
 *      On struct (or classes) returns a new struct (or class)
 *      with the minimum of each field. For example, if a = {2, 3} and b = {3, 1}
 *      it returns {2, 1}.
 *   - 'T maximum(const T& a, const T& b)' returns the maximum between a and b
 *      On struct (or classes) returns a new struct (or class)
 *      with the maximum of each field. For example, if a = {2, 3} and b = {3, 1}
 *      it returns {3, 3}.
 */
template <typename T> class Smoother{
    template<typename V>
    friend std::ostream& operator<<(std::ostream& os,
                                    const Smoother<V>& obj);

    template<typename V>
    friend std::ofstream& operator<<(std::ofstream& os,
                                     const Smoother<V>& obj);
private:
    T _min, _max;
protected:
    /**
     * Adds a sample to the smoother.
     * @param value The sample to be added.
     */
    virtual void addImpl(const T& value) = 0;
public:
    virtual ~Smoother(){;}

    /**
     * Returns the smoothing factor.
     * @return The smoothing factor.
     */
    virtual double getSmoothingFactor() const = 0;

    /**
     * Sets the smoothing factor.
     * @param s The smoothing factor.
     */
    virtual void setSmoothingFactor(double s) = 0;

    /**
     * Adds a sample to the smoother.
     * @param value The sample to be added.
     */
    void add(const T& value){
        _min = minimum(_min, value);
        _max = maximum(_max, value);
        addImpl(value);
    }

    /**
     * Returns the minimum values added so far.
     * @return The minimum values added so far.
     */
    T min(){return _min;}

    /**
     * Returns the maximum values added so far.
     * @return The maximum values added so far.
     */
    T max(){return _max;}

    /**
     * Gets the last stored sample.
     * @return The last stored sample.
     */
    virtual T getLastSample() const = 0;

    /**
     * Resets the smoother.
     */
    virtual void reset() = 0;

    /**
     * Returns the number of samples stored by the smoother.
     * @return The number of samples stored by the smoother.
     */
    virtual size_t size() const = 0;

    /**
     * Returns the average of the stored samples.
     * @return The average of the stored samples.
     */
    virtual T average() const = 0;

    /**
     * Returns the variance of the stored samples.
     * @return The variance of the stored samples.
     */
    virtual T variance() const = 0;

    /**
     * Returns the standard deviation of the stored samples.
     * @return The standard deviation of the stored samples.
     */
    virtual T standardDeviation() const = 0;

    /**
     * Returns the coefficient of variation of the stored samples [0, 100].
     * @return The coefficient of variation of the stored samples [0, 100].
     */
    virtual T coefficientVariation() const = 0;
};

/**
 * Smoothing technique: Moving average
 */
template<typename T> class MovingAverageSimple: public Smoother<T>{
private:
    std::vector<T> _windowImpl;
    size_t _span;
    size_t _nextIndex;
    size_t _storedValues;
    T _lastSample;
    T _oldAverage, _average;
    T _oldTmpVariance, _tmpVariance;
    T _oldVariance, _variance;
    T _standardDeviation;
    T _coefficientVariation;
public:
    explicit MovingAverageSimple(size_t span):_span(span), _nextIndex(0),
                                     _storedValues(0){
        _windowImpl.resize(_span);
        zero(_lastSample);
        zero(_oldAverage);
        zero(_average);
        zero(_oldTmpVariance);
        zero(_tmpVariance);
        zero(_oldVariance);
        zero(_variance);
        zero(_standardDeviation);
        zero(_coefficientVariation);
    }

    double getSmoothingFactor() const{
        return _span;
    }

    void setSmoothingFactor(double s){
        reset();
        _span = s;
        _windowImpl.resize(_span);
    }

    void addImpl(const T& value){
        _lastSample = value;
        if(_storedValues == 0){
            ++_storedValues;
            _oldAverage = _average = value;
        }else if(_storedValues < _span){
            ++_storedValues;
            _average = _oldAverage + (value - _oldAverage) / _storedValues;
            _tmpVariance = _oldTmpVariance + (value - _oldAverage)*
                                             (value - _average);
            _variance = _tmpVariance / (_storedValues - 1);
        }else{
            T removedValue = _windowImpl.at(_nextIndex);
            _average = _oldAverage + (value - removedValue)/
                                      _span;
            _variance = _oldVariance + (value - _average +
                                        removedValue - _oldAverage)*
                                        (value - removedValue)/(_span - 1);
        }

        _oldAverage = _average;
        _oldTmpVariance = _tmpVariance;
        _oldVariance = _variance;

        regularize(_variance);

        _standardDeviation = squareRoot(_variance);
        _coefficientVariation = (_standardDeviation / _average) * 100.0;

        _windowImpl.at(_nextIndex) = value;
        _nextIndex = (_nextIndex + 1) % _span;
    }

    T getLastSample() const{
        return _lastSample;
    }

    void reset(){
        _windowImpl.clear();
        _windowImpl.resize(_span);
        _nextIndex = 0;
        _storedValues = 0;
        _lastSample = T();
        _oldAverage = T();
        _average = T();
        _oldTmpVariance = T();
        _tmpVariance = T();
        _oldVariance = T();
        _variance = T();
        _standardDeviation = T();
        _coefficientVariation = T();
    }

    size_t size() const{
        return _storedValues;
    }

    T average() const{
        return _average;
    }

    T variance() const{
        return _variance;
    }

    T standardDeviation() const{
        return _standardDeviation;
    }

    T coefficientVariation() const{
        return _coefficientVariation;
    }
};


template<typename T> class MovingAverageExponential: public Smoother<T>{
private:
    double _alpha;
    size_t _storedValues;
    T _lastSample;
    T _average;
    T _variance;
    T _standardDeviation;
    T _coefficientVariation;
public:
    explicit MovingAverageExponential(double alpha):_alpha(alpha),
                                                     _storedValues(0){
        if(_alpha < 0 || _alpha > 1.0){
            throw std::runtime_error("Alpha must be between 0 and 1 "
                                     "(included)");
        }
        zero(_lastSample);
        zero(_average);
        zero(_variance);
        zero(_standardDeviation);
        zero(_coefficientVariation);
    }

    double getSmoothingFactor() const{
        return _alpha;
    }

    void setSmoothingFactor(double s){
        _alpha = s;
    }

    void addImpl(const T& value){
        ++_storedValues;
        _lastSample = value;
        if(_storedValues == 1){
            _average = value;
            _variance = T();
        }else{
            T diff = value - _average;
            T incr = diff * _alpha;
            _average += incr;
            _variance = (_variance + diff * incr) * (1 - _alpha);
        }
        regularize(_variance);
        _standardDeviation = squareRoot(_variance);
        _coefficientVariation = (_standardDeviation / _average) * 100.0;
    }

    T getLastSample() const{
        return _lastSample;
    }

    void reset(){
        _storedValues = 0;
        _lastSample = T();
        _average = T();
        _variance = T();
        _standardDeviation = T();
        _coefficientVariation = T();
    }

    size_t size() const{
        return _storedValues;
    }

    T average() const{
        return _average;
    }

    T variance() const{
        return _variance;
    }

    T standardDeviation() const{
        return _standardDeviation;
    }

    T coefficientVariation() const{
        return _coefficientVariation;
    }
};

template<class T>
std::ostream& operator<<(std::ostream& os, const Smoother<T>& obj){
    os << "==============================" << std::endl;
    os << "Last sample: " << obj.getLastSample() << std::endl;
    os << "Average: " << obj.average() << std::endl;
    os << "CoeffVar: " << obj.coefficientVariation() << std::endl;
    os << "Variance: " << obj.variance() << std::endl;
    os << "StdDev: " << obj.standardDeviation() << std::endl;
    os << "==============================" << std::endl;
    return os;
}

template<class T>
std::ofstream& operator<<(std::ofstream& os, const Smoother<T>& obj){
    os << obj.getLastSample() << "\t";
    os << obj.average() << "\t";
    os << obj.coefficientVariation() << "\t";
    os << obj.variance() << "\t";
    os << obj.standardDeviation() << "\t";
    return os;
}

template<typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v){
    if(!v.empty()){
        out << '[';
        std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
        out << "\b\b]";
    }
    return out;
}

/*****************************************************************************************/
/* To disable compiler warnings.                                                         */
/* Code obtained from: https://github.com/facebook/folly/blob/master/folly/Portability.h */
/*****************************************************************************************/

// Generalize warning push/pop.
#if defined(_MSC_VER)
# define PUSH_WARNING __pragma(warning(push))
# define POP_WARNING __pragma(warning(pop))
// Disable the GCC warnings.
# define GCC_DISABLE_WARNING(warningName)
# define MSVC_DISABLE_WARNING(warningNumber) __pragma(warning(disable: warningNumber))
#elif defined(__clang__) || __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
# define PUSH_WARNING _Pragma("GCC diagnostic push")
# define POP_WARNING _Pragma("GCC diagnostic pop")
#define GCC_DISABLE_WARNING_INTERNAL3(warningName) #warningName
#define GCC_DISABLE_WARNING_INTERNAL2(warningName) \
  GCC_DISABLE_WARNING_INTERNAL3(warningName)
#define GCC_DISABLE_WARNING(warningName)                       \
  _Pragma(GCC_DISABLE_WARNING_INTERNAL2(GCC diagnostic ignored \
          GCC_DISABLE_WARNING_INTERNAL3(-W##warningName)))
// Disable the MSVC warnings.
# define MSVC_DISABLE_WARNING(warningNumber)
#else
# define PUSH_WARNING
# define POP_WARNING
# define GCC_DISABLE_WARNING(warningName)
# define MSVC_DISABLE_WARNING(warningNumber)
#endif

}

#endif /* NORNIR_UTILS_HPP_ */

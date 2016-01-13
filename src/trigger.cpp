/*
 * trigger.cpp
 *
 * Created on: 02/01/2016
 *
 * =========================================================================
 *  Copyright (C) 2015-, Daniele De Sensi (d.desensi.software@gmail.com)
 *
 *  This file is part of AdaptiveFastFlow.
 *
 *  AdaptiveFastFlow is free software: you can redistribute it and/or
 *  modify it under the terms of the Lesser GNU General Public
 *  License as published by the Free Software Foundation, either
 *  version 3 of the License, or (at your option) any later version.

 *  AdaptiveFastFlow is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  Lesser GNU General Public License for more details.
 *
 *  You should have received a copy of the Lesser GNU General Public
 *  License along with AdaptiveFastFlow.
 *  If not, see <http://www.gnu.org/licenses/>.
 *
 * =========================================================================
 */

#include "trigger.hpp"

namespace adpff{

TriggerQBlocking::TriggerQBlocking(TriggerConfQBlocking confQBlocking,
                                   double thresholdQBlocking,
                                   Smoother<MonitoredSample> const* samples,
                                   AdaptiveNode* emitter):
        _confQBlocking(confQBlocking),
        _thresholdQBlocking(thresholdQBlocking),
        _samples(samples),
        _emitter(emitter),
        _blocking(false){
    ;
}

double TriggerQBlocking::getIdleTime() const{
    MonitoredSample avg = _samples->average();

    /**
     * Latency = Utilization * Interarrival
     * Idle = Interarrival - Latency
     *
     * Interarrival = Latency/Utilization
     * Idle = Latency/Utilization - Latency
     **/

    // We need to convert to microsec since the threshold is specified in
    // microsec.
    double latencyMicroSec = avg.latency / 1000.0;

    // We convert from [0, 100] to [0, 1]
    double utilization = avg.utilization / 100.0;
    if(!utilization){utilization = 1;} // Should never happen

    return latencyMicroSec/utilization - latencyMicroSec;
}

bool TriggerQBlocking::setBlocking(){
    if(!_blocking){
        _emitter->setQBlocking();
        _blocking = true;
        return true;
    }
    return false;
}

bool TriggerQBlocking::setNonBlocking(){
    if(_blocking){
        _emitter->setQNonblocking();
        _blocking = false;
        return true;
    }
    return false;
}

bool TriggerQBlocking::trigger(){
    switch(_confQBlocking){
        case TRIGGER_Q_BLOCKING_YES:{
            return setBlocking();
        }break;
        case TRIGGER_Q_BLOCKING_NO:{
            return setNonBlocking();
        }break;
        case TRIGGER_Q_BLOCKING_AUTO:{
            double idleTime = getIdleTime();
            if(idleTime > _thresholdQBlocking){
                return setBlocking();
            }else if(idleTime < _thresholdQBlocking){
                return setNonBlocking();
            }
        }break;
    }
    return false;
}
}

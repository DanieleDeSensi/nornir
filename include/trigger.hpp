/*
 * trigger.hpp
 *
 * Created on: 02/01/2016
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

#ifndef NORNIR_TRIGGER_HPP_
#define NORNIR_TRIGGER_HPP_

#include "node.hpp"
#include "utils.hpp"

namespace nornir{

typedef enum{
    TRIGGER_TYPE_Q_BLOCKING = 0,
    TRIGGER_TYPE_NUM // <---- This must always be the last value
}TriggerType;

class Trigger: public mammut::utils::NonCopyable{
public:
    /**
     * Triggers the trigger, if needed.
     * @return True if the trigger has been triggered, false otherwise.
     */
    virtual bool trigger() = 0;

    virtual ~Trigger(){;}
};

/**
 * When the idle time is == 0, blocking reaches the same performances of non
 * blocking only for high grain computation. However, since it has no advantages
 * with respect to nonblocking, in this case we always use a nonblocking support.
 *
 * If the idle time is > 0, blocking support is more convenient than nonblocking
 * from a power consumption point of view for very long idle times (specified
 * through the 'thresholdQBlocking', independently from the computation grain.
 */
class TriggerQBlocking: public Trigger{
public:
    TriggerQBlocking(TriggerConfQBlocking confQBlocking,
                     double thresholdQBlocking,
                     double thresholdQBlockingBelt,
                     Smoother<MonitoredSample> const* samples,
                     AdaptiveNode* emitter);
    bool trigger();
private:
    double getIdleTime() const;
    bool setBlocking();
    bool setNonBlocking();

    TriggerConfQBlocking _confQBlocking;
    const double _thresholdQBlocking;
    const double _thresholdQBlockingBelt;
    Smoother<MonitoredSample> const* _samples;
    AdaptiveNode* _emitter;
    bool _blocking;
};

}

#endif

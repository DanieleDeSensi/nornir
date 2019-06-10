/**
 *  Different tests on triggers.
 **/
#include "parametersLoader.hpp"
#include <algorithm>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include "../src/nornir.hpp"
#include "gtest/gtest.h"

using namespace nornir;
using namespace ff;
using namespace mammut;
using namespace mammut::cpufreq;
using namespace mammut::energy;
using namespace mammut::task;
using namespace mammut::topology;
using namespace mammut::utils;

int dummyMsg;

class DummyEmitter: public AdaptiveNode{
    void* svc(void*){
        return (void*) &dummyMsg;
    }
};

TEST(TriggersTest, Blocking){
    getParameters("repara"); // Only to force setting XDG_CONFIG_DIRS environment variables
    MovingAverageSimple<MonitoredSample> samples(4);
    DummyEmitter emitter;
    TriggerQBlocking trigger(TRIGGER_Q_BLOCKING_YES, 0, 0, &samples, &emitter);
    EXPECT_TRUE(trigger.trigger());
    ManagementRequest* mr = emitter.checkManagementRequest();
    EXPECT_EQ(mr->type, MGMT_REQ_SWITCH_BLOCKING);
    EXPECT_EQ(mr->mark, (void*) FF_BLK);
}

TEST(TriggersTest, NonBlocking){
    getParameters("repara"); // Only to force setting XDG_CONFIG_DIRS environment variables
    MovingAverageSimple<MonitoredSample> samples(4);
    DummyEmitter emitter;
    TriggerQBlocking trigger(TRIGGER_Q_BLOCKING_NO, 0, 0, &samples, &emitter);
    // We expect false because it should be non blocking by default so nothing should be changed
    EXPECT_FALSE(trigger.trigger());
}


TEST(TriggersTest, Auto){
    getParameters("repara"); // Only to force setting XDG_CONFIG_DIRS environment variables
    MovingAverageSimple<MonitoredSample> samples(4);
    DummyEmitter emitter;
    // < 90  -> NonBlocking
    // > 110 -> Blocking
    TriggerQBlocking trigger(TRIGGER_Q_BLOCKING_AUTO, 100, 0.1, &samples, &emitter);
    // First try to trigger blocking
    MonitoredSample ms;
    ms.latency = 100*1000.0;
    ms.loadPercentage = 10;
    // Latency 100 and loadPercentage 10 -> Idle = 900
    samples.add(ms);
    EXPECT_TRUE(trigger.trigger());
    ManagementRequest* mr = emitter.checkManagementRequest();
    EXPECT_EQ(mr->type, MGMT_REQ_SWITCH_BLOCKING);
    EXPECT_EQ(mr->mark, (void*) FF_BLK);
    // Then try to trigger nonblocking
    samples.reset();
    ms.latency = 100*1000.0;
    ms.loadPercentage = 90;
    // Latency 100 and loadPercentage 90 -> Idle = 11.11
    samples.add(ms);
    EXPECT_TRUE(trigger.trigger());
    mr = emitter.checkManagementRequest();
    EXPECT_EQ(mr->type, MGMT_REQ_SWITCH_BLOCKING);
    EXPECT_EQ(mr->mark, (void*) FF_NBLK);
}

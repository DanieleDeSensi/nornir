/**
 *  Different tests on configuration.
 **/
#include "parametersLoader.hpp"
#include <algorithm>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include "../src/nornir.hpp"
#include "gtest/gtest.h"

using namespace nornir;
using namespace mammut;
using namespace mammut::cpufreq;
using namespace mammut::energy;
using namespace mammut::task;
using namespace mammut::topology;
using namespace mammut::utils;

TEST(ConfigurationTest, SimpleTest) {
    Parameters  p = getParameters("repara");
    p.strategyUnusedVirtualCores = STRATEGY_UNUSED_VC_OFF;
    p.knobHyperthreadingEnabled = true;
    ConfigurationExternal configuration(p);
    dynamic_cast<KnobMappingExternal*>(configuration.getKnob(KNOB_MAPPING))->setPid(getpid());
    dynamic_cast<KnobClkModEmulated*>(configuration.getKnob(KNOB_CLKMOD))->setPid(getpid());
    configuration.createAllRealCombinations();
    EXPECT_EQ(configuration.getNumServiceNodes(), (uint) 0);
    EXPECT_TRUE(configuration.knobsChangeNeeded());
    configuration.maxAllKnobs();
    KnobsValues kv = configuration.getRealValues();
    EXPECT_EQ(kv[KNOB_VIRTUAL_CORES], 48);
    EXPECT_EQ(kv[KNOB_HYPERTHREADING], 2);
    EXPECT_EQ(kv[KNOB_MAPPING], MAPPING_TYPE_INTERLEAVED);
    EXPECT_EQ(kv[KNOB_FREQUENCY], 2400000);
    EXPECT_EQ(kv[KNOB_CLKMOD], 100.0);

    // Test equality and correct frequency set.
    ConfigurationExternal configuration2(p);
    dynamic_cast<KnobMappingExternal*>(configuration2.getKnob(KNOB_MAPPING))->setPid(getpid());
    dynamic_cast<KnobClkModEmulated*>(configuration2.getKnob(KNOB_CLKMOD))->setPid(getpid());
    KnobsValues kv2(KNOB_VALUE_REAL);
    kv2[KNOB_VIRTUAL_CORES] = 12;
    kv2[KNOB_HYPERTHREADING] = 1;
    kv2[KNOB_MAPPING] = MAPPING_TYPE_LINEAR;
    kv2[KNOB_FREQUENCY] = 1600000;
    kv2[KNOB_CLKMOD] = 100.0;
    EXPECT_FALSE(configuration.equal(kv2));
    configuration2.setValues(kv2);
    EXPECT_EQ(p.mammut.getInstanceCpuFreq()->getDomains()[0]->getCurrentGovernor(), GOVERNOR_USERSPACE);
    EXPECT_EQ(p.mammut.getInstanceCpuFreq()->getDomains()[0]->getCurrentFrequencyUserspace(), (Frequency) 1600000);
    // Second domain should be off
    EXPECT_EQ(p.mammut.getInstanceCpuFreq()->getDomains()[1]->getCurrentGovernor(), GOVERNOR_USERSPACE);
    for(mammut::topology::VirtualCore* vc : p.mammut.getInstanceCpuFreq()->getDomains()[1]->getVirtualCores()){
        EXPECT_FALSE(vc->isHotPlugged());
    }

    // Test unneded change configuration.
    ConfigurationExternal configuration3(p);
    dynamic_cast<KnobMappingExternal*>(configuration3.getKnob(KNOB_MAPPING))->setPid(getpid());
    dynamic_cast<KnobClkModEmulated*>(configuration3.getKnob(KNOB_CLKMOD))->setPid(getpid());
    for(size_t i = 0; i < KNOB_NUM; i++){
        configuration3.getKnob((KnobType) i)->lockToMax();
    }
    EXPECT_FALSE(configuration3.knobsChangeNeeded());
}

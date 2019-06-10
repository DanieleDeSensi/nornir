/**
 *  Different tests on knobs.
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
using namespace mammut;
using namespace mammut::cpufreq;
using namespace mammut::energy;
using namespace mammut::task;
using namespace mammut::topology;
using namespace mammut::utils;

TEST(KnobsTest, WriteRead){
    KnobsValues kv(KNOB_VALUE_REAL), validation(KNOB_VALUE_REAL);
    kv[KNOB_VIRTUAL_CORES] = 1;
    kv[KNOB_HYPERTHREADING] = 10;
    kv[KNOB_MAPPING] = MAPPING_TYPE_INTERLEAVED;
    kv[KNOB_FREQUENCY] = 2.4;
    kv[KNOB_CLKMOD] = 20.0;
    std::ofstream outstream;
    outstream.open("/tmp/testknobs.dat");
    outstream << kv;
    outstream.close();
    std::ifstream instream("/tmp/testknobs.dat");
    instream >> validation;
    instream.close();
    for(size_t i = 0; i < (size_t) KNOB_NUM; i++){
        EXPECT_EQ(kv[(KnobType) i], validation[(KnobType) i]);
    }
}

TEST(KnobsTest, KnobsVirtualCores) {
    Parameters  p = getParameters("repara");
    p.knobHyperthreadingEnabled = false;
    KnobVirtualCores knob(p);

    // Check values.
    std::vector<double> values = knob.getAllowedValues();
    EXPECT_EQ(values.size(), (size_t) 24);
    for(size_t i = 0; i < values.size(); i++){
        EXPECT_EQ(values[i], i + 1);
    }

    // Check max/min
    knob.setRelativeValue(100);
    EXPECT_EQ(knob.getRealValue(), 24);
    knob.setRelativeValue(0);
    EXPECT_EQ(knob.getRealValue(), 1);

    // Knob-specific calls
    // Check max/min after changing the max
    knob.changeMax(12);
    knob.setRelativeValue(0);
    EXPECT_EQ(knob.getRealValue(), 1);
    knob.setRelativeValue(100);
    EXPECT_EQ(knob.getRealValue(), 12);
    knob.changeMax(24);

    // Lock
    knob.lockToMax();
    EXPECT_TRUE(knob.isLocked());
    EXPECT_EQ(knob.getRealValue(), 24);
    // After lock we lost the values. So to
    // try lockToMin we need to build a new knob.
    KnobVirtualCores knob2(p);
    knob2.lockToMin();
    EXPECT_TRUE(knob2.isLocked());
    EXPECT_EQ(knob2.getRealValue(), 1);
}

// Test trying to disallow some number of cores (used in Nornir dataflow runtime).
TEST(KnobsTest, KnobsVirtualCores2) {
    Parameters  p = getParameters("repara");
    for(size_t i = 0; i < 48; i++){
        if(i % 3 == 0){
            p.disallowedNumCores.push_back(i);
        }
    }
    KnobVirtualCores knob(p);
    // Check values.
    std::vector<double> values = knob.getAllowedValues();
    // With % 3 we skipped 8 values.
    EXPECT_EQ(values.size(), (size_t) 24 - 8);
    size_t expected = 1;
    for(size_t i = 0; i < values.size(); i++){
        EXPECT_EQ(values[i], expected);
        if((expected + 1) % 3 == 0){
            expected += 2;
        }else{
            expected += 1;
        }
    }
}

// Test trying to enable/disable hyperthreading
TEST(KnobsTest, KnobsVirtualCores3) {
    Parameters  p = getParameters("repara");
    p.knobHyperthreadingEnabled = true;
    KnobVirtualCores knob(p);
    // Check values.
    std::vector<double> values = knob.getAllowedValues();
    EXPECT_EQ(values.size(), (size_t) 48);
    for(size_t i = 0; i < values.size(); i++){
        EXPECT_EQ(values[i], i + 1);
    }

    p.knobHyperthreadingEnabled = false;
    KnobVirtualCores knob2(p);
    // Check values.
    values = knob2.getAllowedValues();
    EXPECT_EQ(values.size(), (size_t) 24);
    for(size_t i = 0; i < values.size(); i++){
        EXPECT_EQ(values[i], i + 1);
    }
}

TEST(KnobsTest, KnobsHyperThreading){
    Parameters  p = getParameters("repara");
    KnobHyperThreading knob(p);

    // Check values.
    std::vector<double> values = knob.getAllowedValues();
    EXPECT_EQ(values.size(), (size_t) 2);
    for(size_t i = 0; i < values.size(); i++){
        EXPECT_EQ(values[i], i + 1);
    }

    // Check max/min
    knob.setRelativeValue(100);
    EXPECT_EQ(knob.getRealValue(), 2);
    knob.setRelativeValue(0);
    EXPECT_EQ(knob.getRealValue(), 1);

    // Lock
    knob.lockToMax();
    EXPECT_TRUE(knob.isLocked());
    EXPECT_EQ(knob.getRealValue(), 2);
    // After lock we lost the values. So to
    // try lockToMin we need to build a new knob.
    KnobHyperThreading knob2(p);
    knob2.lockToMin();
    EXPECT_TRUE(knob2.isLocked());
    EXPECT_EQ(knob2.getRealValue(), 1);
}

TEST(KnobsTest, KnobsMapping){
    Parameters  p = getParameters("repara");
    p.knobHyperthreadingEnabled = true;
    KnobVirtualCores knobCores(p);
    knobCores.setToMax();
    KnobHyperThreading knobHT(p);
    knobHT.setToMax();
    KnobMappingExternal knob(p, knobCores, knobHT);
    knob.setPid(getpid());

    // Check values.
    std::vector<double> values = knob.getAllowedValues();
    EXPECT_EQ(values.size(), MAPPING_TYPE_NUM);
    for(size_t i = 0; i < values.size(); i++){
        EXPECT_EQ(values[i], i);
    }

    // Knob-specific calls
    // Testing linear mapping
    knob.setRealValue(MAPPING_TYPE_LINEAR);
    std::vector<mammut::topology::VirtualCore*> vcs = knob.getActiveVirtualCores();
    EXPECT_EQ(vcs.size(), (size_t) 48);
    for(size_t i = 0; i < vcs.size(); i++){
        EXPECT_EQ(vcs[i]->getVirtualCoreId(), i);
    }

    // Testing interleaved mapping
    knob.setRealValue(MAPPING_TYPE_INTERLEAVED);
    vcs = knob.getActiveVirtualCores();
    EXPECT_EQ(vcs.size(), (size_t) 48);
    for(size_t i = 0; i < vcs.size(); i++){
        int offset = 0;
        if((i % 2) == 0){
            if(i < 24){
                offset = - (i / 2);
            }else{
                offset = - ((i - 24) / 2);
            }
            EXPECT_EQ(vcs[i]->getVirtualCoreId(), i + offset);
        }else{
            if(i < 24){
                offset = 12 - (i+1)/2;
            }else{
                offset = 12 - ((i-24)+1)/2;
            }
            EXPECT_EQ(vcs[i]->getVirtualCoreId(), i + offset);
        }
    }

    // Testing linear mapping and only 1 level of hyperthreading
    knobHT.setRealValue(1);
    // Testing linear mapping
    knob.setRealValue(MAPPING_TYPE_LINEAR);
    vcs = knob.getActiveVirtualCores();
    EXPECT_EQ(vcs.size(), (size_t) 24);
    for(size_t i = 0; i < vcs.size(); i++){
        EXPECT_EQ(vcs[i]->getVirtualCoreId(), i);
    }
    // Testing interleaved mapping
    knob.setRealValue(MAPPING_TYPE_INTERLEAVED);
    vcs = knob.getActiveVirtualCores();
    EXPECT_EQ(vcs.size(), (size_t) 24);
    for(size_t i = 0; i < vcs.size(); i++){
        int offset = 0;
        if((i % 2) == 0){
            offset = - (i / 2);
            EXPECT_EQ(vcs[i]->getVirtualCoreId(), i + offset);
        }else{
            offset = 12 - (i+1)/2;
            EXPECT_EQ(vcs[i]->getVirtualCoreId(), i + offset);
        }
    }

    // Isolate manager
    p.isolateManager = true;
    knobHT.setRealValue(2);
    KnobMappingExternal knob2(p, knobCores, knobHT);
    knob2.setPid(getpid());
    knob2.setRealValue(MAPPING_TYPE_LINEAR);
    vcs = knob2.getActiveVirtualCores();
    vcs = knob2.getActiveVirtualCores();
    EXPECT_FALSE(mammut::utils::contains(vcs, p.mammut.getInstanceTopology()->getVirtualCore(NORNIR_MANAGER_VIRTUAL_CORE)));
}

TEST(KnobsTest, KnobsFrequency) {
    Parameters  p = getParameters("repara");
    p.knobHyperthreadingEnabled = false;
    KnobVirtualCores knobCores(p);
    KnobHyperThreading knobHT(p);
    KnobMappingExternal knobMapping(p, knobCores, knobHT);
    knobMapping.setPid(getpid());
    KnobFrequency knob(p, knobMapping);

    // Check values.
    std::vector<double> values = knob.getAllowedValues();
    EXPECT_EQ(values.size(), (size_t) 13);
    for(size_t i = 0; i < values.size(); i++){
        EXPECT_EQ(values[i], 1200000 + (i*100000));
    }

    // Check max/min
    knob.setRelativeValue(100);
    EXPECT_EQ(knob.getRealValue(), (Frequency) 2400000);
    knob.setRelativeValue(0);
    EXPECT_EQ(knob.getRealValue(), (Frequency) 1200000);

    // Lock
    knob.lockToMax();
    EXPECT_TRUE(knob.isLocked());
    EXPECT_EQ(knob.getRealValue(), (Frequency) 2400000);
    // After lock we lost the values. So to
    // try lockToMin we need to build a new knob.
    KnobFrequency knob2(p, knobMapping);
    knob2.lockToMin();
    EXPECT_TRUE(knob2.isLocked());
    EXPECT_EQ(knob2.getRealValue(), (Frequency) 1200000);
}

// Global test with strategy for unused virtual cores = NONE
TEST(KnobsTest, GlobalUnusedNone){
    Parameters  p = getParameters("repara");
    p.knobHyperthreadingEnabled = true;
    p.strategyUnusedVirtualCores = STRATEGY_UNUSED_VC_NONE;
    KnobVirtualCores knobCores(p);
    KnobHyperThreading knobHT(p);
    KnobMappingExternal knobMapping(p, knobCores, knobHT);
    knobMapping.setPid(getpid());
    KnobFrequency knobFrequency(p, knobMapping);

    // Both domains - linear
    knobCores.setRealValue(14);
    knobHT.setRealValue(2);
    knobMapping.setRealValue(MAPPING_TYPE_LINEAR);
    knobFrequency.setRealValue(1700000);
    for(mammut::cpufreq::Domain* d : p.mammut.getInstanceCpuFreq()->getDomains()){
        EXPECT_EQ(d->getCurrentGovernor(), GOVERNOR_USERSPACE);
        EXPECT_EQ(d->getCurrentFrequencyUserspace(), (Frequency) 1700000);
    }

    // One domain - linear
    knobCores.setRealValue(4);
    knobHT.setRealValue(2);
    knobMapping.setRealValue(MAPPING_TYPE_LINEAR);
    knobFrequency.setRealValue(1800000);
    EXPECT_EQ(p.mammut.getInstanceCpuFreq()->getDomains()[0]->getCurrentGovernor(), GOVERNOR_USERSPACE);
    EXPECT_EQ(p.mammut.getInstanceCpuFreq()->getDomains()[0]->getCurrentFrequencyUserspace(), (Frequency) 1800000);
    // Second domain should be left unchanged
    EXPECT_EQ(p.mammut.getInstanceCpuFreq()->getDomains()[1]->getCurrentGovernor(), GOVERNOR_USERSPACE);
    EXPECT_EQ(p.mammut.getInstanceCpuFreq()->getDomains()[1]->getCurrentFrequencyUserspace(), (Frequency) 1700000);


    // Both domains - interleaved
    knobCores.setRealValue(14);
    knobHT.setRealValue(2);
    knobMapping.setRealValue(MAPPING_TYPE_INTERLEAVED);
    knobFrequency.setRealValue(1700000);
    for(mammut::cpufreq::Domain* d : p.mammut.getInstanceCpuFreq()->getDomains()){
        EXPECT_EQ(d->getCurrentGovernor(), GOVERNOR_USERSPACE);
        EXPECT_EQ(d->getCurrentFrequencyUserspace(), (Frequency) 1700000);
    }

    // One domain - interleaved
    knobCores.setRealValue(4);
    knobHT.setRealValue(2);
    knobMapping.setRealValue(MAPPING_TYPE_INTERLEAVED);
    knobFrequency.setRealValue(1800000);
    for(mammut::cpufreq::Domain* d : p.mammut.getInstanceCpuFreq()->getDomains()){
        EXPECT_EQ(d->getCurrentGovernor(), GOVERNOR_USERSPACE);
        EXPECT_EQ(d->getCurrentFrequencyUserspace(), (Frequency) 1800000);
    }
}

// Global test with strategy for unused virtual cores = LOWEST_FREQUENCY
TEST(KnobsTest, GlobalUnusedLowest){
    Parameters  p = getParameters("repara");
    p.knobHyperthreadingEnabled = true;
    p.strategyUnusedVirtualCores = STRATEGY_UNUSED_VC_LOWEST_FREQUENCY;
    KnobVirtualCores knobCores(p);
    KnobHyperThreading knobHT(p);
    KnobMappingExternal knobMapping(p, knobCores, knobHT);
    knobMapping.setPid(getpid());
    KnobFrequency knobFrequency(p, knobMapping);

    // Both domains - linear
    knobCores.setRealValue(14);
    knobHT.setRealValue(2);
    knobMapping.setRealValue(MAPPING_TYPE_LINEAR);
    knobFrequency.setRealValue(1700000);
    for(mammut::cpufreq::Domain* d : p.mammut.getInstanceCpuFreq()->getDomains()){
        EXPECT_EQ(d->getCurrentGovernor(), GOVERNOR_USERSPACE);
        EXPECT_EQ(d->getCurrentFrequencyUserspace(), (Frequency) 1700000);
    }

    // One domain - linear
    knobCores.setRealValue(4);
    knobHT.setRealValue(2);
    knobMapping.setRealValue(MAPPING_TYPE_LINEAR);
    knobFrequency.setRealValue(1800000);
    EXPECT_EQ(p.mammut.getInstanceCpuFreq()->getDomains()[0]->getCurrentGovernor(), GOVERNOR_USERSPACE);
    EXPECT_EQ(p.mammut.getInstanceCpuFreq()->getDomains()[0]->getCurrentFrequencyUserspace(), (Frequency) 1800000);
    // Second domain should be at lowest frequency
    EXPECT_EQ(p.mammut.getInstanceCpuFreq()->getDomains()[1]->getCurrentGovernor(), GOVERNOR_USERSPACE);
    EXPECT_EQ(p.mammut.getInstanceCpuFreq()->getDomains()[1]->getCurrentFrequencyUserspace(), (Frequency) 1200000);


    // Both domains - interleaved
    knobCores.setRealValue(14);
    knobHT.setRealValue(2);
    knobMapping.setRealValue(MAPPING_TYPE_INTERLEAVED);
    knobFrequency.setRealValue(1700000);
    for(mammut::cpufreq::Domain* d : p.mammut.getInstanceCpuFreq()->getDomains()){
        EXPECT_EQ(d->getCurrentGovernor(), GOVERNOR_USERSPACE);
        EXPECT_EQ(d->getCurrentFrequencyUserspace(), (Frequency) 1700000);
    }

    // One domain - interleaved
    knobCores.setRealValue(4);
    knobHT.setRealValue(2);
    knobMapping.setRealValue(MAPPING_TYPE_INTERLEAVED);
    knobFrequency.setRealValue(1800000);
    for(mammut::cpufreq::Domain* d : p.mammut.getInstanceCpuFreq()->getDomains()){
        EXPECT_EQ(d->getCurrentGovernor(), GOVERNOR_USERSPACE);
        EXPECT_EQ(d->getCurrentFrequencyUserspace(), (Frequency) 1800000);
    }
}

// Global test with strategy for unused virtual cores = OFF
TEST(KnobsTest, GlobalUnusedOff){
    Parameters  p = getParameters("repara");
    p.knobHyperthreadingEnabled = true;
    p.strategyUnusedVirtualCores = STRATEGY_UNUSED_VC_OFF;
    KnobVirtualCores knobCores(p);
    KnobHyperThreading knobHT(p);
    KnobMappingExternal knobMapping(p, knobCores, knobHT);
    knobMapping.setPid(getpid());
    KnobFrequency knobFrequency(p, knobMapping);

    // Both domains - linear
    knobCores.setRealValue(14);
    knobHT.setRealValue(2);
    knobMapping.setRealValue(MAPPING_TYPE_LINEAR);
    knobFrequency.setRealValue(1700000);
    for(mammut::cpufreq::Domain* d : p.mammut.getInstanceCpuFreq()->getDomains()){
        EXPECT_EQ(d->getCurrentGovernor(), GOVERNOR_USERSPACE);
        EXPECT_EQ(d->getCurrentFrequencyUserspace(), (Frequency) 1700000);
    }

    // One domain - linear
    knobCores.setRealValue(4);
    knobHT.setRealValue(2);
    knobMapping.setRealValue(MAPPING_TYPE_LINEAR);
    knobFrequency.setRealValue(1800000);
    EXPECT_EQ(p.mammut.getInstanceCpuFreq()->getDomains()[0]->getCurrentGovernor(), GOVERNOR_USERSPACE);
    EXPECT_EQ(p.mammut.getInstanceCpuFreq()->getDomains()[0]->getCurrentFrequencyUserspace(), (Frequency) 1800000);
    // Second domain should be off
    EXPECT_EQ(p.mammut.getInstanceCpuFreq()->getDomains()[1]->getCurrentGovernor(), GOVERNOR_USERSPACE);
    for(mammut::topology::VirtualCore* vc : p.mammut.getInstanceCpuFreq()->getDomains()[1]->getVirtualCores()){
        EXPECT_FALSE(vc->isHotPlugged());
    }

    // Both domains - interleaved
    knobCores.setRealValue(14);
    knobHT.setRealValue(2);
    knobMapping.setRealValue(MAPPING_TYPE_INTERLEAVED);
    knobFrequency.setRealValue(1700000);
    for(mammut::cpufreq::Domain* d : p.mammut.getInstanceCpuFreq()->getDomains()){
        EXPECT_EQ(d->getCurrentGovernor(), GOVERNOR_USERSPACE);
        EXPECT_EQ(d->getCurrentFrequencyUserspace(), (Frequency) 1700000);
    }

    // One domain - interleaved
    knobCores.setRealValue(4);
    knobHT.setRealValue(2);
    knobMapping.setRealValue(MAPPING_TYPE_INTERLEAVED);
    knobFrequency.setRealValue(1800000);
    for(mammut::cpufreq::Domain* d : p.mammut.getInstanceCpuFreq()->getDomains()){
        EXPECT_EQ(d->getCurrentGovernor(), GOVERNOR_USERSPACE);
        EXPECT_EQ(d->getCurrentFrequencyUserspace(), (Frequency) 1800000);
    }
}

// Global test with strategy for unused virtual cores = SAME
TEST(KnobsTest, GlobalUnusedSame){
    Parameters  p = getParameters("repara");
    p.knobHyperthreadingEnabled = true;
    p.strategyUnusedVirtualCores = STRATEGY_UNUSED_VC_SAME;
    KnobVirtualCores knobCores(p);
    KnobHyperThreading knobHT(p);
    KnobMappingExternal knobMapping(p, knobCores, knobHT);
    knobMapping.setPid(getpid());
    KnobFrequency knobFrequency(p, knobMapping);

    // Both domains - linear
    knobCores.setRealValue(14);
    knobHT.setRealValue(2);
    knobMapping.setRealValue(MAPPING_TYPE_LINEAR);
    knobFrequency.setRealValue(1700000);
    for(mammut::cpufreq::Domain* d : p.mammut.getInstanceCpuFreq()->getDomains()){
        EXPECT_EQ(d->getCurrentGovernor(), GOVERNOR_USERSPACE);
        EXPECT_EQ(d->getCurrentFrequencyUserspace(), (Frequency) 1700000);
    }

    // One domain - linear
    knobCores.setRealValue(4);
    knobHT.setRealValue(2);
    knobMapping.setRealValue(MAPPING_TYPE_LINEAR);
    knobFrequency.setRealValue(1800000);
    EXPECT_EQ(p.mammut.getInstanceCpuFreq()->getDomains()[0]->getCurrentGovernor(), GOVERNOR_USERSPACE);
    EXPECT_EQ(p.mammut.getInstanceCpuFreq()->getDomains()[0]->getCurrentFrequencyUserspace(), (Frequency) 1800000);
    // Second domain should be at lowest frequency
    EXPECT_EQ(p.mammut.getInstanceCpuFreq()->getDomains()[1]->getCurrentGovernor(), GOVERNOR_USERSPACE);
    EXPECT_EQ(p.mammut.getInstanceCpuFreq()->getDomains()[1]->getCurrentFrequencyUserspace(), (Frequency) 1800000);


    // Both domains - interleaved
    knobCores.setRealValue(14);
    knobHT.setRealValue(2);
    knobMapping.setRealValue(MAPPING_TYPE_INTERLEAVED);
    knobFrequency.setRealValue(1700000);
    for(mammut::cpufreq::Domain* d : p.mammut.getInstanceCpuFreq()->getDomains()){
        EXPECT_EQ(d->getCurrentGovernor(), GOVERNOR_USERSPACE);
        EXPECT_EQ(d->getCurrentFrequencyUserspace(), (Frequency) 1700000);
    }

    // One domain - interleaved
    knobCores.setRealValue(4);
    knobHT.setRealValue(2);
    knobMapping.setRealValue(MAPPING_TYPE_INTERLEAVED);
    knobFrequency.setRealValue(1800000);
    for(mammut::cpufreq::Domain* d : p.mammut.getInstanceCpuFreq()->getDomains()){
        EXPECT_EQ(d->getCurrentGovernor(), GOVERNOR_USERSPACE);
        EXPECT_EQ(d->getCurrentFrequencyUserspace(), (Frequency) 1800000);
    }
}

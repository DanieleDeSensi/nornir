/**
 *  Different tests on parallel for.
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

void runTest(int startloop, int endloop, int step, int chunksize, uint loopDuration = 0){
    int nworkers = 4;
    std::vector<uint> v;
    nornir::Parameters p = getParameters("repara");
    p.requirements.throughput = NORNIR_REQUIREMENT_MAX;
    nornir::ParallelFor pf(nworkers, &p);
    size_t iterations = 10;
    if(loopDuration){
        iterations = 1;
    }
    for(size_t i = 0; i < iterations; i++){
        v.resize(endloop - startloop, 0);

        pf.parallel_for(startloop, endloop, step, chunksize,
        [&](long long int idx, long long int id){
            if(loopDuration){
                usleep((loopDuration*1000000) / ((endloop - startloop) * nworkers));
            }
            v[idx] = idx;
        });

        for(int j = startloop; j < endloop; j += step){
            if((j - startloop) % step == 0){
                EXPECT_EQ(v[j], (uint) j);
            }else{
                EXPECT_EQ(v[j], (uint) 0);
            }
        }
        v.clear();
    }
}

TEST(ParallelForTest, Simple){
    runTest(0, 100, 1, 0);
}

TEST(ParallelForTest, Step3){
    runTest(0, 100, 3, 0);
}

TEST(ParallelForTest, DynamicChunk1){
    runTest(0, 100, 3, 1);
}

TEST(ParallelForTest, DynamicChunk3){
    runTest(0, 100, 3, 3);
}

TEST(ParallelForTest, LongLoop10Seconds){
    runTest(0, 100, 1, 0, 10);
}
/**
 *  Different tests on Manager.
 **/
#include "parametersLoader.hpp"
#include <algorithm>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <linux/limits.h>
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
using std::vector;
using std::string;

string getResultsPath(const string& archName){
    return "./validationdata/" + archName + "/gcc-pthreads-nornir/"; // TODO Different confs??
}

string getBenchPath(const string& archName, const string& benchName){
    return getResultsPath(archName) + benchName;
}

vector<string> getBenchmarks(const string& archName){
    return mammut::utils::getFilesNamesInDir(getResultsPath(archName), false, true);
}

vector<string> getTestCases(const string& archName, const string& benchName){
    vector <string> r;
    for(std::string contract : mammut::utils::getFilesNamesInDir(getBenchPath(archName, benchName), false, true)){
        for(std::string selector : mammut::utils::getFilesNamesInDir(getBenchPath(archName, benchName) + "/" + contract+ "/", false, true)){
            for(std::string testcase : mammut::utils::getFilesNamesInDir(getBenchPath(archName, benchName) + "/" + contract + "/" + selector + "/", true, false)){
                if(testcase.find("parameters.xml") != string::npos){
                    std::vector<std::string> tokens = split(testcase, '.');
                    std::string testcasePrefix = tokens.at(0) + "." + tokens.at(1) + ".";
                    // It will insert something like results/..../10.0.
                    r.push_back(getBenchPath(archName, benchName) + "/" + contract + "/" + selector + "/" + testcasePrefix);
                }
            }
        }
    }
    return r;
}

uint getNumThreads(const string& archName, const string& benchName){
    vector<string> lines = mammut::utils::readFile(getBenchPath(archName, benchName) + "/numthreads.csv");
    return mammut::utils::stringToUint(lines[0]);
}

inline vector<string> getFieldsToCompare(){
    return {"[VirtualCores]", "Workers", "Frequency"};
}

size_t getPosition(const string& fieldName, const string& headerLine){
    vector<string> headerFields = mammut::utils::split(headerLine, '\t');
    auto it = std::find(headerFields.begin(), headerFields.end(), fieldName);
    if (it == headerFields.end()){
        throw std::runtime_error("Field " + fieldName + " not present.");
    }else{
        return std::distance(headerFields.begin(), it);
    }
}

/**
 * We get as input the monitored data, and we verify that the
 * algorithm take the same reconfiguration decision we have in
 * the stored sample files.
 */
TEST(ManagerTest, GlobalTest) {
    // Stored samples were taken with 4 custom fileds.
    EXPECT_EQ(RIFF_MAX_CUSTOM_FIELDS, 4);
    for(string arch : getTestingArchitectures()){ // For all architectures
        for(string bench : getBenchmarks(arch)){ // For all benchmarks
            for(string testcase : getTestCases(arch, bench)){ // And for all testcases on that benchmark
                std::cout << "Running test: " << testcase << std::endl;
                Parameters p = getParameters(arch, testcase + "parameters.xml");
                if(testcase.find("leo") != std::string::npos){
                    // For leo, fix the paths.
                    std::string tmpPath = testcase;
                    // Remove the '10.0.' part
                    while(tmpPath.back() != '/'){
                        tmpPath.pop_back();
                    }
                    p.leo.powerData = tmpPath + "leo_power_raw.csv";
                    p.leo.throughputData = tmpPath + "leo_performance_raw.csv";
                    p.leo.namesData = tmpPath + "leo_names.csv";
                    // p.leo.applicationName doesn't need to be modified
                }
                // Sampling can be done faster since it is just simulated.
                p.samplingIntervalCalibration = 1;
                p.samplingIntervalSteady = 1;

                // Simulate the execution
                ManagerTest* m = new ManagerTest(p, getNumThreads(arch, bench));
                m->setSimulationParameters(testcase + "samples.csv");
                m->start();
                m->join();
                delete m; // To force destructor call for rollback
                std::cout << "Test completed: " << std::endl;

                /*********************************************************/
                /* Now compare output of the test with the expected one. */
                /*********************************************************/
                vector<string> expectedLines = mammut::utils::readFile(testcase + "stats.csv");
                vector<string> realLines = mammut::utils::readFile("stats.csv");
                // Remove empty lines at end of files
                while(expectedLines.back().compare("") == 0){expectedLines.pop_back();}
                while(realLines.back().compare("") == 0){realLines.pop_back();}
                EXPECT_EQ(expectedLines.size(), realLines.size());
                // Get Fields positions (may be different in different logs).
                // We try to be independent from the position so we keep a map where
                // the first element is the field name and the second is its column
                // position inside the log.
                std::map<string, size_t> expectedFields, realFields;
                for(string field : getFieldsToCompare()){
                    expectedFields.insert(std::pair<string, size_t>(field, getPosition(field, expectedLines[0])));
                    realFields.insert(std::pair<string, size_t>(field, getPosition(field, realLines[0])));
                }
                // Now compare the output line by line.
                for(size_t i = 1; i < expectedLines.size(); i++){
                    for(string field : getFieldsToCompare()){
                        // Get the fields by name (to be independent from the column position).
                        vector<string> currentExpectedFields = mammut::utils::split(expectedLines[i], '\t');
                        vector<string> currentRealFields = mammut::utils::split(realLines[i], '\t');
                        string currentExpectedField = currentExpectedFields[expectedFields[field]];
                        string currentRealField = currentRealFields[realFields[field]];
                        // For the cores string we need to remove the positions related to service nodes
                        // (the test may be executed with a farm but the validation is done with managertest that
                        // doesn't have service nodes).
                        if(field.compare("[VirtualCores]") == 0){
                            // Remove square brackets
                            currentExpectedField = currentExpectedField.substr(1, currentExpectedField.size() - 2);
                            currentRealField = currentRealField.substr(1, currentRealField.size() - 2);
                            vector<string> expectedIdentifiers = mammut::utils::split(currentExpectedField, ',');
                            vector<string> realIdentifiers = mammut::utils::split(currentRealField, ',');
                            EXPECT_EQ(expectedIdentifiers.size(), realIdentifiers.size());
                            // Compare all the identifiers except those of service nodes.
                            for(size_t i = 0; i < expectedIdentifiers.size(); i++){
                                EXPECT_EQ(expectedIdentifiers[i], realIdentifiers[i]);
                            }
                        }else{
                            // Compare a generic field.
                            EXPECT_EQ(currentExpectedField, currentRealField);
                        }
                    }
                }
            }
        }
    }
}

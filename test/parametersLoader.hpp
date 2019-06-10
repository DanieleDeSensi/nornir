/**
 * Includes routine to load parameters for different architectures.
 **/
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <linux/limits.h>
#include <time.h>
#include "../src/nornir.hpp"

inline std::vector<std::string> getTestingArchitectures(){
    return {"repara"};
}

inline nornir::Parameters getParameters(const std::string& archName,
                                        const std::string& fileName = ""){
    //p.strategyUnusedVirtualCores // For knobFrequency
    //p.isolateManager //For knobMapping
    //p.knobHyperthreadingEnabled // For knobCores
    //p.disallowedNumCores // For knobCores
    char resolved_path[PATH_MAX];
    std::string relativePath = "./archconfig/" + archName + "/";
    if(realpath(relativePath.c_str(), resolved_path) == NULL){
        throw std::runtime_error("realpath not working for path " + relativePath);
    }
    // We set this env variable to load the nornir config for this architecture
    setenv("XDG_CONFIG_DIRS", resolved_path, 1);
    nornir::Parameters* p;
    if(fileName.compare("") == 0){
        p = new nornir::Parameters();
    }else{
        p = new nornir::Parameters(fileName);
    }
    mammut::SimulationParameters simulationParameters;
    simulationParameters.sysfsRootPrefix = "../src/external/mammut/test/archs/" + archName;
    p->mammut.setSimulationParameters(simulationParameters);
    return *p;
}

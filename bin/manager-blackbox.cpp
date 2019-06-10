/*
 * manager-external.cpp
 *
 * Created on: 21/06/2016
 *
 * This executable starts a manager which can monitor applications not written
 * with the Nornir framework.
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

#include <vector>
#include <iostream>
#include <fstream>
#include <sys/mman.h>
#include <sys/ptrace.h>
#include <sys/wait.h>
#include <stdlib.h>
#include "../src/nornir.hpp"
#include <tclap/CmdLine.h>

using namespace nornir;
using namespace std;

#undef DEBUG
#undef DEBUGB

#define DEBUG_MBB

#ifdef DEBUG_MBB
#define DEBUG(x) do { cerr << "[BlackBox Manager] " << x << endl; } while (0)
#define DEBUGB(x) do {x;} while(0)
#else
#define DEBUG(x)
#define DEBUGB(x)
#endif

struct ScheduledProgram;

typedef struct ScheduledProgram{
    friend std::ostream& operator <<(std::ostream& os, const ScheduledProgram& instance);
    double start;
    double minPerfRequired;
    std::vector<std::string> program;
    pid_t pid;
}ScheduledProgram;

std::ostream& operator <<(std::ostream& os, const ScheduledProgram& instance){
    os << "Starts at: [" << instance.start << "] ";
    os << "Minimum performance required: [" << instance.minPerfRequired << "] ";
    os << "Program: [";
    for(size_t i = 0; i < instance.program.size(); i++){
        os << instance.program[i] << " ";
    }
    os << "]";
    return os;
}

void initializeParameters(Parameters& p){
    p.requirements.throughput = NORNIR_REQUIREMENT_MAX;
    p.samplingIntervalCalibration = 100;
    p.samplingIntervalSteady = 1000;
    p.strategySelection = STRATEGY_SELECTION_LEARNING;
    p.strategyPredictionPerformance = STRATEGY_PREDICTION_PERFORMANCE_USLP;
    p.knobMappingEnabled = false;
    p.knobHyperthreadingEnabled = false;
}

static std::string pidToString(pid_t pid){
    stringstream out;
    out << pid;
    return out.str();
}

int main(int argc, char * argv[]){
    std::vector<long> pids;
    std::vector<double> pidPerfs;
    std::vector<ScheduledProgram> scheduledPrograms;
    std::string logDir;
    std::string parametersFile;
    double powerCap;
    double startTime = mammut::utils::getMillisecondsTime();
    try {
        TCLAP::CmdLine cmd("Runs Nornir on already existing processes, coordinating "
                           "their requirements. It is also possible to specify new "
                           "programs to be ran. One parameter between --pid or "
                           "--schedule must be specified.", ' ', "1.0");
        TCLAP::ValueArg<std::string> logArg("l", "logdir", "Directory where the log files will be stored. (default = ./logs)", false, "./logs/", "string", cmd);
        TCLAP::ValueArg<std::string> parametersArg("r", "parameters", "Nornir parameters file. Can be specified only when only one application needs "
                                                   "to be controlled.", false, "", "string", cmd);
        TCLAP::ValueArg<double> capArg("c", "cap", "Power cap, expressed in watts. If 0, no power cap will "
                                       "be applied and we will just maximise performance for each submitted "
                                       "application. (default = 0).", false, 0, "double", cmd);
        TCLAP::MultiArg<long> pidArg("p", "pid", "PID of a running process.", false, "long", cmd);
        TCLAP::MultiArg<double> perfArg("m", "minperf", "Minimum performance required (one value for each PID, in the same "
                                     "order. It is expressed as a percentage of the maximum achievable performance. For "
                                     "example, a value of 80 means that we would like to do not decrease the performance "
                                     "of the application below the 80% of its maximum. 0 means no specific requirement.",
                                     false, "double", cmd);
        TCLAP::ValueArg<std::string> scheduleArg("s", "schedule", "Name of the file containing the schedule of "
                                                  "the programs to be run. The file has multiple lines, each one with "
                                                  "the following syntax:\n"
                                                  "StartTime MinPerformanceRequired ProgramPath ProgramArguments\n"
                                                  "The start time is a relative offset starting from the start of this "
                                                  "executable.\n"
                                                  "MinPerformanceRequired is a value in the range [0, 100] "
                                                  "representing the minimum required performance, expressed as a percentage "
                                                  "of the maximum achievable performance. For example, a value of 80 means "
                                                  "that we would like to do not decrease the performance of the application "
                                                  "below the 80% of its maximum. 0 means no specific requirement.\n"
                                                  "The lines must be sorted increasingly by start time.",
                                                  false, "", "string", cmd);
        cmd.parse(argc, argv);

        /* Validate parameters. */
        if(!pidArg.getValue().size() &&
           !scheduleArg.getValue().compare("")){
            std::cerr << "[ERROR] One between --pid or --schedule must be specified." << std::endl;
            return -1;
        }

        if(pidArg.getValue().size() !=
           perfArg.getValue().size()){
            std::cerr << "[ERROR] You must specify one minimum performance requirement for each pid." << std::endl;
            return -1;
        }

        logDir = logArg.getValue();
        powerCap = capArg.getValue();
        pids = pidArg.getValue();
        pidPerfs = perfArg.getValue();
        parametersFile = parametersArg.getValue();

        if(pids.size()){
            std::cerr << "Sorry, --pid is still not supported. This is mainly due to difficulties in attaching to an already running "
                         "(non child) process to read the counters. Moreover, even if we could attach, we would not be able to monitor "
                         "the already created childrens/threads." << std::endl;
        }
        
        if(!mammut::utils::existsDirectory(logDir)){
            std::cerr << "Specified logging directory doesn't exist. Please create it before running." << std::endl;
            return -1;
        }

        /* Load schedule. */
        std::ifstream scheduleFile(scheduleArg.getValue());
        std::string line;
        while(std::getline(scheduleFile, line)){
            std::vector<std::string> fields;
            fields = mammut::utils::split(line, ' ');
            if(fields.size() < 3){
                std::cerr << "Invalid scheduled line: " << line << ". You must at least "
                             "specify the scheduled time, minimum performance required "
                             "and program name." << std::endl;
                return -1;
            }
            ScheduledProgram sp;
            sp.start = mammut::utils::stringToDouble(fields[0]);
            sp.minPerfRequired = mammut::utils::stringToDouble(fields[1]);
            sp.program = std::vector<std::string>(fields.begin() + 2, fields.end());
            DEBUG("Read schedule: " << sp);
            scheduledPrograms.push_back(sp);
        }
        scheduleFile.close();
    }catch (TCLAP::ArgException &e){
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
        return -1;
    }

    if(parametersFile.compare("") && scheduledPrograms.size() > 1){
        std::cerr << "Parameter file can only be specified when you want to control only one application." << std::endl;
        return -1;
    }

    ManagerMultiConfiguration mmc;
    mmc.powerCap = powerCap;
    ManagerMulti mm(mmc);
    bool multiManagerNeeded = false;
    if(pids.size() + scheduledPrograms.size() > 1){
        multiManagerNeeded = true;
        mm.start();
        DEBUG("Multi manager started.");
    }

    /* First we add the already running pid. */
    for(size_t i = 0; i < pids.size(); i++){
        Parameters p;
        if(parametersFile.compare("")){
            p.load(parametersFile);
        }else{
            initializeParameters(p);
        }
        std::string logPrefix = logDir + "/" + pidToString(pids.at(i));
        p.loggers.push_back(new LoggerFile(logPrefix + "_stats.csv",
                                  logPrefix + "_calibration.csv",
                                  logPrefix + "_summary.csv",
                                  mammut::utils::getMillisecondsTime() - startTime));
        ManagerBlackBox* m = new ManagerBlackBox(pids.at(i), p);
        if(multiManagerNeeded){
            mm.addManager(m, pidPerfs.at(i));
            DEBUG("Added PID: " << pids.at(i) <<
                  " to the multimanager with minPerf " << pidPerfs.at(i));
        }else{
            m->start();
            DEBUG("Started PID: " << pids.at(i));
            m->join();
            delete m;
            delete p.loggers.back();
            p.loggers.pop_back();
        }
    }

    /* Then we start following the schedule. */
    volatile bool* started = (bool*) mmap(NULL, sizeof(volatile bool), PROT_READ | PROT_WRITE,
                           MAP_SHARED | MAP_ANONYMOUS, -1, 0);
    volatile bool* handlerCreated = (bool*) mmap(NULL, sizeof(volatile bool), PROT_READ | PROT_WRITE,
                                  MAP_SHARED | MAP_ANONYMOUS, -1, 0);

    if(started == MAP_FAILED || handlerCreated == MAP_FAILED){
        std::cerr << "mmap failed."<< std::endl;
        return -1;
    }

    *started = false;
    *handlerCreated = false;
    double lastStart = 0;

    for(size_t i = 0; i < scheduledPrograms.size(); i++){
        ScheduledProgram& sp = scheduledPrograms[i];
        sleep(sp.start - lastStart);
        lastStart = sp.start;

        void* mmem = (void*) mmap(NULL, sizeof(ManagerBlackBox), PROT_READ | PROT_WRITE,
                                  MAP_SHARED | MAP_ANONYMOUS, -1, 0);
        if(mmem == MAP_FAILED){
            std::cerr << "mmap failed." << std::endl;
            return -1;
        }

        pid_t pid = fork();
        if(pid == -1){
            std::cerr << "Fork failed." << std::endl;
            return -1;
        }else if(pid){
            /* Father - Manager. */
            sp.pid = pid;
            Parameters p;
            if(parametersFile.compare("")){
                p.load(parametersFile);
            }else{
                initializeParameters(p);
            }
            stringstream out;
            out << sp.start;
            std::string logPrefix = logDir + "/" + out.str() + "_" + mammut::utils::split(sp.program.at(0), '/').back();
            p.loggers.push_back(new LoggerFile(logPrefix + "_stats.csv",
                                      logPrefix + "_calibration.csv",
                                      logPrefix + "_summary.csv",
                                      mammut::utils::getMillisecondsTime() - startTime));
            while(!*started){;}

            ManagerBlackBox* m = new (mmem) ManagerBlackBox(pid, p);
            *handlerCreated = true;
            DEBUG("Started scheduled: " << sp.program[0] << " with pid " << pid << " on manager " << m);

            if(multiManagerNeeded){
                mm.addManager(m, sp.minPerfRequired);
                // Wait for the child to reset the flags.
                while(*handlerCreated){;}
            }else{
                m->start();
                waitpid(pid, NULL, 0);
                m->join();
                m->~ManagerBlackBox(); // Because created with placement new
                delete p.loggers.back();
                p.loggers.pop_back();
            }
        }else{
            /* Child - Application */
            //ptrace(PTRACE_TRACEME, 0, 0, NULL);
            *started = true;
            // If the handler is not created, we would not catch the counters of
            // the threads/processes created by this process
            while(!*handlerCreated){;}

            extern char** environ;
            char** arguments = new char*[sp.program.size() + 1];
            for(size_t i = 0; i < sp.program.size(); i++){
                arguments[i] = new char[sp.program[i].size() + 1];
                strcpy(arguments[i], sp.program[i].c_str());
            }
            arguments[sp.program.size()] = NULL;
            DEBUG("Running program " << sp.program[0]);
            /* Resets flag for next scheduled application. */
            *started = false;
            *handlerCreated = false;
            if(execve(arguments[0], arguments, environ) == -1){
                std::cerr << "Impossible to run the specified executable. Terminating the manager." << std::endl;
                static_cast<ManagerBlackBox*>(mmem)->terminate();
                // Useless since is never executed, inserted only to avoid cppcheck errors.
                for(size_t i = 0; i < sp.program.size(); i++){
                    delete[] arguments[i];
                }
                delete[] arguments;
                return -1;
            }
            // Useless since is never executed, inserted only to avoid cppcheck errors.
            for(size_t i = 0; i < sp.program.size(); i++){
                delete[] arguments[i];
            }
            delete[] arguments;
            return -1; // execve never returns (except in the error case above).
        }
    }

    DEBUG("All specified programs have been added to the monitoring system. Waiting for their termination...");
    /* All programs scheduled, now wait for termination and clean. */
    if(multiManagerNeeded){
        for(size_t i = 0; i < scheduledPrograms.size(); i++){
            waitpid(scheduledPrograms[i].pid, NULL, 0);
            DEBUG(scheduledPrograms[i].pid << " PID terminated (" << scheduledPrograms[i].program << ")");
        }

        for(size_t i = 0; i < pids.size() + scheduledPrograms.size(); i++){
            ManagerBlackBox* m;
            while((m = static_cast<ManagerBlackBox*>(mm.getTerminatedManager())) == NULL){
                sleep(1);
            }
            DEBUG("Manager (" << m << ") terminated.");
            if(mammut::utils::contains(pids, (long) m->getPid())){
                delete m;
            }else{
                m->~ManagerBlackBox(); // Because we used placement new
            }
            //delete m->_p.observer;
        }
    }
}


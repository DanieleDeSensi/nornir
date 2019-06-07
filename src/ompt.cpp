/*
 * ompt.cpp
 *
 * Created on: 29/05/2019
 *
 * Contains hooks to OpenMP to let it communicate with a nornir external manager.
 * ./bin/manger-external needs to be started before starting the OpenMP application.
 * Before starting the OpenMP application, simply LD_PRELOAD nornir.
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
#define __STDC_FORMAT_MACROS 1

#include <nornir/instrumenter.hpp>

#include <omp.h>
#include <ompt.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <sys/types.h>
#include <unistd.h>

#undef DEBUG
#undef DEBUGB

#ifdef DEBUG_OMPT
//#define DEBUG(x) do { cerr << "[Nornir-OMPT] " << x << endl; } while (0)
#define DEBUG(...) do{ fprintf( stderr, __VA_ARGS__ ); } while( false )
#define DEBUGB(x) do {x;} while(0)
#else
#define DEBUG(...) do {} while(0)
#define DEBUGB(x) do {} while(0)
#endif

static nornir::Instrumenter* instr = NULL;
static bool initialized = false;

static inline void init(){
    if(omp_get_thread_num() == 0 && !initialized){
        const char* parameters_file = getenv("NORNIR_OMP_PARAMETERS");
        nornir::Instrumenter* tmpinstr = NULL;
        if(parameters_file){
            tmpinstr = new nornir::Instrumenter(std::string(parameters_file), omp_get_num_threads());
            riff::ApplicationConfiguration ac;
            ac.samplingLengthMs = 0;
            ac.consistencyThreshold = std::numeric_limits<double>::max();
            tmpinstr->setConfiguration(ac);
        }
        initialized = true;
        instr = tmpinstr;
    }
}

static inline void instrument(unsigned long long times = 1){
    init();
    if(instr){
        unsigned int tid = omp_get_thread_num();
        instr->begin(tid);
        instr->end(tid, times);
    }
}

#define register_callback_t(name, type)                       \
    do{                                                           \
    type f_##name = &on_##name;                                 \
    if (ompt_set_callback(name, (ompt_callback_t)f_##name) ==   \
        ompt_set_never)                                         \
        printf("0: Could not register callback '" #name "'\n");   \
    }while(0)

#define register_callback(name) register_callback_t(name, name##_t)

#ifdef DEBUG_OMPT
static const char* ompt_task_status_t_values[] = {
    NULL,
    "ompt_task_complete",
    "ompt_task_yield",
    "ompt_task_cancel",
    "ompt_task_others"
};
#endif

static ompt_get_thread_data_t ompt_get_thread_data;
static ompt_get_unique_id_t ompt_get_unique_id;

static void
on_ompt_callback_task_schedule(
                               ompt_data_t *first_task_data,
                               ompt_task_status_t prior_task_status,
                               ompt_data_t *second_task_data)
{
    DEBUG("[Nornir] %" PRIu64 ": ompt_event_task_schedule: first_task_id=%" PRIu64 ", second_task_id=%" PRIu64 ", prior_task_status=%s=%d\n", ompt_get_thread_data()->value, first_task_data->value, second_task_data->value, ompt_task_status_t_values[prior_task_status], prior_task_status);
    if(prior_task_status == ompt_task_complete){
        DEBUG("[Nornir] %" PRIu64 ": ompt_event_task_end: task_id=%" PRIu64 "\n", ompt_get_thread_data()->value, first_task_data->value);
        instrument();
    }
}

static void
on_ompt_callback_implicit_task(
                               ompt_scope_endpoint_t endpoint,
                               ompt_data_t *parallel_data,
                               ompt_data_t *task_data,
                               unsigned int team_size,
                               unsigned int thread_num)
{
    switch(endpoint)
        {
        case ompt_scope_begin:
            if(task_data->ptr)
                DEBUG("[Nornir] %s\n", "0: task_data initially not null");
            task_data->value = ompt_get_unique_id();
            DEBUG("[Nornir] %" PRIu64 ": ompt_event_implicit_task_begin: parallel_id=%" PRIu64 ", task_id=%" PRIu64 ", team_size=%" PRIu32 ", thread_num=%" PRIu32 "\n", ompt_get_thread_data()->value, parallel_data->value, task_data->value, team_size, thread_num);
            break;
        case ompt_scope_end:
            DEBUG("[Nornir] %" PRIu64 ": ompt_event_implicit_task_end: parallel_id=%" PRIu64 ", task_id=%" PRIu64 ", team_size=%" PRIu32 ", thread_num=%" PRIu32 "\n", ompt_get_thread_data()->value, (parallel_data)?parallel_data->value:0, task_data->value, team_size, thread_num);
            instrument();
            break;
        }
}

static void
on_ompt_callback_chunk(unsigned long long chunk_size){
    DEBUG("[Nornir] Chunk size: %llu\n", chunk_size);
    instrument(chunk_size);
}

int ompt_initialize(
                    ompt_function_lookup_t lookup,
                    ompt_data_t* data)
{
    DEBUG("[Nornir] libomp init time: %f\n", omp_get_wtime() - *(double*)(data->ptr));
    *(double*)(data->ptr) = omp_get_wtime();

    ompt_set_callback_t ompt_set_callback = (ompt_set_callback_t) lookup("ompt_set_callback");
    ompt_get_thread_data = (ompt_get_thread_data_t) lookup("ompt_get_thread_data");
    ompt_get_unique_id = (ompt_get_unique_id_t) lookup("ompt_get_unique_id");

    register_callback(ompt_callback_task_schedule);
    register_callback(ompt_callback_implicit_task);
    register_callback(ompt_callback_chunk);

    return 1; //success
}

static nornir::InstrumenterServer* is;

void ompt_finalize(ompt_data_t* data)
{
    DEBUG("[Nornir] application runtime: %f\n", omp_get_wtime()-*(double*)(data->ptr));
    if(omp_get_thread_num() == 0 && instr){
        instr->terminate();        
        delete instr;
    }
}

ompt_start_tool_result_t* ompt_start_tool(
                                          unsigned int omp_version,
                                          const char *runtime_version)
{
    static double time = 0;
    time = omp_get_wtime();
    static ompt_start_tool_result_t ompt_start_tool_result = {&ompt_initialize,&ompt_finalize,{.ptr=&time}};

    pid_t serverPid = fork();
    if(!serverPid){
      is = new nornir::InstrumenterServer(true);
      is->start();
      is->join();
      exit(0);
    }

    return &ompt_start_tool_result;
}
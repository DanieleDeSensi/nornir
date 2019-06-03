/*
 * instrumenter.h
 *
 * Created on: 11/08/2016
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

#ifndef NORNIR_INSTRUMENTER_H_
#define NORNIR_INSTRUMENTER_H_

#ifdef __cplusplus
extern "C"{
#endif

struct NornirInstrumenter;
typedef struct NornirInstrumenter NornirInstrumenter;

NornirInstrumenter* nornir_instrumenter_create(const char* parametersFile);

NornirInstrumenter* nornir_instrumenter_create_with_threads(const char* parametersFile, size_t numThreads);

void nornir_instrumenter_destroy(NornirInstrumenter* instrumenter);

void nornir_instrumenter_begin(NornirInstrumenter* instrumenter);

void nornir_instrumenter_begin_with_threads(NornirInstrumenter* instrumenter, size_t threadId);

void nornir_instrumenter_end(NornirInstrumenter* instrumenter);

void nornir_instrumenter_end_with_threads(NornirInstrumenter* instrumenter, size_t threadId);

void nornir_instrumenter_terminate(NornirInstrumenter* instrumenter);

unsigned long nornir_instrumenter_get_execution_time(NornirInstrumenter* instrumenter);

unsigned long long nornir_instrumenter_get_total_tasks(NornirInstrumenter* instrumenter);

void nornir_instrumenter_set_total_threads(NornirInstrumenter* instrumenter, unsigned int totalThreads);

void nornir_instrumenter_set_phase_id(NornirInstrumenter* instrumenter, unsigned int phaseId);

void nornir_instrumenter_mark_inconsistent_samples(NornirInstrumenter* instrumenter);

#ifdef __cplusplus
}
#endif

#endif /* NORNIR_INSTRUMENTER_H_ */

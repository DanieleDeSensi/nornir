/*
 * ffincs.hpp
 *
 * Created on: 23/03/2016
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

/**
 * This file includes all the includes to FastFlow files.
 * They are centralized in this file since we want to ensure
 * that they are included after some appropriate macros have been defined.
 * All the nornir files that need to use some specific FastFlow should include
 * this file instead of directly including FastFlow files.
 */
#ifndef NORNIR_FFINCS_HPP_
#define NORNIR_FFINCS_HPP_

#if defined(FF_NODE_HPP) || defined(FF_CONFIG_HPP)
#error "Please include nornir headers BEFORE including FastFlow files."
#endif

#define FF_TASK_CALLBACK
#define TRACE_FASTFLOW
#undef BLOCKING_MODE // Rethreading does not work with blocking mode
#define FF_ESAVER
#define NORNIR_NEW_FF_API_SUPPORT

// This 2 following includes are needed since ff/mapping_utils.hpp does not include them.
//#include <unistd.h>
//#include <sys/syscall.h>

#include "external/fastflow/ff/farm.hpp"

#endif /* NORNIR_FFINCS_HPP_ */

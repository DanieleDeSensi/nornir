/*
 * nornir.hpp
 *
 * Created on: 20/07/2017
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
#ifndef NORNIR_HPP_
#define NORNIR_HPP_

#ifdef __x86_64__
#include "./dataflow/interpreter.hpp"
#endif

#include "interface.hpp"
#include "instrumenter.hpp"
#include "manager.hpp"
#include "manager-multi.hpp"
#include "stats.hpp"

#endif // NORNIR_HPP_

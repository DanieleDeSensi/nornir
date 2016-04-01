/*
 * interpreter.hpp
 *
 * Created on: 26/03/2016
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

#ifndef NORNIR_DF_INTERPRETER_HPP_
#define NORNIR_DF_INTERPRETER_HPP_

#include "../interface.hpp"
#include "interpreter.hpp"
#include "tokens.hpp"
#include "stream.hpp"
#include "mdfi.hpp"
#include "mdfg.hpp"
#include "hashMap.hpp"
#include "skeleton/computable.hpp"
#include "skeleton/farm.hpp"
#include "skeleton/pipeline.hpp"
#include "skeleton/ewc.hpp"
#include "../external/fastflow/ff/ubuffer.hpp"
#include "../external/fastflow/ff/squeue.hpp"
#include <vector>
#include <queue>
#include <typeinfo>
#include <errno.h>

#ifndef MAXPOOLSIZE
#define MAXPOOLSIZE 500000
#endif

namespace nornir{
namespace dataflow{

/**
 * This is the worker of the fastflow's farm.
 */
class WorkerMdf: public nornir::Worker<Mdfi>{
private:
    ff::uSWSR_Ptr_Buffer* _buffer;
public:
    WorkerMdf(ff::uSWSR_Ptr_Buffer* buffer);

    /**
     * Computes a macro data flow instruction.
     * \param task A macro data flow instruction.
     * \return The result of the computation.
     */
    void compute(Mdfi* t);
};

/**
 * This is the interpreter of the macro data flow instructions.
 */
class Interpreter{
private:
    ff::uSWSR_Ptr_Buffer** _buffers;
    nornir::Scheduler<Mdfi>* _s;
    nornir::Farm<Mdfi>* _farm;
    nornir::Parameters* _p;
    nornir::Observer* _o;
    Mdfg* _compiledGraph;
    size_t _maxWorkers;
public:
    /**
     * Constructor of the interpreter.
     * \param c The skeleton to be executed.
     * \param i The input stream.
     * \param o The output stream.
     * \param parDegree Number of workers to activate.
     * \param groupSize Number of instructions to send to the interpreter
     *                  before waiting the results.
     * \param orderedTasks If true, we will guarantee that a macro data flow
     *                     instruction generated by the stream element N + 1, will
     *                     always be executed before after the same instruction
     *                     generated by the stream element N terminated.
     */
    Interpreter(Computable *c, InputStream *i, OutputStream *o, size_t parDegree,
                unsigned long int groupSize, bool orderedTasks);

    /**
     * Constructor of the interpreter.
     * \param graph The macro data flow graph to be executed.
     * \param i The input stream.
     * \param o The output stream.
     * \param parDegree Number of workers to activate.
     * \param groupSize Number of instructions to send to the interpreter
     *                  before waiting the results.
     * \param orderedTasks If true, we will guarantee that a macro data flow
     *                     instruction generated by the stream element N + 1, will
     *                     always be executed before after the same instruction
     *                     generated by the stream element N terminated.
     */
    Interpreter(Mdfg *graph, InputStream *i, OutputStream *o, size_t parDegree,
                unsigned long int groupSize, bool orderedTasks);

    /**
     * Destructor of the interpreter.
     */
    ~Interpreter();

    /**
     * Prints some stats.
     * \param out The stream where print the stats.
     */
    inline void stats(std::ostream& out=std::cout){
        out << "Stats not available." << std::endl;
    }

    inline void start(){
        _farm->start();
    }

    inline void wait(){
        _farm->wait();
    }
};

}
}

#endif /* NORNIR_DF_INTERPRETER_HPP_ */

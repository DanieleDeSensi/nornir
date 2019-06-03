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
#include "tokens.hpp"
#include "stream.hpp"
#include "mdfi.hpp"
#include "mdfg.hpp"
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

PUSH_WARNING
GCC_DISABLE_WARNING(vla)
#include "../external/queues/MSqueue.hpp"
POP_WARNING

#define QUEUE MSQueue

namespace nornir{
namespace dataflow{

// TODO: ArrayWrapper<T> -> std::vector<T>
// TODO: Use arma::vec and arma::mat for Map and Reduce computations

/**
 * This is the worker of the fastflow's farm.
 */
class WorkerMdf: public nornir::Worker<Mdfi>{
private:
    QUEUE& _q;
    bool _init;
    size_t _qId;
    size_t _processedTasks;
public:
    explicit WorkerMdf(QUEUE& q);

    /**
     * Computes a macro data flow instruction.
     * \param task A macro data flow instruction.
     * \return The result of the computation.
     */
    void compute(Mdfi* t);

    size_t getProcessedTasks() const;
};

/**
 * This is the interpreter of the macro data flow instructions.
 */
class Interpreter: mammut::utils::NonCopyable{
private:
    QUEUE _q;
    nornir::Scheduler<Mdfi>* _s;
    std::vector<WorkerMdf*> _workers;
    nornir::Farm<Mdfi>* _farm;
    nornir::Parameters* _p;
    Mdfg* _compiledGraph;
    size_t _maxWorkers;
public:
    /**
     * Constructor of the interpreter.
     * \param p The nornir parameters.
     * \param c The skeleton to be executed.
     * \param i The input stream.
     * \param o The output stream. If NULL, results are not sent on the output
     *          stream.
     */
    Interpreter(Parameters* p, Computable *c, InputStream *i, OutputStream *o = NULL);

    /**
     * Constructor of the interpreter.
     * \param p The nornir parmeters.
     * \param graph The macro data flow graph to be executed.
     * \param i The input stream.
     * \param o The output stream. If NULL, results are not sent on the output
     *          stream.
     */
    Interpreter(Parameters* p, Mdfg *graph, InputStream *i, OutputStream *o = NULL);

    /**
     * Destructor of the interpreter.
     */
    ~Interpreter();

    /**
     * Prints some stats.
     * \param out The stream where print the stats.
     */
    inline void stats(std::ostream& out = std::cout){
        _farm->stats(out);
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

/*
 * computable.hpp
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

#ifndef NORNIR_DF_COMPUTABLE_HPP_
#define NORNIR_DF_COMPUTABLE_HPP_

#include "../stream.hpp"
#include <map>

namespace nornir{
namespace dataflow{

/**
 * \class Computable
 * This class represents a generic computation.
 * If a class extends \e Computable, its "compute" function can be computed sequentially in a Macro Data Flow Instruction.
 * \e Compute is a function that takes an array of Task* and returns an array of Task*. The array must be dynamically
 * allocated using new[].
 *
 * This is an example of a function that takes \e n values as input and returns \e m values.
 *
 * \code
 * #include "nornir/dataflow/skeleton/computable.hpp"
 *
 * TODO: FIX documentation
 * class myCompute: public skel::Computable{
 * public:
 * Task** compute(Task** input){
 *      IntTask* x;
 *      //Creation of the array to return.
 *      //(To note that the array is created using new[].
 *      Task** toReturn=new Task*[n];
 *      for(int i=0; i<n; i++){
 *          //Gets the input task (to note the cast).
 *          x=(IntTask*) input[i];
 *          //Insertion of the elements to return.
 *          //toReturn[i] will be sent to the MDFi linked
 *          //to the i-th output.
 *          toReturn[i]=new IntArrayTask(x[i]->get()+10);
 *      }
 *      return toReturn;
 *  }
 *};
 * \endcode
 */

class Mdfi;
class Farm;
class Pipeline;
class EmitterWorkerCollector;
class Computable;

class Data{
private:
    friend class Mdfi;
    friend class Farm;
    friend class Pipeline;
    friend class EmitterWorkerCollector;
    std::map<Computable*, void*> _sources;
    std::map<Computable*, void**> _destinations;

    /**
     * NULL = Input stream.
     */
    void setSource(void* s, Computable* c = NULL){
        //cppcheck-suppress nullPointerDefaultArg
        _sources[c] = s;
    }

    /**
     * NULL = Output stream.
     */
    void setDestination(void** d, Computable* c = NULL){
        //cppcheck-suppress nullPointerDefaultArg
        _destinations[c] = d;
    }
public:
    /**
     * Retreive the data received from a specific computable.
     * @param c The computable from which the data should be received.
     *          If NULL, the data is received from the input stream.
     *          For map computations, c is the index of the worker
     *          (cast from int to Computable*).
     */
    void* getInput(Computable* c = NULL){
        if(!c){
            if(_sources.size() != 1){
                throw std::runtime_error("You can avoid to specify the instruction from" \
                                         "which you want to receive the data only if this" \
                                         "instruction has only 1 input.");
            }else{
                return _sources.begin()->second;
            }
        }else{
            if(_sources.find(c) == _sources.end()){
                throw std::runtime_error("Impossible to receive from computable.");
            }
            return _sources[c];
        }
    }

    /**
     * Sends a result to a specific computable.
     * @param c The computable to which the data should be sent. If NULL,
     *          the data is sent to the output stream. For map computations,
     *          c is the index of the worker (cast from int to Computable*).
     */
    void setOutput(void* x, Computable* c = NULL){
        if(!c){
            if(_destinations.size() != 1){
                throw std::runtime_error("You can avoid to specify the instruction to" \
                                         "which you want to send the data only if this" \
                                         "instruction has only 1 output.");
            }else{
                *(_destinations.begin()->second) = x;
            }
        }else{
            if(_destinations.find(c) == _destinations.end()){
                throw std::runtime_error("Impossible to send to computable.");
            }
            *(_destinations[c]) = x;
        }
    }

};

class Computable{

public:
    /**
     * Destructor of the computable.
     */
    virtual ~Computable(){;}

    /**
     * This method computes the result.
     */
    virtual void compute(Data* d) = 0;
};

}
}

#endif /* NORNIR_DF_COMPUTABLE_HPP_ */

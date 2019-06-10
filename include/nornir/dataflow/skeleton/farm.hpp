/*
 * farm.hpp
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

#ifndef NORNIR_DF_FARM_HPP_
#define NORNIR_DF_FARM_HPP_

namespace nornir{
namespace dataflow{

/**
 * \class Farm
 * This class is the farm skeleton class.
 */
class Farm: public Computable{
private:
    Computable* worker;
    bool deleteWorker;
public:
    /**
     * The constructor of the farm.
     * \param w Is the skeleton used to implement the farm workers.
     *          This computation is indipendently performed each time that a
     *          task is submitted to the farm.
     * \param deleteWorker If true, the destructor of the farm deletes the worker.
     */
    Farm(Computable* w, bool deleteWorker = false):worker(w), deleteWorker(deleteWorker){;}

    /**
     * Destructor of the farm.
     * If deleteAll is true, deletes the worker of the farm.
     */
    ~Farm(){
        if(deleteWorker){
            delete worker;
        }
    }
    /**
     * This method computes the result sequentially.
     */
    void compute(Data* d){
        void* result;
        Data dw;
        dw.setSource(d->getInput());
        dw.setDestination(&result);
        worker->compute(&dw);
        d->setOutput(result);
    }

    /**
     * Returns a pointer to the worker.
     * \return A pointer to the worker.
     */
    Computable* getWorker(){
        return worker;
    }
};


template <typename T, typename V, V*(*fun)(T*)> class StandardFarmWorker: public Computable{
public:
    void compute(Data* d){
        V* result = fun(static_cast<T*>(d->getInput()));
        d->setOutput(result);
    }
};

/**
 * Creates a standard farm.
 *
 * \tparam T* is the type of the input elements.
 * \tparam V* is the type of the output elements.
 * \tparam fun is the function to compute over the elements.
 *
 * \return A pointer to a standard farm.
 */

template<typename T, typename V, V*(*fun)(T*)> Farm* createStandardFarm(){
    return new Farm(new StandardFarmWorker<T, V, fun>, true);
}

}
}


#endif /* NORNIR_DF_FARM_HPP_ */

/*
 * ewc.hpp
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

#ifndef NORNIR_DF_EWC_HPP_
#define NORNIR_DF_EWC_HPP_

#include "computable.hpp"
#include "../stream.hpp"

namespace nornir{
namespace dataflow{

class EmitterWorkerCollector;

class Scatterer: public Computable{
    friend Mdfg* compile(Computable* c);
private:
    std::vector<Computable*> _workers;
protected:
    size_t _numPartitions;
public:
    explicit Scatterer(size_t numPartitions):_numPartitions(numPartitions){;}

    virtual std::vector<void*> compute(void* in) = 0;

    void compute(Data* d){
        std::vector<void*> r = compute(d->getInput());
        for(size_t i = 0; i < _numPartitions; i++){
            d->setOutput(r.at(i), _workers.at(i));
        }
    }
};

class Gatherer: public Computable{
    friend Mdfg* compile(Computable* c);
private:
    std::vector<Computable*> _workers;
protected:
    size_t _numPartitions;
public:
    explicit Gatherer(size_t numPartitions):_numPartitions(numPartitions){;}

    virtual void* compute(std::vector<void*> in) = 0;

    void compute(Data* d){
        std::vector<void*> inputs;
        for(size_t i = 0; i < _numPartitions; i++){
            inputs.push_back(d->getInput(_workers.at(i)));
        }
        d->setOutput(compute(inputs));
    }
};

/**
 * A generic emitter of a reduce.
 * \tparam T* is the type of the elements of the array (The input stream must produce ArrayWrapper<T*>).\n
 *         T must be a subclass of Task.
 */
template <typename T> class ReduceScatterer: public Scatterer{
private:
    bool _autoDelete;
public:
    /**
     * Constructor of the emitter of the reduce.
     * \param numPartitions The number of partitions.
     * \param autoDelete If true, the elements received from the input stream (ArrayWrapper<T*>) are automatically deleted.
     */
    ReduceScatterer(size_t numPartitions, bool autoDelete=true);

    std::vector<void*> compute(void* in);
};

/**
 * A generic worker of a reduce.
 *
 * \tparam T* is the type of the elements of the array.
 * \tparam fun Is a binary, associative and commutative function to compute over the elements.
 */
template <typename T, T*(*fun)(T*,T*)> class ReduceWorker: public Computable{
public:
    void compute(Data* d);
};


/**
 * A generic collector of the reduce.
 *
 * \tparam T* is the type of the elements of the array.
 * \tparam fun Is a binary, associative and commutative function to compute over the elements.
 * The map returns a T*.
 */
template <typename T, T*(*fun)(T*,T*)> class ReduceGatherer: public Gatherer{
public:
    /**
     * Constructor of the collector of the reduce.
     * \param numPartitions The number of partitions.
     */
    explicit ReduceGatherer(size_t numPartitions);

    void* compute(std::vector<void*> in);
};

/**
 * A generic emitter of a map.
 *
 * \tparam T* is the type of the elements of the array.(The input stream must produce ArrayWrapper<T*>).
 */
template <typename T> class MapScatterer: public Scatterer{
private:
    bool _autoDelete;
public:
    /**
     * Constructor of the emitter of the map.
     * \param numPartitions The number of partitions.
     * \param autoDelete If true, the elements received from the input stream (ArrayWrapper<T*>) are automatically deleted.
     */
    MapScatterer(size_t numPartitions, bool autoDelete = true);

    std::vector<void*> compute(void* in);
};

/**
 * Generic worker of a map.
 *
 * \tparam T* is the type of the input elements. T must be subclass of Task.
 * \tparam V* is the type of the output elements. V must be subclass of Task.
 * \tparam fun is the function to compute over the elements.
 *
 */
template <typename T, typename V, V*(*fun)(T*) > class MapWorker: public Computable{
public:
    void compute(Data* d);
};

/**
 * Generic collector of a map.
 *
 * \tparam V* is the type of the elements received from the workers.
 *
 * The \e map returns an ArrayWrapper<V*>
 */
template <typename V> class MapGatherer: public Gatherer{
public:
    explicit MapGatherer(size_t numPartitions);

    void* compute(std::vector<void*> in);
};


/**
 * This class represents a standard pattern with an Emitter, N workers and a Collector.
 * This pattern is used in Map and Reduce skeletons.
 */
class EmitterWorkerCollector: public Computable{
private:
    size_t _nWorkers;
    Scatterer *_scatterer;
    std::vector<Computable*> _workers;
    Gatherer *_gatherer;
    bool deleteAll;
public:
    /**
     * Constructor skeleton.
     * \param scatterer The skeleton used to implement the splitting of the data\n.
     *                Scatterer can be an instance of the generic scatterers
     *                (MapScatterer or ReduceScatterer) or can be redefined for the
     *                purpose simply extending the Scatterer
     *                class.
     * \param workers The computables used to implement the map workers.\n
      *               Each worker can be an instance of the generic worker
     *                (MapWorker or ReduceWorker) or can be redefined for the
     *                purpose simply extending the \e Computable class.
     * \param gatherer The skeleton used to implement the merging of the data.\n
     *                  Gatherer can be an instance of the generic collector
     *                  (MapGatherer or ReduceGatherer) or can be redefined
     *                  for the purpose simply extending the Gatherer class.
     * \param deleteAll If true, the destructor deletes the emitter, the worker
     *                  and the collector.
     */
    inline EmitterWorkerCollector(Scatterer* scatterer, std::vector<Computable*> workers,
                                  Gatherer* gatherer, bool deleteAll = false):
                               _nWorkers(workers.size()), _scatterer(scatterer),
                               _workers(workers), _gatherer(gatherer), deleteAll(deleteAll){
        ;
    }

    /**
     * Destructor of the EmitterWorkerCollector.
     * If deleteAll is true, deletes emitter, workers and collector.
     */
    inline ~EmitterWorkerCollector(){
        if(deleteAll){
            delete _scatterer;
            for(size_t i = 0; i < _workers.size(); i++){
                delete _workers.at(i);
            }
            delete _gatherer;
        }
    }

    /**
     * This method computes the result sequentially.
     */
    inline void compute(Data* d){
        void *fromWorker;

        std::vector<void*> eRes = _scatterer->compute(d->getInput());
        std::vector<void*> cInput;

        Data dw;
        for(uint i = 0; i < _nWorkers; i++){
            dw.setSource(eRes.at(i));
            dw.setDestination(&fromWorker);
            _workers.at(i)->compute(&dw);
            cInput.push_back(fromWorker);
        }

        d->setOutput(_gatherer->compute(cInput));
    }

    /**
     * Returns the number of workers.
     * \return Number of workers.
     */
    inline uint getNWorkers(){
        return _nWorkers;
    }

    /**
     * Returns a pointer to the scatterer.
     * \return A pointer to the scatterer.
     */
    inline Scatterer* getScatterer(){
        return _scatterer;
    }

    /**
     * Returns a pointer to the worker.
     * \return A pointer to the worker.
     */
    inline std::vector<Computable*>& getWorkers(){
        return _workers;
    }

    /**
     * Returns a pointer to the gatherer.
     * \return A pointer to the gatherer.
     */
    inline Gatherer* getGatherer(){
        return _gatherer;
    }
};

/**
 * Creates a standard reduce.
 *
 * \tparam T is the type of the elements of the array (The input stream must produce ArrayWrapper<T*>)
 * \tparam fun Is a binary, associative and commutative function to compute over the elements.
 *
 * \param numPartitions Number of partitions.
 * \param autoDelete If true, the elements received from the input stream (ArrayWrapper<T*>) are automatically deleted by the emitter.
 *
 * \return A pointer to a standard reduce.
 */

template<typename T, T*(*fun)(T*,T*)>
EmitterWorkerCollector* createStandardReduce(size_t numPartitions, bool autoDelete = true){
    std::vector<Computable*> workers;
    for(size_t i = 0; i < numPartitions; i++){
        workers.push_back(new ReduceWorker<T,fun>);
    }
    return new EmitterWorkerCollector(new ReduceScatterer<T>(numPartitions, autoDelete),
                                      workers,
                                      new ReduceGatherer<T,fun>(numPartitions), true);
}

/**
 * Creates a standard map.
 *
 * \tparam T is the type of the input elements (The input stream must produce ArrayWrapper<T*>).
 * \tparam V is the type of the output elements (The map returns ArrayWrapper<V*>).
 * \tparam fun is the function to compute over the elements.
 *
 * \param numPartitions Number of partitions.
 * \param autoDelete If true, the elements received from the input stream (ArrayWrapper<T*>) are automatically deleted by the emitter.
 *
 * Input array (single element of the input stream): [x1,x2,...,xn] where typeof(x1)==typeof(x2)==...==typeof(xn)==T\n
 * Output array (single element of the output stream): [y1=fun(x1),y2=fun(x2),...,yn=fun(xn)] where typeof(y1)==typeof(y2)==...==typeof(yn)==V\n
 *
 *
 * \return A pointer to a standard map.
 */

template<typename T, typename V, V*(*fun)(T*)>
EmitterWorkerCollector* createStandardMap(size_t numPartitions, bool autoDelete=true){
    std::vector<Computable*> workers;
    for(size_t i = 0; i < numPartitions; i++){
        workers.push_back(new MapWorker<T, V, fun>);
    }
    return new EmitterWorkerCollector(new MapScatterer<T>(numPartitions, autoDelete),
                                      workers,
                                      new MapGatherer<T>(numPartitions), true);
}

}
}

#include "ewc.tpp"

#endif /* NORNIR_DF_EWC_HPP_ */

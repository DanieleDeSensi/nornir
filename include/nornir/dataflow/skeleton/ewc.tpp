/*
 * ewc.tpp
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

namespace nornir{
namespace dataflow{

template <typename T>
ReduceScatterer<T>::ReduceScatterer(size_t numPartitions, bool autoDelete):
    Scatterer(numPartitions), _autoDelete(autoDelete){;}

template <typename T>
std::vector<void*> ReduceScatterer<T>::compute(void* in){
    ArrayWrapper<T*>* task = (ArrayWrapper<T*>*) in;
    std::vector<void*> r;
    int dim = task->size();
    int mod = dim%_numPartitions;
    int size = dim/_numPartitions;
#ifdef NOCOPY
    int h = 0;
    for(uint i = 0; i < chunkNum; i++){
        int l = h;
        if(mod && i == chunkNum - mod){
            size += 1;
            mod = 0;
        }
        h = l + size;
        toRet[i] = new ArrayIndexes<T*>(task, l, h);
    }
#else
    int k=0;
    for(uint i=0; i<_numPartitions; i++){
        if(mod && i==_numPartitions-mod){
            size+=1;
            mod=0;
        }
        ArrayWrapper<T*>* toAdd=new ArrayWrapper<T*>(size);
        for(int j=0; j<size; j++){
            toAdd->set(j,task->get(k));
            k++;
        }
        r.push_back(toAdd);
    }
    if(_autoDelete) delete task;
#endif
    return r;
}


template <typename T, T*(*fun)(T*,T*)>
void ReduceWorker<T, fun>::compute(Data* d){
    void* t = d->getInput();
    ArrayWrapper<T*>* w1;
    int nElem,first;
#ifdef NOCOPY
    ArrayIndexes<T*>* ai=(ArrayIndexes<T*>*)t;
    w1=ai->getArray();
    nElem=ai->getj();
    first=ai->geti();
#else
    void* toDelete = t;
    w1=(ArrayWrapper<T*>*)t;
    nElem=w1->size();
    first=0;
#endif
    if(nElem-first==1){
#ifndef NOCOPY
        t = w1->get(first);
        delete toDelete;
#endif
    }else{
        T *a=w1->get(first),*b=w1->get(first+1);
        T *x=fun(a,b);
        for(int i=first+2; i<nElem;i++){
            a=w1->get(i);
            x=fun(x,a);
        }
#ifdef NOCOPY
    w1->set(first,x);
#else
    t = x;
    delete toDelete;
#endif
    }
    d->setOutput(t);
}


template<typename T, T*(*fun)(T*,T*)>
ReduceGatherer<T, fun>::ReduceGatherer(size_t numPartitions):
    Gatherer(numPartitions){;}

template<typename T, T*(*fun)(T*,T*)>
void* ReduceGatherer<T, fun>::compute(std::vector<void*> in){
    assert(in.size() > 1);
    T *p, *q;
#ifdef NOCOPY
    ArrayIndexes<T*>* ai;
    ai = dynamic_cast<ArrayIndexes<T*>*>(d->getInput(dynamic_cast<Computable*>(0)));
    p = ai->getArray()->get(ai->geti());
    delete ai;
    ai = dynamic_cast<ArrayIndexes<T*>*>(d->getInput(dynamic_cast<Computable*>(1)));
    q = ai->getArray()->get(ai->geti());
    delete ai;
#else
    p = static_cast<T*>(in.at(0));
    q = static_cast<T*>(in.at(1));
#endif
    T* x = fun(p,q);
    for(int i = 2; i < _numPartitions; i++){
#ifdef NOCOPY
        ArrayIndexes<T*>* ai;
        ai = (ArrayIndexes<T*>*)t[i];
        T* a = ai->getArray()->get(ai->geti());
        if(i == chunkNum - 1){
            delete ai->getArray();
        }
        delete ai;
#else
        T* a = static_cast<T*>(in.at(i));
#endif
        x = fun(x, a);
    }
    return x;
}

template <typename T>
MapScatterer<T>::MapScatterer(size_t numPartitions, bool autoDelete):
    Scatterer(numPartitions), _autoDelete(autoDelete){;}

template <typename T>
std::vector<void*> MapScatterer<T>::compute(void* in){
    ArrayWrapper<T*>* task = (ArrayWrapper<T*>*) in;
    std::vector<void*> r;
    uint dim = task->size();
    uint mod = dim % _numPartitions;
    uint size = dim / _numPartitions;
#ifdef NOCOPY
    int h = 0;
    for(uint i = 0; i < numWorkers; i++){
        int l = h;
        if(mod && i == numWorkers-mod){
            size += 1;
            mod = 0;
        }
        h = l+size;
        toRet[i] = new ArrayIndexes<T*>(task,l,h);
    }
#else
    int k = 0;
    for(uint i = 0; i < _numPartitions; i++){
        if(mod && i == _numPartitions - mod){
            size += 1;
            mod = 0;
        }
        ArrayWrapper<void*>* toAdd = new ArrayWrapper<void*>(size);
        for(uint j = 0; j < size; j++){
            toAdd->set(j,task->get(k));
            k++;
        }
        r.push_back(toAdd);
    }
    if(_autoDelete){
        delete task;
    }
#endif
    return r;
}

template <typename T, typename V, V*(*fun)(T*) >
void MapWorker<T, V, fun>::compute(Data* d){
    ArrayWrapper<void*> *w1;
    int nElem,first;
    void* t = d->getInput();
#ifdef NOCOPY
    ArrayIndexes<void*>* ai=(ArrayIndexes<void*>*) t;
    w1=ai->getArray();
    nElem=ai->getj();
    first=ai->geti();
#else
    w1=(ArrayWrapper<void*>*) t;
    nElem=w1->size();
    first=0;
#endif
    for(int i=first; i<nElem;i++){
        T* x = static_cast<T*>(w1->get(i));
        V* y = fun(x);
        w1->set(i, y);
    }
    d->setOutput(t);
}

template <typename V>
MapGatherer<V>::MapGatherer(size_t numPartitions):Gatherer(numPartitions){;}

template <typename V>
void* MapGatherer<V>::compute(std::vector<void*> in){
#ifdef NOCOPY
    ArrayIndexes<void*>* ai;
    for(uint i=0; i<_numPartitions-1; i++){
        ai = dynamic_cast<ArrayIndexes<void*>*>(receiveData(dynamic_cast<Computable*>(i)));
        delete ai;
    }
    ai = dynamic_cast<ArrayIndexes<void*>*>(receiveData(dynamic_cast<Computable*>(_numPartitions-1)));
    sendData(ai->getArray());
    delete ai;
#else
    uint size = 0;
    for(uint i = 0; i < _numPartitions; i++){
        size += ((ArrayWrapper<void*>*) in.at(i))->size();
    }


    ArrayWrapper<V*> *aw = new ArrayWrapper<V*>(size), *tempAw;
    uint tempSize, k = 0;
    for(uint i = 0; i < _numPartitions; i++){
        tempAw = ((ArrayWrapper<V*>*) in.at(i));
        tempSize = tempAw->size();
        for(uint j = 0; j < tempSize; j++){
            aw->set(k,tempAw->get(j));
            ++k;
        }
        delete tempAw;
    }
    return aw;
#endif
}

}
}

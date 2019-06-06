/*
 * map.cpp
 *
 * Created on: 04/05/2016
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

#include <iostream>
#include <stdlib.h>
#include <nornir/nornir.hpp>

using namespace nornir::dataflow;

#define MAX_SIZE 10;

class DemoInputStream: public nornir::dataflow::InputStream{
private:
    size_t _currentElem;
    size_t _streamSize;
    bool _eos;
public:
    explicit inline DemoInputStream(size_t streamSize):
            _currentElem(0), _streamSize(streamSize), _eos(false){
        srand(time(NULL));
    }

    inline void* next(){
        if(_currentElem < _streamSize){
            ++_currentElem;
            size_t size = rand() % MAX_SIZE;
            ArrayWrapper<float*>* aw = new ArrayWrapper<float*>(size);
            std::cout << "Generated: [";
            for(size_t i = 0; i < size; i++){
                float* x = new float();
                *x = (rand() % 100)* 0.1;
                aw->set(i, x);
                std::cout << *x << ", ";
            }
            std::cout << "]" << std::endl;
            return (void*) aw;
        }else{
            _eos = true;
            return NULL;
        }
    }

    inline bool hasNext(){
        return !_eos;
    }

};

class DemoOutputStream: public nornir::dataflow::OutputStream{
public:
    void put(void* a){
        ArrayWrapper<float*>* aw = (ArrayWrapper<float*>*) a;
        std::cout << "Result: ";
        for(size_t i = 0; i < aw->size(); i++){
            std::cout << *(aw->get(i)) << ", ";
        }
        std::cout << std::endl;
    }
};

int* fun(float* x){
    *x = floor(*x);
    //cppcheck-suppress invalidPointerCast
    return (int*) x;
}

int main(int argc, char** argv){
    if(argc < 2){
        std::cerr << "Usage: " << argv[0] << " streamSize [numWorkers]" << std::endl;
        return -1;
    }

    /* Create streams. */
    DemoInputStream inp(atoi(argv[1]));
    DemoOutputStream out;

    uint nWorkers = 4;
    if(argc >= 3){
        nWorkers = atoi(argv[2]);
    }

    Computable* map = createStandardMap<float, int, fun>(nWorkers);

    nornir::Parameters p("parameters.xml");
    nornir::dataflow::Interpreter m(&p, map, &inp, &out);
    m.start();
    m.wait();
}


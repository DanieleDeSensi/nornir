/*
 * pipeline.cpp
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
            int* x;
            ++_currentElem;
            x = new int();
            *x = rand() % 1000;
            std::cout << "Generated: " << *x << std::endl;
            std::cout << "ExpectedOutput: " << ((*x + 3) * 4) + 1 << std::endl;
            return (void*) x;
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
        int* x = (int*) a;
        std::cout << "Result: " << *x << std::endl;
        delete x;
    }
};

int* fun1(int* x){
    *x = *x + 3;
    return (int*) x;
}

int* fun2(int* x){
    *x = *x * 4;
    return (int*) x;
}

class LastStage: public Computable{
public:
    void compute(Data* d){
        int* in = (int*) d->getInput();
        *in = *in + 1;
        d->setOutput((void*) in);
    }
};

int main(int argc, char** argv){
    if(argc < 2){
        std::cerr << "Usage: " << argv[0] << " streamSize" << std::endl;
        return -1;
    }

    /* Create streams. */
    DemoInputStream inp(atoi(argv[1]));
    DemoOutputStream out;


    Computable* pipeTmp = createStandardPipeline<int, int, int, fun1, fun2>();
    Computable* lastStage = new LastStage();
    Pipeline* pipe = new Pipeline(pipeTmp, lastStage);
    nornir::Parameters p("parameters.xml");
    nornir::dataflow::Interpreter m(&p, pipe, &inp, &out);
    m.start();
    m.wait();
    delete pipe;
    delete lastStage;
}

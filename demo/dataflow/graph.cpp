/*
 * graph.cpp
 *
 * Created on: 30/04/2016
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
            int* s = new int();
            *s = rand();
            std::cout << "Generated " << *s << std::endl;
            std::cout << "Expected output " << ((*s + 1)*2  + (*s + 2) + 10 + 2) * ( (*s + 2) + 20) << std::endl;
            return (void*) s;
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
        std::cout << "Result: " << *((int*) a) << std::endl;
        delete (int*) a;
    }
};

Computable *A, *B, *C, *D, *E;

class ACode: public Computable{
public:
    void compute(Data* d){
        int* x = (int*) d->getInput();
        int* toB = new int(*x + 1);
        int* toC = new int(*x + 2);
        delete x;
        d->setOutput(toB, B);
        d->setOutput(toC, C);
    }
};

class BCode: public Computable{
    void compute(Data* d){
        int* x = (int*) d->getInput(A);
        (*x) = (*x) * 2;
        d->setOutput(x, D);
    }
};

class CCode: public Computable{
    void compute(Data* d){
        int* x = (int*) d->getInput(A);
        int* toD = new int((*x) + 10);
        int* toE = new int((*x) + 20);
        delete x;
        d->setOutput(toD, D);
        d->setOutput(toE, E);
    }
};

class DCode: public Computable{
    void compute(Data* d){
        int* fromB = (int*) d->getInput(B);
        int* fromC = (int*) d->getInput(C);
        *fromC = *fromC + *fromB + 2;
        delete fromB;
        d->setOutput(fromC, E);
    }
};

class ECode: public Computable{
    void compute(Data* d){
        int* fromC = (int*) d->getInput(C);
        int* fromD = (int*) d->getInput(D);
        *fromC = *fromC * *fromD;
        delete fromD;
        d->setOutput(fromC);
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

    /* Create instructions. */
    A = new ACode();
    B = new BCode();
    C = new CCode();
    D = new DCode();
    E = new ECode();

    /* Link instructions. */
    Mdfg graph;
    graph.link(A, B);
    graph.link(A, C);
    graph.link(B, D);
    graph.link(C, D);
    graph.link(C, E);
    graph.link(D, E);

    nornir::Parameters p("parameters.xml");
    nornir::dataflow::Interpreter m(&p, &graph, &inp, &out);
    m.start();
    m.wait();
}


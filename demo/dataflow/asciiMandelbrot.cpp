/*
 * asciiMandelbrot.cpp
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

class MandelInputStream: public nornir::dataflow::InputStream{
private:
    int x, y;
    bool eos;
    clock_t t;
public:
    inline MandelInputStream():x(-39), y(-39), eos(false),
                               t(clock()){;}

    inline void* next(){
        if(x == 39 && y == 38){
            eos = true;
        }
        if(x == 39){
            x =- 39;
            ++y;
            ArrayWrapper<float> *t = new ArrayWrapper<float>(2);
            t->set(0, -100000);
            return t;
        }
        ArrayWrapper<float> *t = new ArrayWrapper<float>(2);
        t->set(0, x/40.0f);
        t->set(1, y/40.0f);
        ++x;
        return t;
    }

    inline bool hasNext(){return !eos;}

};

class MandelOutputStream: public nornir::dataflow::OutputStream{
public:
    void put(void* a){
        float* dt = (float*) a;
        float x = *dt;
        if(x == 0)
            std::cout << "*";
        else if(x == -1)
            std::cout << std::endl;
        else
            std::cout << " ";
        delete dt;
    }
};

int BAILOUT = 16;
int MAX_ITERATIONS;

/**Function taken from http://www.timestretch.com/FractalBenchmark.html#67b4f5a3200c7b7e900c38ff21321741 **/
float* iterateTask(ArrayWrapper<float>* t){
    float f1 = t->get(1);
    float f2 = t->get(0);
    delete t;
    if(f2 == -100000){
        return new float(-1);
    }

    float cr = f1-0.5f;
    float ci = f2;
    float zi = 0.0f;
    float zr = 0.0f;
    int i = 0;
    while (true){
        i++;
        float temp = zr * zi;
        float zr2 = zr * zr;
        float zi2 = zi * zi;
        zr = zr2 - zi2 + cr;
        zi = temp + temp + ci;
        if (zi2 + zr2 > BAILOUT)
            return new float(i);

        if (i > MAX_ITERATIONS)
            return new float(0);

    }
}

int main(int argc, char** argv){
    if(argc < 2){
        std::cerr << "Usage:" << std::endl;
        std::cerr << argv[0] << " maxIterations" << std::endl;
        return -1;
    }else{
        MAX_ITERATIONS = atoi(argv[1]);
        MandelInputStream inp;
        MandelOutputStream out;
        nornir::dataflow::Farm* f = createStandardFarm<ArrayWrapper<float>, float, iterateTask>();
        nornir::Parameters p("parameters.xml");
        Interpreter m(&p, f, &inp, &out);
        m.start();
        m.wait();
        delete f;
        return 0;
    }
}


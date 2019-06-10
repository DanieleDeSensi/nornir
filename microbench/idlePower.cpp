#include <mammut/mammut.hpp>

#include <cassert>
#include <iostream>
#include <unistd.h>

#define PERIOD 10

using namespace mammut;
using namespace mammut::energy;
using namespace std;

int main(int argc, char** argv){
    Mammut m;
    Energy* energy = m.getInstanceEnergy();
    Joules j;

    /** Gets the energy counters (one per CPU). **/
    Counter* counter = energy->getCounter();
    if(!counter){
        cout << "Power counters not available on this machine." << endl;
        return -1;
    }

    counter->reset();
    sleep(PERIOD);
    j = counter->getJoules();
    cout << "<idlePower>" << j / (double) PERIOD << "</idlePower>" << endl;
}

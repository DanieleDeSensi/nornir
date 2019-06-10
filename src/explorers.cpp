/*
 * explorers.cpp
 *
 * Created on: 10/07/2015
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

#include "explorers.hpp"
#include "configuration.hpp"

#undef DEBUG
#undef DEBUGB

#ifdef DEBUG_EXPLORERS
#define DEBUG(x) do { std::cerr << "[Explorers] " << x << std::endl; } while (0)
#define DEBUGB(x) do {x} while (0)
#else
#define DEBUG(x)
#define DEBUGB(x)
#endif

namespace nornir{
    Explorer::Explorer(std::vector<bool> knobs,
                       std::vector<KnobsValues> additionalPoints):
        _knobs(knobs),
        _additionalPoints(additionalPoints){;}

    KnobsValues Explorer::getNextAdditionalPoint() const{
        KnobsValues r = _additionalPoints.front();
        if(!r.areRelative()){
            throw std::runtime_error("Additional calibration points must be relative.");
        }
        _additionalPoints.erase(_additionalPoints.begin());
        return r;
    }

    ExplorerRandom::ExplorerRandom(std::vector<bool> knobs, std::vector<KnobsValues> additionalPoints):
            Explorer(knobs, additionalPoints){
        srand(time(NULL));
    }

    KnobsValues ExplorerRandom::nextRelativeKnobsValues() const{
        KnobsValues r(KNOB_VALUE_RELATIVE);
        if(_additionalPoints.size()){
            r = getNextAdditionalPoint();
        }else{
            for(size_t i = 0; i < KNOB_NUM; i++){
                if(generate((KnobType) i)){
                    r[(KnobType)i] = rand() % 100;
                }else{
                    /**
                     * If we do not need to automatically find the value for this knob,
                     * then it has only 0 or 1 possible value. Accordingly, we can set
                     * it to any value.
                     */
                    r[(KnobType)i] = 0;
                }
            }
        }
        return r;
    }

    void ExplorerRandom::reset(){;}

    ExplorerLowDiscrepancy::ExplorerLowDiscrepancy(
                std::vector<bool> knobs,
                StrategyExploration explorationStrategy,
                std::vector<KnobsValues> additionalPoints):
            Explorer(knobs, additionalPoints), _explorationStrategy(explorationStrategy){
        uint d = 0;
        for(size_t i = 0; i < KNOB_NUM; i++){
            if(generate((KnobType) i)){
                ++d;
            }
        }
        DEBUG("[Explorator] Will generate low discrepancy points in " << d <<
              " dimensions");
        const gsl_qrng_type* generatorType;
        switch(_explorationStrategy){
            case STRATEGY_EXPLORATION_NIEDERREITER:{
                generatorType = gsl_qrng_niederreiter_2;
            }break;
            case STRATEGY_EXPLORATION_SOBOL:{
                generatorType = gsl_qrng_sobol;
            }break;
            case STRATEGY_EXPLORATION_HALTON:{
                generatorType = gsl_qrng_halton;
            }break;
            case STRATEGY_EXPLORATION_HALTON_REVERSE:{
                generatorType = gsl_qrng_reversehalton;
            }break;
            default:{
                throw std::runtime_error("ExplorerLowDiscrepancy: Unknown "
                                         "exploration strategy: " +
                                         _explorationStrategy);
            }break;
        }
        _generator = gsl_qrng_alloc(generatorType, d);
        _normalizedPoint = new double[d];
    }

    ExplorerLowDiscrepancy::~ExplorerLowDiscrepancy(){
        gsl_qrng_free(_generator);
        delete[] _normalizedPoint;
    }


    KnobsValues ExplorerLowDiscrepancy::nextRelativeKnobsValues() const{
        KnobsValues r(KNOB_VALUE_RELATIVE);
        if(_additionalPoints.size()){
            r = getNextAdditionalPoint();
        }else{
            gsl_qrng_get(_generator, _normalizedPoint);
            size_t nextCoordinate = 0;
            for(size_t i = 0; i < KNOB_NUM; i++){
                if(generate((KnobType) i)){
                    r[(KnobType)i] = _normalizedPoint[nextCoordinate]*100.0;
                    ++nextCoordinate;
                }else{
                    /**
                     * If we do not need to automatically find the value for this knob,
                     * then it has only 0 or 1 possible value. Accordingly, we can set
                     * it to any value.
                     */
                    r[(KnobType)i] = 0;
                }
            }
        }
        return r;
    }

    void ExplorerLowDiscrepancy::reset(){
        gsl_qrng_init(_generator);
    }

    ExplorerMultiple::ExplorerMultiple(std::vector<bool> knobs,
                                       Explorer* explorer,
                                       KnobType kt,
                                       size_t numValues):
                Explorer(knobs), _explorer(explorer), _kt(kt),
                _numValues(numValues), _nextValue(0), _lastkv(explorer->nextRelativeKnobsValues()){
        ;
    }

    ExplorerMultiple::~ExplorerMultiple(){
        ;
    }

    void ExplorerMultiple::reset(){
        _explorer->reset();
    }

    KnobsValues ExplorerMultiple::nextRelativeKnobsValues() const{
        if(_nextValue == _numValues){
            _nextValue = 0;
            _lastkv = _explorer->nextRelativeKnobsValues();
        }
        KnobsValues kv = _lastkv;
        kv[_kt] = _nextValue * (100.0 / (double) _numValues);
        ++_nextValue;
        return kv;
    }
}

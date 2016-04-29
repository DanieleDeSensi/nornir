/*
 * mdfi.hpp
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

#ifndef NORNIR_DF_MDFI_HPP_
#define NORNIR_DF_MDFI_HPP_

#include "skeleton/computable.hpp"
#include "tokens.hpp"

#include <cassert>
#include <iostream>
#include <vector>


namespace nornir{
namespace dataflow{

const uint sourceInput = std::numeric_limits<uint>::max();

#define MAX_INSTRUCTION_INPUTS 1

/**Macro data flow instruction.**/
class Mdfi{
private:
    std::map<Computable*, uint> sourcesMap; // Map between computable and its position in the sources vector
    std::map<Computable*, uint> destinationsMap;
    uint sources[MAX_INSTRUCTION_INPUTS];
    InputToken tInput[MAX_INSTRUCTION_INPUTS]; ///<Input tokens. Are numbered from 0 to \e dInput -1.
    TokenId dest[MAX_INSTRUCTION_INPUTS]; ///<Destinations of the output. Are numbered from 0 to \e dOutput -1.
    OutputToken tOutput[MAX_INSTRUCTION_INPUTS];
    Computable* comp; ///<\e Computable to compute.
    unsigned int dInput, ///<Input size.
        dOutput, ///<Output size.
        id; ///<Instruction's identifier.
    unsigned long int graphId; ///<Id of the graph to which the instruction belongs.
public:
    /**
     * Constructor of the instruction.
     * \param c A pointer to the \e Computable that the instruction have to compute.
     * \param i Instruction's identifier.
     */
    inline Mdfi(Computable* c, unsigned int i):
        comp(c),dInput(0),dOutput(0),id(i),graphId(0){
        ;
    }

    /**
     * Copy constructor.
     * \param ins Reference to instruction to copy.
     */
    inline Mdfi(const Mdfi& ins):
            sourcesMap(ins.sourcesMap), destinationsMap(ins.destinationsMap),
            comp(ins.comp), dInput(ins.dInput),
            dOutput(ins.dOutput), id(ins.id), graphId(ins.graphId){
        for(uint i = 0; i < dInput; i++){
            sources[i] = ins.sources[i];
        }
        for(uint i = 0; i < dOutput; i++){
            dest[i] = ins.dest[i];
        }
        assert(dInput <= MAX_INSTRUCTION_INPUTS);
    }

    /**
     * Destructor of the instruction.
     */
    inline ~Mdfi(){
        ;
    }

    /**
     * Computes the result of the instruction.
     */
    void compute();

    inline uint getNumOutTokens() const{
        return dOutput;
    }

    inline OutputToken* getOutToken(uint i){
        assert(i < dOutput);
        return &(tOutput[i]);
    }

    /**
     * Checks if the instruction is fireable.
     * \return \e true if the instruction is fireable, otherwise returns \e false.
     */
    inline bool isFireable() const{
        uint i = 0;
        while(i < dInput){
            if(!tInput[i].isPresent()) return false;
            i++;
        }
        return true;
    }

    /**
     * Sets an input task.
     * \param t The task to set.
     * \param i The position where to set the task.
     * \return \e false if \e i is higher than \e dInput, otherwise returns false.
     */
    inline bool setInput(void* t, uint i){
        if(i < dInput) {
            tInput[i].setTask(t);
            return true;
        }else
            return false;
    }

    inline void setSourceInStream(){
        sources[dInput] = sourceInput;
        ++dInput;
        assert(dInput <= MAX_INSTRUCTION_INPUTS);
    }

    inline void setSource(uint s, Computable* c = NULL){
        sources[dInput] = s;
        if(c){
            sourcesMap.emplace(c, dInput);
        }
        ++dInput;
        assert(dInput <= MAX_INSTRUCTION_INPUTS);
    }

    inline void setDestination(uint d, Computable* c = NULL){
        dest[dOutput] = TokenId(d, dOutput);
        if(c){
            sourcesMap.emplace(c, dOutput);
        }
        ++dOutput;
        assert(dOutput <= MAX_INSTRUCTION_INPUTS);
    }

    inline void setDestinationOutStream(){
        dest[dOutput].setOutputStream();
        ++dOutput;
        assert(dOutput <= MAX_INSTRUCTION_INPUTS);
    }

    /**
     * Returns the id of the instruction.
     * \return The id of the instruction.
     */
    inline uint getId(){
        return id;
    }

    /**
     * Sets the id of the instruction.
     * \param i The id of the instruction.
     */
    inline void setId(uint i){
        id = i;
    }

    /**
     * Returns the size of the input.
     * \return The size of the input.
     */
    inline unsigned int getInputSize(){
        return dInput;
    }

    /**
     * Returns the size of the output.
     * \return The size of the output.
     */
    inline unsigned int getOutputSize(){
        return dOutput;
    }

    /**
     * Returns the id of the graph.
     * \return The id of the graph.
     */
    inline unsigned long int getGid(){
        return graphId;
    }

    /**
     * Sets the id of the graph.
     * \param newd The new id of the graph.
     */
    inline void setGid(ulong newd){
        graphId = newd;
        for(uint i = 0; i < dOutput; i++)
            dest[i].setGraphId(newd);
    }

    /**
     * Returns a pointer to the \e Computable that the instruction have to compute.
     * \return A pointer to the \e Computable that the instruction have to compute.
     */
    inline Computable* getComputable(){
        return comp;
    }

    /**
     * Updates the destinations of the instructions.
     * \param v \e v[i] is the new identifier of the instruction that previously has \e i as identifier.
     */
    void updateDestinations(int* v);

    /**
     * Resets the instruction.
     * \param newId The new id of the graph wich belongs the instruction.
     */
    void reset(ulong newId);
};

}
}

#endif /* NORNIR_DF_MDFI_HPP_ */


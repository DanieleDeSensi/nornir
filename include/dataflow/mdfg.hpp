/*
 * mdfg.hpp
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

#ifndef NORNIR_DF_MDFG_HPP_
#define NORNIR_DF_MDFG_HPP_

#include "mdfi.hpp"

namespace nornir{
namespace dataflow{

/**
 * Macro data flow graph.
 * Accessory methods are provided to interact with and use the graph..
 *
 * Usage example:
 * \code
 *      using namespace nornir;
 *      using namespace nornir::dataflow;
 *
 *      Computable *A, *B, *C, *D, *E;
 *      A = new ACode(...);
 *      B = new BCode(...);
 *      C = new CCode(...);
 *      D = new DCode(...);
 *      E = new ECode(...);
 * \endcode
 *
 * To link them together:
 *
 * \code
 *      Mdfg graph;
 *      graph.link(A, B);
 *      graph.link(A, C);
 *      graph.link(B, D);
 *      graph.link(C, D);
 *      graph.link(C, E);
 *      graph.link(D, E);
 *  //TODO AGGIUNGERE IMG GRAFO.
 * \endcode
 */
class Mdfg{
private:
    /**Next free index**/
    uint _nextId;
    /**Graph's id.**/
    ulong _id;
    /**Instructions of the graph.**/
    std::vector<Mdfi> _instructions;
    ulong _firstId;
    ulong _lastId;
    bool _init;

    friend Mdfg* compile(Computable* c);

    inline void clearInputInstruction(){
        _instructions.at(_firstId).clearInput();
        _firstId = std::numeric_limits<ulong>::max();
    }

    inline void clearOutputInstruction(){
        _instructions.at(_lastId).clearOutput();
        _lastId = std::numeric_limits<ulong>::max();
    }
public:
    /**
     * Constructor of the graph.
     */
    Mdfg();

    /**
     * Builds a graph with only one instruction.
     * \param c A pointer to the \e Computable that the instruction have to
     *        compute.
     */
    explicit Mdfg(Computable* c);

    /**
     * Copy constructor.
     * \param g Reference to graph to copy
     * \param gid The new id of the graph's copy.
     */
    Mdfg(const Mdfg& g, ulong gid);

    /**
     * Graph's destructor.
     */
    ~Mdfg();


    inline ulong getId(){
        return _id;
    }


    /**
     * Returns the number of instructions.
     * \return Number of instructions.
     */
    inline uint getNumMdfi(){
        return _instructions.size();
    }

    /**
     * Creates an instruction of the graph.
     * \param c A pointer to the \e Computable that the instruction have to
     *        compute.
     * \return The id of the created instruction.
     */
    inline uint createMdfi(Computable* c){
        _instructions.emplace_back(c, _nextId);
        ++_nextId;
        return _instructions.size() - 1;
    }

    /**
     * Creates an instruction of the graph.
     * \param in The instruction to be included into the graph.
     * \return The id of the created instruction.
     */
    inline uint createMdfi(Mdfi* in){
        _instructions.emplace_back(*in);
        _instructions.back().setId(_nextId);
        ++_nextId;
        return _instructions.size() - 1;
    }

    /**
     * Creates an instruction of the graph starting from an
     * instruction of another graph.
     * \param g The other graph.
     * \param id The id of the mdfi to be copied.
     * \return The id of the created instruction.
     */
    inline uint createMdfi(Mdfg* g, uint id){
        _instructions.emplace_back(g->_instructions.at(id));
        _instructions.back().setId(_nextId);
        ++_nextId;
        return _instructions.size() - 1;
    }

    inline void addOffset(size_t id, size_t offset){
        _instructions.at(id).addOffset(offset);
    }

    /**
     * Returns the id of the first instruction of the graph.
     * \return The id of the created instruction.
     */
    inline Computable* getFirst(){
        if(_firstId != std::numeric_limits<ulong>::max()){
            return _instructions.at(_firstId).getComputable();
        }else{
            throw std::runtime_error("No Mdfi in the graph.");
        }
        return 0;
    }

    /**
     * Returns the id of the last instruction of the graph.
     * \return The id of the created instruction.
     */
    inline Computable* getLast(){
        if(_lastId != std::numeric_limits<ulong>::max()){
            return _instructions.at(_lastId).getComputable();
        }else{
            throw std::runtime_error("No Mdfi in the graph.");
        }
    }

    inline uint getFirstId(){
        if(_firstId != std::numeric_limits<ulong>::max()){
            return _firstId;
        }else{
            throw std::runtime_error("No Mdfi in the graph.");
        }
        return 0;
    }

    /**
     * Returns the id of the last instruction of the graph.
     * \return The id of the created instruction.
     */
    inline uint getLastId(){
        if(_lastId != std::numeric_limits<ulong>::max()){
            return _lastId;
        }else{
            throw std::runtime_error("No Mdfi in the graph.");
        }
    }

    inline void link(Computable* a, Computable* b){
        long idA = -1;
        long idB = -1;
        size_t i = 0;
        while((idA == -1 || idB == -1) && i < _instructions.size()){
            if(_instructions.at(i).getComputable() == a){
                idA = i;
            }
            if(_instructions.at(i).getComputable() == b){
                idB = i;
            }
            i++;
        }

        if(idA == -1){
            /* A not already present in the graph. */
            idA = createMdfi(a);
        }

        if(idB == -1){
            /* B not already present in the graph. */
            idB = createMdfi(b);
        }
        _instructions.at(idA).setDestination(idB, b);
        _instructions.at(idB).setSource(a);
    }

    /**
     * MUST be called before using the graph.
     */
    inline void init(){
        if(!_init){
            bool inFound = false, outFound = false;
            for(size_t i = 0; i < _instructions.size(); i++){
                if(!_instructions.at(i).getInputSize()){
                    if(!inFound){
                        inFound = true;
                        _instructions.at(i).setSourceInStream();
                        _firstId = i;
                    }else{
                        throw std::runtime_error("More than 1 instruction have "
                                                 "no input links.");
                    }
                }

                if(!_instructions.at(i).getOutputSize()){
                    if(!outFound){
                        outFound = true;
                        _instructions.at(i).setDestinationOutStream();
                        _lastId = i;
                    }else{
                        throw std::runtime_error("More than 1 instruction have "
                                                 "no output links.");
                    }
                }
            }
            _init = true;
            if(!inFound){
                throw std::runtime_error("Impossible to find the first instruction "
                                         "of the graph (Must have only 1 input).");
            }

            if(!outFound){
                throw std::runtime_error("Impossible to find the last instruction "
                                         "of the graph (Must have only 1 output).");
            }
        }
    }

    inline bool isFireable(uint id) const{
        return _instructions.at(id).isFireable();
    }

    /**
     * ATTENTION: SHOULD ONLY BE CALLED AFTER THE GRAPH HAS BEEN CREATED.
     *            If instructions are added/removed, this pointer will be
     *            invalidated.
     */
    inline Mdfi* getMdfi(uint id){
        return &(_instructions.at(id));
    }

    inline Computable* getComputable(uint id){
        return _instructions.at(id).getComputable();
    }

    /**
     * Resets all the instructions.
     * \param newId The new id of the graph.
     */
    void reset(ulong newId);

    inline void deinit(){
        clearInputInstruction();
        clearOutputInstruction();
    }

};
}
}

#endif /* MDFG_HPP_ */

/*
 * pipeline.hpp
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

#ifndef NORNIR_DF_PIPELINE_HPP_
#define NORNIR_DF_PIPELINE_HPP_

namespace nornir{
namespace dataflow{

/**
 * \class Pipeline
 * This class represents a generic 2-stage pipeline.
 * Multiple stage pipelines can be defined nesting
 * two stage pipelines.
 * The following example shows how to create a 3-stage pipeline.
 *
 * \code
 * Pipeline x(Stage1,Stage2);
 * Pipeline pipe(x,Stage3);
 * \endcode
 *
 * Where Stage1,Stage2 and Stage3 are other skeletons.
 *
 */
class Pipeline: public Computable{
private:
    Computable *s1,
        *s2;
    bool deleteStages;
public:
    /**
     * The constructor of the pipeline.
     * \param stage1 The skeleton used to implement the first stage of the pipeline.
     * \param stage2 The skeleton used to implement the second stage of the pipeline.
     * \param deleteStages If true, the destructor of the pipeline deletes the stages.
     */
    Pipeline(Computable *stage1, Computable *stage2, bool deleteStages=false):s1(stage1),s2(stage2),deleteStages(deleteStages){;}

    /**
     * Destructor of the pipeline.
     * If deleteStages is true, deletes the stages of the pipeline.
     */
    ~Pipeline(){
        if(deleteStages){
            delete s1;
            delete s2;
        }
    }

    /**
     * This method computes the result sequentially.
     * \param t Array of Task*.
     * \return Array of Task* as result. The array of Task* and its elements must be deallocated using delete[] and delete.
     *
     */
    void compute(Data* d){
        void* stage1Result;
        void* stage2Result;

        Data d1, d2;
        d1.setSource(d->getInput());
        d1.setDestination(&stage1Result);
        s1->compute(&d1);

        d2.setSource(stage1Result);
        d2.setDestination(&stage2Result);
        s2->compute(&d2);

        d->setOutput(stage2Result);
    }

    /**
     * Returns a pointer to the first stage of the pipeline.
     * \return A pointer to the second stage of the pipeline.
     */
    Computable* getFirstStage(){
        return s1;
    }

    /**
     * Returns a pointer to the second stage of the pipeline.
     * \return A pointer to the second stage of the pipeline.
     */
    Computable* getSecondStage(){
        return s2;
    }
};


template <typename In, typename Out, Out*(*fun)(In*)> class PipelineStage: public Computable{
public:
    void compute(Data* d){
        In* x = (In*) d->getInput();
        d->setOutput((void*) fun(x));
    }
};


/**
 * Creates a standard pipeline.
 *
 * \tparam T is the type of the input elements of the first stage.
 * \tparam U is the type of the output elements of the first stage (and also of the input elements of the second stage).
 * \tparam V is the type of the output elements of the second stage.
 * \tparam stage1 Is the function that the first stage have to compute.
 * \tparam stage2 Is the function that the second stage have to compute.
 * \return A pointer to a standard pipeline.
 */

template<typename T, typename U, typename V, U*(*stage1)(T*), V*(*stage2)(U*) > Pipeline* createStandardPipeline(){
    return new Pipeline(new PipelineStage<T,U,stage1>,new PipelineStage<U,V,stage2>,true);
}


}
}

#endif /* NORNIR_DF_PIPELINE_HPP_ */

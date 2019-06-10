/*
 * manager-ff.hpp
 *
 * Created on: 23/03/2015
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

/*!
 * @file manager-ff.hpp
 * @brief Implementation of adaptive fastflow patterns.
 */
#ifndef NORNIR_MANAGER_FF_HPP_
#define NORNIR_MANAGER_FF_HPP_
#include <nornir/ffincs.hpp>
#include <nornir/manager.hpp>

#include <nornir/external/fastflow/ff/config.hpp>
#include <nornir/external/fastflow/ff/pipeline.hpp>

namespace nornir{
    template<typename IN_t, typename OUT_t = IN_t>
struct nrnr_node_t: public AdaptiveNode {
    typedef IN_t  in_type;
    typedef OUT_t out_type;

    nrnr_node_t():
        GO_ON((OUT_t*)ff::FF_GO_ON),
        EOS((OUT_t*)ff::FF_EOS),
        EOSW((OUT_t*)ff::FF_EOSW),
        GO_OUT((OUT_t*)ff::FF_GO_OUT),
        EOS_NOFREEZE((OUT_t*) ff::FF_EOS_NOFREEZE),
        BLK((OUT_t*)ff::FF_BLK), NBLK((OUT_t*)ff::FF_NBLK) {
    }
    OUT_t * const GO_ON,  *const EOS, *const EOSW, *const GO_OUT, *const EOS_NOFREEZE, *const BLK, *const NBLK;
    virtual ~nrnr_node_t()  {}
    virtual OUT_t* svc(IN_t*)=0;
    inline  void *svc(void *task) {
        void* r = svc(reinterpret_cast<IN_t*>(task));
        if(!r){
            TERMINATE_APPLICATION;
        }else{
            return r;
        }
    };
};

template<typename TIN, typename TOUT, typename FUNC>
struct nrnr_node_F: public nrnr_node_t<TIN, TOUT> {
   nrnr_node_F(FUNC f):F(f) {}
   TOUT* svc(TIN* task) {
       TOUT* r = F(task, this);
       if(!r){
           TERMINATE_APPLICATION_TYPED(TOUT);
       }else{
           return r;
       }
   }
   FUNC F;
};

/*!
 * \class ManagerFastFlow
 * \brief This class manages the adaptivity in applications written
 * with FastFlow programming framework.
 *
 * This class manages the adaptivity in applications written
 * with FastFlow programming framework.
 */
class ManagerFastFlow: public Manager{
    template <typename I, typename O> friend class FarmBase;
public:
    /**
     * Creates a farm adaptivity manager.
     * @param farm The farm to be managed.
     * @param nornirParameters The parameters to be used for
     * adaptivity decisions.
     */
    ManagerFastFlow(ff::ff_farm<>* farm, Parameters nornirParameters);

    /**
     * Destroyes this adaptivity manager.
     */
    ~ManagerFastFlow();
private:
    // The managed farm.
    ff::ff_farm<>* _farm;

    // The emitter (if present).
    AdaptiveNode* _emitter;

    // The collector (if present).
    AdaptiveNode* _collector;

    // The vector of active workers.
    std::vector<AdaptiveNode*> _activeWorkers;

    void waitForStart();
    MonitoredSample getSample();
    void postConfigurationManagement();
    void terminationManagement();
    ulong getExecutionTime();
    void shrinkPause();
    void stretchPause();
};

/*!
 * \class ManagerFastFlowPipeline
 * \brief This class manages the adaptivity in applications written
 * with pipeline pattern in the FastFlow programming framework.
 *
 * This class manages the adaptivity in applications written
 * with pipeline pattern in the FastFlow programming framework.
 */
class ManagerFastFlowPipeline: public Manager{
    template <typename I, typename O> friend class FarmBase;
public:
    /**
     * Creates a pipeline adaptivity manager.
     * @param pipe The pipeline to be managed.
     * @param nornirParameters The parameters to be used for
     * adaptivity decisions.
     */
    ManagerFastFlowPipeline(ff::ff_pipeline* pipe, std::vector<bool> farmsFlags, Parameters nornirParameters);

    /**
     * Destroyes this adaptivity manager.
     */
    ~ManagerFastFlowPipeline();
private:
    // The managed farm.
    ff::ff_pipeline* _pipe;
    std::vector<bool> _farmsFlags;
    std::vector<KnobVirtualCoresFarm*> _farmsKnobs;
    std::vector<AdaptiveNode*> _activeWorkers;
    std::vector<std::vector<AdaptiveNode*>> _allWorkers;
    std::vector<std::vector<double>> _allowedValues;

    
    void waitForStart();
    MonitoredSample getSample();
    void postConfigurationManagement();
    ulong getExecutionTime();
    void shrinkPause();
    void stretchPause();
};
}

#endif // NORNIR_MANAGER_FF_HPP_
/*
 * interface.hpp
 *
 * Created on: 27/02/2016
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

#ifndef NORNIR_INTERFACE_HPP_
#define NORNIR_INTERFACE_HPP_

#include "manager.hpp"

#include <cstddef>
#include <queue>

namespace nornir{

//TODO Gestire svc_end svc_init etc...

template <typename I, typename O> class FarmBase;

// For internal use
class OrderedTask{
private:
    std::pair<unsigned long long, void*> _pair;
public:
    OrderedTask(unsigned long long id, void* task):_pair(id, task){
        ;
    }

    unsigned long long getId() const{
        return _pair.first;
    }

    void* getTask() const{
        return _pair.second;
    }

    void setTask(void* task){
        _pair.second = task;
    }
};

class CompareOTasks{
public:
    bool operator()(const OrderedTask& a, const OrderedTask& b){
        return a.getId() > b.getId();
    }
};

/**
 * @class SchedulerBase
 * Common base class for scheduler (e.g. scheduler, accelerator scheduler, etc...).
 */
template <typename O> class SchedulerBase: public AdaptiveNode{
    template <typename T, typename V> friend class FarmBase;
private:
    ff::ff_loadbalancer* _lb;
    bool _ondemand;
    bool _preserveOrdering;
    unsigned long long _nextTaskId;

    void setLb(ff::ff_loadbalancer* lb){
        _lb = lb;
    }

    void setOndemand(){
        _ondemand = true;
    }

    void preserveOrdering(){
        _preserveOrdering = true;
    }
protected:
    void* transformTaskForOrdering(O* task){
        if(!task){
            return NULL;
        }
        if(_preserveOrdering){
            OrderedTask* ot = new OrderedTask(_nextTaskId, (void*) task);
            ++_nextTaskId;
            return (void*) ot;
        }else{
            return (void*) task;
        }
    }

    size_t getCurrentNumWorkers() const{
        return _lb->getnworkers();
    }

public:
    SchedulerBase():_lb(NULL), _ondemand(false), _preserveOrdering(false), _nextTaskId(0){;}

    /**
     * Sends a task to one of the workers.
     * @param task The task to be sent.
     */
    void send(O* task) CX11_KEYWORD(final){
        void* realTask = transformTaskForOrdering(task);
        while(!ff_send_out(realTask)){;}
    }

    /**
     * Sends a task to one of the workers (non blocking).
     * @param task The task to be sent.
     * @return true if the task has been sent, false otherwise.
     */
    bool sendNonBlocking(O* task) CX11_KEYWORD(final){
        void* realTask = transformTaskForOrdering(task);
        return ff_send_out(realTask);
    }

    /**
     * Sends a task to a specific worker.
     * @param task The task to be sent.
     * @param id The identifier of the worker. Starts from 0.
     */
    void sendTo(O* task, uint id) CX11_KEYWORD(final){
        if(id >= _lb->getnworkers()){
            throw new std::runtime_error("FATAL: Trying to send to a non running worker."
                                         "Please ensure that your application is "
                                         "notified when a rethreading occur.");
        }
        void* realTask = transformTaskForOrdering(task);
        while(!_lb->ff_send_out_to(realTask, id)){;}
    }

    /**
     * Sends a task to a specific worker (non blocking).
     * @param task The task to be sent.
     * @param id The identifier of the worker. Starts from 0.
     * @return true if the task has been sent, false otherwise.
     */
    bool sendToNonBlocking(O* task, uint id) CX11_KEYWORD(final){
        if(id >= _lb->getnworkers()){
            throw new std::runtime_error("FATAL: Trying to send to a non running worker."
                                         "Please ensure that your application is "
                                         "notified when a rethreading occur.");
        }
        void* realTask = transformTaskForOrdering(task);
        return _lb->ff_send_out_to(realTask, id);
    }

    /**
     * Sends a task to all the running workers.
     * @param task The task to be sent.
     */
    void broadcast(O* task) CX11_KEYWORD(final){
        // We can't rely on FastFlow broadcast
        // otherwise our ordering implementation
        // would not work.
        disableRethreading(); // Avoid doing rethreading during broadcast
        size_t nworkers = _lb->getnworkers();
        for(size_t i = 0; i < nworkers; i++){
            sendTo(task, i);
        }
        enableRethreading();
    }

    /**
     * Returns a special element which will be
     * ignored by the runtime.
     * @return A special element which will be
     * ignored by the runtime.
     **/
    O* nothing(){
        return (O*) GO_ON;
    }

    /**
     * Returns a special element which will be
     * considered as the last of the stream.
     * @return A special element which will be
     * considered as the last of the stream.
     **/
    O* lastElement(){
        return (O*) NULL;
    }
};

/**
 * @class Scheduler
 * @brief The scheduler of the farm (accelerator case).
 * @tparam O The type of the tasks produced by the scheduler.
 */
template <typename I, typename O = std::nullptr_t> class Scheduler: public SchedulerBase<O>{
private:
    // We force it to be private to avoid misuse by the user.
    using SchedulerBase<O>::transformTaskForOrdering;

    void* svc(void* task) CX11_KEYWORD(final){
        O* r = schedule((I*)(task));
        void* outTask = transformTaskForOrdering(r);
        if(outTask){
           return outTask;
        }else{
            AdaptiveNode::terminate();
            return (void*) ff::FF_EOS;
        }
    }
public:
    virtual ~Scheduler(){;}

    virtual O* schedule(I* task) = 0;
};


/**
 * @class Scheduler
 * @brief The scheduler of the farm.
 * @tparam O The type of the tasks produced by the scheduler.
 */
/** ATTENTION: Semantic of targs inverted!!! This is a dirty trick
 *             but we can avoid to have different class names according
 *             to the presence of an input channel.
 *             This a partial template specialisation. From the standard:
 *               The members of the class template partial specialization are
 *               unrelated to the members of the primary template. Class
 *               template partial specialization members that are used in a way
 *               that requires a definition shall be defined; the definitions of
 *               members of the primary template are never used as definitions
 *               for members of a class template partial specialization. [...]
 */
template <typename O>
class Scheduler<O, std::nullptr_t>: public SchedulerBase<O>{
private:
    // We force it to be private to avoid misuse by the user.
    using SchedulerBase<O>::transformTaskForOrdering;

    void* svc(void* task) CX11_KEYWORD(final){
        O* r = schedule();
        void* outTask = transformTaskForOrdering(r);
        if(outTask){
           return outTask;
        }else{
            AdaptiveNode::terminate();
            return (void*) ff::FF_EOS;
        }
    }
public:
    virtual ~Scheduler(){;}

    virtual O* schedule() = 0;
};

template <typename I, typename O> class SchedulerDummy: public Scheduler<I, O>{
public:
    O* schedule(I* task){
        return (O*) task;
    }
};

//! @cond
/**
 * @class WorkerBase
 * Common base class for workers (e.g. worker, worker with no output, etc...).
 */
template <typename I, typename O>
class WorkerBase: public AdaptiveNode{
    template <typename T, typename V> friend class FarmBase;
protected:
    bool _ordering;

    I* getComputeInput(void* t){
        if(_ordering){
            void* tvoid = reinterpret_cast<OrderedTask*>(t)->getTask();
            return reinterpret_cast<I*>(tvoid);
        }else{
            return reinterpret_cast<I*>(t);
        }
    }

    /**
     * Returns a special element which will be
     * ignored by the runtime.
     * @return A special element which will be
     * ignored by the runtime.
     **/
    O* nothing(){
        return (O*) GO_ON;
    }
private:
    void preserveOrdering(){
        _ordering = true;
    }

    using AdaptiveNode::enableRethreading; // Can only be used on schedueler
    using AdaptiveNode::disableRethreading; // Can only be used on schedueler
public:
    WorkerBase():_ordering(false){;}

    virtual ~WorkerBase(){;}

    /**
     * Returns the identifier of this worker.
     * @return The identifier of this worker.
     */
    uint getId() const{
        return get_my_id();
    }
};
//! @endcond

/**
 * @class Worker
 * @brief The worker of the farm.
 * @tparam I The type of the tasks received from the scheduler.
 * @tparam O The type of the tasks produced from the worker and sent to the
 *           gatherer. If not present, it means that the worker does not produce
 *           data.
 */
template <typename I, typename O = std::nullptr_t> class Worker: public WorkerBase<I, O>{
private:
    // We force it to be private to avoid misuse by the user.
    using WorkerBase<I, O>::getComputeInput;
    using WorkerBase<I, O>::_ordering;

    void* svc(void* t) CX11_KEYWORD(final){
        I* computeInput = getComputeInput(t);
        void* computeOutput = (void*) compute(computeInput);
        if(_ordering){
            OrderedTask* ot = reinterpret_cast<OrderedTask*>(t);
            ot->setTask((void*) computeOutput);
            computeOutput = (void*) ot;
        }
        return computeOutput;
    }
public:
    virtual ~Worker(){;}

    /**
     * Computes a function over a received task.
     * @param task The task received from the scheduler.
     * @return The processed task. It will be sent to the gatherer.
     */
    virtual O* compute(I* task) = 0;
};


//! @cond
/**
 * Template specialisation for the worker with no output.
 */
template <typename I>
class Worker<I, std::nullptr_t>: public WorkerBase<I, std::nullptr_t>{
private:
    // We force it to be private to avoid misuse by the user.
    using WorkerBase<I, std::nullptr_t>::getComputeInput;

    void* svc(void* t) CX11_KEYWORD(final){
        I* computeInput = getComputeInput(t);
        compute(computeInput);
        return (void*) GO_ON;
    }
public:
    virtual ~Worker(){;}

    /**
     * Computes a function over a received task.
     * No output will be produced.
     * @param task The task received from the scheduler.
     */
    virtual void compute(I* task) = 0;
};
//! @endcond


//! @cond
/**
 * @class GathererBase
 * @tparam I The type of data received from the workers.
 * Common base class for gatherers (e.g. gatherer, accelerator gatherer, etc...).
 */
template <typename I>
class GathererBase: public AdaptiveNode{
    template <typename T, typename V> friend class FarmBase;
private:
    bool _ordering;
    unsigned long long _nextTaskId;
    std::priority_queue<OrderedTask, std::vector<OrderedTask>, CompareOTasks> _priorityQueue;

    void preserveOrdering(){
        _ordering = true;
    }

    using AdaptiveNode::enableRethreading; // Can only be used on schedueler
    using AdaptiveNode::disableRethreading; // Can only be used on schedueler
protected:
    // Gets the gather inputs.
    // The vector can contain more than one
    // element only if ordering is required.
    void getGatherInputs(void* t, std::vector<I*>& toReturn){
        if(_ordering){
            OrderedTask* ot = reinterpret_cast<OrderedTask*>(t);
            if(ot->getId() == _nextTaskId){
                // Lucky. This is the correct task.
                toReturn.push_back(reinterpret_cast<I*>(ot->getTask()));
                ++_nextTaskId;
            }else{
                // Put in queue
                _priorityQueue.push(*ot);
            }
            while(!_priorityQueue.empty() && _priorityQueue.top().getId() == _nextTaskId){
                toReturn.push_back(reinterpret_cast<I*>(_priorityQueue.top().getTask()));
                _priorityQueue.pop();
                ++_nextTaskId;
            }
            delete ot;
        }else{
            toReturn.push_back(reinterpret_cast<I*>(t));
        }
    }
public:
    GathererBase():_ordering(false), _nextTaskId(0){;}

    //TODO: Receivefrom?
};
//! @endcond


/**
 * @class Gatherer
 * @brief The gatherer of the farm.
 * @tparam I The type of the tasks received from the workers.
 */
template <typename I, typename O = std::nullptr_t> class Gatherer: public GathererBase<I>{
private:
    using GathererBase<I>::getGatherInputs;

    void* svc(void* t) CX11_KEYWORD(final){
        std::vector<I*> realTasks;
        getGatherInputs(t, realTasks);
        for(I* task : realTasks){
            this->ff_send_out((void*) task);
        }
        return GO_ON;
    }
protected:
    /**
     * Returns a special element which will be
     * ignored by the runtime.
     * @return A special element which will be
     * ignored by the runtime.
     **/
    O* nothing(){
        return (O*) GO_ON;
    }
public:
    Gatherer(){;}

    virtual ~Gatherer(){;}

    /**
     * Computes a function over a task received from a worker.
     * @param task The task received from a worker.
     */
    virtual O* gather(I* task) = 0;
};

/**
 * Template specialisation for the worker with no output.
 */
template <typename I>
class Gatherer<I, std::nullptr_t>: public GathererBase<I>{
private:
    using GathererBase<I>::getGatherInputs;

    void* svc(void* t) CX11_KEYWORD(final){
        std::vector<I*> realTasks;
        getGatherInputs(t, realTasks);
        for(I* task : realTasks){
            gather(task);
        }
        return GO_ON;
    }
public:
    virtual ~Gatherer(){;}

    /**
     * Computes a function over a task received from a worker.
     * @param task The task received from a worker.
     */
    virtual void gather(I* task) = 0;
};

template <typename I> class GathererDummy: public Gatherer<I, I>{
public:
    I* gather(I* task){
        return task;
    }
};

/**
 * @class FarmBase
 * @brief The nornir farm.
 * @tparam I The type of the tasks produced by the scheduler and received
 *           by the workers.
 * @tparam O The type of the tasks produced by the workers and received by the
 *           gatherer
 */
template <typename I, typename O> class FarmBase: public mammut::utils::NonCopyable{
protected:
    ff::ff_farm<>* _farm;
    std::vector<ff::ff_node*> _rawWorkers;
    bool _schedulerHasInput;
private:
    bool _nodesCreated;
    const Parameters* _params;
    bool _paramsCreated;
    ManagerFastFlow* _manager;
    SchedulerBase<I>* _scheduler;
    GathererBase<O>* _gatherer;
    bool _feedback;
protected:
    std::vector<WorkerBase<I, O>* > _workers;

    template <class S, class W>
    void init(size_t numWorkers){
        _nodesCreated = true;
        setScheduler(new S());
        for(size_t i = 0; i < numWorkers; i++){
            setWorker(new W());
        }
    }

protected:
    /**
     * Adds the scheduler to the farm.
     * @param s The scheduler of the farm.
     */
    void setScheduler(SchedulerBase<I>* s){
        if(!s){
            throw std::runtime_error("setScheduler: Scheduler must be != NULL.");
        }
        if(!_farm){
            createFarm();
        }
        _scheduler = s;
        _farm->add_emitter(s);
    }

    /**
     * Adds a worker to the farm.
     * @param w A worker of the farm.
     */
    void setWorker(WorkerBase<I, O>* w){
        if(!w){
            throw std::runtime_error("setWorker: Worker must be != NULL.");
        }
        if(!_farm){
            createFarm();
        }
        _rawWorkers.push_back(dynamic_cast<ff::ff_node*>(w));
        _workers.push_back(w);
    }

    /**
     * Adds the gatherer to the farm.
     * @param g The gatherer of the farm.
     */
    void setGatherer(GathererBase<O>* g){
        if(!g){
            throw std::runtime_error("setGatherer: Gatherer must be != NULL.");
        }
        if(!_farm){
            createFarm();
        }
        _gatherer = g;
        _farm->add_collector(g);
    }

    virtual void createFarm(){
        _farm = new ff::ff_farm<>();
    }

    virtual void preStart(){
        ;
    }
public:
    /**
     * The constructor of the farm.
     * @param parameters The configuration parameters.
     */
    // TODO Parameters pointer or Copy?
    explicit FarmBase(const Parameters* parameters):_farm(NULL),
                                       _params(NULL), _manager(NULL),
                                       _scheduler(NULL), _gatherer(NULL){
        _nodesCreated = false;
        _params = parameters;
        _paramsCreated = false;
        _feedback = false;
        _schedulerHasInput = false;
    }

    /**
     * The constructor of the farm.
     * @param paramFileName The filename of the XML file containing the
     *        configuration parameters.
     */
    explicit FarmBase(const std::string& paramFileName):_farm(NULL),
                                          _params(NULL), _manager(NULL),
                                          _scheduler(NULL), _gatherer(NULL){
        _nodesCreated = false;
        _params = new Parameters(paramFileName);
        _paramsCreated = true;
        _feedback = false;
        _schedulerHasInput = false;
    }

    /**
     * The destructor of the farm.
     */
    virtual ~FarmBase(){
        if(_paramsCreated){
            delete _params;
        }

        if(_farm){
            delete _farm;
        }
    }

    /**
     * Starts the farm.
     * @tparam S The type of the scheduler. The scheduler must have a
     *           constructor with no parameters.
     * @tparam W The type of the workers. The worker must have a constructor with
     *           no parameters.
     * @param numWorkers The maximum number of workers to be used.
     */
    template <class S, class W>
    void start(size_t numWorkers){
        // TODO Metter un valore di default per i workers (e.g. 0). IN tal caso girarlo con un numero di workers uguali al numero di cores.
        init<S, W>(numWorkers);
        start();
    }

    /**
     * Starts the farm.
     * @tparam S The type of the scheduler. The scheduler must have a
     *           constructor with no parameters.
     * @tparam W The type of the workers. The worker must have a constructor with
     *           no parameters.
     * @tparam G The type of the gatherer. The gatherer must have a constructor
     *           with no parameters.
     * @param numWorkers The number of workers.
     */
    template <class S, class W, class G>
    void start(size_t numWorkers){
        init<S, W>(numWorkers);
        setGatherer(new G());
        start();
    }

    /**
     * Starts the farm.
     */
    void start(){
        if(!_farm){
            createFarm();
        }
        _farm->add_workers(_rawWorkers);
        if(_farm->getEmitter()){
            (dynamic_cast<SchedulerBase<I>*>(_farm->getEmitter()))->setLb(_farm->getlb());
        }
        if(_feedback){
            if(_schedulerHasInput){
                if(!_gatherer){
                    // We need to explicitely remove it otherwise it won't work.
                    _farm->remove_collector();
                }
                _farm->wrap_around();
            }else{
                throw std::runtime_error("If you want to use feedback, you need to specify as scheduler "
                                         "extending from nornir::Scheduler<I, O>, i.e. a scheduler with "
                                         "some input data, corresponding to the data received back on the feedback "
                                         "channel from workers or from gatherer.");
            }
        }else if(_schedulerHasInput){
            throw std::runtime_error("Unless you set a feedback communication, there is no need to have a scheduler "
                                     "accepting some input. Accordingly, if you need a feedback call "
                                     "the setFeedback() call, otherwise let your scheduler extend the "
                                     "nornir::Scheduler<I> class instead of the nornir::Scheduler<I, O> class.");
        }
        preStart();
        _manager = new ManagerFastFlow(_farm, *_params);
        _manager->start();
    }

    /**
     * Waits for the farm end.
     */
    void wait(){
        if(_manager){
            _manager->join();
        }
        if(_nodesCreated){
            if(_farm->getEmitter()){
                delete _farm->getEmitter();
            }
            for(size_t i = 0; i < _rawWorkers.size(); i++){
                delete _rawWorkers.at(i);
            }
            if(_farm->getCollector()){
                delete _farm->getCollector();
            }
        }
        if(_manager){
            delete _manager;
        }
    }

    /**
     * Returns the number of workers currently present in the farm.
     */
    size_t getCurrentNumWorkers() const{
        return _scheduler->getCurrentNumWorkers();
    }

    void stats(std::ostream& o) const{
        _farm->ffStats(o);
    }

    /**
     * Sets on demand scheduling.
     * If you use on demand scheduling, DON't use
     * sendTo(...) call if you have a user
     * defined scheduler.
     * This function must be called after the scheduler
     * has been set.
     **/
    void setOndemandScheduling(){
        _farm->set_scheduling_ondemand();
        if(_scheduler){
            // This is done to allow the scheduler to check
            // that sendTo is not used.
            _scheduler->setOndemand();
        }
    }

    /**
     * By default, the farm does not preserve the order of the
     * input elements. I.e. corresponding output elements may
     * be produced in a different order.
     * By calling this function it is possible to force
     * the farm to preserve the order of the elements.
     * This function must be called after the scheduler
     * and the gatherer have been set.
     **/
    void preserveOrdering(){
        if(!_gatherer){
            setGatherer(new GathererDummy<O>());
        }
        _scheduler->preserveOrdering();
        for(auto w : _workers){
            w->preserveOrdering();
        }
        _gatherer->preserveOrdering();
    }

    /**
     * Connects the gatherer to the scheduler.
     * If a gatherer is not present, connects each
     * worker to the scheduler.
     * ATTENTION: You must provide a scheduler that
     * extends nornir::Scheduler<O, I>, since the scheduler
     * needs to have an input (i.e. the output of the
     * gatherer (or the output of workers if gatherer not present)).
     **/
    void setFeedback(){
#ifdef BLOCKING_MODE
        throw std::runtime_error("setFeedback cannot be used with FastFlow BLOCKING_MODE macro defined.");
        // TODO Set nonblocking at runtime
#endif
        _feedback = true;
    }
};

/**
 * @class Farm
 * @brief The nornir farm.
 * @tparam I The type of the tasks produced by the scheduler and received
 *           by the workers.
 * @tparam O The type of the tasks produced by the workers and received by the
 *           gatherer (or back by the scheduler if feedback is used).
 */
template <typename I, typename O = std::nullptr_t> class Farm: public FarmBase<I, O>{
public:
    /**
     * The constructor of the farm.
     * @param parameter The configuration parameters.
     */
    explicit Farm(const Parameters* parameters):FarmBase<I,O>(parameters){;}

    /**
     * The constructor of the farm.
     * @param paramFileName The filename of the XML file containing the
     *        configuration parameters.
     */
    explicit Farm(const std::string& paramFileName):FarmBase<I,O>(paramFileName){;}

    /**
     * Adds the scheduler to the farm.
     * @param s The scheduler of the farm.
     */
    void addScheduler(Scheduler<I>* s){
        FarmBase<I,O>::setScheduler(s);
    }

    /**
     * Adds the scheduler to the farm.
     * @param s The scheduler of the farm.
     */
    virtual void addScheduler(Scheduler<O, I>* s){
        FarmBase<I,O>::setScheduler(s);
        FarmBase<I,O>::_schedulerHasInput = true;
    }

    /**
     * Adds a worker to the farm.
     * @param w A worker of the farm.
     */
    void addWorker(Worker<I, O>* w){
        FarmBase<I,O>::setWorker(w);
    }

    /**
     * Adds the gatherer to the farm.
     * @param g The gatherer of the farm.
     */
    void addGatherer(Gatherer<O>* g){
        FarmBase<I,O>::setGatherer(g);
    }
};

/**
 * @class FarmAcceleratorBase
 * @brief The base class for farm accelerators.
 * @tparam S The type of data sent from the application to the scheduler.
 * @tparam I The type of data sent by the scheduler to the workers.
 * @tparam O The type of data sent by the workers to the gatherer.
 * @tparam G The type of data sent by the gatherer to the application.
 */
template <typename S, typename I, typename O, typename G>
class FarmAcceleratorBase: public FarmBase<I, O>{
private:
    SchedulerDummy<S, I>* _schedulerDummy;
protected:
    void createFarm(){
        FarmBase<I, O>::_farm = new ff::ff_farm<>(true);
    }

    void preStart(){
        // Adds default gatherer and scheduler.
        if(!FarmBase<I, O>::_farm->getEmitter()){
            _schedulerDummy = new SchedulerDummy<S, I>();
            FarmBase<I, O>::setScheduler(_schedulerDummy);
        }
    }
public:
    /**
     * The constructor of the farm.
     * @param parameters The configuration parameters.
     */
    explicit FarmAcceleratorBase(const Parameters* parameters):
            FarmBase<I,O>::FarmBase(parameters),
            _schedulerDummy(NULL){;}

    /**
     * The constructor of the farm.
     * @param paramFileName The filename of the XML file containing the
     *        configuration parameters.
     */
    explicit FarmAcceleratorBase(const std::string& paramFileName):
             FarmBase<I,O>::FarmBase(paramFileName),
             _schedulerDummy(NULL){;}

    /**
     * Denstructor of the accelerator.
     */
    ~FarmAcceleratorBase(){
        if(_schedulerDummy){
            delete _schedulerDummy;
        }
    }

    /**
     * Starts the accelerator.
     * @tparam W The type of the workers. The worker must have a constructor with
     *           no parameters.
     * @param numWorkers The maximum number of workers to be used.
     */
#if 0
    template <class W>
    void start(size_t numWorkers){
        FarmBase<I, O>::start<SchedulerDummy<S, I>, W>(numWorkers); //TODO
    }
#endif

    /**
     * Adds the scheduler to the farm.
     * @param s The scheduler of the farm.
     */
    void addScheduler(Scheduler<S, I>* s){
        FarmBase<I,O>::setScheduler(s);
    }

    /**
     * Adds a worker to the farm.
     * @param w A worker of the farm.
     */
    void addWorker(Worker<I, O>* w){
        FarmBase<I,O>::setWorker(w);
    }

    /**
     * Offloads a task to the accelerator.
     * Don't use NULL since it is reserved for internal use.
     * This call is blocking and only returns when the task has been sent.
     * @param task The task to be offloaded.
     */
    void offload(S* task){
        void* realTask = (void*) task;
        while(!FarmBase<I,O>::_farm->offload(realTask?realTask:EOS));
        if(realTask == EOS){
            ((Scheduler<S, I>*) FarmBase<I, O>::_farm->getEmitter())->terminate();
        }
    }

    /**
     * Offloads a task to the accelerator.
     * Don't use NULL since it is reserved for internal use.
     * This call is non blocking.
     * @param task The task to be offloaded.
     * @return True if the task has been sent, false otherwise.
     */
    bool offloadNonBlocking(S* task){
        void* realTask = (void*) task;
        bool res = FarmBase<I,O>::_farm->offload(realTask?realTask:EOS, 1);
        if(res && realTask == EOS){
            ((Scheduler<S, I>*) FarmBase<I, O>::_farm->getEmitter())->terminate();
        }
        return res;
    }

    /**
     * Shutdown the farm. You need to wait for its termination with the
     * wait() call.
     * This call blocks until it is possible to send a shutdown request.
     */
    void shutdown(){
        offload(static_cast<S*>(EOS));
    }

    /**
     * Shutdown the farm. You need to wait for its termination with the
     * wait() call.
     * @return True if the shutdown request has been sent, false otherwise.
     */
    bool shutdownNonBlocking(){
        return offloadNonBlocking(EOS);
    }
};

/**
 * @class FarmAccelerator
 * @brief Accelerator with collector and output channel.
 *
 * @tparam S The type of data sent from the application to the scheduler.
 * @tparam I The type of data sent by the scheduler to the workers.
 * @tparam O The type of data sent by the workers to the gatherer.
 * @tparam G The type of data sent by the gatherer to the application.
 *
 * - 1 tparams specified: Accelerator with no collector. Same type between
 *     application and scheduler and between scheduler and workers.
 * - 2 tparams specified: Accelerator with no collector. One type between
 *     application and scheduler and the other between scheduler and workers.
 *     You MUST explicitly provide a Scheduler.
 * - 3 tparams specified: Accelerator with collector. One type between application
 *     and scheduler, one type between scheduler and workers, and one type between
 *     workers and collector.
 * - 4 tparams specified: Accelerator with collector. One type for each channel(s).
 *
 */
template <typename S, typename I = S, typename O = std::nullptr_t, typename G = std::nullptr_t>
class FarmAccelerator: public FarmAcceleratorBase<S, I, O, G>{
private:
    GathererDummy<O>* _gathererDummy;

    bool isManagementTask(void* task){
        return task == EOSW   ||
               task == GO_OUT ||
               task == BLK    ||
               task == NBLK;
    }
protected:
    void preStart(){
        FarmAcceleratorBase<S, I, O, G>::preStart();
        if(!FarmBase<I, O>::_farm->getCollector()){
            _gathererDummy = new GathererDummy<O>();
            FarmBase<I, O>::setGatherer(_gathererDummy);
        }
    }
public:
    /**
     * The constructor of the farm.
     * @param parameters The configuration parameters.
     */
    explicit FarmAccelerator(const Parameters* parameters):
            FarmAcceleratorBase<S, I, O, G>(parameters),
            _gathererDummy(NULL){;}

    /**
     * The constructor of the farm.
     * @param paramFileName The filename of the XML file containing the
     *        configuration parameters.
     */
    explicit FarmAccelerator(const std::string& paramFileName):
            FarmAcceleratorBase<S, I, O, G>(paramFileName),
            _gathererDummy(NULL){;}

    /**
     * Denstructor of the accelerator.
     */
    ~FarmAccelerator(){
        if(_gathererDummy){
            delete _gathererDummy;
        }
    }


    /**
     * Adds the gatherer to the farm.
     * @param g The gatherer of the farm.
     */
    void addGatherer(Gatherer<O, G>* g){
        FarmBase<I, O>::setGatherer(g);
    }

    /**
     * Gets a result from the accelerator.
     * This call is blocking and only returns when the task has been received.
     * @return The retrieved result. If NULL, then this is the last produced
     *         result.
     */
    G* getResult(){
        G* r = NULL;
        G** pr = &r;
        while((!FarmBase<I,O>::_farm->load_result((void**)pr) && *pr != EOS) ||
              isManagementTask((void*)*pr)){
            ;
        }

        if(*pr == EOS){
            *pr = NULL;
        }
        return *pr;
    }

    /**
     * Gets a result from the accelerator.
     * This call is non blocking.
     * @param validResult Will be set to true if a result has been retrieved,
     *        to false otherwise.
     * @return The retrieved result. If NULL, then this is the last produced
     *         result.
     */
    G* getResultNonBlocking(bool& validResult){
        G* r = NULL;
        G** pr = &r;
        if(!FarmBase<I,O>::_farm->load_result_nb((void**)pr) ||
           isManagementTask((void*)*pr)){
            validResult = false;
            return NULL;
        }else{
            validResult = true;
            if(*pr == EOS){
                *pr = NULL;
            }
            return *pr;
        }
    }

};

/**
 * @brief Accelerator with no collector.
 */
template <typename S, typename I>
class FarmAccelerator<S, I, std::nullptr_t, std::nullptr_t>: public FarmAcceleratorBase<S, I, std::nullptr_t, std::nullptr_t>{
public:
    /**
     * The constructor of the farm.
     * @param parameters The configuration parameters.
     */
    explicit FarmAccelerator(const Parameters* parameters):
            FarmAcceleratorBase<S, I, std::nullptr_t, std::nullptr_t>(parameters){;}

    /**
     * The constructor of the farm.
     * @param paramFileName The filename of the XML file containing the
     *        configuration parameters.
     */
    explicit FarmAccelerator(const std::string& paramFileName):
            FarmAcceleratorBase<S, I, std::nullptr_t, std::nullptr_t>(paramFileName){;}
};

/**
 * @brief Accelerator with collector but no output channel.
 */
template <typename S, typename I, typename O>
class FarmAccelerator<S, I, O, std::nullptr_t>: public FarmAcceleratorBase<S, I, O, std::nullptr_t>{
public:
    /**
     * The constructor of the farm.
     * @param parameters The configuration parameters.
     */
    explicit FarmAccelerator(const Parameters* parameters):
            FarmAcceleratorBase<S, I, O, std::nullptr_t>(parameters){;}

    /**
     * The constructor of the farm.
     * @param paramFileName The filename of the XML file containing the
     *        configuration parameters.
     */
    explicit FarmAccelerator(const std::string& paramFileName):
            FarmAcceleratorBase<S, I, O, std::nullptr_t>(paramFileName){;}

    /**
     * Adds the gatherer to the farm.
     * @param g The gatherer of the farm.
     */
    void addGatherer(Gatherer<O>* g){
        FarmBase<I, O>::setGatherer(g);
    }
};


typedef struct ParallelForRange{
    long long int start;
    long long int end;
    long long int step;
}ParallelForRange;

extern ParallelForRange terminationRange; // Just a dummy value to signal termination

class ParallelForScheduler: public nornir::Scheduler<ParallelForRange, ParallelForRange>{
public:
    ParallelForRange* schedule(ParallelForRange* r){
        if(r == &terminationRange){
            disableRethreading();
            // We record how many workers were active when termination
            // range was sent.
            terminationRange.start = getCurrentNumWorkers();
            broadcast(&terminationRange);
            enableRethreading();
            return nothing();
        }else{
            return r;
        }
    }
};

class ParallelForWorker: public nornir::Worker<ParallelForRange, ParallelForRange>{
private:
    std::function<void(unsigned long long, unsigned long long)> _function;
public:
    explicit ParallelForWorker(){;}

    void setFunction(const std::function<void(unsigned long long, unsigned long long)>& function){
        _function = function;
    }

    ParallelForRange* compute(ParallelForRange* range) {
        if(range == &terminationRange){
            // We only forward the termination range so the gatherer can sleep
            // while the workers are working.
            return range;
        }else{
            for(long long int i = range->start; i < range->end; i += range->step){
                _function(i, getId());
            }
            return nothing();
        }
    }
};

class ParallelForGatherer: public nornir::Gatherer<ParallelForRange, ParallelForRange>{
public:
    ParallelForRange* gather(ParallelForRange* r){
        if(r == &terminationRange){
            return r;
        }else{
            throw std::runtime_error("ParallelForGatherer is supposed to receive only terminationRange.");
        }
    }
};

class ParallelFor: public mammut::utils::NonCopyable{
private:
    FarmAccelerator<ParallelForRange, ParallelForRange, ParallelForRange, ParallelForRange>* _acc;
    std::vector<ParallelForWorker*> _workers;
    size_t _numThreads;

    void pause(){
        long long int receivedTerminations = 0;
        _acc->offload(&terminationRange);
        ParallelForRange* r;
        do{
            r = _acc->getResult();
            if(r != &terminationRange){
                throw std::runtime_error("ParallelFor getResult is supposed to receive only terminationRange.");
            }
            ++receivedTerminations;
        }while(receivedTerminations < r->start); // r->start is the number of expected ranges
        // TODO Find a way to freeze the threads of the farm without destroying
        // everything (it should be enough to call the freeze and run methods of knobworkersfastflow)
    }

    void resume(){;}

public:
    // TODO: Dire quand'è che possiamo prendere un sample: piu sample per ogni tipo di loop, un sample per ogni tipo di loop, un sample per il blocco di loop
    ParallelFor(unsigned long int numThreads,
                nornir::Parameters* parameters){
        _acc = new FarmAccelerator<ParallelForRange, ParallelForRange, ParallelForRange, ParallelForRange>(parameters);
        for(unsigned long int i = 0; i < numThreads; i++){
            _workers.push_back(new ParallelForWorker());
            _acc->addWorker(_workers.back());
        }
        _acc->addScheduler(new ParallelForScheduler());
        _acc->addGatherer(new ParallelForGatherer());
        _acc->setOndemandScheduling();
        _numThreads = numThreads;
        _acc->start();
        pause();
    }

    ~ParallelFor(){
        _acc->shutdown();
        _acc->wait();
        delete _acc;
        for(auto w : _workers){
            delete w;
        }
    }

    inline void parallel_for(long long int start, long long int end, long long int step,
                             long int chunkSize, const std::function<void(unsigned long long, unsigned long long)>& function){
        // Allocate ranges here so pointers will be valid for all the function duration.
        // We need list because with vector we invalidate pointers when doing push_back
        std::list<ParallelForRange> ranges;
        if(step > (end - start)){
            throw std::runtime_error("parallel_for: step cannot be greater than iteration space.");
        }
        for(auto w : _workers){
            w->setFunction(function);
        }

        if(!chunkSize){
            unsigned long long numIterations = std::ceil((end - start)/(double) step);
            chunkSize = std::ceil(numIterations / (double) _numThreads);
        }

        resume();

        ParallelForRange pfr;
        bool setStart = true;
        long long numIterations = 0;
        long long lastValidId = 0;
        for(long long int i = start; i < end; i += step){
            lastValidId = i;
            if(setStart){
                pfr.start = i;
                setStart = false;
            }
            if(++numIterations == chunkSize){
                pfr.end = i + 1; // +1 because this iteration needs to be computed
                pfr.step = step;
                ranges.push_back(pfr);
                _acc->offload(&(ranges.back()));

                setStart = true;
                numIterations = 0;
            }
        }
        // Spurious element
        if(!setStart){
            //TODO Avoid duplicated code)
            pfr.end = lastValidId + 1; // +1 because this iteration needs to be computed
            pfr.step = step;
            ranges.push_back(pfr);
            _acc->offload(&(ranges.back()));
        }
        pause();
    }
};


/**
 * @param chunkSize If 0, iteration space is statically divided among threads,
 * i.e. each thread gets numIterations/numThreads iterations
 **/
template <typename Function>
inline void parallel_for(long long int start, long long int end, long long int step,
                         long int chunkSize, unsigned long int numThreads,
                         nornir::Parameters* parameters, const Function& function){
    ParallelFor pf(numThreads, parameters);
    pf.parallel_for(start, end, step, chunkSize, function);

}


template <typename Function>
inline void parallel_for(long long int start, long long int end, long long int step,
                         long int chunkSize, unsigned long int numThreads,
                         std::string parametersFile, const Function& function){
    Parameters p(parametersFile);
    ParallelFor pf(numThreads, &p);
    pf.parallel_for(start, end, step, chunkSize, function);
}

}

#endif /* NORNIR_INTERFACE_HPP_ */

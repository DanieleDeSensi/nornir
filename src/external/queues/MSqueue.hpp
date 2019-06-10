#ifndef MSQUEUE_H
#define MSQUEUE_H
#include "align.h"
#include "hzdptr.h"
#include "primitives.h"

// cppcheck-suppress noConstructor
class MSQueue {
private:
#define EMPTY (void *) -1
    
    typedef struct _node_t {
	struct _node_t * volatile next DOUBLE_CACHE_ALIGNED;
	void * data DOUBLE_CACHE_ALIGNED;
    } node_t DOUBLE_CACHE_ALIGNED;
    
    typedef struct _queue_t {
	struct _node_t * volatile head DOUBLE_CACHE_ALIGNED;
	struct _node_t * volatile tail DOUBLE_CACHE_ALIGNED;
	int nprocs;
    } queue_t DOUBLE_CACHE_ALIGNED;
    
    typedef struct _handle_t {
	hzdptr_t hzd;
    } handle_t DOUBLE_CACHE_ALIGNED;
    
    queue_t * q = NULL;
    handle_t ** hds;

private:
        
    void enqueue(handle_t * handle, void * data) {
	node_t * node = (node_t*)malloc(sizeof(node_t));
	
	node->data = data;
	node->next = NULL;
	
	node_t * tail;
	node_t * next;
	
	while (1) {
	    tail = (node_t*)hzdptr_setv(&q->tail, &handle->hzd, 0);
	    next = tail->next;
	    
	    if (tail != q->tail) {
		continue;
	    }
	    
	    if (next != NULL) {
		CAS(&q->tail, &tail, next);
		continue;
	    }
	    
	    if (CAS(&tail->next, &next, node)) break;
	}
	
	CAS(&q->tail, &tail, node);
    }
    
    void * dequeue(handle_t * handle)  {
	void * data;
	
	node_t * head;
	node_t * tail;
	node_t * next;
	
	while (1) {
	    head = (node_t*)hzdptr_setv(&q->head, &handle->hzd, 0);
	    tail = q->tail;
	    next = (node_t*)hzdptr_set(&head->next, &handle->hzd, 1);
	    
	    if (head != q->head) {
		continue;
	    }
	    
	    if (next == NULL) {
		return (void *) -1;
	    }
	    
	    if (head == tail) {
		CAS(&q->tail, &tail, next);
		continue;
	    }
	    
	    data = next->data;
	    if (CAS(&q->head, &head, next)) break;
	}
	
	hzdptr_retire(&handle->hzd, head);
	return data;
    }

    void queue_init(int nprocs) {
	node_t * node = (node_t*)malloc(sizeof(node_t));
	node->next = NULL;
	
	q->head = node;
	q->tail = node;
	q->nprocs = nprocs;
    }
    
    void queue_register(handle_t * th, int id)  {
	hzdptr_init(&th->hzd, q->nprocs, 2);
    }
      
public:

    // must be called once before starting using the queue
    bool init(int nprocs) {
	if (q != NULL) return false;
	if (nprocs > MAX_PROCS) return false;

	q = (queue_t*)align_malloc(PAGE_SIZE, sizeof(queue_t));
	if (!q) return false;
	queue_init(nprocs);
	hds = (handle_t**)align_malloc(PAGE_SIZE, sizeof(handle_t * [nprocs]));
	if (!hds) {
	    return false;
	}
	return true;
    }

    // called once by a thread before calling push/pop
    void registerq(int pid) {
	hds[pid] = (handle_t*)align_malloc(PAGE_SIZE, sizeof(handle_t));
	queue_register(hds[pid], pid);
    }

    void deregisterq(int pid) {}

    bool push(void *data, int pid) {
	enqueue(hds[pid], data);
	//P[pid] += 1;
	return true;
    }

    bool pop(void **val, int pid) {
	*val = dequeue(hds[pid]);
	//if (val>=0) C[pid] += 1;
	return (*val == EMPTY) ? false: true;
    }

    void printStats(int nthreads) {
	size_t sum = 0;
	for(int i=0;i<nthreads;++i) {
	    printf("TH%d n push = %ld, n pop = %ld\n", i, P[i], C[i]);
	    sum += C[i];
	}
	printf("Total elements read: %zu\n", sum);
    }

private:
    // used only for statistics
    long P[MAX_PROCS] CACHE_ALIGNED;
    long C[MAX_PROCS] CACHE_ALIGNED;
};



#endif /* end of include guard: MSQUEUE_H */

#ifndef AOS_LIB_THREAD_POOL_H
#define AOS_LIB_THREAD_POOL_H
typedef void*(*thread_fun)(void*);
typedef struct thread_pool_t thread_pool_t;
void thread_pool_init(int nthread);
void thread_pool_queue(thread_fun fun, void *arg);
void thread_pool_wait(void);
void thread_pool_destroy(void);
#endif

#ifndef AOS_LIB_STFUN
#define AOS_LIB_STFUN
typedef struct stfun_t stfun_t;
stfun_t* stfun_init(long nx, long ny, double *amp);
void stfun_push(stfun_t *A, dmat *opd);
dmat *stfun_finalize(stfun_t *A);
dmat* stfun_kolmogorov(loc_t *loc, double r0);
#endif

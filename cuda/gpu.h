#ifndef AOS_CUDA_GPU_H
#define AOS_CUDA_GPU_H
#include "../lib/aos.h"
#include "../maos/parms.h"
#include "../maos/types.h"
/* Data for propagation to WFS */
typedef struct{
    dmat *phi;//return.
    int iwfs;
    int isim;
    const PARMS_T *parms;
    const POWFS_T *powfs;
    float atmalpha;
    float dmalpha;
}gpu_wfs_t;

void gpu_info(void);
void gpu_atm2gpu(map_t **atm, int nps);
void gpu_dm2gpu(map_t **dmreal, int ndm, DM_CFG_T *dmcfg);
void gpu_saloc2gpu(int npowfs, int ipowfs, int nwfs, loc_t *loc);
void gpu_plocs2gpu(loc_t *plocs, dmat *amp);
void gpu_print_cc(char *msg);
void gpu_wfs(gpu_wfs_t *info);
void gpu_perfevl(thread_t *info);
#endif

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
extern struct cumap_t *cudmreal;
extern struct cumap_t *cudmproj;
void gpu_info(void);
int gpu_init(int igpu);
void gpu_cleanup(void);
void gpu_atm2gpu(map_t **atm, int nps);
void gpu_dm2gpu(struct cumap_t **cudm, map_t **dmreal, int ndm, DM_CFG_T *dmcfg);
void gpu_telws2gpu(dmat *telws);
void gpu_wfsgrad_init(const PARMS_T *parms, const POWFS_T *powfs);
void gpu_wfsgrad_seeding(const PARMS_T *parms, const POWFS_T *powfs, rand_t *rstat);
void gpu_wfssurf2gpu(const PARMS_T *parms, POWFS_T *powfs);
void gpu_evlsurf2gpu(APER_T *aper);
void gpu_perfevl_init(const PARMS_T *parms, APER_T *aper);
void gpu_wfsgrad(thread_t *info);
void gpu_perfevl(thread_t *info);
void gpu_perfevl_save(SIM_T *simu);
void gpu_setup_recon(const PARMS_T *parms, POWFS_T *powfs, RECON_T *recon);
void gpu_tomofit(SIM_T *simu);
#endif


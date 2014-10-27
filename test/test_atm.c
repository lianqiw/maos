/*
  Test the difference in generated atmosphere in taurus and opteron.
*/
#include "../lib/aos.h"
int main(){
    rand_t rstat;
    rand_t bstat;
    seed_rand(&bstat,1);
    seed_rand(&rstat,lrand(&bstat));
    int nlayer=2;
    double wt[7]={0.2887, 0.17795, 0.06602, 0.07833, 0.1405, 0.1216, 0.1269};
    double r0=0.1987; 
    double L0=30;
    double dx=1./64.;
    int m=4096*2;
    NTHREAD=6;
    THREAD_POOL_INIT(NTHREAD);
    GENATM_T data;
    memset(&data, 0, sizeof(GENATM_T));
    data.rstat=&rstat;
    data.nx=m;
    data.ny=m;
    data.dx=dx;
    data.r0=r0;
    data.l0=L0;
    data.wt=wt;
    data.nlayer=nlayer;
    mapcell *map=vonkarman_screen(&data);
    writebin(map, "atm.bin");
}

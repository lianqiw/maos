/*
  Test the difference in generated atmosphere in taurus and opteron.
*/
#include "../lib/aos.h"
int main(){
    rand_t rstat;
    rand_t bstat;
    seed_rand(&bstat,1);
    seed_rand(&rstat,lrand(&bstat));
    int nthread=6;
    int nlayer=2;
    double wt[7]={0.2887, 0.17795, 0.06602, 0.07833, 0.1405, 0.1216, 0.1269};
    double r0=0.1987; 
    double L0=30;
    double dx=1./64.;
    int m=4096*2;
    dmat *spect=turbpsd(m, m, dx, r0, L0,0.5);
    dwrite(spect, "spect");
    THREAD_POOL_INIT(nthread);
    map_t **map=genscreen_from_spect(&rstat, spect, r0,L0,dx, wt, nlayer, nthread);
    maparrwrite(map, nlayer, "atm.bin");
}

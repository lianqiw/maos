/*
  Test the difference in generated atmosphere in taurus and opteron.
*/
#include "../lib/aos.h"
int main(){
    struct_rand rstat;
    struct_rand bstat;
    seed_rand(&bstat,1);
    seed_rand(&rstat,lrand(&bstat));
    int nthread=6;
    int nlayer=2;
    double wt[7]={0.2887, 0.17795, 0.06602, 0.07833, 0.1405, 0.1216, 0.1269};
    double r0=0.1987; 
    double L0=30;
    double dx=1./64.;
    int m=4096*2;
    dmat *spect=vonkarman_spect(m, m, dx, r0, L0);
    dwrite(spect, "spect");
    if(nthread>1){
#if USE_PTHREAD == 2
	default_pool=thr_pool_create(1,nthread,3600,NULL);
#endif
    }
    MAP_T **map=genscreen_from_spect(&rstat, spect, dx, wt, nlayer, nthread);
    sqmaparrwrite(map, nlayer, "atm.bin");
}

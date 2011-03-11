#include "../lib/aos.h"
int main(){
    ccell *psf;
    int nsim=5000;
    for(int iwfs=6; iwfs<64; iwfs++){
	cellarr *ca=cellarr_init(nsim,1,"psfout_wfs%d.bin",iwfs);
	for(int isim=0; isim<nsim; isim++){
	    if(isim%100==0)
		info("iwfs=%d, isim=%d\n",iwfs,isim);
	    psf=ccellread("psfout/psfout_wfs%d_isim%d.bin.gz",iwfs,isim);
	    cellarr_ccell(ca,psf);
	    ccellfree(psf);
	}
	cellarr_close(ca);
    }
}

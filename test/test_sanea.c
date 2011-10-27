#include "../aos.h"
char *dirout;
static void test_ints(){
    rand_t init;
    seed_rand(&init,1);
    const int nwfs=6;
    lrand(&init);/*atm */
    rand_t wfs_rand[nwfs];
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	seed_rand(wfs_rand+iwfs,lrand(&init));
    }
    dcell *mtche=dcellread("powfs0_mtche.bin");
    int nsim=500;
    int nsa=2582;
    dmat *nea=dnew(nsim,nsa*2);
    double(*pnea)[nsim]=(void*)nea->p;
    double rne=3;
    double bkgrnd=0;
    double siglev=1000;
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	for(int isim=0; isim<nsim; isim++){
	    info("iwfs %d isim=%d\n",iwfs,isim);
	    /*dcell *ints=dcellread("ints_%d_wfs%d.bin",isim,iwfs); */
	    dcell *ints=dcellread("ints_%d_wfs%d.bin",isim,iwfs);
	    /*dcell *i0=dcellread("powfs0_i0.bin"); */
	    dmat *im=NULL, *imy=NULL;
	    double gnf[2], gny[2];

	    for(int isa=0; isa<nsa; isa++){
		dcp(&im,ints->p[isa]);
		dscale(im,siglev);
		gnf[0]=0; gnf[1]=0;
		dmulvec(gnf, mtche->p[isa], im->p,1.);
		gny[0]=0; gny[1]=0;
		dcp(&imy, im);
		addnoise(imy, &wfs_rand[iwfs], bkgrnd,1.,rne);
		dmulvec(gny, mtche->p[isa], imy->p,1.);
		pnea[isa][isim]=gny[0]-gnf[0];
		pnea[isa+nsa][isim]=gny[1]-gnf[1];
	    }
	}
	dwrite(nea,"test_sanea_wfs%d.bin",iwfs);
    }
}
/*
static void test_i0(){
    rand_t i0rand;
    seed_rand(&i0rand,100);
    dcell *mtche=dcellread("powfs0_mtche.bin");
    int nsim=500;
    int nsa=2582;
    dmat *nea=dnew(nsim,nsa*2);
    double(*pnea)[nsim]=(void*)nea->p;
    double rne=3;
    double bkgrnd=0;
    double siglev=1000;
    dcell *i0=dcellread("powfs0_i0.bin");
    for(int isim=0; isim<nsim; isim++){
	info("isim=%d\n",isim);
	dmat *im=NULL, *imy=NULL;
	double gnf[2], gny[2];
	for(int isa=0; isa<nsa; isa++){
	    dcp(&im,i0->p[isa]);
	    dscale(im,siglev);
	    gnf[0]=0; gnf[1]=0;
	    dmulvec(gnf, mtche->p[isa], im->p,1.);
	    gny[0]=0; gny[1]=0;
	    dcp(&imy, im);
	    addnoise(imy, &i0rand, bkgrnd,1.,rne);
	    dmulvec(gny, mtche->p[isa], imy->p,1.);
	    pnea[isa][isim]=gny[0]-gnf[0];
	    pnea[isa+nsa][isim]=gny[1]-gnf[1];
	}
    }
    dwrite(nea,"test_sanea_i0.bin");
    }*/
int main(){
    test_ints();
}

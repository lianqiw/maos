/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
  This file is part of Multithreaded Adaptive Optics Simulator (MAOS).

  MAOS is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  MAOS is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  MAOS.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
\file genotf.c
   Routines to generate short exposure OTFs of an aperture in present of
   atmosphere turbulence.
*/
#include "common.h"
#include "thread.h"
#include "daemonize.h"
#include "dmat.h"
#include "cmat.h"
#include "mathmisc.h"
#include "loc.h"
#include "genotf.h"
/**
private data struct to mark valid pairs of points.  */
typedef struct T_VALID{
    int n;
    int (*loc)[2];
}T_VALID;
/**
 Wrap the data to genotf to have multi-thread capability.*/
typedef struct GENOTF_T{
    long isa;
#if USE_PTHREAD > 0
    pthread_mutex_t mutex_isa;
#endif
    cmat **otf;
    loc_t *loc;     /**<the common aperture grid*/
    const double *amp;    /**<The amplitude map of all the (sub)apertures*/
    const double *opdbias;/**<The static OPD bias. */
    const double *area;   /**<area of a (sub)aperture*/
    double thres;/**<threshold to consider an (sub)aperture as full*/
    double wvl;  /**<The wavelength. only needef if opdbias is not null*/
    long ncompx; /**<Size of OTF*/
    long ncompy; /**<Size of OTF*/
    long nsa;   /**<Number of (sub)apertures*/
    long pttr;  /**<Remove piston/tip/tilt*/

    double *B;
    const T_VALID *pval;
    long isafull;
    const cmat *otffull;
}GENOTF_T;
/**
   Remove tip/tilt from the B matrix
*/
static double* pttr_B(const double *B0,   /**<The B matrix. */
		      loc_t *loc,  /**<The aperture grid*/
		      const double *amp /**<The amplitude map*/
		   ){
    double *locx=loc->locx;
    double *locy=loc->locy;
    int nloc=loc->nloc;
    double (*B)[nloc]=(double(*)[nloc]) B0;
  
    char transn='N';
    char transt='T';
    int nmod=3;
    double *mod[3];
    dmat *mcc=dnew(3,3);
    double (*cc)[3]=(double(*)[3])mcc->p;
    double (*restrict BP)[nloc]=malloc(sizeof(double)*nloc*nloc);
    
    double *restrict M, * restrict MW;
    double (*restrict MCCT)[3],(*restrict Mtmp)[3];
    double dpone=1;
    double dpnone=-1;
    double dpzero=0;
    mod[0]=NULL;
    mod[1]=locx;
    mod[2]=locy;
    for(int im=0; im<3;im++){
	for(int jm=im;jm<3;jm++){
	    cc[jm][im]=dotdbl(mod[im], mod[jm], amp, nloc);
	    if(im!=jm)
		cc[im][jm]=cc[jm][im];
	}
    }
    dinvspd_inplace(mcc);
    M   = malloc(sizeof(double)*nloc*3);
    MW  = malloc(sizeof(double)*nloc*3);
    MCCT= malloc(sizeof(double)*nloc*3);
    Mtmp= malloc(sizeof(double)*nloc*3);
    for(int iloc=0; iloc<nloc; iloc++){
	M[iloc]=1;
    }
    memcpy(M+nloc, locx, nloc*sizeof(double));
    memcpy(M+nloc*2, locy, nloc*sizeof(double));
    for(int iloc=0; iloc<nloc; iloc++){
	MW[iloc]=amp[iloc];
	MW[iloc+nloc]=amp[iloc]*locx[iloc];
	MW[iloc+nloc*2]=amp[iloc]*locy[iloc];
    }
    /* MCCT = - cci'*M' */
    dgemm_(&transt, &transt, 
	   &nmod, &nloc, &nmod, 
	   &dpnone, 
	   (double*)cc,&nmod,
	   M, &nloc, 
	   &dpzero,
	   (double*)MCCT,&nmod);
    /* Mtmp = MW'*B */
    dgemm_(&transt, &transn,
	   &nmod, &nloc, &nloc,
	   &dpone, 
	   MW, &nloc,
	   (double*)B, &nloc,
	   &dpzero, 
	   (double*)Mtmp, &nmod);
    for(int iloc=0; iloc<nloc; iloc++){
	double tmp1=Mtmp[iloc][0];
	double tmp2=Mtmp[iloc][1];
	double tmp3=Mtmp[iloc][2];
	for(int im=0; im<nloc; im++){
	    BP[iloc][im]=B[iloc][im]+
		(MCCT[im][0]*tmp1+MCCT[im][1]*tmp2+MCCT[im][2]*tmp3);
	}
    }
    /* Mtmp = MW'*BP' */
    dgemm_(&transt, &transt,
	   &nmod, &nloc, &nloc,
	   &dpone, 
	   MW, &nloc,
	   (double*)BP, &nloc,
	   &dpzero, 
	   (double*)Mtmp, &nmod);

    //double *restrict BPD=malloc(sizeof(double)*nloc);
    for(int iloc=0; iloc<nloc; iloc++){
	double tmp1=MCCT[iloc][0];
	double tmp2=MCCT[iloc][1];
	double tmp3=MCCT[iloc][2];
	/*BPD[iloc]=exp(BP[iloc][iloc]
		      +tmp1*Mtmp[iloc][0]
		      +tmp2*Mtmp[iloc][1]
		      +tmp3*Mtmp[iloc][2]);*/
	for(int jloc=0; jloc<nloc; jloc++){
	    double tmp=BP[iloc][jloc]+tmp1*Mtmp[jloc][0]
		+tmp2*Mtmp[jloc][1]+tmp3*Mtmp[jloc][2];
	    BP[iloc][jloc]=exp(-2.*tmp);
	}
    }
    dfree(mcc);
    free(M);
    free(MW);
    free(MCCT);
    free(Mtmp);
    return (double*)BP;
}
/**
   Generate OTF from the B or tip/tilted removed B matrix.
*/
static void genotf_do(cmat **otf, int pttr, long notfx, long notfy, 
		      loc_t *loc, const double *amp, const double *opdbias, double wvl,
		      const double* B,  const T_VALID *pval){
    int nloc=loc->nloc;
    double (*BP)[nloc];
    if(pttr){//remove p/t/t from the B matrix
	BP=(void*)pttr_B(B,loc,amp);
    }else{
	BP=(void*)B;
    }
    if(!*otf){
	*otf=cnew(notfx,notfy);
    }
 
    PCMAT(*otf,OTF);
 

    double *restrict BPD=malloc(sizeof(double)*nloc);
    for(int iloc=0; iloc<nloc; iloc++){
	BPD[iloc]=pow(BP[iloc][iloc],-0.5);
    }

    double otfnorm;
    otfnorm=0;
    for(int iloc=0; iloc<nloc; iloc++){
	otfnorm+=amp[iloc]*amp[iloc];
    }

    otfnorm=1./otfnorm;
    struct T_VALID (*qval)[notfx]=(struct T_VALID (*)[notfx])pval;
    if(opdbias){
	dcomplex wvk=2.*M_PI/wvl*I;
	for(int jm=0; jm<notfy; jm++){
	    for(int im=0; im<notfx; im++){
		int (*jloc)[2]=qval[jm][im].loc;
		dcomplex tmp1,tmp2;
		register dcomplex tmp=0.;
		for(int iloc=0; iloc<qval[jm][im].n; iloc++){
		    int iloc1=jloc[iloc][0];//iloc1 is continuous.
		    int iloc2=jloc[iloc][1];//iloc2 is not continuous.
		    tmp1=amp[iloc1]*cexp(wvk*opdbias[iloc1])*BPD[iloc1]*BP[iloc1][iloc2];
		    tmp2=amp[iloc2]*cexp(-wvk*opdbias[iloc2])*BPD[iloc2];
		    tmp+=tmp1*tmp2;
		}
		OTF[jm][im]=tmp*otfnorm;
	    }
	}
    }else{
	for(int jm=0; jm<notfy; jm++){
	    for(int im=0; im<notfx; im++){
		int (*jloc)[2]=qval[jm][im].loc;
		double tmp1,tmp2;
		register double tmp=0.;
		for(int iloc=0; iloc<qval[jm][im].n; iloc++){
		    int iloc1=jloc[iloc][0];//iloc1 is continuous.
		    int iloc2=jloc[iloc][1];//iloc2 is not continuous.
		    tmp1=amp[iloc1]*BPD[iloc1]*BP[iloc1][iloc2];
		    tmp2=amp[iloc2]*BPD[iloc2];
		    tmp+=tmp1*tmp2;
		}
		OTF[jm][im]=tmp*otfnorm;
	    }
	}
    }
    free(BPD);
    if((void*)BP!=(void*)B){
	free(BP);
    }
}
/**
   A wrapper to execute pttr parallel in pthreads
 */
static void *genotf_wrap(GENOTF_T *data){
    long isa=0;
    const long nsa=data->nsa;
    cmat**otf=(cmat**)data->otf;
    loc_t *loc=data->loc;
    const long nxsa=loc->nloc;
    const double wvl=data->wvl;
    const long ncompx=data->ncompx;
    const long ncompy=data->ncompy;
    const double *area=data->area;
    const double thres=data->thres;
    const cmat *otffull=data->otffull;
    const double *amp=data->amp;
    const int pttr=data->pttr;
    double *B=data->B;
    const T_VALID *pval=data->pval;
    while(LOCK(data->mutex_isa),isa=data->isa++,UNLOCK(data->mutex_isa),isa<nsa){
	if(!detached){
	    fprintf(stderr,"%6ld of %6ld\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", isa,nsa);
	}
	const double *opdbiasi=NULL;
	if(data->opdbias){
	    opdbiasi=data->opdbias+isa*nxsa;
	}else{
	    opdbiasi=NULL;
	}
	if(otffull && area[isa]>thres){
	    ccp(&otf[isa],otffull);//just copy the full array
	}else{ 
	    genotf_do(&otf[isa],pttr,ncompx,ncompy,loc,amp+isa*nxsa,opdbiasi,wvl,B,pval);
	}
    }
    return NULL;
}
/**
   Generate pairs of overlapping points for structure function. 
   Changelog:

   2010-11-08: removed amp. It caused wrong otf because it uses the amp of the
   first subaperture to build pval, but this one is not fully illuminated. 
 */
static T_VALID *gen_pval(int notfx, int notfy, loc_t *loc,
			 double dtheta, double wvl){
    double dux=1./(dtheta*notfx);
    double duy=1./(dtheta*notfy);
    int nloc=loc->nloc;
    double *locx=loc->locx;
    double *locy=loc->locy;
    const int pvaltot=notfx*notfy*nloc*2;
    int (*pval0)[2]=malloc(sizeof(int)*pvaltot);
    T_VALID *pval=malloc(sizeof(T_VALID)*notfx*notfy);
    T_VALID (*restrict qval)[notfx]=(T_VALID (*)[notfx])(pval);
    int count=0,count2;
    loc_create_map(loc);
    locmap_t *map=loc->map;
    int notfx2=notfx/2;
    int notfy2=notfy/2;
    double duxwvl=dux*wvl;
    double duywvl=duy*wvl;
    double dx1=1./loc->dx;
    long (*mapp)[map->nx]=(long(*)[map->nx])map->p;
    for(int jm=0; jm<notfy; jm++){
	int jm2=(jm-notfy2);//peak in the center
	//int jm2=jm<notfy2?jm:jm-notfy;//peak in the corner
	for(int im=0; im<notfx; im++){
	    int im2=(im-notfx2);
	    //int im2=im<notfx2?im:im-notfx;
	    count2=count;
	    for(int iloc=0; iloc<loc->nloc; iloc++){
		int iy=(int)round((locy[iloc]+jm2*duywvl-map->oy)*dx1);
		int ix=(int)round((locx[iloc]+im2*duxwvl-map->ox)*dx1);
		if (ix>=0 && ix<map->nx && iy>=0 && iy<map->ny) {
		    long iloc2=mapp[iy][ix];
		    if(iloc2--){
			pval0[count][0]=iloc;
			pval0[count][1]=iloc2;
			count++;
		    }
		}
	    }
	    qval[jm][im].loc=pval0+count2;
	    qval[jm][im].n=count-count2;
	}
    }
    loc_free_map(loc);
    //pval0=realloc(pval0, sizeof(int)*count*2);
    //can not realloc. will change position.
    return pval;
}
/**
   Generate the B matrix.
*/
static double* genotfB(loc_t *loc, double wvl, double r0, double L0){
    (void)L0;
    long nloc=loc->nloc;
    double *locx=loc->locx;
    double *locy=loc->locy;
    double(*B)[nloc]=(double(*)[nloc])malloc(sizeof(double)*nloc*nloc);
    const double coeff=6.88*pow(0.5e-6/wvl,2)*pow(r0,-5./3.)*0.25;
    for(int i=0; i<nloc; i++){
	for(int j=i; j<nloc; j++){
	    double rdiff2=pow(locx[i]-locx[j],2)+pow(locy[i]-locy[j],2);
	    B[j][i]=B[i][j]=coeff*pow(rdiff2,5./6.);
	}
    }
    return (double*)B;
}
/**
   Generate OTFs for multiple (sub)apertures. ALl these apertures must share the
   same geometry, but may come with different amplitude map and OPD biasas. if
   pttr is 1, the OTF will have tip/tilt removed. make r0 to infinity to build
   diffraction limited OTF. make r0 to infinity and opdbias to none null to
   build OTF for a static map.*/
void genotf(cmat **otf,    /**<The otf array for output*/
	    loc_t *loc,    /**<the aperture grid (same for all apertures)*/
	    const double *amp,   /**<The amplitude map of all the (sub)apertures*/
	    const double *opdbias,  /**<The static OPD bias. */
	    const double *area, /**<normalized area of the (sub)apertures*/
	    double thres,  /**<The threshold to consider a (sub)aperture as full*/
	    double wvl,    /**<The wavelength. only needef if opdbias is not null*/
	    double dtheta, /**<Sampling of PSF.*/
	    double r0,     /**<Fried parameter*/
	    double l0,     /**<Outer scale*/
	    long ncompx,    /**<Size of OTF*/
	    long ncompy,    /**<Size of OTF*/
	    long nsa,       /**<Number of (sub)apertures*/
	    long pttr,      /**<Remove piston/tip/tilt*/
	    long nthread    /**<Number of threads*/
	     ){
   /*creating pairs of points that both exist with given separation for
      computing structure function*/
    T_VALID *pval=gen_pval(ncompx, ncompy, loc, dtheta, wvl);//returns T_VALID array.
    /* Generate the B matrix. */
    double *B=genotfB(loc, wvl, r0, l0);
    cmat *otffull=NULL;
    const long nloc=loc->nloc;
    int isafull=-1;
    if(!opdbias && nsa>1){
	double maxarea=0;
	for(int isa=0; isa<nsa; isa++){
	    if(area[isa]>maxarea){
		maxarea=area[isa];
		isafull=isa;
	    }
	}
	if(isafull>0){
	    genotf_do(&otffull,pttr,ncompx,ncompy,loc,amp+isafull*nloc,NULL,wvl,B,pval);
	}
    }
    
    GENOTF_T data;
    memset(&data, 0, sizeof(GENOTF_T));
    data.isa=0;
    PINIT(data.mutex_isa);
    data.otf=otf;
    data.loc=loc;
    data.amp=amp;
    data.opdbias=opdbias;
    data.area=area;
    data.thres=thres;
    data.wvl=wvl;
    data.ncompx=ncompx;
    data.ncompy=ncompy;
    data.nsa=nsa;
    data.pttr=pttr;//was missing.
    data.B=B;
    data.pval=pval;
    data.isafull=isafull;
    data.otffull=otffull;

    CALL(genotf_wrap, &data, nthread);
    cfree(otffull);
    free(B);
    free(pval[0].loc);
    free(pval);
}

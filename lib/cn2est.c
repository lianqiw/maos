/*
  Copyright 2009-2016 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "../math/mathdef.h"
#include "cn2est.h"
#include "mkh.h"
#include "turbulence.h"
#define INTERP_NEAREST 0 /*set to 0 after debugging */
#define MIN_SA_OVERLAP 5 /*minimum of subaperture overlap at this separation*/
#define COV_ROTATE 0     /*1: rotate the covariance and then cut along x. (old method).*/
cn2est_t *cn2est_new(const dmat *wfspair, /**<2n*1 vector for n pair of WFS indices.*/
		     const dmat *wfstheta,/**<nwfs*2: angular direction of each WFS.*/
		     const loc_t *saloc,  /**<nsa*2: Subaperture low left corner coordinates*/
		     const dmat *saa,     /**<nsa*1: Normalized subaperture area*/
		     const double saat,   /**<Threshold for keeping subapertures*/
		     const dmat *hs,      /**<nwfs*1: altitude of guide star*/
		     const dmat *htrecon, /**<Layers height intended for tomography*/
		     int keepht,          /**<2: slodar directly to htrecon,
					   * otherwise: interpolate onto htrecon
					   * from native slodar heights*/
		     double L0            /**<The Outer scale*/
    ){
    info2("Cn2 estimation:");
    /*We need at least a pair */
    if(!wfspair) return 0;
    int npair=wfspair->nx*wfspair->ny;
    if(!npair) return 0;
    if(npair%2==1){
	error("pair must have even number of entries\n");
    }
    /*>>1 is a short cut for /2 */
    int nwfspair=npair>>1;
    cn2est_t *cn2est=calloc(1, sizeof(cn2est_t));
    cn2est->nwfspair=nwfspair;
    int nwfs=wfstheta->nx;
    cn2est->nwfs=nwfs;
    cn2est->nsa=saloc->nloc;
    /*wfscov is a flag for wfs showing whether this wfs participates in
     * covariance computation. */
    cn2est->wfscov=calloc(nwfs, sizeof(int));
    for(int ind=0; ind<npair; ind++){
	int iwfs=(int)wfspair->p[ind];
	if(iwfs<0 || iwfs>=nwfs){
	    error("wfspair has invalid values: %d\n", iwfs);
	}
	cn2est->wfscov[iwfs]=1;
    }
    /*embed has dimension nembed*nembed. It is used to embed 1-d gradient vector
      to a 2-d map for each curvature, so that we can use FFT to compute the covariance*/
    cn2est->nembed=0;
    cn2est->embed=loc_create_embed(&cn2est->nembed, saloc, 2, 0);
    const long nx=cn2est->nembed;
    const long nxnx=nx*nx;
    /* mask is a mask defined on square grid for subapertures that have
     * normalized area above the threshold.*/
    lmat *mask=lnew(nx,nx);
    if(saa && saa->nx!=cn2est->nsa){
	error("saa and saloc mismatch\n");
    }
    
    double saat2=dmax(saa)*saat;
    for(int isa=0; isa<cn2est->nsa; isa++){
	if(!saa || saa->p[isa]>saat2){
	    mask->p[cn2est->embed->p[isa]]=1;/*use this subaperture */
	}
    }
    /* this mask for for subapertures that we can compute curvature.*/
    cn2est->mask=lnew(nx,nx);
    cmat *overlap=cnew(nx, nx);
    cmat*  pover=overlap;
    int iymin=nx,iymax=0,ixmin=nx,ixmax=0;
    for(int iy=0; iy<nx; iy++){
	for(int ix=0; ix<nx; ix++){
	    /*Only use a subaperture if we are able to make curvature. */
	    if(IND(mask,ix,iy) && IND(mask,ix+1,iy) && IND(mask,ix-1,iy) &&
	       IND(mask,ix,iy+1) && IND(mask,ix,iy-1)){
		IND(cn2est->mask,ix,iy)=1;
		IND(pover,ix,iy)=1;
		if(ix>ixmax) ixmax=ix;
		if(ix<ixmin) ixmin=ix;
		if(iy>iymax) iymax=iy;
		if(iy<iymin) iymin=iy;
	    }
	}
    }
    int maxsep=MIN((iymax-iymin), (ixmax-ixmin));
    lfree(mask);
    cfft2(overlap, -1);
    for(long i=0; i<nxnx; i++){
	overlap->p[i]=overlap->p[i]*conj(overlap->p[i]);
    }
    cfft2(overlap, 1);
    cfftshift(overlap);
    /*cn2est->overlap: number of overlaps for each separation. peak value at (nx/2, nx/2)*/
    cn2est->overlapi=dnew(nx,nx);
    for(long i=0; i<nxnx; i++){
	double over=creal(overlap->p[i]);
	cn2est->overlapi->p[i]=(over>MIN_SA_OVERLAP)?(nxnx/over):0;
    }
    cfree(overlap);
    /*2-d arrays to store x y gradient, and "curvature" */
    cn2est->gxs=cellnew(nwfs, 1);/*stores gradient in 2-d map */
    cn2est->gys=cellnew(nwfs, 1);
    cn2est->curi=cellnew(nwfs, 1);
    for(int ix=0; ix<nwfs; ix++){
	if(cn2est->wfscov[ix]){
	    cn2est->gxs->p[ix]=dnew(nx, nx);
	    cn2est->gys->p[ix]=dnew(nx, nx);
	    cn2est->curi->p[ix]=cnew(nx, nx);
	}
    }
    /*
      Now we prepare indexes that can be used to compute gradient(curvature)
      cross-covariances conveniently during simulation.
    */
    /*first get a few constants */
    const double dsa=saloc->dx;
    /*determine the layer height used for tomography. */
    cn2est->htrecon=ddup(htrecon);
    cn2est->wtrecon=cellnew(1,1);
    cn2est->wtrecon->p[0]=dnew(cn2est->htrecon->nx,1);
    {
	info2("htrecon=[");
	for(int iht=0; iht<cn2est->htrecon->nx; iht++){
	    info2("%.2f ", cn2est->htrecon->p[iht]*0.001);
	}
	info2("]km\n");
    }
    /*stores cross-covariance data during simulation */
    cn2est->covc=cellnew(nwfspair,1);
    cn2est->cov1=cellnew(nwfspair,1);
    cn2est->cov2=cellnew(nwfspair,1);
    /*height of layers for each wfs pair of slodar output. */
    cn2est->ht=cellnew(nwfspair,1);
    /*record sapair to use for each separation */
    cn2est->pair=calloc(nwfspair, sizeof(CN2PAIR_T));
    long nhtsx[nwfspair]; 
    long nhtsy[nwfspair];
    double hmin, hmax;
    dmaxmin(cn2est->htrecon->p, cn2est->htrecon->nx, &hmax, &hmin);
    /*ovs is the over sampling factor in mc. need to be at least 2 to
      cover the non-zero PSD (sinc functions.)*/
      
    /*deprecated:for ovs bigger than 1, the
      results changes even with wfs pairs aligned along x/y direction. The
      reason is that we are sampling more PSDs with bigger ovs. This is probably
      the aliasing effect. With larger ovs, we are using larger
      frequencies. That gets aliased into the results. Then how about padding?*/
    cn2est->ovs=2;
    /*Pnk stores the opeator from layer weight to curvature covariance
      matrix. iPnk is the psuedo inverse of Pnk. Each diagonal cell is for each
      WFS pair.*/
    cn2est->Pnk=cellnew(nwfspair, nwfspair);
    cn2est->iPnk=cellnew(nwfspair,nwfspair);
    /*wtconvert is the matrix to down/up sample the CN2 estimates to layers
      used for tomography*/
    cn2est->wtconvert=cellnew(1,nwfspair);
    cn2est->L0=L0;
    for(int iwfspair=0; iwfspair<nwfspair; iwfspair++){
	/*get pointer for this pair */
	CN2PAIR_T *pair=cn2est->pair+iwfspair;
	/*The WFS in this pair. */
	const int wfs0=wfspair->p[iwfspair*2];
	const int wfs1=wfspair->p[iwfspair*2+1];
	pair->wfs0=wfs0;
	pair->wfs1=wfs1;
	/*The separation between the stars */
	const double dthetax=IND(wfstheta,wfs0,0)-IND(wfstheta,wfs1,0);
	const double dthetay=IND(wfstheta,wfs0,1)-IND(wfstheta,wfs1,1);
	/*the angular distance between WFS */
	double dtheta=sqrt(dthetax*dthetax+dthetay*dthetay);
	/*The direction of the WFS pair baseline vector */
	double beta=atan2(dthetay, dthetax);
	if(wfs0==wfs1 || dtheta<1e-14){
	    dtheta=0;
	    beta=0;
	}
	/*The average wfs guide star height of the two*/
	const double cb=cos(beta);
	const double sb=sin(beta);
#if COV_ROTATE
	const double slang=1;
#else
	const double slang=MAX(fabs(cb), fabs(sb));
#endif
	double hsm=0;
	if(hs->nx*hs->ny>=nwfs){
	    hsm=(hs->p[wfs0]+hs->p[wfs1])*0.5;
	}else if(hs->nx*hs->ny==1){
	    hsm=hs->p[0];
	}else{
	    error("hs in in wrong format, should have 1 or %d elements\n", nwfs);
	}
	if(dtheta==0){//seeing only for self-correlation.
	    pair->iht0=0;
	    pair->iht1=1;
	    warning("Can only obtain r0 using auto-correlation.\n");
	}else{
	    /*The number of layers is determined */
	    if(keepht==2){
		/*slodar output directly to layers used for tomography. may not work well */
		pair->iht0=0;
		pair->iht1=cn2est->htrecon->nx;
	    }else{
		/*output to nature slodar heights. and then bin to tomography layers */
		pair->iht0=(int)floor(hmin*dtheta*slang/(dsa*(1-hmin/hsm)));
		pair->iht1=(int)floor(hmax*dtheta*slang/(dsa*(1-hmax/hsm)));
		if(pair->iht0<-maxsep){
		    pair->iht0=-maxsep;
		}
		if(pair->iht1>maxsep){
		    pair->iht1=maxsep;
		}
		pair->iht1++;//exclusive
	    }
	}
	pair->nht=(pair->iht1-pair->iht0);
	nhtsx[iwfspair]=pair->nht;
	nhtsy[iwfspair]=1;
#if COV_ROTATE
	pair->nsep=pair->nht;
#else
	pair->nsep=nxnx; //use all combinations.
#endif
	/*initialize ht array to store layer heights */
	cn2est->ht->p[iwfspair]=dnew(pair->nht,1);
	/*stores the cross-covariance of the curvature for x, and y grads. */
	cn2est->covc->p[iwfspair]=cnew(nx, nx);
        cn2est->cov2->p[iwfspair]=dnew(nx, nx);
#if COV_ROTATE
	cn2est->cov1->p[iwfspair]=dnew(pair->nsep,1);//1d cut
#else
        cn2est->cov1->p[iwfspair]=dref_reshape(cn2est->cov2->p[iwfspair], nxnx, 1);
#endif
	/*info2("Pair %d: wfs %d and %d. dtheta=%4f\" iht=[%d, %d)\n", 
	  iwfspair, wfs0, wfs1, pair->dtheta*206265, pair->iht0, pair->iht1);*/
 

	/*
	  Now we start to build the model that will be used to estimate Cn2 from the
	  cross-covariance of curvature.
	*/
	const int nsep=pair->nsep;
	/*nm is the size of the mc array */
	const int nm=4*(MAX(abs(pair->iht0), abs(pair->iht1-1))+1)*cn2est->ovs;
	const int nm2=nm>>1;
	/*sampling in mc after FFT */
	const double dx=dsa/cn2est->ovs;
	const double df=1./(nm*dx);
	const double L02=pow(cn2est->L0,-2);
	/*initialize */
	cmat *mc=cnew(nm,nm);
	/*create 2-d pointers */
	cmat* pmc=mc;
	/*the forward operator from layer weights to cross-covariance */
 	dmat *Pnk=dnew(nsep, pair->nht);
	info2("Pair %d: hk=[", iwfspair);
	for(int iht=pair->iht0; iht<pair->iht1; iht++){
	    /*the height of each layer */
	    double hk;
	    if(keepht==2){//do not use.
		hk=cn2est->htrecon->p[iht-pair->iht0];
	    }else if(dtheta==0){
		hk=0;
	    }else{
		hk=dsa*iht/(dtheta*slang+dsa*iht/hsm);
	    }
	    info2("%.2f ",hk*0.001);
	    cn2est->ht->p[iwfspair]->p[iht-pair->iht0]=hk;
	    /*the cone effect */
	    const double zeta=1.-hk/hsm;
	    const double zetan2=pow(zeta,-2);
	    /*
	      coefficients for the PSD. zeta is cone effect. df is the frequency
	      bin. df*df is from the descretization of the PSD. we use df
	      instead of 1/dx like in the matlab code because our inverse fft
	      does not apply (1/N)^ scaling.
	    */
	    /*the coefficients in fron of the PSD. r0 is left out. */
	    const double psd_coef=pow(zeta,-2)*4*M_PI*M_PI
		*0.0229*pow(2*M_PI/0.5e-6,-2)*(df*df);
	    for(int iy=0; iy<nm; iy++){
		/*0 freq at corner so no need to do fftshift before fft */
		const double fy=(iy-nm2)*df;
		/*the sinc subaperture fnction */
		const double sincfy=sinc(fy*dsa);
		for(int ix=0; ix<nm; ix++){
		    const double fx=(ix<nm2?ix:ix-nm)*df;
		    const double sincfx=sinc(fx*dsa);
		    /*the turbulence PSD with outerscale */
		    const double psd=psd_coef*pow((fx*fx+fy*fy)*zetan2+L02,-11./6.);
		    /*gx diff is along x, gy diff is along y to form real curvature */
		    const double cur=pow(2*fx*(cos(2*M_PI*dsa*fx)-1)+2*fy*(cos(2*M_PI*dsa*fy)-1),2);
		    IND(pmc,ix,iy)=pow(sincfy*sincfx,2)*psd*cur;
		}/*ix */
	    }/*iy */
	    /*doing fft */
	    cfft2(mc,1);/*inverse fft */
	    /*shift 0 freq to center. */
	    cfftshift(mc);
	    for(int isep=0; isep<nsep; isep++){
		/*Different from matlab prototype because I used interpolation. */
#if COV_ROTATE //subaperture separation align with guide star separation baseline
		int xsep=isep*cb+pair->iht0;
		int ysep=isep*sb+pair->iht0;
#else
		int ysep=(isep/nx);
		int xsep=(isep-ysep*nx)-nx/2;
		ysep-=nx/2;
#endif
		double xx=(iht*cb/slang-xsep)*cn2est->ovs + nm2;
		double yy=(iht*sb/slang-ysep)*cn2est->ovs + nm2;

		if(xx<0 || yy<0 || xx+1>=nm || yy+1>=nm){
		    //warning("Out of range: xx=%g, yy=%g\n", xx, yy);
		}else{
#if INTERP_NEAREST
		    /*Do interpolation using nearest neighbor */
		    int ixx=(int)round(xx);
		    int iyy=(int)round(yy);
		    double imc=creal(IND(pmc,ixx,iyy));
#else
		    /*Do interpolation using bilinear spline interp. */
		    int ixx=(int)floor(xx); xx=xx-ixx;
		    int iyy=(int)floor(yy); yy=yy-iyy;
		    double imc=creal((IND(pmc,ixx,iyy)*(1-xx)+IND(pmc,ixx+1,iyy)*(xx))*(1-yy)
				     +(IND(pmc,ixx,iyy+1)*(1-xx)+IND(pmc,ixx+1,iyy+1)*(xx))*yy);
#endif
		    IND(Pnk, isep, iht-pair->iht0)=imc;
		}
	    }
	}
	info2("]km\n");
	/*
	  iPnk is a block diagonal matrix for Cn2 Estimation.
	*/
	cn2est->Pnk->p[iwfspair+iwfspair*nwfspair]=dref(Pnk);
	cn2est->iPnk->p[iwfspair+iwfspair*nwfspair]=dpinv(Pnk,0);
	cn2est->wtconvert->p[iwfspair]=mkhbin1d(cn2est->ht->p[iwfspair],cn2est->htrecon);
	dfree(Pnk);
	cfree(mc);
    }/*iwfspair */
/*stores estimated weight of layers during simulation and output to file finally. */
    cn2est->wt=dcellnew3(nwfspair, 1, nhtsx, nhtsy);
/*stores estimated r0 during simulation */
    cn2est->r0=dnew(nwfspair,1);
    return cn2est;
}/*cn2est_prepare */

/**
   Embed gradent vector to gradient map 
*/
static void cn2est_embed(cn2est_t *cn2est, dcell *gradol, int icol){
    long *embed=cn2est->embed->p;
    for(int iwfs=0; iwfs<gradol->nx; iwfs++){
	if(!cn2est->wfscov[iwfs]) continue;
	const int nsa=cn2est->nsa;
	dmat *grad=gradol->p[iwfs];
	if(!grad){
	    error("wfs %d: PSOL grads is required to do cn2 estimation\n", iwfs);
	}
	if(icol>=grad->ny){
	    error("icol=%d is invalid\n", icol);
	}
	if(grad->nx!=nsa*2){
	    error("grad and saloc does not match\n");
	}
	double *pgrad=grad->p+grad->nx*icol;
	/*Embed gradients in a 2-d array */
	for(int isa=0; isa<nsa; isa++){
	    cn2est->gxs->p[iwfs]->p[embed[isa]]=pgrad[isa];
	    cn2est->gys->p[iwfs]->p[embed[isa]]=pgrad[isa+nsa];
	}
	/*Compute curvature of wavefront from gradients. */
	cmat*  cur=cn2est->curi->p[iwfs];
	dmat*  gx=cn2est->gxs->p[iwfs];
	dmat*  gy=cn2est->gys->p[iwfs];
	const int ny=cn2est->curi->p[iwfs]->ny;
	const int nx=cn2est->curi->p[iwfs]->nx;
	for(int iy=0; iy<ny; iy++){
	    for(int ix=0; ix<nx; ix++){
		if(IND(cn2est->mask,ix,iy)){
		    IND(cur,ix,iy)=IND(gx,ix+1,iy)+IND(gx,ix-1,iy)-2*IND(gx,ix,iy)/*gx along x */
			+IND(gy,ix,iy+1)+IND(gy,ix,iy-1)-2*IND(gy,ix,iy);/*gy along y */
		}else{
		    IND(cur,ix,iy)=0;//must set to zero.
		}
	    }
	}
	cfft2(cn2est->curi->p[iwfs], -1);
    }
}
/**
   Compute cross-covairance from gradient curvature
*/
static void cn2est_cov(cn2est_t *cn2est){
    /*accumulate cross-covariance of gradient curvature. */
    const int nwfspair=cn2est->nwfspair;
    cn2est->count++;
    /*for each wfs pair */
    for(int iwfspair=0; iwfspair<nwfspair; iwfspair++){
	const int wfs0=cn2est->pair[iwfspair].wfs0;
	const int wfs1=cn2est->pair[iwfspair].wfs1;
	const dcomplex *cur1=cn2est->curi->p[wfs0]->p;
	const dcomplex *cur2=cn2est->curi->p[wfs1]->p;
	dcomplex *cov=cn2est->covc->p[iwfspair]->p;
	for(long i=0; i<cn2est->nembed*cn2est->nembed; i++){
	    cov[i]+=conj(cur1[i])*(cur2[i]);
	}
    }/*iwfspair */
}
void cn2est_push(cn2est_t *cn2est, dcell *gradol){
    int ncol=0;
    if(gradol->nx<cn2est->nwfs){
	error("Grad has less number of wfs than required %d\n", cn2est->nwfs);
    }
    for(int iwfs=0; iwfs<cn2est->nwfs; iwfs++){
	if(!cn2est->wfscov[iwfs]) continue;
	if(ncol==0) {
	    ncol=gradol->p[iwfs]->ny;
	}else if(ncol!=gradol->p[iwfs]->ny){
	    error("different wfs has differnt number of columns\n");
	}
    }
    for(int icol=0; icol<ncol; icol++){
	cn2est_embed(cn2est, gradol, icol);
	cn2est_cov(cn2est);
    }
}
DEF_ENV_FLAG(CN2EST_NO_NEGATIVE, 1);

/**
   Do the Cn2 Estimation.
*/
void cn2est_est(cn2est_t *cn2est, int verbose, int reset){
    info2("cn2est from %d measurements\n", cn2est->count);
    cmat *covi=cnew(cn2est->nembed, cn2est->nembed);
#if COV_ROTATE
    dmat *covr=dnew(cn2est->nembed, cn2est->nembed);
#endif
    const int nwfspair=cn2est->wt->nx;
    for(int iwfspair=0; iwfspair<nwfspair; iwfspair++){
	cadd(&covi, 0, cn2est->covc->p[iwfspair], 1./(cn2est->count*cn2est->nembed*cn2est->nembed));
	cfft2(covi, 1);
	cfftshift(covi);
	creal2d(&cn2est->cov2->p[iwfspair], 0, covi, 1);
#if COV_ROTATE
	//roate and embed;
	dembed(covr, cn2est->cov2->p[iwfspair], -cn2est->pair[iwfspair].beta);
	double *cc=cn2est->cov1->p[iwfspair]->p;
	CN2PAIR_T *pair=cn2est->pair+iwfspair;
	int off=cn2est->nembed/2;
	for(long isep=pair->iht0; isep<pair->iht1; isep++){
	    cc[isep-pair->iht0]=IND(covr, off+isep, off)*IND(cn2est->overlapi,off+isep,off);
	}
#else
	for(long isep=0; isep<covi->nx*covi->ny; isep++){
	    cn2est->cov2->p[iwfspair]->p[isep]*=cn2est->overlapi->p[isep];//temporary. apply to mc instead.
	}
#endif
    }
    cfree(covi);
#if COV_ROTATE
    dfree(covr);
#endif
    dcellzero(cn2est->wt);
    dcellmm(&cn2est->wt, cn2est->iPnk, cn2est->cov1, "nn", 1);
  
    double wtsumsum=0;
    for(int iwfspair=0; iwfspair<nwfspair; iwfspair++){
	double wtsum=0;
	dmat *wt=cn2est->wt->p[iwfspair];//the layer weights. 
	dmat *ht=cn2est->ht->p[iwfspair];
	int nlayer=wt->nx;
	if(CN2EST_NO_NEGATIVE){
	    /*The following tries to remove negative weights.  For small
	      negatives, we just zero them. For significant ones, we remove the
	      layer from being estimated to avoid biasing the result.
	    */
	    dmat *Pnk=ddup(cn2est->Pnk->p[iwfspair+iwfspair*nwfspair]);
	    //number of negative layers found. 
	    int nfd=0,nfdi=0;
	    do{
		nfdi=0;
		for(int ix=0; ix<nlayer; ix++){
		    const double val=wt->p[ix];
		    if(val<0){
			if(val<-1e-2){
			    //negatives found. Remove columns in the forward matrix to disable this point.
			    //Pnk and redo the estimation matrix iPnk
			    nfd++;
			    nfdi++;
			    dmat *tmp=drefcols(Pnk, ix, 1);
			    dzero(tmp);
			    dfree(tmp);
			}else{
			    wt->p[ix]=0;
			}
		    }
		}
		if(nfdi>0){
		    warning_once("Ignore %d negative weights. set MAOS_CN2EST_NO_NEGATIVE=0 in your shell to disable this feature.\n", nfd);
		    dmat *iPnk2=dpinv(Pnk, NULL);
		    //compute the new result 
		    dmm(&wt, 0, iPnk2, cn2est->cov1->p[iwfspair], "nn", 1);
		    dfree(iPnk2);
		}
	    }while(nfdi);
	}
	for(int ix=0; ix<nlayer; ix++){
	    wtsum+=wt->p[ix];
	}
	wtsumsum+=wtsum;
	double r0=1.;
	if(wtsum>0){
	    r0=pow(wtsum,-3./5.);
	    cn2est->r0->p[iwfspair]=r0;
	}
	dscale(wt, 1./wtsum);
	if(verbose){
	    info2("r0=%.4fm theta0=%6f\" ",r0,calc_aniso(r0,wt->nx,ht->p,wt->p)*206265);
	    if(cn2est->dmht && cn2est->dmht->nx==2){
		info2("theta2=%6f\" ", calc_aniso2(r0,wt->nx,ht->p,wt->p,
						   cn2est->dmht->p[0], cn2est->dmht->p[1])*206265);
	    }
	    info2("wt=[");
	    for(int iht=0; iht<wt->nx; iht++){
		info2("%.4f ", wt->p[iht]);
	    }
	    info2("]\n");
	}
    }
    cn2est->r0m=pow(wtsumsum/cn2est->wt->nx, -3./5.);
    dcellzero(cn2est->wtrecon);
    dcellmm(&cn2est->wtrecon, cn2est->wtconvert, cn2est->wt, "nn", 1);
/*only 1 cell. norm to sum to 1. */
    normalize_sum(cn2est->wtrecon->p[0]->p, cn2est->wtrecon->p[0]->nx, 1);
    if(verbose){
	info2("r0m=%.4f theta0=%.4f\" ",cn2est->r0m, 
	      calc_aniso(cn2est->r0m,cn2est->wtrecon->p[0]->nx,
			 cn2est->htrecon->p,cn2est->wtrecon->p[0]->p)*206265);
	if(cn2est->dmht && cn2est->dmht->nx==2){
	    info2("theta2=%6f\" ", calc_aniso2(cn2est->r0m,cn2est->wtrecon->p[0]->nx,
					       cn2est->htrecon->p,cn2est->wtrecon->p[0]->p,
					       cn2est->dmht->p[0], cn2est->dmht->p[1])*206265);
	}
	info2("wt=[");
	for(int iht=0; iht<cn2est->wtrecon->p[0]->nx; iht++){
	    info2("%.4f ", cn2est->wtrecon->p[0]->p[iht]);
	}
	info2("]\n");
    }
/*divide by the number of accumulated frames. */
    if(reset){
	if(verbose) info2("reset the covariance");
	cn2est->count=0;/*reset the counter; */
	ccellzero(cn2est->covc);/*reset the numbers. */
    }
}

/**
   Free all the data.
*/
void cn2est_free(cn2est_t *cn2est){
    if(!cn2est) return;
    free(cn2est->pair);
    free(cn2est->wfscov);
    lfree(cn2est->embed);
    lfree(cn2est->mask);
    cellfree(cn2est->gxs);
    cellfree(cn2est->gys);
    cellfree(cn2est->curi);
    cellfree(cn2est->covc);
    cellfree(cn2est->cov1);
    cellfree(cn2est->cov2);
    dfree(cn2est->overlapi);
    cellfree(cn2est->Pnk);
    cellfree(cn2est->iPnk);
    cellfree(cn2est->wt);
    cellfree(cn2est->ht);
    dfree(cn2est->dmht);
    dfree(cn2est->r0);
    cellfree(cn2est->htrecon);
    cellfree(cn2est->os);
    cellfree(cn2est->dx);
    cellfree(cn2est->dmht);
    cellfree(cn2est->wtrecon);
    cellfree(cn2est->wtconvert);

    free(cn2est);
}

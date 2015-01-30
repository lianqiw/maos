/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
CN2EST_T *cn2est_new(dmat *wfspair, dmat *wfstheta, loc_t *saloc, dmat *saa, const double saat,
		     const dmat *hs, dmat *htrecon, int keepht, double L0){
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
    CN2EST_T *cn2est=calloc(1, sizeof(CN2EST_T));
    cn2est->nwfspair=nwfspair;
    int nwfs=wfstheta->nx;
    cn2est->nwfs=nwfs;
    cn2est->nsa=saloc->nloc;
    /*wfscov is a flag for wfs show whether wfs participates in covariance. */
    cn2est->wfscov=calloc(nwfs, sizeof(int));
    for(int ind=0; ind<npair; ind++){
	int iwfs=(int)wfspair->p[ind];
	if(iwfs<0 || iwfs>=nwfs){
	    error("wfspair has invalid values: %d\n", iwfs);
	}
	cn2est->wfscov[iwfs]=1;
    }
    /*embed has dimension nembed*nembed. It is used to embed 1-d gradient vector
      to a 2-d map for each curvature and covariance compuation later.*/
    cn2est->nembed=0;
    cn2est->embed=loc_create_embed(&cn2est->nembed, saloc, 2, 0);
    const int nx=cn2est->nembed;
    /*mask is a mask for valid subapertures that have enough area */
    int *mask=calloc(nx*nx,sizeof(int));
    if(saa && saa->nx!=cn2est->nsa){
	error("saa and saloc mismatch\n");
    }
    
    double saat2=dmax(saa)*saat;
    for(int isa=0; isa<cn2est->nsa; isa++){
	if(!saa || saa->p[isa]>saat2){
	    mask[cn2est->embed->p[isa]]=1;/*use this subaperture */
	}
    }
    int (*pmask)[nx]=(void*)mask;
    cn2est->mask=lnew(nx,nx);
    int (*pmask2)[nx]=(void*)cn2est->mask->p;
    cmat *overlap=cnew(nx, nx);
    //cfft2plan(overlap, -1);
    //cfft2plan(overlap, 1);
    PCMAT(overlap, pover);
    int iymin=nx,iymax=0,ixmin=nx,ixmax=0;
    for(int iy=0; iy<nx; iy++){
	for(int ix=0; ix<nx; ix++){
	    /*Only use a subaperture if we are able to make curvature. */
	    if(pmask[iy][ix] && pmask[iy][ix+1] && pmask[iy][ix-1] &&
	       pmask[iy+1][ix] && pmask[iy-1][ix]){
		pmask2[iy][ix]=1;
		pover[iy][ix]=1;
		if(ix>ixmax) ixmax=ix;
		if(ix<ixmin) ixmin=ix;
		if(iy>iymax) iymax=iy;
		if(iy<iymin) iymin=iy;
	    }
	}
    }
    int maxsep=MIN((iymax-iymin), (ixmax-ixmin));
    free(mask);
    cfft2(overlap, -1);
    for(long i=0; i<nx*nx; i++){
	overlap->p[i]=overlap->p[i]*conj(overlap->p[i]);
    }
    cfft2(overlap, 1);
    cfftshift(overlap);
    creal2d(&cn2est->overlap, 0, overlap, 1./(nx*nx));
    cfree(overlap);
    /*2-d arrays to store x y gradient  and x y "curvature" */
    cn2est->gxs=cellnew(nwfs, 1);/*stores gradient in 2-d map */
    cn2est->gys=cellnew(nwfs, 1);
    cn2est->curi=cellnew(nwfs, 1);
    for(int ix=0; ix<nwfs; ix++){
	if(cn2est->wfscov[ix]){
	    cn2est->gxs->p[ix]=dnew(nx, nx);
	    cn2est->gys->p[ix]=dnew(nx, nx);
	    cn2est->curi->p[ix]=cnew(nx, nx);
	    //cfft2plan(cn2est->curi->p[ix], -1);
	    //cfft2plan(cn2est->curi->p[ix], 1);
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
    PDMAT(wfstheta, pwfstheta);
    double hmin, hmax;
    dmaxmin(cn2est->htrecon->p, cn2est->htrecon->nx, &hmax, &hmin);
    for(int iwfspair=0; iwfspair<nwfspair; iwfspair++){
	/*get pointer for this pair */
	CN2PAIR_T *pair=cn2est->pair+iwfspair;
	/*The WFS in this pair. */
	int wfs0=wfspair->p[iwfspair*2];
	int wfs1=wfspair->p[iwfspair*2+1];
	if(wfs0==wfs1){
	    error("must pair different wfs.\n");
	}
	pair->wfs0=wfs0;
	pair->wfs1=wfs1;
	/*The separation between the stars */
	double dthetax=pwfstheta[0][wfs0]-pwfstheta[0][wfs1];
	double dthetay=pwfstheta[1][wfs0]-pwfstheta[1][wfs1];
	/*The direction of the WFS pair baseline vector */
	double ang=atan2(dthetay, dthetax);
	pair->beta=ang;
	/*the angular distance between WFS */
	pair->dtheta=sqrt(dthetax*dthetax+dthetay*dthetay);
	/*The average wfs guide star height of the two*/
	if(hs->nx*hs->ny>=nwfs){
	    pair->hsm=(hs->p[wfs0]+hs->p[wfs1])*0.5;
	}else if(hs->nx*hs->ny==1){
	    pair->hsm=hs->p[0];
	}else{
	    error("hs in in wrong format, should have 1 or %d elements\n", nwfs);
	}
	/*The number of layers is determined */
	if(keepht==2){
	    /*slodar output directly to layers used for tomography. may not work well */
	    pair->iht0=0;
	    pair->iht1=cn2est->htrecon->nx;
	}else{
	    /*output to nature slodar heights. and then bin to tomography layers */
	    pair->iht0=(int)floor(hmin*pair->dtheta/(dsa*(1-hmin/pair->hsm)));
	    pair->iht1=(int)floor(hmax*pair->dtheta/(dsa*(1-hmax/pair->hsm)));
	    if(pair->iht0<-maxsep){
		pair->iht0=-maxsep;
	    }
	    if(pair->iht1>maxsep){
		pair->iht1=maxsep;
	    }
	    pair->iht1++;//exclusive
	}
	pair->nsep=(pair->iht1-pair->iht0);
	nhtsx[iwfspair]=pair->nsep;
	nhtsy[iwfspair]=1;
	/*initialize ht array to store layer heights */
	cn2est->ht->p[iwfspair]=dnew(pair->nsep,1);
	/*stores the cross-covariance of the curvature for x, and y grads. */
	cn2est->cov1->p[iwfspair]=dnew(pair->nsep,1);
	cn2est->covc->p[iwfspair]=cnew(nx, nx);

	/*info2("Pair %d: wfs %d and %d. dtheta=%4f\" iht=[%d, %d)\n", 
	  iwfspair, wfs0, wfs1, pair->dtheta*206265, pair->iht0, pair->iht1);*/
    }/*iwfspair */
    /*stores estimated weight of layers during simulation and output to file finally. */
    cn2est->wt=dcellnew3(nwfspair, 1, nhtsx, nhtsy);
    /*stores estimated r0 during simulation */
    cn2est->r0=dnew(nwfspair,1);
    /*
      Now we start to build the model that will be used to estimate Cn2 from the
      cross-covariance.
    */

    /*ovs is the over sampling factor in mxx, myy. need to be at least 2 to
      cover the non-zero PSD (sinc functions.)*/
      
      /*deprecated:for ovs bigger than 1, the
      results changes even with wfs pairs aligned along x/y direction. The
      reason is that we are sampling more PSDs with bigger ovs. This is probably
      the aliasing effect. With larger ovs, we are using larger
      frequencies. That gets aliased into the results. Then how about padding?*/
    cn2est->ovs=2;
    /*Pnk stores the opeator from layer weight to gradient curvature covariance
      matrix. iPnk is the covariance of it. Each diagonal cell is for each WFS
      pair.*/
    cn2est->Pnk=cellnew(nwfspair, nwfspair);
    cn2est->iPnk=cellnew(nwfspair,nwfspair);
    /*wtconvert is the matrix to down/up sample the CN2 estimates to layers
      used for tomography*/
    cn2est->wtconvert=cellnew(1,nwfspair);
    cn2est->L0=L0;
    for(int iwfspair=0; iwfspair<nwfspair; iwfspair++){
	CN2PAIR_T *pair=cn2est->pair+iwfspair;
	const int nsep=pair->nsep;
	/*nm is the size of the mxx myy array */
	const int nm=4*(MAX(abs(pair->iht0), abs(pair->iht1-1))+1)*cn2est->ovs;
	const int nm2=nm>>1;
	/*sampling in mxx, myy after FFT */
	const double dx=dsa/cn2est->ovs;
	const double df=1./(nm*dx);
	const double L02=pow(cn2est->L0,-2);
	/*initialize */
	cmat *mc=cnew(nm,nm);
	/*create 2-d pointers */
	PCMAT(mc,pmc);
	/*create FFT plans */
	//cfft2plan(mc,1);
	/*the forward operator from layer weights to cross-covariance */
 	dmat *Pnk=dnew(nsep, nsep);
	PDMAT(Pnk, pPnk);
	info2("Pair %d: hk=[", iwfspair);
	double cb=cos(pair->beta);
	double sb=sin(pair->beta);
	for(int iht=pair->iht0; iht<pair->iht1; iht++){
	    /*the height of each layer */
	     double hk;
	    if(keepht==2){
		hk=cn2est->htrecon->p[iht-pair->iht0];
	    }else{
		hk=dsa*iht/(pair->dtheta+dsa*iht/pair->hsm);
	    }
	    info2("%.2f ",hk*0.001);
	    cn2est->ht->p[iwfspair]->p[iht-pair->iht0]=hk;
	    /*the cone effect */
	    const double zeta=1.-hk/pair->hsm;
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
		    pmc[iy][ix]=pow(sincfy*sincfx,2)*psd*cur;
		}/*ix */
	    }/*iy */
	    /*doing fft */
	    cfft2(mc,1);/*fft */
	    /*shift 0 freq to center. */
	    cfftshift(mc);
	    for(int isep=pair->iht0; isep<pair->iht1; isep++){
		/*Different from matlab prototype because I used interpolation. */
		double xx=(iht-isep)*cb*cn2est->ovs + nm2;
		double yy=(iht-isep)*sb*cn2est->ovs + nm2;
#if INTERP_NEAREST
		/*Do interpolation using nearest neighbor */
		int ixx=(int)round(xx);
		int iyy=(int)round(yy);
		double imc=creal(pmc[iyy][ixx]);
#else
		/*Do interpolation using bilinear spline interp. */
		int ixx=(int)floor(xx); xx=xx-ixx;
		int iyy=(int)floor(yy); yy=yy-iyy;
		double imc=creal((pmc[iyy][ixx]*(1-xx)+pmc[iyy][ixx+1]*(xx))*(1-yy)
				 +(pmc[iyy+1][ixx]*(1-xx)+pmc[iyy+1][ixx+1]*(xx))*yy);
#endif
		pPnk[iht-pair->iht0][isep-pair->iht0]=imc;
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
    return cn2est;
}/*cn2est_prepare */

/**
   Embed gradent vector to gradient map 
*/
static void cn2est_embed(CN2EST_T *cn2est, dcell *gradol, int icol){
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
	PCMAT(cn2est->curi->p[iwfs], cur);
	PDMAT(cn2est->gxs->p[iwfs], gx);
	PDMAT(cn2est->gys->p[iwfs], gy);
	int (*mask)[cn2est->nembed]=(void*)cn2est->mask->p;
	const int ny=cn2est->curi->p[iwfs]->ny;
	const int nx=cn2est->curi->p[iwfs]->nx;
	for(int iy=0; iy<ny; iy++){
	    for(int ix=0; ix<nx; ix++){
		if(mask[iy][ix]){
		    cur[iy][ix]=gx[iy][ix+1]+gx[iy][ix-1]-2*gx[iy][ix]/*gx along x */
			+gy[iy+1][ix]+gy[iy-1][ix]-2*gy[iy][ix];/*gy along y */
		}else{
		    cur[iy][ix]=0;//must set to zero.
		}
	    }
	}
	cfft2(cn2est->curi->p[iwfs], -1);
    }
}
/**
   Compute cross-covairance from gradient curvature
*/
static void cn2est_cov(CN2EST_T *cn2est){
    /*accumulate cross-covariance of gradient curvature. */
    const int nwfspair=cn2est->nwfspair;
    cn2est->count++;
    /*for each wfs pair */
    for(int iwfspair=0; iwfspair<nwfspair; iwfspair++){
	const int wfs0=cn2est->pair[iwfspair].wfs0;
	const int wfs1=cn2est->pair[iwfspair].wfs1;
	dcomplex *cur1=cn2est->curi->p[wfs0]->p;
	dcomplex *cur2=cn2est->curi->p[wfs1]->p;
	dcomplex *cov=cn2est->covc->p[iwfspair]->p;
	for(long i=0; i<cn2est->nembed*cn2est->nembed; i++){
	    cov[i]+=cur1[i]*conj(cur2[i]);
	}
    }/*iwfspair */
}
void cn2est_push(CN2EST_T *cn2est, dcell *gradol){
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

/**
   Do the Cn2 Estimation.
 */
void cn2est_est(CN2EST_T *cn2est, int verbose, int reset){
    info2("cn2est from %d measurements\n", cn2est->count);
    cmat *covi=cnew(cn2est->nembed, cn2est->nembed);
    dmat *covr=dnew(cn2est->nembed, cn2est->nembed);
    PDMAT(covr, pcovr);
    PDMAT(cn2est->overlap, pover);
    const int nwfspair=cn2est->wt->nx;
    //cfft2plan(covi, 1);
    for(int iwfspair=0; iwfspair<nwfspair; iwfspair++){
	//cfft2plan(covi, 1);
	cadd(&covi, 0, cn2est->covc->p[iwfspair], 1./(cn2est->count*cn2est->nembed*cn2est->nembed));
	cfft2(covi, 1);
	cfftshift(covi);
	creal2d(&cn2est->cov2->p[iwfspair], 0, covi, 1);
	//roate and embed;
	dembed(covr, cn2est->cov2->p[iwfspair], M_PI-cn2est->pair[iwfspair].beta);
	double *cc=cn2est->cov1->p[iwfspair]->p;
	CN2PAIR_T *pair=cn2est->pair+iwfspair;
	int off=cn2est->nembed/2;
	for(int isep=pair->iht0; isep<pair->iht1; isep++){
	    if(pover[off][off+isep]>=MIN_SA_OVERLAP){
		cc[isep-pair->iht0]=(double)pcovr[off][off+isep]/pover[off][off+isep];
	    }else{
		cc[isep-pair->iht0]=0;
	    }
	}
    }
    cfree(covi);
    dfree(covr);
    dcellzero(cn2est->wt);
    dcellmm(&cn2est->wt, cn2est->iPnk, cn2est->cov1, "nn", 1);
    double wtsumsum=0;
    DEF_ENV_FLAG(CN2EST_NO_NEGATIVE, 1);
    for(int iwfspair=0; iwfspair<nwfspair; iwfspair++){
	double wtsum=0;
	dmat *wt=cn2est->wt->p[iwfspair];//the layer weights. 
	dmat *ht=cn2est->ht->p[iwfspair];
	int nlayer=wt->nx;
	if(CN2EST_NO_NEGATIVE){
	    /*The following tries to remove negative heights. 
	      This should not be done because noise causing small negative and postive values. 
	      If we simply remove negative values, it will bias the r0.*/
	    dmat *Pnk=cn2est->Pnk->p[iwfspair+iwfspair*nwfspair];
	    dmat *Pnkwt=dnew(nlayer, 1);//to remove columns in Pnk if negative weights are returned 
	    //number of negative layers found. 
	    int nfd=0;
	    dset(Pnkwt, 1);
	  repeat:
	    nfd=0;
	    for(int ix=0; ix<nlayer; ix++){
		const double val=wt->p[ix];
		if(val<0 && Pnkwt->p[ix]>0){
		    //negatives found. Remove columns in the forward matrix to disable this point.
		    //Pnk and redo the estimation matrix iPnk
		    nfd++;
		    //disable diagonal elements. 
		    Pnkwt->p[ix]=0;
		}
	    }
	    if(nfd>0){
		if(nfd) warning_once("Ignore %d negative points. set MAOS_CN2EST_NO_NEGATIVE=0 to disable.\n", nfd);

		dmat *iPnk2=dpinv(Pnk, Pnkwt);
		//compute the new result 
		dmm(&wt, 0, iPnk2, cn2est->cov1->p[iwfspair], "nn", 1);
		for(int ix=0; ix<nlayer; ix++){
		    if(Pnkwt->p[ix]==0){
			wt->p[ix]=0;
		    }
		}
		dfree(iPnk2);
		goto repeat;
	    }
	    dfree(Pnkwt);
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
void cn2est_free(CN2EST_T *cn2est){
    if(!cn2est) return;
    lfree(cn2est->embed);
    lfree(cn2est->mask);
    free(cn2est->wfscov);
    free(cn2est->pair);
    dcellfree(cn2est->gxs);
    dcellfree(cn2est->gys);
    ccellfree(cn2est->curi);
    ccellfree(cn2est->covc);
    dcellfree(cn2est->cov1);
    dcellfree(cn2est->cov2);
    dcellfree(cn2est->iPnk);
    dcellfree(cn2est->wt);
    dcellfree(cn2est->ht);
    dfree(cn2est->dmht);
    dfree(cn2est->r0);
    free(cn2est);
}

/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "maos.h"
#include "cn2est.h"
#include "setup_recon.h"


/**
   Prepare arrays for cn2 estimation. Multiple pairs can be used to do Cn2
Estimation. The result will be an average of them.  */
#define INTERP_NEAREST 0 /*set to 0 after debugging */

CN2EST_T *cn2est_prepare(const PARMS_T *parms, const POWFS_T *powfs){
    info2("Cn2 estimation:");
    /*We need at least a pair */
    if(!parms->cn2.pair || parms->cn2.npair==0) return NULL;
    if(parms->cn2.npair%2==1){
	error("parms->cn2.pair must have even number of entries\n");
    }
    /*>>1 is a short cut for /2 */
    int nwfspair=parms->cn2.npair>>1;
    CN2EST_T *cn2est=calloc(1, sizeof(CN2EST_T));
    cn2est->nwfspair=nwfspair;
    cn2est->reset=parms->cn2.reset;
    /*wfscov is a flag for wfs show whether wfs participates in covariance. */
    cn2est->wfscov=calloc(parms->nwfs, sizeof(int));
    int ipowfs=-1;
    for(int ind=0; ind<parms->cn2.npair; ind++){
	int iwfs=parms->cn2.pair[ind];
	cn2est->wfscov[iwfs]=1;
	if(ipowfs==-1){
	    ipowfs=parms->wfs[iwfs].powfs;
	}else if(ipowfs!=parms->wfs[iwfs].powfs){
	    error("All wfs in parms->cn2.pair do not belong to the same powfs\n");
	}
    }
    cn2est->ipowfs=ipowfs;
    /*embed has dimension nembed*nembed. It is used to embed 1-d gradient vector
      to a 2-d map for each curvature and covariance compuation later.*/
    cn2est->nembed=parms->powfs[ipowfs].order*2;
    cn2est->embed=loc_create_embed(&cn2est->nembed, powfs[ipowfs].saloc);
    const int nx=cn2est->nembed;
    /*mask is a mask for valid subapertures that have enough area */
    int *mask=calloc(nx*nx,sizeof(int));
    for(int isa=0; isa<powfs[ipowfs].saloc->nloc; isa++){
	if(powfs[ipowfs].saa->p[isa]>parms->cn2.saat){
	    mask[cn2est->embed[isa]]=1;/*use this subaperture */
	}
    }
    int (*pmask)[nx]=(void*)mask;
    int (*pmask2)[nx]=calloc(nx*nx,sizeof(int));
    for(int iy=0; iy<nx; iy++){
	for(int ix=0; ix<nx; ix++){
	    /*Only use a subaperture if we are able to make curvature. */
	    if(pmask[iy][ix] && pmask[iy][ix+1] && pmask[iy][ix-1] &&
	       pmask[iy+1][ix] && pmask[iy-1][ix]){
		pmask2[iy][ix]=1;
	    }
	}
    }
    free(mask);
    /*2-d arrays to store x y gradient  and x y "curvature" */
    cn2est->gxs=dcellnew(parms->nwfs, 1);/*stores gradient in 2-d map */
    cn2est->gys=dcellnew(parms->nwfs, 1);
    cn2est->cxs=dcellnew(parms->nwfs, 1);/*stores curvature in 2-d map. */
    cn2est->cys=dcellnew(parms->nwfs, 1);
    for(int ix=0; ix<parms->nwfs; ix++){
	if(cn2est->wfscov[ix]){
	    cn2est->gxs->p[ix]=dnew(nx, nx);
	    cn2est->gys->p[ix]=dnew(nx, nx);
	    cn2est->cxs->p[ix]=dnew(nx, nx);
	    cn2est->cys->p[ix]=dnew(nx, nx);
	}
    }
    /*
      Now we prepare indexes that can be used to compute gradient(curvature)
      cross-covariances conveniently during simulation.
     */
    /*first get a few constants */
    const double hmax=parms->cn2.hmax;/*maximum height to estimate. */
    const double dsa=powfs[ipowfs].pts->dsa;
    const double hs=parms->powfs[ipowfs].hs;
    const int nsa=powfs[ipowfs].pts->nsa;
    /*determine the layer height used for tomography. */
    if(parms->cn2.keepht){/*use atmr.ht */
	const int nht=parms->atmr.nps;
	cn2est->htrecon=dnew(nht,1);
	cn2est->os=dnew(nht,1);
	for(int iht=0; iht<nht; iht++){
	    cn2est->htrecon->p[iht]=parms->atmr.ht[iht];
	    cn2est->os->p[iht]=parms->atmr.os[iht];
	}
    }else{/*use linearly spaced number of layers between ground and hmax */
	int nht=parms->cn2.nhtomo;
	cn2est->htrecon=dnew(nht,1);
	cn2est->os=dnew(nht,1);
	if(nht<1) error("invalid nhtomo");
	/*calculate the number of over sampled layers. */
	int osc=0;
	for(int ips=0; ips<parms->atmr.nps; ips++){
	    if(parms->atmr.os[ips]>1){
		osc++;
		if(parms->atmr.os[ips]!=2){
		    error("os is not 2. adept this code to it.\n");
		}
	    }
	}
	double dht=parms->cn2.hmax/(double)(nht-1);
	for(int iht=0; iht<nht; iht++){
	    cn2est->htrecon->p[iht]=dht*iht;
	    if(iht<osc){
		cn2est->os->p[iht]=2;
	    }else{
		cn2est->os->p[iht]=1;
	    }
	}
    }
    cn2est->wtrecon=dcellnew(1,1);
    cn2est->wtrecon->p[0]=dnew_mmap(cn2est->htrecon->nx,1,NULL,"Res_Cn2_wtrecon");
    dwrite(cn2est->htrecon,"Res_Cn2_htrecon");
    {
	info2("htrecon=[");
	for(int iht=0; iht<cn2est->htrecon->nx; iht++){
	    info2("%.2f ", cn2est->htrecon->p[iht]*0.001);
	}
	info2("]km\n");
    }
    /*stores cross-covariance data during simulation */
    cn2est->cc=dcellnew(nwfspair,1);
    /*height of layers for each wfs pair of slodar output. */
    cn2est->ht=dcellnew(nwfspair,1);
    /*record sapair to use for each separation */
    cn2est->pair=calloc(nwfspair, sizeof(CN2PAIR_T));
    long nhtsx[nwfspair]; 
    long nhtsy[nwfspair];
    for(int iwfspair=0; iwfspair<nwfspair; iwfspair++){
	/*get pointer for this pair */
	CN2PAIR_T *pair=cn2est->pair+iwfspair;
	/*The WFS in this pair. */
	int wfs0=parms->cn2.pair[iwfspair*2];
	int wfs1=parms->cn2.pair[iwfspair*2+1];
	if(wfs0==wfs1){
	    error("must pair different wfs.\n");
	}
	pair->wfs0=wfs0;
	pair->wfs1=wfs1;
	/*The separation between the stars */
	double dthetax=parms->wfs[wfs0].thetax-parms->wfs[wfs1].thetax;
	double dthetay=parms->wfs[wfs0].thetay-parms->wfs[wfs1].thetay;
	/*The direction of the WFS pair baseline vector */
	double ang=atan2(dthetay, dthetax);
	pair->beta=ang;
	/*the angular distance between WFS */
	pair->dtheta=sqrt(dthetax*dthetax+dthetay*dthetay);
	/*The number of layers is determined */
	int nht;
	if(parms->cn2.keepht!=2){/*output to nature slodar heights.  */
	    nht=(int)ceil(hmax*pair->dtheta/(dsa*(1-hmax/hs)))+1;
	}else{ /*slodar output directly to layers used for tomography. may not work well */
	    nht=cn2est->htrecon->nx;
	}
	pair->nht=nht;
	nhtsx[iwfspair]=nht;
	nhtsy[iwfspair]=1;
	/*initialize ht array to store layer heights */
	cn2est->ht->p[iwfspair]=dnew(nht,1);
	/*number of subaperture separations to choose. */
	const int nsep=pair->nsep=nht;
	/*round the WFS pair baseline vector to along x or y axis so covariances
	  along that direction can be easily computed.  We choose the
	  subaperture separation direction as close to the star separation as
	  possible, which maximizes sensitivity. Can further improve by
	  considering the directions along PI/4 where we can do xstep=1 and
	  ystep=1;*/
	double ang2=round(ang*2/M_PI)*M_PI/2;
	/*either xstep or ystep is 0. the other one is 1 or -1 */
	int xstep=(int)round(cos(ang2));/*integers */
	int ystep=(int)round(sin(ang2));/*integers */
	pair->xstep=xstep;
	pair->ystep=ystep;
	/*records the subaperture pairs to use for each subaperture
	  separation. Used in wfsgrad.c.*/
	pair->sapair  = calloc(nsep, sizeof(void*));
	pair->nsapair = calloc(nsep, sizeof(int));
	/*stores the cross-covariance of the curvature for x, and y grads. */
	cn2est->cc->p[iwfspair]=dnew(nsep*2,1);
	info2("Pair %d: wfs %d and %d. dtheta=%4f\" nht=%d xstep=%d ystep=%d\n", 
	     iwfspair, wfs0, wfs1, pair->dtheta*206265, nht, xstep, ystep);
	for(int isep=0; isep<nsep; isep++){
	    pair->sapair[isep]=calloc(nsa*2, sizeof(int));
	    int count=0;
	    for(int iy=0; iy<nx; iy++){
		for(int ix=0; ix<nx; ix++){
		    int xsep=isep*xstep;
		    int ysep=isep*ystep;
		    /*if both subapertures exist for the pair, record it. */
		    if(pmask2[iy][ix] && pmask2[iy+ysep][ix+xsep]){
			pair->sapair[isep][count][0]=ix+iy*nx;
			pair->sapair[isep][count][1]=(ix+xsep)+(iy+ysep)*nx;
			count++;
		    }/*mask */
		}/*ix */
	    }/*iy */
	    pair->sapair[isep]=realloc(pair->sapair[isep], sizeof(int)*2*count);
	    if(count==0){
		error("there is no overlapping for this sep: %d\n", isep);
	    }
	    pair->nsapair[isep]=count;
	}/*isep */
    }/*iwfspair */
    free(pmask2);
    /*stores estimated weight of layers during simulation and output to file finally. */
    cn2est->wt=dcellnew_mmap(nwfspair, 1, nhtsx, nhtsy, NULL, NULL, "Res_Cn2_wt");
    /*stores estimated r0 during simulation */
    cn2est->r0=dnew_mmap(nwfspair,1,NULL,"Res_Cn2_r0");
    dcellwrite(cn2est->ht, "Res_Cn2_ht");
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
    /*Pkn stores the opeator from layer weight to gradient curvature covariance
      matrix. iPkn is the covariance of it. Each diagonal cell is for each WFS
      pair.*/
    cn2est->Pkn=dcellnew(nwfspair, nwfspair);
    cn2est->iPkn=dcellnew(nwfspair,nwfspair);
    /*wtconvert is the matrix to down/up sample the CN2 estimates to layers
      used for tomography*/
    cn2est->wtconvert=spcellnew(1,nwfspair);
    cn2est->l0=parms->atmr.l0;
    for(int iwfspair=0; iwfspair<nwfspair; iwfspair++){
	CN2PAIR_T *pair=cn2est->pair+iwfspair;
	const int nsep=pair->nsep;
	const int nht=pair->nht;
	/*nm is the size of the mxx myy array */
	const int nm=2*nsep*cn2est->ovs;
	const int nm2=nm>>1;
	/*sampling in mxx, myy after FFT */
	const double dx=dsa/cn2est->ovs;
	const double df=1./(nm*dx);
	const double l02=pow(cn2est->l0,-2);
	/*initialize */
	cmat *mxx=cnew(nm,nm);
	cmat *myy=cnew(nm,nm);
	/*create 2-d pointers */
	PCMAT(mxx,pmxx);
	PCMAT(myy,pmyy);
	/*create FFT plans */
	cfft2plan(mxx,1);
	cfft2plan(myy,1);
	/*the forward operator from layer weights to cross-covariance */
 	dmat *Pkn=dnew(nsep*2,nht);
	PDMAT(Pkn, pPkn);
	info2("Pair %d: hk=[", iwfspair);
	for(int iht=0; iht<nht; iht++){
	    /*the height of each layer */
	     double hk;
	    if(parms->cn2.keepht!=2){
		hk=dsa*iht/(pair->dtheta+dsa*iht/hs);
	    }else{
		hk=cn2est->htrecon->p[iht];
	    }
	    info2("%.2f ",hk*0.001);
	    cn2est->ht->p[iwfspair]->p[iht]=hk;
	    /*the cone effect */
	    const double zeta=1.-hk/hs;
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
		const double fy=(iy<nm2?iy:iy-nm)*df;
		/*the sinc subaperture fnction */
		const double sincfy=sinc(fy*dsa);
		for(int ix=0; ix<nm; ix++){
		    const double fx=(ix<nm2?ix:ix-nm)*df;
		    const double sincfx=sinc(fx*dsa);
		    /*the turbulence PSD with outerscale */
		    const double psd=psd_coef*pow((fx*fx+fy*fy)*zetan2+l02,-11./6.);
		    /*gx diff is along x, gy diff is along y to form real curvature */
		    const double curx=pow(2*(cos(2*M_PI*dsa*fx)-1),2);
		    const double cury=pow(2*(cos(2*M_PI*dsa*fy)-1),2);
		    const double pmmm=pow(sincfy*sincfx,2)*psd;
		    pmxx[iy][ix]=fx*fx*pmmm*curx;/*x gradient */
		    pmyy[iy][ix]=fy*fy*pmmm*cury;/*y gradient */
		}/*ix */
	    }/*iy */
	    /*doing fft */
	    cfft2(mxx,1);/*fft */
	    cfft2(myy,1);
	    /*shift 0 freq to center. */
	    cfftshift(mxx);
	    cfftshift(myy);
	    
	    for(int isep=0; isep<nsep; isep++){
		/*Different from matlab prototype because I used interpolation. */
		double xx=(iht*cos(pair->beta)-isep*pair->xstep)*cn2est->ovs + nm2;
		double yy=(iht*sin(pair->beta)-isep*pair->ystep)*cn2est->ovs + nm2;
#if INTERP_NEAREST
		/*Do interpolation using nearest neighbor */
		int ixx=(int)round(xx);
		int iyy=(int)round(yy);
		double imxx=creal(pmxx[iyy][ixx]);
		double imyy=creal(pmyy[iyy][ixx]);
#else
		/*Do interpolation using bilinear spline interp. */
		int ixx=(int)floor(xx); xx=xx-ixx;
		int iyy=(int)floor(yy); yy=yy-iyy;
		double imxx=creal((pmxx[iyy][ixx]*(1-xx)+pmxx[iyy][ixx+1]*(xx))*(1-yy)
				  +(pmxx[iyy+1][ixx]*(1-xx)+pmxx[iyy+1][ixx+1]*(xx))*yy);
		double imyy=creal((pmyy[iyy][ixx]*(1-xx)+pmyy[iyy][ixx+1]*(xx))*(1-yy)
				  +(pmyy[iyy+1][ixx]*(1-xx)+pmyy[iyy+1][ixx+1]*(xx))*yy);
#endif
		pPkn[iht][isep]=imxx;
		pPkn[iht][isep+nsep]=imyy;
	    }
	}
	info2("]km\n");
	/*
	  iPkn is a block diagonal matrix for Cn2 Estimation.
	*/
	cn2est->Pkn->p[iwfspair+iwfspair*nwfspair]=dref(Pkn);
	cn2est->iPkn->p[iwfspair+iwfspair*nwfspair]=dpinv(Pkn,NULL,NULL);
	cn2est->wtconvert->p[iwfspair]=mkhbin1d(cn2est->ht->p[iwfspair],cn2est->htrecon);
	dfree(Pkn);
	cfree(mxx);
	cfree(myy);
    }/*iwfspair */
    if(parms->save.setup){
	dcellwrite(cn2est->iPkn,"%s/cn2_iPkn",dirsetup);
	dcellwrite(cn2est->Pkn,"%s/cn2_Pkn",dirsetup);
	dcellwrite(cn2est->ht,"%s/cn2_ht",dirsetup);
	spcellwrite(cn2est->wtconvert,"%s/cn2_wtconvert",dirsetup);
    }
    return cn2est;
}/*cn2est_prepare */

/**
   Embed gradent vector to gradient map.
*/
void cn2est_embed(CN2EST_T *cn2est, dmat *grad, int iwfs){
    if(!grad){
	error("wfs %d: PSOL grads is required to do cn2 estimation\n", iwfs);
    }
    long *embed=cn2est->embed;
    const int nsa=grad->nx>>1;
    /*Embed gradients in a 2-d array */
    for(int isa=0; isa<nsa; isa++){
	cn2est->gxs->p[iwfs]->p[embed[isa]]=grad->p[isa];
	cn2est->gys->p[iwfs]->p[embed[isa]]=grad->p[isa+nsa];
    }
    /*Compute curvature of wavefront from gradients. */
    PDMAT(cn2est->cxs->p[iwfs], curx);
    PDMAT(cn2est->cys->p[iwfs], cury);
    PDMAT(cn2est->gxs->p[iwfs], gx);
    PDMAT(cn2est->gys->p[iwfs], gy);
    const int ny=cn2est->cxs->p[iwfs]->ny;
    const int nx=cn2est->cxs->p[iwfs]->nx;
    for(int iy=0; iy<ny; iy++){
	for(int ix=0; ix<nx; ix++){
	    curx[iy][ix]=gx[iy][ix+1]+gx[iy][ix-1]-2*gx[iy][ix];/*gx along x */
	    cury[iy][ix]=gy[iy+1][ix]+gy[iy-1][ix]-2*gy[iy][ix];/*gy along y */
	}
    }
}

/**
   Compute cross-covairance from gradient curvature
*/
void cn2est_cov(CN2EST_T *cn2est){
    /*accumulate cross-covariance of gradient curvature. */
    const int nwfspair=cn2est->nwfspair;
    cn2est->nstep++;
    /*for each wfs pair */
    for(int iwfspair=0; iwfspair<nwfspair; iwfspair++){
	const int nsep=cn2est->pair[iwfspair].nsep;
	double *ccx=cn2est->cc->p[iwfspair]->p;
	double *ccy=cn2est->cc->p[iwfspair]->p+nsep;
	const int wfs0=cn2est->pair[iwfspair].wfs0;
	const int wfs1=cn2est->pair[iwfspair].wfs1;
	const double *cxs1=cn2est->cxs->p[wfs0]->p;
	const double *cxs2=cn2est->cxs->p[wfs1]->p;
	const double *cys1=cn2est->cys->p[wfs0]->p;
	const double *cys2=cn2est->cys->p[wfs1]->p;
	
	/*for each subaperture separation */
	for(int isep=0; isep<nsep; isep++){
	    const int nsapair=cn2est->pair[iwfspair].nsapair[isep];
	    int (*const sapair)[2]=cn2est->pair[iwfspair].sapair[isep];
	    double accx=0;
	    double accy=0;
	    /*for each pair of the subapertures */
	    for(int isa=0; isa<nsapair; isa++){
		accx+=cxs1[sapair[isa][0]]*cxs2[sapair[isa][1]];
		accy+=cys1[sapair[isa][0]]*cys2[sapair[isa][1]];
	    }/*isa */
	    ccx[isep]+=accx/nsapair;
	    ccy[isep]+=accy/nsapair;
	}/*isep */
    }/*iwfspair */
}
/**
   Do the Cn2 Estimation.
 */
void cn2est_est(CN2EST_T *cn2est, const PARMS_T *parms){
    if(parms->cn2.verbose){
	info2("Cn2 is estimated with %d averages\n", cn2est->nstep);
    }
    dcellzero(cn2est->wt);
    dcellmm(&cn2est->wt, cn2est->iPkn, cn2est->cc, "nn", 1./cn2est->nstep);
    double wtsumsum=0;
    const int nwfspair=cn2est->wt->nx;
    for(int iwfspair=0; iwfspair<nwfspair; iwfspair++){
	dmat *Pkn=cn2est->Pkn->p[iwfspair+iwfspair*nwfspair];
	dmat *wt=cn2est->wt->p[iwfspair];/*the layer weights. */
	dmat *ht=cn2est->ht->p[iwfspair];
	int nlayer=wt->nx;
	dmat *Pknwt=dnew(nlayer,nlayer);/*to remove columns in Pkn if zero weights are returned */
	/*number of negative layers found. */
 	int nfd=0;
	int nfdlast=0;
	double wtsum;
	dzero(Pknwt);
	daddI(Pknwt,1);/*start with identity matrix */
    repeat:
	nfd=0;
	wtsum=0; 

	for(int ix=0; ix<nlayer; ix++){
	    const double val=wt->p[ix];
	    if(val<0){
		/*negatives found. need to remove columns in the forward matrix
		  Pkn and redo the estimation matrix iPkn*/
		nfd++;
		/*disable diagonal elements. */
		Pknwt->p[ix*(1+nlayer)]=0;
	    }else{
		wtsum+=val;
	    }
	}
	if(nfd>nfdlast){
	    /*info("nfd=%d\n",nfd); */
	    if(fabs(wtsum)<1e-10){
		error("No valid layers found\n");
	    }
	    dmat *Pkn2=NULL;
	    /*this will zero out columns in the forward matrx */
	    dmm(&Pkn2, Pkn, Pknwt, "nn", 1);
	    /*will redo the estimation matrix. */
	    dmat *iPkn2=dpinv(Pkn2, NULL, NULL);
	    /*dwrite(wt,"wt");
	    dwrite(Pknwt,"Pknwt");
	    dwrite(Pkn,"Pkn");
	    dwrite(Pkn2,"Pkn2");
	    dwrite(iPkn2,"iPkn2");
	    exit(0);*/
	    dzero(wt);
	    /*compute the new result */
	    dmm(&wt, iPkn2, cn2est->cc->p[iwfspair], "nn", 1./cn2est->nstep);
	    dfree(iPkn2);
	    dfree(Pkn2);
	    nfdlast=nfd;
	    /*repeat the computation of wtsum */
	    goto repeat;
	}
	dfree(Pknwt);
	wtsumsum+=wtsum;
	double r0=pow(wtsum,-3./5.);
	cn2est->r0->p[iwfspair]=r0;
	dscale(wt, 1./wtsum);
	if(parms->cn2.verbose){
	    info2("r0=%.4fm theta0=%6f\" ",r0,calc_aniso(r0,wt->nx,ht->p,wt->p)*206265);
	    if(parms->ndm==2){
		info2("theta2=%6f\" ", calc_aniso2(r0,wt->nx,ht->p,wt->p,
						   parms->dm[0].ht,parms->dm[1].ht)*206265);
	    }
	    info2("wt=[");
	    for(int iht=0; iht<wt->nx; iht++){
		info2("%5f ", wt->p[iht]);
	    }
	    info2("]\n");
	}
    }
    cn2est->r0m=pow(wtsumsum/cn2est->wt->nx, -3./5.);
    dcellzero(cn2est->wtrecon);
    spcellmulmat(&cn2est->wtrecon, cn2est->wtconvert, cn2est->wt, 1);
    /*only 1 cell. norm to sum to 1. */
    normalize(cn2est->wtrecon->p[0]->p, cn2est->wtrecon->p[0]->nx, 1);
    if(parms->cn2.verbose){
	info2("r0m=%.4f theta0=%.4f\" ",cn2est->r0m, 
	      calc_aniso(cn2est->r0m,cn2est->wtrecon->p[0]->nx,
			 cn2est->htrecon->p,cn2est->wtrecon->p[0]->p)*206265);
	if(parms->ndm==2){
	    info2("theta2=%6f\" ", calc_aniso2(cn2est->r0m,cn2est->wtrecon->p[0]->nx,
					       cn2est->htrecon->p,cn2est->wtrecon->p[0]->p,
					       parms->dm[0].ht, parms->dm[1].ht)*206265);
	}
	info2("wt=[");
	for(int iht=0; iht<cn2est->wtrecon->p[0]->nx; iht++){
	    info2("%5f ", cn2est->wtrecon->p[0]->p[iht]);
	}
	info2("]\n");
    }
    /*divide by the number of accumulated frames. */
    if(cn2est->reset){
	info("reset the covariance");
	cn2est->nstep=0;/*reset the counter; */
	dcellzero(cn2est->cc);/*reset the numbers. */
    }
}
/**
   Implemented mechanism to move height of layers.
 */
void cn2est_moveht(RECON_T *recon){
    (void)recon;
    /*CN2EST_T *cn2est=recon->cn2est; */
    /*
      Implemented mechanism to move height of layers. Need to redo HXF, GX, etc.
    */
    error("moveht not yet implemented");
}
/**
   Update the tomographic reconstructor by updating L2.
 */
void cn2est_updatetomo(RECON_T *recon, const PARMS_T *parms){
    CN2EST_T *cn2est=recon->cn2est;
    /*wtrecon is referenced so should be updated automaticaly. */
    if(recon->wt->p!=cn2est->wtrecon->p[0]->p){
	warning("wtrecon is not referenced\n");
	dfree(recon->wt);
	recon->wt=dref(cn2est->wtrecon->p[0]);
    }
    recon->r0=cn2est->r0m;
    recon->l0=cn2est->l0;
    setup_recon_tomo_update(recon, parms);
}
/**
   Wrapper of Cn2 Estimation operations in recon.c
*/
void cn2est_isim(RECON_T *recon, const PARMS_T *parms, dcell *gradol, int isim){
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	if(recon->cn2est->wfscov[iwfs]){
	    cn2est_embed(recon->cn2est, gradol->p[iwfs], iwfs);
	}
    }
    CN2EST_T *cn2est=recon->cn2est;
    cn2est_cov(cn2est);/*convert gradients to cross covariance. */
    if((isim+1-parms->sim.start)%parms->cn2.step == 0){
	/*dcellswrite(cn2est->cc, 1./cn2est->nstep, "cc_%d",isim+1);*/
	cn2est_est(cn2est, parms);/*do the CN2 estimation */
	if(parms->cn2.moveht){
	    cn2est_moveht(recon);
	}
	if(parms->cn2.tomo){
	    if(parms->cn2.verbose){
		info2("Updating tomography weights\n");
	    }
	    cn2est_updatetomo(recon,parms);/*notice, cannot be parallel with tomofit(). */
	}
    }
}
/**
   Free all the data.
 */
void cn2est_free(CN2EST_T *cn2est){
    if(!cn2est) return;
    free(cn2est->embed);
    free(cn2est->wfscov);
    for(int iwfspair=0;iwfspair<cn2est->nwfspair; iwfspair++){
	for(int isep=0; isep<cn2est->pair[iwfspair].nsep; isep++){
	    free(cn2est->pair[iwfspair].sapair[isep]);
	}
	free(cn2est->pair[iwfspair].sapair);
	free(cn2est->pair[iwfspair].nsapair);
    }
    free(cn2est->pair);
    dcellfree(cn2est->gxs);
    dcellfree(cn2est->gys);
    dcellfree(cn2est->cxs);
    dcellfree(cn2est->cys);
    dcellfree(cn2est->cc);
    dcellfree(cn2est->iPkn);
    dcellfree(cn2est->wt);
    dcellfree(cn2est->ht);
    dfree(cn2est->r0);
    free(cn2est);
}

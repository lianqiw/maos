/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
/**
   Initialize cn2est_t
*/
cn2est_t* cn2est_new(const dmat* wfspair, /**<2n*1 vector for n pair of WFS indices.*/
	const dmat* wfstheta,/**<nwfs*2: angular direction of each WFS.*/
	const loc_t* saloc,  /**<nsa*2: Subaperture low left corner coordinates*/
	const dmat* saa,     /**<nsa*1: Normalized subaperture area*/
	const real saat,     /**<Threshold for keeping subapertures*/
	const dmat* hs,      /**<nwfs*1: altitude of guide star*/
	const dmat* htrecon, /**<Layers height intended for tomography*/
	int keepht,          /**<2: slodar directly to htrecon,
			  * otherwise: interpolate onto htrecon
			  * from native slodar heights*/
	real L0            /**<The Outer scale*/
){
	info("Cn2 estimation: ");
	/*We need at least a pair */
	if(!wfspair) return 0;
	int npair=NX(wfspair)*NY(wfspair);
	if(!npair) return 0;
	if(npair%2==1){
		error("pair must have even number of entries\n");
	}
	/*>>1 is a short cut for /2 */
	int nwfspair=npair>>1;
	cn2est_t* cn2est=mycalloc(1, cn2est_t);
	cn2est->nwfspair=nwfspair;
	int nwfs=NX(wfstheta);
	cn2est->nwfs=nwfs;
	cn2est->nsa=saloc->nloc;
	/*wfscov is a flag for wfs showing whether this wfs participates in
	 * covariance computation. */
	cn2est->wfscov=mycalloc(nwfs, int);
	for(int ind=0; ind<npair; ind++){
		int iwfs=(int)P(wfspair,ind);
		if(iwfs<0||iwfs>=nwfs){
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
	lmat* mask=lnew(nx, nx);
	if(saa&&NX(saa)!=cn2est->nsa){
		error("saa and saloc mismatch\n");
	}

	real saat2=dmax(saa)*saat;
	for(int isa=0; isa<cn2est->nsa; isa++){
		if(!saa||P(saa,isa)>saat2){
			P(mask,P(cn2est->embed,isa))=1;/*use this subaperture */
		}
	}
	/* this mask for for subapertures that we can compute curvature.*/
	cn2est->mask=lnew(nx, nx);
	cmat* overlap=cnew(nx, nx);
	int iymin=nx, iymax=0, ixmin=nx, ixmax=0;
	for(int iy=0; iy<nx; iy++){
		for(int ix=0; ix<nx; ix++){
			/*Only use a subaperture if we are able to make curvature. */
			if(P(mask, ix, iy)&&P(mask, ix+1, iy)&&P(mask, ix-1, iy)&&
				P(mask, ix, iy+1)&&P(mask, ix, iy-1)){
				P(cn2est->mask, ix, iy)=1;
				P(overlap, ix, iy)=1;
				if(ix>ixmax) ixmax=ix;
				if(ix<ixmin) ixmin=ix;
				if(iy>iymax) iymax=iy;
				if(iy<iymin) iymin=iy;
			}
		}
	}
	lfree(mask);
	int maxsep=MIN((iymax-iymin), (ixmax-ixmin));
	if(maxsep<0){
		warning("slodar is not possible with not enough subapertures\n");
		cn2est_free(cn2est);
		free(overlap);
		return 0;
	}

	cfft2(overlap, -1);
	for(long i=0; i<nxnx; i++){
		P(overlap,i)=P(overlap,i)*conj(P(overlap,i));
	}
	cfft2(overlap, 1);
	cfftshift(overlap);
	/*cn2est->overlap: number of overlaps for each separation. peak value at (nx/2, nx/2)*/
	cn2est->overlapi=dnew(nx, nx);
	for(long i=0; i<nxnx; i++){
		real over=creal(P(overlap,i));
		P(cn2est->overlapi,i)=(over>MIN_SA_OVERLAP)?(nxnx/over):0;
	}
	cfree(overlap);
	/*2-d arrays to store x y gradient, and "curvature" */
	cn2est->gxs=dcellnew(nwfs, 1);/*stores gradient in 2-d map */
	cn2est->gys=dcellnew(nwfs, 1);
	cn2est->curi=ccellnew(nwfs, 1);
	for(int ix=0; ix<nwfs; ix++){
		if(cn2est->wfscov[ix]){
			P(cn2est->gxs,ix)=dnew(nx, nx);
			P(cn2est->gys,ix)=dnew(nx, nx);
			P(cn2est->curi,ix)=cnew(nx, nx);
		}
	}
	/*
	  Now we prepare indexes that can be used to compute gradient(curvature)
	  cross-covariances conveniently during simulation.
	*/
	/*first get a few constants */
	const real dsa=saloc->dx;
	/*determine the layer height used for tomography. */
	cn2est->htrecon=ddup(htrecon);
	cn2est->wtrecon=dcellnew(1, 1);
	P(cn2est->wtrecon,0)=dnew(NX(cn2est->htrecon), 1);
	{
		info("htrecon=[");
		for(int iht=0; iht<NX(cn2est->htrecon); iht++){
			info("%.2f ", P(cn2est->htrecon,iht)*0.001);
		}
		info("]km\n");
	}
	/*stores cross-covariance data during simulation */
	cn2est->covc=ccellnew(nwfspair, 1);
	cn2est->cov1=dcellnew(nwfspair, 1);
	cn2est->cov2=dcellnew(nwfspair, 1);
	/*height of layers for each wfs pair of slodar output. */
	cn2est->ht=dcellnew(nwfspair, 1);
	/*record sapair to use for each separation */
	cn2est->pair=mycalloc(nwfspair, cn2est_pair_t);
	long nhtsx[nwfspair];
	long nhtsy[nwfspair];
	real hmin, hmax;
	dmaxmin(P(cn2est->htrecon), NX(cn2est->htrecon), &hmax, &hmin);
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
	cn2est->Pnk=dcellnew(nwfspair, nwfspair);
	cn2est->iPnk=dcellnew(nwfspair, nwfspair);
	/*wtconvert is the matrix to down/up sample the CN2 estimates to layers
	  used for tomography*/
	cn2est->wtconvert=dspcellnew(1, nwfspair);
	cn2est->L0=L0;
	for(int iwfspair=0; iwfspair<nwfspair; iwfspair++){
	/*get pointer for this pair */
		cn2est_pair_t* pair=cn2est->pair+iwfspair;
		/*The WFS in this pair. */
		const int wfs0=P(wfspair,iwfspair*2);
		const int wfs1=P(wfspair,iwfspair*2+1);
		pair->wfs0=wfs0;
		pair->wfs1=wfs1;
		/*The separation between the stars */
		const real dthetax=P(wfstheta, wfs0, 0)-P(wfstheta, wfs1, 0);
		const real dthetay=P(wfstheta, wfs0, 1)-P(wfstheta, wfs1, 1);
		/*the angular distance between WFS */
		real dtheta=sqrt(dthetax*dthetax+dthetay*dthetay);
		/*The direction of the WFS pair baseline vector */
		real beta=atan2(dthetay, dthetax);
		if(wfs0==wfs1||dtheta<1e-14){
			dtheta=0;
			beta=0;
		}
		/*The average wfs guide star height of the two*/
		const real cb=cos(beta);
		const real sb=sin(beta);
#if COV_ROTATE
		const real slang=1;
#else
		const real slang=MAX(fabs(cb), fabs(sb));
#endif
		real hsm=0;
		if(NX(hs)*NY(hs)>=nwfs){
			hsm=(P(hs,wfs0)+P(hs,wfs1))*0.5;
		} else if(NX(hs)*NY(hs)==1){
			hsm=P(hs,0);
		} else{
			error("hs in in wrong format, should have 1 or %d elements\n", nwfs);
		}
		if(dtheta==0){//seeing only for self-correlation.
			pair->iht0=0;
			pair->iht1=1;
			warning("Can only obtain r0 using auto-correlation.\n");
		} else{
			/*The number of layers is determined */
			if(keepht==2){
			/*slodar output directly to layers used for tomography. may not work well */
				pair->iht0=0;
				pair->iht1=NX(cn2est->htrecon);
			} else{
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
		P(cn2est->ht,iwfspair)=dnew(pair->nht, 1);
		/*stores the cross-covariance of the curvature for x, and y grads. */
		P(cn2est->covc,iwfspair)=cnew(nx, nx);
		P(cn2est->cov2,iwfspair)=dnew(nx, nx);
#if COV_ROTATE
		P(cn2est->cov1,iwfspair)=dnew(pair->nsep, 1);//1d cut
#else
		P(cn2est->cov1,iwfspair)=dref_reshape(P(cn2est->cov2,iwfspair), nxnx, 1);
#endif
	/*info2("Pair %d: wfs %d and %d. dtheta=%4f\" iht=[%d, %d)\n",
	  iwfspair, wfs0, wfs1, pair->dtheta*RAD2AS, pair->iht0, pair->iht1);*/


	/*
	  Now we start to build the model that will be used to estimate Cn2 from the
	  cross-covariance of curvature.
	*/
		const int nsep=pair->nsep;
		/*nm is the size of the mc array */
		const int nm=4*(MAX(abs(pair->iht0), abs(pair->iht1-1))+1)*cn2est->ovs;
		const int nm2=nm>>1;
		/*sampling in mc after FFT */
		const real dx=dsa/cn2est->ovs;
		const real df=1./(nm*dx);
		const real L02=pow(cn2est->L0, -2);
		/*initialize */
		cmat* mc=cnew(nm, nm);
		/*create 2-d pointers */
		cmat* pmc=mc;
		/*the forward operator from layer weights to cross-covariance */
		dmat* Pnk=dnew(nsep, pair->nht);
		info("Pair %d: hk=[", iwfspair);
		for(int iht=pair->iht0; iht<pair->iht1; iht++){
			/*the height of each layer */
			real hk;
			if(keepht==2){//do not use.
				hk=P(cn2est->htrecon,iht-pair->iht0);
			} else if(dtheta==0){
				hk=0;
			} else{
				hk=dsa*iht/(dtheta*slang+dsa*iht/hsm);
			}
			info("%.2f ", hk*0.001);
			P(P(cn2est->ht,iwfspair),iht-pair->iht0)=hk;
			/*the cone effect */
			const real zeta=1.-hk/hsm;
			const real zetan2=pow(zeta, -2);
			/*
			  coefficients for the PSD. zeta is cone effect. df is the frequency
			  bin. df*df is from the descretization of the PSD. we use df
			  instead of 1/dx like in the matlab code because our inverse fft
			  does not apply (1/N)^ scaling.
			*/
			/*the coefficients in fron of the PSD. r0 is left out. */
			const real psd_coef=pow(zeta, -2)*4*M_PI*M_PI
				*0.0229*pow(2*M_PI/0.5e-6, -2)*(df*df);
			for(int iy=0; iy<nm; iy++){
			/*0 freq at corner so no need to do fftshift before fft */
				const real fy=(iy-nm2)*df;
				/*the sinc subaperture fnction */
				const real sincfy=sinc(fy*dsa);
				for(int ix=0; ix<nm; ix++){
					const real fx=(ix<nm2?ix:ix-nm)*df;
					const real sincfx=sinc(fx*dsa);
					/*the turbulence PSD with outerscale */
					const real psd=psd_coef*pow((fx*fx+fy*fy)*zetan2+L02, -11./6.);
					/*gx diff is along x, gy diff is along y to form real curvature */
					const real cur=pow(2*fx*(cos(2*M_PI*dsa*fx)-1)+2*fy*(cos(2*M_PI*dsa*fy)-1), 2);
					P(pmc, ix, iy)=pow(sincfy*sincfx, 2)*psd*cur;
				}/*ix */
			}/*iy */
			/*doing fft */
			cfft2(mc, 1);/*inverse fft */
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
				real xx=(iht*cb/slang-xsep)*cn2est->ovs+nm2;
				real yy=(iht*sb/slang-ysep)*cn2est->ovs+nm2;

				if(xx<0||yy<0||xx+1>=nm||yy+1>=nm){
					//warning("Out of range: xx=%g, yy=%g\n", xx, yy);
				} else{
#if INTERP_NEAREST
			/*Do interpolation using nearest neighbor */
					int ixx=(int)round(xx);
					int iyy=(int)round(yy);
					real imc=creal(P(pmc, ixx, iyy));
#else
			/*Do interpolation using bilinear spline interp. */
					int ixx=(int)floor(xx); xx=xx-ixx;
					int iyy=(int)floor(yy); yy=yy-iyy;
					real imc=creal((P(pmc, ixx, iyy)*(1-xx)+P(pmc, ixx+1, iyy)*(xx))*(1-yy)
						+(P(pmc, ixx, iyy+1)*(1-xx)+P(pmc, ixx+1, iyy+1)*(xx))*yy);
#endif
					P(Pnk, isep, iht-pair->iht0)=imc;
				}
			}
		}
		info("] km\n");
		/*
		  iPnk is a block diagonal matrix for Cn2 Estimation.
		*/
		P(cn2est->Pnk, iwfspair, iwfspair)=dref(Pnk);
		P(cn2est->iPnk, iwfspair, iwfspair)=dpinv(Pnk, NULL);
		P(cn2est->wtconvert,iwfspair)=mkhbin1d(P(cn2est->ht,iwfspair), cn2est->htrecon);
		dfree(Pnk);
		cfree(mc);
	}/*iwfspair */
/*stores estimated weight of layers during simulation and output to file finally. */
	cn2est->wt=dcellnew3(nwfspair, 1, nhtsx, nhtsy);
/*stores estimated r0 during simulation */
	cn2est->r0=dnew(nwfspair, 1);
	return cn2est;
}/*cn2est_new */

/**
   Compute covariance from gradients and embed into a 2d map cn2est->curi.
*/
static void cn2est_embed(cn2est_t* cn2est, const dcell* gradol, int icol){
	long* embed=P(cn2est->embed);
	for(int iwfs=0; iwfs<NX(gradol); iwfs++){
		if(!cn2est->wfscov[iwfs]) continue;
		const int nsa=cn2est->nsa;
		dmat* grad=P(gradol,iwfs);
		if(!grad){
			error("wfs %d: PSOL grads is required to do cn2 estimation\n", iwfs);
		}
		if(icol>=NY(grad)){
			error("icol=%d is invalid\n", icol);
		}
		if(NX(grad)!=nsa*2){
			error("grad and saloc does not match\n");
		}
		dmat* gx=P(cn2est->gxs,iwfs);
		dmat* gy=P(cn2est->gys,iwfs);
		const real* pgrad=PCOL(grad,icol);
		/*Embed gradients in a 2-d array */
		for(int isa=0; isa<nsa; isa++){
			P(gx,embed[isa])=pgrad[isa];
			P(gy,embed[isa])=pgrad[isa+nsa];
		}
		/*Compute curvature of wavefront from gradients. */
		cmat* cur=P(cn2est->curi,iwfs);
		const int ny=P(cn2est->curi,iwfs)->ny;
		const int nx=P(cn2est->curi,iwfs)->nx;
		for(int iy=0; iy<ny; iy++){
			for(int ix=0; ix<nx; ix++){
				if(P(cn2est->mask, ix, iy)){
					P(cur, ix, iy)=P(gx, ix+1, iy)+P(gx, ix-1, iy)-2*P(gx, ix, iy)/*gx along x */
						+P(gy, ix, iy+1)+P(gy, ix, iy-1)-2*P(gy, ix, iy);/*gy along y */
				} else{
					P(cur, ix, iy)=0;//must set to zero.
				}
			}
		}
		cfft2(P(cn2est->curi,iwfs), -1);
	}
}
/**
   Compute cross-covairance from gradient curvature
*/
static void cn2est_cov(cn2est_t* cn2est){
	/*accumulate cross-covariance of gradient curvature. */
	const int nwfspair=cn2est->nwfspair;
	cn2est->count++;
	/*for each wfs pair */
	for(int iwfspair=0; iwfspair<nwfspair; iwfspair++){
		const int wfs0=cn2est->pair[iwfspair].wfs0;
		const int wfs1=cn2est->pair[iwfspair].wfs1;
		const comp* cur1=P(P(cn2est->curi,wfs0));
		const comp* cur2=P(P(cn2est->curi,wfs1));
		comp* cov=P(P(cn2est->covc,iwfspair));
		for(long i=0; i<cn2est->nembed*cn2est->nembed; i++){
			cov[i]+=conj(cur1[i])*(cur2[i]);
		}
	}/*iwfspair */
}
/**
   Accumulate the coveriance.
 */
void cn2est_push(cn2est_t* cn2est, const dcell* gradol){
	int ncol=0;
	if(NX(gradol)<cn2est->nwfs){
		error("Grad has less number of wfs than required %d\n", cn2est->nwfs);
	}
	for(int iwfs=0; iwfs<cn2est->nwfs; iwfs++){
		if(!cn2est->wfscov[iwfs]) continue;
		if(ncol==0){
			ncol=P(gradol,iwfs)->ny;
		} else if(ncol!=P(gradol,iwfs)->ny){
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
void cn2est_est(cn2est_t* cn2est, int verbose){
	info("cn2est from %d measurements\n", cn2est->count);
	cmat* covi=cnew(cn2est->nembed, cn2est->nembed);
#if COV_ROTATE
	dmat* covr=dnew(cn2est->nembed, cn2est->nembed);
#endif
	const int nwfspair=NX(cn2est->wt);
	for(int iwfspair=0; iwfspair<nwfspair; iwfspair++){
		cadd(&covi, 0, P(cn2est->covc,iwfspair), 1./(cn2est->count*cn2est->nembed*cn2est->nembed));
		cfft2(covi, 1);
		cfftshift(covi);
		creal2d(&P(cn2est->cov2,iwfspair), 0, covi, 1);
#if COV_ROTATE
	//roate and embed;
		dembed(covr, P(cn2est->cov2,iwfspair), -cn2est->pair[iwfspair].beta);
		real* cc=P(P(cn2est->cov1,iwfspair));
		cn2est_pair_t* pair=cn2est->pair+iwfspair;
		int off=cn2est->nembed/2;
		for(long isep=pair->iht0; isep<pair->iht1; isep++){
			cc[isep-pair->iht0]=P(covr, off+isep, off)*P(cn2est->overlapi, off+isep, off);
		}
#else
		for(long isep=0; isep<NX(covi)*NY(covi); isep++){
			P(P(cn2est->cov2,iwfspair),isep)*=P(cn2est->overlapi,isep);//temporary. apply to mc instead.
		}
#endif
	}
	cfree(covi);
#if COV_ROTATE
	dfree(covr);
#endif
	dcellzero(cn2est->wt);
	dcellmm(&cn2est->wt, cn2est->iPnk, cn2est->cov1, "nn", 1);

	real wtsumsum=0;
	for(int iwfspair=0; iwfspair<nwfspair; iwfspair++){
		real wtsum=0;
		dmat* wt=P(cn2est->wt,iwfspair);//the layer weights. 
		dmat* ht=P(cn2est->ht,iwfspair);
		int nlayer=NX(wt);
		if(CN2EST_NO_NEGATIVE){
			/*The following tries to remove negative weights.  For small
			  negatives, we just zero them. For significant ones, we remove the
			  layer from being estimated to avoid biasing the result.
			*/
			dmat* Pnk=ddup(P(cn2est->Pnk, iwfspair, iwfspair));
			//number of negative layers found. 
			int nfd=0, nfdi=0;
			do{
				nfdi=0;
				for(int ix=0; ix<nlayer; ix++){
					const real val=P(wt,ix);
					if(val<0){
						if(val<-1e-2){
							//negatives found. Remove columns in the forward matrix to disable this point.
							//Pnk and redo the estimation matrix iPnk
							nfd++;
							nfdi++;
							dmat* tmp=drefcols(Pnk, ix, 1);
							dzero(tmp);
							dfree(tmp);
						} else{
							P(wt,ix)=0;
						}
					}
				}
				if(nfdi>0){
					info_once("Ignore %d negative weights. set MAOS_CN2EST_NO_NEGATIVE=0 in your shell to disable this feature.\n", nfd);
					dmat* iPnk2=dpinv(Pnk, NULL);
					//compute the new result 
					dmm(&wt, 0, iPnk2, P(cn2est->cov1,iwfspair), "nn", 1);
					dfree(iPnk2);
				}
			} while(nfdi);
		}
		for(int ix=0; ix<nlayer; ix++){
			wtsum+=P(wt,ix);
		}
		wtsumsum+=wtsum;
		real r0=1.;
		if(wtsum>0){
			r0=pow(wtsum, -3./5.);
			P(cn2est->r0,iwfspair)=r0;
		}
		dscale(wt, 1./wtsum);
		if(verbose){
			info2("Pair%d: r0=%.4fm theta0=%.2f\" ", iwfspair, r0, calc_aniso(r0, NX(wt), P(ht), P(wt))*RAD2AS);
			if(cn2est->dmht&&NX(cn2est->dmht)==2){
				info2("theta2=%.2f\" ", calc_aniso2(r0, NX(wt), P(ht), P(wt),
					P(cn2est->dmht,0), P(cn2est->dmht,1))*RAD2AS);
			}
			info2("wt=[");
			for(int iht=0; iht<NX(wt); iht++){
				info2("%.4f ", P(wt,iht));
			}
			info2("]\n");
		}
	}
	cn2est->r0m=pow(wtsumsum/NX(cn2est->wt), -3./5.);
	dcellzero(cn2est->wtrecon);
	dcellmm(&cn2est->wtrecon, cn2est->wtconvert, cn2est->wt, "nn", 1);
	/*only 1 cell. norm to sum to 1. */
	dnormalize_sumabs(P(P(cn2est->wtrecon,0)), P(cn2est->wtrecon,0)->nx, 1);
	if(verbose){
		info2("Mean : r0=%.4fm theta0=%.2f\" ", cn2est->r0m,
			calc_aniso(cn2est->r0m, P(cn2est->wtrecon,0)->nx,
				P(cn2est->htrecon), P(P(cn2est->wtrecon,0)))*RAD2AS);
		if(cn2est->dmht&&NX(cn2est->dmht)==2){
			info2("theta2=%.2f\" ", calc_aniso2(cn2est->r0m, P(cn2est->wtrecon,0)->nx,
				P(cn2est->htrecon), P(P(cn2est->wtrecon,0)),
				P(cn2est->dmht,0), P(cn2est->dmht,1))*RAD2AS);
		}
		info2("wt=[");
		for(int iht=0; iht<P(cn2est->wtrecon,0)->nx; iht++){
			info2("%.4f ", P(P(cn2est->wtrecon,0),iht));
		}
		info2("]\n");
	}
}
/**
   Reset the accumulation
 */
void cn2est_reset(cn2est_t* cn2est){
	if(cn2est->count){
		info("cn2est: reset the accumulation\n");
		cn2est->count=0;/*reset the counter; */
		ccellzero(cn2est->covc);/*reset the numbers. */
	}
}

/**
   Free all the data.
*/
void cn2est_free(cn2est_t* cn2est){
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
/**
   One stop wrapper for #cn2est_new, #cn2est_push, and #cn2est_est.
 */
cn2est_t* cn2est_all(const dmat* wfspair, dmat* wfstheta, const loc_t* saloc,
	const dmat* saa, const real saat,
	const dmat* hs, const dmat* htrecon, int keepht, real l0, dcell* grad){
	if(dmaxabs(wfstheta)>1){
		dmat* tmp=wfstheta;
		wfstheta=ddup(tmp);
		dfree(tmp);
		//Don't scale the matlab one.
		dscale(wfstheta, 1.*AS2RAD);
	}
	if(NX(grad)==1&&NY(grad)>1){
		reshape(grad, NY(grad), 1);
	}
	struct cn2est_t* cn2est=cn2est_new(wfspair, wfstheta, saloc, saa, saat, hs, htrecon, keepht, l0);
	cn2est_push(cn2est, grad);
	cn2est_est(cn2est, 1);
	return cn2est;
}

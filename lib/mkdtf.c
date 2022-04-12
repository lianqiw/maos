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
#include "mkdtf.h"
#include "mkh.h"
/**
   mkdtf() computes the parameters used to sample PSFs onto detectors in the
   Fourier domain. It incorporates size of the detector pixels, and charge
   diffusion. For polar coordinate detectors, it also rotates the pixels and
   pixel coordinates so that they are aranged along the radial and azimuthal
   direction.
 */
dtf_t* mkdtf(const dmat* wvls, /**<List of wavelength*/
	real dxsa,        /**<Subaperture size*/
	real embfac,      /**<Embedding factor (2)*/
	long notfx,       /**<FFT size along x*/
	long notfy,       /**<FFT size along y*/
	long pixpsax,     /**<Number of pixels along x(r)*/
	long pixpsay,     /**<Number of pixels along y(a)*/
	real pixthetax,   /**<Pixel size along x (r)*/
	real pixthetay,   /**<Pixel size along y (a)*/
	const dmat* pixoffx,  /**<offset of image center from center of pixel array, along x or radial*/
	const dmat* pixoffy,  /**<offset of image center from center of pixel array, along y or azimuthal*/
	real pixblur,     /**<Pixel blur sigma(fraction of pixel)*/
	const dcell* pixrot /**<Rotation angle of pixels islands in each subaperture. for polar coordinate only*/
){
	int nwvl=NX(wvls)*NY(wvls);
	dtf_t* dtfs=mycalloc(nwvl, dtf_t);
	dtfs->nwvl=nwvl;
	const real blurx=pixblur*pixthetax;
	const real blury=pixblur*pixthetay;
	const real e0x=-2*M_PI*M_PI*blurx*blurx;//blurring factors
	const real e0y=-2*M_PI*M_PI*blury*blury;
	const int do_blur=fabs(blurx)>EPS&&fabs(blury)>EPS;
	const long notfx2=notfx>>1;
	const long notfy2=notfy>>1;
	const real pxo=-(pixpsax*0.5-0.5)*pixthetax;
	const real pyo=-(pixpsay*0.5-0.5)*pixthetay;
	
	const int nwfs=MAX(pixoffx?NY(pixoffx):1, pixrot?NX(pixrot):1);
	const int nsa=MAX(pixoffx?NX(pixoffx):1, pixrot?NX(P(pixrot, 0)):1);

	/*both nominal and si depends on wavelength.*/
	for(int iwvl=0; iwvl<nwvl; iwvl++){
		const real wvl=P(wvls,iwvl);
		const real dtheta=wvl/(dxsa*embfac);/*PSF sampling. */
		const real dux=1./(dtheta*notfx);
		const real duy=1./(dtheta*notfy);
		const real dux2=dux*dux;
		const real duy2=duy*duy;
		const real pdtheta=pixthetax*pixthetay/(dtheta*dtheta);//scaling factor due to binning into detectors
		const real duxp=dux*pixthetax;
		const real duyp=duy*pixthetay;
		/*For LGS WFS (srot!=NULL), there is elongation effect. For radial pixel
		  CCD without using rotating psf/otf method, we need to create
		  nominal/si for each subaperture. The PSF/OTF and hence DTF are on x/y
		  coordinate, so pixels and pixel coordinates need to rotated to r/a
		  direction. */
		//When DTF is along x/y coord, pixel is r/a coord, or when there is pixoffx, we need multiple DTFs.
		
		dtfs[iwvl].dtheta=dtheta;
		dtfs[iwvl].dxsa=dxsa;
		dtfs[iwvl].radpix=pixrot?1:0;
		dtfs[iwvl].wvl=wvl;
		dtfs[iwvl].notfx=notfx;
		dtfs[iwvl].notfy=notfy;
		dtfs[iwvl].pixpsax=pixpsax;
		dtfs[iwvl].pixpsay=pixpsay;
		dtfs[iwvl].pixthetax=pixthetax;
		dtfs[iwvl].pixthetay=pixthetay;
		dtfs[iwvl].nominal=ccellnew_same(nsa, nwfs, notfx, notfy);
		dtfs[iwvl].si=dspcellnew(nsa, nwfs);
		ccell* nominals=dtfs[iwvl].nominal;
		dspcell* sis=dtfs[iwvl].si;
		cmat* nominal=cnew(notfx, notfy);
		//Coordinate of PSF points
		loc_t* loc_psf=mksqloc(notfx, notfy, dtheta, dtheta, -notfx2*dtheta, -notfy2*dtheta);
		real theta=0;
		real ct=cos(theta);
		real st=sin(theta);
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			for(int isa=0; isa<nsa; isa++){
				if(pixrot){
					theta=P(PR(pixrot,iwfs,0),isa);
					ct=cos(theta);
					st=sin(theta);
				}

				for(int iy=0; iy<notfy; iy++){
					int jy=iy-notfy2;
					for(int ix=0; ix<notfx; ix++){
						int jx=ix-notfx2;
						real ir=ct*jx+st*jy;
						real ia=-st*jx+ct*jy;
						//Pixel function
						P(nominal, ix, iy)=sinc(ir*duxp)*sinc(ia*duyp)*pdtheta;
						if(do_blur){//Charge diffusion.
							P(nominal, ix, iy)*=exp(e0x*(ir*ir*dux2)+e0y*(ia*ia*duy2));
						}
					}
				}
				/*put peak in corner. */
				cfftshift(nominal);
				/*Embed a fftshift operation in nominal so that we avoid fftshift when forming i0.*/
				cfft2(nominal, -1);
				cfftshift(nominal);
				cfft2(nominal, 1);
				//cancel FFT scaling effect.
				cscale(nominal, 1./(real)(NX(nominal)*NY(nominal)));
				ccp(&P(nominals, isa, iwfs), nominal);
				//Coordinate of pixels
				
				real dx=pixoffx?(PR(pixoffx, isa, iwfs)*pixthetax):0;
				real dy=pixoffy?(PR(pixoffy, isa, iwfs)*pixthetay):0;
			
				loc_t* loc_ccd=mksqloc(pixpsax, pixpsay, pixthetax, pixthetay, pxo-dx, pyo-dy);
				locrot(loc_ccd, theta);
				P(sis, isa, iwfs)=mkh(loc_psf, loc_ccd, 0, 0, 1);
				locfree(loc_ccd);
			}/*isa */
		}/*iwfs */
		cfree(nominal);
		locfree(loc_psf);

		/*Create an excessive high frequency in nominal so that
		  we don't have to do fftshift later.*/
		dtfs[iwvl].Ux=cnew(notfx, 1);
		dtfs[iwvl].Uy=cnew(notfy, 1);
		comp* Ux=P(dtfs[iwvl].Ux);
		comp* Uy=P(dtfs[iwvl].Uy);

		/*The following is used in genseotf to compute shifted i0.*/
		for(int ix=0; ix<notfx; ix++){
			int jx=ix<notfx2?ix:(ix-notfx);
			Ux[ix]=COMPLEX(0, -2.*M_PI*jx*dux);
		}
		for(int iy=0; iy<notfy; iy++){
			int jy=iy<notfy2?iy:(iy-notfy);
			Uy[iy]=COMPLEX(0, -2.*M_PI*jy*duy);
		}

	}/*iwvl */
	return dtfs;
}

/**
   Wrap the index for dataset with total of n frames for continuity. The actual data becomes
   0, 1, 2, ..., n-2, n-1, n-2, ..., 0, 1
*/
static inline int wrap_seq(long index, long n){
	long m=n*2-1;
	index=index%m;
	if(index<0) index+=m;
	if(index>=n) index=m-1-index;
	return index;
}

etf_t* mketf(const dtf_t* dtfs,  /**<The dtfs*/
	const cell* sodium,/**<The sodium profile. In each cell First column is coordinate. Different cells must have the same coordinate. Can also be a dmat*/
	int icol,     /**<Which sodium profile to use*/
	const dcell* srot,  /**<Rotation angle of each subaperture. NULL for NGS WFS*/
	const dcell* srsa,  /**<Subaperture to LLT distance*/
	real hs,      /**<LGS WFS focus altitude*/
	real htel,    /**<Telescope altitude*/
	real za_rad,   /**<Zenith angle in radian*/
	int no_interp /**<Use direct sum instead of interpolation + FFT. Slower */
){
	int nwvl=dtfs[0].nwvl;
	etf_t* etfs=mycalloc(nwvl, etf_t);
	etfs->nwvl=nwvl;
	etfs->hs=hs;
	/*setup elongation along radial direction. don't care azimuthal. */
	if(!srot) error("srot is required");
	const int nllt=MAX(NX(srot), iscell(sodium)?NX(sodium):1);
	const int nsa=P(srot,0)->nx;
	dmat* sodium0=dmat_cast(iscell(sodium)?P(sodium,0):sodium);
	const int ncol=sodium0->ny-1;
	const int nhp=NX(sodium0);
	const real* px=P(sodium0);
	const real cosza=cos(za_rad);
	const int icolwrap=wrap_seq(icol, ncol);
	//info("Na using column %d.\n",icol);

	//adjusting sodium height for the zenith angle;
	real hpmin=0, dhp1=0;
	if(nhp==1){
		no_interp=1;
	}
	if(!no_interp){/**Requires linear spacing*/
		hpmin=(px[0]-htel)/cosza;
		dhp1=cosza/(px[1]-px[0]);
		/*assume linear spacing. check the assumption valid */
		real diff;
		if((diff=fabs((px[nhp-1]-px[0])/((nhp-1)*(px[1]-px[0]))-1.))>1.e-5){
			error("sodium profile is not evenly spaced: %g\n", diff);
		}
	}//else: Any spacing is ok.

	/*the sum of pp determines the scaling of the pixel intensity. */
	real* pna[nllt];
	real i0scale[nllt];
	for(int illt=0; illt<nllt; illt++){
		dmat* sodiumi=dmat_cast(iscell(sodium)?PR(sodium, illt, 0):sodium);
		pna[illt]=PCOL(sodiumi, 1+icolwrap);
		i0scale[illt]=dvecsum(pna[illt], nhp);
		if(fabs(i0scale[illt]-1)>0.01){
			warning("Siglev is scaled by %g by sodium profile\n", i0scale[illt]);
		}
	}
	for(int iwvl=0; iwvl<nwvl; iwvl++){
		const real dtheta=dtfs[iwvl].dtheta;
		const long notfx=dtfs[iwvl].notfx;
		const long notfy=dtfs[iwvl].notfy;
		const long notfx2=notfx>>1;
		const long notfy2=notfy>>1;
		const real dux=1./(dtheta*notfx);
		const real duy=1./(dtheta*notfy);
		
		/*
		  PSF/OTF/ETF is defined in x/y coordinate.
		  2010-01-04: Fuse dtf nominal into etf for this case.
		*/
		ccell* petf=etfs[iwvl].etf=ccellnew_same(nsa, nllt, notfx, notfy);

		etfs[iwvl].icol=icol;
		if(no_interp){ /* No interpolation, no fft. intensity scaling is automatically taken care of */
			OMP_TASK_FOR_COLLAPSE(2)		
			for(int illt=0; illt<nllt; illt++){
				for(int isa=0; isa<nsa; isa++){
					const real rsa=P(P(srsa,illt),isa);
					const real theta=P(PR(srot, illt, 0), isa);
					const real ct=cos(theta);
					const real st=sin(theta);
					cmat* etf2d=P(petf, isa, illt);
					for(int icompy=0; icompy<notfy; icompy++){
						const real ky=duy*(icompy>=notfy2?(icompy-notfy):icompy);
						for(int icompx=0; icompx<notfx; icompx++){
							const real kx=dux*(icompx>=notfx2?(icompx-notfx):icompx);
							const real kr=(ct*kx+st*ky);/*along radial*/
							for(int ih=0; ih<nhp; ih++){
								if(pna[illt][ih]>0){
									const real tmp=(-2*M_PI*(kr*(rsa*cosza/(px[ih]-htel)-rsa/hs)));
									P(etf2d, icompx, icompy)+=COMPLEX(pna[illt][ih]*cos(tmp), pna[illt][ih]*sin(tmp));
								}
							}
						}
					}
				}//isa
			}//illt
		} else{
			/*
			  The ETF is computed as DFT
			  ETF(k_i)=\sum_j exp(-2*\pi*I*k_i*{\theta}_j)P({\theta}_j)
			  Where \theta_j is the radial coord of pixel j, and
			  P({\theta}_j}=sodium(h_j) with h_j=rsa*hs/(rsa+hs*{\theta}_j)
			  where hs is GS distance and rsa is subaperture to LLT distance
			  We replace the above interpolation (and FFT) by doing summation directly
			  ETF(k_i)=\sum_j exp(-2*\pi*I*k_i*(rsa/h_j-rsa/hs))P(h_j)
			*/
			const int npad=2;/*zero padding to reduce aliasing */
			const int nover=2;/*enough size for rotation*/
			const int netf=notfx*nover*npad;
			const real dtetf=dtheta/nover;
			const real dusc=(netf*dtetf)/(dtheta*notfx);
			dmat* thetas=dnew(netf, 1);
			const int netf2=netf>>1;
			/*Only interpolating the center part. the rest is padding. */
			const int etf0=netf2-(int)round(notfx2*(dtheta/dtetf));
			const int etf1=etf0+(int)round(notfx*(dtheta/dtetf));
			if(etf0<0) error("Invalid configuration\n");
			for(int it=etf0; it<etf1; it++){
				P(thetas,it)=(it-netf2)*dtetf;
			}
			int ncrop=0;
			real max_crop=0;
			ccell *etf_cache=ccellnew_same(MAXTHREAD,1,netf, 1);
			OMP_TASK_FOR_COLLAPSE(2)
			for(int illt=0; illt<nllt; illt++){
				for(long isa=0; isa<nsa; isa++){
					/*1d ETF along radius. */
					int ithread=0;
#ifdef _OPENMP
					ithread=omp_get_thread_num();
#endif									
					cmat* etf=P(etf_cache, ithread);
					real rsa=P(P(srsa,illt),isa);
					real etf2sum=0;
					czero(etf);
					for(int icomp=etf0; icomp<etf1; icomp++){
					/*peak in center */
						const real itheta=P(thetas,icomp);
						/*non linear mapping. changed from - to + on 2014-05-07.*/
						const real ih=hs*rsa/(rsa+itheta*hs);
						/*interpolating to get Na profile strenght. */
						/*this is bilinear interpolation. not good. need to do averaging */
						const real iih=(ih-hpmin)*dhp1;
						const int iihf=ifloor(iih);
						const real iihw=iih-iihf;
						if(iihf<0||iihf>nhp-2){
							P(etf,icomp)=0.;
						} else{
							real tmp=pna[illt][iihf]*(1.-iihw)+pna[illt][iihf+1]*iihw;
							/*neglected rsa1 due to renormalization. */
							P(etf,icomp)=tmp;
							etf2sum+=tmp;
						}
					}
					if(fabs(etf2sum)>1.e-20){
						/*
							Changed: We normalize the etf by sum(profile), so we can model
							the variation of the intensity and meteor trails.
						*/
						cscale(etf, i0scale[illt]/etf2sum);

						//check for truncation
						const real ratio_edge=0.5*creal(P(etf,etf0)+P(etf,etf1-1))/creal(P(etf,netf2));
						if(ratio_edge>0.1){
							ncrop++;
							if(ratio_edge>max_crop) max_crop=ratio_edge;
						}
							
						cfftshift(etf);/*put peak in corner; */
						cfft2(etf, -1);

						/*Rotate the ETF. */
						/*need to put peak in center. */
						cfftshift(etf);
						const real theta=P(PR(srot, illt, 0), isa);
						real ct=cos(theta);
						real st=sin(theta);
						cmat* etf2d=P(petf, isa, illt);
						for(int icompy=0; icompy<notfy; icompy++){
							real iy=(icompy-notfy2);
							for(int icompx=0; icompx<notfx; icompx++){
								real ix=(icompx-notfx2);
								real ir=(dusc*(ct*ix+st*iy))+netf2;/*index in etf */
								int iir=ifloor(ir);
								ir=ir-iir;
								if(iir>=0&&iir<netf-1){
									/*bilinear interpolation. */
									P(etf2d, icompx, icompy)=P(etf,iir)*(1.-ir)+P(etf,iir+1)*ir;
								}/*else{P(etf2d,icompx,icompy)=0;}*/
							}
						}
						cfftshift(P(petf, isa, illt));/*peak in corner; */
					} else{
						warning_once("Wrong focus!\n");
						cset(P(petf, isa, illt), 1);
					}
					
				}//for isa
			}//for illt.
			ccellfree(etf_cache);
			dfree(thetas);
			if(ncrop){
				warning("Sodium profile is cropped by up to %.2f%% for %d subaps. Increase powfs.pixpsa or powfs.notf.\n", max_crop*100, ncrop);
			}
		}//if na_interp else
		//fuse nominal to etf to avoid multiply again.
		for(int illt=0; illt<nllt; illt++){
			for(int isa=0; isa<nsa; isa++){
				ccwm(P(petf, isa, illt), PR(dtfs[iwvl].nominal, isa, illt));
			}
		}
	}//for iwvl
	return etfs;
}

void dtf_free_do(dtf_t* dtfs){
	if(!dtfs) return;
	for(int iwvl=0;iwvl<dtfs->nwvl;iwvl++){
		ccellfree(dtfs[iwvl].nominal);
		dspcellfree(dtfs[iwvl].si);
		cfree(dtfs[iwvl].Ux);
		cfree(dtfs[iwvl].Uy);
	}
	free(dtfs);
}

void etf_free_do(etf_t* etfs){
	if(!etfs) return;
	for(int iwvl=0;iwvl<etfs->nwvl;iwvl++){
		ccellfree(etfs[iwvl].etf);
	}
	free(etfs);
}
/**
   profile is nxm array with m>=2. First column is coordinate and the remaining
   is function of coordinate. This routine smoothes the profile on new
   coordinate with spacing of dxnew.
 */
dmat* smooth(const dmat* prof, real dxnew){
	const long nxin=NX(prof);
	const real x0in=P(prof,0);
	const real dxin=(P(prof,nxin-1)-x0in)/(nxin-1);
	dmat* out;
	if(dxnew>dxin*2){
		dbg("Smoothing sodium profile from %g to %g sampling\n", dxin, dxnew);
		const long nxnew=ceil((P(prof,nxin-1)-x0in)/dxnew);
		loc_t* loc_in=mk1dloc_vec(P(prof), nxin);
		loc_t* loc_out=mk1dloc(x0in, dxnew, nxnew);
		dsp* ht=mkhb(loc_out, loc_in, 0, 0, 1);
		out=dnew(nxnew, NY(prof));
		memcpy(P(out), loc_out->locx, sizeof(real)*nxnew);
		OMP_TASK_FOR(4)
		for(long icol=1; icol<NY(prof); icol++){
			/*input profile */
			const real* pin=PCOL(prof, icol);
			/*output profile */
			real* pout=PCOL(out, icol);
			/*preserve sum of input profile */
			real Nasum=dvecsum(pin, nxin);
			dspmulvec(pout, ht, pin, 'n', 1);
			dnormalize_sumabs(pout, nxnew, Nasum);
		}
		dspfree(ht);
		locfree(loc_in);
		locfree(loc_out);
	} else{
		dbg("Not smoothing sodium profile from %g to %g\n", dxin, dxnew);
		out=dref(prof);
	}
	return out;
}
dcell *smooth_cell(const dcell *profs, real dxnew){
	dcell *out=dcellnew(NX(profs), NY(profs));
	for(int ix=0; ix<PN(profs); ix++){
		P(out,ix)=smooth(P(profs, ix), dxnew);
	}
	return out;
}

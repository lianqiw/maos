/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

DTF_T *mkdtf(const dmat *wvls, /**<List of wavelength*/
	     double dxsa,/**<Subaperture size*/
	     double embfac,/**<Embedding factor (2)*/
	     long ncompx,/**<FFT size along x*/
	     long ncompy,/**<FFT size along y*/
	     long pixpsax,/**<Number of pixels along x(r)*/
	     long pixpsay,/**<Number of pixels along y(a)*/
	     double pixthetax,/**<Pixel size along x (r)*/
	     double pixthetay,/**<Pixel size along y (a)*/
	     const dmat* pixoffx,  /**<offset of image center from center of pixel array, along x or radial*/
	     const dmat* pixoffy,  /**<offset of image center from center of pixel array, along y or azimuthal*/
	     double pixblur,  /**<Pixel blur sigma(fraction of pixel)*/
	     const dcell *srot, /**<Rotation angle of each subaperture. NULL for NGS WFS*/
	     int radpix,  /**<1: Pixels are along radial/azimuthal direction*/
	     int radrot  /**<For radial format CCD, rotate PSF/OTF into r/a coord. uses less memory*/
    ){
    int nwvl=wvls->nx;
    DTF_T *dtfs=mycalloc(nwvl,DTF_T);
    const double blurx=pixblur*pixthetax;
    const double blury=pixblur*pixthetay;
    const double e0x=-2*M_PI*M_PI*blurx*blurx;//blurring factors
    const double e0y=-2*M_PI*M_PI*blury*blury;
    const int do_blur=fabs(blurx)>EPS && fabs(blury)>EPS;
    const long ncompx2=ncompx>>1;
    const long ncompy2=ncompy>>1;
    const double pxo=-(pixpsax*0.5-0.5)*pixthetax;
    const double pyo=-(pixpsay*0.5-0.5)*pixthetay;
    double pxo2=pxo;
    double pyo2=pyo;
    if(pixoffx || pixoffy){
	int nwfs=1; int nsa=1;
	if(srot){
	    nwfs=srot->nx;
	    nsa=srot->p[0]->nx;
	}
	if(!pixoffx || !pixoffy){
	    error("pixoffx and pixoffy must simultaneously be supplied\n");
	}
	if(pixoffx->nx!=1 && pixoffx->nx!=nsa && pixoffx->ny!=1 && pixoffx->ny!=nwfs){
	    error("pixoffx has wrong format\n");
	}
	if(pixoffy->nx!=pixoffx->nx || pixoffy->ny!=pixoffx->ny){
	    error("pixoffy must have the same format as pixoffx\n");
	}
    }

    /*both nominal and si depends on wavelength.*/
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	const double wvl=wvls->p[iwvl];
	const double dtheta=wvl/(dxsa*embfac);/*PSF sampling. */
	const double dux=1./(dtheta*ncompx);
	const double duy=1./(dtheta*ncompy);
	const double dux2=dux*dux;
	const double duy2=duy*duy;
	const double pdtheta=pixthetax*pixthetay/(dtheta*dtheta);//scaling factor due to binning into detectors
	const double duxp=dux*pixthetax;
	const double duyp=duy*pixthetay;
	/*For LGS WFS (srot!=NULL), there is elongation effect. For radial pixel
	  CCD without using rotating psf/otf method, we need to create
	  nominal/si for each subaperture. The PSF/OTF and hence DTF are on x/y
	  coordinate, so pixels and pixel coordinates need to rotated to r/a
	  direction. */
	const int multi_dtf=(srot && radpix && !radrot);//DTF is along x/y coord, pixel is r/a coord.
	const int ndtf=(multi_dtf || (pixoffx && pixoffx->nx>1))?srot->p[0]->nx:1;
	const int multi_wfs=((srot && srot->nx>1 && radpix && !radrot) || (pixoffx && pixoffx->ny>1));
	const int nwfs=multi_wfs?MAX(srot->nx, pixoffx->ny):1;

	dtfs[iwvl].dtheta=dtheta;
	dtfs[iwvl].dxsa=dxsa;
	dtfs[iwvl].radpix=radpix;
	dtfs[iwvl].radrot=radrot;
	dtfs[iwvl].wvl=wvl;
	dtfs[iwvl].ncompx=ncompx;
	dtfs[iwvl].ncompy=ncompy;
	dtfs[iwvl].nominal=ccellnew(ndtf,nwfs);
	dtfs[iwvl].si=dspcellnew(ndtf,nwfs);
	ccell*  nominals=dtfs[iwvl].nominal;
	dspcell*  sis=dtfs[iwvl].si;
	cmat *nominal=cnew(ncompx,ncompy);
	//Coordinate of PSF points
	loc_t *loc_psf=mksqloc(ncompx,ncompy,dtheta,dtheta,-ncompx2*dtheta, -ncompy2*dtheta);
	double theta=0;
	double ct=cos(theta);
	double st=sin(theta);
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    for(int isa=0; isa<ndtf; isa++){
		if(multi_dtf){
		    theta=PR(srot, iwfs, 1)->p[isa];
		    ct=cos(theta);
		    st=sin(theta);
		}
		
		for(int iy=0; iy<ncompy; iy++){
		    int jy=iy-ncompy2;
		    for(int ix=0; ix<ncompx; ix++){
			int jx=ix-ncompx2;
			double ir=ct*jx+st*jy;
			double ia=-st*jx+ct*jy;
			//Pixel function
			P(nominal,ix,iy)=sinc(ir*duxp)*sinc(ia*duyp)*pdtheta;
			if(do_blur){//Charge diffusion.
			    P(nominal,ix,iy)*=exp(e0x*(ir*ir*dux2)+e0y*(ia*ia*duy2));
			}
		    }
		}
		/*put peak in corner. */
		cfftshift(nominal);
		/*Embed a fftshift operation in nominal so that we avoid fftshift when forming i0.*/
		cfft2(nominal,-1);
		cfftshift(nominal);
		cfft2(nominal,1);
		//cancel FFT scaling effect.
		cscale(nominal,1./(double)(nominal->nx*nominal->ny));
		ccp(PP(nominals,isa,iwfs), nominal);
		//Coordinate of pixels
		if(pixoffx){
		    double dx=PR(pixoffx, isa, iwfs)*pixthetax;
		    double dy=PR(pixoffy, isa, iwfs)*pixthetay;
		    pxo2=pxo-dx;
		    pyo2=pyo-dy;
		}
		loc_t *loc_ccd=mksqlocrot(pixpsax,pixpsay, pixthetax,pixthetay,pxo2,pyo2,theta);
		P(sis,isa,iwfs)=mkh(loc_psf,loc_ccd,0,0,1);
		locfree(loc_ccd);
	    }/*isa */
	}/*iwfs */
	cfree(nominal);
	locfree(loc_psf);

	/*Create an excessive high frequency in nominal so that
	  we don't have to do fftshift later.*/
	dtfs[iwvl].Ux=cnew(ncompx,1);
	dtfs[iwvl].Uy=cnew(ncompy,1);
	dcomplex *Ux=dtfs[iwvl].Ux->p;
	dcomplex *Uy=dtfs[iwvl].Uy->p;

	/*The following is used in genseotf to compute shifted i0.*/
	for(int ix=0; ix<ncompx; ix++){
	    int jx=ix<ncompx2?ix:(ix-ncompx);
	    Ux[ix]=COMPLEX(0, -2.*M_PI*jx*dux);
	}
	for(int iy=0; iy<ncompy; iy++){
	    int jy=iy<ncompy2?iy:(iy-ncompy);
	    Uy[iy]=COMPLEX(0, -2.*M_PI*jy*duy);
	}

    }/*iwvl */
    return dtfs;
}

ETF_T *mketf(DTF_T *dtfs,  /**<The dtfs*/
	     double hs,    /**<Guide star focus range*/
	     const dcell *sodium,/**<The sodium profile. In each cell First column is coordinate.*/
	     int icol,     /**<Which sodium profile to use*/
	     int nwvl,     /**<Number of wavelength*/
	     const dcell *srot,  /**<Rotation angle of each subaperture. NULL for NGS WFS*/
	     const dcell *srsa,  /**<Subaperture to LLT distance*/
	     double za,    /**<Zenith angle*/
	     int no_interp /**<Use direct sum instead of interpolation + FFT. Slower */
    ){
    ETF_T *etfs=mycalloc(nwvl,ETF_T);
    /*setup elongation along radial direction. don't care azimuthal. */
    if(!srot) error("srot is required");
    const int nllt=MAX(srot->nx, sodium->nx);
    const int nsa=srot->p[0]->nx;
    dmat *sodium0=sodium->p[0];
    const int ncol=(sodium0->ny-1);
    icol=wrap(icol, ncol);
    info("Na using column %d.\n",icol);
    const int nhp=sodium0->nx; 
    //adjusting sodium height for the zenith angle;
    double hpmin=0, dhp1=0;
    if(!no_interp){/**Requires linear spacing*/
	hpmin=sodium0->p[0]/cos(za);
	dhp1=cos(za)/(sodium0->p[1]-sodium0->p[0]);
	/*assume linear spacing. check the assumption valid */
	if(fabs(sodium0->p[nhp-1]-sodium0->p[0]-(nhp-1)*(sodium0->p[1]-sodium0->p[0]))>1.e-7){
	    error("llt profile is not evenly spaced\n");
	}
    }//else: Any spacing is ok.
    const double* px=sodium0->p;
    /*the sum of pp determines the scaling of the pixel intensity. */
    double *psrot[nllt];
    double *pna[nllt];
    double i0scale[nllt];
    for(int illt=0; illt<nllt; illt++){
	psrot[illt]=srot->p[srot->nx>1?illt:0]->p;
	pna[illt]=sodium->p[sodium->nx>1?illt:0]->p+nhp*(1+icol);
	i0scale[illt]=dblsum(pna[illt], nhp);
	if(fabs(i0scale[illt]-1)>0.01){
	    warning("Siglev is scaled by %g by sodium profile\n", i0scale[illt]);
	}
    }
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	const double dtheta=dtfs[iwvl].dtheta;
	const long ncompx=dtfs[iwvl].ncompx;
	const long ncompy=dtfs[iwvl].ncompy;
	const long ncompx2=ncompx>>1;
	const long ncompy2=ncompy>>1;
	const double dux=1./(dtheta*ncompx);
	const double duy=1./(dtheta*ncompy);
	ccell *petf=0;
	int use1d;
	if(dtfs[iwvl].radrot){
	    if(!dtfs[iwvl].radpix){
		error("radrot can only be used with radpix\n");
	    }
	    /*
	      PSF/OTF is defined in a/r coordinate. ETF is 1D only, along r.
	      Good for off-axis launch, otherwise, ETF and DTF takes a lot of
	      space for many LGS.
	    */
	    warning("Rotate PSF to do radial format detector (preferred)\n");
	    etfs[iwvl].p1=ccellnew(nsa,nllt);
	    petf=etfs[iwvl].p1;
	    use1d=1;
	}else{
	    /*
	      PSF/OTF/ETF is defined in x/y coordinate. 
	      2010-01-04: Fuse dtf nominal into etf for this case.
	    */
	    if(dtfs[iwvl].radpix){
		info_once("2D ETF for Radial CCD\n");
	    }else{
		info_once("Non-Radial CCD\n");
	    }
	    etfs[iwvl].p2=ccellnew(nsa,nllt);
	    petf=etfs[iwvl].p2;
	    use1d=0;
	}

	if(no_interp){ /* No interpolation, no fft. intensity scaling is automatically taken care of */
	    TIC;tic;
	    for(int illt=0; illt<nllt; illt++){
		for(int isa=0; isa<nsa; isa++){
		    double rsa=srsa->p[illt]->p[isa];
		    double rsa_za=rsa*cos(za);
		    if(use1d){ /*1d ETF along radius. */
			P(petf,isa,illt)=cnew(ncompx,1);
			dcomplex *etf1d=P(petf,isa,illt)->p;
#pragma omp parallel for default(shared)
			for(int icompx=0; icompx<ncompx; icompx++){
			    const double kr=dux*(icompx>=ncompx2?(icompx-ncompx):icompx);
			    for(int ih=0; ih<nhp; ih++){
				const double tmp=(-2*M_PI*(kr*(rsa_za/sodium0->p[ih]-rsa/hs)));
				etf1d[icompx]+=COMPLEX(pna[illt][ih]*cos(tmp), pna[illt][ih]*sin(tmp));
			    }
			}
		    }else{
			const double theta=psrot[illt][isa];
			const double ct=cos(theta);
			const double st=sin(theta);
			P(petf,isa,illt)=cnew(ncompx,ncompy);
			cmat *etf2d=P(petf,isa,illt);
#pragma omp parallel for default(shared)
			for(int icompy=0; icompy<ncompy; icompy++){
			    const double ky=duy*(icompy>=ncompy2?(icompy-ncompy):icompy);
			    for(int icompx=0; icompx<ncompx; icompx++){
				const double kx=dux*(icompx>=ncompx2?(icompx-ncompx):icompx);
				const double kr=(ct*kx+st*ky);/*along radial*/
				for(int ih=0; ih<nhp; ih++){
				    const double tmp=(-2*M_PI*(kr*(rsa_za/px[ih]-rsa/hs)));
				    P(etf2d,icompx,icompy)+=COMPLEX(pna[illt][ih]*cos(tmp), pna[illt][ih]*sin(tmp));
				}
			    }
			}
		    }
		}//isa
	    }//illt
	    toc2("ETF");
	}else{
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
	    const int netf=ncompx*nover*npad;
	    const double dtetf=dtheta/nover;
	    const double dusc=(netf*dtetf)/(dtheta*ncompx);
	    cmat *etf=cnew(netf,1);
	    double *thetas=mycalloc(netf,double);
	    const int netf2=netf>>1;
	    /*Only interpolating the center part. the rest is padding. */
	    const int etf0=netf2-(int)round(ncompx2*(dtheta/dtetf));
	    const int etf1=etf0+(int)round(ncompx*(dtheta/dtetf));
	    if(etf0<0) error("Invalid configuration\n");
	    for(int it=etf0; it<etf1; it++){
		thetas[it]=(it-netf2)*dtetf;
	    }
	    //cfft2plan(etf, -1);

	    for(int illt=0; illt<nllt; illt++){
		for(int isa=0; isa<nsa; isa++){
		    /*1d ETF along radius. */
		    double rsa=srsa->p[illt]->p[isa];
		    double etf2sum=0;
		    czero(etf);
		    for(int icomp=etf0; icomp<etf1; icomp++){
			/*peak in center */
			const double itheta=thetas[icomp];
			/*non linear mapping. changed from - to + on 2014-05-07.*/
			const double ih=hs*rsa/(rsa+itheta*hs);
			/*interpolating to get Na profile strenght. */
			/*this is bilinear interpolation. not good. need to do averaging */
			const double iih=(ih-hpmin)*dhp1;
			const int iihf=ifloor(iih);
			const double iihw=iih-iihf;
			if(iihf<0 || iihf>nhp-2){
			    etf->p[icomp]=0.;
			}else{
			    double tmp=pna[illt][iihf]*(1.-iihw)+pna[illt][iihf+1]*iihw;
			    /*neglected rsa1 due to renormalization. */
			    etf->p[icomp]=tmp;
			    etf2sum+=tmp;
			}
		    }
		    if(fabs(etf2sum)>1.e-20){
			/*2010-11-09:

			  We used to normalize the etf before fft so that after fft
			  it max to 1. The strength of original profile doesn't
			  matter.
		    
			  Changed: We no longer normalize the etf, so we can model
			  the variation of the intensity and meteor trails.
		      
			*/
			cscale(etf,i0scale[illt]/etf2sum);
			cfftshift(etf);/*put peak in corner; */
			cfft2(etf, -1);
			if(use1d){
			    if(npad==1 && nover==1){
				ccp(PP(petf,isa,illt),etf);
			    }else{
				cfftshift(etf);
				P(petf,isa,illt)=cnew(ncompx,1);
				dcomplex *etf1d=P(petf,isa,illt)->p;
#pragma omp parallel for default(shared)
				for(int icompx=0; icompx<ncompx; icompx++){
				    double ir=dusc*(icompx-ncompx2)+netf2;
				    int iir=ifloor(ir);
				    ir=ir-iir;
				    if(iir>=0 && iir<netf-1){
					etf1d[icompx]=etf->p[iir]*(1.-ir)
					    +etf->p[iir+1]*ir;
				    }/*else{etf1d[icompx]=0;}*/
				}
				cfftshift(P(petf,isa,illt));
			    }
			}else{
			    /*Rotate the ETF. */
			    /*need to put peak in center. */
			    cfftshift(etf);
			    double theta=psrot[illt][isa];
			    double ct=cos(theta);
			    double st=sin(theta);
			    P(petf,isa,illt)=cnew(ncompx,ncompy);
			    cmat *etf2d=P(petf,isa,illt);
#pragma omp parallel for default(shared)
			    for(int icompy=0; icompy<ncompy; icompy++){
				double iy=(icompy-ncompy2);
				for(int icompx=0; icompx<ncompx; icompx++){
				    double ix=(icompx-ncompx2);
				    double ir=(dusc*(ct*ix+st*iy))+netf2;/*index in etf */
				    int iir=ifloor(ir);
				    ir=ir-iir;
				    if(iir>=0 && iir<netf-1){
					/*bilinear interpolation. */
					P(etf2d,icompx,icompy)=etf->p[iir]*(1.-ir)
					    +etf->p[iir+1]*ir;
				    }/*else{P(etf2d,icompx,icompy)=0;}*/
				}
			    }
			    cfftshift(P(petf,isa,illt));/*peak in corner; */
			}
		    }else{
			warning_once("Wrong focus!\n");
			if(use1d){
			    P(petf,isa,illt)=cnew(ncompx,1);
			}else{
			    P(petf,isa,illt)=cnew(ncompx,ncompy);
			}
			cset(P(petf,isa,illt),1);
		    }
		}//for isa
	    }//for illt.
	    cfree(etf);
	    free(thetas);
	}//if na_interp
	if(!use1d){//fuse nominal to etf to avoid multiply again.
	    ccell*  pnominal=dtfs[iwvl].nominal;
	    int mnominal=0;/*multiply with isa to get index into pnominal. */
	    if(dtfs[iwvl].nominal->nx>1){
		mnominal=1;//polar ccd
	    }
	    int mllt=0;
	    if(dtfs[iwvl].nominal->ny>1){
		mllt=1;
	    }
	    for(int illt=0; illt<nllt; illt++){
		for(int isa=0; isa<nsa; isa++){
		    ccwm(P(petf,isa,illt), P(pnominal,isa*mnominal,illt*mllt));
		}
	    }
	    dtfs[iwvl].fused=1;
	}
    }//for iwvl
    return etfs;
}

void dtf_free_do(DTF_T *dtfs, int nwvl){
    if(!dtfs) return;
    for(int iwvl=0;iwvl<nwvl;iwvl++){
	ccellfree(dtfs[iwvl].nominal);
	dspcellfree(dtfs[iwvl].si);
	cfree(dtfs[iwvl].Ux);
	cfree(dtfs[iwvl].Uy);
    }
    free(dtfs);
}

void etf_free_do(ETF_T *etfs, int nwvl){
    if(!etfs) return;
    for(int iwvl=0;iwvl<nwvl;iwvl++){
	ccellfree(etfs[iwvl].p1);
	ccellfree(etfs[iwvl].p2);
    }
    free(etfs);
}
/**
   profile is nxm array with m>=2. First column is coordinate and the remaining
   is function of coordinate. This routine smoothes the profile on new
   coordinate with spacing of dxnew.
 */
dmat* smooth(const dmat *prof, double dxnew){
    const long nxin=prof->nx;
    const double x0in=prof->p[0];
    const double dxin=(prof->p[nxin-1]-x0in)/(nxin-1);
    dmat *out;
    if(dxnew > dxin * 2){
	const long nxnew=ceil((prof->p[nxin-1]-x0in)/dxnew);
	loc_t *loc_in=mk1dloc_vec(prof->p, nxin);
	loc_t *loc_out=mk1dloc(x0in, dxnew, nxnew);
	dsp *ht=mkhb(loc_out, loc_in, 0, 0, 1);
	out=dnew(nxnew, prof->ny);
	memcpy(out->p, loc_out->locx, sizeof(double)*nxnew);
#pragma omp parallel for default(shared)
	for(long icol=1; icol<prof->ny; icol++){
	    /*input profile */
	    const double *pin=PCOL(prof,icol);
	    /*output profile */
	    double *pout=PCOL(out, icol);
	    /*preserve sum of input profile */
	    double Nasum=dblsum(pin, nxin);
	    dspmulvec(pout, ht, pin, 'n', 1);
	    normalize_sumabs(pout, nxnew, Nasum);
	}
	dspfree(ht);
	locfree(loc_in);
	locfree(loc_out);
    }else{
	out=dref(prof);
    }
    return out;
}

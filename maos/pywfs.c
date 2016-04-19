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

/**
   \file pywfs.c
   Routines related to pyramid WFS
*/
#include "common.h"
#include "setup_powfs.h"
#include "pywfs.h"
#include "../cuda/gpu.h"
#define PYWFS_GUARD 1.5 //separate the pupil by 1.1 times more
/**
   Setup pyramid WFS based on configuration.
*/
void pywfs_setup(POWFS_T *powfs, const PARMS_T *parms, APER_T *aper, int ipowfs){
    pywfs_free(powfs[ipowfs].pywfs);
    PYWFS_T *pywfs=powfs[ipowfs].pywfs=calloc(1, sizeof(PYWFS_T));
    map_t *map=0;
    pywfs->hs=parms->powfs[ipowfs].hs;
    pywfs->poke=parms->recon.poke;//How many meters to poke
    if(pywfs->poke>1e-5 || pywfs->poke<1e-10){
	warning("poke=%g m is out of range\n", pywfs->poke);
    }
    pywfs->iwfs0=parms->powfs[ipowfs].wfs->p[0];
    double dx=parms->powfs[ipowfs].dx; 
    create_metapupil(&map, 0, 0, parms->dirs, parms->aper.d, 0, dx, dx, 0, 0, 0, 0, 0, 0);
    pywfs->loc=map2loc(map, 0);
    powfs[ipowfs].loc=pywfs->loc;//do not free here.
    mapfree(map);
    pywfs->amp=mkamp(powfs[ipowfs].loc, aper->ampground, 
		     parms->misreg.pupil->p[0],parms->misreg.pupil->p[1], 
		     parms->aper.d, parms->aper.din);
    powfs[ipowfs].amp=dref(pywfs->amp);
    setup_powfs_misreg(powfs, parms, aper, ipowfs);
    powfs[ipowfs].realamp=dcellnew(parms->powfs[ipowfs].nwfs,1);
    for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
	if(powfs[ipowfs].amp_tel){
	    powfs[ipowfs].realamp->p[jwfs]=dref(powfs[ipowfs].amp_tel->p[jwfs]);
	}else{
	    powfs[ipowfs].realamp->p[jwfs]=dref(powfs[ipowfs].amp); 
	}
    }
    int nwvl=parms->powfs[ipowfs].wvl->nx;
    double oversize=2*PYWFS_GUARD;
    pywfs->locfft=locfft_init(powfs[ipowfs].loc, pywfs->amp, parms->powfs[ipowfs].wvl, 0, oversize, 0);
    pywfs->wvlwts=ddup(parms->powfs[ipowfs].wvlwts);
    pywfs->modulate=parms->powfs[ipowfs].modulate;
    pywfs->modulpos=pywfs->modulate>0?parms->powfs[ipowfs].modulpos:1;
    long nembed=pywfs->locfft->nembed->p[0];
    double wvlmin, wvlmax;
    dmaxmin(parms->powfs[ipowfs].wvl->p, nwvl, &wvlmax, &wvlmin);
    double dtheta_min=wvlmin/(dx*nembed);
    pywfs->gain=1;
    //size of the part of the PSF captured by pyramid
    long ncomp=nembed;
    if(parms->powfs[ipowfs].fieldstop){
	ncomp=nextfftsize(ceil(parms->powfs[ipowfs].fieldstop/dtheta_min));
	if(ncomp>nembed){
	    ncomp=nembed;
	}
    }

    long ncomp2=ncomp/2;
    pywfs->pyramid=cellnew(nwvl, 1);
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	pywfs->pyramid->p[iwvl]=cnew(ncomp, ncomp);
	PCMAT(pywfs->pyramid->p[iwvl], pp);
	dcomplex coeff=COMPLEX(0, M_PI*0.5);
	long skip=0;
	if(parms->powfs[ipowfs].fieldstop){//Limit fov per wvl
	    double dtheta=parms->powfs[ipowfs].wvl->p[iwvl]/(dx*nembed);
	    int nstop=ceil(parms->powfs[ipowfs].fieldstop/dtheta*0.5)*2;
	    if(nstop>ncomp) nstop=ncomp;
	    skip=(ncomp-nstop)/2;
	}
	for(long iy=skip; iy<ncomp-skip; iy++){
	    for(long ix=skip; ix<ncomp-skip; ix++){
		pp[iy][ix]=cexp((labs(iy-ncomp2)+labs(ix-ncomp2))*coeff);
	    }
	}
    }

    cmat *nominal=pywfs->nominal=cnew(ncomp, ncomp);
    PCMAT(nominal, pn);
    long order=parms->powfs[ipowfs].order;
    double dsa=parms->aper.d/order;//size of detector pixel mapped on pupil
    double dx2=dx*nembed/ncomp;//sampling of pupil after inverse fft
    double du=1./(dx2*ncomp);
    double dupix=dsa*du;
    double pdmeter=pow(dsa/dx2, 2);
    double pixblur=parms->powfs[ipowfs].pixblur;
    double e0b=-2*pow(M_PI*pixblur*dsa*du, 2);
    for(int iy=0; iy<ncomp; iy++){
	int jy=iy-ncomp2;
	for(int ix=0; ix<ncomp; ix++){
	    int jx=ix-ncomp2; 
	    pn[iy][ix]=sinc(jy*dupix)*sinc(jx*dupix)*pdmeter;
	    if(pixblur){
		pn[iy][ix]*=exp(e0b*(jx*jx+jy*jy));
	    }
	}
    }
    cfftshift(nominal);
    cfft2(nominal, -1);
    cfftshift(nominal);
    cfft2(nominal, 1);
    cscale(nominal, 1./(nominal->nx*nominal->ny));
    if(parms->dbg.pwfs_psx){
	if(parms->dbg.pwfs_psx->nx!=4){
	    error("dbg.pwfs_psx has wrong format: expected 4x1, got %ldx%ld\n",
		  parms->dbg.pwfs_psx->nx, parms->dbg.pwfs_psx->ny);
	}
	if(parms->dbg.pwfs_psy->nx!=4){
	    error("dbg.pwfs_psy has wrong format: expected 4x1, got %ldx%ld\n",
		  parms->dbg.pwfs_psy->nx, parms->dbg.pwfs_psy->ny);
	}
	pywfs->pupilshift=dnew(4,2);
	for(int i=0; i<4; i++){
	    double tmp;
	    tmp=IND(parms->dbg.pwfs_psx, i);
	    IND(pywfs->pupilshift, i, 0)=tmp-round(tmp);
	    tmp=IND(parms->dbg.pwfs_psy, i);
	    IND(pywfs->pupilshift, i, 1)=tmp-round(tmp);

	}
	if(parms->save.setup){
	    writebin(pywfs->pupilshift, "powfs%d_pupilshift", ipowfs);
	}
    }
    pywfs->si=cellnew(4,1);//for each quadrant.
    //Make loc_t symmetric to ensure proper sampling onto detector. Center of subaperture
    powfs[ipowfs].saloc=mksqloc(order, order, dsa, dsa, 
				(-order*0.5+0.5)*dsa, (-order*0.5+0.5)*dsa);
    loc_t *loc_fft=mksqloc(ncomp, ncomp, dx2, dx2, (-ncomp2+0.5)*dx2, (-ncomp2+0.5)*dx2);
    for(int iy=0; iy<2; iy++){
	for(int ix=0; ix<2; ix++){
	    const int ind=ix+iy*2;
	    double shx=0, shy=0;
	    if(pywfs->pupilshift){
		shx=IND(pywfs->pupilshift, ind, 0)*dsa;
		shy=IND(pywfs->pupilshift, ind, 1)*dsa;
	    }
	    pywfs->si->p[ind]=mkh(loc_fft, powfs[ipowfs].saloc, 
				  ((ix-0.5)*ncomp2)*dx2+shx, 
				  ((iy-0.5)*ncomp2)*dx2+shy,
				  1); 
	}
    }
    if(parms->save.setup){
	writebin(pywfs->si, "powfs%d_si0", ipowfs);
	locwrite(pywfs->locfft->loc, "powfs%d_locfft", ipowfs);
    }
    {
	//Determine subapertures area
	dmat *opd=dnew(pywfs->locfft->loc->nloc, 1);
	dmat *ints=0;
	pywfs_fft(&ints, powfs[ipowfs].pywfs, opd);
	if(parms->save.setup){
	    writebin(ints, "powfs%d_ints0", ipowfs);
	}
	int nints=ints->nx;
	pywfs->saa=dnew(nints, 1);
	for(int i=0; i<ints->nx; i++){
	    pywfs->saa->p[i]=ints->p[i]+ints->p[i+nints]+ints->p[i+nints*2]+ints->p[i+nints*3];
	}
	cellfree(ints);
	dfree(opd);
    }
    if(parms->powfs[ipowfs].saat>0){
	dmat *saa=pywfs->saa;
	double samax=dmaxabs(saa);
	loc_reduce(powfs[ipowfs].saloc, saa, parms->powfs[ipowfs].saat*samax, 0, 0);
	dscale(saa, saa->nx/dsum(saa));//saa average to one.
	powfs[ipowfs].saa=dref(pywfs->saa);
	for(int iy=0; iy<2; iy++){
	    for(int ix=0; ix<2; ix++){
		dspfree(pywfs->si->p[ix+iy*2]);
		pywfs->si->p[ix+iy*2]=mkh(loc_fft, powfs[ipowfs].saloc, 
					  ((ix-0.5)*dx2*ncomp2), 
					  ((iy-0.5)*dx2*ncomp2),
					  1);	    
	    }
	}
    }
    const int nsa=powfs[ipowfs].saloc->nloc;
    dscale(pywfs->saa, pywfs->saa->nx/dsum(pywfs->saa));//saa average to one.
    locfree(loc_fft);
  
    //Determine the gain and offset of PyWFS
    {
	//offset: grad of a flat wavefront
	dmat *opd=dnew(pywfs->locfft->loc->nloc, 1);
	dmat *ints=0;
	dmat *goff=0;
	pywfs_fft(&ints, pywfs, opd);//writebin(ints, "ints_0");
	pywfs_grad(&goff, pywfs, ints);//writebin(goff, "goff_0");
	dadd(&pywfs->gradoff, 1, goff, 1);
	dfree(opd);
	dfree(ints);
	//gain
	dmat *TT=pywfs_tt(pywfs);
	double gxm=0, gym=0;
	for(int isa=0; isa<nsa; isa++){
	    gxm+=TT->p[isa];
	    gym+=TT->p[isa+nsa*3];
	}
	double gainscl=2.*nsa/(gxm+gym);
	//gainscl*=2;//inject an error;
	pywfs->gain*=gainscl;
	dscale(pywfs->gradoff, gainscl);
	dscale(TT, gainscl);
	pywfs->GTT=TT;
	info("pywfs_gain=%g\n", pywfs->gain);
    }

    //Determine the NEA. It will be changed by powfs.gradscale as dithering converges    
    {
	powfs[ipowfs].saneaxy=cellnew(nsa,1);
	double rne=parms->powfs[ipowfs].rne;
	for(int isa=0; isa<nsa; isa++){
	    dmat *tmp=dnew(2,2);
	    double ogi=pywfs->gain*parms->powfs[ipowfs].gradscale;
	    double sig=pywfs->saa->p[isa]*parms->powfs[ipowfs].siglev;//siglev of subaperture
	    IND(tmp,0,0)=IND(tmp,1,1)=pow(ogi/sig,2)*(sig+4*rne*rne);
	    powfs[ipowfs].saneaxy->p[isa]=tmp;
	}
    }
    if(parms->save.setup){
	writebin(powfs[ipowfs].loc, "powfs%d_loc", ipowfs);
	writebin(powfs[ipowfs].saloc, "powfs%d_saloc", ipowfs);
	writebin(pywfs->amp, "powfs%d_amp", ipowfs);
	writebin(pywfs->locfft->embed, "powfs%d_embed", ipowfs);
	writebin(pywfs->pyramid, "powfs%d_pyramid", ipowfs);
	writebin(nominal, "powfs%d_nominal", ipowfs);
	writebin(pywfs->si, "powfs%d_si", ipowfs);
	writebin(pywfs->gradoff, "powfs%d_gradoff", ipowfs);
	writebin(powfs[ipowfs].saneaxy, "powfs%d_sanea", ipowfs);
    }
    if(0){//Test implementation using zernikes
	dmat *ints=0;
	int nn=1;
	double wve=1e-9*160;
	dmat *opds=zernike(pywfs->locfft->loc, parms->aper.d, 0, 3, 0);
	cellarr *pupsave=cellarr_init(nn,opds->ny,"ints");
	cellarr *grads=cellarr_init(nn,opds->ny,"grads");
	dmat *opd=0;
	dmat *grad=0;
	for(int im=0; im<opds->ny; im++){
	    for(int j=0; j<nn; j++){
		info2("im=%d, j=%d\n", im, j);
		opd=dsub(opds, 0, 0, im, 1);
		dscale(opd, pow(2, j)*wve);
		dzero(ints);
		pywfs_fft(&ints, powfs[ipowfs].pywfs, opd);
		pywfs_grad(&grad, powfs[ipowfs].pywfs, ints);
		cellarr_push(pupsave, j+im*nn, ints);
		cellarr_push(grads, j+im*nn, grad);
		dfree(opd);
	    }
	}
	cellarr_close(pupsave);
	cellarr_close(grads);
	cellfree(ints);
	dfree(opds);
	exit(0);
    }
    if(0){//Test linearity of a zenike mode with noise
	double wve=1e-9*20;
	dmat *opds=zernike(pywfs->locfft->loc, parms->aper.d, 0, 0, -parms->powfs[ipowfs].dither);
	dmat *opdi=0;
	dmat *ints=0, *grad=0;
	dadd(&opdi, 0, opds, wve);
	pywfs_fft(&ints, powfs[ipowfs].pywfs, opdi);
	pywfs_grad(&grad, powfs[ipowfs].pywfs, ints);
	dmat *reg=dpinv(grad, 0);
	writebin(grad, "dither_grad");
	writebin(reg, "dither_reg");
	writebin(ints, "dither_ints");
	rand_t rstat;
	seed_rand(&rstat, 1);
	dmat *tmp=0;
	int nj=10, nn=10;
	dmat *res=dnew(nj, nn);
	dmat *ints2=0;
	for(int j=0; j<nj; j++){
	    dzero(ints);
	    dadd(&opdi, 0, opds, wve*(j+1));
	    pywfs_fft(&ints, powfs[ipowfs].pywfs, opdi);
	    for (int in=0; in<nn; in++){
		dadd(&ints2, 0, ints, 100);
		addnoise(ints2, &rstat, 0, 0, 0, 0, in);
		pywfs_grad(&grad, powfs[ipowfs].pywfs, ints2);
		dmm(&tmp, 0, reg, grad, "nn", 1);
		IND(res, j, in)=tmp->p[0];
		info2("%d of %d, %d of %d: %g\n", j, nj, in, nn, tmp->p[0]);
	    }
	}
	writebin(opds, "dither_opd");
	writebin(res, "dither_res");
	dfree(opds); dfree(opdi); dfree(ints); dfree(grad); dfree(reg), dfree(tmp); dfree(res); dfree(ints2);
	exit(0);
    }
    //Test NCPA calibration
    int PYWFS_NCPA=0;
    READ_ENV_INT(PYWFS_NCPA,0,1);
    if(PYWFS_NCPA){
	dmat *opdatm=dread("opdatm");
	dmat *opdbias_full=dread("opdbias_full");
	dmat *opdbias_astigx=dread("opdbias_astigx");
	dmat *opdbias_polish=dread("opdbias_polish");
	const double atmscale=1;
#pragma omp parallel for
	for(int i=0; i<100; i++){
	    info2("%d ", i);
	    dmat *ints=0;
	    dmat *grad=0;
	    dmat *opd=0;

	    dadd(&opd, 0, opdatm, (i+1)*0.02);
	    dzero(ints);
	    pywfs_fft(&ints, powfs[ipowfs].pywfs, opd);
	    pywfs_grad(&grad, powfs[ipowfs].pywfs, ints);
	    writebin(grad, "grad_atm_%d", i);

	    dadd(&opd, 0, opdbias_full, (i+1)*0.02);
	    dzero(ints);
	    pywfs_fft(&ints, powfs[ipowfs].pywfs, opd);
	    pywfs_grad(&grad, powfs[ipowfs].pywfs, ints);
	    writebin(grad, "gradbias_full_%d", i);

	    dadd(&opd, 0, opdbias_astigx, (i+1)*0.02);
	    dzero(ints);
	    pywfs_fft(&ints, powfs[ipowfs].pywfs, opd);
	    pywfs_grad(&grad, powfs[ipowfs].pywfs, ints);
	    writebin(grad, "gradbias_astigx_%d", i);

	    dadd(&opd, 0, opdbias_polish, (i+1)*0.02);
	    dzero(ints);
	    pywfs_fft(&ints, powfs[ipowfs].pywfs, opd);
	    pywfs_grad(&grad, powfs[ipowfs].pywfs, ints);
	    writebin(grad, "gradbias_polish_%d", i);

	    dadd(&opd, 0, opdbias_full, (i+1)*0.02);
	    dadd(&opd, 1, opdatm, atmscale);
	    dzero(ints);
	    pywfs_fft(&ints, powfs[ipowfs].pywfs, opd);
	    pywfs_grad(&grad, powfs[ipowfs].pywfs, ints);
	    writebin(grad, "gradboth_full_%d", i);

	    
	    dadd(&opd, 0, opdbias_astigx, (i+1)*0.02);
	    dadd(&opd, 1, opdatm, atmscale);
	    dzero(ints);
	    pywfs_fft(&ints, powfs[ipowfs].pywfs, opd);
	    pywfs_grad(&grad, powfs[ipowfs].pywfs, ints);
	    writebin(grad, "gradboth_astigx_%d", i);

	    dadd(&opd, 0, opdbias_polish, (i+1)*0.02);
	    dadd(&opd, 1, opdatm, atmscale);
	    dzero(ints);
	    pywfs_fft(&ints, powfs[ipowfs].pywfs, opd);
	    pywfs_grad(&grad, powfs[ipowfs].pywfs, ints);
	    writebin(grad, "gradboth_polish_%d", i);
	    dfree(opd);
	    dfree(ints);
	    dfree(grad);

	    
	    dadd(&opd, 0, opdbias_full, (i+1)*0.02);
	    dadd(&opd, 1, opdatm, (i+1)*0.02);
	    dzero(ints);
	    pywfs_fft(&ints, powfs[ipowfs].pywfs, opd);
	    pywfs_grad(&grad, powfs[ipowfs].pywfs, ints);
	    writebin(grad, "gradall_%d", i);
	    dfree(opd);
	    dfree(ints);
	    dfree(grad);

	}

	exit(0);	
    }
}
/**
   Perform FFT over the complex PSF with additional phases caused by the
   pyramid. FFT on each quadrant of the PSF creates diffraction effects.
*/
void pywfs_fft(dmat **ints, const PYWFS_T *pywfs, const dmat *opd){
    locfft_t *locfft=pywfs->locfft;
    ccell *psfs=0;
    locfft_psf(&psfs, locfft, opd, NULL, 1);//psfs.^2 sum to 1. peak in center
    int nwvl=locfft->wvl->nx;
    double dx=locfft->loc->dx;
    long nembed=locfft->nembed->p[0];
    long nembed2=nembed/2;
    dmat *wvlwts=pywfs->wvlwts;
    //position of pyramid for modulation
    int pos_n=pywfs->modulpos;
    double pos_r=pywfs->modulate;
    if(pos_r<=0) pos_n=1;
    long ncomp=pywfs->nominal->nx;
    long ncomp2=ncomp/2;
    cmat *otf=cnew(ncomp, ncomp);
    dmat *pupraw=dnew(ncomp, ncomp);
    //writebin(psfs, "cpu_wvf0");
    for(int ipos=0; ipos<pos_n; ipos++){
	//whether the first point falls on the edge or not makes little difference
	double theta=2*M_PI*(ipos+0.)/pos_n;
	double posx=cos(theta)*pos_r;
	double posy=sin(theta)*pos_r;
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    double dtheta=locfft->wvl->p[iwvl]/(dx*nembed);
	    long offy=(long)round(posy/dtheta);
	    long offy2=nembed2+offy-ncomp2;
	    long iy0=MAX(-offy2, 0);
	    long ny2=MIN(ncomp, nembed-offy2)-iy0;

	    long offx=(long)round(posx/dtheta);
	    long offx2=nembed/2+offx-ncomp2;
	    long ix0=MAX(-offx2, 0);
	    long nx2=MIN(ncomp, nembed-offx2)-ix0;

	    czero(otf);
	    dcomplex *pyramid=pywfs->pyramid->p[iwvl]->p;
	    for(long iy=iy0; iy<ny2; iy++){
		for(long ix=ix0; ix<nx2; ix++){
		    long indin=ix+offx2+(iy+offy2)*nembed;
		    long indout=ix+iy*ncomp;
		    otf->p[indout]=psfs->p[iwvl]->p[indin]*pyramid[indout];
		}
	    }
	    //writebin(otf, "cpu_wvf1");
	    cfft2(otf, 1);
	    //writebin(otf, "cpu_wvf2");
	    cabs22d(&pupraw, 1., otf, wvlwts->p[iwvl]/(ncomp*ncomp*pos_n));
	    //writebin(pupraw, "cpu_wvf3");
	}//for iwvl
    }//for ipos
    //writebin(pupraw, "cpu_psf"); exit(0);
    ccpd(&otf, pupraw);//pupraw sum to one.
    //writebin(otf, "cpu_wvf4");
    dfree(pupraw);
    cfft2(otf, -1);
    //writebin(otf, "cpu_wvf5");
    ccwm(otf, pywfs->nominal);
    cfft2(otf, 1);
    //writebin(otf, "cpu_wvf6");
    const int nsa=pywfs->si->p[0]->nx;
    if(!(*ints)){
	(*ints)=dnew(nsa, 4);
    }
    for(int i=0; i<4; i++){
	//normalized so that each "subaperture" sum to 1.
	dspmulcreal((*ints)->p+nsa*i, pywfs->si->p[i], otf->p, (double)nsa/(ncomp*ncomp));
    }
    //writebin(*ints, "cpu_ints"); exit(0);
    ccellfree(psfs);
    cfree(otf);
}
/**
   Compute gradients. It replaces the result, not accumulate.
 */
void pywfs_grad(dmat **pgrad, const PYWFS_T *pywfs, const dmat *ints){
    const long nsa=ints->nx;
    if(!*pgrad){
	*pgrad=dnew(nsa*2,1);
    }
    double *pgx=(*pgrad)->p;
    double *pgy=(*pgrad)->p+nsa;
    PDMAT(ints, pi);
    double gain=pywfs->gain;
    if(0){//Use standard Quad Cell algorithm
	warning_once("Do not use mean i0\n");
	for(int isa=0; isa<nsa; isa++){
	    double alpha2=gain/(pi[0][isa]+pi[1][isa]+pi[2][isa]+pi[3][isa]);
	    pgx[isa]=(pi[1][isa]-pi[0][isa]
		      +pi[3][isa]-pi[2][isa])*alpha2;
	    pgy[isa]=(pi[2][isa]+pi[3][isa]
		      -pi[0][isa]-pi[1][isa])*alpha2;
	}
    }else{//Denominator is replaced by SAA*mean(i0).
	double imean=dsum(ints)/nsa;
	double alpha0=gain/imean;
	for(int isa=0; isa<nsa; isa++){
	    double alpha2=alpha0/pywfs->saa->p[isa];
	    pgx[isa]=(pi[1][isa]-pi[0][isa]
		      +pi[3][isa]-pi[2][isa])*alpha2;
	    pgy[isa]=(pi[2][isa]+pi[3][isa]
		      -pi[0][isa]-pi[1][isa])*alpha2;
	}
    }
    if(pywfs->gradoff){
	dadd(pgrad, 1, pywfs->gradoff, -1);
    }
}
/**
   Return measurement of T/T mode, normalized for 1 unit of input.
*/
dmat *pywfs_tt(const PYWFS_T *pywfs){
    if(pywfs->GTT) {
	warning("Reusing cached pywfs->GTT\n");
	return dref(pywfs->GTT);
    }
    TIC;tic;info2("Computing pywfs_tt...");
    const loc_t *loc=pywfs->locfft->loc;
    dmat *opd=dnew(loc->nloc,1);
    dmat *ints=0;
    long nsa=pywfs->si->p[0]->nx;
    dmat *out=dnew(nsa*2,2);
    dmat *gradx=drefcols(out, 0, 1);
    dmat *grady=drefcols(out, 1, 1);
    dmat *gradx2=dnew(nsa*2,1);
    dmat *grady2=dnew(nsa*2,1);

    double ptt[3]={0,0,0};
    double alpha=0.005/206265.;

    //+x
    ptt[1]=alpha;  ptt[2]=0;
    loc_add_ptt(opd->p, ptt, loc);
    dzero(ints);
    pywfs_fft(&ints, pywfs, opd);
    pywfs_grad(&gradx, pywfs, ints);
    //+y
    ptt[1]=-alpha; ptt[2]=alpha;
    loc_add_ptt(opd->p, ptt, loc);
    dzero(ints);
    pywfs_fft(&ints, pywfs, opd);
    pywfs_grad(&grady, pywfs, ints);
    //-x
    ptt[1]=-alpha; ptt[2]=-alpha;
    loc_add_ptt(opd->p, ptt, loc);
    dzero(ints);
    pywfs_fft(&ints, pywfs, opd);
    pywfs_grad(&gradx2, pywfs, ints);
    //-y
    ptt[1]=+alpha; ptt[2]=-alpha;
    loc_add_ptt(opd->p, ptt, loc);
    dzero(ints);
    pywfs_fft(&ints, pywfs, opd);
    pywfs_grad(&grady2, pywfs, ints);

    dadd(&gradx, 1, gradx2, -1);
    dadd(&grady, 1, grady2, -1);
    dscale(out, 0.5/alpha);
    dfree(gradx);
    dfree(grady);
    dfree(gradx2);
    dfree(grady2);
    dfree(opd);
    dfree(ints);
    toc2("done");
    return out;
}
static uint32_t pywfs_hash(const PYWFS_T *pywfs, uint32_t key){
    key=lochash(pywfs->loc, key);
    key=dhash(pywfs->amp, key);
    key=dhash(pywfs->saa, key);
    key=dhash(pywfs->wvlwts, key);
    key=chash(pywfs->pyramid->p[0], key);
    if(pywfs->pupilshift){
	key=dhash(pywfs->pupilshift, key);
    }
    return key;
}
/**
   There is no need to simulate turbulence to fake optical gain. Any optical
   gain can be used as long as it "correct", i.e., a radian of tilt produces one
   radian of tilt, which is gauranteed when pywfs->gain is computed under the
   same conditions.
 */
static dmat *pywfs_mkg_do(const PYWFS_T *pywfs, const loc_t* locin, const dmat *mod, double displacex, double displacey, double scale){
    const loc_t *locfft=pywfs->locfft->loc;
    const int nsa=pywfs->si->p[0]->nx;
    dmat *grad0=dnew(nsa*2,1);
    {
	dmat *opd=dnew(locfft->nloc, 1);
	dmat *ints=0;
	pywfs_fft(&ints, pywfs, opd);
	pywfs_grad(&grad0, pywfs, ints);
	//writebin(grad0, "grad0_cpu");
	//writebin(ints, "ints0_cpu");
	dfree(opd);
	dfree(ints);
    }
    int count=0;
    int nmod=mod?mod->ny:locin->nloc;
    dmat *ggd=dnew(nsa*2, nmod);
    if(mod && mod->nx!=locin->nloc){
	error("mod->nx must equal to %ld", locin->nloc);
    }
    TIC;tic;
#pragma omp parallel for shared(count)
    for(int imod=0; imod<nmod; imod++){
	dmat *opdin=dnew(locin->nloc, 1);
	dmat *opdfft=dnew(locfft->nloc, 1);
	dmat *ints=0;
	dmat *grad=drefcols(ggd, imod, 1);
	double poke=pywfs->poke;
	if(mod){
	    dmat *tmp=drefcols(mod, imod, 1);
	    //the equivalent radimodl order of zernike.
	    //double radial=ceil((sqrt(8.*(imod+1)+1)-3)*0.5)+1;
	    //double std=dstd(tmp);
	    double tmax,tmin;
	    dmaxmin(tmp->p, tmp->nx, &tmax, &tmin);
	    poke/=(tmax-tmin);//sqrt(radial);
	    dadd(&opdin, 0, tmp, poke);
	    dfree(tmp);
	}else{
	    opdin->p[imod]=poke;
	}
	prop_nongrid((loc_t*)locin, opdin->p, locfft, opdfft->p, 1, displacex, displacey, scale, 0, 0);
	//writebin(opdfft, "phiout_cpu_%d", imod);
	pywfs_fft(&ints, pywfs, opdfft);
	//writebin(ints, "ints_cpu_%d", imod);
	pywfs_grad(&grad, pywfs, ints);
	dadd(&grad, 1, grad0, -1);
	dscale(grad, 1./poke);	
	atomicadd(&count, 1);
	if(count%10==0){
	    double ts=myclockd()-tk;
	    info2("%d of %ld. %.2f of %.2f seconds. std(grad)=%g.\n", count, locin->nloc, ts, ts/count*locin->nloc, dstd(grad));
	}
	dfree(opdfft);
	dfree(opdin);
	dfree(grad);
	dfree(ints);
    }
    info2("\n");
    dfree(grad0);	
    return ggd;
}
/**
   ploc is on pupil.
 */
dmat* pywfs_mkg(const PYWFS_T *pywfs, const loc_t* locin, const dmat *mod, double displacex, double displacey, double scale){
    if(mod && mod->ny<=6){
	return pywfs_mkg_do(pywfs, locin, mod, displacex, displacey, scale);
    }
    uint32_t key=0;
    key=lochash(locin, key);
    key=pywfs_hash(pywfs, key);
    if(mod) key=dhash(mod, key);
    char fn[PATH_MAX];
    char fnlock[PATH_MAX];
    mymkdir("%s/.aos/cache/", HOME);
    snprintf(fn, PATH_MAX, "%s/.aos/cache/G_%u_%ld_%ld_%g_%d_%g_%g_%g_%g_%g.bin", HOME, 
	     key, pywfs->locfft->nembed->p[0], locin->nloc, pywfs->modulate, pywfs->modulpos,
	     locin->iac, displacex, displacey, scale, pywfs->poke);
    snprintf(fnlock, PATH_MAX, "%s.lock", fn);
    info2("Using G in %s\n", fn);
    dmat *gg=0;
    if(0){//test amount of poke 
	dmat *mod1=dnew(locin->nloc, 3);
	if(mod){
	    memcpy(PCOL(mod1, 0), PCOL(mod, 0), sizeof(double)*mod->nx);
	    memcpy(PCOL(mod1, 1), PCOL(mod, 1000), sizeof(double)*mod->nx);
	    memcpy(PCOL(mod1, 2), PCOL(mod, 3200), sizeof(double)*mod->nx);
	}else{
	    IND(mod1, 30, 0)=1;//sqrt(locin->nloc);
	    IND(mod1, 800, 1)=1;//sqrt(locin->nloc);
	    IND(mod1, 3000, 2)=1;//sqrt(locin->nloc);
	}
	writebin(mod1, "mod1");
	dcell *gg1=dcellnew(13,1);
	dcell *gg2=dcellnew(13,1);
	double poke=1e-9;
	double step=pow(10,0.25);
	for(int ig=0; ig<gg1->nx; ig++){
	    ((PYWFS_T*)pywfs)->poke=poke;
	    gg1->p[ig]=gpu_pywfs_mkg(pywfs, locin, mod1, displacex, displacey);
	    gg2->p[ig]=pywfs_mkg_do(pywfs, locin, mod1, displacex, displacey, scale);
	    poke=poke*step;
	}
	writebin(gg1, "gg1g");
	cellfree(gg1);
	writebin(gg2, "gg1c");
	cellfree(gg2);
	exit(0);
    }
  retry:
    if(exist(fnlock) || !zfexist(fn)){
	int fd=lock_file(fnlock, 0, 0);//non blocking, exclusive
	if(fd>0){//succeed
	    info2("Generating PYWFS poke matrix\n");
#if USE_CUDA
	    if(global->parms->gpu.wfs){
		gg=gpu_pywfs_mkg(pywfs, locin, mod, displacex, displacey);
	    }else
#endif
		gg=pywfs_mkg_do(pywfs, locin, mod, displacex, displacey, scale);
	    writebin(gg, "%s", fn);
	    snprintf(fn, PATH_MAX, "%s/.aos/cache/", HOME);
	    remove_file_older(fn, 365*24*3600);//one year
	    close(fd); remove(fnlock);
	}else{
	    info2("Trying to lock %s\n", fnlock);
	    fd=lock_file(fnlock, 1, 0);
	    close(fd); remove(fnlock);
	    goto retry;
	}
    }else{
	gg=dread("%s", fn);
    }
    return gg;
}

void pywfs_free(PYWFS_T *pywfs){
    if(!pywfs) return;
    dfree(pywfs->amp);
    pywfs->loc=0;
    locfft_free(pywfs->locfft);
    cellfree(pywfs->pyramid);
    cfree(pywfs->nominal);
    dspcellfree(pywfs->si);
    free(pywfs);
}

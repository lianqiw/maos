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

/**
   \file pywfs.c
   Routines related to pyramid WFS
*/
#include "common.h"
#include "setup_powfs.h"
#include "pywfs.h"
#include "../cuda/gpu.h"
#define PYWFS_GUARD 1.5 //separate the pupil by 1.1 times more
#define PYWFS_POKE 1e-6 //How many meters to poke
/**
   Setup pyramid WFS based on configuration.
*/
void pywfs_setup(POWFS_T *powfs, const PARMS_T *parms, APER_T *aper, int ipowfs){
    pywfs_free(powfs[ipowfs].pywfs);
    PYWFS_T *pywfs=powfs[ipowfs].pywfs=calloc(1, sizeof(PYWFS_T));
    map_t *map=0;
    double dx=parms->powfs[ipowfs].dx; 
    create_metapupil(&map, 0, 0, parms->dirs, parms->aper.d, 0, dx, dx, 0, 0, 0, 0, 0, 0);
    powfs[ipowfs].loc=map2loc(map);
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
	ncomp=ceil(parms->powfs[ipowfs].fieldstop/dtheta_min*0.5)*2;
	if(ncomp>nembed){
	    error("nembed=%ld is smaller than ncomp=%ld\n", nembed, ncomp);
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
    pywfs->si=cellnew(4,1);//for each quadrant.
    //Make loc_t symmetric to ensure proper sampling onto detector. Center of subaperture
    powfs[ipowfs].saloc=mksqloc(order, order, dsa, dsa, 
				(-order*0.5+0.5)*dsa, (-order*0.5+0.5)*dsa);
    loc_t *loc_fft=mksqloc(ncomp, ncomp, dx2, dx2, (-ncomp2+0.5)*dx2, (-ncomp2+0.5)*dx2);
    for(int iy=0; iy<2; iy++){
	for(int ix=0; ix<2; ix++){
	    pywfs->si->p[ix+iy*2]=mkh(loc_fft, powfs[ipowfs].saloc, NULL, 
				      ((ix-0.5)*dx2*ncomp2), 
				      ((iy-0.5)*dx2*ncomp2),
				      1, 0, 0);	    
	}
    }
    if(parms->save.setup){
	writebin(pywfs->si, "pywfs_si0");
    }

    {
	//Determine subapertures area
	dmat *opd=dnew(pywfs->locfft->loc->nloc, 1);
	dmat *ints=0;
	pywfs_fft(&ints, powfs[ipowfs].pywfs, opd);
	if(parms->save.setup){
	    writebin(ints, "pywfs_ints0");
	}
	const int nsa=ints->nx;
	pywfs->saa=dnew(nsa, 1);
	for(int i=0; i<ints->nx; i++){
	    pywfs->saa->p[i]=ints->p[i]+ints->p[i+nsa]+ints->p[i+nsa*2]+ints->p[i+nsa*3];
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
		pywfs->si->p[ix+iy*2]=mkh(loc_fft, powfs[ipowfs].saloc, NULL, 
					  ((ix-0.5)*dx2*ncomp2), 
					  ((iy-0.5)*dx2*ncomp2),
					  1, 0, 0);	    
	    }
	}
    }
    dscale(pywfs->saa, pywfs->saa->nx/dsum(pywfs->saa));//saa average to one.
    locfree(loc_fft);
  
    if(parms->dbg.pywfs_atm){
	pywfs->atm=genatm_loc(pywfs->locfft->loc, parms->atm.r0, dsa);
	writebin(pywfs->atm, "pywfs_loc_atm");
    }

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
	const int nsa=powfs[ipowfs].saloc->nloc;
	dmat *TT=pywfs_tt(pywfs);
	double gxm=0, gym=0;
	for(int isa=0; isa<nsa; isa++){
	    gxm+=TT->p[isa];
	    gym+=TT->p[isa+nsa*3];
	}
	double gainscl=2.*nsa/(gxm+gym);
	//gainscl*=2.3; //2.3 is imperical adjustment.
	pywfs->gain*=gainscl;
	dscale(pywfs->gradoff, gainscl);
	info("pywfs_gain=%g\n", pywfs->gain);
	dfree(TT);
    }
    if(parms->save.setup){
	writebin(powfs[ipowfs].loc, "powfs%d_loc", ipowfs);
	writebin(powfs[ipowfs].saloc, "powfs%d_saloc", ipowfs);
	writebin(pywfs->amp, "pywfs_amp");
	writebin(pywfs->locfft->embed, "pywfs_embed");
	writebin(pywfs->pyramid, "pywfs_pyramid");
	writebin(nominal, "pywfs_nominal");
	writebin(pywfs->si, "pywfs_si");
	writebin(pywfs->gradoff, "pywfs_gradoff");
    }
    if(0){//Test implementation
	dmat *ints=0;
	int nn=1;
	double wve=1e-9*10;
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
    }
   
  
 
    /*
    {//test GA
	double dsa=parms->aper.d/order;
	loc_t *ploc=mkcirloc(parms->aper.d+dsa*2, dsa/2);
	dsp *gg=pywfs_mkg(powfs[ipowfs].pywfs, ploc);
	writebin(ploc, "pywfs_ploc");
	writebin(gg, "pywfs_gp");
	locfree(ploc);
	dspfree(gg);
	}*/
 
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
    double pos_r=pywfs->modulate*0.5;
    long ncomp=pywfs->nominal->nx;
    long ncomp2=ncomp/2;
    cmat *otf=cnew(ncomp, ncomp);
    dmat *pupraw=dnew(ncomp, ncomp);
    //writebin(psfs, "cpu_wvf");
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
	    long ny2=MIN(ncomp+offy2, nembed)-offy2;

	    long offx=(long)round(posx/dtheta);
	    long offx2=nembed/2+offx-ncomp2;
	    long ix0=MAX(-offx2, 0);
	    long nx2=MIN(ncomp+offx2, nembed)-offx2;

	    czero(otf);
	    dcomplex *pyramid=pywfs->pyramid->p[iwvl]->p;
	    for(long iy=iy0; iy<ny2; iy++){
		for(long ix=ix0; ix<nx2; ix++){
		    long indin=ix+offx2+(iy+offy2)*nembed;
		    long indout=ix+iy*ncomp;
		    otf->p[indout]=psfs->p[iwvl]->p[indin]*pyramid[indout];
		}
	    }
	    cfft2(otf, 1);
	    cabs22d(&pupraw, 1., otf, wvlwts->p[iwvl]/(ncomp*ncomp*pos_n));
	}//for iwvl
    }//for ipos
    //writebin(pupraw, "cpu_psf"); exit(0);
    ccpd(&otf, pupraw);//pupraw sum to one.
    dfree(pupraw);
    //cfftshift(otf);
    cfft2(otf, -1);
    ccwm(otf, pywfs->nominal);
    cfft2(otf, 1);
    //cfftshift(otf);
    const int nsa=pywfs->si->p[0]->nx;
    if(!(*ints)){
	(*ints)=dnew(nsa, 4);
    }
    for(int i=0; i<4; i++){
	//normalized so that each "subaperture" sum to 1.
	dspmulcreal((*ints)->p+nsa*i, pywfs->si->p[i], otf->p, (double)nsa/(ncomp*ncomp));
    }
    //writebin(*ints, "cpu_ints"); exit(0);
    /*{
	static int count=-1; count++;
	writebin(otf, "otf_%d", count);
	writebin(opd, "opd_%d", count);
	writebin(*ints, "ints_%d", count);
	exit(0);
	}*/
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
    double isum=dsum(ints);
    double alpha0=gain*nsa/isum;
    for(int isa=0; isa<nsa; isa++){
	double alpha2=alpha0/pywfs->saa->p[isa];
	pgx[isa]=(pi[1][isa]-pi[0][isa]
		  +pi[3][isa]-pi[2][isa])*alpha2;
	pgy[isa]=(pi[2][isa]+pi[3][isa]
		  -pi[0][isa]-pi[1][isa])*alpha2;
    }
    if(pywfs->gradoff){
	dadd(pgrad, 1, pywfs->gradoff, -1);
    }
}
/**
   Return T/T mode
*/
dmat *pywfs_tt(const PYWFS_T *pywfs){
    const loc_t *loc=pywfs->locfft->loc;
    dmat *opd=dnew(loc->nloc,1);
    dmat *ints=0;
    long nsa=pywfs->si->p[0]->nx;
    dmat *out=dnew(nsa*2,2);
    dmat *gradx=drefcols(out, 0, 1);
    dmat *grady=drefcols(out, 1, 1);
    dmat *gradx2=dnew(nsa*2,1);
    dmat *grady2=dnew(nsa*2,1);
    //dmat *grad0=dnew(nsa*2,1);
    double ptt[3]={0,0,0};
    double alpha=0.005/206265.;
    if(pywfs->atm){
	dcp(&opd, pywfs->atm);
    }
    /*dzero(ints);
    pywfs_fft(&ints, pywfs, opd);
    pywfs_grad(&grad0, pywfs, ints);*/
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

    /*writebin(gradx, "pywfs_gradx");
    writebin(grady, "pywfs_grady");
    writebin(gradx2, "pywfs_gradx2");
    writebin(grady2, "pywfs_grady2");
    writebin(grad0, "pywfs_grad0");
    exit(0);*/
    dadd(&gradx, 1, gradx2, -1);
    dadd(&grady, 1, grady2, -1);
    dscale(out, 0.5/alpha);
    dfree(gradx);
    dfree(grady);
    dfree(gradx2);
    dfree(grady2);
    dfree(opd);
    dfree(ints);
    //dfree(grad0);
    return out;
}
/**
   There is no need to simulate turbulence to fake optical gain. Any optical
   gain can be used as long as it "correct", i.e., a radian of tilt produces one
   radian of tilt, which is gauranteed when pywfs->gain is computed under the
   same conditions.
 */
static dsp *pywfs_mkg_do(const PYWFS_T *pywfs, const loc_t* ploc, int cubic, double iac){
    const loc_t *loc=pywfs->locfft->loc;
    dmat *opd=dnew(loc->nloc, 1);
    const double poke1=1./PYWFS_POKE;
    dmat *ints=0;
    const int nsa=pywfs->si->p[0]->nx;
    dmat *grad=dnew(nsa*2,1);
    dmat *grad0=dnew(nsa*2,1);
    dmat *opdin=dnew(ploc->nloc, 1);
    dsp *gg=dspnew(nsa*2, ploc->nloc, nsa*2*ploc->nloc);
    long count=0;
    if(pywfs->atm){
	dcp(&opd, pywfs->atm);
    }
    pywfs_fft(&ints, pywfs, opd);
    pywfs_grad(&grad0, pywfs, ints);
    for(int ia=0; ia<ploc->nloc; ia++){
	info2("%d of %ld\n", ia, ploc->nloc);
	if(pywfs->atm){
	    dcp(&opd, pywfs->atm);
	}else{
	    dzero(opd);
	}
	dzero(ints);
	dzero(grad);
	opdin->p[ia]=PYWFS_POKE;
	if(ia>0){
	    opdin->p[ia-1]=0;//remove old value from opd
	}
	if(cubic){
	    prop_nongrid_cubic((loc_t*)ploc, opdin->p, loc, NULL, opd->p, 1, 0, 0, 1, iac, 0, 0);
	}else{
	    prop_nongrid((loc_t*)ploc, opdin->p, loc, NULL, opd->p, 1, 0, 0, 1, 0, 0);
	}
	pywfs_fft(&ints, pywfs, opd);
	pywfs_grad(&grad, pywfs, ints);
	dadd(&grad, 1, grad0, -1);
	gg->p[ia]=count;
	const double thres=dmaxabs(grad)*EPS;
	for(int ig=0; ig<grad->nx; ig++){
	    if(fabs(grad->p[ig])>thres){
		gg->x[count]=grad->p[ig]*poke1;
		gg->i[count]=ig;
		count++;
	    }
	}
	/*
	    writebin(opd, "opd_%d", ia);
	    writebin(ints, "ints_%d", ia);
	    writebin(grad, "grad_%d", ia);
	    if(ia>50){
		exit(0);
	    }
	*/
	
    }
    gg->p[ploc->nloc]=count;
    dspsetnzmax(gg, count);
    cellfree(ints);
    dfree(opd);
    dfree(opdin);
    dfree(grad);
    dfree(grad0);
    return gg;
}
/**
   ploc is on pupil.
 */
dsp* pywfs_mkg(const PYWFS_T *pywfs, const loc_t* ploc, int cubic, double iac){
    uint32_t key=0;
    key=hashlittle(ploc->locx, ploc->nloc*sizeof(double), key);
    key=hashlittle(ploc->locy, ploc->nloc*sizeof(double), key);
    key=dhash(pywfs->wvlwts, key);
    key=dhash(pywfs->amp, key);
    key=chash(pywfs->nominal, key);
    key=chash(pywfs->pyramid->p[0], key);
    key=hashlittle(pywfs->si->p[0]->i, pywfs->si->p[0]->nzmax*sizeof(spint), key);
    key=hashlittle(pywfs->si->p[0]->x, pywfs->si->p[0]->nzmax*sizeof(double), key);
    if(pywfs->atm){
	key=dhash(pywfs->atm, key);
    }
    char fn[PATH_MAX];
    char fnlock[PATH_MAX];
    mymkdir("%s/.aos/pywfs/", HOME);
    snprintf(fn, PATH_MAX, "%s/.aos/pywfs/G_%u_%ld_%g.bin", HOME, 
	     key, pywfs->si->p[0]->nx, PYWFS_POKE);
    snprintf(fnlock, PATH_MAX, "%s.lock", fn);
    dsp *gg=0;
  retry:
    if(exist(fnlock) || !zfexist(fn)){
	int fd=lock_file(fnlock, 0, 0);//non blocking, exclusive
	if(fd>0){//succeed
	    info2("Generating PYWFS poke matrix\n");
	    gg=pywfs_mkg_do(pywfs, ploc, cubic, iac);
	    writebin(gg, "%s", fn);
	    snprintf(fn, PATH_MAX, "%s/.aos/pywfs/", HOME);
	    remove_file_older(fn, 30*24*3600);
	    close(fd); remove(fnlock);
	}else{
	    info2("Trying to lock %s\n", fnlock);
	    fd=lock_file(fnlock, 1, 0);
	    close(fd); remove(fnlock);
	    goto retry;
	}
    }else{
	gg=dspread(fn);
    }
    return gg;
}
/*
  Temporary wrapping
*/
dsp *pywfs_mkg_ga(const PARMS_T *parms, const POWFS_T *powfs, loc_t *aloc,
		  int iwfs, int idm){
    const int ipowfs=parms->wfs[iwfs].powfs;
    PYWFS_T *pywfs=powfs[ipowfs].pywfs;
    uint32_t key=0;
    key=hashlittle(aloc->locx, aloc->nloc*sizeof(double), key);
    key=hashlittle(aloc->locy, aloc->nloc*sizeof(double), key);
    key=dhash(pywfs->wvlwts, key);
    key=dhash(pywfs->amp, key);
    key=chash(pywfs->nominal, key);
    key=chash(pywfs->pyramid->p[0], key);
    key=hashlittle(pywfs->si->p[0]->i, pywfs->si->p[0]->nzmax*sizeof(spint), key);
    key=hashlittle(pywfs->si->p[0]->x, pywfs->si->p[0]->nzmax*sizeof(double), key);
    if(pywfs->atm){
	key=dhash(pywfs->atm, key);
    }
    char fn[PATH_MAX];
    char fnlock[PATH_MAX];
    mymkdir("%s/.aos/pywfs/", HOME);
    snprintf(fn, PATH_MAX, "%s/.aos/pywfs/G_%u_%ld_%g.bin", HOME, 
	     key, pywfs->si->p[0]->nx, PYWFS_POKE);
    snprintf(fnlock, PATH_MAX, "%s.lock", fn);
    dsp *gg=0;
  retry:
    if(exist(fnlock) || !zfexist(fn)){
	int fd=lock_file(fnlock, 0, 0);//non blocking, exclusive
	if(fd>0){//succeed
	    info2("Generating PYWFS poke matrix\n");
#if USE_CUDA
	    if(parms->gpu.wfs){
		gg=gpu_pywfs_mkg(parms, powfs, aloc, iwfs, idm);
	    }else
#endif
		gg=pywfs_mkg_do(pywfs, aloc, 
				parms->dm[idm].cubic, parms->dm[idm].iac);
	    writebin(gg, "%s", fn);
	    snprintf(fn, PATH_MAX, "%s/.aos/pywfs/", HOME);
	    remove_file_older(fn, 30*24*3600);
	    close(fd); remove(fnlock);
	}else{
	    info2("Trying to lock %s\n", fnlock);
	    fd=lock_file(fnlock, 1, 0);
	    close(fd); remove(fnlock);
	    goto retry;
	}
    }else{
	gg=dspread(fn);
    }
    return gg;
}
void pywfs_free(PYWFS_T *pywfs){
    if(!pywfs) return;
    dfree(pywfs->amp);
    locfft_free(pywfs->locfft);
    cellfree(pywfs->pyramid);
    cfree(pywfs->nominal);
    dspcellfree(pywfs->si);
    free(pywfs);
}

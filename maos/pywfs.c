/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include <unistd.h>
#include "common.h"
#include "setup_powfs.h"
#include "pywfs.h"
#include "../cuda/gpu.h"
#define PYWFS_GUARD 1.5 //separate the pupil by this amount in relative
#define PWFS_DEBUG 0 //For testing
/**
   \file pywfs.h

   Setup pyramid WFS and do simulation.
*/

static void pywfs_mksi(pywfs_t* pywfs, loc_t* loc_fft, loc_t* saloc0, real dx2, real pupelong){
	dspcellfree(pywfs->si);
	cellfree(pywfs->msaloc);
	const int pyside=pywfs->nside;
	pywfs->si=dspcellnew(pyside, 1);
	if(!pywfs->sioff){
		pywfs->sioff=dnew(pyside, 2);
	}
	const real dsa=saloc0->dx;
	const long ncomp2=NX(pywfs->nominal)/2;
	for(int ind=0; ind<pyside; ind++){
		const int iy=ind/2;//4-sided
		const int ix=ind%2;//4-sided
		loc_t* saloc=0;
		if(pupelong){//pupil elongation (along radial direction)
			if(pyside!=4){
				error("Revise implementation\n");
			} else{
				if(!pywfs->msaloc){
					pywfs->msaloc=loccellnew(pyside, 1);
				}
				saloc=locdup(saloc0);
				real angle=atan2(iy-0.5, ix-0.5);
				//squeeze the detector pixel coordinate radially to simulate pupil elongation
				real frac=1.-pupelong;
				saloc=P(pywfs->msaloc,ind)=locdup(saloc0);
				locstretch(saloc, angle, frac);
			}
		} else{
			saloc=saloc0;
		}
		real shx=0, shy=0;
		if(pywfs->pupilshift){
			shx=P(pywfs->pupilshift, ind, 0)*dsa;
			shy=P(pywfs->pupilshift, ind, 1)*dsa;
		}
		real offx=0, offy=0;
		switch(pyside){
		case 2:
			offx=ix-0.5;
			offy=0.5;
			break;
		case 3:
		{
			if(ind==0){
				offx=0;
				offy=-0.5;
			} else{
				offx=(ind-1.5)*sqrt(3.)*0.5;
				offy=0.25;
			}
		}
		break;
		case 4:
			offx=ix-0.5;
			offy=iy-0.5;
			break;
		default:
			error("Invalid dbg.pwfs_side=%d\n", pyside);
		}
		P(pywfs->sioff, ind, 0)=offx;
		P(pywfs->sioff, ind, 1)=offy;
		P(pywfs->si,ind)=mkh(loc_fft, saloc,
			(offx*ncomp2)*dx2+shx,
			(offy*ncomp2)*dx2+shy,
			1.);
	}
}
/**
   Setup pyramid WFS based on configuration.
*/
void pywfs_setup(powfs_t* powfs, const parms_t* parms, aper_t* aper, int ipowfs){
	pywfs_free(powfs[ipowfs].pywfs);
	pywfs_t* pywfs=powfs[ipowfs].pywfs=mycalloc(1, pywfs_t);
	const int pyside=pywfs->nside=parms->dbg.pwfs_side;
	map_t* map=0;
	pywfs->hs=parms->powfs[ipowfs].hs;
	pywfs->hc=parms->powfs[ipowfs].hc;
	pywfs->sigmatch=parms->powfs[ipowfs].sigmatch;
	pywfs->siglev=parms->powfs[ipowfs].siglev;
	pywfs->poke=parms->recon.poke;//How many meters to poke
	if(pywfs->poke>1e-5||pywfs->poke<1e-10){
		warning("poke=%g m is out of range\n", pywfs->poke);
	}
	pywfs->iwfs0=P(parms->powfs[ipowfs].wfs,0);
	real dx=parms->powfs[ipowfs].dx;
	create_metapupil(&map, 0, 0, parms->dirs, parms->aper.d, 0, dx, dx, 0, 0, 0, 0, 0, 0);
	powfs[ipowfs].loc=map2loc(map, 0);
	mapfree(map);
	powfs[ipowfs].amp=mkamp(powfs[ipowfs].loc, aper->ampground,
		P(parms->misreg.pupil,0), P(parms->misreg.pupil,1),
		parms->aper.d, parms->aper.din);
	loc_reduce(powfs[ipowfs].loc, powfs[ipowfs].amp, EPS, 0, NULL);
	//For convenience.
	pywfs->loc=locref(powfs[ipowfs].loc);
	pywfs->amp=dref(powfs[ipowfs].amp);
	setup_powfs_misreg_tel(powfs, parms, aper, ipowfs);
	setup_powfs_misreg_dm(powfs, parms, aper, ipowfs);
	if(powfs[ipowfs].realamp){
		warning("realamp is already set\n");
		dcellfree(powfs[ipowfs].realamp);
	}
	powfs[ipowfs].realamp=dcellnew(parms->powfs[ipowfs].nwfs, 1);
	for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
		int iwfs=P(parms->powfs[ipowfs].wfs,jwfs);
		if(parms->misreg.tel2wfs&&parms->misreg.tel2wfs[iwfs]){
			P(powfs[ipowfs].realamp,jwfs)=dref(P(powfs[ipowfs].amp_tel,jwfs));
		} else{
			P(powfs[ipowfs].realamp,jwfs)=dref(powfs[ipowfs].amp);
		}
	}

	int nwvl=NX(parms->powfs[ipowfs].wvl);
	real oversize=2*PYWFS_GUARD;
	pywfs->locfft=locfft_init(powfs[ipowfs].loc, pywfs->amp, parms->powfs[ipowfs].wvl, 0, oversize, 0);
	pywfs->wvlwts=ddup(parms->powfs[ipowfs].wvlwts);
	pywfs->modulate=parms->powfs[ipowfs].modulate;
	pywfs->modulpos=pywfs->modulate>0?(parms->powfs[ipowfs].modulpos/pyside*pyside):1;
	pywfs->modulring=pywfs->modulate>0?MAX(1, parms->powfs[ipowfs].modulring):1;
	long nembed=P(pywfs->locfft->nembed,0);
	real wvlmin, wvlmax;
	dmaxmin(P(parms->powfs[ipowfs].wvl), nwvl, &wvlmax, &wvlmin);
	real dtheta_min=wvlmin/(dx*nembed);
	pywfs->gain=1;
	//size of the part of the PSF captured by pyramid
	long ncomp=nembed;
	if(parms->powfs[ipowfs].fieldstop){
		ncomp=nextfftsize(ceil(parms->powfs[ipowfs].fieldstop/dtheta_min));
		if(ncomp>nembed){
			ncomp=nembed;
			warning("pywfs sampling is too coarse. Need %ld but has %ld pixels\n",
				ncomp, nembed);
		}
	}

	const long ncomp2=ncomp/2;
	pywfs->pyramid=ccellnew(nwvl, 1);
	dmat* pyramid=dnew(ncomp, ncomp);
	for(int iwvl=0; iwvl<nwvl; iwvl++){
		P(pywfs->pyramid,iwvl)=cnew(ncomp, ncomp);
		cmat* pp=P(pywfs->pyramid,iwvl)/*PCMAT*/;
		comp coeff=COMPLEX(0, M_PI*0.5);
		long skip=0;
		real dtheta=P(parms->powfs[ipowfs].wvl,iwvl)/(dx*nembed);//PSF sampling
		int nstop=ncomp;
		if(parms->powfs[ipowfs].fieldstop){//Limit fov per wvl
			nstop=ceil(parms->powfs[ipowfs].fieldstop/dtheta*0.5)*2;
			skip=(ncomp-nstop)/2;
			if(skip<0) skip=0;
		}
		real radius2=nstop*nstop*0.25;
		//Make pyramid edge or vertax flat within certain range
		real eskip=(parms->dbg.pwfs_flate/dtheta/2);
		real vskip=(parms->dbg.pwfs_flatv/dtheta/2);
		const real sqrt3=sqrt(3.);
		//const real ratio2=acos(sqrt(0.5))/acos(0.5);
		for(long iy=skip; iy<ncomp-skip; iy++){
			for(long ix=skip; ix<ncomp-skip; ix++){
				real xd=fabs(ix-ncomp2);
				real yy=iy-ncomp2;
				real yd=fabs(yy);
				real opd=0;
				if(!(pyside==4&&(xd<eskip||yd<eskip||(xd<vskip&&yd<vskip)))
					&&(xd*xd+yd*yd)<radius2){
					 //not on flat edge
					switch(pyside){
					case 2://roof with slope 
						opd=(xd+yy);
						break;
					case 3://3 sided pyrmaid
					/*
					  We keep one side of the pyramid the same as slope as a
					  roof in order to shift pupil by half.
					*/
						if(yy*sqrt3>xd){
							opd=yy;
						} else{
							//Rotate the coordinate by 120 and apply above surface formula.
							opd=fabs(yy-xd*sqrt3)*0.5;
						}
						break;
					case 4://4-sided pyramid
						opd=(xd+yd);
						break;
					default:
						error("Invalid pwfs_side=%d\n", pyside);
					}
				}
				P(pp, ix, iy)=cexp(opd*coeff);
				P(pyramid, ix, iy)=opd;//saving only.
			}
		}
	}
	if(parms->save.setup){
		writebin(pyramid, "powfs%d_pyramid", ipowfs);
	}
	dfree(pyramid);
	//Detector transfer function (sampling onto pixels).
	cmat* nominal=pywfs->nominal=cnew(ncomp, ncomp);
	cmat* pn=nominal/*PCMAT*/;
	long order=parms->powfs[ipowfs].order;
	real dsa=parms->powfs[ipowfs].dsa;//size of detector pixel mapped on pupil
	real dx2=dx*nembed/ncomp;//sampling of pupil after inverse fft
	real du=1./(dx2*ncomp);
	real dupix=dsa*du;
	real pdmeter=pow(dsa/dx2, 2);
	real pixblur=parms->powfs[ipowfs].pixblur*dsa;
	real e0b=-2*pow(M_PI*pixblur*du, 2);
	for(int iy=0; iy<ncomp; iy++){
		int jy=iy-ncomp2;
		for(int ix=0; ix<ncomp; ix++){
			int jx=ix-ncomp2;
			P(pn, ix, iy)=sinc(jy*dupix)*sinc(jx*dupix)*pdmeter;
			if(pixblur){
				P(pn, ix, iy)*=exp(e0b*(jx*jx+jy*jy));
			}
		}
	}
	cfftshift(nominal);
	/*
	cfft2(nominal, -1);
	cfftshift(nominal);
	cfft2(nominal, 1);
	cscale(nominal, 1./(NX(nominal)*NY(nominal)));
	*/
	if(parms->dbg.pwfs_psx){
		if(NX(parms->dbg.pwfs_psx)!=4){
			error("dbg.pwfs_psx has wrong format: expected 4x1, got %ldx%ld\n",
				NX(parms->dbg.pwfs_psx), NY(parms->dbg.pwfs_psx));
		}
		if(NX(parms->dbg.pwfs_psy)!=4){
			error("dbg.pwfs_psy has wrong format: expected 4x1, got %ldx%ld\n",
				NX(parms->dbg.pwfs_psy), NY(parms->dbg.pwfs_psy));
		}
		int redefine=parms->powfs[ipowfs].saloc?0:1;
		pywfs->pupilshift=dnew(4, 2);
		for(int i=0; i<4; i++){
			real tmp;
			tmp=P(parms->dbg.pwfs_psx, i);
			P(pywfs->pupilshift, i, 0)=tmp-redefine?round(tmp):0;
			tmp=P(parms->dbg.pwfs_psy, i);
			P(pywfs->pupilshift, i, 1)=tmp-redefine?round(tmp):0;

		}
		if(parms->save.setup){
			writebin(pywfs->pupilshift, "powfs%d_pupilshift", ipowfs);
		}
	}

	//Make loc_t symmetric to ensure proper sampling onto detector. Center of subaperture
	if(parms->powfs[ipowfs].saloc){
		powfs[ipowfs].saloc=locread("%s", parms->powfs[ipowfs].saloc);
		if(fabs(powfs[ipowfs].saloc->dx-dsa)>1e-6*fabs(dsa)){
			warning("loaded saloc has dx=%g, while powfs.dsa=%g\n",
				powfs[ipowfs].saloc->dx, dsa);
		}
	} else{
	//Pad the grid to avoid missing significant pixels (subapertures).
		long order2=ceil(order)+2*MAX(0, ceil(parms->dbg.pwfs_pupelong));
		if(order2>order){
			warning("order=%ld, order2=%ld.\n", order, order2);
		}
		powfs[ipowfs].saloc=mksqloc(order2, order2, dsa, dsa,
			(-order2*0.5+0.5)*dsa, (-order2*0.5+0.5)*dsa);
	}

	loc_t* loc_fft=mksqloc(ncomp, ncomp, dx2, dx2, (-ncomp2+0.5)*dx2, (-ncomp2+0.5)*dx2);
	const real pupelong=parms->dbg.pwfs_pupelong*sqrt(2)/(order*0.5);
	pywfs_mksi(pywfs, loc_fft, powfs[ipowfs].saloc, dx2, pupelong);//for each quadrant.

	if(parms->save.setup){
		writebin(pywfs->si, "powfs%d_si0", ipowfs);
		locwrite(pywfs->locfft->loc, "powfs%d_locfft", ipowfs);
		writebin(powfs[ipowfs].saloc, "powfs%d_saloc0", ipowfs);
	}
	{
	//Determine subapertures area
		dmat* opd=dnew(pywfs->locfft->loc->nloc, 1);
		dmat* ints=0;
		pywfs_fft(&ints, powfs[ipowfs].pywfs, opd);
		if(parms->save.setup){
			writebin(ints, "powfs%d_ints0", ipowfs);
		}
		int nints=NX(ints);
		pywfs->saa=dnew(nints, 1);
		for(int i=0; i<NX(ints); i++){
			for(int j=0; j<NY(ints); j++){
				P(pywfs->saa,i)+=P(ints, i, j);
			}
		}
		cellfree(ints);
		dfree(opd);
	}
	if(parms->powfs[ipowfs].saat>0&&!parms->powfs[ipowfs].saloc){
	//Get ride of insufficiently illuminated subapertures.
		dmat* saa=pywfs->saa;
		real samax=dmaxabs(saa);
		loc_reduce(powfs[ipowfs].saloc, saa, parms->powfs[ipowfs].saat*samax, 0, 0);
		dscale(saa, NX(saa)/dsum(saa));//saa average to one.
		pywfs_mksi(pywfs, loc_fft, powfs[ipowfs].saloc, dx2, pupelong);
	}
	powfs[ipowfs].saa=dref(pywfs->saa);
	const int nsa=powfs[ipowfs].saloc->nloc;
	dscale(pywfs->saa, NX(pywfs->saa)/dsum(pywfs->saa));//saa average to one.
	locfree(loc_fft);
	if(parms->save.setup){
		writebin(powfs[ipowfs].loc, "powfs%d_loc", ipowfs);
		writebin(powfs[ipowfs].saloc, "powfs%d_saloc", ipowfs);
		writebin(powfs[ipowfs].saa, "powfs%d_saa", ipowfs);
		writebin(powfs[ipowfs].amp, "powfs%d_amp", ipowfs);
		writebin(pywfs->locfft->embed, "powfs%d_embed", ipowfs);
		writebin(nominal, "powfs%d_nominal", ipowfs);
		writebin(pywfs->si, "powfs%d_si", ipowfs);
	}
	//Determine the gain and offset of PyWFS
	{
	//offset: grad of a flat wavefront
		dmat* opd=dnew(pywfs->locfft->loc->nloc, 1);
		dmat* ints=0;
		dmat* goff=0;
		pywfs_fft(&ints, pywfs, opd);
		if(parms->save.setup){
			writebin(ints, "powfs%d_ints1", ipowfs);
		}
		pywfs_grad(&goff, pywfs, ints);
		if(parms->save.setup){
			writebin(goff, "powfs%d_goff1", ipowfs);
		}
		dadd(&pywfs->gradoff, 1, goff, 1);
		if(0){//test TT response
			real ptt[3]={0,0.001/206265,0};
			loc_add_ptt(opd, ptt, pywfs->locfft->loc);
			dzero(ints); dzero(goff);
			pywfs_fft(&ints, pywfs, opd);
			writebin(ints, "powfs%d_ttx_ints", ipowfs);
			pywfs_grad(&goff, pywfs, ints);
			writebin(goff, "powfs%d_ttx_grad", ipowfs);
			exit(0);
		}
		dfree(goff);
		dfree(opd);
		dfree(ints);
		//Determine optical gain.
		dmat* TT=pywfs_tt(pywfs);
		real gxm=0, gym=0;
		for(int isa=0; isa<nsa; isa++){
			gxm+=P(TT, isa, 0);
			gym+=P(TT, isa+nsa, 1);
		}
		//info2("gxm=%g, gym=%g.\n", gxm/nsa, gym/nsa);
		real gainscl=2.*nsa/(gxm+gym);
		/*
		  pywfs->gain is inverse of optical gain, to insure 1rad of input
		  tip/tilt wavefront gives 1rad of gradient output.*/
		pywfs->gain*=gainscl;
		dscale(pywfs->gradoff, gainscl);
		dscale(TT, gainscl);
		pywfs->GTT=TT;
		//dbg("pywfs_gain=%g\n", pywfs->gain);
	}
	//Determine the NEA. It will be changed by powfs.gradscale as dithering converges    
	{
		powfs[ipowfs].sanea=dcellnew(1, 1);
		dmat* sanea=P(powfs[ipowfs].sanea,0)=dnew(nsa, 3);
		real rne=parms->powfs[ipowfs].rne;
		for(int isa=0; isa<nsa; isa++){
			real ogi=pywfs->gain*parms->powfs[ipowfs].gradscale;
			real sig=P(pywfs->saa,isa)*parms->powfs[ipowfs].siglev;//siglev of subaperture
			P(sanea, isa, 0)=P(sanea, isa, 1)=pow(ogi/sig, 2)*(sig+4*rne*rne);
		}
	}
	if(parms->save.setup){
		writebin(pywfs->gradoff, "powfs%d_gradoff", ipowfs);
		writebin(powfs[ipowfs].sanea, "powfs%d_sanea", ipowfs);
		writebin(pywfs->GTT, "powfs%d_GTT", ipowfs);
	}

	if(0){//Test implementation using zernikes
		dmat* ints=0;
		int nn=1;
		real wve=1e-9*160;
		dmat* opds=zernike(pywfs->locfft->loc, parms->aper.d, 0, 3, 0);
		zfarr* pupsave=zfarr_init(nn, NY(opds), "ints");
		zfarr* grads=zfarr_init(nn, NY(opds), "grads");
		dmat* opd=0;
		dmat* grad=0;
		for(int im=0; im<NY(opds); im++){
			for(int j=0; j<nn; j++){
				info2("im=%d, j=%d\n", im, j);
				opd=dsub(opds, 0, 0, im, 1);
				dscale(opd, pow(2, j)*wve);
				dzero(ints);
				pywfs_fft(&ints, powfs[ipowfs].pywfs, opd);
				pywfs_grad(&grad, powfs[ipowfs].pywfs, ints);
				zfarr_push(pupsave, j+im*nn, ints);
				zfarr_push(grads, j+im*nn, grad);
				dfree(opd);
			}
		}
		zfarr_close(pupsave);
		zfarr_close(grads);
		cellfree(ints);
		dfree(opds);
		exit(0);
	}
	if(0){//Test linearity of a zenike mode with noise
		real wve=1e-9*20;
		dmat* opds=zernike(pywfs->locfft->loc, parms->aper.d, 0, 0, -parms->powfs[ipowfs].dither);
		dmat* opdi=0;
		dmat* ints=0, * grad=0;
		dadd(&opdi, 0, opds, wve);
		pywfs_fft(&ints, powfs[ipowfs].pywfs, opdi);
		pywfs_grad(&grad, powfs[ipowfs].pywfs, ints);
		dmat* reg=dpinv(grad, 0);
		writebin(grad, "dither_grad");
		writebin(reg, "dither_reg");
		writebin(ints, "dither_ints");
		rand_t rstat;
		seed_rand(&rstat, 1);
		dmat* tmp=0;
		int nj=10, nn=10;
		dmat* res=dnew(nj, nn);
		dmat* ints2=0;
		for(int j=0; j<nj; j++){
			dzero(ints);
			dadd(&opdi, 0, opds, wve*(j+1));
			pywfs_fft(&ints, powfs[ipowfs].pywfs, opdi);
			for(int in=0; in<nn; in++){
				dadd(&ints2, 0, ints, 100);
				addnoise(ints2, &rstat, 0, 0, 0, 0, 0, in, 1);
				pywfs_grad(&grad, powfs[ipowfs].pywfs, ints2);
				dmm(&tmp, 0, reg, grad, "nn", 1);
				P(res, j, in)=P(tmp,0);
				info2("%d of %d, %d of %d: %g\n", j, nj, in, nn, P(tmp,0));
			}
		}
		writebin(opds, "dither_opd");
		writebin(res, "dither_res");
		dfree(opds); dfree(opdi); dfree(ints); dfree(grad); dfree(reg); dfree(tmp); dfree(res); dfree(ints2);
		exit(0);
	}
	//Test NCPA calibration
	int PYWFS_NCPA=0;
	READ_ENV_INT(PYWFS_NCPA, 0, 1);
	if(PYWFS_NCPA){
		dmat* opdatm=dread("opdatm");
		dmat* opdbias_full=dread("opdbias_full");
		dmat* opdbias_astigx=dread("opdbias_astigx");
		dmat* opdbias_polish=dread("opdbias_polish");
		const real atmscale=1;
		OMP_TASK_FOR(4)
		for(int i=0; i<100; i++){
			info2("%d ", i);
			dmat* ints=0;
			dmat* grad=0;
			dmat* opd=0;

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
#if PWFS_DEBUG
	exit(0);
#endif
}
/**
   Perform FFT over the complex PSF with additional phases caused by the
   pyramid. FFT on each quadrant of the PSF creates diffraction effects.
*/
void pywfs_fft(dmat** ints, const pywfs_t* pywfs, const dmat* opd){
	locfft_t* locfft=pywfs->locfft;
	ccell* psfs=0;
	locfft_psf(&psfs, locfft, opd, NULL, 1);//psfs.^2 sum to 1. peak in center
	int nwvl=NX(locfft->wvl);
	real dx=locfft->loc->dx;
	long nembed=P(locfft->nembed,0);
	long nembed2=nembed/2;
	dmat* wvlwts=pywfs->wvlwts;
	//position of pyramid for modulation
	int pos_n=pywfs->modulpos;
	int pos_nr=pywfs->modulring;
	real pos_r=pywfs->modulate;
	long ncomp=NX(pywfs->nominal);
	long ncomp2=ncomp/2;
	cmat* otf=cnew(ncomp, ncomp);
	dmat* pupraw=dnew(ncomp, ncomp);
#if PWFS_DEBUG
	static int savec=-1; savec++;
	writebin(psfs, "pwfs_fft_cpu_psf_%d", savec);
#endif
	for(int ir=0; ir<pos_nr; ir++){
	//Radius of the current ring
		real pos_ri=pos_r*(ir+1)/pos_nr;
		//Scale number of points by ring size to have even surface brightness
		int pos_ni=pos_n*(ir+1)/pos_nr;
		for(int ipos=0; ipos<pos_ni; ipos++){
				//whether the first point falls on the edge or not makes little difference
			real theta=2*M_PI*((real)ipos/pos_ni);
			real posx=cos(theta)*pos_ri;
			real posy=sin(theta)*pos_ri;
			for(int iwvl=0; iwvl<nwvl; iwvl++){
				real dtheta=P(locfft->wvl,iwvl)/(dx*nembed);
				long offy=(long)round(posy/dtheta);
				long offy2=nembed2+offy-ncomp2;
				long iy0=MAX(-offy2, 0);
				long ny2=MIN(ncomp, nembed-offy2)-iy0;

				long offx=(long)round(posx/dtheta);
				long offx2=nembed/2+offx-ncomp2;
				long ix0=MAX(-offx2, 0);
				long nx2=MIN(ncomp, nembed-offx2)-ix0;

				czero(otf);
				comp* pyramid=P(P(pywfs->pyramid,iwvl));
				for(long iy=iy0; iy<ny2; iy++){
					for(long ix=ix0; ix<nx2; ix++){
						long indin=ix+offx2+(iy+offy2)*nembed;
						long indout=ix+iy*ncomp;
						P(otf,indout)=P(P(psfs,iwvl),indin)*pyramid[indout];
					}
				}
				cfft2(otf, 1);
				cabs22d(&pupraw, 1., otf, P(wvlwts,iwvl)/(ncomp*ncomp*pos_ni*pos_nr));
			}//for iwvl
		}//for ipos
	}//for ir
#if PWFS_DEBUG
	writebin(pupraw, "pwfs_fft_cpu_pupil_%d", savec);
#endif
	//writebin(pupraw, "cpu_psf"); exit(0);
	ccpd(&otf, pupraw);//pupraw sum to one.
	//writebin(otf, "cpu_wvf4");
	dfree(pupraw);
	cfft2(otf, -1);
	//writebin(otf, "cpu_wvf5");
	ccwm(otf, pywfs->nominal);
	cfft2(otf, 1);
	const int nsa=P(pywfs->si,0)->nx;
	if(!(*ints)){
		(*ints)=dnew(nsa, NX(pywfs->si));
	}
	for(int i=0; i<NX(pywfs->si); i++){
		//normalized so that each "subaperture" sum to 1.
		dspmulcreal(PCOL(*ints,i), P(pywfs->si,i), P(otf), (real)nsa/(ncomp*ncomp));
	}
	//writebin(*ints, "cpu_ints"); exit(0);
	ccellfree(psfs);
	cfree(otf);
}
/**
   Compute gradients. It replaces the result, not accumulate.
 */
void pywfs_grad(dmat** pgrad, const pywfs_t* pywfs, const dmat* ints){
	const long nsa=NX(ints);
	const int pyside=pywfs->nside;
	if(!*pgrad){
		*pgrad=dnew(nsa*2, 1);
	}
	real* pgx=P(*pgrad);
	real* pgy=P(*pgrad)+nsa;
	real gain=pywfs->gain;
	real triscalex=sqrt(3.)/2;
	real imean=0;
	if(pywfs->sigmatch==2){
		imean=dsum(ints)/nsa;
	}
	for(int isa=0; isa<nsa; isa++){
		real isum=0;
		switch(pywfs->sigmatch){
		case 0:
			info_once("PWFS: No siglev correction.\n");
			isum=pywfs->siglev*P(pywfs->saa,isa);
			break;
		case 1:
			info_once("PWFS: Individual siglev correction.\n");
			for(int i=0; i<pyside; i++){
				isum+=P(ints, isa, i);
			}
			break;
		case 2:
			info_once("PWFS: Global siglev correction.\n");//preferred.
			isum=imean*P(pywfs->saa,isa);
			break;
		}
		real alpha2=gain/isum;
		switch(pyside){
		case 3:
			pgx[isa]=(P(ints, isa, 1)-P(ints, isa, 2))*alpha2*triscalex;
			pgy[isa]=(P(ints, isa, 0)-0.5*(P(ints, isa, 1)+P(ints, isa, 2)))*alpha2;
			break;
		case 4:
			pgx[isa]=(P(ints, isa, 0)-P(ints, isa, 1)
				+P(ints, isa, 2)-P(ints, isa, 3))*alpha2;
			pgy[isa]=(P(ints, isa, 0)+P(ints, isa, 1)
				-P(ints, isa, 2)-P(ints, isa, 3))*alpha2;
			break;
		}
	}

	if(pywfs->gradoff){
		dadd(pgrad, 1, pywfs->gradoff, -1);
	}
}
/**
   Return measurement of T/T mode, normalized for 1 unit of input.
*/
dmat* pywfs_tt(const pywfs_t* pywfs){
	TIC;tic;info2("Computing pywfs_tt...");
	const loc_t* loc=pywfs->locfft->loc;
	dmat* opd=dnew(loc->nloc, 1);
	dmat* ints=0;
	long nsa=P(pywfs->si,0)->nx;
	dmat* out=dnew(nsa*2, 2);
	dmat* gradx=drefcols(out, 0, 1);
	dmat* grady=drefcols(out, 1, 1);

	real ptt[3]={0,0,0};
	real alpha=0.005/206265.;

	//+x
	ptt[1]=alpha;  ptt[2]=0;
	loc_add_ptt(opd, ptt, loc);
	dzero(ints);
	pywfs_fft(&ints, pywfs, opd);
	pywfs_grad(&gradx, pywfs, ints);
#if PWFS_DEBUG
	writebin(ints, "pwfs_ttx");
#endif
	//+y
	ptt[1]=-alpha; ptt[2]=alpha;
	loc_add_ptt(opd, ptt, loc);
	dzero(ints);
	pywfs_fft(&ints, pywfs, opd);
	pywfs_grad(&grady, pywfs, ints);
#if PWFS_DEBUG
	writebin(ints, "pwfs_tty");
#endif
#if PWFS_DEBUG
#define pywfs_tT_DUAL 1
#else
#define pywfs_tT_DUAL 0
#endif
#if pywfs_tT_DUAL
	dmat* gradx2=dnew(nsa*2, 1);
	dmat* grady2=dnew(nsa*2, 1);
	//-x
	ptt[1]=-alpha; ptt[2]=-alpha;
	loc_add_ptt(opd, ptt, loc);
	dzero(ints);
	pywfs_fft(&ints, pywfs, opd);
	pywfs_grad(&gradx2, pywfs, ints);
#if PWFS_DEBUG
	writebin(ints, "pwfs_ttx2");
#endif
	//-y
	ptt[1]=+alpha; ptt[2]=-alpha;
	loc_add_ptt(opd, ptt, loc);
	dzero(ints);
	pywfs_fft(&ints, pywfs, opd);
	pywfs_grad(&grady2, pywfs, ints);
#if PWFS_DEBUG
	writebin(ints, "pwfs_tty2");
#endif
	dadd(&gradx, 1, gradx2, -1);
	dadd(&grady, 1, grady2, -1);
	dscale(out, 0.5/alpha);
	dfree(gradx2);
	dfree(grady2);
#else
	dscale(out, 1./alpha);
#endif
	dfree(gradx);
	dfree(grady);
	dfree(opd);
	dfree(ints);
	toc2("done");
	return out;
}
static uint32_t pywfs_hash(const pywfs_t* pywfs, uint32_t key){
	key=lochash(pywfs->loc, key);
	key=dhash(pywfs->amp, key);
	key=dhash(pywfs->saa, key);
	key=dhash(pywfs->wvlwts, key);
	key=chash(P(pywfs->pyramid,0), key);
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
static dmat* pywfs_mkg_do(const pywfs_t* pywfs, const loc_t* locin, const loc_t* locfft, const dmat* mod,
	real displacex, real displacey){
#if USE_CUDA
	if(pywfs->gpu){
		return gpu_pywfs_mkg(pywfs, locin, locfft, mod, displacex, displacey);
	}
#endif
	const int nsa=P(pywfs->si,0)->nx;
	dmat* grad0=dnew(nsa*2, 1);
	dmat* opd0;
	if(pywfs->opdadd){
		opd0=dref(pywfs->opdadd);
	} else{
		opd0=dnew(locfft->nloc, 1);
	}
	{
		dmat* ints=0;
		pywfs_fft(&ints, pywfs, opd0);
		pywfs_grad(&grad0, pywfs, ints);
		//writebin(grad0, "grad0_cpu");
		//writebin(ints, "ints0_cpu");
		dfree(ints);
	}
	int count=0;
	int nmod=mod?NY(mod):locin->nloc;
	dmat* ggd=dnew(nsa*2, nmod);
	if(mod&&NX(mod)!=locin->nloc){
		error("NX(mod) must equal to %ld", locin->nloc);
	}

	const real scale=1.-(locin->ht-pywfs->hc)/pywfs->hs;
	TIC;tic;
	OMP_TASK_FOR(4)
	for(int imod=0; imod<nmod; imod++){
		dmat* opdin=dnew(locin->nloc, 1);
		dmat* opdfft=ddup(opd0);
		dmat* ints=0;
		dmat* grad=drefcols(ggd, imod, 1);
		real poke=pywfs->poke;
		if(mod){
			dmat* tmp=drefcols(mod, imod, 1);
			//the equivalent radimodl order of zernike.
			//real radial=ceil((sqrt(8.*(imod+1)+1)-3)*0.5)+1;
			//real std=dstd(tmp);
			real tmax, tmin;
			dmaxmin(P(tmp), NX(tmp), &tmax, &tmin);
			poke/=(tmax-tmin);//sqrt(radial);
			dadd(&opdin, 0, tmp, poke);
			dfree(tmp);
		} else{
			P(opdin,imod)=poke;
		}
		prop_nongrid((loc_t*)locin, P(opdin), locfft, P(opdfft), 1, displacex, displacey, scale, 0, 0);
		//writebin(opdfft, "phiout_cpu_%d", imod);
		pywfs_fft(&ints, pywfs, opdfft);
		//writebin(ints, "ints_cpu_%d", imod);
		pywfs_grad(&grad, pywfs, ints);
		dadd(&grad, 1, grad0, -1);
		dscale(grad, 1./poke);
		atomicadd(&count, 1);
		if(count%10==0){
			real ts=myclockd()-tk;
			info2("%d of %ld. %.2f of %.2f seconds. std(grad)=%g.\n", count, locin->nloc, ts, ts/count*locin->nloc, dstd(grad));
		}
		dfree(opdfft);
		dfree(opdin);
		dfree(grad);
		dfree(ints);
	}
	info2("\n");
	dfree(grad0);
	dfree(opd0);
	return ggd;
}
/**
   locin is on pupil.
 */
dmat* pywfs_mkg(pywfs_t* pywfs, const loc_t* locin, const char* distortion, const dmat* mod, const dmat* opdadd,
	real displacex, real displacey){
	if(opdadd){
		dfree(pywfs->opdadd);
		pywfs->opdadd=dnew(pywfs->locfft->loc->nloc, 1);
		prop_nongrid(pywfs->loc, P(opdadd), pywfs->locfft->loc, P(pywfs->opdadd), 1, 0, 0, 1, 0, 0);
	}
	loc_t* locfft=pywfs->locfft->loc;
	if(distortion){
		locfft=loctransform(locfft, distortion);
	}
	/*if(mod && NY(mod)<=6){
	return pywfs_mkg_do(pywfs, locin, locfft, mod, displacex, displacey);
	}*/
	uint32_t key=0;
	key=lochash(locin, key);
	key=pywfs_hash(pywfs, key);
	if(mod) key=dhash(mod, key);
	if(opdadd) key=dhash(opdadd, key);
	char fn[PATH_MAX-10];
	char fnlock[PATH_MAX];
	mymkdir("%s/G/", CACHE);
	snprintf(fn, sizeof(fn), "%s/G/G_%u_%ld_%ld_%g_%d_%g_%g_%g_%g_%g_v2.bin", CACHE,
		key, P(pywfs->locfft->nembed,0), locin->nloc, pywfs->modulate, pywfs->modulpos,
		locin->iac, displacex, displacey, 1., pywfs->poke);
	snprintf(fnlock, sizeof(fnlock), "%s.lock", fn);
	info2("Using G in %s\n", fn);
	dmat* gg=0;
	/*
	if(0){//test amount of poke
	dmat *mod1=dnew(locin->nloc, 3);
	if(mod){
		memcpy(PCOL(mod1, 0), PCOL(mod, 0), sizeof(real)*NX(mod));
		memcpy(PCOL(mod1, 1), PCOL(mod, 1000), sizeof(real)*NX(mod));
		memcpy(PCOL(mod1, 2), PCOL(mod, 3200), sizeof(real)*NX(mod));
	}else{
		P(mod1, 30, 0)=1;//sqrt(locin->nloc);
		P(mod1, 800, 1)=1;//sqrt(locin->nloc);
		P(mod1, 3000, 2)=1;//sqrt(locin->nloc);
	}
	writebin(mod1, "mod1");
	dcell *gg1=dcellnew(13,1);
	real poke=1e-9;
	real step=pow(10,0.25);
	for(int ig=0; ig<NX(gg1); ig++){
		((pywfs_t*)pywfs)->poke=poke;
		P(gg1,ig)=pywfs_mkg_do(pywfs, locin, locfft, mod1, displacex, displacey);
		poke=poke*step;
	}
	writebin(gg1, "gg1");
	cellfree(gg1);
	exit(0);
	}*/
retry:
	if(exist(fnlock)||!zfexist("%s",fn)){
		int fd=lock_file(fnlock, 0, 0);//non blocking, exclusive
		if(fd>0){//succeed
			info2("Generating PYWFS poke matrix\n");
			gg=pywfs_mkg_do(pywfs, locin, locfft, mod, displacex, displacey);
			writebin(gg, "%s", fn);
			close(fd); remove(fnlock);
		} else{
			info2("Trying to lock %s\n", fnlock);
			fd=lock_file(fnlock, 1, 0);
			close(fd); remove(fnlock);
			goto retry;
		}
	} else{
		gg=dread("%s", fn);
	}
	if(distortion){
		locfree(locfft);
	}
	return gg;
}
/**
   frees pywfs_t
*/
void pywfs_free(pywfs_t* pywfs){
	if(!pywfs) return;
	cellfree(pywfs->amp);
	locfree(pywfs->loc);

	locfft_free(pywfs->locfft);
	cellfree(pywfs->pyramid);
	cellfree(pywfs->nominal);
	cellfree(pywfs->si);
	cellfree(pywfs->sioff);
	cellfree(pywfs->opdadd);
	cellfree(pywfs->gradoff);
	cellfree(pywfs->GTT);
	cellfree(pywfs->pupilshift);
	cellfree(pywfs->msaloc);
	cellfree(pywfs->wvlwts);
	cellfree(pywfs->saa);
	free(pywfs);
}

/*
  Copyright 2009-2026 Lianqi Wang
  
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
#include "phase.h"
/**
   Embed amp*exp(i*phi) into the center of wvf. Uses the same centering as
   dembed() so that it is consistent for odd sizes.

   @param[out] wvf The complex field. Zeroed and then filled.
   @param[in,out] tmp Temporary array of the same size as amp. Allocated if NULL.
   @param[in] phi The phase (radian). Can be NULL.
   @param[in] amp The amplitude.
*/
static void gs_embed(cmat *wvf, cmat **tmp, const dmat *phi, const dmat *amp){
	if(phi){
		if(!*tmp) *tmp=cnew(NX(amp), NY(amp));
		cembed_wvf(*tmp, P(phi), P(amp), NX(amp), NY(amp), TWOPI, 0);
		cembed(wvf, *tmp, 0);
	}else{
		cembedd(wvf, amp, 0);
	}
}
/**
   Replace the amplitude of a complex number by amp while keeping its phase.
   The phase is set to 0 if the input is 0.
*/
static inline comp gs_setamp(comp val, real amp){
	real a=cabs(val);
	return a>0?(val*(amp/a)):COMPLEX(amp, 0);
}
/**
   Estimate the pupil phase from the pupil and focal plane amplitude using the
   Gerchberg-Saxton algorithm.

   @param[out] pphi1 The final pupil phase estimate in radian. Same dimension as amp1. Can be NULL.
   @param[out] mphi1 The phase update in modes. Only valid if mode and rmod are set. Can be NULL.
   @param[in] amp1 Pupil amplitude map. Size is n*n.
   @param[in] amp2 Field amplitude map. Size is m*m, usually square root of the PSF, with peak in center.
   @param[in] phi1b The initial pupil phase estimate in radian. Can be NULL.
   @param[in] mode Modal matrix. If not null, project pupil phase update onto the mode span.
   @param[in] rmod Modal reconstructor. Used with mode.
   @param[in] fembed Embedding factor for amp1.
   @param[in] nrep Number of repetitions.
*/
void gerchberg_saxton(dmat **pphi1, dmat **mphi1, const dmat *amp1, const dmat *amp2, const dmat *phi1b,
	const dmat *mode, const dmat *rmod, const real fembed, const int nrep){
	if(!amp1 || !amp2){
		error("amp1, and am2 must be all set\n");
		return;
	}
	if(mode&&!rmod){
		error("rmod must be set when mode is set\n");
		return;
	}

	long nx=(long)round(NX(amp1)*fembed);
	long ny=(long)round(NY(amp1)*fembed);
	cmat *wvf1=cnew(nx, ny);
	cmat *wvf0=NULL;//temporary of the same size as amp1
	dmat *amp1e=NULL;//same size as wvf1
	dmat *amp2e=NULL;//same size as wvf1
	dmat *phi1e=NULL;//phase estimate.
	dmat *phi1v=NULL;//phi1e as vector
	dmat *phi1u=NULL;//phase estimated unembeded.
	dmat *mv=NULL;//modal estimate
	if(nx!=NX(amp1)||ny!=NY(amp1)){
		amp1e=dnew(nx, ny);
		dembed(amp1e, amp1, 0);
	} else{
		amp1e=dref(amp1);
	}

	if(nx!=NX(amp2) || ny!=NY(amp2)){
		amp2e=dnew(nx,ny);
		dembed(amp2e,amp2,0);
	}else{
		amp2e=ddup(amp2);//duplicate so that we can do fftshift
	}
	gs_embed(wvf1, &wvf0, phi1b, amp1);
	if(pphi1){
		dinit(pphi1, NX(amp1), NY(amp1));
		phi1u=*pphi1;
	}else{
		phi1u=dnew(NX(amp1),NY(amp1));
	}
	if(mphi1){
		dinit(mphi1, NX(rmod), 1);
		mv=*mphi1;
	}
	//info("nrep=%d, sum(amp1e)=%g, sum(amp2e)=%g, max(amp2e)=%g\n", nrep, dsum(amp1e), dsum(amp2e), dmax(amp2e));
	/*{
		static int count=0; count++;
		writebin(amp2e, "amp2e_%d", count);
	}*/
	dfftshift(amp2e);
	if(P(amp2e,0)>1 ||P(amp2e,0)<0.1){
		warning("Amplitude at the focal plane is %g, outside of range [0.1,1]\n", P(amp2e,0));
	}
	for(int irep=0; irep<nrep; irep++){
		cfft2(wvf1, -1);//to image plane
		for(long i=0; i<PN(wvf1); i++){
			P(wvf1,i)=gs_setamp(P(wvf1,i), P(amp2e,i));//fix image plane amplitude
		}
		cfft2(wvf1, 1);//to pupil plane
		if(mode || (irep+1)==nrep){
			carg2d(&phi1e, 0, wvf1, 1);
			dembed(phi1u, phi1e, 0);
		}
		if(mode){
			if(phi1b) dadd(&phi1u, 1, phi1b, -1);//remove initial phase
			if(!phi1v) phi1v=dref_reshape(phi1u, PN(phi1u), 1);
			dmm(&mv, 0, rmod, phi1v, "nn", 1);
			dmm(&phi1v, 0, mode, mv, "nn", 1);
			if(phi1b) dadd(&phi1u, 1, phi1b, 1);//add back initial phase
			gs_embed(wvf1, &wvf0, phi1u, amp1);
		}else{
			for(long i=0; i<PN(wvf1); i++){
				P(wvf1, i)=gs_setamp(P(wvf1, i), P(amp1e, i));//fix pupil plane amplitude
			}
		}
	}
	dfree(amp1e);
	dfree(amp2e);
	dfree(phi1e);
	dfree(phi1v);
	cfree(wvf1);
	cfree(wvf0);
	if(!mphi1) dfree(mv);
	if(!pphi1) dfree(phi1u);
}

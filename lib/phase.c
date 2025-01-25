/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
/// @brief 
/// @param pphi1 
/// @param mphi1 
/// @param phi1b 
/// @param amp1 
/// @param amp2 
/// @param mode 
/// @param rmod 
/// @param fembed 
/// @param nrep 
void gerchberg_saxton(
	dmat **pphi1, /**<[Out]The final pupil phase estimate in radian. same dimension as amp1.*/
	dmat **mphi1, /**<[Out]The phase update in modes only if mode and rmode are set.*/
  	const dmat *amp1, /**<[In]Pupil amplitude map. size is n*n. embeded. */
  	const dmat *amp2, /**<[In]Field amplitude map. size is m*m. usually square root of the PSF. */
	const dmat *phi1b,/**<[In] The intial final pupil phase estimate in radian.*/
  	const dmat *mode, /**[In]modal matrix. If not null, project pupil phase update onto the mode span*/
  	const dmat *rmod, /**[In]modeal reconstructor. If not null, project pupil phase update onto the mode span*/
  	const real fembed,/**<[in]Embeded factor for amp1*/
  	const int nrep    /**<[In]Number of repeatitions*/
  ){
	if(!amp1 || !amp2){
		error("amp1, and am2 must be all set\n");
	}

	long nx=(long)round(NX(amp1)*fembed);
	long ny=(long)round(NY(amp1)*fembed);
	cmat *wvf1=cnew(nx, ny);
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
	if(phi1b){
		cembed_wvf(wvf1, P(phi1b), P(amp1), NX(amp1), NY(amp1), TWOPI, 0);
	}else{
		cembedd(wvf1, amp1, 0);
	}
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
			P(wvf1,i)*=P(amp2e,i)/cabs(P(wvf1,i));//fix image plane amplitude
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
			cembed_wvf(wvf1, P(phi1u), P(amp1), NX(amp1), NY(amp1), TWOPI, 0);
		}else{
			for(long i=0; i<PN(wvf1); i++){
				P(wvf1, i)*=P(amp1e, i)/cabs(P(wvf1, i));//fix pupil plane amplitude
			}
		}
	}
	dfree(amp1e);
	dfree(amp2e);
	dfree(phi1e);
	dfree(phi1v);
	cfree(wvf1);
	if(!mphi1) dfree(mv);
	if(!pphi1) dfree(phi1u);
}

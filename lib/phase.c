/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>

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
dmat* gerchberg_saxton(
  const dmat *phi1b,      /**<[In/Out] The (intial and) final pupil phase estimate in radian.*/
  const dmat *amp1, /**<[In]Pupil amplitude map. size is n*n. embeded. */
  const dmat *amp2, /**<[In]Field amplitude map. size is n*n. usually square root of the PSF.*/
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
	dmat *phi1u=NULL;//unembeded phi1e
	dmat *phi1e=NULL;//phase estimate.
	dmat *phi1v=NULL;//phi1e as vector
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
	if(mode){
		phi1u=dnew(NX(amp1), NY(amp1));
	}
	dfftshift(amp2e);
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
	dfree(mv);
	return phi1u;
}

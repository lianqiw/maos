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
#include "hyst.h"

/**
    contains data related to DM hysterisis modeling for all the common DMs (not
	MOAO). let input be x, and output of each mode be y.  
*/
struct hyst_t{
	real alpha;       /**<Material constant*/
	real beta;        /**<constants*/
	real beta_m_u;    /**<constants, beta-u*/
	dmat* xlast;      /**<Record x from last time step*/
	dmat* ylast;      /**<Record y from last time step*/

};

/**
  Hysteresis modeling. Create the hysteresis model. It assumes that the hysteresis is measured with a stroke.
  For details, see the documentation 50_hyst.md
*/
hyst_t* hyst_new(
	real hysteresis,   /**<The ratio of the maximum difference of ascending/descending curve over 2*stroke. */
	real alpha,        /**<The alpha parameter*/
	real stroke,       /**<The DM command is from -stroke to stroke*/
	long nact          /**<Number of actuators*/
){
	hyst_t* hyst=mycalloc(1, hyst_t);
	if(!hyst) return NULL;
	hyst->xlast=dnew(nact, 1);
	hyst->ylast=dnew(nact, 1);
	if(hysteresis<0||hysteresis>1){
		error("Hysteresis should be between 0 and 1\n");
	}
	if(isinf(stroke)){
		stroke=10e-6;
	}
	if(alpha<0) alpha=2; //a good default
	real eas=exp(alpha);
	real enas=1./eas;
	real u=hysteresis*alpha/(1.-2./(eas+enas));
	hyst->beta=u/alpha*(1.-2.*enas/(eas+enas))+1.;
	hyst->beta_m_u=hyst->beta-u;
	hyst->alpha=alpha/stroke;
	dbg("hyst_new: stroke=%g, alpha=%g, beta=%g, u=%g, nact=%ld\n", stroke, alpha, hyst->beta, u, nact);
	return hyst;
}

/**
   Reset hysteresis state
*/
void hyst_reset(hyst_t* hyst){
	if(!hyst) return;
	dzero(hyst->xlast);
	dzero(hyst->ylast);
}

/**
   Free hysteresis model.
*/
void hyst_free(hyst_t* hyst){
	if(!hyst) return;
	dfree(hyst->xlast);
	dfree(hyst->ylast);
	free(hyst);
}

/**
   Apply hysteresis to DM vector. dmreal and dmcmd may be the same
*/
void hyst_dmat(hyst_t* hyst, dmat* dmreal, const dmat* dmcmd){
	if(!hyst){
		dcp(&dmreal, dmcmd);
		return;
	}
	real* restrict xlast=P(hyst->xlast);
	real* restrict ylast=P(hyst->ylast);
	assert(NX(dmcmd)==NX(hyst->xlast));
	for(int it=0; it<NY(dmcmd); it++){
		for(int ia=0; ia<NX(dmcmd); ia++){
			real xia=P(dmcmd, ia, it);
			real dx=xia-xlast[ia];
			real dy=fabs(dx)*hyst->alpha*(hyst->beta*xlast[ia]-ylast[ia])+dx*(hyst->beta_m_u);
			xlast[ia]=xia;
			ylast[ia]+=dy;
			P(dmreal, ia, it)=ylast[ia];
		}
	}
}

/**
   Apply hysteresis to set of DM vectors
*/
void hyst_dcell(hyst_t** hyst, dcell* dmreal, const dcell* dmcmd){
	if(!hyst){
		dcellcp(&dmreal, dmcmd);
		return;
	}
	for(int idm=0; idm<NX(dmcmd)*NY(dmcmd); idm++){
		hyst_dmat(hyst[idm], P(dmreal,idm), P(dmcmd,idm));
	}
}

/**
 * Test hysteresis model. Dmcmd has dimensions of nact*ntime
 * */
dmat* hyst_test(real hysteresis, real hyst_alpha, real hyst_stroke, const dmat* dmcmd){
	hyst_t* hyst=hyst_new(hysteresis, hyst_alpha, hyst_stroke, NX(dmcmd));
	dmat* dmreal=dnew(NX(dmcmd), NY(dmcmd));
	hyst_dmat(hyst, dmreal, dmcmd);
	hyst_free(hyst);
	return dmreal;
}

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

   @param[in] hysteresis The ratio of the maximum difference of ascending/descending curve over 2*stroke. Between 0 and 1.
   @param[in] alpha The alpha parameter. Negative to use the default of 2.
   @param[in] stroke The DM command is from -stroke to stroke (m). Infinity is replaced by 10e-6.
   @param[in] nact Number of actuators
   @return Newly allocated hysteresis model. Free with hyst_free().
*/
hyst_t* hyst_new(real hysteresis, real alpha, real stroke, long nact){
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

   @param[in,out] hyst Hysteresis model. May be NULL.
*/
void hyst_reset(hyst_t* hyst){
	if(!hyst) return;
	dzero(hyst->xlast);
	dzero(hyst->ylast);
}

/**
   Free hysteresis model.

   @param[in] hyst Hysteresis model to free. May be NULL.
*/
void hyst_free(hyst_t* hyst){
	if(!hyst) return;
	dfree(hyst->xlast);
	dfree(hyst->ylast);
	free(hyst);
}

/**
   Apply hysteresis to DM vector. dmreal and dmcmd may be the same

   @param[in,out] hyst Hysteresis model. Its state is updated. If NULL, dmcmd is copied to dmreal.
   @param[out] dmreal Actual DM position of nact*ntime. Must be preallocated.
   @param[in] dmcmd DM commands of nact*ntime. Each column is a time step.
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

   @param[in,out] hyst Array of hysteresis models, one per DM. If NULL, dmcmd is copied to dmreal.
   @param[out] dmreal Actual DM positions. Must be preallocated.
   @param[in] dmcmd DM commands.
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
   Test hysteresis model. Dmcmd has dimensions of nact*ntime

   @param[in] hysteresis The ratio of the maximum difference of ascending/descending curve over 2*stroke.
   @param[in] hyst_alpha The alpha parameter.
   @param[in] hyst_stroke The DM command is from -stroke to stroke.
   @param[in] dmcmd DM commands of nact*ntime.
   @return Newly allocated actual DM position of nact*ntime. The caller frees it.
*/
dmat* hyst_test(real hysteresis, real hyst_alpha, real hyst_stroke, const dmat* dmcmd){
	hyst_t* hyst=hyst_new(hysteresis, hyst_alpha, hyst_stroke, NX(dmcmd));
	dmat* dmreal=dnew(NX(dmcmd), NY(dmcmd));
	hyst_dmat(hyst, dmreal, dmcmd);
	hyst_free(hyst);
	return dmreal;
}

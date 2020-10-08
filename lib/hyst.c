/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>

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
	\page algorithm

	\section hysteresis DM Hysteresis

   The DM hysteresis modeling is based on "Modeling the hysteresis of a scanning probe microscope"
   J. Vac. Sci. Technol. B 18 (2), Mar/Apr 2000

   Formula (2) is used with  V replaced by command x, and x replaced by actual DM position y
	\f[
	   \frac{dy}{dx}=\alpha sign(dx) (\beta x-y)+(\beta-u)
	\f]


	For closed hyteresis loop within the stroke limit \f$S\f$, the linear ordinary differential equation can be solved analytically:
	\f{eqnarray*}{
		y_{t/r}(x)=\beta x \mp \frac{u}{\alpha}\left(1-\frac{2\exp(\mp\alpha x)}{\exp(-\alpha S)+\exp(\alpha S)}\right)
	\f}

	where t and r represents trace and retrace directions. It is desirable to have y be
	within the same range as x to simulate calibrated command:

	\f$y_{t/r}(S)=S\f$. Let the hysteresis ratio be \f$h\f$, we have

	\f{eqnarray*}{
		\frac{y_r(0)}{S}&=&\frac{u}{\alpha S}\left(1-\frac{2}{\exp(-\alpha S)+exp(\alpha S)}\right)=h \\
		\frac{y_{t/r}(S)}{S}&=&\beta-\frac{u}{\alpha S}\left(1-\frac{2\exp(-\alpha S)}{\exp(-\alpha S)+\exp(\alpha S)}\right) = 1
	\f}

	The equation can be solved uniquely when \f$\alpha S\f$ is fixed.
	The hysteresis curve is a weak function of \f$\alpha S\f$ and having \f$\alpha S=2\f$ produces a good looking curve.
	Smaller \f$\alpha S\f$ has slower convergence of the hysteresis curve.

	The hysteresis model is simulated using simple integration. For each step with a new x, the position y is updated as
	\f{eqnarray*}{
	   dx &=& x-x_{last} \\
	   dy &=& \alpha*abs(dx)*(\beta*x_{last}-y_{last})+dx*(\beta-u) \\
	   y_{last}&\mathrel{+}=&dy \\
	   x_{last}&=&x \\
	   y&=&y_{last}
	\f}

*/

/**
   contains data related to DM hysterisis modeling for all the common DMs (not
MOAO). let input be x, and output of each mode be y.  */
struct HYST_T{
	real alpha;       /**<Material constant*/
	real beta;        /**<constants*/
	real beta_m_u;    /**<constants, beta-u*/
	dmat* xlast;      /**<Record x from last time step*/
	dmat* ylast;      /**<Record y from last time step*/

};

/**
  Hysteresis modeling. Create the hysteresis model. It assumes that the hysteresis is measured with a stroke.
*/
HYST_T* hyst_new(
	real hysteresis,   /**<The ratio of the maximum difference of ascending/descending curve over 2*stroke. */
	real alpha,        /**<The alpha parameter*/
	real stroke,       /**<The DM command is from -stroke to stroke*/
	long nact          /**<Number of actuators*/
){
	HYST_T* hyst=mycalloc(1, HYST_T);
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
	info("hyst_new: stroke=%g, alpha=%g, beta=%g, u=%g, nact=%ld\n", stroke, alpha, hyst->beta, u, nact);
	return hyst;
}

/**
   Reset hysteresis state
*/
void hyst_reset(HYST_T* hyst){
	dzero(hyst->xlast);
	dzero(hyst->ylast);
}

/**
   Free hysteresis model.
*/
void hyst_free(HYST_T* hyst){
	dfree(hyst->xlast);
	dfree(hyst->ylast);
	free(hyst);
}

/**
   Apply hysteresis to DM vector. dmreal and dmcmd may be the same
*/
void hyst_dmat(HYST_T* hyst, dmat* dmreal, const dmat* dmcmd){
	real* restrict xlast=hyst->xlast->p;
	real* restrict ylast=hyst->ylast->p;
	assert(dmcmd->nx==hyst->xlast->nx);
	for(int it=0; it<dmcmd->ny; it++){
		for(int ia=0; ia<dmcmd->nx; ia++){
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
void hyst_dcell(HYST_T** hyst, dcell* dmreal, const dcell* dmcmd){
	if(!hyst) return;
	for(int idm=0; idm<dmcmd->nx*dmcmd->ny; idm++){
		hyst_dmat(hyst[idm], dmreal->p[idm], dmcmd->p[idm]);
	}
}

/**
 * Test hysteresis model. Dmcmd has dimensions of nact*ntime
 * */
dmat* hyst_test(real hysteresis, real hyst_alpha, real hyst_stroke, const dmat* dmcmd){
	HYST_T* hyst=hyst_new(hysteresis, hyst_alpha, hyst_stroke, dmcmd->nx);
	dmat* dmreal=dnew(dmcmd->nx, dmcmd->ny);
	hyst_dmat(hyst, dmreal, dmcmd);
	hyst_free(hyst);
	return dmreal;
}

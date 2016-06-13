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

#include "hyst.h"

/**
   contains data related to DM hysterisis modeling for all the common DMs (not
MOAO). let input be x, and output of each mode be y.  */
struct HYST_T{
    dmat *coeff;      /**<contains data from parms->dm.hyst*/
    dmat *xlast;      /**<Record x from last time step*/
    dmat *ylast;      /**<Record y from last time step*/
    dmat *dxlast;     /**<Change of x*/
    dmat *x0;         /**<Initial x (before dx changes sign)*/
    dmat *y0;         /**<Initial y*/
    double stroke;    /**<Surface stroke of DM. coeff->alpha is scaled by this*/
    double power;     /**<Slope of changing alpha*/
    double xscale;    /**<Scaling of x to fix overal slope error.*/
};

/**
  Hysteresis modeling. Create the hysteresis model
*/
HYST_T *hyst_new(dmat *coeff, int naloc){
    HYST_T *hyst=calloc(1, sizeof(HYST_T));
    int nhmod=coeff->ny;
    if(coeff->nx!=3 || nhmod<1){
	error("DM hystereis has wrong format. Expect 3 rows, but has %ldx%ld\n", coeff->nx, coeff->ny);
    }
    hyst->coeff=dref(coeff);
    hyst->xlast=dnew(naloc,1);
    hyst->ylast=dnew(nhmod,naloc);
    hyst->dxlast=dnew(naloc,1);
    hyst->x0=dnew(naloc,1);
    hyst->y0=dnew(nhmod,naloc);
    hyst->xscale=1;
    if(coeff->header){
	double stroke=search_header_num(coeff->header, "stroke");
	if(isfinite(stroke)){
	    warning("stroke=%g\n", stroke);
	    for(int i=0; i<coeff->ny; i++){
		coeff->p[i*coeff->nx+1]*=stroke;
	    }
	    hyst->stroke=stroke;
	}
	double alpha_power=search_header_num(coeff->header, "alpha_power");
	if(isfinite(alpha_power)){
	    warning("alpha_power=%g\n", alpha_power);
	    hyst->power=alpha_power;
	}
    }
    return hyst;
}

/**
   Reset hysteresis state
*/
void hyst_reset(HYST_T *hyst){
    dzero(hyst->xlast);
    dzero(hyst->ylast);
    dzero(hyst->dxlast);
    dzero(hyst->x0);
    dzero(hyst->y0);
}

/**
   Free hysteresis model.
*/
void hyst_free(HYST_T *hyst){
    dfree(hyst->coeff);
    dfree(hyst->xlast);
    dfree(hyst->ylast);
    dfree(hyst->dxlast);
    dfree(hyst->x0);
    dfree(hyst->y0);
    free(hyst);
}

/**
   Apply hysteresis to DM vector.
*/
void hyst_dmat(HYST_T *hyst, dmat *dmreal, const dmat *dmcmd){
    double *restrict x=dmcmd->p;
    double *restrict xout=dmreal->p;
    double *restrict xlast=hyst->xlast->p;
    double *restrict dxlast=hyst->dxlast->p;
    double *restrict x0=hyst->x0->p;
    dmat*  ylast=hyst->ylast;
    dmat*  py0=hyst->y0;
    dmat*  coeff=hyst->coeff;
    int nmod=hyst->coeff->ny;
    int naloc=dmcmd->nx;
    for(int ia=0; ia<naloc; ia++){
	double xia=x[ia]*hyst->xscale;
	double dx=xia-xlast[ia];
	if(fabs(dx)>1e-14){/*There is change in command */
	    if(dx*dxlast[ia]<0){
		/*Changes in moving direction, change the initial condition */
		x0[ia]=xlast[ia];
		for(int imod=0; imod<nmod; imod++){
		    IND(py0,imod,ia)=IND(ylast,imod,ia);
		}
	    }
	    double alphasc=dx>0?1:-1;/*To revert the sign of alpha when dx<0 */
	    if(hyst->power){
		alphasc*=pow((hyst->stroke-x[ia])/(hyst->stroke*2.), hyst->power);
	    }
	    for(int imod=0; imod<nmod; imod++){
		const double alpha=alphasc*IND(coeff,1,imod);
		const double alphabeta=alpha*IND(coeff,2,imod);
		IND(ylast,imod,ia)=xia-alphabeta+(IND(py0,imod,ia)-x0[ia]+alphabeta)*exp(-(xia-x0[ia])/alpha);
	    }
	    xlast[ia]=xia;
	    dxlast[ia]=dx;
	}/*else: no change in voltage, no change in output. */
	/*update output. */
	double y=0;
	for(int imod=0; imod<nmod; imod++){
	    y+=IND(ylast,imod,ia)*IND(coeff,0,imod);
	}
	xout[ia]=y;
    }
}

/**
   Apply hysteresis to set of DM vectors
*/
void hyst_dcell(HYST_T **hyst, dcell *dmreal, const dcell *dmcmd){
    if(!hyst) return;
    for(int idm=0; idm<dmcmd->nx*dmcmd->ny; idm++){
	hyst_dmat(hyst[idm], dmreal->p[idm], dmcmd->p[idm]);
    }
}

/*Calibrate the slope of hysteresis curve.*/
void hyst_calib(HYST_T *hyst, int ih){
    int N=1000;
    dmat *cmd=dnew(N,1);
    dmat *real=dnew(N,1);
    dmat *cmd0=dnew(1,1);
    dmat *real0=dnew(1,1);
    double stroke=hyst->stroke;
    if(!hyst->stroke){
	stroke=10e-6;
    }
    for(int i=0; i<N; i++){
	cmd0->p[0]=cmd->p[i]=sin(i*M_PI*0.01)*stroke;
	hyst_dmat(hyst, real0, cmd0);
	real->p[i]=real0->p[0];
    }
    hyst_reset(hyst);
    hyst->xscale*=(dmax(cmd)-dmin(cmd))/(dmax(real)-dmin(real));
    warning("x would be multiplied by %g\n", hyst->xscale);
    writebin(real, "hyst%d_real", ih);
    writebin(cmd,  "hyst%d_cmd", ih);
    for(int i=0; i<N; i++){
	cmd0->p[0]=cmd->p[i];
	hyst_dmat(hyst, real0, cmd0);
	real->p[i]=real0->p[0];
    }
    hyst_reset(hyst);
    writebin(real, "hyst%d_realcalib", ih);
    dfree(cmd);
    dfree(real);
    dfree(cmd0);
    dfree(real0);
}

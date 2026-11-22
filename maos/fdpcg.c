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
/**
   \file maos/fdpcg.c
   Fourier Domain Preconditioner for the tomography. The implementation is in
   lib/fdpcg.c. This file collects its inputs from the parameters and the
   reconstructor.
*/
#include "common.h"
#include "fdpcg.h"
/**
   Prepare data for Tomography Fourier Domain Preconditioner.

   @param[in] parms Simulation parameters.
   @param[in] recon Reconstructor.
   @param[in] powfs WFS type data.
   @param[in] atm Atmosphere screens. Used to provide wind velocity information for delay compensation. May be NULL.
   @return Newly allocated preconditioner data. Free with fdpcg_free().
*/
fdpcg_t* fdpcg_prepare(const parms_t* parms, const recon_t* recon, const powfs_t* powfs, mapcell* atm){
	TIC;tic;
	int hipowfs=-1;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(!parms->powfs[ipowfs].lo&&!parms->powfs[ipowfs].skip){
			if(hipowfs==-1){
				hipowfs=ipowfs;
			} else{
				error("Multiple High order powfs found\n");
			}
		}
	}
	loc_t* saloc=powfs[hipowfs].saloc;
	const long nps=recon->npsr;
	long pos=parms->tomo.pos;
	const real* os=P(parms->atmr.os);
	if(pos!=os[0]){
		dbg("ploc over sampling does not equal to ground layer over sampling.\n");
	}
	long* nx=P(recon->xnx);
	long* ny=P(recon->xny);
	const long nxp=nx[0]/os[0]*parms->tomo.pos;
	const long nyp=ny[0]/os[0]*parms->tomo.pos;
	const real dxp=recon->ploc->dx;
	const real* ht=P(parms->atmr.ht);
	long nxtot=0;
	for(long ips=0; ips<nps; ips++){
		nxtot+=nx[ips]*ny[ips];
		if(os[ips]>pos){
			dbg("Layer %ld os is greater than ground\n", ips);
		}
	}
	/*Concatenate invpsd; */
	cmat* invpsd=cnew(nxtot, 1);
	long offset=0;
	switch(parms->tomo.cxxalg){
	case 0:/*forward matrix uses biharmonic approx. We use here also. */
		for(long ips=0; ips<nps; ips++){
			cmat* psd=cnew(nx[ips], ny[ips]);
			//cfft2plan(psd,-1);
			dsp* L2;
			if(parms->tomo.square){
				L2=dspref(P(recon->L2,ips,ips));
			} else{/*L2 is for non square xloc. need to build L2 for square xloc. */
				L2=mklaplacian_map(nx[ips], ny[ips],
					P(recon->xloc,ips)->dx, recon->r0,
					P(recon->wt,ips));
				dspscale(L2, sqrt(parms->tomo.cxxscale*TOMOSCALE));
			}
			dsp* tmp=dspmulsp(L2, L2, "tn");
			dspfree(L2);
			for(long irow=tmp->pp[0]; irow<tmp->pp[1]; irow++){/*first column of tmp to psf. */
				P(psd,tmp->pi[irow])=tmp->px[irow];
			}
			dspfree(tmp);
			cfft2(psd, -1);
			cfftshift(psd);
			/*look for a way to obtain this automatically. */
			const real eps=2.220446049250313e-16;
			real max=cmax(psd)*sqrt(eps);
			for(long i=0; i<nx[ips]*ny[ips]; i++){
				P(invpsd, offset+i)=creal(P(psd,i))+max;
			}
			offset+=nx[ips]*ny[ips];
			cfree(psd);
		}
		break;
	case 1:
	case 2:
	/*forward matrix uses inverse PSD or fractal. we use PSD here. */
		for(long ips=0; ips<nps; ips++){
			dmat* tmp=ddup(P(recon->invpsd->invpsd,ips));
			dfftshift(tmp);
			/*cancel the scaling applied in invpsd routine. */
			dscale(tmp, (real)(nx[ips]*ny[ips]));
			for(long i=0; i<nx[ips]*ny[ips]; i++){
				P(invpsd, offset+i)=P(tmp,i);
			}
			offset+=nx[ips]*ny[ips];
			dfree(tmp);
		}
		break;
	}
	//toc2("fdpcg_prepare:invspd");
	const real delay=parms->sim.dt*(parms->powfs[hipowfs].dtrat+1+parms->sim.alhi);
	const long nwfs=parms->powfs[hipowfs].nwfsr;
	dmat* dispx=dnew(nps, nwfs);
	dmat* dispy=dnew(nps, nwfs);
	real* neai=mycalloc(nwfs, real);
	long* iwfss=mycalloc(nwfs, long);
	for(int jwfs=0; jwfs<nwfs; jwfs++){
		int iwfs=(int)P(parms->powfs[hipowfs].wfsr,jwfs);
		iwfss[jwfs]=iwfs;
		neai[jwfs]=P(recon->neam,iwfs);
		for(long ips=0; ips<nps; ips++){
			/*
			  2010-05-28: The cone effect cancels with the cone
			  coordinate. so we are like doing parallel beam
			  propagation. Removed scaling by 1/(1-ht[ips]/hs);
			*/
			P(dispx, ips, jwfs)=ht[ips]*parms->wfsr[iwfs].thetax;
			P(dispy, ips, jwfs)=ht[ips]*parms->wfsr[iwfs].thetay;
			if(atm){
				int ips0=(int)P(parms->atmr.indps,ips);
				P(dispx, ips, jwfs)+=P(atm,ips0)->vx*delay;
				P(dispy, ips, jwfs)+=P(atm,ips0)->vy*delay;
			}
		}
	}
	csp* Mhat=fdpcg_mhat(nps, nx, ny, nxp, nyp, dxp, saloc, P(powfs[hipowfs].saamax), invpsd, nwfs, dispx, dispy, neai, iwfss, parms->save.fdpcg);
	cfree(invpsd);
	free(neai);
	free(iwfss);
	dfree(dispx);
	dfree(dispy);
	fdpcg_t* fdpcg=fdpcg_new(Mhat, nps, nx, ny, os, parms->gpu.tomo, recon->xloc, parms->tomo.square, parms->save.fdpcg);
	toc2("fdpcg_prepare");
	return fdpcg;
}
/**
   Apply the Fourier domain preconditioner with fdpcg_apply(). Used as the
   preconditioner of the tomography CG.

   @param[out] xout Output. Allocated if NULL. The content is replaced.
   @param[in] A The reconstructor (recon_t).
   @param[in] xin Input.
*/
void fdpcg_precond(dcell** xout, const void* A, const dcell* xin){
	const recon_t* recon=(const recon_t*)A;
	fdpcg_apply(xout, recon->fdpcg, xin, recon->nthread);
}

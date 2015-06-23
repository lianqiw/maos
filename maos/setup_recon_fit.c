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

#include "common.h"
#include "setup_recon.h"
#include "recon_utils.h"

/**
   Setup ray tracing operator HXF from xloc to aperture ploc along DM fiting directions*/
static void
setup_recon_HXF(RECON_T *recon, const PARMS_T *parms){
    if(parms->load.HXF && zfexist(parms->load.HXF)){
	warning("Loading saved HXF\n");
	recon->HXF=dspcellread("%s",parms->load.HXF);
    }else{
	info2("Generating HXF");TIC;tic;
	const int nfit=parms->fit.nfit;
	const int npsr=recon->npsr;
	recon->HXF=cellnew(nfit, npsr);
	PDSPCELL(recon->HXF,HXF);
	for(int ifit=0; ifit<nfit; ifit++){
	    double hs=parms->fit.hs->p[ifit];
	    for(int ips=0; ips<npsr; ips++){
		const double ht = recon->ht->p[ips];
		const double scale=1.-ht/hs;
		double displace[2];
		displace[0]=parms->fit.thetax->p[ifit]*ht;
		displace[1]=parms->fit.thetay->p[ifit]*ht;
		HXF[ips][ifit]=mkh(recon->xloc->p[ips], recon->floc, 
				   displace[0], displace[1], scale,
				   parms->tomo.iac);
	    }
	}
	if(parms->save.setup){
	    writebin(recon->HXF,"HXF");
	}
	toc2(" ");
    }
}

/**
   Assemble the DM fitting matrix

   The fitting is done by minimizing \f$||H_X x - H_A a||^2_W\f$ where \f$H_X,
   H_A\f$ are ray tracing operator from tomography grid xloc, and deformable
   mirror grid aloc to pupil grid ploc. The norm is weighted using bilinear
   influence functions within the telescope aperture. We have
   
   \f$a=\left[H_A^T(W_0-W_1 W_1^T)H_A\right]^{-1} H_A^T (W_0-W_1) H_X x\f$

   For details see www.opticsinfobase.org/abstract.cfm?URI=josaa-19-9-1803 
*/
static void
setup_recon_fit_matrix(RECON_T *recon, const PARMS_T *parms){
    const int nfit=parms->fit.nfit;
    const int ndm=parms->ndm;
    if(ndm==0) return;
    dspcell *HATc=dspcelltrans(recon->HA);
    PDSPCELL(HATc, HAT);
    PDSPCELL(recon->HA,HA);

    info2("Before assembling fit matrix:\t%.2f MiB\n",get_job_mem()/1024.);
    /*Assemble Fit matrix. */
    int npsr=recon->npsr;
    if(parms->load.fit){
	if(!(zfexist("FRM") && zfexist("FRU") && zfexist("FRV"))){
	    error("FRM, FRU, FRV (.bin) not all exist\n");
	}
	warning("Loading saved recon->FR\n");
	recon->FR.M=readbin("FRM");
	recon->FR.U=dcellread("FRU");
	recon->FR.V=dcellread("FRV");
    }else{
	if(recon->HXF){
	    info2("Building recon->FR\n");
	    recon->FR.M=cellnew(ndm, npsr);
	    PDSPCELL(recon->FR.M, FRM);
	    PDSPCELL(recon->HXF, HXF);

	    for(int ips=0; ips<npsr; ips++){
		for(int ifit=0; ifit<nfit; ifit++){
		    if(fabs(recon->fitwt->p[ifit])<1.e-12) continue;
		    dsp *tmp=dspmulsp(recon->W0, HXF[ips][ifit],"nn");
		    for(int idm=0; idm<ndm; idm++){
			dspmulsp2(&FRM[ips][idm],HAT[ifit][idm], tmp, "nn",
				  recon->fitwt->p[ifit]);
		    }
		    dspfree(tmp);
		}
	    }
	    recon->FR.V=cellnew(npsr, 1);
	    dmat **FRV=recon->FR.V->p;  
	
	    for(int ips=0; ips<npsr; ips++){
		int nloc=recon->xloc->p[ips]->nloc;
		FRV[ips]=dnew(nloc,nfit);
		for(int ifit=0; ifit<nfit; ifit++){
		    /*notice the sqrt. */
		    if(fabs(recon->fitwt->p[ifit])<1.e-12) continue;
		    dspmulvec(FRV[ips]->p+ifit*nloc, 
			      HXF[ips][ifit], recon->W1->p, 't',
			      sqrt(recon->fitwt->p[ifit]));
		}
	    }
	    if(parms->save.recon){
		writebin(recon->FR.M,"FRM");
		writebin(recon->FR.V,"FRV");
	    }
	}else{
	    info("Avoid building recon->FR.M\n");
	    recon->FR.M=NULL;
	    recon->FR.V=NULL;
	}
	/*Always need FR.U as it is used to do FL.U, FL.V */
	recon->FR.U=cellnew(ndm, 1);
	dmat **FRU=recon->FR.U->p;
	
	for(int idm=0; idm<ndm; idm++){    
	    int nloc=recon->aloc->p[idm]->nloc;
	    FRU[idm]=dnew(nloc, nfit);
	    for(int ifit=0; ifit<nfit; ifit++){
		/*notice the sart. */
		if(fabs(recon->fitwt->p[ifit])<1.e-12) continue;
		dspmulvec(FRU[idm]->p+ifit*nloc, 
			  HA[idm][ifit], recon->W1->p,'t',
			  sqrt(recon->fitwt->p[ifit]));
	    }
	}
	if(parms->save.recon){
	    writebin(recon->FR.U,"FRU");
	}
    }

    if(parms->load.fit){
	if(!(zfexist("FLM") && zfexist("FLU") && zfexist("FLV"))){
	    error("FLM, FLU, FLV (.bin) not all exist\n");
	}
	warning("Loading saved recon->FL\n");
	recon->FL.M=readbin("FLM");
	recon->FL.U=dcellread("FLU");
	recon->FL.V=dcellread("FLV");
    }else{
	info2("Building recon->FL\n");
	recon->FL.M=cellnew(ndm, ndm);
	dsp *(*FLM)[ndm]=(dsp*(*)[ndm])recon->FL.M->p;
	for(int idm=0; idm<ndm; idm++){
	    for(int ifit=0; ifit<nfit; ifit++){
		if(fabs(recon->fitwt->p[ifit])<1.e-12) continue;
		dsp *tmp=dspmulsp(recon->W0, HA[idm][ifit],"nn");
		for(int jdm=0; jdm<ndm; jdm++){
		    dspmulsp2(&FLM[idm][jdm],HAT[ifit][jdm], tmp,"nn",
			      recon->fitwt->p[ifit]);
		}
		dspfree(tmp);
	    }
	}
	dspcellfree(HATc);
	if(fabs(parms->fit.tikcr)>1.e-15){
	    double tikcr=parms->fit.tikcr;
	    /*Estimated from the formula.  1/nloc is due to W0, the other
	      scaling is due to ray tracing between different sampling freq.*/
	    int nact=0;
	    for(int idm=0; idm<parms->ndm; idm++){
		nact+=recon->aloc->p[idm]->nloc;
	    }
	    double maxeig=4./nact;
	    info2("Adding tikhonov constraint of %g to FLM\n", tikcr);
	    info2("The maximum eigen value is estimated to be around %e\n", maxeig);
	    dcelladdI(recon->FL.M,tikcr*maxeig);
	}

	{/*Low rank terms. */
	    recon->FL.U=dcellcat_each(recon->FR.U, recon->fitNW, 2);
	    dcell *tmp=NULL;/*negative NW. */
	    dcelladd(&tmp, 1, recon->fitNW, -1);
	    recon->FL.V=dcellcat_each(recon->FR.U, tmp, 2);
	    dcellfree(tmp);
	}
	if(recon->actslave){
	    dcelladd(&recon->FL.M, 1, recon->actslave, 1);
	}
	/*dspcellsym(recon->FL.M); */
	info2("DM Fit number of Low rank terms: %ld in RHS, %ld in LHS\n",
	      recon->FR.U->p[0]->ny, recon->FL.U->p[0]->ny);
	if(parms->save.recon){
	    writebin(recon->FL.M,"FLM.bin");
	    writebin(recon->FL.U,"FLU");
	    writebin(recon->FL.V,"FLV");
	}
    }
    if((parms->fit.alg==0 || parms->fit.alg==2) && parms->fit.bgs){
	muv_direct_diag_prep(&(recon->FL),(parms->fit.alg==2)*parms->fit.svdthres);
    }
    if((parms->fit.alg==0 || parms->fit.alg==2) && !parms->fit.bgs){
	if(fabs(parms->fit.tikcr)<1.e-14){
	    warning("tickcr=%g is too small, chol may fail.\n", parms->fit.tikcr);
	}
	muv_direct_prep(&(recon->FL),(parms->fit.alg==2)*parms->fit.svdthres);
	info2("After cholesky/svd on matrix:\t%.2f MiB\n",get_job_mem()/1024.);
    }
 
    if(parms->save.recon){
       	if(recon->FL.C){
	    chol_convert(recon->FL.C, 1);
	    chol_save(recon->FL.C,"FLC.bin");
	}
	if(recon->FL.MI)
	    writebin(recon->FL.MI,"FLMI");
	if(recon->FL.Up)
	    writebin(recon->FL.Up, "FLUp");
	if(recon->FL.Vp)
	    writebin(recon->FL.Vp, "FLVp");
	if(recon->FL.CB){
	    for(int ib=0; ib<recon->FL.nb; ib++){
		chol_save(recon->FL.CB[ib],"FLCB_%d.bin", ib);
	    }
	}
	if(recon->FL.MIB){
	    writebin(recon->FL.MIB,"FLMIB");
	}
    }
    info2("After assemble fit matrix:\t%.2f MiB\n",get_job_mem()/1024.);
}

void setup_recon_fit(RECON_T *recon, const PARMS_T *parms){
    TIC;tic;
    if(!parms->sim.idealfit){
	/*In idealfit, xloc has high sampling. We avoid HXF. */
	setup_recon_HXF(recon,parms);
    }
    /*copy over fitwt since we need a dmat */
    int nfit=parms->fit.nfit;
    recon->fitwt=dnew(nfit,1);
    dcp(&recon->fitwt,parms->fit.wt);
    
    /*always assemble fit matrix, faster if many directions */
    if(parms->fit.assemble){
	setup_recon_fit_matrix(recon,parms);
    }
    toc2("Generating fit matrix ");
    /*Fall back function method if FR.M is NULL (!HXF<-idealfit) */
    recon->FR.Mfun  = FitR;
    recon->FR.Mdata = recon;
    /*Fall back function method if FL.M is NULL */
    recon->FL.Mfun  = FitL;
    recon->FL.Mdata = recon;
    recon->FL.alg = parms->fit.alg;
    recon->FL.bgs = parms->fit.bgs;
    recon->FL.warm  = parms->recon.warm_restart;
    recon->FL.maxit = parms->fit.maxit;
}


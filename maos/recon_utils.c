/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#include "recon_utils.h"

/**
   \file recon_utils.c
   Reusable utilities for wavefront reconstruction and DM fitting.
*/

/**
   Apply Laplacian2 to xin and accumulates to xout.
*/
void apply_L2(dcell **xout, spcell *L2, const dcell *xin, 
	      double alpha, int nthread){
    dcell *xx=NULL;
    spcellmulmat_thread(&xx, L2, xin, 1.,nthread);
    sptcellmulmat_thread(xout, L2, xx, alpha, nthread);
    dcellfree(xx);
}
/**
   Apply turbulence invpsd to x in Fourier space, scaled by alpha.
*/
void apply_invpsd(dcell **xout, INVPSD_T *extra,
		  //long **xembed, dcell *invpsd, ccell *fftxopd, 
		  const dcell *xin, double alpha){
    long **xembed=extra->xembed;
    dcell *invpsd=extra->invpsd;
    ccell *fftxopd=extra->fftxopd;
    static int icall=0;
    icall++;
    for(int ips=0; ips<xin->nx*xin->ny; ips++){
	dmat *xini;
	long nx=fftxopd->p[ips]->nx;
	long ny=fftxopd->p[ips]->ny;
	if(xembed){
	    xini=dnew(nx, ny);
	    embed_in(xini->p, xin->p[ips]->p, xin->p[ips]->nx, xembed[ips]);
	}else{
	    xini=dref_reshape(xin->p[ips], nx, ny);
	}
	ccpd(&fftxopd->p[ips], xini);
	dfree(xini);
	cfft2(fftxopd->p[ips],-1);
	ccwmd(fftxopd->p[ips], invpsd->p[ips], 1);
	cfft2(fftxopd->p[ips],1);
	dmat *xouti=NULL;
	if(xembed){
	    creal2d(&xouti,1,fftxopd->p[ips],alpha);
	    embed_out(xouti->p, (*xout)->p[ips]->p, 
		      (*xout)->p[ips]->nx, xembed[ips]);
	}else{
	    xouti=dref_reshape((*xout)->p[ips], nx, ny);
	    creal2d(&xouti,1,fftxopd->p[ips],alpha);
	}
	dfree(xouti);
    }
}
/**
   Removing Tip/Tilt/Focus from LGS grads. TTF is the Tip/tilt/focus modes, and
PTTF is the pseudo inverse of it, weighted by subaperture noise.  */
void TTFR(dcell* x, const dcell *TTF, const dcell *PTTF){
    if(!TTF || !PTTF){
	return;
    }
    dcell *junk=NULL;
    dcellmm(&junk, PTTF, x, "nn", 1);
    dcellmm(&x, TTF, junk, "nn", -1);
    dcellfree(junk);
}
/**
   Apply weighting W0/W1 to a vector. W0*x-W1*(W1'*x)
 */
static void applyWeach(dmat *xin, const dsp *W0, const dmat *W1, const double wt){
    if(!W0 || !W1) {
	warning("W0 or W1 is NULL\n");
	return;
    }
    dmat *xout=NULL;
    dmat *tmp=NULL;
    spmulmat(&xout, W0, xin, wt);
    dmm(&tmp, W1, xin, "tn", -1);
    dmm(&xout,W1, tmp, "nn", wt);
    dcp(&xin, xout);
    dfree(xout); dfree(tmp);
}
/**
   apply weighting W0/W1 with weighting wt for each block.
   W0*xin-W1*(W1'*xin);
*/
void applyW(dcell *xin, const dsp *W0, const dmat *W1, const double *wt){
    if(!W0 || !W1) {
	warning("W0 or W1 is NULL\n");
	return;
    }
    const int nevl=xin->nx;
    for(int iy=0; iy<xin->ny; iy++){
	for(int ievl=0; ievl<nevl; ievl++){
	    int ind=iy*nevl+ievl;
	    applyWeach(xin->p[ind], W0, W1, wt[ievl]);
	}
    }
}

/**
   Compute W0/W1 weighting dot product: \f$A^T(W0 B-W1 (W1^T B))\f$
*/
dcell* calcWmcc(const dcell *A, const dcell *B, const dsp *W0, 
		const dmat *W1, const dmat *wt){

    assert(wt->nx==B->nx && wt->ny==1 && A->nx == B->nx);
    const int nevl=B->nx;
    dcell *res=dcellnew(A->ny, B->ny);
    for(int iy=0; iy<B->ny; iy++){
	for(int ievl=0; ievl<nevl; ievl++){
	    int ind=iy*nevl+ievl;
	    dmat *xout=NULL;
	    dmat *tmp=NULL;
	    spmulmat(&xout, W0, B->p[ind], wt->p[ievl]);
	    dmm(&tmp, W1, B->p[ind], "tn", -1);
	    dmm(&xout, W1, tmp, "nn", wt->p[ievl]);
	    for(int ix=0; ix<A->ny; ix++){
		dmm(&res->p[ix+iy*res->nx], A->p[ix*nevl+ievl], xout, "tn", 1);
	    }
	    dfree(xout);
	    dfree(tmp);
	}
    }
    return res;
}

/**
   Compute slaving actuator regularization. HA is used to compute active
   actuators. If NW is non NULL, orthogonalize it with the slaving
   regularization.  When the actuators are in the NULL space of HA, we want to
   contraint their values to be close to the ones that are active. We put an
   additional term in the fitting matrix to force this. Becaure with it when
   using tip/tilt constraint and cholesky back solve.*/
spcell *act_slaving(loc_t **aloc, spcell *HA, dmat *W1, dcell *NW){
    if(!HA) {
	error("HA is not supplied\n");
    }
    PSPCELL(HA,pHA);
    int ndm=HA->ny;
    dcell *actcplc=dcellnew(ndm, 1);
    for(int idm=0; idm<ndm; idm++){
	int nact=aloc[idm]->nloc;
	for(int ifit=0; ifit<HA->ny; ifit++){
	    sptmulmat(&actcplc->p[idm], pHA[idm][ifit], W1, 1);
	}
	normalize_max(actcplc->p[idm]->p, nact, 1);//bring max to 1;
    }
    double scl=1./pHA[0][0]->m;
    
    int nslavetot=0;
    const double slave_thres=0.1;
    spcell *actslavec=spcellnew(ndm, ndm);//block diagonal.
    PSPCELL(actslavec, actslave);
    
    for(int idm=0; idm<ndm; idm++){
	int  nact     = aloc[idm]->nloc;
	double *actcpl= actcplc->p[idm]->p;
	double *actcpl0 = actcpl-1;
	int  nslave   = 0;
	for(int iact=0; iact<nact; iact++){
	    if(actcpl[iact]<slave_thres){
		nslave++;
	    }
	}
	nslavetot+=nslave;
	info("dm %d: there are %d slave actuators\n", idm, nslave);
	if(nslave==0) {
	    continue;
	}
	loc_create_map_npad(aloc[idm],1);
	long (*map)[aloc[idm]->map->nx]=(void*)aloc[idm]->map->p;
	double ox=aloc[idm]->map->ox;
	double oy=aloc[idm]->map->oy;
	double dx1=1./aloc[idm]->dx;
	dsp *slavet=spnew(nact,nslave,nslave*5);
	spint *pp=slavet->p;
	spint *pi=slavet->i;
	double *px=slavet->x;
	const double *locx=aloc[idm]->locx;
	const double *locy=aloc[idm]->locy;
	long count=0;
	long icol=0;
	for(int iact=0; iact<nact; iact++){
	    if(actcpl[iact]<slave_thres){//slave actuators
		pp[icol]=count;
		long mapx=(long)round((locx[iact]-ox)*dx1);
		long mapy=(long)round((locy[iact]-oy)*dx1);
		if(map[mapy][mapx]-1!=iact){
		    error("mapping is used incorrectly\n");
		}
		int near_active=0;
		int near_exist=0;
		for(int idy=-1; idy<2; idy++){
		    for(int idx=-1; idx<2; idx++){
			if((idx!=0 && idy!=0) || (idx==0 && idy==0)){
			    continue;//skip center and corner
			}
			if(map[mapy+idy][mapx+idx]){
			    near_exist++;
			    if(actcpl0[map[mapy+idy][mapx+idx]]>0.1){
				near_active++;
			    }
			}
		
		    }
		}
		if(near_exist){
		    pi[count]=iact;
		    px[count]=scl;
		    //also limits the strength
		    if(actcpl[iact]<0.1){
			px[count]+=scl*1;
		    }
		   
		    /*
		    if(actcpl[iact]<0.1 && actcpl[iact]>0){
			double scale=-0.001*log(actcpl[iact]);
			if(scale>0.01){
			    scale=0.01;
			}
			px[count]=scl*(1+scale);
			}*/
		    count++;
		}else{
		    error("This is an isolated actuator\n");
		}
		if(near_active>0){
		    double value=-scl/near_active;
		    for(int idy=-1; idy<2; idy++){
			for(int idx=-1; idx<2; idx++){
			    if((idx!=0 && idy!=0) || (idx==0 && idy==0)){
				continue;//skip center and corner
			    }
			    if(map[mapy+idy][mapx+idx] && actcpl0[map[mapy+idy][mapx+idx]]>0.1){
				pi[count]=map[mapy+idy][mapx+idx]-1;
				px[count]=value;
				count++;
			    }
			}
		    }
		}else{
		    double value=-scl/near_exist;
		    for(int idy=-1; idy<2; idy++){
			for(int idx=-1; idx<2; idx++){
			    if((idx!=0 && idy!=0) || (idx==0 && idy==0)){
				continue;//skip center and corner
			    }
			    if(map[mapy+idy][mapx+idx]){
				pi[count]=map[mapy+idy][mapx+idx]-1;
				px[count]=value;
				count++;
			    }
			}
		    }
	
		}
		icol++;
	    }
	}
	pp[icol]=count;
	if(icol!=nslave){
	    error("Doesnot match\n");
	}
	spsetnzmax(slavet, count);

	loc_free_map(aloc[idm]);
	dsp *slave=sptrans(slavet);
	actslave[idm][idm]=spmulsp(slavet, slave);

	if(NW){
	    /*Now we need to make sure NW is in the NULL
	      space of the slaving regularization, especiall
	      the tip/tilt constraints.*/
	    if(NW->p[idm]){
		dmat *H=NULL;
		spfull(&H, slavet, 1);
		dmat *Hinv=dpinv(H,NULL,NULL);
		dmat *mod=NULL;
		dmm(&mod, Hinv, NW->p[idm], "nn", 1);
		dmm(&NW->p[idm], H, mod,"nn", -1);
		dfree(H);
		dfree(Hinv);
		dfree(mod);
	    }
	}
	spfree(slave);
	spfree(slavet);
    }//idm
    dcellfree(actcplc);
    if(nslavetot==0){
	spcellfree(actslavec);
	actslavec=NULL;
    }
    return actslavec;
}

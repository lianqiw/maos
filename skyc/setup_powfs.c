/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#include "skyc.h"
#include "parms.h"
#include "types.h"
#include "setup_powfs.h"
/**
   \file skyc/setup_powfs.c
*/
/**
   Setup the detector transfer functions. See maos/setup_powfs.c
 */
static void setup_powfs_dtf(POWFS_S *powfs, const PARMS_S* parms){
    const int npowfs=parms->maos.npowfs;
    for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
	const int ncomp=parms->maos.ncomp[ipowfs];
	const int ncomp2=ncomp>>1;
	const int embfac=parms->maos.embfac[ipowfs];
	const int pixpsa=parms->skyc.pixpsa[ipowfs];
	const double pixtheta=parms->skyc.pixtheta[ipowfs];
	const double blur=parms->skyc.pixblur[ipowfs]*pixtheta;
	const double e0=exp(-2*M_PI*M_PI*blur*blur);
	const double dxsa=parms->maos.dxsa[ipowfs];
	const double pixoffx=parms->skyc.pixoffx[ipowfs];
	const double pixoffy=parms->skyc.pixoffy[ipowfs];
	const double pxo=-(pixpsa/2-0.5+pixoffx)*pixtheta;
	const double pyo=-(pixpsa/2-0.5+pixoffy)*pixtheta;
	loc_t *loc_ccd=mksqloc(pixpsa, pixpsa, pixtheta, pxo, pyo);
	powfs[ipowfs].dtf=calloc(parms->maos.nwvl, sizeof(DTF_S));
	for(int iwvl=0; iwvl<parms->maos.nwvl; iwvl++){
	    const double wvl=parms->maos.wvl[iwvl];
	    const double dtheta=wvl/(dxsa*embfac);
	    const double pdtheta=pixtheta*pixtheta/(dtheta*dtheta);
	    const double du=1./(dtheta*ncomp);
	    const double du2=du*du;
	    const double dupth=du*pixtheta;
	    cmat *nominal=cnew(ncomp,ncomp);
	    cfft2plan(nominal,-1);
	    cfft2plan(nominal,1);
	    PCMAT(nominal,pn);
	    const double theta=0;
	    const double ct=cos(theta);
	    const double st=sin(theta);
	    for(int iy=0; iy<ncomp; iy++){
		int jy=iy-ncomp2;
		for(int ix=0; ix<ncomp; ix++){
		    int jx=ix-ncomp2;
		    double ir=ct*jx+st*jy;
		    double ia=-st*jx+ct*jy;
		    pn[iy][ix]=sinc(ir*dupth)*sinc(ia*dupth)
			*pow(e0,ir*ir*du2)*pow(e0,ia*ia*du2)
			*pdtheta;
		}
	    }
	    if(parms->skyc.fnpsf1[ipowfs]){
		warning("powfs %d has additional otf to be multiplied\n", ipowfs);
		dcell *psf1c=dcellread(parms->skyc.fnpsf1[ipowfs]);
		dmat *psf1=NULL;
		if(psf1c->nx == 1){
		    psf1=dref(psf1c->p[0]);
		}else if(psf1c->nx==parms->maos.nwvl){
		    psf1=dref(psf1c->p[iwvl]);
		}else{
		    error("skyc.fnpsf1 has wrong dimension\n");
		}
		dcellfree(psf1c);
		if(psf1->ny!=2){
		    error("skyc.fnpsf1 has wrong dimension\n");
		}
		dmat *psf1x=dnew_ref(psf1->nx, 1, psf1->p);
		dmat *psf1y=dnew_ref(psf1->nx, 1, psf1->p+psf1->nx);
		dmat *psf2x=dnew(ncomp*ncomp, 1);
		for(int iy=0; iy<ncomp; iy++){
		    int jy=iy-ncomp2;
		    for(int ix=0; ix<ncomp; ix++){
			int jx=ix-ncomp2;
			psf2x->p[ix+iy*ncomp]=sqrt(jx*jx+jy*jy)*dtheta;
		    }
		}
		info("powfs %d, iwvl=%d, dtheta=%g\n", ipowfs, iwvl, dtheta*206265000);
		dwrite(psf2x, "powfs%d_psf2x_%d", ipowfs,iwvl);
		dmat *psf2=dinterp1(psf1x, psf1y, psf2x);
		normalize_sum(psf2->p, psf2->nx*psf2->ny, 1);
		psf2->nx=ncomp; psf2->ny=ncomp;
		dwrite(psf2, "powfs%d_psf2_%d", ipowfs,iwvl);
		cmat *otf2=cnew(ncomp, ncomp);
		cfft2plan(otf2, -1);
		ccpd(&otf2, psf2);
		cfftshift(otf2);
		cfft2(otf2, -1);
		cfftshift(otf2);
		cwrite(otf2, "powfs%d_otf2_%d", ipowfs, iwvl);
		cwrite(nominal, "powfs%d_dtf%d_nominal_0",ipowfs,iwvl);
		for(int i=0; i<ncomp*ncomp; i++){
		    nominal->p[i]*=otf2->p[i];
		}
		cwrite(nominal, "powfs%d_dtf%d_nominal_1",ipowfs,iwvl);
		dfree(psf1x);
		dfree(psf1y);
		dfree(psf2x);
		dfree(psf1);
		cfree(otf2);
		dfree(psf2);
	    }
	    cfftshift(nominal);
	    cfft2(nominal,-1);
	    cfftshift(nominal);
	    cfft2i(nominal,1);
	    ccp(&powfs[ipowfs].dtf[iwvl].nominal, nominal);
	    cfree(nominal);

	    loc_t *loc_psf=mksqloc(ncomp, ncomp, dtheta, -ncomp2*dtheta, -ncomp2*dtheta);
	    powfs[ipowfs].dtf[iwvl].si=mkh(loc_psf,loc_ccd,NULL,0,0,1,0,0);
	    locfree(loc_psf);
	    if(parms->skyc.dbg){
		cwrite(powfs[ipowfs].dtf[iwvl].nominal,
		       "%s/powfs%d_dtf%d_nominal",dirsetup,ipowfs,iwvl);
		spwrite(powfs[ipowfs].dtf[iwvl].si,
			"%s/powfs%d_dtf%d_si",dirsetup,ipowfs,iwvl);
	    }
	    powfs[ipowfs].dtf[iwvl].U=cnew(ncomp,1);
	    dcomplex *U=powfs[ipowfs].dtf[iwvl].U->p;

	    for(int ix=0; ix<ncomp; ix++){
		int jx=ix<ncomp2?ix:(ix-ncomp);
		U[ix]=-2.*I*M_PI*jx*du;
	    }
	}/*iwvl */
	locfree(loc_ccd);
    }/*ipowfs */
}
/**
   Read in the WFS grid and amplitude map.
 */
static void read_powfs_locamp(POWFS_S *powfs, const PARMS_S *parms){
    const int npowfs=parms->maos.npowfs;
    for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
	powfs[ipowfs].loc=locread("%s",parms->maos.fnwfsloc[ipowfs]);
	powfs[ipowfs].amp=dread("%s",parms->maos.fnwfsamp[ipowfs]);
	powfs[ipowfs].saloc=locread("%s",parms->maos.fnsaloc[ipowfs]);
	powfs[ipowfs].saloc->dx=parms->maos.dxsa[ipowfs];
	const loc_t *loc=powfs[ipowfs].loc;
	const dmat *amp=powfs[ipowfs].amp;
	const long nsa=parms->maos.nsa[ipowfs];
	const long ptspsa=loc->nloc/nsa;
	if(loc->nloc!=nsa*ptspsa){
	    error("loc %ld does not divide to %ld sa\n",loc->nloc, nsa);
	}
	powfs[ipowfs].locxamp=calloc(nsa,sizeof(double));
	powfs[ipowfs].locyamp=calloc(nsa,sizeof(double));
	for(long isa=0; isa<nsa; isa++){
	    const double* iamp=amp->p+isa*ptspsa;
	    const double* locx=loc->locx+isa*ptspsa;
	    const double* locy=loc->locy+isa*ptspsa;
	    double ampsum=0, locxamp=0, locyamp=0;
	    for(long iloc=0; iloc<ptspsa; iloc++){
		ampsum+=iamp[iloc];
		locxamp+=iamp[iloc]*locx[iloc];
		locyamp+=iamp[iloc]*locy[iloc];
	    }
	    powfs[ipowfs].locxamp[isa]=locxamp/ampsum;
	    powfs[ipowfs].locyamp[isa]=locyamp/ampsum;
	}
    }
}
/**
   Setup the corase loc grid corresponding to coarse sampled complex pupil
function.  */
static void setup_powfs_coarseloc(POWFS_S *powfs, const PARMS_S *parms){
    for(long ipowfs=0; ipowfs<parms->maos.npowfs; ipowfs++){
	const long nsa=parms->maos.nsa[ipowfs];
	powfs[ipowfs].cloc=calloc(nsa, sizeof(loc_t*));
	powfs[ipowfs].fpc=calloc(nsa,sizeof(double));

	const long ptspsa=parms->maos.ncomp[ipowfs]/2;
	const double dxsa=parms->maos.dxsa[ipowfs];
	double dx=dxsa/ptspsa;
	for(long isa=0; isa<nsa; isa++){
	    double ox=powfs[ipowfs].saloc->locx[isa]+dx*0.5;
	    double oy=powfs[ipowfs].saloc->locy[isa]+dx*0.5;
	    powfs[ipowfs].cloc[isa]=mksqloc(ptspsa,ptspsa,dx,ox,oy);
	    if(parms->skyc.dbg){
		locwrite(powfs[ipowfs].cloc[isa],"%s/powfs%ld_cloc%ld",
			 dirsetup, ipowfs, isa);
	    }
	    double fpc=0;
	    const double *locx=powfs[ipowfs].cloc[isa]->locx;
	    const double *locy=powfs[ipowfs].cloc[isa]->locy;
	    const long nloc=powfs[ipowfs].cloc[isa]->nloc;
	    for(long iloc=0; iloc<nloc; iloc++){
		fpc+=locx[iloc]*locx[iloc]+locy[iloc]*locy[iloc];
	    }
	    powfs[ipowfs].fpc[isa]=fpc/nloc;
	    /*info("powfs[%ld].fpc[%ld]=%g\n",ipowfs,isa,powfs[ipowfs].fpc[isa]); */
	}
    }
}
/**
   Method to demote TTF. deprecated. 
 */
static void setup_powfs_dettf(POWFS_S *powfs, const PARMS_S *parms){
    for(long ipowfs=0; ipowfs<parms->maos.npowfs; ipowfs++){
	const long nsa=parms->maos.nsa[ipowfs];
	powfs[ipowfs].dettf=dnew(2,nsa*2);	
	if(nsa==1){
	    daddI(powfs[ipowfs].dettf, 1);
	}else{
	    PDMAT(powfs[ipowfs].dettf, ppdettf);
	    double insa=1./nsa;
	    for(long ig=0; ig<2; ig++){
		for(long isa=0; isa<nsa; isa++){
		    ppdettf[isa+nsa*ig][ig]=insa;
		}
	    }
	}
    }
}
/**
   Setup POWFS_S struct.
*/    
void setup_powfs(POWFS_S *powfs, const PARMS_S *parms){
    for(int ipowfs=0; ipowfs<parms->maos.npowfs; ipowfs++){
	powfs[ipowfs].ipowfs=ipowfs;
    }
    setup_powfs_dtf(powfs, parms);
    read_powfs_locamp(powfs, parms);
    setup_powfs_coarseloc(powfs,parms);
    setup_powfs_dettf(powfs,parms);
}
/**
   Release memory.
 */
void free_powfs(POWFS_S *powfs, const PARMS_S *parms){
    const int npowfs=parms->maos.npowfs;
    for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
	for(int iwvl=0; iwvl<parms->maos.nwvl; iwvl++){
	    cfree(powfs[ipowfs].dtf[iwvl].U);
	    cfree(powfs[ipowfs].dtf[iwvl].nominal);
	    spfree(powfs[ipowfs].dtf[iwvl].si);
	}
	free(powfs[ipowfs].dtf);
	locfree(powfs[ipowfs].loc);
	locfree(powfs[ipowfs].saloc);
	locarrfree(powfs[ipowfs].cloc, parms->maos.nsa[ipowfs]);
	free(powfs[ipowfs].fpc);
	dfree(powfs[ipowfs].amp);
	free(powfs[ipowfs].locxamp);
	free(powfs[ipowfs].locyamp);
	dfree(powfs[ipowfs].dettf);
    }
}

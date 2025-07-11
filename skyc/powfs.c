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
#include "skyc.h"
#include "parms.h"
#include "types.h"
#include "powfs.h"
/**
   \file skyc/setup_powfs.c
*/
/**
   Setup the detector transfer functions. See maos/setup_powfs.c
 */
static void setup_powfs_dtf(powfs_s* powfs, const parms_s* parms){
	const int npowfs=parms->maos.npowfs;
	for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
		const int ncomp=parms->maos.ncomp[ipowfs];
		const int ncomp2=ncomp>>1;
		const int embfac=parms->maos.embfac[ipowfs];
		const int pixpsa=parms->skyc.pixpsa[ipowfs];
		const real pixtheta=parms->skyc.pixtheta[ipowfs];
		const real blur=parms->skyc.pixblur[ipowfs]*pixtheta;
		const real e0=exp(-2*M_PI*M_PI*blur*blur);
		const real dxsa=parms->maos.dxsa[ipowfs];
		const real pixoffx=parms->skyc.pixoffx[ipowfs];
		const real pixoffy=parms->skyc.pixoffy[ipowfs];
		const real pxo=-(pixpsa/2-0.5+pixoffx)*pixtheta;
		const real pyo=-(pixpsa/2-0.5+pixoffy)*pixtheta;
		loc_t* loc_ccd=mksqloc(pixpsa, pixpsa, pixtheta, pixtheta, pxo, pyo);
		powfs[ipowfs].dtf=mycalloc(parms->maos.nwvl, dtf_s);
		for(int iwvl=0; iwvl<parms->maos.nwvl; iwvl++){
			const real wvl=parms->maos.wvl[iwvl];
			const real dtheta=wvl/(dxsa*embfac);
			const real pdtheta=pixtheta*pixtheta/(dtheta*dtheta);
			const real du=1./(dtheta*ncomp);
			const real du2=du*du;
			const real dupth=du*pixtheta;
			cmat* nominal=cnew(ncomp, ncomp);
			//cfft2plan(nominal,-1);
			//cfft2plan(nominal,1);
			cmat* pn=nominal;
			const real theta=0;
			const real ct=cos(theta);
			const real st=sin(theta);
			for(int iy=0; iy<ncomp; iy++){
				int jy=iy-ncomp2;
				for(int ix=0; ix<ncomp; ix++){
					int jx=ix-ncomp2;
					real ir=ct*jx+st*jy;
					real ia=-st*jx+ct*jy;
					P(pn, ix, iy)=sinc(ir*dupth)*sinc(ia*dupth)
						*pow(e0, ir*ir*du2)*pow(e0, ia*ia*du2)
						*pdtheta;
				}
			}
			if(parms->skyc.fnpsf1[ipowfs]){
				warning("powfs %d has additional otf to be multiplied\n", ipowfs);
				dcell* psf1c=dcellread("%s", parms->skyc.fnpsf1[ipowfs]);
				dmat* psf1=NULL;
				if(psf1c->nx==1){
					psf1=dref(P(psf1c,0));
				} else if(psf1c->nx==parms->maos.nwvl){
					psf1=dref(P(psf1c,iwvl));
				} else{
					error("skyc.fnpsf1 has wrong dimension\n");
				}
				dcellfree(psf1c);
				if(psf1->ny!=2){
					error("skyc.fnpsf1 has wrong dimension\n");
				}
				dmat* psf1x=drefcols(psf1, 0, 1);
				dmat* psf1y=drefcols(psf1, 1, 1);
				dmat* psf2x=dnew(ncomp*ncomp, 1);
				for(int iy=0; iy<ncomp; iy++){
					int jy=iy-ncomp2;
					for(int ix=0; ix<ncomp; ix++){
						int jx=ix-ncomp2;
						P(psf2x,ix,iy)=sqrt(jx*jx+jy*jy)*dtheta;
					}
				}
				dbg("powfs %d, iwvl=%d, dtheta=%g\n", ipowfs, iwvl, dtheta*RAD2MAS);
				writebin(psf2x, "powfs%d_psf2x_%d", ipowfs, iwvl);
				dmat* psf2=dinterp1(psf1x, psf1y, psf2x, 0);
				dnormalize_sumabs(psf2, 1);
				psf2->nx=ncomp; psf2->ny=ncomp;
				writebin(psf2, "powfs%d_psf2_%d", ipowfs, iwvl);
				cmat* otf2=cnew(ncomp, ncomp);
				//cfft2plan(otf2, -1);
				ccpd(&otf2, psf2);//peak in center
				cfftshift(otf2);//peak in corner
				cfft2(otf2, -1);
				cfftshift(otf2);//peak in center
				writebin(otf2, "powfs%d_otf2_%d", ipowfs, iwvl);
				writebin(nominal, "powfs%d_dtf%d_nominal_0", ipowfs, iwvl);
				for(int i=0; i<ncomp*ncomp; i++){
					P(nominal,i)*=P(otf2,i);
				}
				writebin(nominal, "powfs%d_dtf%d_nominal_1", ipowfs, iwvl);
				dfree(psf1x);
				dfree(psf1y);
				dfree(psf2x);
				dfree(psf1);
				cfree(otf2);
				dfree(psf2);
			}
			cfftshift(nominal);//peak in corner
			cfft2(nominal, -1);
			cfftshift(nominal);//peak in center
			cfft2i(nominal, 1);
			/*This nominal will multiply to OTF with peak in corner. But after
			  inverse fft, peak will be in center*/
			ccp(&powfs[ipowfs].dtf[iwvl].nominal, nominal);
			cfree(nominal);

			loc_t* loc_psf=mksqloc(ncomp, ncomp, dtheta, dtheta, -ncomp2*dtheta, -ncomp2*dtheta);
			powfs[ipowfs].dtf[iwvl].si=mkh(loc_psf, loc_ccd, 0, 0, 1, 0);
			locfree(loc_psf);
			if(parms->skyc.dbg>3){
				writebin(powfs[ipowfs].dtf[iwvl].nominal,
					"%s/powfs%d_dtf%d_nominal", dirsetup, ipowfs, iwvl);
				writebin(powfs[ipowfs].dtf[iwvl].si,
					"%s/powfs%d_dtf%d_si", dirsetup, ipowfs, iwvl);
			}
			powfs[ipowfs].dtf[iwvl].U=cnew(ncomp, 1);
			comp* U=P(powfs[ipowfs].dtf[iwvl].U);

			for(int ix=0; ix<ncomp; ix++){
				int jx=ix<ncomp2?ix:(ix-ncomp);
				U[ix]=COMPLEX(0, -2.*M_PI*jx*du);
			}
		}/*iwvl */
		locfree(loc_ccd);
	}/*ipowfs */
}
/**
   Read in the WFS grid and amplitude map.
 */
static void read_powfs_locamp(powfs_s* powfs, const parms_s* parms){
	const int npowfs=parms->maos.npowfs;
	for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
		powfs[ipowfs].loc=locread("%s", parms->maos.fnwfsloc[ipowfs]);
		powfs[ipowfs].amp=dread("%s", parms->maos.fnwfsamp[ipowfs]);
		powfs[ipowfs].saloc=locread("%s", parms->maos.fnsaloc[ipowfs]);
		powfs[ipowfs].saloc->dx=parms->maos.dxsa[ipowfs];
		powfs[ipowfs].saloc->dy=parms->maos.dxsa[ipowfs];
		const loc_t* loc=powfs[ipowfs].loc;
		const dmat* amp=powfs[ipowfs].amp;
		const long nsa=parms->maos.nsa[ipowfs];
		const long ptspsa=loc->nloc/nsa;
		if(loc->nloc!=nsa*ptspsa){
			error("loc %ld does not divide to %ld sa\n", loc->nloc, nsa);
		}
		powfs[ipowfs].sacent=dnew(nsa, 2);
		for(long isa=0; isa<nsa; isa++){
			const real* iamp=P(amp)+isa*ptspsa;
			const real* locx=loc->locx+isa*ptspsa;
			const real* locy=loc->locy+isa*ptspsa;
			real ampsum=0, locxamp=0, locyamp=0;
			for(long iloc=0; iloc<ptspsa; iloc++){
				ampsum+=iamp[iloc];
				locxamp+=iamp[iloc]*locx[iloc];
				locyamp+=iamp[iloc]*locy[iloc];
			}
			P(powfs[ipowfs].sacent, isa, 0)=locxamp/ampsum;
			P(powfs[ipowfs].sacent, isa, 1)=locyamp/ampsum;
		}
	}
}
/**
   Setup powfs_s struct.
*/
void setup_powfs(powfs_s* powfs, const parms_s* parms){
	for(int ipowfs=0; ipowfs<parms->maos.npowfs; ipowfs++){
		powfs[ipowfs].ipowfs=ipowfs;
		powfs[ipowfs].nxwvf=parms->maos.ncomp[ipowfs]/parms->maos.embfac[ipowfs];
		powfs[ipowfs].dxwvf=parms->maos.dxsa[ipowfs]/powfs[ipowfs].nxwvf;
	}
	setup_powfs_dtf(powfs, parms);
	read_powfs_locamp(powfs, parms);
}
/**
   Release memory.
 */
void free_powfs(powfs_s* powfs, const parms_s* parms){
	const int npowfs=parms->maos.npowfs;
	for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
		for(int iwvl=0; iwvl<parms->maos.nwvl; iwvl++){
			cfree(powfs[ipowfs].dtf[iwvl].U);
			cfree(powfs[ipowfs].dtf[iwvl].nominal);
			dspfree(powfs[ipowfs].dtf[iwvl].si);
		}
		free(powfs[ipowfs].dtf);
		locfree(powfs[ipowfs].loc);
		locfree(powfs[ipowfs].saloc);
		dfree(powfs[ipowfs].amp);
		dfree(powfs[ipowfs].sacent);
	}
}

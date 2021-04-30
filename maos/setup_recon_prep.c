/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>

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
#include "setup_powfs.h"
#include "pywfs.h"
#include "ahst.h"
/*

   Setup grid and ray tracing operators regarding DM. This is independent of
   1) WFS geometry or noise parameters
   2) Tomography
*/

/**
   Setting up PLOC grid, which is a coarse sampled (usually halves the
   subaperture spacing) grid that defines the circular aperture for tomography.*/
static void
setup_recon_ploc(recon_t* recon, const parms_t* parms){
	real dxr=parms->atmr.dx/parms->tomo.pos;/*sampling of ploc */
	if(parms->load.ploc){/*optionally load ploc from the file. see dbg.conf */
		warning("Loading ploc from %s\n", parms->load.ploc);
		recon->ploc=locread("%s", parms->load.ploc);
		if(fabs(recon->ploc->dx-dxr)>dxr*1e-6){
			warning("Loaded ploc has unexpected sampling of %g, should be %g\n",
				recon->ploc->dx, dxr);
		}
	} else{
	/*
	  Create a circular PLOC with telescope diameter by calling
	  create_metapupil with height of 0. We don't add any guard points. PLOC
	  does not need to follow XLOC in FDPCG.*/
		real guard=parms->tomo.guard*dxr;
		map_t* pmap=0;
		create_metapupil(&pmap, 0, 0, parms->dirs, parms->aper.d, 0, dxr, dxr, 0, guard, 0, 0, 0, parms->tomo.square);
		info("PLOC is %ldx%ld, with sampling of %.2fm\n", pmap->nx, pmap->ny, dxr);
		recon->ploc=map2loc(pmap, 0);/*convert map_t to loc_t */
		mapfree(pmap);
	}
	if(parms->save.setup){
		locwrite(recon->ploc, "ploc");
	}
	loc_create_map_npad(recon->ploc, parms->tomo.square?0:1, 0, 0);
	recon->pmap=recon->ploc->map;
	loc_create_stat(recon->ploc);
	if(parms->recon.misreg_tel2wfs){
		real ploc_xm=0, ploc_ym=0;
		for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
			if(parms->recon.misreg_tel2wfs[iwfs]){
				if(!recon->ploc_tel){
					recon->ploc_tel=loccellnew(parms->nwfsr, 1);
					locmean(&ploc_xm, &ploc_ym, recon->ploc);
				}
				recon->ploc_tel->p[iwfs]=loctransform(recon->ploc, parms->recon.misreg_tel2wfs[iwfs]);
				loc_t* ploc2=recon->ploc_tel->p[iwfs];
				real xm, ym;
				locmean(&xm, &ym, ploc2);
				parms->wfsr[iwfs].misreg_x=xm-ploc_xm;
				parms->wfsr[iwfs].misreg_y=ym-ploc_ym;
				parms->wfsr[iwfs].misreg_r=loc_angle(recon->ploc, ploc2);
				warning("ploc for wfs%d shifted by (%g, %g), and rotated by %g.\n", iwfs,
					parms->wfsr[iwfs].misreg_x, parms->wfsr[iwfs].misreg_y,
					parms->wfsr[iwfs].misreg_r);
			}
		}
		if(parms->save.setup){
			writebin(recon->ploc_tel, "ploc_tel");
		}
	}
	if(0){
		for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
			parms->wfsr[iwfs].misreg_r=0.5/15;
			warning("ploc for wfs%d is rotated by %g.\n", iwfs, parms->wfsr[iwfs].misreg_r);
		}
	}
}
/**
   Create loc/amp that can be used to build GP. It has points on edge of subapertures. The amplitude depends on the distortion.
*/
static loc_t*
make_gloc(dmat** gamp, const parms_t* parms, const aper_t* aper, int iwfsr){
	const int ipowfs=parms->wfsr[iwfsr].powfs;
	real dx=parms->powfs[ipowfs].dx;
	loc_t* gloc=mkannloc(parms->aper.d+1, 0, dx, 0);
	if(gamp){
		real dout=parms->aper.d;
		real din=parms->aper.din;
		map_t* ampground=parms->dbg.gp_noamp>0?0:aper->ampground;
		if(parms->dbg.gp_noamp==2){
			dout+=1;
			din=0;
		}
		loc_t* gloc2=0;
		if(ampground&&parms->recon.misreg_tel2wfs&&parms->recon.misreg_tel2wfs[iwfsr]){
			gloc2=loctransform(gloc, parms->recon.misreg_tel2wfs[iwfsr]);
		}
		*gamp=mkamp(gloc2?gloc2:gloc, ampground, 0, 0, dout, din);
		locfree(gloc2);
	}
	return gloc;
}

/**
   Like ploc, but for DM fitting
*/
static void
setup_recon_floc(recon_t* recon, const parms_t* parms){
	real dxr=parms->atmr.dx/parms->fit.pos;/*sampling of floc */
	if(parms->load.floc){
		warning("Loading floc from %s\n", parms->load.floc);
		recon->floc=locread("%s", parms->load.floc);
		if(fabs(recon->floc->dx-dxr)>dxr*1e-6){
			warning("Loaded floc has unexpected sampling of %g, should be %g\n",
				recon->floc->dx, dxr);
		}
	} else{
		real guard=parms->tomo.guard*dxr;
		map_t* fmap=0;
		create_metapupil(&fmap, 0, 0, parms->dirs, parms->aper.d, 0, dxr, dxr, 0, guard, 0, 0, 0, parms->fit.square);
		info("FLOC is %ldx%ld, with sampling of %.2fm\n", fmap->nx, fmap->ny, dxr);
		recon->floc=map2loc(fmap, 0);/*convert map_t to loc_t */
		mapfree(fmap);
		/*Do not restrict fmap to within active pupil. */
	}
	loc_create_map_npad(recon->floc, parms->fit.square?0:1, 0, 0);
	recon->fmap=recon->floc->map;
	/*create the weighting W for bilinear influence function. See [Ellerbroek 2002] */
	if(parms->load.W){
		if(!(zfexist("W0")&&zfexist("W1"))){
			error("W0 or W1 not exist\n");
		}
		warning("Loading W0, W1");
		recon->W0=dspread("W0");
		recon->W1=dread("W1");
	} else{
	/*
	  Compute W0,W1 weighting matrix that can be used to compute piston
	  removed wavefront variance of OPD: RMS=OPD'*(W0-W1*W1')*OPD; W0 is
	  sparse. W1 is a vector. These matrices are used for weighting in DM
	  fitting.
	*/
		real rin=0;
		if(parms->dbg.annular_W&&parms->aper.din>0){
			warning("Define the W0/W1 on annular aperture instead of circular.\n");
			rin=parms->aper.din/2;
		}
		mkw_annular(recon->floc, 0, 0, rin, parms->aper.d/2,
			&(recon->W0), &(recon->W1));
	}
	if(parms->save.setup){
		writebin(recon->W0, "W0");
		writebin(recon->W1, "W1");
	}
	if(parms->save.setup){
		locwrite(recon->floc, "floc");
	}
	loc_create_stat(recon->floc);
}

/**
   Setup the tomography grids xloc which is used for Tomography.
*/
static void
setup_recon_xloc(recon_t* recon, const parms_t* parms){
	const int npsr=recon->npsr;
	long nin0=0;
	if(parms->load.xloc){
		char* fn=parms->load.xloc;
		warning("Loading xloc from %s\n", fn);
		recon->xloc=loccellread("%s", fn);
		int nxloc=recon->xloc->nx;
		if(nxloc!=npsr)
			error("Invalid saved file. npsr=%d, nxloc=%d\n", npsr, nxloc);
		for(int ips=0; ips<npsr; ips++){
			real dxr=recon->dx->p[ips];
			if(fabs(recon->xloc->p[ips]->dx-dxr)>0.01*dxr){
				warning("xloc[%d]: sampling is %g, expected %g\n", ips, recon->xloc->p[ips]->dx, dxr);
			}
		}
	} else{
		recon->xloc=loccellnew(npsr, 1);
		info("Tomography grid is %ssquare:\n", parms->tomo.square?"":"not ");
		/*FFT in FDPCG prefers power of 2 dimensions. for embeding and fast FFT*/
		if(parms->tomo.nxbase){
			nin0=parms->tomo.nxbase;
		} else if(!parms->sim.idealfit&&(parms->tomo.precond==1||parms->tomo.square==2)){
			/*same square grid dimension in meter on all layers.*/
			long nxmin=LONG_MAX, nymin=LONG_MAX;
			long nxmax=0, nymax=0;
			for(int ips=0; ips<npsr; ips++){
				long nxi, nyi;
				real dxr=recon->dx->p[ips];
				create_metapupil(0, &nxi, &nyi, parms->dirs, parms->aper.d, recon->ht->p[ips], dxr, dxr, 0,
					dxr*parms->tomo.guard, 0, 0, 0, 1);
				nxi/=recon->os->p[ips];
				nyi/=recon->os->p[ips];
				if(nxmax<nxi) nxmax=nxi;
				if(nymax<nyi) nymax=nyi;
				if(nxmin>nxi) nxmin=nxi;
				if(nymin>nyi) nymin=nyi;
			}
			/*FFT grid Must be at least 1.5 times the smallest (on pupil) to avoid
			 * severe aliasing penalty*/
			long nx=MAX(nxmax, nxmin*1.5);
			long ny=MAX(nymax, nymin*1.5);
			nin0=nextfftsize(MAX(nx, ny));

			while(parms->tomo.precond==1&&(nin0&1)){//FFT need even number
				nin0=nextfftsize(nin0+1);
			}
		}
		for(int ips=0; ips<npsr; ips++){
			const real ht=recon->ht->p[ips];
			real dxr=(parms->sim.idealfit)?parms->atm.dx:recon->dx->p[ips];
			const real guard=parms->tomo.guard*dxr;
			long nin=nin0*recon->os->p[ips];
			map_t* map=0;
			create_metapupil(&map, 0, 0, parms->dirs, parms->aper.d, ht, dxr, dxr, 0, guard, nin, nin, 0, parms->tomo.square);
			recon->xloc->p[ips]=map2loc(map, 0);
			loc_create_stat(recon->xloc->p[ips]);
			info("    layer %d: xloc grid is %3ld x %3ld, sampling is %.3f m, %5ld points\n",
				ips, map->nx, map->ny, dxr, recon->xloc->p[ips]->nloc);
			mapfree(map);
		}
	}
	if(parms->gpu.fit==2&&parms->fit.cachex){//to cache x on grid matching floc.
		recon->xcmap=mapcellnew(npsr, 1);
		for(int ips=0; ips<npsr; ips++){
			const real ht=recon->ht->p[ips];
			real dxr=parms->atmr.dx/parms->fit.pos;
			const real guard=parms->tomo.guard*dxr;
			create_metapupil(&recon->xcmap->p[ips], 0, 0, parms->dirs, parms->aper.d, ht, dxr, dxr, 0, guard, 0, 0, 0, parms->fit.square);
			mem_unref(&recon->xcmap->p[ips]->mem);
			recon->xcmap->p[ips]->mem=0;
			recon->xcmap->p[ips]->p=NULL;
		}
	}
	recon->xmap=mapcellnew(npsr, 1);
	recon->xnx=lnew(recon->npsr, 1);
	recon->xny=lnew(recon->npsr, 1);
	recon->xnloc=lnew(recon->npsr, 1);
	for(long i=0; i<recon->npsr; i++){
		recon->xloc->p[i]->iac=parms->tomo.iac;
		loc_create_map_npad(recon->xloc->p[i], (nin0||parms->tomo.square)?0:1,
			nin0*recon->os->p[i], nin0*recon->os->p[i]);
		recon->xmap->p[i]=mapref(recon->xloc->p[i]->map);
		recon->xmap->p[i]->h=recon->ht->p[i];
		recon->xnx->p[i]=recon->xmap->p[i]->nx;
		recon->xny->p[i]=recon->xmap->p[i]->ny;
		recon->xnloc->p[i]=recon->xloc->p[i]->nloc;
	}
	recon->xmcc=dcellnew(npsr, 1);
	for(int ipsr=0; ipsr<npsr; ipsr++){
		recon->xmcc->p[ipsr]=loc_mcc_ptt(recon->xloc->p[ipsr], NULL);
		dinvspd_inplace(recon->xmcc->p[ipsr]);
	}
	if(parms->save.setup){
		writebin(recon->xloc, "xloc");
		writebin(recon->xmap, "xmap");
	}
}
long count_nonzero(const lmat* in){
	long count=0;
	if(in){
		for(long i=0; i<in->nx*in->ny; i++){
			if(in->p[i]){
				count++;
			}
		}
	}
	return count;
}
/**
   Setup the deformable mirrors grid aloc. This is used for DM fitting.
*/
static void
setup_recon_aloc(recon_t* recon, const parms_t* parms){
	const int ndm=parms->ndm;
	if(ndm==0) return;
	if(parms->fit.cachedm){
		recon->acmap=mapcellnew(ndm, 1);
	}
	if(parms->load.aloc){
		char* fn=parms->load.aloc;
		warning("Loading aloc from %s\n", fn);
		recon->aloc=loccellread("%s", fn);
		if(recon->aloc->nx!=ndm||recon->aloc->ny!=1){
			error("Loaded aloc should have %dx1 cells but has %ldx%ld.\n",
				ndm, recon->aloc->nx, recon->aloc->ny);
		}
		for(int idm=0; idm<ndm; idm++){
			if(fabs(parms->dm[idm].dx-recon->aloc->p[idm]->dx)>1e-7){
				error("DM[%d]: loaded aloc has dx=%g while dm.dx=%g\n", idm,
					recon->aloc->p[idm]->dx, parms->dm[idm].dx);
			}
			real max, min;
			dmaxmin(recon->aloc->p[idm]->locx, recon->aloc->p[idm]->nloc, &max, &min);
			if(max-min<parms->aper.d){
				warning("DM[%d]: loaded aloc is too small: diameter is %g while aper.d is %g\n",
					idm, max-min, parms->aper.d);
			}
		}
	} else{
		recon->aloc=loccellnew(ndm, 1);
		/*int nxmax=0, nymax=0; */
		for(int idm=0; idm<ndm; idm++){
			real ht=parms->dm[idm].ht;
			real dx=parms->dm[idm].dx;
			real dy=parms->dm[idm].dy;
			real offset=parms->dm[idm].offset+((int)round(parms->dm[idm].order)%2)*0.5;
			real guard=parms->dm[idm].guard*MAX(dx, dy);
			map_t* map;
			if(parms->dbg.dmfullfov&&!parms->fit.square){//DM covers full fov
				real D=(parms->sim.fov*fabs(ht)+parms->aper.d+guard*2);
				long nx=D/dx+1;
				long ny=D/dy+1;
				map=mapnew(nx, ny, dx, dy);
				map->h=ht;
				map->ox+=offset*dx;
				mapcircle_symbolic(map, D*0.5);
			} else{
				create_metapupil(&map, 0, 0, parms->dirs, parms->aper.d, ht, dx, dy, offset, guard, 0, 0, 0, parms->fit.square);
			}
			info("    DM %d: grid is %ld x %ld\n", idm, map->nx, map->ny);
			recon->aloc->p[idm]=map2loc(map, 0);
			mapfree(map);
		}
	}
	recon->amap=mapcellnew(parms->ndm, 1);
	for(int idm=0; idm<parms->ndm; idm++){
		real ht=parms->dm[idm].ht;
		real offset=parms->dm[idm].offset+((int)round(parms->dm[idm].order)%2)*0.5;
		real dx=parms->dm[idm].dx;
		recon->aloc->p[idm]->iac=parms->dm[idm].iac;
		loc_create_map_npad(recon->aloc->p[idm], parms->fit.square?0:1, 0, 0);
		recon->amap->p[idm]=recon->aloc->p[idm]->map;
		recon->amap->p[idm]->h=ht;
		if(parms->fit.cachedm){
			const real dx2=parms->atmr.dx/parms->fit.pos;
			const real dy2=dx2;
			create_metapupil(&recon->acmap->p[idm], 0, 0, parms->dirs, parms->aper.d,
				ht, dx2, dy2, offset*dx/dx2, dx2, 0, 0, 0, parms->fit.square);
		}
	}

	recon->aimcc=dcellnew(ndm, 1);
	for(int idm=0; idm<ndm; idm++){
		recon->aimcc->p[idm]=loc_mcc_ptt(recon->aloc->p[idm], NULL);
		dinvspd_inplace(recon->aimcc->p[idm]);
	}
	recon->anx=lnew(ndm, 1);
	recon->any=lnew(ndm, 1);
	recon->anloc=lnew(ndm, 1);
	for(int idm=0; idm<ndm; idm++){
		recon->anx->p[idm]=recon->amap->p[idm]->nx;
		recon->any->p[idm]=recon->amap->p[idm]->ny;
		recon->anloc->p[idm]=recon->aloc->p[idm]->nloc;
	}
	/*Dealing with stuck/floating actuators. */
	int anyfloat=0, anystuck=0;
	for(int idm=0; idm<ndm; idm++){
		if(parms->dm[idm].actstuck) anystuck=1;
		if(parms->dm[idm].actfloat) anyfloat=1;
	}
	if(anystuck){
	/**
	   2019-04-17: Stuck actautor is implemented as follows:

	   MVR (recon.alg=0):

	   1) HA is treated so that stuck actuator influence is zeroed. The DM
	   fitting output will therefore set 0 to stuck actuators.

	   2) GA is also treated so that stuck actuators influence is
	   zeroed. This is saying that we treat stuck actautors as input
	   wavefront that the PSOL reconstructor is supposed to control. Without
	   this treatment, the PSOL gradients would not measure those stuck
	   actuators and then cannot use other actuators to reduce their
	   effects.

	   3) For dmreal the stuck actuators are set to their stuck value in
	   clipdm(). The clipped value is fed into the integrator.

	 */
		recon->actstuck=lcellnew(parms->ndm, 1);
		for(int idm=0; idm<ndm; idm++){
			if(!parms->dm[idm].actstuck) continue;
			recon->actstuck->p[idm]=loc_coord2ind(recon->aloc->p[idm], parms->dm[idm].actstuck);
		}
	}
	if(anyfloat){
		recon->actfloat=lcellnew(parms->ndm, 1);
		for(int idm=0; idm<ndm; idm++){
			if(!parms->dm[idm].actfloat) continue;
			recon->actfloat->p[idm]=loc_coord2ind(recon->aloc->p[idm], parms->dm[idm].actfloat);
		}
	}
	if(anystuck||anyfloat){
		for(int idm=0; idm<ndm; idm++){
			int nstuck=recon->actstuck?count_nonzero(recon->actstuck->p[idm]):0;
			int nfloat=recon->actfloat?count_nonzero(recon->actfloat->p[idm]):0;
			info2("DM %d has %d stuck and %d floating actuators\n", idm, nstuck, nfloat);
		}
	}
	if(parms->save.setup){
		writebin(recon->aloc, "aloc");
		writebin(recon->amap, "amap");
	}
}
/**
   Setup ray tracing operator from xloc to ploc for guide stars
*/

static void
setup_recon_HXW(recon_t* recon, const parms_t* parms){
	loc_t* ploc=recon->ploc;
	const int nwfs=parms->nwfsr;
	const int npsr=recon->npsr;
	if(parms->load.HXW){
		warning("Loading saved HXW\n");
		recon->HXW=dspcellread("%s", parms->load.HXW);
		if(recon->HXW->nx!=nwfs||recon->HXW->ny!=npsr){
			error("Wrong saved HXW\n");
		}
		dspcell* HXW=recon->HXW/*PDSPCELL*/;
		int nploc=ploc->nloc;
		for(int ips=0; ips<npsr; ips++){
			int nloc=recon->xloc->p[ips]->nloc;
			for(int iwfs=0; iwfs<nwfs; iwfs++){
				if(!P(HXW, iwfs, ips)
					||P(HXW, iwfs, ips)->nx!=nploc
					||P(HXW, iwfs, ips)->ny!=nloc){
					error("Wrong saved HXW\n");
				}
			}
		}
	} else{
		info("Generating HXW");TIC;tic;
		recon->HXW=dspcellnew(nwfs, npsr);
		dspcell* HXW=recon->HXW/*PDSPCELL*/;
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			int ipowfs=parms->wfsr[iwfs].powfs;

			if(parms->recon.split!=2&&parms->powfs[ipowfs].skip){
			//don't need HXW for low order wfs that does not participate in tomography. 
				continue;
			}
			const real  hs=parms->wfs[iwfs].hs;
			const real  hc=parms->wfs[iwfs].hc;
			loc_t* loc=recon->ploc;
			if(recon->ploc_tel&&recon->ploc_tel->p[iwfs]){
				loc=recon->ploc_tel->p[iwfs];
			}
			for(int ips=0; ips<npsr; ips++){
				const real  ht=recon->ht->p[ips];
				const real  scale=1.-(ht-hc)/hs;
				const real dispx=parms->wfsr[iwfs].thetax*ht;
				const real dispy=parms->wfsr[iwfs].thetay*ht;
				P(HXW, iwfs, ips)=mkh(recon->xloc->p[ips], loc,
					dispx, dispy, scale);
			}
		}
		toc(" ");
	}
	if(parms->save.setup){
		writebin(recon->HXW, "HXW");
	}
	recon->HXWtomo=dspcellnew(recon->HXW->nx, recon->HXW->ny);
	dspcell* HXWtomo=recon->HXWtomo/*PDSPCELL*/;
	dspcell* HXW=recon->HXW/*PDSPCELL*/;
	for(int iwfs=0; iwfs<nwfs; iwfs++){
		int ipowfs=parms->wfsr[iwfs].powfs;
		if(!parms->powfs[ipowfs].skip){/*for tomography */
			for(int ips=0; ips<npsr; ips++){
				P(HXWtomo, iwfs, ips)=dspref(P(HXW, iwfs, ips));
			}
		}
	}
}

/**
   Setup gradient operator from ploc to wavefront sensors.
*/
static void
setup_recon_GP(recon_t* recon, const parms_t* parms, const aper_t* aper){
	loc_t* ploc=recon->ploc;
	const int nwfs=parms->nwfsr;
	recon->GP=dspcellnew(nwfs, 1);
	if(parms->load.GP){
		warning("Loading saved GP\n");
		dspcell* GPload=dspcellread("%s", parms->load.GP);
		int assign=0;
		if(GPload->nx==nwfs){
			assign=1;
		} else if(GPload->nx==parms->npowfs){
			assign=2;
		} else{
			error("GP loaded from %s has wrong size.\n", parms->load.GP);
		}
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			const int ipowfs=parms->wfs[iwfs].powfs;
			const long nsa=recon->saloc->p[ipowfs]->nloc;
			recon->GP->p[iwfs]=dspref(GPload->p[assign==1?iwfs:ipowfs]);
			if(recon->GP->p[iwfs]->nx!=nsa*2||recon->GP->p[iwfs]->ny!=ploc->nloc){
				error("Wrong saved GP[%d]: size is %ldx%ld, need %ldx%ld\n", iwfs,
					recon->GP->p[iwfs]->nx, recon->GP->p[iwfs]->ny, nsa*2, ploc->nloc);
			}
		}
		dspcellfree(GPload);
	} else{
		int share_gp=1;
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			if(parms->recon.misreg_tel2wfs&&parms->recon.misreg_tel2wfs[iwfs]){
				share_gp=0;
			}
		}
		info("Generating GP with ");TIC;tic;
#pragma omp parallel for
		for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
			const int ipowfs=parms->wfsr[iwfs].powfs;
			const int iwfs0=parms->powfs[ipowfs].wfsr->p[0];
			if(parms->powfs[ipowfs].type==1){
				info(" PWFS (skip)");
			} else{
				dmat* gamp=0;
				loc_t* gloc=make_gloc(&gamp, parms, aper, iwfs);
				if(!share_gp||iwfs==iwfs0){
					dsp* gp=0;
					if(parms->powfs[ipowfs].gtype_recon==0){//Average tilt
						info(" Gploc");
						gp=mkg(ploc, gloc, gamp, recon->saloc->p[ipowfs], 1, 0, 0, 1);
					} else if(parms->powfs[ipowfs].gtype_recon==1){//Zernike fit
						info(" Zploc");
						dsp* ZS0=mkz(gloc, gamp->p, recon->saloc->p[ipowfs], 1, 1, 0, 0);
						dsp* H=mkh(ploc, gloc, 0, 0, 1);
						gp=dspmulsp(ZS0, H, "nn");
						dspfree(H);
						dspfree(ZS0);
					} else{
						error("Invalid gtype_recon\n");
					}
					recon->GP->p[iwfs]=gp;
				}
				locfree(gloc);
				dfree(gamp);
			}
		}
		if(share_gp){
			for(int iwfs=0; iwfs<nwfs; iwfs++){
				const int ipowfs=parms->wfsr[iwfs].powfs;
				const int iwfs0=parms->powfs[ipowfs].wfsr->p[0];
				if(iwfs!=iwfs0){
					recon->GP->p[iwfs]=dspref(recon->GP->p[iwfs0]);
				}
			}
		}
		toc(" ");
		if(parms->save.setup){
			writebin(recon->GP, "GP");
		}
	}
}

/**
   Setup gradient operator form aloc for wfs by using GP.
*/
static void
setup_recon_GA(recon_t* recon, const parms_t* parms, const powfs_t* powfs){
	if(parms->nwfs==0) return;
	loc_t* ploc=recon->ploc;
	const int nwfs=parms->nwfsr;
	const int ndm=parms->ndm;
	if(parms->load.GA){
		warning("Loading GA from %s\n", parms->load.GA);
		recon->GA=dspcellread("%s", parms->load.GA);
		if(recon->GA->nx!=nwfs||recon->GA->ny!=ndm)
			error("Wrong saved GA (%ldx%ld). Need (%dx%d)\n",
				recon->GA->nx, recon->GA->ny, nwfs, ndm);
		for(int idm=0; idm<ndm; idm++){
			int nloc=recon->aloc->p[idm]->nloc;
			for(int iwfs=0; iwfs<nwfs; iwfs++){
				int ipowfs=parms->wfsr[iwfs].powfs;
				if(parms->tomo.ahst_idealngs&&parms->powfs[ipowfs].lo){
					continue;
				}
				int nsa=recon->saloc->p[ipowfs]->nloc;
				if(P(recon->GA, iwfs, idm)->nx!=nsa*2||(!parms->recon.modal&&P(recon->GA, iwfs, idm)->ny!=nloc)
					||(parms->recon.modal&&P(recon->GA, iwfs, idm)->ny!=recon->amod->p[idm]->ny)){
					error("Wrong saved GA\n");
				}
			}
		}
	} else{
		info("Generating GA ");TIC;tic;
		recon->GA=dspcellnew(nwfs, ndm);
		if(parms->recon.modal){
			recon->GM=dcellnew(nwfs, ndm);

			recon->amod=dcellnew(ndm, 1);
			recon->anmod=lnew(ndm, 1);
			for(int idm=0; idm<ndm; idm++){
				int nmod=parms->recon.nmod;
				const long nloc=recon->aloc->p[idm]->nloc;
				switch(parms->recon.modal){
				case -2: {//dummy modal control, emulating zonal mode with identity modal matrix
					if(nmod&&nmod!=nloc){
						warning("recon.mod should be 0 or %ld when recon.modal=2 \n", nloc);
					}
					recon->amod->p[idm]=dnew(nloc, nloc);
					real val=sqrt(nloc);
					daddI(recon->amod->p[idm], val);
					dadds(recon->amod->p[idm], -val/nloc);
				}
					   break;
				case -1://zernike
				{
					if(!nmod) nmod=nloc;
					int rmax=floor((sqrt(1+8*nmod)-3)*0.5);
					recon->amod->p[idm]=zernike(recon->aloc->p[idm], 0, 0, rmax, 0);
				}
				break;
				case 1://Karhunen loeve. Don't limit number of modes here to make caching of G_M right.
					recon->amod->p[idm]=KL_vonkarman(recon->aloc->p[idm], 0, parms->atmr.L0);
					break;
				default:
					error("Invalid recon.modal");
				}
				recon->anmod->p[idm]=recon->amod->p[idm]->ny;
			}
		}
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			int ipowfs=parms->wfsr[iwfs].powfs;
			if(parms->tomo.ahst_idealngs&&parms->powfs[ipowfs].lo){
				continue;
			}
			if(parms->powfs[ipowfs].skip==2){//no need for TWFS
				continue;
			}
			const real hs=parms->wfs[iwfs].hs;
			const real hc=parms->wfs[iwfs].hc;
			const loc_t* saloc=recon->saloc->p[ipowfs];
			for(int idm=0; idm<ndm; idm++){
				const real  ht=parms->dm[idm].ht;
				const real  scale=1.-(ht-hc)/hs;
				const loc_t* aloc=recon->aloc->p[idm];
				const real dispx=parms->wfsr[iwfs].thetax*ht;
				const real dispy=parms->wfsr[iwfs].thetay*ht;

				if(parms->powfs[ipowfs].type==1){//PWFS
					if(!parms->powfs[ipowfs].lo){
						dmat* opdadd=0;
						if(!parms->recon.glao&&powfs[ipowfs].opdadd&&0){
							int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
							opdadd=powfs[ipowfs].opdadd->p[wfsind];
						}
						if(!parms->recon.modal){
							info("\nPyWFS from aloc to saloc directly\n");
							dmat* tmp=pywfs_mkg(powfs[ipowfs].pywfs, aloc, parms->recon.misreg_dm2wfs[iwfs+idm*nwfs],
								0, opdadd, dispx, dispy);
							P(recon->GA, iwfs, idm)=d2sp(tmp, dmaxabs(tmp)*1e-6);
							dfree(tmp);
						} else{
							info("\nPyWFS from amod to saloc directly\n");
							//We compute the GM for full set of modes so that it is cached only once.
							P(recon->GM, iwfs, idm)=pywfs_mkg(powfs[ipowfs].pywfs, aloc, parms->recon.misreg_dm2wfs[iwfs+idm*nwfs],
								recon->amod->p[idm], opdadd, dispx, dispy);
						}
					}
				} else{//SHWFS
					char* input=parms->misreg.dm2wfs?parms->misreg.dm2wfs[iwfs+idm*nwfs]:0;
					char* calib=parms->recon.misreg_dm2wfs?parms->recon.misreg_dm2wfs[iwfs+idm*nwfs]:0;
					loc_t* loc=ploc;

					if(input&&saloc->nloc>4){//there is distortion input			
						if(calib){//there is magical calibration input
							loc=loctransform(ploc, calib);
						} else{//determine distortion by fitting
							warning_once("dm2wfs: determine distortion by fitting \n");
							//First, simulate the poke matrix measurement process.
							loc_t* loc2=loctransform(ploc, input);
							dsp* H2=mkh(aloc, loc2, dispx, dispy, scale);
							locfree(loc2);
							dsp* GA2=dspmulsp(recon->GP->p[iwfs], H2, "nn"); //measured GA.
							dspfree(H2);
							dmat* calib2=loc_calib(GA2, aloc, saloc, dispx, dispy, scale, 2);
							dspfree(GA2);
							loc=loctransform2(ploc, calib2);
							dfree(calib2);
						}
					}
					dsp* H=mkh(aloc, loc, dispx, dispy, scale);
					P(recon->GA, iwfs, idm)=dspmulsp(recon->GP->p[iwfs], H, "nn");
					dspfree(H);
					if(loc!=ploc){
						locfree(loc);
					}
					if(parms->recon.modal){
						dspmm(PP(recon->GM, iwfs, idm), P(recon->GA, iwfs, idm), recon->amod->p[idm], "nn", 1);
					}
				}
			}/*idm */
		}//iwfs
		toc(" ");
	}
	if(parms->recon.modal&&parms->recon.nmod>0){
		for(int idm=0; idm<ndm; idm++){
			if(parms->recon.nmod<recon->amod->p[idm]->ny){
				info("DM %d:Reduce number of controlled modes from %ld to %d\n",
					idm, recon->amod->p[idm]->ny, parms->recon.nmod);
				dresize(recon->amod->p[idm], 0, parms->recon.nmod);
				recon->anmod->p[idm]=recon->amod->p[idm]->ny;
				for(int iwfs=0; iwfs<nwfs; iwfs++){
					if(!P(recon->GM, iwfs, idm)) continue;
					dresize(P(recon->GM, iwfs, idm), 0, parms->recon.nmod);
				}
			}
		}
	}
	if(!parms->recon.modal){
		if(recon->actstuck&&parms->recon.alg==1&&parms->dbg.recon_stuck){
		/*This is need for LSR reconstructor to skip stuck actuators.  GA is
		  also used to form PSOL gradients, but that one doesn't need this
		  modification because actuator extropolation was already applied.*/
		/*
		  For MVR we want to treat stuck actuators as input OPD. Therefore its
		  effect on GA is also removed
		 */
			warning("Apply stuck actuators to GA\n");
			act_stuck(recon->aloc, recon->GA, recon->actstuck);

		}
		if(parms->recon.alg==1){//LSR.
			recon->actcpl=genactcpl(recon->GA, 0);
			act_stuck(recon->aloc, recon->actcpl, recon->actfloat);
			if(parms->lsr.actinterp){
				recon->actinterp=act_extrap(recon->aloc, recon->actcpl, parms->lsr.actthres);
			} else if(recon->actfloat){
				warning("There are float actuators, but fit.actinterp is off\n");
			}
			if(recon->actinterp){
				dspcell* GA2=0;
				dcellmm(&GA2, recon->GA, recon->actinterp, "nn", 1);
				dspcellfree(recon->GA);
				recon->GA=GA2;
			}
			if(parms->save.setup){
				if(recon->actinterp){
					writebin(recon->actinterp, "lsr_actinterp");
				}
				if(recon->actcpl){
					writebin(recon->actcpl, "lsr_actcpl");
				}
			}
		}
	}
	/*Create GAlo that only contains GA for low order wfs */
	recon->GAlo=cellnew(recon->GA->nx, recon->GA->ny);
	recon->GAhi=dspcellnew(recon->GA->nx, recon->GA->ny);
	if(parms->recon.modal) recon->GMhi=dcellnew(nwfs, ndm);

	for(int idm=0; idm<ndm; idm++){
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			int ipowfs=parms->wfsr[iwfs].powfs;
			if(parms->powfs[ipowfs].lo
				||(parms->recon.split&&parms->nlopowfs==0&&!parms->powfs[ipowfs].trs)){/*for low order wfs */
				if(parms->recon.modal){
					P(recon->GAlo, iwfs, idm)=(cell*)dref(P(recon->GM, iwfs, idm));
				}else{
					P(recon->GAlo, iwfs, idm)=(cell*)dspref(P(recon->GA, iwfs, idm));
				}
			}
			if(!parms->powfs[ipowfs].skip){
				P(recon->GAhi, iwfs, idm)=dspref(P(recon->GA, iwfs, idm));
				if(parms->recon.modal){
					P(recon->GMhi, iwfs, idm)=dref(P(recon->GM, iwfs, idm));
				}
			}
		}
	}
	if(parms->save.setup){
		writebin(recon->GA, "GA");
		if(parms->recon.modal){
			writebin(recon->amod, "amod");
			writebin(recon->GM, "GM");
		}
	}
}
/**
   Crate the xloc to wfs gradient operator.
*/
static void
setup_recon_GX(recon_t* recon, const parms_t* parms){
	const int nwfs=parms->nwfsr;
	const int npsr=recon->npsr;
	recon->GX=dspcellnew(nwfs, npsr);
	dspcell* GX=recon->GX/*PDSPCELL*/;
	dspcell* HXW=recon->HXW/*PDSPCELL*/;
	info("Generating GX ");TIC;tic;
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	/*gradient from xloc. Also useful for lo WFS in MVST mode. */
		for(int ips=0; ips<npsr; ips++){
			P(GX, iwfs, ips)=dspmulsp(recon->GP->p[iwfs], P(HXW, iwfs, ips), "nn");
		}/*ips */
	}
	toc(" ");
	recon->GXtomo=dspcellnew(recon->GX->nx, recon->GX->ny);
	dspcell* GXtomo=recon->GXtomo/*PDSPCELL*/;

	recon->GXlo=dspcellnew(recon->GX->nx, recon->GX->ny);
	dspcell* GXlo=recon->GXlo/*PDSPCELL*/;

	int nlo=parms->nlopowfs;
	for(int iwfs=0; iwfs<nwfs; iwfs++){
		int ipowfs=parms->wfsr[iwfs].powfs;
		for(int ips=0; ips<npsr; ips++){
			if(!parms->powfs[ipowfs].skip){/*for tomography */
				P(GXtomo, iwfs, ips)=dspref(P(GX, iwfs, ips));
			}
			if(parms->powfs[ipowfs].lo
				||(parms->recon.split&&nlo==0&&!parms->powfs[ipowfs].trs)){
			 /*for low order wfs or extracted t/t for high order ngs wfs.*/
				P(GXlo, iwfs, ips)=dspref(P(GX, iwfs, ips));
			}
		}
	}/*iwfs */
}

/**
   From focus mode to gradients. This acts on WFS, not WFSR.
 */
static void
setup_recon_GF(recon_t* recon, const parms_t* parms){
	/*Create GFall: Focus mode -> WFS grad. This is model*/
	recon->GFall=dcellnew(parms->nwfs, 1);
	recon->GFngs=dcellnew(parms->nwfs, 1);
	{
		dmat* opd=dnew(recon->ploc->nloc, 1);
		loc_add_focus(opd, recon->ploc, 1);
		for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
			const int ipowfs=parms->wfs[iwfs].powfs;
			dspmm(&recon->GFall->p[iwfs], recon->GP->p[parms->recon.glao?ipowfs:iwfs], opd, "nn", 1);
			if(parms->powfs[ipowfs].lo&&!parms->powfs[ipowfs].llt){
				recon->GFngs->p[iwfs]=dref(recon->GFall->p[iwfs]);
			}
		}
		dfree(opd);
	}
	if(parms->save.setup){
		writebin(recon->GFall, "GFall");
	}
}
/**
   From radial order modes to gradients.
 */
static void
setup_recon_GR(recon_t* recon, const powfs_t* powfs, const parms_t* parms){
	recon->GRall=dcellnew(parms->nwfs, 1);
	dmat* opd=0;
	real reduce=recon->ploc->dx*2;//to reduce the edge effect.
	int zmax=parms->dbg.twfsrmax?parms->dbg.twfsrmax:(parms->powfs[parms->itpowfs].order/2);
	int zmin=parms->dbg.twfsflag>2?2:3;
	int zradonly=parms->dbg.twfsflag==1||parms->dbg.twfsflag==3;
	dbg("twfs mode is %s from order %d to %d.\n", zradonly?"radial":"all modes", zmin, zmax);
	opd=zernike(recon->ploc, parms->aper.d-reduce, zmin, zmax, zradonly);

	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		const int ipowfs=parms->wfs[iwfs].powfs;
		if(parms->powfs[ipowfs].skip==2||parms->powfs[ipowfs].llt){
			if(parms->powfs[ipowfs].type==1){//PWFS
				recon->GRall->p[iwfs]=pywfs_mkg(powfs[ipowfs].pywfs, recon->ploc, NULL, opd, 0, 0, 0);
			} else{//SHWFS
				dspmm(&recon->GRall->p[iwfs], recon->GP->p[parms->recon.glao?ipowfs:iwfs], opd, "nn", 1);
			}
		}
	}
	if(parms->save.setup){
		writebin(recon->GRall, "twfs_GR");
		writebin(opd, "twfs_opd");
	}
	dfree(opd);
}
/**
   Tilt removal from DM command. Used by filter.c
 */
void setup_recon_dmttr(recon_t* recon, const parms_t* parms){
	recon->DMTT=dcellnew(parms->ndm, 1);
	recon->DMPTT=dcellnew(parms->ndm, 1);
	/*if(!recon->actcpl && parms->nwfs>0){
	error("actcpl must not be null\n");
	}*/
	for(int idm=0; idm<parms->ndm; idm++){
		recon->DMTT->p[idm]=loc2mat(recon->aloc->p[idm], 0);
	}
	if(parms->dbg.recon_stuck){
		act_zero(recon->aloc, recon->DMTT, recon->actstuck);
	}
	for(int idm=0; idm<parms->ndm; idm++){
		recon->DMPTT->p[idm]=dpinv(recon->DMTT->p[idm], 0);
	}
	if(parms->save.setup){
		writebin(recon->DMTT, "DMTT");
		writebin(recon->DMPTT, "DMPTT");
	}
}
/**
   setting up global tip/tilt remove operator from WFS gradients.
*/

static void
setup_recon_TT(recon_t* recon, const parms_t* parms, const powfs_t* powfs){
	int nwfs=parms->nwfsr;
	recon->TT=dcellnew(nwfs, nwfs);
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(parms->powfs[ipowfs].nwfs==0) continue;
		if(parms->powfs[ipowfs].trs
			||(parms->recon.split&&!parms->powfs[ipowfs].lo&&!parms->powfs[ipowfs].skip)
			||parms->powfs[ipowfs].dither==1){
			int nsa=recon->saloc->p[ipowfs]->nloc;
			dmat* TT=0;
			if(parms->powfs[ipowfs].type==0){//SHWFS
				TT=dnew(nsa*2, 2);
				real* TTx=TT->p;
				real* TTy=TT->p+nsa*2;
				for(int isa=0; isa<nsa; isa++){
					TTx[isa]=1;
					TTy[isa]=0;
				}
				for(int isa=nsa; isa<nsa*2; isa++){
					TTx[isa]=0;
					TTy[isa]=1;
				}
			} else if(parms->powfs[ipowfs].type==1){//PYWFS
				TT=dref(powfs[ipowfs].pywfs->GTT); //pywfs_tt(powfs[ipowfs].pywfs);
			} else{
				error("Invalid powfs.type\n");
			}
			if(parms->recon.glao){
				recon->TT->p[ipowfs*(parms->npowfs+1)]=ddup(TT);
			} else{
				for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
					int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
					recon->TT->p[iwfs*(1+parms->nwfsr)]=ddup(TT);
				}
			}
			dfree(TT);
		}
	}
	if(parms->save.setup){
		writebin(recon->TT, "TT");
	}
}
/**
   operator to remove diffrential focus modes that might be caused by sodium layer
   horizontal structure.
*/

static void
setup_recon_DF(recon_t* recon, const parms_t* parms){
	if(!recon->has_dfr) return;

	int nwfs=parms->nwfsr;
	recon->DF=dcellnew(nwfs, nwfs);
	/*Then differential focus modes. */
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(parms->powfs[ipowfs].nwfs==0) continue;
		if(parms->powfs[ipowfs].dfrs){
			if(parms->powfs[ipowfs].nwfs<2){
				error("This powfs group has only 1 wfs. Could not remove diff focus\n");
			}
			int nsa=recon->saloc->p[ipowfs]->nloc;
			dmat* DF=dnew(nsa*2, 1);
			/*saloc is lower left corner of subaperture. don't have to be the center. */
			memcpy(DF->p, recon->saloc->p[ipowfs]->locx, sizeof(real)*nsa);
			memcpy(DF->p+nsa, recon->saloc->p[ipowfs]->locy, sizeof(real)*nsa);
			/**
			   postive focus on first wfs. negative focus on diagnonal wfs.
			*/
			for(int jwfs=1; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
				int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
				if(parms->powfs[ipowfs].skip){
					error("This WFS %d should be included in Tomo.\n", iwfs);
				}
				dcp(&recon->DF->p[iwfs*nwfs], DF);
				dadd(&recon->DF->p[iwfs+iwfs*nwfs], 0, DF, -1);
			}
			dfree(DF);
		}
	}
	if(parms->save.setup){
		writebin(recon->DF, "DF");
	}
}
/**
   Dither using command path (DM) aberration
 */
void setup_recon_dither_dm(recon_t* recon, const powfs_t* powfs, const parms_t* parms){
	int any=0;
	int dither_mode=0;
	real dither_amp=0;
	int dither_npoint=0;
	int dither_dtrat=0;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(parms->powfs[ipowfs].dither>1){//common path dithering
			if(any){//already found, check consistency
				if(dither_mode!=-parms->powfs[ipowfs].dither||
					fabs(dither_amp-parms->powfs[ipowfs].dither_amp)>dither_amp*1e-5
					||dither_npoint!=parms->powfs[ipowfs].dither_npoint
					||dither_dtrat!=parms->powfs[ipowfs].dtrat){
					error("Multiple dither with different configuration is not supported\n");
				}
			}
			dither_mode=-parms->powfs[ipowfs].dither;
			dither_amp=parms->powfs[ipowfs].dither_amp;
			dither_npoint=parms->powfs[ipowfs].dither_npoint;
			dither_dtrat=parms->powfs[ipowfs].dtrat;
			any=1;
		}
	}
	if(any){
		const int idm=parms->idmground;
		recon->dither_npoint=dither_npoint;
		recon->dither_dtrat=dither_dtrat;
		dcellfree(recon->dither_m);
		dcellfree(recon->dither_ra);
		dcellfree(recon->dither_rg);
		recon->dither_m=dcellnew(parms->ndm, 1);
		recon->dither_m->p[idm]=zernike(recon->aloc->p[idm], parms->aper.d, 0, 0, dither_mode);
		dscale(recon->dither_m->p[idm], dither_amp);
		recon->dither_ra=dcellnew(parms->ndm, parms->ndm);
		P(recon->dither_ra, idm, idm)=dpinv(recon->dither_m->p[idm], 0);
		recon->dither_rg=dcellnew(parms->nwfsr, parms->nwfsr);
		for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
			int ipowfs=parms->wfsr[iwfs].powfs;
			const real hc=parms->wfsr[iwfs].hc;
			if(parms->powfs[ipowfs].dither>1){
				dmat* opd=dnew(powfs[ipowfs].loc->nloc, 1);
				real ht=parms->dm[idm].ht+parms->dm[idm].vmisreg;
				real scale=1.-(ht-hc)/parms->powfs[ipowfs].hs;
				real dispx=ht*parms->wfsr[iwfs].thetax;
				real dispy=ht*parms->wfsr[iwfs].thetay;
				prop_nongrid(recon->aloc->p[idm], recon->dither_m->p[idm]->p,
					powfs[ipowfs].loc, opd->p,
					-1, dispx, dispy, scale, 0, 0);
				dmat* ints=0;
				dmat* grad=0;
				pywfs_fft(&ints, powfs[ipowfs].pywfs, opd);
				pywfs_grad(&grad, powfs[ipowfs].pywfs, ints);
				P(recon->dither_rg, iwfs, iwfs)=dpinv(grad, P(recon->saneai, iwfs, iwfs));
				if(0){//test linearity
					dscale(opd, 1./4.);
					dmat* tmp=0;
					dmat* res=dnew(10, 1);
					for(int i=0; i<10; i++){
						dscale(opd, 2);
						dzero(ints);
						pywfs_fft(&ints, powfs[ipowfs].pywfs, opd);
						pywfs_grad(&grad, powfs[ipowfs].pywfs, ints);
						dmm(&tmp, 0, P(recon->dither_rg, iwfs, iwfs), grad, "nn", 1);
						res->p[i]=tmp->p[0];
					}
					writebin(res, "linearity");
					dfree(tmp);
					dfree(res);
					exit(0);
				}
				dfree(ints);
				dfree(opd);
				dfree(grad);
			}
		}
		if(parms->save.setup){
			writebin(recon->dither_m, "dither_m");
			writebin(recon->dither_ra, "dither_ra");
			writebin(recon->dither_rg, "dither_rg");
		}
	}
}

/**
   Create reconstruction parameters that are related to the geometry only, and
   will not be updated when estimated WFS measurement noise changes.

   This can be used to do NCPA calibration.
 */
recon_t* setup_recon_prep(const parms_t* parms, const aper_t* aper, const powfs_t* powfs){
	info("\n%sSetting up reconstructor geometry.%s\n\n", GREEN, BLACK);
	recon_t* recon=mycalloc(1, recon_t);
	if(parms->recon.warm_restart){
		info("Wavefront reconstruction uses warm restart.\n");
	} else{
		warning("Wavefront reconstruction does not use warm restart.\n");
	}
	if(parms->cn2.pair&&parms->cn2.pair->nx>0&&!recon->cn2est){
	/*setup CN2 Estimator. It determines the reconstructed layer heigh can be fed to the tomography */
		recon->cn2est=cn2est_prepare(parms, powfs);
	}
	recon->saloc=loccellnew(parms->npowfs, 1);
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		P(recon->saloc, ipowfs)=locref(powfs[ipowfs].saloc);
	}
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(parms->powfs[ipowfs].nwfs==0) continue;
		/*we will remove tip/tilt from the high order NGS wfs in split tomo mode.*/
		if(parms->powfs[ipowfs].trs||(parms->recon.split&&!parms->powfs[ipowfs].lo)){
			recon->has_ttr=1;
			break;
		}
	}
	if(!parms->recon.glao){
		for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
			if(parms->powfs[ipowfs].nwfs<=1) continue;
			if(parms->powfs[ipowfs].dfrs){
				recon->has_dfr=1;
				break;
			}
		}
	}
	recon->ngrad=lnew(parms->nwfsr, 1);
	for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
		const int ipowfs=parms->wfsr[iwfs].powfs;
		if(!parms->powfs[ipowfs].skip){
			recon->ngrad->p[iwfs]=recon->saloc->p[ipowfs]->nloc*2;
		}
	}
	/*to be used in tomography. */
	recon->nthread=NTHREAD;
	/*setup pupil coarse grid for gradient operator*/
	setup_recon_ploc(recon, parms);
	/*setup DM actuator grid */
	setup_recon_aloc(recon, parms);
	/*Grid for DM fitting*/
	setup_recon_floc(recon, parms);
	/*Gradient operators*/
	setup_recon_GP(recon, parms, aper);
	//TT Removal
	setup_recon_TT(recon, parms, powfs);
	//DF removal. //Deprecated
	if(recon->has_dfr){
		setup_recon_DF(recon, parms);
	}
	if(recon->DF){
		recon->TTF=dcellcat(recon->TT, recon->DF, 2);
	} else{
		recon->TTF=dcellref(recon->TT);
	}
	if(parms->recon.alg==0&&!parms->sim.idealfit){//tomography parameters
		if(parms->cn2.tomo&&recon->cn2est){
			/*Use cn2 estimation results for tomography. Use its ht to build
			  reconstructor matrices.*/
			cn2est_t* cn2est=recon->cn2est;
			recon->ht=dref(cn2est->htrecon);
			recon->os=dref(cn2est->os);
			recon->wt=dref(cn2est->wtrecon->p[0]);
			/*the following will be updated later in simulation. */
			if(parms->cn2.keepht){
				for(int ips=0; ips<recon->wt->nx; ips++){
					recon->wt->p[ips]=parms->atmr.wt->p[ips];
				}
			} else{
				dset(recon->wt, 1./recon->wt->nx);/*evenly distributed.  */
			}
		} else{/*use input information from atmr */
			recon->wt=dnew(parms->atmr.nps, 1);
			recon->ht=dnew(parms->atmr.nps, 1);
			recon->os=dnew(parms->atmr.nps, 1);
			for(int ips=0; ips<recon->wt->nx; ips++){
				recon->wt->p[ips]=parms->atmr.wt->p[ips];
				recon->ht->p[ips]=parms->atmr.ht->p[ips];
				recon->os->p[ips]=parms->atmr.os->p[ips];
			}
		}
		recon->r0=parms->atmr.r0;
		recon->L0=parms->atmr.L0;

		/*sampling of xloc */
		recon->dx=dnew(recon->ht->nx, 1);
		for(int iht=0; iht<recon->ht->nx; iht++){
			real scale=1.0-recon->ht->p[iht]/parms->atmr.hs;
			recon->dx->p[iht]=(parms->atmr.dx/recon->os->p[iht])*scale;
		}
		/*number of reconstruction layers */
		recon->npsr=recon->ht->nx;
		/*setup atm reconstruction layer grid */
		setup_recon_xloc(recon, parms);
		/*setup xloc/aloc to WFS grad */
		setup_recon_HXW(recon, parms);
		setup_recon_GX(recon, parms);//uses HXW.
		dspcellfree(recon->HXW);/*only keep HXWtomo for tomography */
	}

	return recon;
}
/**
   That may depend on GPU data.
 */
void setup_recon_prep2(recon_t* recon, const parms_t* parms, const aper_t* aper, const powfs_t* powfs){
	info2("\n%sSetting up reconstructor%s\n\n", GREEN, BLACK);
	setup_recon_GA(recon, parms, powfs);//PWFS uses GPU data.
	setup_recon_GF(recon, parms);//GF depends on GA.
	if(parms->itpowfs!=-1){
		setup_recon_GR(recon, powfs, parms);
	}
	if(parms->recon.split){
		setup_ngsmod_prep(parms, recon, aper, powfs);
	}
	setup_recon_dmttr(recon, parms);
}

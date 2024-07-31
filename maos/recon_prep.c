/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>

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
#include "recon.h"
#include "recon_utils.h"
#include "powfs.h"
#include "ahst.h"
#undef GREEN
#define GREEN BLACK
/**
   \file recon_prep.c

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
	int square=parms->tomo.square;
	if(parms->load.ploc){/*optionally load ploc from the file. see dbg.conf */
		warning("Loading ploc from %s\n", parms->load.ploc);
		recon->ploc=locread("%s", parms->load.ploc);
		if(fabs(recon->ploc->dx-dxr)>dxr*1e-6){
			warning("Loaded ploc has unexpected sampling of %g, should be %g\n",
				recon->ploc->dx, dxr);
		}
	} else{
		map_t* pmap=0;
		real rmax=0;//additional guard band to accomodate misregistration.
		for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
			if(rmax<parms->powfs[ipowfs].misregrmax) rmax=parms->powfs[ipowfs].misregrmax;
		}
		real guard=MAX(1, parms->tomo.guard)*dxr;
		if(rmax) guard=MAX(guard,rmax); 
		//an additional ring of guard points is needed if square is used.
		create_metapupil(&pmap, 0, 0, parms->dirs, parms->aper.d, 0, dxr, dxr, 0, guard, 0, 0, 0, square);
		info("PLOC is %ldx%ld, with sampling of %.2fm (%ssquare)\n", NX(pmap), NY(pmap), dxr, square?"":"not ");
		recon->ploc=map2loc(pmap, 0);/*convert map_t to loc_t */
		mapfree(pmap);
	}
	if(parms->save.setup){
		locwrite(recon->ploc, "ploc");
	}
	loc_create_map_npad(recon->ploc, square?0:1, 0, 0);
	recon->pmap=recon->ploc->map;
	loc_create_stat(recon->ploc);
	if(parms->recon.distortion_tel2wfs){
		//determine the misregistration in shift and rotation from distortion_tel2wfs.
		real ploc_xm=0, ploc_ym=0;
		locmean(&ploc_xm, &ploc_ym, recon->ploc);
		for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
			if(parms->recon.distortion_tel2wfs[iwfs]){
				loc_t *ploc2=loctransform(recon->ploc, parms->recon.distortion_tel2wfs[iwfs]);
				real xm=0, ym=0;
				locmean(&xm, &ym, ploc2);
				locshift(ploc2, -xm, -ym);
				parms->wfsr[iwfs].misregx+=xm-ploc_xm;
				parms->wfsr[iwfs].misregy+=ym-ploc_ym;
				parms->wfsr[iwfs].misregc+=loc_angle(recon->ploc, ploc2);
				dbg("ploc for wfs%d shifted by (%g, %g), and rotated by %g CCW.\n", iwfs,
					parms->wfsr[iwfs].misregx, parms->wfsr[iwfs].misregy, parms->wfsr[iwfs].misregc);
				locfree(ploc2);
			}
		}
	}
}
/**
   Create loc/amp that can be used to build GP. It has points on edge of
   subapertures which is different from powfs.lo. The amplitude depends on the
   distortion as we assume that the distortion can be extracted from dm to wfs
   interaction matrix calibration.

   gloc and the rotated version are centered around misregx and misregy so mkamp does not use misreg.
   the returned gloc is defined at the WFS reference frame
*/

void shwfs_gloc_gamp(loc_t** pgloc, dmat** gamp, const parms_t* parms, const aper_t* aper, int iwfsr){
	const int ipowfs=parms->wfsr[iwfsr].powfs;
	real dx=parms->powfs[ipowfs].dx;
	if(!*pgloc){
		real rmax=parms->powfs[ipowfs].misregrmax;
		*pgloc=mkannloc(parms->aper.d+rmax*2+2*dx, 0, dx, 0);
	}
	
	loc_t *gloc=*pgloc;
	if(gamp){
		if((parms->distortion.tel2wfs&&parms->distortion.tel2wfs[iwfsr])||parms->wfsr[iwfsr].misregc){
			if(parms->distortion.tel2wfs&&parms->distortion.tel2wfs[iwfsr]){
				gloc=loctransform(gloc, parms->distortion.tel2wfs[iwfsr]);
			} else{//misregc is set. need to rotate grid, so duplicate it.
				gloc=locdup(gloc);
			}
			if(parms->wfsr[iwfsr].misregc){
				locrot(gloc, -parms->wfsr[iwfsr].misregc);
			}
		}
		//now, gloc is transformed to telescoe pupil reference frame.
		real dout=parms->aper.d;
		real din=parms->aper.din;
		map_t* ampground=parms->dbg.gp_noamp>0?0:aper->ampground;
		if(parms->dbg.gp_noamp==2){//make telescope aperture larger.
			dout+=1;
			din=0;
		}
		if(*gamp) dfree(*gamp);
		*gamp=mkamp(gloc, ampground, parms->wfsr[iwfsr].misregx-P(parms->aper.misreg, 0), 
			parms->wfsr[iwfsr].misregy-P(parms->aper.misreg, 1), dout, din);
	
		if(gloc!=*pgloc){
			locfree(gloc);//we do not return rotated or distorted gloc. 
		}
	}
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
		real guard=MAX(1,parms->tomo.guard)*dxr;
		map_t* fmap=0;
		create_metapupil(&fmap, 0, 0, parms->dirs, parms->aper.d, 0, dxr, dxr, 0, guard, 0, 0, 0, parms->fit.square);
		info("FLOC is %ldx%ld, with sampling of %.2fm (%ssquare)\n", NX(fmap), NY(fmap), dxr, parms->fit.square?"":"not ");
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
		real rout=parms->aper.d/2;
		if(parms->dbg.annular_W&&parms->aper.din>0){
			dbg("Define the W0/W1 on annular aperture instead of circular.\n");
			rin=parms->aper.din/2;
		}
		{
			rout=loc_diam(recon->floc)/2;
			dbg("rout is set to floc radius %g\n", rout);
		}
		mkw_annular(recon->floc, 0, 0, rin, rout, &(recon->W0), &(recon->W1));
	}
	if(parms->save.recon){
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
	if(parms->sim.idealtomo) return;
	const int npsr=recon->npsr;
	long nin0=0;
	if(parms->load.xloc){
		char* fn=parms->load.xloc;
		warning("Loading xloc from %s\n", fn);
		recon->xloc=loccellread("%s", fn);
		int nxloc=NX(recon->xloc);
		if(nxloc!=npsr)
			error("Invalid saved file. npsr=%d, nxloc=%d\n", npsr, nxloc);
		for(int ips=0; ips<npsr; ips++){
			real dxr=P(recon->dx, ips);
			if(fabs(P(recon->xloc, ips)->dx-dxr)>0.01*dxr){
				warning("xloc[%d]: sampling is %g, expected %g\n", ips, P(recon->xloc, ips)->dx, dxr);
			}
		}
	} else{
		recon->xloc=loccellnew(npsr, 1);
		info("Tomography grid xloc is %ssquare:\n", parms->tomo.square?"":"not ");
		const real guard=P(recon->dx, 0)*MAX(parms->tomo.guard, 1);//use ground layer info as guard
		/*FFT in FDPCG prefers power of 2 dimensions. for embeding and fast FFT*/
		if(parms->tomo.nxbase){
			nin0=parms->tomo.nxbase;
		} else if(parms->tomo.precond==1||parms->tomo.square==2){
			/*same square grid dimension in meter on all layers.*/
			long nxmin=LONG_MAX, nymin=LONG_MAX;
			long nxmax=0, nymax=0;
			for(int ips=0; ips<npsr; ips++){
				long nxi, nyi;
				const real ht=P(recon->ht, ips);
				const real dxr=P(recon->dx, ips);
				const real offset=ips==0?0:1./(ips+1);
				create_metapupil(0, &nxi, &nyi, parms->dirs, parms->aper.d, ht, dxr, dxr, offset, guard, 0, 0, 0, 1);
				nxi/=P(recon->os, ips);
				nyi/=P(recon->os, ips);
				if(nxmax<nxi) nxmax=nxi;
				if(nymax<nyi) nymax=nyi;
				if(nxmin>nxi) nxmin=nxi;
				if(nymin>nyi) nymin=nyi;
			}
			if(parms->tomo.precond==1){
				/*FFT grid Must be at least 1.5 times the smallest (on pupil) to avoid
			 		* severe aliasing penalty*/
				long nx=MAX(nxmax, nxmin*1.5);
				long ny=MAX(nymax, nymin*1.5);
				long nn=MAX(nx,ny);
				if(nn&1){
					nn=nn+1;
				}
				nin0=nextfftsize(nn);
			}else{
				nin0=MAX(nxmax, nymax);
			}
		}

		for(int ips=0; ips<npsr; ips++){
			const real ht=P(recon->ht, ips);
			const real dxr=P(recon->dx, ips);
			const long nin=nin0*P(recon->os, ips);
			map_t* map=0;
			const real offset=ips==0?0:1./(ips+1);
			create_metapupil(&map, 0, 0, parms->dirs, parms->aper.d, ht, dxr, dxr, offset, guard, nin, nin, 0, parms->tomo.square);
			P(recon->xloc, ips)=map2loc(map, 0);
			loc_create_stat(P(recon->xloc, ips));
			info("    layer %d: grid is %4ldx%4ld @ %.3f m, offset is %.2f, %5ld points\n",
				ips, NX(map), NY(map), dxr, offset, P(recon->xloc, ips)->nloc);
			mapfree(map);
		}
	}
	if(parms->gpu.fit==2&&parms->fit.cachex){//to cache x on grid matching floc.
		recon->xcmap=mapcellnew(npsr, 1);
		for(int ips=0; ips<npsr; ips++){
			const real ht=P(recon->ht, ips);
			real dxr=parms->atmr.dx/parms->fit.pos;
			const real guard=MAX(1, parms->tomo.guard)*dxr;
			create_metapupil(&P(recon->xcmap, ips), 0, 0, parms->dirs, parms->aper.d, ht, dxr, dxr, 0, guard, 0, 0, 0, parms->fit.square);
			mem_unref(&P(recon->xcmap, ips)->mem);
			P(recon->xcmap, ips)->mem=0;
			P(P(recon->xcmap, ips))=NULL;
		}
	}
	recon->xmap=mapcellnew(npsr, 1);
	recon->xnx=lnew(recon->npsr, 1);
	recon->xny=lnew(recon->npsr, 1);
	recon->xnloc=lnew(recon->npsr, 1);
	for(long i=0; i<recon->npsr; i++){
		P(recon->xloc, i)->iac=parms->tomo.iac;
		loc_create_map_npad(P(recon->xloc, i), (nin0||parms->tomo.square)?0:1,
			nin0*P(recon->os, i), nin0*P(recon->os, i));
		P(recon->xmap, i)=mapref(P(recon->xloc, i)->map);
		P(recon->xmap, i)->h=P(recon->ht, i);
		P(recon->xnx, i)=P(recon->xmap, i)->nx;
		P(recon->xny, i)=P(recon->xmap, i)->ny;
		P(recon->xnloc, i)=P(recon->xloc, i)->nloc;
	}
	recon->xmcc=dcellnew(npsr, 1);
	for(int ipsr=0; ipsr<npsr; ipsr++){
		P(recon->xmcc, ipsr)=loc_mcc_ptt(P(recon->xloc, ipsr), NULL);
		dinvspd_inplace(P(recon->xmcc, ipsr));
	}
	if(parms->save.setup){
		writebin(recon->xloc, "xloc");
	}
}
long count_nonzero(const lmat* in){
	long count=0;
	if(in){
		for(long i=0; i<NX(in)*NY(in); i++){
			if(P(in, i)){
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
		if(NX(recon->aloc)!=ndm||NY(recon->aloc)!=1){
			error("Loaded aloc should have %dx1 cells but has %ldx%ld.\n",
				ndm, NX(recon->aloc), NY(recon->aloc));
		}
		for(int idm=0; idm<ndm; idm++){
			if(fabs(parms->dm[idm].dx-P(recon->aloc, idm)->dx)>1e-7){
				error("DM[%d]: loaded aloc has dx=%g while dm.dx=%g\n", idm,
					P(recon->aloc, idm)->dx, parms->dm[idm].dx);
			}
			real max, min;
			dvecmaxmin(P(recon->aloc, idm)->locx, P(recon->aloc, idm)->nloc, &max, &min);
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
			real offset=wrap2range(parms->dm[idm].offset,-0.5,0.5)+((int)round(parms->dm[idm].order)%2)*0.5;
			real guard=parms->dm[idm].guard*MAX(dx, dy);
			map_t* map;
			if(parms->dbg.dmfullfov&&!parms->fit.square){//DM covers full fov
				real D=(parms->sim.fov*fabs(ht)+parms->aper.d+guard*2);
				long nx=D/dx+1;
				long ny=D/dy+1;
				map=mapnew(nx, ny, dx, dy);
				map->h=ht;
				map->ox+=offset*dx;
				mapcircle(map, D*0.5, 1);
			} else{
				create_metapupil(&map, 0, 0, parms->dirs, parms->aper.d, ht, dx, dy, offset, guard, 0, 0, 0, parms->fit.square);
			}
			info("    DM %d: grid is %ld x %ld\n", idm, NX(map), NY(map));
			P(recon->aloc, idm)=map2loc(map, 0);
			mapfree(map);
		}
	}
	recon->amap=mapcellnew(parms->ndm, 1);
	for(int idm=0; idm<parms->ndm; idm++){
		real ht=parms->dm[idm].ht;
		real offset=wrap2range(parms->dm[idm].offset,-0.5,0.5)+((int)round(parms->dm[idm].order)%2)*0.5;
		real dx=parms->dm[idm].dx;
		P(recon->aloc, idm)->iac=parms->dm[idm].iac;
		loc_create_map_npad(P(recon->aloc, idm), parms->fit.square?0:1, 0, 0);
		P(recon->amap, idm)=P(recon->aloc, idm)->map;
		P(recon->amap, idm)->h=ht;
		if(parms->fit.cachedm){
			const real dx2=parms->atmr.dx/parms->fit.pos;
			const real dy2=dx2;
			create_metapupil(&P(recon->acmap, idm), 0, 0, parms->dirs, parms->aper.d,
				ht, dx2, dy2, offset*dx/dx2, dx2, 0, 0, 0, parms->fit.square);
		}
	}

	recon->aimcc=dcellnew(ndm, 1);
	for(int idm=0; idm<ndm; idm++){
		P(recon->aimcc, idm)=loc_mcc_ptt(P(recon->aloc, idm), NULL);
		dinvspd_inplace(P(recon->aimcc, idm));
	}
	recon->anx=lnew(ndm, 1);
	recon->any=lnew(ndm, 1);
	recon->anloc=lnew(ndm, 1);
	for(int idm=0; idm<ndm; idm++){
		P(recon->anx, idm)=P(recon->amap, idm)->nx;
		P(recon->any, idm)=P(recon->amap, idm)->ny;
		P(recon->anloc, idm)=P(recon->aloc, idm)->nloc;
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
			P(recon->actstuck, idm)=loc_coord2ind(P(recon->aloc, idm), parms->dm[idm].actstuck);
		}
	}
	if(anyfloat){
		recon->actfloat=lcellnew(parms->ndm, 1);
		for(int idm=0; idm<ndm; idm++){
			if(!parms->dm[idm].actfloat) continue;
			P(recon->actfloat, idm)=loc_coord2ind(P(recon->aloc, idm), parms->dm[idm].actfloat);
		}
	}
	if(anystuck||anyfloat){
		for(int idm=0; idm<ndm; idm++){
			int nstuck=recon->actstuck?count_nonzero(P(recon->actstuck, idm)):0;
			int nfloat=recon->actfloat?count_nonzero(P(recon->actfloat, idm)):0;
			info2("    DM %d has %d stuck and %d floating actuators\n", idm, nstuck, nfloat);
		}
	}
	if(parms->save.setup){
		writebin(recon->aloc, "aloc");
		writebin(recon->amap, "amap");
	}
}
static void
setup_recon_amod(recon_t* recon, const parms_t* parms){
	const int ndm=parms->ndm;
	if(!ndm) return;
	recon->amod=dcellnew(ndm, 1);
	recon->anmod=lnew(ndm, 1);
	for(int idm=0; idm<ndm; idm++){
		int nmod=parms->recon.nmod;
		const long nloc=P(recon->aloc, idm)->nloc;
		switch(parms->recon.modal){
		case -2: {//dummy modal control, emulating zonal mode with identity modal matrix
			if(nmod&&nmod!=nloc){
				warning("recon.mod should be 0 or %ld when recon.modal=2 \n", nloc);
			}
			P(recon->amod, idm)=dnew(nloc, nloc);
			real val=sqrt(nloc);
			daddI(P(recon->amod, idm), val);
			dadds(P(recon->amod, idm), -val/nloc);
		}
			   break;
		case -1://zernike
		{
			if(!nmod) nmod=nloc;
			int rmax=floor((sqrt(1+8*nmod)-3)*0.5);
			P(recon->amod, idm)=zernike(P(recon->aloc, idm), 0, 0, rmax, 0);
		}
		break;
		case 1://Karhunen loeve. Don't limit number of modes here to make caching of G_M right.
			P(recon->amod, idm)=KL_vonkarman(P(recon->aloc, idm), 0, parms->atmr.L0);
			break;
		default:
			error("Invalid recon.modal");
		}
		P(recon->anmod, idm)=P(recon->amod, idm)->ny;
	}
}
/**
   Setup ray tracing operator from xloc to ploc for guide stars.
   @param recon
   @param parms
   @param atm		If set, atm[].vx,vy will be used for predictive offsetting.
*/

void setup_recon_HXW(recon_t* recon, const parms_t* parms, mapcell *atm){
	loc_t* ploc=recon->ploc;
	const int nwfsr=parms->nwfsr;
	const int npsr=recon->npsr;
	if(parms->load.HXW && !atm){
		warning("Loading saved HXW\n");
		if(recon->HXW){
			warning("recon->HXW already exists.\n");
			cellfree(recon->HXW);
		}
		recon->HXW=dspcellread("%s", parms->load.HXW);
		if(NX(recon->HXW)!=nwfsr||NY(recon->HXW)!=npsr){
			error("Loaded HXW has wrong cell size.\n");
		}
		dspcell* HXW=recon->HXW/*PDSPCELL*/;
		int nploc=ploc->nloc;
		for(int ips=0; ips<npsr; ips++){
			int nloc=P(recon->xloc, ips)->nloc;
			for(int iwfs=0; iwfs<nwfsr; iwfs++){
				if(!P(HXW, iwfs, ips)
					||P(HXW, iwfs, ips)->nx!=nploc
					||P(HXW, iwfs, ips)->ny!=nloc){
					error("Loaded HXW has wrong size.\n");
				}
			}
		}
	} else{
		TIC;tic;
		if(!recon->HXW) recon->HXW=dspcellnew(nwfsr, npsr);
		dspcell* HXW=recon->HXW/*PDSPCELL*/;
		for(int iwfsr=0; iwfsr<nwfsr; iwfsr++){
			int ipowfs=parms->wfsr[iwfsr].powfs;
			int shwfs=parms->powfs[ipowfs].type==WFS_SH;
			if(parms->recon.split!=2&&parms->powfs[ipowfs].skip){
				//don't need HXW for low order wfs that does not participate in tomography.
				continue;
			}
			const real delay=parms->sim.dt*(parms->powfs[ipowfs].dtrat+1+parms->sim.alhi);
			const real hs=parms->wfsr[iwfsr].hs;
			loc_t* loc=recon->ploc; //do not use ploc_tel. the rotation is handled in recon_utils.
			for(int ips=0; ips<npsr; ips++){
				const real  ht=P(recon->ht, ips);
				const real  scale=1.-ht/hs;
				real dispx=parms->wfsr[iwfsr].thetax*ht+(shwfs?scale*parms->wfsr[iwfsr].misregx:0);
				real dispy=parms->wfsr[iwfsr].thetay*ht+(shwfs?scale*parms->wfsr[iwfsr].misregy:0);
				if(atm){//point ahead to reduce lag.
					int ips2=P(parms->atmr.indps, ips);
					dispx+=P(atm, ips2)->vx*delay;
					dispy+=P(atm, ips2)->vy*delay;
				}
				cellfree(P(HXW, iwfsr, ips));
				P(HXW, iwfsr, ips)=mkh(P(recon->xloc, ips), loc, dispx, dispy, scale, 0);
			}
		}
		toc2("HXW");
	}
	if(parms->save.recon){
		writebin(recon->HXW, "HXW");
	}
	if(!recon->HXWtomo) recon->HXWtomo=dspcellnew(NX(recon->HXW), NY(recon->HXW));
	dspcell* HXWtomo=recon->HXWtomo/*PDSPCELL*/;
	dspcell* HXW=recon->HXW/*PDSPCELL*/;
	for(int iwfs=0; iwfs<nwfsr; iwfs++){
		int ipowfs=parms->wfsr[iwfs].powfs;
		if(!parms->powfs[ipowfs].skip){/*for tomography */
			for(int ips=0; ips<npsr; ips++){
				cellfree(P(HXWtomo, iwfs, ips));
				P(HXWtomo, iwfs, ips)=dspref(P(HXW, iwfs, ips));
			}
		}
	}
}

/**
   Setup gradient operator from ploc to wavefront sensors.
*/
static void
setup_recon_GP(recon_t* recon, const parms_t* parms, const aper_t *aper, const powfs_t *powfs){
	loc_t* ploc=recon->ploc;
	const int nwfsr=parms->nwfsr;
	recon->GP=dspcellnew(nwfsr, 1);
	if(parms->load.GP){
		warning("Loading saved GP\n");
		dspcell* GPload=dspcellread("%s", parms->load.GP);
		int assign=0;
		if(NX(GPload)==nwfsr){
			assign=1;
		} else if(NX(GPload)==parms->npowfs){
			assign=2;
		} else{
			error("GP loaded from %s has wrong size.\n", parms->load.GP);
		}
		for(int iwfs=0; iwfs<nwfsr; iwfs++){
			const int ipowfs=parms->wfsr[iwfs].powfs;
			const long nsa=P(recon->saloc, ipowfs)->nloc;
			P(recon->GP, iwfs)=dspref(P(GPload, assign==1?iwfs:ipowfs));
			if(P(recon->GP, iwfs)->nx!=nsa*2||P(recon->GP, iwfs)->ny!=ploc->nloc){
				error("Wrong saved GP[%d]: size is %ldx%ld, need %ldx%ld\n", iwfs,
					P(recon->GP, iwfs)->nx, P(recon->GP, iwfs)->ny, nsa*2, ploc->nloc);
			}
		}
		dspcellfree(GPload);
	} else{
		TIC;tic;
		for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
			if(parms->powfs[ipowfs].type!=WFS_SH) continue;
			const int iwfs0=P(parms->powfs[ipowfs].wfsr, 0);
			int share_gp=NX(powfs[ipowfs].amp)==1;
			/*
				During tomography operation, raytracing to ploc handles shift,
				an rotation step handles clocking. Higher order distortion is
				ignored. GP maps from ploc to saloc with no misregistration
				needed. This mimimizes non-zeros in GP and enables matrix free
				approach in the GPU.
			*/
	OMP_FOR(parms->powfs[ipowfs].nwfsr) 
			for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfsr; jwfs++){
				const int iwfs=P(parms->powfs[ipowfs].wfsr, jwfs);
				if(!share_gp||iwfs==iwfs0){
					//create gloc that aligns with saloc edges. The GP created this way has fewer elements.
					dmat* gamp=0;
					loc_t *gloc=0;
					shwfs_gloc_gamp(&gloc, &gamp, parms, aper, iwfs);
					if(parms->plot.setup>1){
						drawopd("Aperture", gloc, gamp, 0, "GP amplitude map", "x (m)", "y (m)", "gamp %d", iwfs);
					}
					dsp* gp=0;
					if(parms->powfs[ipowfs].gtype_recon==GTYPE_G){//Average tilt
						gp=mkg(ploc, gloc, gamp, P(recon->saloc, ipowfs), 
							PR(powfs[ipowfs].saa, jwfs), parms->powfs[ipowfs].saat,
							1, 0, 0, 1);
					} else if(parms->powfs[ipowfs].gtype_recon==GTYPE_Z){//Zernike fit
						dsp* ZS0=mkz(gloc, P(gamp), P(recon->saloc, ipowfs), 1, 1, 0, 0);
						dsp* H=mkh(ploc, gloc, 0, 0, 1, 0);
						gp=dspmulsp(ZS0, H, "nn");
						dspfree(H);
						dspfree(ZS0);
					} else{
						error("Invalid gtype_recon\n");
					}
					P(recon->GP, iwfs)=gp;
					locfree(gloc);
					dfree(gamp);
				}
			}
			//outside of OMP_FOR
			if(share_gp){
				for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfsr; jwfs++){
					const int iwfs=P(parms->powfs[ipowfs].wfsr, jwfs);
					P(recon->GP, iwfs)=dspref(P(recon->GP, iwfs0));
				}
			}
		}
		toc2("GP");
		if(parms->save.recon){
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
	const int nwfsr=parms->nwfsr;
	const int ndm=parms->ndm;
	if(parms->load.GA){
		warning("Loading GA from %s\n", parms->load.GA);
		recon->GA=dspcellread("%s", parms->load.GA);
		if(NX(recon->GA)!=nwfsr||NY(recon->GA)!=ndm)
			error("Wrong saved GA (%ldx%ld). Need (%dx%d)\n",
				NX(recon->GA), NY(recon->GA), nwfsr, ndm);
		for(int idm=0; idm<ndm; idm++){
			int nloc=P(recon->aloc, idm)->nloc;
			for(int iwfs=0; iwfs<nwfsr; iwfs++){
				int ipowfs=parms->wfsr[iwfs].powfs;
				if(parms->tomo.ahst_idealngs&&parms->powfs[ipowfs].lo){
					continue;
				}
				int nsa=P(recon->saloc, ipowfs)->nloc;
				if(P(recon->GA, iwfs, idm)->nx!=nsa*2||(!parms->recon.modal&&P(recon->GA, iwfs, idm)->ny!=nloc)
					||(parms->recon.modal&&P(recon->GA, iwfs, idm)->ny!=P(recon->amod, idm)->ny)){
					error("Wrong saved GA\n");
				}
			}
		}
	} else{
		TIC;tic;
		recon->GA=dspcellnew(nwfsr, ndm);
		if(parms->recon.modal){
			recon->GM=dcellnew(nwfsr, ndm);
		}
		for(int iwfsr=0; iwfsr<nwfsr; iwfsr++){
			int ipowfs=parms->wfsr[iwfsr].powfs;
			if(parms->tomo.ahst_idealngs&&parms->powfs[ipowfs].lo){
				continue;
			}
			/*if(parms->powfs[ipowfs].skip==2){//no need for TWFS [use for TWFS mode]
				continue;
			}*/
			const int shwfs=parms->powfs[ipowfs].type==WFS_SH;
			const real hs=parms->wfsr[iwfsr].hs;
			const loc_t* saloc=P(recon->saloc, ipowfs);
			for(int idm=0; idm<ndm; idm++){
				const real  ht=parms->dm[idm].ht;
				const real  scale=1.-ht/hs;
				const loc_t* aloc=P(recon->aloc, idm);
				const real dispx=parms->wfsr[iwfsr].thetax*ht+(shwfs?scale*parms->wfsr[iwfsr].misregx:0);
				const real dispy=parms->wfsr[iwfsr].thetay*ht+(shwfs?scale*parms->wfsr[iwfsr].misregy:0);

				if(parms->powfs[ipowfs].type==WFS_PY){//PWFS
					/*if(!parms->powfs[ipowfs].lo)*/
					{
						dmat* opdadd=0;
						/*if(!parms->recon.glao&&powfs[ipowfs].opdadd&&0){
								int wfsind=P(parms->powfs[ipowfs].wfsind, iwfs);
							opdadd=P(powfs[ipowfs].opdadd, wfsind);
						}*/
						if(!parms->recon.modal){
							info("\nPyWFS from aloc to saloc directly\n");
							dmat* tmp=pywfs_mkg(powfs[ipowfs].pywfs, aloc, parms->recon.distortion_dm2wfs[iwfsr+idm*nwfsr],
								0, opdadd, dispx, dispy);
							P(recon->GA, iwfsr, idm)=d2sp(tmp, dmaxabs(tmp)*1e-6);
							dfree(tmp);
						} else{
							info("\nPyWFS from amod to saloc directly\n");
							//We compute the GM for full set of modes so that it is cached only once.
							P(recon->GM, iwfsr, idm)=pywfs_mkg(powfs[ipowfs].pywfs, aloc, parms->recon.distortion_dm2wfs[iwfsr+idm*nwfsr],
								P(recon->amod, idm), opdadd, dispx, dispy);
						}
					}
				} else{//SHWFS
					char* input=parms->distortion.dm2wfs?parms->distortion.dm2wfs[iwfsr+idm*nwfsr]:0;
					char* calib=parms->recon.distortion_dm2wfs?parms->recon.distortion_dm2wfs[iwfsr+idm*nwfsr]:0;
					loc_t* loc=ploc;

					if(input&&saloc->nloc>4){//there is distortion input
						if(calib){//there is magical calibration input
							loc=loctransform(ploc, calib);
						} else{//determine distortion by fitting
							warning_once("dm2wfs: determine distortion by fitting \n");
							//First, simulate the poke matrix measurement process.
							loc_t* loc2=loctransform(ploc, input);
							dsp* H2=mkh(aloc, loc2, dispx, dispy, scale, 0);
							locfree(loc2);
							dsp* GA2=dspmulsp(P(recon->GP, iwfsr), H2, "nn"); //measured GA.
							dspfree(H2);
							dmat* calib2=loc_calib(GA2, aloc, saloc, dispx, dispy, scale, 2);
							dshow(calib2, "dm2wfs calib fitting");
							dspfree(GA2);
							loc=loctransform2(ploc, calib2);
							dfree(calib2);
						}
					}
					/*if(parms->wfsr[iwfsr].misregc){//absorbed in mkht
						if(loc==ploc){
							loc=locdup(ploc);
						}
						locrot(loc, -parms->wfsr[iwfsr].misregc);
					}*/
					dsp *Ht=mkht(aloc, loc, dispx, dispy, scale, -parms->wfsr[iwfsr].misregc);
					P(recon->GA, iwfsr, idm)=dspmulsp(P(recon->GP, iwfsr), Ht, "nt");
					dspfree(Ht);
					if(loc!=ploc){
						locfree(loc);
					}
					if(parms->recon.modal){
						dspmm(&P(recon->GM, iwfsr, idm), P(recon->GA, iwfsr, idm), P(recon->amod, idm), "nn", 1);
					}
				}
			}/*idm */
		}//iwfs
		toc2("GA");
	}
	if(parms->recon.modal&&parms->recon.nmod>0){
		for(int idm=0; idm<ndm; idm++){
			if(parms->recon.nmod<P(recon->amod, idm)->ny){
				info("DM %d:Reduce number of controlled modes from %ld to %d\n",
					idm, P(recon->amod, idm)->ny, parms->recon.nmod);
				dresize(P(recon->amod, idm), 0, parms->recon.nmod);
				P(recon->anmod, idm)=P(recon->amod, idm)->ny;
				for(int iwfs=0; iwfs<nwfsr; iwfs++){
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
			if(parms->save.setup){
				writebin(recon->actcpl, "lsr_actcpl");
			}
			act_stuck(recon->aloc, recon->actcpl, recon->actfloat);
			if(parms->lsr.actextrap){
				//when lor is enabled, the resulting matrix is much less sparse.
				recon->actextrap=act_extrap(recon->aloc, recon->actcpl, parms->lsr.actthres, 0);
				dspcell *GA2=0;
				dcellmm(&GA2, recon->GA, recon->actextrap, "nn", 1);
				dspcellfree(recon->GA);
				recon->GA=GA2;
				if(parms->save.setup){
					writebin(recon->actextrap, "lsr_actextrap");
				}
			} else if(recon->actfloat){//avoid commanding floating actuators
				act_float(recon->aloc, &recon->GA, NULL, recon->actfloat);
			}
		}
	}
	/*Create GAlo that only contains GA for low order wfs */
	recon->GAlo=cellnew(NX(recon->GA), NY(recon->GA));
	recon->GAhi=dspcellnew(NX(recon->GA), NY(recon->GA));
	if(parms->recon.modal) recon->GMhi=dcellnew(nwfsr, ndm);

	for(int idm=0; idm<ndm; idm++){
		for(int iwfs=0; iwfs<nwfsr; iwfs++){
			int ipowfs=parms->wfsr[iwfs].powfs;
			if(parms->powfs[ipowfs].lo
				||(parms->recon.split&&parms->nlopowfs==0&&!parms->powfs[ipowfs].trs)){/*for low order wfs */
				P(recon->GAlo, iwfs, idm)=(cell*)dspref(P(recon->GA, iwfs, idm));
			}
			if(!parms->powfs[ipowfs].skip){
				P(recon->GAhi, iwfs, idm)=dspref(P(recon->GA, iwfs, idm));
				if(parms->recon.modal){
					P(recon->GMhi, iwfs, idm)=dref(P(recon->GM, iwfs, idm));
				}
			}
		}
	}
	if(parms->save.recon){
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
	TIC;tic;
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	/*gradient from xloc. Also useful for lo WFS in MVST mode. */
		for(int ips=0; ips<npsr; ips++){
			P(GX, iwfs, ips)=dspmulsp(P(recon->GP, iwfs), P(HXW, iwfs, ips), "nn");
		}/*ips */
	}
	toc2("GX");
	recon->GXtomo=dspcellnew(NX(recon->GX), NY(recon->GX));
	dspcell* GXtomo=recon->GXtomo/*PDSPCELL*/;

	recon->GXlo=dspcellnew(NX(recon->GX), NY(recon->GX));
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
			dspmm(&P(recon->GFall, iwfs), P(recon->GP, parms->recon.glao?ipowfs:iwfs), opd, "nn", 1);
			if(parms->powfs[ipowfs].lo&&!parms->powfs[ipowfs].llt){
				P(recon->GFngs, iwfs)=dref(P(recon->GFall, iwfs));
			}
		}
		dfree(opd);
	}
	if(parms->save.recon){
		writebin(recon->GFall, "GFall");
	}
}
/**
   From radial order modes to gradients.
 */
static void
setup_recon_GR(recon_t* recon, const parms_t* parms){
	if(parms->itpowfs==-1&&!(parms->ilgspowfs!=-1&&parms->powfs[parms->ilgspowfs].dither==-1&&parms->powfs[parms->ilgspowfs].phytype_sim2==PTYPE_COG)){
		return;
	}
	int nlayer=1;
	//Need two layers when there are multiple TWFS or LGS WFS gradient offset needs projection adjustment
	if((parms->itpowfs!=-1 && parms->powfs[parms->itpowfs].nwfs>1)
		||(parms->ilgspowfs!=-1&&parms->powfs[parms->ilgspowfs].dither==-1&&parms->powfs[parms->ilgspowfs].phytype_sim2==PTYPE_COG&&parms->powfs[parms->ilgspowfs].nwfs>1)){
		nlayer=parms->ndm;
	}
	recon->GRall=dcellnew(parms->nwfs, nlayer);

	const int rmax=parms->recon.twfs_rmax?parms->recon.twfs_rmax:(parms->powfs[parms->itpowfs].order/2);
	const int rmin=parms->recon.twfs_rmin?parms->recon.twfs_rmin:3;
	const int zradonly=parms->recon.twfs_radonly;
	dbg("twfs mode is %s from order %d to %d on %d layers.\n", zradonly?"radial":"all modes", rmin, rmax, nlayer);
	for(int ilayer=0; ilayer<nlayer; ilayer++){
		const loc_t* loc=P(recon->aloc, ilayer);
		int rmin2=rmin;
		if(ilayer>0 && rmin2<3){//don't place those on upper layer.
			rmin2=3;
		}
		//must use aper.d here to make sure mode in different layers match in strength for TWFS.
		dmat* opd=zernike(loc, -parms->aper.d, rmin2, rmax, zradonly);

		for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
			const int ipowfs=parms->wfs[iwfs].powfs;
			if(parms->powfs[ipowfs].skip==2||parms->powfs[ipowfs].llt){
				if(P(recon->GA, iwfs, ilayer)){
					dspmm(&P(recon->GRall, iwfs, ilayer), P(recon->GA, iwfs, ilayer), opd, "nn", 1);
				}else{
					error("Please implement without GA\n");
				}
			}
		}
		dfree(opd);
	}
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		int ipowfs=parms->wfs[iwfs].powfs;
		if(parms->powfs[ipowfs].skip==2){//twfs
			int nlayer2=MIN(parms->powfs[ipowfs].nwfs, nlayer);
			if(parms->powfs[ipowfs].nwfs>1&&nlayer==1){
				error("recon.GRwfs should have more than 1 layer when there are multiple twfs.\n");
			}
			if(rmin<3){
				warning("rmin should be 3 for truth wfs.\n");
			}
			if(!recon->GRtwfs){
				recon->GRtwfs=dcellnew(parms->nwfs, nlayer2);
			}
			for(int ilayer=0; ilayer<nlayer2; ilayer++){
				P(recon->GRtwfs, iwfs, ilayer)=dref(P(recon->GRall, iwfs, ilayer));
			}
		}
		if(parms->powfs[ipowfs].llt&&parms->powfs[ipowfs].dither==-1&&parms->powfs[ipowfs].phytype_sim2==PTYPE_COG){
			int nlayer2=MIN(parms->powfs[ipowfs].nwfs, nlayer);
			if(parms->powfs[ipowfs].nwfs>1&&nlayer==1){
				error("recon.GRwfs should have more than 1 layer for sodium fitting projection.\n");
			}
			if(rmin>2){
				error("rmin should be 1 or 2 for sodium fitting projection\n");
			}
			if(!recon->GRlgs){
				recon->GRlgs=dcellnew(parms->powfs[ipowfs].nwfs, nlayer2);
			}
			int jwfs=P(parms->powfs[ipowfs].wfsind, iwfs);
			for(int ilayer=0; ilayer<nlayer2; ilayer++){
				P(recon->GRlgs, jwfs, ilayer)=dref(P(recon->GRall, iwfs, ilayer));
			}
		}
	}
	if(recon->GRlgs){
		//2021-10-15: Since we are not selecting modes, there is no need for high threshold
		//to high threshold makes the filtering ill formed
		const real thres=1e-14;
		info("RRlgs svd thres is %g\n", thres);
		recon->RRlgs=dcellpinv2(recon->GRlgs, NULL, thres);
	}
	if(parms->save.recon){
		writebin(recon->GRall, "twfs_GR");
		if(recon->RRlgs) writebin(recon->RRlgs, "twfs_RRlgs");
	}
}
/**
   Tilt removal from DM command. Used by filter.c
 */
static void
setup_recon_dmttr(recon_t* recon, const parms_t* parms){
	recon->DMTT=dcellnew(parms->ndm, 1);
	recon->DMPTT=dcellnew(parms->ndm, 1);
	/*if(!recon->actcpl && parms->nwfs>0){
	error("actcpl must not be null\n");
	}*/
	for(int idm=0; idm<parms->ndm; idm++){
		P(recon->DMTT, idm)=loc2mat(P(recon->aloc, idm), 0);
	}
	if(parms->dbg.recon_stuck){
		act_zero(recon->aloc, recon->DMTT, recon->actstuck);
	}
	for(int idm=0; idm<parms->ndm; idm++){
		P(recon->DMPTT, idm)=dpinv(P(recon->DMTT, idm), 0);
	}
	if(parms->save.recon){
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
		if(parms->powfs[ipowfs].trs
			||(parms->recon.split&&!parms->powfs[ipowfs].lo&&!parms->powfs[ipowfs].skip)
			||parms->powfs[ipowfs].dither==1){
			int nsa=P(recon->saloc, ipowfs)->nloc;
			dmat* TT=0;
			if(parms->powfs[ipowfs].type==WFS_SH){//SHWFS
				TT=dnew(nsa*2, 2);
				real* TTx=PCOL(TT, 0);
				real* TTy=PCOL(TT, 1);
				for(int isa=0; isa<nsa; isa++){
					TTx[isa]=1;
					TTy[isa]=0;
				}
				for(int isa=nsa; isa<nsa*2; isa++){
					TTx[isa]=0;
					TTy[isa]=1;
				}
			} else if(parms->powfs[ipowfs].type==WFS_PY){//PYWFS
				TT=dref(powfs[ipowfs].pywfs->GTT); //pywfs_tt(powfs[ipowfs].pywfs);
			} else{
				error("Invalid powfs.type\n");
			}
			if(parms->recon.glao){
				P(recon->TT, ipowfs, ipowfs)=ddup(TT);
			} else{
				for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
					int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
					P(recon->TT, iwfs, iwfs)=ddup(TT);
				}
			}
			dfree(TT);
		}
	}
	if(parms->save.recon){
		writebin(recon->TT, "TT");
	}
}
/**
   Operator to remove focus modes.
*/

static void
setup_recon_FF(recon_t* recon, const parms_t* parms){
	int nwfs=parms->nwfsr;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(parms->powfs[ipowfs].nwfs>0
			&& !parms->powfs[ipowfs].skip
			&& parms->powfs[ipowfs].llt
			&& parms->powfs[ipowfs].frs){
			const int nsa=P(recon->saloc, ipowfs)->nloc;
			dmat* FF=dnew(nsa*2, 1);
			/*saloc is lower left corner of subaperture. don't have to be the center. */
			memcpy(P(FF), P(recon->saloc, ipowfs)->locx, sizeof(real)*nsa);
			memcpy(P(FF)+nsa, P(recon->saloc, ipowfs)->locy, sizeof(real)*nsa);
			if(!recon->FF) recon->FF=dcellnew(nwfs, nwfs);
			dbg("powfs %d has focus mode in recon->FF\n", ipowfs);
			for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
				int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
				P(recon->FF, iwfs, iwfs)=dref(FF);
			}
			dfree(FF);
		}
	}
	if(parms->save.recon && recon->FF){
		writebin(recon->FF, "FF");
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
				if(dither_mode!=parms->powfs[ipowfs].dither||
					fabs(dither_amp-parms->powfs[ipowfs].dither_amp)>dither_amp*1e-5
					||dither_npoint!=parms->powfs[ipowfs].dither_npoint
					||dither_dtrat!=parms->powfs[ipowfs].dtrat){
					error("Multiple dither with different configuration is not supported\n");
				}
			}
			dither_mode=parms->powfs[ipowfs].dither;
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
		recon->dither_ra=dcellnew(parms->nwfs, parms->ndm);
		//recon->dither_rm=dcellnew(parms->nwfs, parms->ndm);
		recon->dither_rg=dcellnew(parms->nwfs, parms->nwfs);
		int DITHER_ND=6;
		READ_ENV_INT(DITHER_ND, 1, 3000); //number of mode bins
		int DITHER_MD2=1;
		READ_ENV_INT(DITHER_MD2, 1, 3000); // average # modes around first bin
		int DITHER_MD3=20;
		READ_ENV_INT(DITHER_MD3, 1, 3000); //average # modes in other bin
		for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
			const int ipowfs=parms->wfs[iwfs].powfs;
			if(parms->powfs[ipowfs].dither>1){
				const real ht=parms->dm[idm].ht+parms->dm[idm].vmisreg;
				const real dispx=ht*parms->wfs[iwfs].thetax;
				const real dispy=ht*parms->wfs[iwfs].thetay;

				const int dither_nd=parms->powfs[ipowfs].dither_mmd?MAX(1, DITHER_ND):1;
				const int dither_md=parms->powfs[ipowfs].dither_mmd?((parms->recon.nmod+dither_nd-1)/dither_nd):1; //number of modes per bin
				recon->dither_md=dither_md;
				const int dither_md2=MIN(dither_md, DITHER_MD2);
				const int dither_md3=MIN(dither_md, DITHER_MD3);

				//(parms->powfs[ipowfs].dither_mode2)?MAX(1,(parms->recon.nmod+md-1)/md):1;//number of dither modes. 2 for modal control
				dbg("dither mds are md=%d, %d, %d. nd=%d. dither_amp=%g\n", dither_md, dither_md2, dither_md3, dither_nd, dither_amp);

				dmat *grad=0;
				for(int id=0; id<dither_nd; id++){
					dmat *dither_rg=NULL;
					dmat *dither_ra=NULL;
					dmat *dither_m=NULL;
					int jm=id*dither_md;
					if (dither_nd==1){//single mode dithering uses zernike mode
						dither_m=zernike(P(recon->aloc, idm), 0, 0, 0, -dither_mode);
						dscale(dither_m, dither_amp);
						if(parms->powfs[ipowfs].type==1){//PWFS
							dfree(grad);
							grad=pywfs_mkg(powfs[ipowfs].pywfs, P(recon->aloc, idm),
								parms->recon.distortion_dm2wfs[iwfs+idm*parms->nwfs],
								dither_m, NULL, dispx, dispy);
						}else{
							error("Please implement for SHWFS\n");
						}
					}else if(dither_md>1){//average every md2 modes
						int md2=(id==0?dither_md2:dither_md3);
						if(md2<2){//a single amod per dithering mode
							if(id==0){
								jm=MAX(0,parms->powfs[ipowfs].dither-2);
								dbg("The first dithering mode uses amod %d directly\n", jm);
							}
							md2=1;
						}
						if(jm+md2>parms->recon.nmod){
							md2=parms->recon.nmod-jm;
						}
						real dither_amp2=dither_amp*(1-id/dither_nd);//higher order modes have higher gain, so reduce its strength.
						real dither_amp3=dither_amp2/sqrt(md2);//split into modes within the bin
						dmat *wt=dnew(md2, 1);
						dset(wt, 1);

						dmat *mtmp=drefcols(P(recon->amod,idm),jm, md2);
						dmat *gtmp=drefcols(P(recon->GM, iwfs, idm), jm, md2);
						dmm(&dither_m, 0, mtmp, wt, "nn", dither_amp3);
						dmm(&grad, 0, gtmp, wt, "nn", dither_amp3);
						dfree(mtmp);
						dfree(gtmp);
						dfree(wt);
					}

					dither_rg=dpinv(grad, P(recon->saneai, iwfs, iwfs));
					dither_ra=dpinv(dither_m, 0);
					if(id==0){
						P(recon->dither_m, idm)=dither_m; dither_m=NULL;
						P(recon->dither_rg, iwfs, iwfs)=dither_rg; dither_rg=NULL;
						P(recon->dither_ra, iwfs, idm)=dither_ra; dither_ra=NULL;
					}else{
						{
							dmat *tmp=P(recon->dither_ra, iwfs, idm);
							P(recon->dither_ra, iwfs, idm)=dcat(tmp, dither_ra, 1);
							dfree(tmp);
							dfree(dither_ra);
						}
						{
							dmat *tmp=P(recon->dither_rg, iwfs, iwfs);
							P(recon->dither_rg, iwfs, iwfs)=dcat(tmp, dither_rg, 1);
							dfree(tmp);
							dfree(dither_rg);
						}
						{
							dadd(&P(recon->dither_m, idm), 1, dither_m, 1);
							dfree(dither_m);
						}
					}
				}

				//dmm(&P(recon->dither_rm, iwfs, idm), 0, P(recon->dither_ra, iwfs, idm), P(recon->amod, idm), "nn", 1.);
				dfree(grad);
			}
		}


		if(parms->save.setup){
			writebin(recon->dither_m, "dither_m");
			writebin(recon->dither_ra, "dither_ra");
			writebin(recon->dither_rg, "dither_rg");
			//writebin(recon->dither_rm, "dither_rm");
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
	if(parms->cn2.pair&&NX(parms->cn2.pair)>0&&!recon->cn2est){
	/*setup CN2 Estimator. It determines the reconstructed layer heigh can be fed to the tomography */
		recon->cn2est=cn2est_prepare(parms, powfs);
	}
	recon->saloc=loccellnew(parms->npowfs, 1);
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		P(recon->saloc, ipowfs)=locref(powfs[ipowfs].saloc);
	}
	recon->wfsr=parms->wfsr;
	recon->nwfsr=parms->nwfsr;
	recon->ngrad=lnew(parms->nwfsr, 1);
	for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
		const int ipowfs=parms->wfsr[iwfs].powfs;
		if(!parms->powfs[ipowfs].skip){
			P(recon->ngrad, iwfs)=P(recon->saloc, ipowfs)->nloc*2;
		}
	}
	/*to be used in tomography. */
	recon->nthread=NTHREAD;
	/*setup pupil coarse grid for gradient operator*/
	setup_recon_ploc(recon, parms);
	/*setup DM actuator grid */
	setup_recon_aloc(recon, parms);
	if(parms->recon.modal){
		//DM modal control matrix.
		setup_recon_amod(recon, parms);
	}
	/*Grid for DM fitting*/
	setup_recon_floc(recon, parms);
	/*Gradient operators*/
	setup_recon_GP(recon, parms, aper, powfs);
	//TT Removal
	setup_recon_TT(recon, parms, powfs);
	//Global or Differential focus removal.
	setup_recon_FF(recon, parms);

	if(recon->FF){
		recon->TTF=dcellcat_each(recon->TT, recon->FF, 2);
	} else{
		recon->TTF=dcellref(recon->TT);
	}
	if(parms->recon.alg==0&&!parms->sim.idealtomo){//tomography parameters
		if(parms->cn2.tomo&&recon->cn2est){
			/*Use cn2 estimation results for tomography. Use its ht to build
			  reconstructor matrices.*/
			cn2est_t* cn2est=recon->cn2est;
			recon->ht=dref(cn2est->htrecon);
			recon->os=dref(cn2est->os);
			recon->wt=dref(P(cn2est->wtrecon, 0));
			/*the following will be updated later in simulation. */
			if(parms->cn2.keepht){
				for(int ips=0; ips<NX(recon->wt); ips++){
					P(recon->wt, ips)=P(parms->atmr.wt, ips);
				}
			} else{
				dset(recon->wt, 1./NX(recon->wt));/*evenly distributed.  */
			}
		} else{/*use input information from atmr */
			recon->wt=dnew(parms->atmr.nps, 1);
			recon->ht=dnew(parms->atmr.nps, 1);
			recon->os=dnew(parms->atmr.nps, 1);
			for(int ips=0; ips<NX(recon->wt); ips++){
				P(recon->wt, ips)=P(parms->atmr.wt, ips);
				P(recon->ht, ips)=P(parms->atmr.ht, ips);
				P(recon->os, ips)=P(parms->atmr.os, ips);
			}
		}
		recon->r0=parms->atmr.r0;
		recon->L0=parms->atmr.L0;

		/*sampling of xloc */
		recon->dx=dnew(NX(recon->ht), 1);
		for(int iht=0; iht<NX(recon->ht); iht++){
			real scale=1.0-P(recon->ht, iht)/parms->atmr.hs;
			P(recon->dx, iht)=(parms->atmr.dx/P(recon->os, iht))*scale;
		}
		/*number of reconstruction layers */
		recon->npsr=NX(recon->ht);
		/*setup atm reconstruction layer grid */
		setup_recon_xloc(recon, parms);
		/*setup xloc/aloc to WFS grad */
		setup_recon_HXW(recon, parms, NULL);
		setup_recon_GX(recon, parms);//uses HXW.
		dspcellfree(recon->HXW);/*only keep HXWtomo for tomography */
	}

	return recon;
}
/**
   That may depend on GPU data.
 */
void setup_recon_prep_ga(recon_t* recon, const parms_t* parms, const aper_t* aper, const powfs_t* powfs){
	info2("\n%sSetting up reconstructor%s\n\n", GREEN, BLACK);
	setup_recon_GA(recon, parms, powfs);//PWFS uses GPU data.
	setup_recon_GF(recon, parms);//GF depends on GA.
	setup_recon_GR(recon, parms);

	if(parms->recon.split){
		ngsmod_prep(parms, recon, aper, powfs);
	}
	setup_recon_dmttr(recon, parms);
}

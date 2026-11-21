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
#include "skyc.h"
#include "genstars.h"
#include "types.h"
#include "photon.h"
/**
   \file genstars.c Generate random stars from besancon star counting model. The
   data file contains J and H magnitudes of all stars within 5 square degree
   FoV. The location of the stars are assumed to be randomly
   distributed. Therefore, the possible number of stars within any sub-field is
   a poisson distribution with mean equal to the average density multiplied with
   area.
 */

/**
   Generate stars for nsky star fields from star catalog.

   @param[in] nsky Number of star fields wanted.
   @param[in] lat Galactic latitude (degree).
   @param[in] lon Galactic longitude (degree).
   @param[in] catscl Scale of the catalog star count. If negative, its absolute value is the maximum number of J<=19
   stars, and the star fields are evenly distributed among 0 to that number of J<=19 stars.
   @param[in] fov Diameter of the patrol field of view (arcsec).
   @param[in] nwvl Number of wavelengths.
   @param[in] wvls Wavelengths (m). Must match the catalog J, H, K bands.
   @param[in] maglimit If nonzero, remove stars dimmer than this in the first wavelength.
   @param[in,out] rstat Random stream.
   @return Cell array of nsky x 1. Each cell contains a (2+nwvl) x nstar array of location (rad) and magnitudes, or
   NULL if the field has no star. Newly allocated; the caller frees it.
*/
dcell* genstars(long nsky, real lat, real lon, real catscl, real fov, int nwvl, real* wvls, real maglimit,
	rand_t* rstat){
	char fn[80];   //file containing list of stars over the entire fov.
	real cat_fov=0;//Catalog FoV in degree squared
	int Jind=-1;   //Index of J band
	const real wvl_cat[3]={1.25e-6, 1.65e-6, 2.2e-6};
	if(nwvl<=3){
		for(int iwvl=0; iwvl<nwvl; iwvl++){
			if(fabs(wvls[iwvl]-wvl_cat[iwvl])>1e-7){
				error("wvl[%d]=%g does not agree with catalog\n", iwvl, wvls[iwvl]);
			}
		}
	}
	/*if(nwvl==2){
		snprintf(fn, 80, "besancon/JH_5sqdeg_lat%g_lon%g_besancon.bin", lat, lon);
		cat_fov=5.0;//5 arc-degree squared. 
		Jind=0;
	} else */
	if(nwvl<=3){
		snprintf(fn, 80, "besancon/JHK_5sqdeg_lat%g_lon%g_besancon.bin", lat, lon);
		cat_fov=5.0;/*5 arc-degree squared. */
		Jind=0;
	} else{
		Jind=-1;
		error("We only have stars for J+H and J+H+K band. Please fill this part\n");
	}
	info("Loading star catalogue from %s\n", fn);
	dmat* catalog=dread("%s", fn);
	if(catalog->ny<nwvl){
		error("Catalogue has not have enough wavelength\n");
	}
	long ntot=catalog->nx;
	long nsky0=0;
	dcell* res=dcellnew(nsky, 1);
	real fov22=pow(fov/2*AS2RAD, 2);


	real navg0=M_PI*pow(fov/2./3600., 2)/cat_fov*ntot;
	if(catscl>0){//regular sky coverage sim
		real navg=navg0*catscl;
		info("Average number of stars in patrol fov: %g, after scaled by %g\n", navg, catscl);
		/*generate nstart && magnitude according to distribution.*/
		for(long isky=0; isky<nsky; isky++){
			long nstar=randp(rstat, navg);
			if(nstar==0) continue;
			P(res,isky)=dnew(nwvl+2, nstar);
			dmat* pres=P(res,isky);
			for(long istar=0; istar<nstar; istar++){
				long ind=round(ntot*randu(rstat));/*randomly draw a star index in the catlog */
				if(ind>=ntot) ind=ntot-1;/*round() can give ntot*/
				for(int iwvl=0; iwvl<nwvl; iwvl++){
					P(pres, 2+iwvl, istar)=P(catalog, ind, iwvl);
				}
			}
		}
	} else{
	/*instead of doing draws on nb of stars, we scan all possibilities and
	  assemble the curve in postprocessing.  catscl is negative, with
	  absolute value indicating the max number of J<=19 stars to consider*/
		long nmax=round(-catscl);//maximum number of J<=19 stars. set to 10
		nsky0=nsky/nmax;//number of random instances for each number of J<=19 stars
		if(nsky0*nmax!=nsky){
			error("nsky=%ld, has to be dividable by max # of stars=%ld", nsky, nmax);
		}
		int counti[nmax];//record count in each bin
		memset(counti, 0, sizeof(int)*nmax);
		int count=0;
		while(count<nsky){
			long nstar=randp(rstat, navg0);
			if(nstar==0) continue;
			dmat* tmp=dnew(nwvl+2, nstar);
			dmat* pres=tmp;
			int J19c=0;
			for(long istar=0; istar<nstar; istar++){
				long ind=round((ntot-1)*randu(rstat));
				for(int iwvl=0; iwvl<nwvl; iwvl++){
					P(pres, 2+iwvl, istar)=P(catalog, ind, iwvl);
				}
				if(P(pres, 2+Jind, istar)<=19){
					J19c++;
				}
			}
			//J19c=0 is ok. Do not skip.
			if(J19c<nmax&&counti[J19c]<nsky0){
				long isky=counti[J19c]+(J19c)*nsky0;
				P(res,isky)=dref(tmp);
				count++;
				counti[J19c]++;
			}
			dfree(tmp);
		}
	}
	long ndrop=0;
	/*Fill in the coordinate*/
	for(long isky=0; isky<nsky; isky++){
		if(!P(res,isky)) continue;
		long nstar=P(res,isky)->ny;
		dmat* pres=P(res,isky);
		long jstar=0;
		for(long istar=0; istar<nstar; istar++){
			/*randomly draw the star location. */
			real r=sqrt(fov22*randu(rstat));
			real th=2*M_PI*randu(rstat);
			P(pres, 0, jstar)=r*cos(th);
			P(pres, 1, jstar)=r*sin(th);
			if(jstar!=istar){
				for(int iwvl=0; iwvl<nwvl; iwvl++){
					P(pres, 2+iwvl, jstar)=P(pres, 2+iwvl, istar);
				}
			}
			if(!maglimit || P(pres, 2, istar)<=maglimit){
				jstar++;
			}
		}
		if(jstar<nstar){
			ndrop+=(nstar-jstar);
			if(jstar==0){
				dfree(P(res, isky));
			}else{
				dresize(P(res, isky), -1, jstar);
			}
		}
	}
	if(ndrop){
		info("Removed %ld stars dimmer than the magnitude limit %g at %g\n", ndrop, maglimit, wvl_cat[0]);
	}
	dfree(catalog);
	return res;
}
/**
   Total J and H band flux of a star. With J band only (nrow==3), use the J band flux.

   @param[in] p    The star (x, y, J, H, ...).
   @param[in] nrow Number of rows (2+nwvl) of the star.
   @return The total flux.
*/
static real star_flux(const real* p, long nrow){
	const real Z_J=3.7666e9;
	const real Z_H=2.7206e9;
	return Z_J*pow(10, -0.4*p[2])+(nrow>3?Z_H*pow(10, -0.4*p[3]):0);
}
/**
   Compare two stars with J and H band magnitudes by total flux.

   @param[in] p1 First star (x, y, J, H, ...).
   @param[in] p2 Second star (x, y, J, H, ...).
   @return -1 if p1 is brighter than p2, 1 if dimmer, 0 if equal.
 */
static int sortfun_jh(const void* p1, const void* p2){
	real tot1=star_flux((const real*)p1, 4);
	real tot2=star_flux((const real*)p2, 4);
	return tot1<tot2?1:(tot1>tot2?-1:0);
}
/**
   Compare two stars with J band magnitude only by flux.

   @param[in] p1 First star (x, y, J).
   @param[in] p2 Second star (x, y, J).
   @return -1 if p1 is brighter than p2, 1 if dimmer, 0 if equal.
 */
static int sortfun_j(const void* p1, const void* p2){
	real tot1=star_flux((const real*)p1, 3);
	real tot2=star_flux((const real*)p2, 3);
	return tot1<tot2?1:(tot1>tot2?-1:0);
}
/**
   Sort the stars in each star field by total J and H band flux from brightest to dimmest.
   When only one band is available, sort by the J band flux.

   @param[in,out] stars Star fields as returned by genstars(). Sorted in place.
*/
void sortstars(dcell* stars){
	for(long isky=0; isky<stars->nx*stars->ny; isky++){
		if(!P(stars,isky)||!P(stars,isky)->ny) continue;
		const long nrow=P(stars,isky)->nx;
		if(nrow<3){
			error("Star field %ld has %ld rows, expect at least 3 (x, y, J).\n", isky, nrow);
		}
		qsort(P(P(stars,isky)), P(stars,isky)->ny, nrow*sizeof(real), nrow>3?sortfun_jh:sortfun_j);
	}
}

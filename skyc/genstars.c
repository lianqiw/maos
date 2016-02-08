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
   The sort function for stars. Sort stars according total flux.
 */
static double Z_J=3.7666e9;
static double Z_H=2.7206e9;
static int sortfun(const double *p1, const double *p2){
    double tot1=Z_J*pow(10,-0.4*p1[2])+Z_H*pow(10,-0.4*p1[3]);/*tot flux */
    double tot2=Z_J*pow(10,-0.4*p2[2])+Z_H*pow(10,-0.4*p2[3]);
    return tot1<tot2?1:-1;
}
/**
   Generate stars for nsky star fields from star catalog.

   \return a cell array of nskyx1, each cell contains (2+nwvl) x nstar array of
location, and magnitudes.  */
dcell *genstars(long nsky,         /**<number of star fields wanted*/
		double lat,        /**<galactic latitude.*/
		double lon,        /**<galactic longitude*/
		double catscl,     /**<Scale the catlog star count.*/
		double fov,        /**<diameter of the patrol field of view in arcsec.*/
		int nwvl,          /**<number of wavelength*/
		double *wvls,      /**<wavelength vector*/
		rand_t *rstat /**<random stream*/
){
    char fn[80];
    double cat_fov=0;/*catalogue fov */
    int Jind=-1;
    if(nwvl==2 && fabs(wvls[0]-1.25e-6)<1.e-10 && fabs(wvls[1]-1.65e-6)<1.e-10){
	snprintf(fn,80,"besancon/JH_5sqdeg_lat%g_lon%g_besancon.bin", lat, lon);
	cat_fov=5.0;/*5 arc-degree squared. */
	Jind=0;
    }else if(nwvl==3 && fabs(wvls[0]-1.25e-6)<1.e-10 && fabs(wvls[1]-1.65e-6)<1.e-10 && fabs(wvls[2]-2.2e-6)<1.e-10){
	snprintf(fn,80,"besancon/JHK_5sqdeg_lat%g_lon%g_besancon.bin", lat, lon);
	cat_fov=5.0;/*5 arc-degree squared. */
	Jind=0;
    }else{
	Jind=-1;
	error("We only have stars for J+H and J+H+K band. Please fill this part\n");
    }
    info2("Loading star catalogue from %s\n",fn);
    dmat *catalog=dread("%s",fn);
    if(catalog->ny!=nwvl){
	error("Catalogue and wanted doesn't match\n");
    }
    long ntot=catalog->nx;
    long nsky0=0;
    dcell *res=cellnew(nsky,1);
    PDMAT(catalog, pcatalog);
    double fov22=pow(fov/2/206265,2);


    double navg0=M_PI*pow(fov/2./3600.,2)/cat_fov * ntot;
    if(catscl>0){//regular sky coverage sim
	double navg=navg0*catscl;
	info2("Average number of stars: %g, after scaled by %g\n", navg, catscl);
	/*generate nstart && magnitude according to distribution.*/
	for(long isky=0; isky<nsky; isky++){
	    long nstar=randp(rstat, navg);
	    if(nstar==0) continue;
	    res->p[isky]=dnew(nwvl+2, nstar);
	    PDMAT(res->p[isky],pres);
	    for(long istar=0; istar<nstar; istar++){
		long ind=round(ntot*randu(rstat));/*randomly draw a star index in the catlog */
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    pres[istar][2+iwvl]=pcatalog[iwvl][ind];
		}
	    }
	}
    }else{
	/*instead of doing draws on nb of stars, we scan all possibilities and
	  assemble the curve in postprocessing.  catscl is negative, with
	  absolute value indicating the max number of J<=19 stars to consider*/
	long nmax=round(-catscl);
	nsky0=nsky/nmax;
	if(nsky0*nmax!=nsky){
	    error("nsky=%ld, has to be dividable by max # of stars=%ld", nsky, nmax);
	}
	int counti[nmax];//record count in each bin
	memset(counti, 0, sizeof(int)*nmax);
	int count=0;
	while(count<nsky){
	    long nstar=randp(rstat, navg0);
	    if(nstar==0) continue;
	    dmat *tmp=dnew(nwvl+2, nstar);
	    PDMAT(tmp, pres);
	    int J19c=0;
	    for(long istar=0; istar<nstar; istar++){
		long ind=round(ntot*randu(rstat));
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    pres[istar][2+iwvl]=pcatalog[iwvl][ind];
		}
		if(pres[istar][2+Jind]<=19){
		    J19c++;
		}
	    }
	    //J19c=0 is ok. Do not skip.
	    if(J19c<nmax && counti[J19c]<nsky0){
		int isky=counti[J19c]+(J19c)*nsky0;
		res->p[isky]=dref(tmp);
		count++;
		counti[J19c]++;
	    }
	    dfree(tmp);
	}
    }
    /*Fill in the coordinate*/
    for(long isky=0; isky<nsky; isky++){
	if(!res->p[isky]) continue;
	long nstar=res->p[isky]->ny;
	PDMAT(res->p[isky],pres);
	for(long istar=0; istar<nstar; istar++){
	    /*randomly draw the star location. */
	    double r=sqrt(fov22*randu(rstat));
	    double th=2*M_PI*randu(rstat);
	    pres[istar][0]=r*cos(th);
	    pres[istar][1]=r*sin(th);
	}
    }
    dfree(catalog);
    return res;
}
/**
   Sort the stars with J band magnitude from brightest to dimmest.
*/
void sortstars(dcell *stars){
    for(long isky=0; isky<stars->nx*stars->ny; isky++){
	if(!stars->p[isky]) continue;
	qsort(stars->p[isky]->p, stars->p[isky]->ny,stars->p[isky]->nx*sizeof(double),
	      (int(*)(const void*, const void *))sortfun);
    }
}

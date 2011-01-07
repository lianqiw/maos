/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
static int sortfun(const double *p1, const double *p2){
    double tot1=Z_J*pow(10,-0.4*p1[2])+Z_H*pow(10,-0.4*p1[3]);//tot flux
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
    char temp[80];
    double cat_fov;//catalogue fov
    if(nwvl==2 && fabs(wvls[0]-1.25e-6)<1.e-10 && fabs(wvls[1]-1.65e-6)<1.e-10){
	snprintf(temp,80,"besancon/JH_5sqdeg_lat%g_lon%g_besancon.bin.gz", lat, lon);
	cat_fov=5.0;//5 arc-degree squared.
    }else{
	error("We only have stars for J+H band. Please fill this part\n");
    }
    char *fn=find_file(temp);
    info("Loading star catalogue from \n%s\n",fn);
    dmat *catalog=dread("%s",fn); free(fn);
    if(catalog->ny!=nwvl){
	error("Catalogue and wanted doesn't match\n");
    }
    long ntot=catalog->nx;
    double navg=M_PI*pow(fov/2./3600.,2)/cat_fov * ntot * catscl;
    info("Average number of stars: %g, after scaled by %g\n", navg, catscl);
    dcell *res=dcellnew(nsky,1);
    PDMAT(catalog, pcatalog);
    double fov22=pow(fov/2/206265,2);
    for(long isky=0; isky<nsky; isky++){
	long nstar=randp(rstat, navg);
	if(nstar>0){
	    res->p[isky]=dnew(nwvl+2, nstar);
	    PDMAT(res->p[isky],pres);
	    for(long istar=0; istar<nstar; istar++){
		long ind=round(ntot*randu(rstat));//randomly draw a star index in the catlog
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    pres[istar][2+iwvl]=pcatalog[iwvl][ind];
		}
		//randomly draw the star location.
		double r=sqrt(fov22*randu(rstat));
		double th=2*M_PI*randu(rstat);
		pres[istar][0]=r*cos(th);
		pres[istar][1]=r*sin(th);
	    }
	    //sort the stars from brigtest to dimmest
	    qsort(res->p[isky]->p, res->p[isky]->ny, res->p[isky]->nx*sizeof(double), 
		  (int(*)(const void*, const void*))sortfun);
	}else{
	    res->p[isky]=NULL;
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

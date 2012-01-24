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

#include "common.h"
#include "misc.h"
#include "loc.h"
#include "dmat.h"
#include "cmat.h"
#include "dsp.h"
#include "bin.h"
#include "mkz.h"
/**
   \file mkg.c
   Contains function that creates ztilt gradient operator
*/
/**
   Returns the transpose of a ztilt gradient operator that converts the OPDs defined
   on xloc to subapertures defines on saloc.
 */
dsp * mkzt(loc_t* xloc, double *amp, loc_t *saloc, 
	      int saorc, double scale, double *displace)
{
    /*compute ztilt influence function from xloc to saloc
      saorc: SALOC is subaperture origin or center. 
      1: origin (lower left corner), 
      0: center.
    */
    long nsa=saloc->nloc;
    double dsa=saloc->dx;
    loc_create_map_npad(xloc,ifloor(dsa/xloc->dx));
    double dx1=1./xloc->dx;
    double dx2=scale*dx1;
    double poffset[2];
    poffset[0]=(displace[0]-xloc->map->ox+saorc*dsa*0.5*scale)*dx1;
    poffset[1]=(displace[1]-xloc->map->oy+saorc*dsa*0.5*scale)*dx1;
    double dsa2=dsa*0.5*dx2;
    long nmax=(dsa2*2+2)*(dsa2*2+2);
    long *ind=calloc(nmax, sizeof(long));
    loc_t *sloc=calloc(1, sizeof(loc_t));
    sloc->dx=xloc->dx;
    sloc->locx=calloc(nmax,sizeof(double));
    sloc->locy=calloc(nmax,sizeof(double));
    double *amploc=NULL;
    if(amp) amploc=calloc(nmax,sizeof(double));
    const long (*xmaps)[xloc->map->nx]=(const long (*)[xloc->map->nx])xloc->map->p;
    dsp*zax=spnew(xloc->nloc,nsa,xloc->nloc);
    dsp*zay=spnew(xloc->nloc,nsa,xloc->nloc);
    long xcount=0,ycount=0;
    spint *xpp=zax->p;
    spint *xpi=zax->i;
    double *xpx=zax->x;
    
    spint *ypp=zay->p;
    spint *ypi=zay->i;
    double *ypx=zay->x;
    const double *locx=xloc->locx;
    const double *locy=xloc->locy;
    double *slocx=sloc->locx;
    double *slocy=sloc->locy;
    for(int isa=0; isa<nsa; isa++){
	/*center of subaperture when mapped onto XLOC*/
	double scx=saloc->locx[isa]*dx2+poffset[0];
	double scy=saloc->locy[isa]*dx2+poffset[1];
	int count=0;
	/*find points that belongs to this subaperture. */
	for(int iy=iceil(scy-dsa2);
	    iy<ifloor(scy+dsa2);iy++){
	    for(int ix=iceil(scx-dsa2);
		ix<ifloor(scx+dsa2);ix++){
		int ii=xmaps[iy][ix];
		if(ii){
		    ii--;
		    ind[count]=ii;
		    slocx[count]=locx[ii];
		    slocy[count]=locy[ii];
		    if(amp) amploc[count]=amp[ii];
		    count++;
		}
	    }
	}
	/*locwrite(sloc,"sloc_isa%d",isa); */
	/*writedbl(amploc,count,1,"amploc_isa%d",isa); */
	sloc->nloc=count;
	dmat *mcc=loc_mcc_ptt(sloc,amploc);
	/*dwrite(mcc,"mcc_isa%d",isa); */
	dinvspd_inplace(mcc);
	/*dwrite(mcc,"imcc_isa%d",isa); */
	double (*MCC)[3]=(double(*)[3])mcc->p;
	xpp[isa]=xcount;
	ypp[isa]=ycount;
	for(int ic=0; ic<count; ic++){
	    double xx=MCC[1][0]+MCC[1][1]*slocx[ic]+MCC[1][2]*slocy[ic];
	    double yy=MCC[2][0]+MCC[2][1]*slocx[ic]+MCC[2][2]*slocy[ic];
	    if(amp){
		xx*=amploc[ic];
		yy*=amploc[ic];
	    }
	    xpi[xcount]=ind[ic];
	    xpx[xcount]=xx;
	    xcount++;
	    ypi[ycount]=ind[ic];
	    ypx[ycount]=yy;
	    ycount++;
	}
	dfree(mcc);
    }
    xpp[nsa]=xcount;
    ypp[nsa]=ycount;
    locfree(sloc);
    free(ind);

    spsetnzmax(zax,xcount);
    spsetnzmax(zay,ycount);
    dsp*ZAT=spcat(zax,zay,1);
    spfree(zax);
    spfree(zay);
    if(amp) free(amploc);
    loc_free_map(xloc);
    return ZAT;
}
/**
   Returns a ztilt gradient operator that converts the OPDs defined
   on xloc to subapertures defines on saloc.
 */
dsp *mkz(loc_t* xloc, double *amp,loc_t *saloc, 
	      int saorc, double scale, double *displace)
{
    dsp *zt=mkzt(xloc,amp,saloc,saorc,scale,displace);
    dsp *z=sptrans(zt);
    spfree(zt);
    return z;
}

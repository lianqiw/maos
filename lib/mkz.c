/*
  Copyright 2009-2018 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#include "../math/mathdef.h"
#include "mkz.h"

/**
   Returns the transpose of a ztilt gradient operator that converts the OPDs defined
   on xloc to subapertures defines on saloc.
 */
dsp * mkzt(loc_t* xloc, double *amp, loc_t *saloc, 
	   int saorc, double scale, double dispx, double dispy)
{
    /*compute ztilt influence function from xloc to saloc
      saorc: SALOC is subaperture origin or center. 
      1: origin (lower left corner), 
      0: center.
    */
    long nsa=saloc->nloc;
    double dsa=saloc->dx;
    double dx1=1./xloc->dx;
    double dx2=scale*dx1;
    double dy1=1./xloc->dy;
    double dy2=scale*dy1;
    loc_create_map(xloc);
    map_t *map=xloc->map;
    dispx=(dispx-map->ox+saorc*dsa*0.5*scale)*dx1;
    dispy=(dispy-map->oy+saorc*dsa*0.5*scale)*dy1;
    double dsa2=dsa*0.5*dx2;
    long nmax=(dsa2*2+2)*(dsa2*2+2);
    long *ind=mycalloc(nmax,long);
    loc_t *sloc=locnew(nmax, xloc->dx, xloc->dy);
    double *amploc=NULL;
    if(amp) amploc=mycalloc(nmax,double);

    dsp*zax=dspnew(xloc->nloc,nsa,xloc->nloc);
    dsp*zay=dspnew(xloc->nloc,nsa,xloc->nloc);
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
	double scx=saloc->locx[isa]*dx2+dispx;
	double scy=saloc->locy[isa]*dy2+dispy;
	int count=0;
	/*find points that belongs to this subaperture. */
	for(int iy=iceil(scy-dsa2); iy<ifloor(scy+dsa2);iy++){
	    for(int ix=iceil(scx-dsa2); ix<ifloor(scx+dsa2);ix++){
		int ii=loc_map_get(map, ix, iy);
		if(ii>0){
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
	/*writebin(mcc,"mcc_isa%d",isa); */
	dinvspd_inplace(mcc);
	/*writebin(mcc,"imcc_isa%d",isa); */
	xpp[isa]=xcount;
	ypp[isa]=ycount;
	for(int ic=0; ic<count; ic++){
	    double xx=IND(mcc,0,1)+IND(mcc,1,1)*slocx[ic]+IND(mcc,2,1)*slocy[ic];
	    double yy=IND(mcc,0,2)+IND(mcc,1,2)*slocx[ic]+IND(mcc,2,2)*slocy[ic];
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

    dspsetnzmax(zax,xcount);
    dspsetnzmax(zay,ycount);
    dsp*ZAT=dspcat(zax,zay,1);
    dspfree(zax);
    dspfree(zay);
    if(amp) free(amploc);
    return ZAT;
}
/**
   Returns a ztilt gradient operator that converts the OPDs defined
   on xloc to subapertures defines on saloc.
 */
dsp *mkz(loc_t* xloc, double *amp,loc_t *saloc, 
	 int saorc, double scale, double dispx, double dispy)
{
    dsp *zt=mkzt(xloc,amp,saloc,saorc,scale,dispx, dispy);
    dsp *z=dsptrans(zt);
    dspfree(zt);
    return z;
}

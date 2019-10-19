/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#include "mkw.h"

#define nres 50
/**
   compute the W0, W1 for bilinear influence function.
   The normalization is very important. We need to have sum(amp)==1;
   W0 and W1 can not be scaled because: 
   1) they scale the same way as amp, 
   2) we are using W0-W1*W1', one is linear, one is quadratic. 
   So the result will be different.
    
   for the result, we must have 
   1) sum(W0,1)==W1, sum(W0,2)==W1, 2) sum(W1)==1
   2009-11-20: bug found: 
   using a circular gray pixel amplitude map, 
   sum(W1) is less than 1 if loc is not enough large.
   This method is not good! Prefer the mkw_circular instead
   which follows closely the loas method. 
   The LAOS W0/W1 gives very good performance in NGS modes.

   This function does not work well. 
*/
void mkw_amp(loc_t *loc,real *amp,dsp **W0,dmat **W1){

    long nloc=loc->nloc;
    real constamp=0;
    if(amp){
	real sumamp=dvecdot(amp,NULL,NULL,nloc);
	if(fabs(sumamp-1)>1.e-6){
	    error("amp must be normalized to sum to 1\n");
	}
    }else{
	constamp=1./(real)nloc;
    }
    *W0=dspnew(nloc,nloc,9*nloc);
    loc_create_map(loc);
    map_t *map=loc->map;
    const real ox=map->ox;
    const real oy=map->oy;
    real idx=1./loc->dx;
    spint *W0p=(*W0)->p;
    spint *W0i=(*W0)->i;
    real *W0x=(*W0)->x;
    long count=0;
    real *amp0=NULL;
    if(amp){
	amp0=amp-1;/*offset by -1 */
    }
    *W1=dnew(nloc,1);
    real *W1p=(*W1)->p;
    for(int iloc=0; iloc<nloc ;iloc++){
	int ix=(int)round((loc->locx[iloc]-ox)*idx);
	int iy=(int)round((loc->locy[iloc]-oy)*idx);
	real camp=0;
	if(amp) camp=amp[iloc];
	real wt2=0;
	W0p[iloc]=count;
	for(int jy=iy-1; jy<iy+2; jy++){
	    for(int jx=ix-1;jx<ix+2; jx++){
		int aloc1=loc_map_get(map, jx, jy);
		if(aloc1<=0) continue; else aloc1--;
		real bamp=0; if(amp) bamp=amp[aloc1]; 
		real wt=0;
		if(abs(ix-jx)==1 && abs(iy-jy)==1){
		    /*corner to corner */
		    int bloc0=loc_map_get(map, ix, jy);
		    int bloc1=loc_map_get(map, jx, iy);
		    if(bloc0>0 && bloc1>0){
			if(amp){
			    wt+=(camp+bamp+amp0[bloc0]+amp0[bloc1])/144.;
			}else{
			    wt+=constamp*4./144.;
			}
		    }
		}else if(abs(iy-jy)==1){
		    /*neighbor, up/down */
		    for(int ioff=-1; ioff<2; ioff+=2){
			int bloc0=loc_map_get(map, ix+ioff, jy);
			int bloc1=loc_map_get(map, jx+ioff, iy);
			if(bloc0>0 && bloc1>0){
			    if(amp){
				wt+=(camp+bamp)/48.+(amp0[bloc0]+amp0[bloc1])/144.;
			    }else{
				wt+=constamp*8./144.;
			    }
			}
		    }
		}else if(abs(ix-jx)==1){
		    /*neighbor, left/right */
		    for(int joff=-1; joff<2; joff+=2){
			int bloc0=loc_map_get(map, ix, jy+joff);
			int bloc1=loc_map_get(map, jx, iy+joff);
			if(bloc0>0 && bloc1>0){
			    if(amp){
				wt+=(camp+bamp)/48.+(amp0[bloc0]+amp0[bloc1])/144.;
			    }else{
				wt+=constamp*8./144.;
			    }
			}
		    }
		}else{
		    /*same point. */
		    /*loop over four corners; */
		    for (int joff=-1; joff<2; joff+=2){
			for(int ioff=-1; ioff<2; ioff+=2){
			    int bloc0=loc_map_get(map, ix, iy+joff);
			    int bloc1=loc_map_get(map, ix+ioff, iy);
			    int cloc0=loc_map_get(map, ix+ioff, iy+joff);
			    if(bloc0>0 && bloc1>0 && cloc0>0){
				if(amp){
				    wt+=camp/16.+(amp0[bloc0]+amp0[bloc1])/48.
					+amp0[cloc0]/144;
				    wt2+=camp/9.+(amp0[bloc0]+amp0[bloc1])/18.
					+amp0[cloc0]/36;
				}else{
				    wt+=constamp*1./9.;
				    wt2+=constamp*0.25;
				}
			    }
			}
		    }
		}
		if(wt>0){
		    W0i[count]=aloc1;
		    W0x[count]=wt;
		    count++;
		}
	    }/*jx */
	}/*jy */
	W1p[iloc]=wt2;
    }/*iloc */
    W0p[nloc]=count;
    dspsetnzmax(*W0,count);
    real sumW1=dsum(*W1);
    if(fabs(sumW1-1)>1.e-12){
	real sc=1./sumW1;
	warning("Sum W1 is not equal to 1. "
		"Probably loc grid is too small. "
		"Rescaling W0, W1\n");
	real *p=(*W1)->p;
	for(int i=0; i<(*W1)->nx; i++){
	    p[i]*=sc;
	}
	p=(*W0)->x;
	for(int i=0; i<(*W0)->nzmax; i++){
	    p[i]*=sc;
	}
    }
}
/**
   calculate the integral of a box with two points on corner to corner.  inside
   a circle at (icx,icy) with circle icr.  */

static real calcwtcorner(int ix, int iy, 
		    int jx, int jy,
		    real icx, real icy, real icr){
    if(ix==jx || iy==jy) error("Invalid\n");
    real icr2=icr*icr;
    if(pow(ix-icx,2)+pow(iy-icy,2)<icr2
       &&pow(ix-icx,2)+pow(jy-icy,2)<icr2
       &&pow(jx-icx,2)+pow(iy-icy,2)<icr2
       &&pow(jx-icx,2)+pow(jy-icy,2)<icr2){
	/*The box is fully inside. */
	return 1./36.;
    }
    if(pow(ix-icx,2)+pow(iy-icy,2)>icr2
       &&pow(ix-icx,2)+pow(jy-icy,2)>icr2
       &&pow(jx-icx,2)+pow(iy-icy,2)>icr2
       &&pow(jx-icx,2)+pow(jy-icy,2)>icr2){
	/*The box is fully outside. */
	return 0;
    }
    real dres=1./(real)nres;
    real ddx=(real)(jx-ix)/(real)nres;
    real ddy=(real)(jy-iy)/(real)nres;
    real wt=0.;
    
    for(int iiy=0; iiy<nres; iiy++){
	real dy=((real)iiy+0.5)*ddy;/*dy may be negative */
	for(int iix=0; iix<nres; iix++){
	    real dx=((real)iix+0.5)*ddx;
	    if(pow(ix+dx-icx,2)+pow(iy+dy-icy,2)<icr2){
		/*inside circle. */
		wt+=(1.-fabs(dx))*(1.-fabs(dy))*fabs(dx*dy);
	    }
	}
    }
    return wt*dres*dres;
}
/**
   calculate the integral of a box with two points on left/right inside a circle
   at (icx,icy) with circle icr.  */
static real calcwtlr(int ix, int iy, 
		int jx, int jy,
		real icx, real icy, real icr){

    if(iy!=jy || ix==jx) error("Invalid\n");
    real icr2=icr*icr;
    if(pow(ix-icx,2)+pow(iy+1-icy,2)<icr2
       &&pow(ix-icx,2)+pow(iy-1-icy,2)<icr2
       &&pow(jx-icx,2)+pow(iy+1-icy,2)<icr2
       &&pow(jx-icx,2)+pow(iy-1-icy,2)<icr2){
	/*The box is fully inside. */
	return 1./9.;
    }
    if(pow(ix-icx,2)+pow(iy+1-icy,2)>icr2
       &&pow(ix-icx,2)+pow(iy-1-icy,2)>icr2
       &&pow(jx-icx,2)+pow(iy+1-icy,2)>icr2
       &&pow(jx-icx,2)+pow(iy-1-icy,2)>icr2){
	/*The box is fully outside. */
	return 0;
    }
    real dres=1./(real)nres;
    real ddx=(real)(jx-ix)/(real)nres;
    /*real ddy=(real)(jy-iy)/(real)nres; */
    real ddy=1./(real)nres;
    real wt=0.;
    
    for(int iiy=0; iiy<nres; iiy++){
	real dy=((real)iiy+0.5)*ddy;
	for(int iix=0; iix<nres; iix++){
	    real dx=((real)iix+0.5)*ddx;
	    /*upper */
	    if(pow(ix+dx-icx,2)+pow(iy+dy-icy,2)<icr2){
		/*inside circle. */
		wt+=(1.-fabs(dx))*pow(1.-fabs(dy),2)*fabs(dx);
	    }
	    /*lower */
	    if(pow(ix+dx-icx,2)+pow(iy-dy-icy,2)<icr2){
		/*inside circle. */
		wt+=(1.-fabs(dx))*pow(1.-fabs(dy),2)*fabs(dx);
	    }
	}
    }
    return wt*dres*dres;
}
/**
   calculate the integral of a box with two points on bottom/top.  inside a
   circle at (icx,icy) with circle icr.
*/
static real calcwtud(int ix, int iy, 
		int jx, int jy,
		real icx, real icy, real icr){
    if(ix!=jx) error("Invalid\n");
    real icr2=icr*icr;
    if(pow(ix-1-icx,2)+pow(iy-icy,2)<icr2
       &&pow(ix-1-icx,2)+pow(jy-icy,2)<icr2
       &&pow(ix+1-icx,2)+pow(iy-icy,2)<icr2
       &&pow(ix+1-icx,2)+pow(jy-icy,2)<icr2){
	/*The box is fully inside. */
	return 1./9.;
    }
    if(pow(ix-1-icx,2)+pow(iy-icy,2)>icr2
       &&pow(ix-1-icx,2)+pow(jy-icy,2)>icr2
       &&pow(ix+1-icx,2)+pow(iy-icy,2)>icr2
       &&pow(ix+1-icx,2)+pow(jy-icy,2)>icr2){
	/*The box is fully outside. */
	return 0;
    }
    
    real dres=1./(real)nres;
    /*real ddx=(real)(jx-ix)/(real)nres; */
    real ddy=(real)(jy-iy)/(real)nres;
    real ddx=1./(real)nres;
    real wt=0.;
    
    for(int iiy=0; iiy<nres; iiy++){
	real dy=((real)iiy+0.5)*ddy;
	for(int iix=0; iix<nres; iix++){
	    real dx=((real)iix+0.5)*ddx;
	    /*upper */
	    if(pow(ix+dx-icx,2)+pow(iy+dy-icy,2)<icr2){
		/*inside circle. */
		wt+=pow(1.-fabs(dx),2)*(1.-fabs(dy))*fabs(dy);
	    }
	    /*lower */
	    if(pow(ix-dx-icx,2)+pow(iy+dy-icy,2)<icr2){
		/*inside circle. */
		wt+=pow(1.-fabs(dx),2)*(1.-fabs(dy))*fabs(dy);
	    }
	}
    }
    return wt*dres*dres;
}
/**
  calculate the integral of a box with the same point.  inside a circle at
  (icx,icy) with circle icr.  */
static real calcwtcenter(int ix, int iy, int jx, int jy,
			   real icx, real icy, real icr){
  
    if(ix!=jx || iy!=jy) error("Invalid\n");
    real icr2=icr*icr;
    if(pow(ix+1-icx,2)+pow(iy+1-icy,2)<icr2
       &&pow(ix+1-icx,2)+pow(iy-1-icy,2)<icr2
       &&pow(ix-1-icx,2)+pow(iy+1-icy,2)<icr2
       &&pow(ix-1-icx,2)+pow(iy-1-icy,2)<icr2){
	/*The box is fully inside. */
	return 4./9.;
    }
    if(pow(ix+1-icx,2)+pow(iy+1-icy,2)>icr2
       &&pow(ix+1-icx,2)+pow(iy-1-icy,2)>icr2
       &&pow(ix-1-icx,2)+pow(iy+1-icy,2)>icr2
       &&pow(ix-1-icx,2)+pow(iy-1-icy,2)>icr2){
	/*The box is fully outside. */
	return 0;
    }
    real dres=1./(real)nres;
    real ddy=1./(real)nres;
    real ddx=ddy;
    real wt=0.;
    for(int iiy=0; iiy<nres; iiy++){
	real dy=((real)iiy+0.5)*ddy;
	for(int iix=0; iix<nres; iix++){
	    real dx=((real)iix+0.5)*ddx;
	    /*lower left */
	    if(pow(ix-dx-icx,2)+pow(iy-dy-icy,2)<icr2){
		/*inside circle. */
		wt+=pow((1.-fabs(dx))*(1.-fabs(dy)),2);
	    }
	    /*lower right */
	    if(pow(ix+dx-icx,2)+pow(iy-dy-icy,2)<icr2){
		/*inside circle. */
		wt+=pow((1.-fabs(dx))*(1.-fabs(dy)),2);
	    }
	    /*upper left */
	    if(pow(ix-dx-icx,2)+pow(iy+dy-icy,2)<icr2){
		/*inside circle. */
		wt+=pow((1.-fabs(dx))*(1.-fabs(dy)),2);
	    }
	    /*upper right */
	    if(pow(ix+dx-icx,2)+pow(iy+dy-icy,2)<icr2){
		/*inside circle. */
		wt+=pow((1.-fabs(dx))*(1.-fabs(dy)),2);
	    }
	}
    }
    return wt*dres*dres;
}
/**
   compute the W0, W1 for bilinear influence function for a circular aperture of
   radius cr, so that for OPD vector A defined on grid loc, the piston removed
   wavefront error, calculated as A'*(W0-W1*W1')*A is equal to the wavefront
   error for a continuous OPD that are interpolated bi-linearly using the OPD on
   the grid.*/
void mkw_circular(loc_t *loc, /**<[in] grid coordinate*/
		  real cx, /**<[in] center of circle, x*/
		  real cy, /**<[in] center of circle, y*/
		  real cr, /**<[in] circle radius*/
		  dsp **W0,  /**<[out] sparse W0*/
		  dmat **W1  /**<[out] dense  W1*/){
    long nloc=loc->nloc;
    *W0=dspnew(nloc,nloc,9*nloc);
    loc_create_map(loc);
    map_t *map=loc->map;
    const real ox=map->ox;
    const real oy=map->oy;
    spint *W0p=(*W0)->p;
    spint *W0i=(*W0)->i;
    real *W0x=(*W0)->x;
    long count=0;
    *W1=dnew(nloc,1);
    real *W1p=(*W1)->p;
    real idx=1./loc->dx;
    /*center of circle in the map. */
    real icx=(cx-ox)*idx;
    real icy=(cy-oy)*idx;
    /*radius of circle in unit of the map. */
    real icr=cr*idx;
    real sc=1./(M_PI*pow(icr,2));
    for(int iloc=0; iloc<nloc ;iloc++){
	real xx=loc->locx[iloc];
	real yy=loc->locy[iloc];
	int ix=(int)round((xx-ox)*idx);
	int iy=(int)round((yy-oy)*idx);
	real wt2=0;
	real wt=0;
	W0p[iloc]=count;
	for(int jy=iy-1; jy<iy+2; jy++){
	    for(int jx=ix-1;jx<ix+2; jx++){
		int aloc1=loc_map_get(map, jx, jy);
		if(aloc1<=0) continue; else aloc1--;
		/*
		  We calculate the integration for different relative position
		  of the points. 
		 */
		if(abs(ix-jx)==1 && abs(iy-jy)==1){
		    /*corner to corner */
		    wt=calcwtcorner(ix,iy,jx,jy,icx,icy,icr);
		}else if(abs(iy-jy)==1){
		    /*up to down */
		    wt=calcwtud(ix,iy,jx,jy,icx,icy,icr);
		}else if(abs(ix-jx)==1){
		    /*left to right */
		    wt=calcwtlr(ix,iy,jx,jy,icx,icy,icr);
		}else{ 
		    /*self */
		    wt=calcwtcenter(ix,iy,jx,jy,icx,icy,icr);
		}
		if(wt>0){
		    W0i[count]=aloc1;
		    W0x[count]=wt*sc;
		    count++;
		}
		wt2+=wt;
	    }/*jx */
	}/*jy */
	W1p[iloc]=wt2*sc;
    }/*iloc */
    W0p[nloc]=count;
    dspsetnzmax(*W0,count);
    real sumW1=dsum(*W1);
    if(fabs(sumW1-1)>1.e-12){
	sc=1./sumW1;
	real *p=(*W1)->p;
	for(int i=0; i<(*W1)->nx; i++){
	    p[i]*=sc;
	}
	p=(*W0)->x;
	for(int i=0; i<(*W0)->nzmax; i++){
	    p[i]*=sc;
	}
    }
    /*loc_free_map(loc);*/
}
/**
   compute the  W0, W1 for bilinear influence function for a annular aperture of
   radius inner radius cri, and outer radius cro. see mkw().
 */
void mkw_annular(loc_t *loc, /**<[in] grid coordinate */
		 real cx,  /**<[in] center of circle, x*/
		 real cy,  /**<[in] center of circle, y*/
		 real cri, /**<[in] inner circular hole radius*/
		 real cro, /**<[in] outer circle radius*/
		 dsp **W0,   /**<[out] sparse W0*/
		 dmat **W1   /**<[out] dense  W1*/
		 ){
    if(cri<loc->dx){/*filled aperture */
	mkw_circular(loc, cx, cy, cro, W0, W1);
    }else{
	dsp *W0o = NULL;
	dmat *W1o= NULL;
	dsp *W0i = NULL;
	dmat *W1i= NULL;
	mkw_circular(loc, cx, cy, cri, &W0i, &W1i);/*inner */
	mkw_circular(loc, cx, cy, cro, &W0o, &W1o);/*outer */

	/*Calculate the factor applied in mkw_circular */
	real idx=1./loc->dx;
	real sco=1/(M_PI*pow(cro*idx, 2));
	real sci=1/(M_PI*pow(cri*idx, 2));

	*W1 = W1o;
	dadd(W1, 1./sco, W1i, -1/sci);/*cancel the factor applied in mkw_circular */
	real sc=1./dsum(*W1);
	dscale(*W1, sc);
	/*cancel the factor applied in mkw_circular and apply the new factor */
	*W0 = dspadd2(W0o, sc/sco, W0i, -sc/sci);

	dspfree(W0o);
	dspfree(W0i);
	dfree(W1i);
    }
}

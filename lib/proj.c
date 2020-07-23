/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

/*
  Compute the projection of a tilted M3 onto pupil.
  input map is located at height ht and tilted at angle
  thetax, thetay
  focus of the guide star is at heigt hs, angle alphax,alphay

  Notice that, for guide stars located at direction alphax, alphay, the input
  parameter betax, betay is betax=alphax*f_tel/f_exit, betay=alphay*f_tel/f_exit
  where f_tel is the telescope focal length. f_exit is the distance between the
  exit pupil to the focal plane. Do not put negative signs. the reflexive is
  taken into account in the function.

 */

#include "accphi.h"
#include "proj.h"
/*const real pi=3.1415926535897932384626433832795; */
static inline real cosangle(real a[3], real b[3]){
    return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])
	/sqrt((a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
	      *(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]));
}
/**
   Ray tracing for tilted surface like M3
 */
void proj_rect_grid(real* phiout,  /**<[out] output OPD*/
                    const rmap_t *mapin, /**<mapin: M3 surface map, NOT OPD. */
                    real thetax,   /**<tilt of surface along x, 0 or pi/2*/
                    real thetay,   /**<tilt of surface along y, pi/2 or 0*/
		    const loc_t *locout, /**<output grid*/
                    const real ratiox,   /**<pupil x scaling.*/
                    const real ratioy,   /**<pupil y scaling.*/
		    const real *ampout,  /**<amplitude of output grid*/
		    real sc,       /**<scaling of OPD*/
                    real hs,       /**<distance from exit pupil to object image*/
                    real ht,       /**<distance from exit pupil to M3 surface*/
		    real betax,    /**<real beam angle*/
                    real betay     /**<real beam angle*/
    ){
    const int wrapx1 = mapin->nx;
    const int wrapy1 = mapin->ny;
    const int wrapx = wrapx1-1;
    const int wrapy = wrapy1-1;
    real offx=hs*betax;
    real offy=hs*betay;
    if(ratiox<0) offx=-offx;
    if(ratioy<0) offy=-offy;
    const real dx_in1 = 1./mapin->dx;
    const real dy_in1 = 1./mapin->dy;
    real a0x=(-offx/sin(thetax)-mapin->ox)*dx_in1;
    real a0y=(-offy/sin(thetay)-mapin->oy)*dy_in1;
    real ddx=(hs-ht)*dx_in1;
    real ddy=(hs-ht)*dy_in1;
  
    int nplocx,nplocy,nplocx1,nplocy1;
    if(fabs(thetax-M_PI*0.5)>M_PI*.45 || fabs(thetay-M_PI*0.5)>M_PI*0.45){
	info("Tilting angle is too much\n");
	return;
    }
    
    real vm3[3];
    vm3[0]=-sin(M_PI/2-thetax);
    vm3[1]=-sin(M_PI/2-thetay);
    vm3[2]=-sqrt(1.-vm3[0]*vm3[0]-vm3[1]*vm3[1]);
    real vi[3];
    vi[2]=-hs;
    real sc2;
    for(int iloc=0; iloc<locout->nloc; iloc++){
	if(ampout && fabs(ampout[iloc])<1.e-10)
	    continue;/*skip points that has zero amplitude */
	real alx=atan2(locout->locx[iloc]*ratiox+offx,hs);
	real aly=atan2(locout->locy[iloc]*ratioy+offy,hs);
	real btx=thetax-alx;
	real bty=thetay-aly;
	real dplocx=ddx*sin(alx)/sin(btx)+a0x;
	real dplocy=ddy*sin(aly)/sin(bty)+a0y;
	vi[0]=locout->locx[iloc]*ratiox+offx;
	vi[1]=locout->locy[iloc]*ratioy+offy;
	SPLIT(dplocx,dplocx,nplocx);
	SPLIT(dplocy,dplocy,nplocy);
	if(nplocx<0||nplocx>=wrapx||nplocy<0||nplocy>=wrapy){
	    continue;
	}else{
	    nplocx1=nplocx+1;
	    nplocy1=nplocy+1;
	}
	sc2=sc*cosangle(vi,vm3);
	/*sc2=sc*0.707; */
	phiout[iloc]+=sc2*(P(mapin,nplocx,nplocy)*(1.-dplocx)*(1.-dplocy)
			  +P(mapin,nplocx1,nplocy)*(dplocx)*(1.-dplocy)
			  +P(mapin,nplocx,nplocy1)*(1.-dplocx)*(dplocy)
			  +P(mapin,nplocx1,nplocy1)*(dplocx)*(dplocy));
    }
}

/**
   Wraps proj_rect_grid for M3.
*/
void m3proj(const rmap_t *tsurf, dmat *opd, const loc_t *locout, real thetax, real thetay, real hs){
    const real alx=tsurf->txdeg/180*M_PI;
    const real aly=tsurf->tydeg/180*M_PI;
    const real ftel=tsurf->ftel;
    const real fexit=tsurf->fexit;
    const real fsurf=tsurf->fsurf;
    const real mag=fexit/ftel;
    const real scalex=-mag;
    const real scaley=mag;
    const real scaleopd=-2;
    const real het=fexit-fsurf;/*distance between exit pupil and M3. */
    //distance from focus to image offinite distance object
    real d_img_focus=1./(1./ftel-1./hs)-ftel;
    
    real d_img_exit=fexit+d_img_focus;
		
    /*2010-04-02: do not put - sign */
    real bx=thetax*(d_img_focus+ftel)/d_img_exit;
    real by=thetay*(d_img_focus+ftel)/d_img_exit;
    proj_rect_grid(opd->p,tsurf,alx,aly,locout,scalex,scaley, NULL,scaleopd, d_img_exit, het, bx, by);
}
/**
   A convenient wrapper for m3proj() to be called from matlab or python
 */
dmat *m3proj2(dmat *mapin_0, const char *header, const loc_t *locout, real thetax, real thetay, real hs){
    free(mapin_0->header);
    mapin_0->header=strdup(header);
    rmap_t *mapin=d2rmap(mapin_0);
    dmat *opd=dnew(locout->nloc, 1);
    m3proj(mapin, opd, locout, thetax, thetay, hs);
    cellfree(mapin);
    cellfree(mapin_0);
    return opd;
}

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

/*
  Compute the projection of a tilted M3 onto pupil.
  input loc locin is located at height ht and tilted at angle
  thetax, thetay
  focus of the guide star is at heigt hs, angle alphax,alphay

  Notice that, for guide stars located at direction alphax, alphay, the input
  parameter betax, betay is betax=alphax*f_tel/f_exit, betay=alphay*f_tel/f_exit
  where f_tel is the telescope focal length. f_exit is the distance between the
  exit pupil to the focal plane. Do not put negative signs. the reflexive is
  taken into account in the function.

 */
#include <math.h>
#include "accphi.h"
#include "proj.h"
/*const double pi=3.1415926535897932384626433832795; */
static inline double cosangle(double a[3], double b[3]){
    return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])
	/sqrt((a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
	      *(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]));
}

void proj_rect_grid(rectmap_t *mapin, double thetax, double thetay,
		    const loc_t *locout,const double ratiox, const double ratioy,
		    const double *ampout, double* phiout, 
		    double sc, double hs, double ht,
		    double betax, double betay){
    /*
      input parameters:
      mapin: M3 surface map. 
      thetax, thetay: tilt of surface map along x or y. one has to be pi/2
      locout: Pupil plan opd grid.
      ratio: scaling of pupil plan to exit pupil plane.
      sc: scaling of opd. don't include projection
      hs: distance from exit pupil to m3.
      betax,betay: beam angle from exit pupil to focal plane.
     */
    const double (*phiin)[mapin->nx]=(void*)mapin->p;
    const int wrapx1 = mapin->nx;
    const int wrapy1 = mapin->ny;
    const int wrapx = wrapx1-1;
    const int wrapy = wrapy1-1;
    double offx=hs*betax;
    double offy=hs*betay;
    if(ratiox<0) offx=-offx;
    if(ratioy<0) offy=-offy;
    const double dx_in1 = 1./mapin->dx;
    const double dy_in1 = 1./mapin->dy;
    double a0x=(-offx/sin(thetax)-mapin->ox)*dx_in1;
    double a0y=(-offy/sin(thetay)-mapin->oy)*dy_in1;
    double ddx=(hs-ht)*dx_in1;
    double ddy=(hs-ht)*dy_in1;
  
    int nplocx,nplocy,nplocx1,nplocy1;
    if(fabs(thetax-M_PI*0.5)>M_PI*.45 || fabs(thetay-M_PI*0.5)>M_PI*0.45){
	info2("Tilting angle is too much\n");
	return;
    }
    
    double vm3[3];
    vm3[0]=-sin(M_PI/2-thetax);
    vm3[1]=-sin(M_PI/2-thetay);
    vm3[2]=-sqrt(1.-vm3[0]*vm3[0]-vm3[1]*vm3[1]);
    double vi[3];
    vi[2]=-hs;
    double sc2;
    for(int iloc=0; iloc<locout->nloc; iloc++){
	if(ampout && fabs(ampout[iloc])<1.e-10)
	    continue;/*skip points that has zero amplitude */
	double alx=atan2(locout->locx[iloc]*ratiox+offx,hs);
	double aly=atan2(locout->locy[iloc]*ratioy+offy,hs);
	double btx=thetax-alx;
	double bty=thetay-aly;
	double dplocx=ddx*sin(alx)/sin(btx)+a0x;
	double dplocy=ddy*sin(aly)/sin(bty)+a0y;
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
	phiout[iloc]+=sc2*(phiin[nplocy][nplocx]*(1.-dplocx)*(1.-dplocy)
			  +phiin[nplocy][nplocx1]*(dplocx)*(1.-dplocy)
			  +phiin[nplocy1][nplocx]*(1.-dplocx)*(dplocy)
			  +phiin[nplocy1][nplocx1]*(dplocx)*(dplocy));
    }
}

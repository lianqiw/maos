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
/**
   Cosine of the angle between two 3-d vectors.

   @param[in] a The first vector.
   @param[in] b The second vector.
   @return The cosine of the angle between a and b.
*/
static inline real cosangle(real a[3], real b[3]){
	return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])
		/sqrt((a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
			*(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]));
}
/**
   Ray tracing for tilted surface like M3

   @param[in,out] phiout Output OPD, accumulated into.
   @param[in] mapin M3 surface map, NOT OPD.
   @param[in] thetax Tilt of surface along x (radian), 0 or pi/2.
   @param[in] thetay Tilt of surface along y (radian), pi/2 or 0.
   @param[in] locout Output grid.
   @param[in] ratiox Pupil x scaling.
   @param[in] ratioy Pupil y scaling.
   @param[in] ampout Amplitude of output grid. Points with zero amplitude are skipped. Can be NULL.
   @param[in] sc Scaling of OPD.
   @param[in] hs Distance from exit pupil to object image (m).
   @param[in] ht Distance from exit pupil to M3 surface (m).
   @param[in] betax Real beam angle along x (radian).
   @param[in] betay Real beam angle along y (radian).
*/
void proj_rect_grid(real* phiout, const rmap_t* mapin, real thetax, real thetay, const loc_t* locout, const real ratiox,
	const real ratioy, const real* ampout, real sc, real hs, real ht, real betax, real betay){
	const long wrapx1=NX(mapin);
	const long wrapy1=NY(mapin);
	const long wrapx=wrapx1-1;
	const long wrapy=wrapy1-1;
	real offx=hs*betax;
	real offy=hs*betay;
	if(ratiox<0) offx=-offx;
	if(ratioy<0) offy=-offy;
	const real dx_in1=1./mapin->dx;
	const real dy_in1=1./mapin->dy;
	real a0x=(-offx/sin(thetax)-mapin->ox)*dx_in1;
	real a0y=(-offy/sin(thetay)-mapin->oy)*dy_in1;
	real ddx=(hs-ht)*dx_in1;
	real ddy=(hs-ht)*dy_in1;

	int nplocx, nplocy, nplocx1, nplocy1;
	if(fabs(thetax-M_PI*0.5)>M_PI*.45||fabs(thetay-M_PI*0.5)>M_PI*0.45){
		info2("Tilting angle is too much\n");
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
		if(ampout&&fabs(ampout[iloc])<1.e-10)
			continue;/*skip points that has zero amplitude */
		real alx=atan2(locout->locx[iloc]*ratiox+offx, hs);
		real aly=atan2(locout->locy[iloc]*ratioy+offy, hs);
		real btx=thetax-alx;
		real bty=thetay-aly;
		real dplocx=ddx*sin(alx)/sin(btx)+a0x;
		real dplocy=ddy*sin(aly)/sin(bty)+a0y;
		vi[0]=locout->locx[iloc]*ratiox+offx;
		vi[1]=locout->locy[iloc]*ratioy+offy;
		SPLIT(dplocx, dplocx, nplocx);
		SPLIT(dplocy, dplocy, nplocy);
		if(nplocx<0||nplocx>=wrapx||nplocy<0||nplocy>=wrapy){
			continue;
		} else{
			nplocx1=nplocx+1;
			nplocy1=nplocy+1;
		}
		sc2=sc*cosangle(vi, vm3);
		/*sc2=sc*0.707; */
		phiout[iloc]+=sc2*(P(mapin, nplocx, nplocy)*(1.-dplocx)*(1.-dplocy)
			+P(mapin, nplocx1, nplocy)*(dplocx)*(1.-dplocy)
			+P(mapin, nplocx, nplocy1)*(1.-dplocx)*(dplocy)
			+P(mapin, nplocx1, nplocy1)*(dplocx)*(dplocy));
	}
}

/**
   Wraps proj_rect_grid for M3.

   @param[in] tsurf The tilted M3 surface map with txdeg, tydeg, ftel, fexit, fsurf set.
   @param[in,out] opd The OPD on locout, accumulated into.
   @param[in] locout The output grid.
   @param[in] thetax Angle of the guide star along x (radian).
   @param[in] thetay Angle of the guide star along y (radian).
   @param[in] hs Range of the guide star (m).
*/
void m3proj(const rmap_t* tsurf, dmat* opd, const loc_t* locout, real thetax, real thetay, real hs){
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
	proj_rect_grid(P(opd), tsurf, alx, aly, locout, scalex, scaley, NULL, scaleopd, d_img_exit, het, bx, by);
}
/**
   A convenient wrapper for m3proj() to be called from matlab or python

   @param[in,out] mapin_0 The surface map. Its keywords are replaced by keywords.
   @param[in] keywords The keywords describing the surface (txdeg, tydeg, ftel, fexit, fsurf, etc.).
   @param[in] locout The output grid.
   @param[in] thetax Angle of the guide star along x (radian).
   @param[in] thetay Angle of the guide star along y (radian).
   @param[in] hs Range of the guide star (m).
   @return Newly allocated OPD on locout; the caller frees it.
*/
dmat* m3proj2(dmat* mapin_0, const char* keywords, const loc_t* locout, real thetax, real thetay, real hs){
	free(mapin_0->keywords);
	mapin_0->keywords=strdup(keywords);
	rmap_t* mapin=rmap_convert(dref(mapin_0));
	dmat* opd=dnew(locout->nloc, 1);
	m3proj(mapin, opd, locout, thetax, thetay, hs);
	cellfree(mapin);
	return opd;
}

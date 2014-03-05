/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#include "../sys/sys.h"
#include "../math/mathdef.h"
#include "loc.h"
#include "mkdtf.h"
#include "mkh.h"

void mkdtf(ccell **pnominal, /**<[out] to be multiplied to the OTF*/
	   spcell **psi,     /**<[out] to be applied after IFFT of the final OTF*/
	   int ncompx,       /**<[in] size of OTF FFT*/
	   int ncompy,       /**<[in] size of OTF FFT*/
	   double dtheta,    /**<[in] sampling of PSF*/
	   int pixpsax,      /**<[in] number of pixels along x dimension*/
	   int pixpsay,      /**<[in] number of pixels along y dimension*/
	   double pixthetax, /**<[in] size of pixel along x dimension*/
	   double pixthetay, /**<[in] size of pixel along y dimension*/
	   double pixoffx,   /**<[in] offset of the image from the center of detector.*/
	   double pixoffy,   /**<[in] offset of the image from the center of detector*/
	   double blurx,     /**<[in] blurring as a percentage of pixel*/
	   double blury,     /**<[in] blurring as a percentage of pixel*/
	   dmat* theta       /**<[in] angle of rotation of each subaps for polar ccd. NULL for  geometry.*/
	   ){
    const int ncompx2=ncompx>>1;
    const int ncompy2=ncompy>>1;
    const double pxo=-(pixpsax*0.5-0.5+pixoffx)*pixthetax;
    const double pyo=-(pixpsay*0.5-0.5+pixoffy)*pixthetay;
    const double dux=1./(dtheta*ncompx);
    const double duy=1./(dtheta*ncompy);
    const double dux2=dux*dux;
    const double duy2=duy*duy;
    const double pdtheta=pixthetax*pixthetay/(dtheta*dtheta);
    const double duxp=dux*pixthetax;
    const double duyp=duy*pixthetay;
    const int do_blur=fabs(blurx)>EPS && fabs(blury)>EPS;
    const int nsa=theta?theta->nx:1;
    if(*pnominal){
	ccellfree(*pnominal);
    }
    *pnominal=ccellnew(nsa,1);
    if(*psi){
	spcellfree(*psi);
    }
    *psi=spcellnew(nsa,1);
    cmat *nominal=cnew(ncompx, ncompy);
    cfft2plan(nominal,-1);	
    cfft2plan(nominal,1);
    PCMAT(nominal,pn);
    loc_t *loc_psf=mksqloc(ncompx, ncompy, dtheta, dtheta, -ncompx2*dtheta, -ncompy2*dtheta);
    const double e0x=-2*M_PI*M_PI*blurx*blurx;/*blurring factors */
    const double e0y=-2*M_PI*M_PI*blury*blury;
    for(int isa=0; isa<nsa; isa++){
	double angle=theta?theta->p[isa]:0;
	double ct=cos(angle);
	double st=sin(angle);
	for(int iy=0; iy<ncompy; iy++){
	    int jy=iy-ncompy2;
	    for(int ix=0; ix<ncompx; ix++){
		int jx=ix-ncompx2;
		double ir=ct*jx+st*jy;
		double ia=-st*jx+ct*jy;
		pn[iy][ix]=sinc(ir*duyp)*sinc(ia*duxp)*pdtheta;
		if(do_blur){
		    pn[iy][ix]*=exp(e0x*(ir*ir*duy2)+e0y*(ia*ia*dux2));
		}
	    }
	}
	/*
	  with the following treatment, the operation to get image in the detector is
	  IM=si*IFFT(otf*nominal) where otf has peak in corner.
	*/
	cfftshift(nominal);
	cfft2(nominal,-1);
	cfftshift(nominal);
	cfft2(nominal,1);
	cscale(nominal,1./(double)(nominal->nx*nominal->ny));
	loc_t *loc_ccd=mksqlocrot(pixpsax,pixpsay, pixthetax,pixthetay,pxo,pyo,angle);
	(*psi)->p[isa]=mkh(loc_psf,loc_ccd,NULL,0,0,1,0,0);
	ccp(&(*pnominal)->p[isa], nominal);
	locfree(loc_ccd);
    }//isa
    cfree(nominal);
    locfree(loc_psf);
}

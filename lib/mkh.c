/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "mkh.h"

/**
   Create ray tracing operator from coordinate locin to locout.  Locin is
   required to be evenly spaced.

   If vector Pin is defined on locin, Pout is defined on locout, H=mkh(locin,
   locout, ...), Pout=H*Pin does the bilinear interpolation.

   If  locin->iac is non zero, will call mkh_cubic to produce a cubical interpolation.

   A cubic influence function that can reproduce piston/tip/tilt is coined by
   Ellerbroek to model the piezostack DM actuator. The influence has the
   form\f[
   h(x;x_{i};\delta)=h_{0}((x-x_{i})/\delta)
   \f]
   where \f$\delta\f$ is the same as the grid spacing, and \f$h_{0}\f$ is the influence
   function defined in the normalized coordinates\f{eqnarray*}{
   h_{0}(x)=\frac{1}{1+2c}\begin{cases}
   1+(4c-\frac{5}{2})|x|^{2}+(\frac{3}{2}-3c)|x|^{3} & |x|\leq1\\
   (2c-\frac{1}{2})(2-|x|)^{2}+(\frac{1}{2}-c)(2-|x|)^{3} & 1<|x|\leq2\\ 0 &
   |x|>2\end{cases}\f} where c is the nearest neighbor coupling
   frequency. The leading coefficient is to normalize the influence function so
   that it sums to 1.

 */
dsp* mkh(const loc_t* locin, const loc_t* locout,
	real displacex, real displacey, real scale){
	dsp* Hb=mkhb(locin, locout, displacex, displacey, scale);
	dsp* H=dsptrans(Hb);
	dspfree(Hb);
	return H;
}
/**
   Create transpose of mkh() result.
*/
dsp* mkhb(const loc_t* locin, const loc_t* locout,
	real displacex, real displacey, real scale){
	if(locin->iac){
		return mkhb_cubic(locin, locout, displacex, displacey, scale, locin->iac);
	}
	loc_create_map(locin);
	dsp* hback;
	int nplocx, nplocy;
	long iloc;
	/*int missing=0; */
	const real dx_in1=1./locin->dx;
	const real dx_in2=scale*dx_in1;
	const real dy_in1=1./locin->dy;
	const real dy_in2=scale*dy_in1;
	displacex=(displacex-locin->map->ox)*dx_in1;
	displacey=(displacey-locin->map->oy)*dy_in1;
	const real* px=locout->locx;
	const real* py=locout->locy;
	/*-1 because we count from 1 in the map. */
	map_t* map=locin->map;
	//const int nxmin=locin->npad;
	//const int nymin=locin->npad;
	//const int nxmax=map->nx-nxmin-1;
	//const int nymax=map->ny-nxmin-1;
	/*transpose of hfor */
	long nzmax=locout->nloc*4;
	hback=dspnew(locin->nloc, locout->nloc, nzmax);
	spint* bp=hback->pp;
	spint* bi=hback->pi;
	real* bx=hback->px;
	long count=0;
	real fx[2], fy[2];
	/*real *phiin0=phiin-1; */
	for(iloc=0; iloc<locout->nloc; iloc++){
		bp[iloc]=count;/*column index */
		if(count+5>nzmax){
			nzmax*=2;
			dspsetnzmax(hback, nzmax);
			bp=hback->pp;
			bi=hback->pi;
			bx=hback->px;
		}
		fx[1]=myfma(px[iloc], dx_in2, displacex);
		fy[1]=myfma(py[iloc], dy_in2, displacey);
		SPLIT(fx[1], fx[1], nplocx);
		SPLIT(fy[1], fy[1], nplocy);
		/*Limit the point to within active region*/
		fx[0]=1.-fx[1];
		fy[0]=1.-fy[1];
		real wtsum=0;
		for(int iy=0; iy<2; iy++){
			for(int ix=0; ix<2; ix++){
				real weight=fx[ix]*fy[iy];
				long iphi;
				/*The test on weight fixes the right/top boundary defect*/
				if(weight>EPS&&(iphi=labs(loc_map_get(map, nplocx+ix, nplocy+iy)))){
					int ic;//look for duplicates (happens when extended)
					for(ic=bp[iloc]; ic<count; ic++){
						if(bi[ic]+1==iphi){
							bx[ic]+=weight;
							break;
						}
					}
					if(ic==count){
						bi[count]=iphi-1;
						bx[count]=weight;
						count++;
					}
					wtsum+=weight;
				}
			}
		}
		if(wtsum>EPS&&wtsum<1-EPS*10){
			wtsum=1./wtsum;
			for(long ip=bp[iloc]; ip<count; ip++){
				bx[count]*=wtsum;
			}
			//warning("Scale weight by 1+%g\n", wtsum-1);
		}
	}
	bp[locout->nloc]=count;
	dspsetnzmax(hback, count);
	dspdroptol(hback, EPS);
	return hback;
}
/**
   Transposes the result from mkhb_cubic.
 */
dsp* mkh_cubic(const loc_t* locin, const loc_t* locout,
	real displacex, real displacey, real scale, real cubic_iac){
	dsp* Hb=mkhb_cubic(locin, locout, displacex, displacey, scale, cubic_iac);
	dsp* H=dsptrans(Hb);
	dspfree(Hb);
	return H;
}
/**
   Create transpose of ray tracing operator from locin to locout using cubic
   influence function that can reproduce piston/tip/tilt.  */
dsp* mkhb_cubic(const loc_t* locin, const loc_t* locout,
	real displacex, real displacey, real scale, real cubic_iac){
	if(!locin||!locout||!locin->nloc||!locout->nloc){
		warning("mkhb_cubic input is empty.\n");
		return NULL;
	}
	dsp* hback;
	real dplocx, dplocy;
	int nplocx, nplocy;
	long iloc;
	int missing=0;
	const real dx_in1=1./locin->dx;
	const real dx_in2=scale*dx_in1;
	const real dy_in1=1./locin->dy;
	const real dy_in2=scale*dy_in1;
	loc_create_map_npad(locin, 0, 0, 0);
	map_t* map=locin->map;
	displacex=(displacex-map->ox)*dx_in1;
	displacey=(displacey-map->oy)*dy_in1;
	const real* px=locout->locx;
	const real* py=locout->locy;
	/*-1 because we count from 1 in the map. */
	/*cubic */
	real fx[4], fy[4];
	//real cubic_iac=locin->iac;
	const real cubicn=1./(1.+2.*cubic_iac);
	const real c0=1*cubicn;
	const real c1=(4*cubic_iac-2.5)*cubicn;
	const real c2=(1.5-3*cubic_iac)*cubicn;
	const real c3=(2*cubic_iac-0.5)*cubicn;
	const real c4=(0.5-cubic_iac)*cubicn;
	real dplocx0, dplocy0;
	long nzmax=locout->nloc*16;
	hback=dspnew(locin->nloc, locout->nloc, nzmax);

	spint* bp=hback->pp;
	spint* bi=hback->pi;
	real* bx=hback->px;
	long count=0;
	const int nxmin=locin->npad;/*allow guarding points*/
	const int nymin=locin->npad;
	const int nxmax=map->nx-nxmin-1;
	const int nymax=map->ny-nxmin-1;
	for(iloc=0; iloc<locout->nloc; iloc++){
		bp[iloc]=count;
		if(count+17>nzmax){
			nzmax*=2;
			dspsetnzmax(hback, nzmax);
			bp=hback->pp;
			bi=hback->pi;
			bx=hback->px;
		}

		dplocx=myfma(px[iloc], dx_in2, displacex);
		dplocy=myfma(py[iloc], dy_in2, displacey);
		//outside of rectangular bounding active region .
		if(dplocy<nymin||dplocy>nymax||dplocx<nxmin||dplocx>nxmax) continue;
		SPLIT(dplocx, dplocx, nplocx);
		SPLIT(dplocy, dplocy, nplocy);
		dplocy0=1.-dplocy;
		dplocx0=1.-dplocx;

		fx[0]=dplocx0*dplocx0*(c3+c4*dplocx0);
		fx[1]=c0+dplocx*dplocx*(c1+c2*dplocx);
		fx[2]=c0+dplocx0*dplocx0*(c1+c2*dplocx0);
		fx[3]=dplocx*dplocx*(c3+c4*dplocx);

		fy[0]=dplocy0*dplocy0*(c3+c4*dplocy0);
		fy[1]=c0+dplocy*dplocy*(c1+c2*dplocy);
		fy[2]=c0+dplocy0*dplocy0*(c1+c2*dplocy0);
		fy[3]=dplocy*dplocy*(c3+c4*dplocy);
		real wtsum=0;
		for(int iy=-1; iy<+3; iy++){
			for(int ix=-1; ix<+3; ix++){
				real weight=fx[ix+1]*fy[iy+1];
				long iphi;
				/*The test on weight fixes the right/top boundary defect*/
				if(weight>EPS&&(iphi=labs(loc_map_get(map, nplocx+ix, nplocy+iy)))>0){
					int ic;//look for duplicates
					for(ic=bp[iloc]; ic<count; ic++){
						if(bi[ic]+1==iphi){
							bx[ic]+=weight;
							break;
						}
					}
					if(ic==count){
						bi[count]=iphi-1;
						bx[count]=weight;
						count++;
					}
				}
			}
		}
		if(wtsum>EPS&&wtsum<1-EPS){
			wtsum=1./wtsum;
			for(long ip=bp[iloc]; ip<count; ip++){
				bx[count]*=wtsum;
			}
			warning("Scale weight by 1+%20g\n", wtsum-1);
		}
	}/*for */
	bp[locout->nloc]=count;
	dspsetnzmax(hback, count);
	dspdroptol(hback, EPS);
	if(missing>0){
		warning("%d points not covered by input screen\n", missing);
	}
	return hback;
}
#undef ONLY_FULL
/**
   Create a matrix to bin from coordinate xin to xout using bilinear interpolation. xin and
   xout should be 1-d arrays of coordinates. We require the coordinates to order
   incrementally monotonically, but do not require them to be evenly spaced.
 */
dsp* mkhbin1d(const dmat* xin, const dmat* xout){
	if(NY(xin)!=1||NY(xout)!=1){
		error("We require both xin and xout to be only one column\n");
	}
	int iout=0;
	dsp* hbin=dspnew(NX(xout), NX(xin), NX(xin)*2);
	int count=0;
	for(int iin=0; iin<NX(xin); iin++){
		hbin->pp[iin]=count;
		real ixin=P(xin,iin);
		while(iout+1<NX(xout)&&P(xout,iout+1)<ixin){
			iout++;//find location in xout to the left of ixin
		}
		//Changes made on Nov 9, 2018 was incorrect. Correct the test
		if((iout==0&&P(xout,iout)>=ixin)||(iout+1==NX(xout))){//outside
			hbin->pi[count]=iout;
			hbin->px[count]=1;
			count++;
		} else{/*within the area */
			real wt=(ixin-P(xout,iout))/(P(xout,iout+1)-P(xout,iout));
			hbin->pi[count]=iout;
			hbin->px[count]=1.-wt;
			count++;
			hbin->pi[count]=iout+1;
			hbin->px[count]=wt;
			count++;
		}
	}
	hbin->pp[NX(xin)]=count;
	return hbin;
}


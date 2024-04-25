/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>

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
#include "petal.h"

/**
 * Create hpetal to map p/t/t per segment to the entire aperture.
*/
static void petal_do(
		dsp **ph, 	/**<[out] If set, the interaction matrix from petal modes to OPD*/
		dmat **popd,/**<out] If set, the opd from petal mode vector*/
		long nx, 	/**<[in]Size of the phase screen grid.*/
		long ny,   	/**<[in]Size of the phase screen grid.*/
		real dx,   	/**<[in]grid spacing.*/
		long nseg,  /**<[in]Number of petals.*/
		real theta0,/**<[in]Petal gap offset angle (in radian). 0 means petal gaps aligns against y axis.*/
		const dmat *mode  /**<[in]If set, output opd with the mode vector.*/
){
	if(!ph&&!popd){
		warning("Either ph or popd should be set. No action\n");
		return;
	}
	long cx=nx/2;
	long cy=ny/2;
	long ncol=nseg*3;
	real dtheta=TWOPI/nseg;//central angle

	if(theta0<0||theta0>TWOPI){
		error("Theta0=%g should be between 0 and 2pi.\n", theta0);
	}
	theta0+=M_PI*1.5;
	real Rq=4*cy*sin(dtheta/2)/(3*dtheta);//CoG of a petal from center
	real CoG[nseg][2];//CoG of each petal
	for(int iseg=0; iseg<nseg; iseg++){
		real theta2=(iseg+0.5)*dtheta-theta0;
		CoG[iseg][0]=cx+Rq*cos(theta2);
		CoG[iseg][1]=cy+Rq*sin(theta2);
	}
	dsp *ht=NULL;
	spint *pp=NULL, *pi=NULL;
	real *px=NULL;
	if(ph){
		ht=dspnew(ncol, nx*ny, nx*ny*ncol/(nseg-1));//transpose of the final result
		pp=ht->pp;
		pi=ht->pi;
		px=ht->px;
	}
	if(popd){
		dinit(popd, nx, ny);
		if(!mode){
			error("Mode must be set in order to output OPD\n");
		}else if(PN(mode)!=nseg*3){
			error("Mode must have length of %ld (is %ld)\n", nseg*3, PN(mode));
		}
	}
	
	spint count=0;
	long icol=0;
	for(long iy=0; iy<ny; iy++){
		for(long ix=0; ix<nx; ix++){
			real theta=atan2(iy-cy, ix-cx)+theta0;//guaranteed to be between [0,2*twopi]
			long iseg=ifloor(theta/dtheta);
			while(iseg>=nseg) iseg-=nseg;
			if(ph){
				pp[icol++]=count;
				px[count]=1; 		   		    pi[count]=iseg;       count++;//piston
				px[count]=(ix-CoG[iseg][0])*dx; pi[count]=iseg+nseg;  count++;//tip
				px[count]=(iy-CoG[iseg][1])*dx; pi[count]=iseg+nseg*2;count++;//tilt
				if(count>ht->nzmax){
					error("hpetal matrix overflows. count=%ld, nzmax=%ld\n", count, ht->nzmax);
				}
			}
			if(popd){
				P(*popd, ix, iy)=P(mode, iseg)
					+P(mode, iseg+nseg)*(ix-CoG[iseg][0])*dx
					+P(mode, iseg+nseg*2)*(iy-CoG[iseg][1])*dx;
			}
		}
	}
	if(ph){
		pp[nx*ny]=count;
		dspsetnzmax(ht, count);
		if(*ph) dspfree(*ph);
		*ph=dsptrans(ht);
		dspfree(ht);
	}
}
dsp *petal_mkh(long nx, long ny, real dx, long nseg, real theta0){
	dsp *h=NULL;
	petal_do(&h, NULL, nx, ny, dx, nseg, theta0, NULL);
	return h;
}
void petal_opd(anydmat opd, real dx, long nseg, real theta0, const dmat *mode){
	petal_do(NULL, &opd.dm, NX(opd.dm), NY(opd.dm), dx, nseg, theta0, mode);
}
dmat *petal_mkopd(long nx, long ny, real dx, long nseg, real theta0, const dmat *mode){
	dmat *opd=NULL;
	petal_do(NULL, &opd, nx, ny, dx, nseg, theta0, mode);
	return opd;
}

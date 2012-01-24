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

#define PRINT_RES 0 /*Print residuals in CG interations. */
#include "pcg.h"

/**
   solve A*x=b using (preconditioned) conjugate gradient method

   ; input x is first guess if not NULL and warm==1; A is symmetric
   positive definite (SPD) matrix plus low rand terms.  x0 contains the initial
   guess. Amul is the operator that applies A to x.  */
void pcg(dcell **px,    /**<[in,out] The output vector. input for warm restart.*/
	 CGFUN Amul,    /**<[in] The function that applies y=y+A*x*alpha: (*Amul)(&y,A,x,alpha)*/
	 const void *A, /**<[in] Contain info about the The left hand side matrix A*/
	 PREFUN Mmul,   /**<[in] The precond func that applies y=M^-1*x;*/
	 const void *M, /**<[in] Contains info about preconditioner matrix M.*/
	 const dcell *b,/**<[in] The right hand side vector to solve*/ 
	 int warm,      /**<[in] Use warm restart (use the value contained in px)*/
	 int maxiter    /**<[in] Max number of iterations*/
	 ){
    
    dcell *r0=NULL, *x0=NULL, *z0=NULL;
    /*computes r0=b-A*x0 */
    dcellcp(&r0, b);
    if(!*px || !warm){/*start from zero guess. */
	x0=dcellnew2(b);
	if(!*px) dcellcp(px, x0);/*initialize the output; */
    }else{
	dcellcp(&x0, *px);
	Amul(&r0, A, x0, -1);/*r0=r0+(-1)*A*x0 */
    }
    double r0z1,r0z2;
    dcell *p0=NULL;
    if(Mmul){
	Mmul(&z0,M,r0);
    }else{
	z0=r0;
    }
    dcellcp(&p0, z0);
    r0z1=dcellinn(r0,z0);
    dcell *Ap=NULL;
    double ak,bk;

    /*double r0zmin=r0z1; */
#if PRINT_RES == 1
    double r0z0=dcellinn(b,b);/*|b| */
    double res[maxiter+1];
    if(Mmul){
	res[0]=sqrt(dcellinn(r0,r0)/r0z0);
    }else{
	res[0]=sqrt(r0z1/r0z0);
    }
    int kres=0;
    info("Step %d, res=%g\n", kres, res[kres]);
#endif
    for(int k=0; k<maxiter; k++){
	if(k!=0 && k%100==0){ /*restart every 100 steps exclude beginning */
	    /*info("Restarting at step %d\n",k); */
	    dcelladd(&r0, 0., b, 1.);/*r0=b; */
	    (*Amul)(&r0, A, x0, -1);/*r0=r0+(-1)*A*x0 */
	    if(Mmul){
		(*Mmul)(&z0,M,r0);
	    }else{
		z0=r0;
	    }
	    dcellcp(&p0, z0);
	}
	if(Ap) dcellzero(Ap);
	(*Amul)(&Ap, A, p0, 1);
	ak=r0z1/dcellinn(p0,Ap);
	dcelladd(&x0, 1, p0, ak);/*x0=x0+ak*p0 */
	dcelladd(&r0, 1, Ap, -ak);/*r0=r0-ak*Ap */
	if(Mmul){
	    (*Mmul)(&z0,M,r0);
	}else{
	    z0=r0;
	}
	r0z2=dcellinn(r0,z0);/*r0z2=r0'*r0 */
	bk=r0z2/r0z1;
	dcelladd(&p0, bk, z0, 1.);/*p0=bk*pi+r0 */

	r0z1=r0z2;
#if PRINT_RES == 1
	if(Mmul){
	    res[k+1]=sqrt(dcellinn(r0,r0)/r0z0);
	}else{
	    res[k+1]=sqrt(r0z2/r0z0);
	}
	info("Step %d, res=%g\n", k+1, res[k+1]);
#endif
    }
    dcellcp(px, x0);
    dcellfree(r0); 
    if(Mmul){
	dcellfree(z0);
    }
    dcellfree(x0);
    dcellfree(Ap);
    dcellfree(p0);
}

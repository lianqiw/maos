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
extern "C"
{
#include <cuda.h>
#include "gpu.h"
}
#include "utils.h"
#include "recon.h"
#include "curmat.h"

void cumuv(curcell **out, float beta, cumuv_t *A, const curcell *in, float alpha){
    if(!A->Mt) error("A->M Can not be empty\n");
    if(A->U && A->U->ny>1 || A->V && A->V->ny>1) error("Not handled yet\n");
    if(!*out){
	*out=curcellnew(A->Mt->ny, 1);
	int nx[A->Mt->ny];
	for(int ix=0; ix<A->Mt->ny; ix++){
	    nx[ix]=A->Mt->p[ix*A->Mt->nx]->ny;
	}
	*out=curcellnew(A->Mt->ny, 1, nx, (int*)NULL);
    }
    for(int idm=0; idm<A->Mt->ny; idm++){
	if(fabs(beta)<EPS) curzero((*out)->p[idm], A->dmstream[idm]);
	else if(fabs(beta-1)>EPS) 
	    curscale((*out)->p[idm], beta, A->dmstream[idm]);
	for(int ifit=0; ifit<A->Mt->nx; ifit++){
	    cusptmul((*out)->p[idm]->p, A->Mt->p[ifit+idm*A->Mt->nx], in->p[ifit]->p,
		     alpha, A->dmsphandle[idm]);
	}
    }
    curmat *tmp=NULL;
    if(A->U && A->V){
	tmp=curnew(A->V->p[0]->ny, 1);
	for(int ifit=0; ifit<A->V->nx; ifit++){
	    curmv(tmp->p, 1.f, A->V->p[ifit], in->p[ifit]->p, 't', 1, A->fithandle[0]);
	}
	cudaStreamSynchronize(A->fitstream[0]);
	for(int idm=0; idm<A->U->nx; idm++){
	    curmv((*out)->p[idm]->p, 1.f, A->U->p[idm], tmp->p, 'n', -alpha, A->dmhandle[idm]);
	}
    }
    for(int idm=0; idm<A->Mt->ny; idm++){
	cudaStreamSynchronize(A->dmstream[idm]);
    }
    curfree(tmp);
}
void cumuv_trans(curcell **out, float beta, cumuv_t *A, const curcell *in, float alpha){
    for(int ifit=0; ifit<A->Mt->nx; ifit++){
	if(fabs(beta)<EPS) curzero((*out)->p[ifit], A->fitstream[ifit]);
	else if(fabs(beta-1)>EPS) 
	    curscale((*out)->p[ifit], beta, A->fitstream[ifit]);
	for(int idm=0; idm<A->Mt->ny; idm++){
	    cuspmul((*out)->p[ifit]->p, A->Mt->p[ifit+idm*A->Mt->nx], in->p[idm]->p, alpha, 
		    A->fitsphandle[ifit]);
	}
    }
    curmat *tmp=NULL;
    if(A->V && A->U){
	tmp=curnew(A->U->p[0]->ny, 1);
	for(int idm=0; idm<A->U->nx; idm++){
	    curmv(tmp->p, 1.f, A->U->p[idm], in->p[idm]->p, 't', 1, A->dmhandle[0]);
	}
	cudaStreamSynchronize(A->dmstream[0]);
	for(int ifit=0; ifit<A->V->nx; ifit++){
	    curmv((*out)->p[ifit]->p, 1.f, A->V->p[ifit], tmp->p, 'n', -alpha, A->fithandle[ifit]);
	}
    }
    for(int ifit=0; ifit<A->Mt->nx; ifit++){
	cudaStreamSynchronize(A->fitstream[ifit]);
    }
    curfree(tmp);
}

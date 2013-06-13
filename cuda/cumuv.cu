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
extern "C"
{
#include <cuda.h>
#include "gpu.h"
}
#include "utils.h"
#include "recon.h"
#include "curmat.h"
/*
  This routine needs to be improved by merging all these operations in loops to a big kernel and only use stream.
*/
void cumuv(curcell **out, float beta, cumuv_t *A, const curcell *in, float alpha, stream_t &stream){
    if(!A->Mt) error("A->M Can not be empty\n");
    if(A->U && A->U->ny>1 || A->V && A->V->ny>1) error("Not implemented yet\n");
    stream.sync();//must sync because it is no used.
    const int ndm=A->Mt->ny;
    curmat *tmp=NULL;
    if(A->U && A->V){
	tmp=curnew(A->V->p[0]->ny, 1);
	//No need to sync dmstream here.
	for(int ifit=0; ifit<A->V->nx; ifit++){
	    curmv(tmp->p, 1.f, A->V->p[ifit], in->p[ifit]->p, 't', 1, A->fitstream[0]);
	}
    }
    if(!*out){
	*out=curcellnew(ndm, 1);
	int nact[ndm];
	for(int idm=0; idm<ndm; idm++){
	    nact[idm]=A->Mt->p[idm*A->Mt->nx]->ny;
	}
	*out=curcellnew(A->Mt->ny, 1, nact, (int*)NULL);
    }
    for(int idm=0; idm<ndm; idm++){
	curscale((*out)->p[idm], beta, A->dmstream[idm]);
	for(int ifit=0; ifit<A->Mt->nx; ifit++){
	    cusptmul((*out)->p[idm]->p, A->Mt->p[ifit+idm*A->Mt->nx], in->p[ifit]->p,
		     alpha, A->dmstream[idm]);
	}
    }
    if(tmp){
	A->fitstream[0].sync();
	for(int idm=0; idm<ndm; idm++){
	    curmv((*out)->p[idm]->p, 1.f, A->U->p[idm], tmp->p, 'n', -alpha, A->dmstream[idm]);
	}
    }
    for(int idm=0; idm<ndm; idm++){
	A->dmstream[idm].sync();
    }
    curfree(tmp);
}
void cumuv_trans(curcell **out, float beta, cumuv_t *A, const curcell *in, float alpha, stream_t &stream){
    stream.sync();//must sync.
    for(int ifit=0; ifit<A->Mt->nx; ifit++){
	if(fabs(beta)<EPS) curzero((*out)->p[ifit], A->fitstream[ifit]);
	else if(fabs(beta-1)>EPS) 
	    curscale((*out)->p[ifit], beta, A->fitstream[ifit]);
	for(int idm=0; idm<A->Mt->ny; idm++){
	    cuspmul((*out)->p[ifit]->p, A->Mt->p[ifit+idm*A->Mt->nx], in->p[idm]->p, alpha, 
		    A->fitstream[ifit]);
	}
    }
    curmat *tmp=NULL;
    if(A->V && A->U){
	tmp=curnew(A->U->p[0]->ny, 1);
	for(int idm=0; idm<A->U->nx; idm++){
	    curmv(tmp->p, 1.f, A->U->p[idm], in->p[idm]->p, 't', 1, A->dmstream[0]);
	}
	cudaStreamSynchronize(A->dmstream[0]);
	for(int ifit=0; ifit<A->V->nx; ifit++){
	    curmv((*out)->p[ifit]->p, 1.f, A->V->p[ifit], tmp->p, 'n', -alpha, A->fitstream[ifit]);
	}
    }
    for(int ifit=0; ifit<A->Mt->nx; ifit++){
	cudaStreamSynchronize(A->fitstream[ifit]);
    }
    curfree(tmp);
}

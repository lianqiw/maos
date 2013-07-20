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
    if(!A->M) error("A->M Can not be empty\n");
    if(!*out){
	*out=curcellnew(A->nx, 1, A->nxs, (int*)NULL);
    }else{
	curscale((*out)->m, beta, stream);
    }
    cuspmul((*out)->m->p, A->M, in->m->p, 1, 'n', alpha, stream);
    if(A->U && A->V){
	curmv(A->Vx->p, 0, A->V, in->m->p, 't', 1, stream);
	curmv((*out)->m->p, 1, A->U, A->Vx->p, 'n', -alpha, stream);
    }
}
void cumuv_trans(curcell **out, float beta, cumuv_t *A, const curcell *in, float alpha, stream_t &stream){ 
    if(!A->M) error("A->M Can not be empty\n");
    if(!*out){
	*out=curcellnew(A->ny, 1, A->nys, (int*)NULL);
    }else{
	curscale((*out)->m, beta, stream);
    }
    
    curscale((*out)->m, beta, stream);
    cuspmul((*out)->m->p, A->M, in->m->p, 1, 't', alpha, stream);
    if(A->U && A->V){
	curmv(A->Vx->p, 0, A->U, in->m->p, 't', 1, stream);
	curmv((*out)->m->p, 1, A->V, A->Vx->p, 'n', -alpha, stream);
    }
}

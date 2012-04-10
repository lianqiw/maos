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
#include "recon.h"

/*solve in place*/
static __global__ void cuchol_solve_lower_do(float *restrict y, float *Cx, int *Cp, int *Ci, int n){
    int id=threadIdx.x;
    int nd=blockDim.x;
    __shared__ float val;
    __shared__ float sum;
    /*first solve L\y*/
    
    for(int icol=0; icol<n; icol++){
	if(id==0){
	    y[icol]/=Cx[Cp[icol]];//divide by diagonal.
	    val=-y[icol];
	}
	__syncthreads();//this is necessary!
	for(int irow=Cp[icol]+1+id; irow<Cp[icol+1]; irow+=nd){
	    y[Ci[irow]]+=val*Cx[irow];
	}
	__syncthreads();
    }
    /*Next solve L'\y*/
    if(id==0) sum=0;
    __syncthreads();
    for(int icol=n-1; icol>-1; icol--){
	for(int irow=Cp[icol]+1+id; irow<Cp[icol+1]; irow+=nd){
	    atomicAdd(&sum, Cx[irow]*y[Ci[irow]]);
	}
	__syncthreads();
	if(id==0){
	    y[icol]=(y[icol]-sum)/Cx[Cp[icol]];
	    sum=0;
	}
	__syncthreads();//this is necessary!
    }
}
/**
   Solve cholesky backsubstitution in GPU. Shared memory is not enough to contain the vector.
*/
void cuchol_solve(float *restrict out, cusp *Cl, int *Cp, const float *restrict in, 
		  cudaStream_t stream){
    if(!Cl || !Cp) error("Invalid\n");
    int n=Cl->nx;
    float *y;
    DO(cudaMalloc(&y, sizeof(float)*n));
    perm_f_do<<<DIM(n, 256),0,stream>>>(y, in, Cp, n);
    cuchol_solve_lower_do<<<1,256,0,stream>>>(y, Cl->x, Cl->p, Cl->i, n); //only 1 block for synchronization.
    perm_i_do<<<DIM(n, 256),0,stream>>>(out, y, Cp, n);
    cudaStreamSynchronize(stream);
    cudaFree(y);
}

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
#include "../lib/aos.h"
#include "gpu.h"
}
#include "kernel.h"
#include "utils.h"
#include "curmat.h"
/**
   \file test.cu
   Routines for testing cuda. Execute from test/test_cuda
   
*/

void test_sum(){
    TIC;tic;
    const int BS=128;
    int N=120*120*6000;
    smat *As=snew(N,1);
    //rand_t srand;
    //seed_rand(&srand, 1);
    //srandu(As, 1, &srand);
    sset(As, 1);
    curmat *Ap;
    cp2gpu(&Ap, As);
    float *res_gpu;
    cudaMalloc(&res_gpu, sizeof(float));
    cudaMemset(res_gpu, 0, sizeof(float));
    float res[4];
    cudaStream_t stream;
    STREAM_NEW(stream);
    toc("malloc");tic;
    sum_wrap(res_gpu, Ap->p, N, stream);
    cudaMemcpyAsync(res, res_gpu, sizeof(float), cudaMemcpyDeviceToHost);
    CUDA_SYNC_STREAM;
    toc("sum_wrap");tic;
    sum_do<<<DIM_REDUCE, BS, DIM_REDUCE*sizeof(float), stream>>>(res_gpu, Ap->p, N);
    cudaMemcpyAsync(res+1, res_gpu, sizeof(float), cudaMemcpyDeviceToHost);
    CUDA_SYNC_STREAM;
    double tim=toc3*1024*1024*1024;
    toc("sum_wrap");tic;
    sum2_wrap(res_gpu, Ap->p, N, stream);
    cudaMemcpyAsync(res+2, res_gpu, sizeof(float), cudaMemcpyDeviceToHost);
    CUDA_SYNC_STREAM;
    toc("sum2_wrap");tic;
    sum2_do<<<DIM_REDUCE, BS, 0, stream>>>(res_gpu, Ap->p, N);
    //sum2_wrap(res_gpu, Ap->p, N, stream);
    cudaMemcpyAsync(res+3, res_gpu, sizeof(float), cudaMemcpyDeviceToHost);
    CUDA_SYNC_STREAM;
    double tim2=toc3*1024*1024*1024;
    toc("sum2_wrap");

    info("sum_wrap  %.2f GB/s\n", N*sizeof(float)/tim);
    info("sum2_wrap %.2f GB/s\n", N*sizeof(float)/tim2);
    info("Result: %g %g %g %g\n", res[0], res[1]-res[0], res[2]-res[1], res[3]-res[2]);
}
void test_gpu(){
    test_sum();
}

/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "../math/cumath.h"
#include "../math/prop_map.h"

/**
   \file test.cu
   Routines for testing cuda. Execute from test/test_cuda

*/

void test_sum(){
	TIC;tic;
	const int BS=128;
	int N=120*120*6000;
	X(mat)* As=X(new)(N, 1);
	//rand_t srand;
	//seed_rand(&srand, 1);
	//srandu(As, 1, &srand);
	X(set)(As, 1);
	curmat Ap;
	cp2gpu(Ap, As);
	Real* res_gpu;
	cudaMalloc(&res_gpu, sizeof(Real));
	cudaMemset(res_gpu, 0, sizeof(Real));
	Real res[4];
	cudaStream_t stream;
	STREAM_NEW(stream);
	toc2("malloc");tic;
	sum_wrap(res_gpu, Ap, N, stream);
	DO(cudaMemcpyAsync(res, res_gpu, sizeof(Real), D2H, stream));
	CUDA_SYNC_STREAM;
	toc2("sum_wrap");tic;
	sum_do<<<DIM_REDUCE, BS, DIM_REDUCE*sizeof(Real), stream>>>(res_gpu, Ap, N);
	DO(cudaMemcpyAsync(res+1, res_gpu, sizeof(Real), D2H,stream));
	CUDA_SYNC_STREAM;
	real tim=toc3*1024*1024*1024;
	toc2("sum_wrap");tic;
	sum2_wrap(res_gpu, Ap, N, stream);
	DO(cudaMemcpyAsync(res+2, res_gpu, sizeof(Real), D2H, stream));
	CUDA_SYNC_STREAM;
	toc2("sum2_wrap");tic;
	sum2_do<<<DIM_REDUCE, BS, 0, stream>>>(res_gpu, Ap, N);
	//sum2_wrap(res_gpu, Ap, N, stream);
	DO(cudaMemcpyAsync(res+3, res_gpu, sizeof(Real), D2H, stream));
	CUDA_SYNC_STREAM;
	real tim2=toc3*1024*1024*1024;
	toc2("sum2_wrap");

	dbg("sum_wrap  %.2f GB/s\n", N*sizeof(Real)/tim);
	dbg("sum2_wrap %.2f GB/s\n", N*sizeof(Real)/tim2);
	dbg("Result: %g %g %g %g\n", res[0], res[1]-res[0], res[2]-res[1], res[3]-res[2]);
}
void gpu_map2map(cumap_t& out, const cumap_t& in, Real dispx, Real dispy, Real alpha, const curmat& cc, char trans){
	map2map_t wrap;
	Array<map2map_t, Gpu> wrap_gpu(1);
	//cudaMalloc(&wrap_gpu, sizeof(PROP_WRAP_T));
	map2map_prep(&wrap, out, in, dispx, dispy, cc);
	wrap.togpu(wrap_gpu());
	Real** p;
	cudaMalloc(&p, sizeof(Real*)*2);
	const Real* tmp[2]={out(), in()};
	DO(cudaMemcpy(p, tmp, sizeof(Real*)*2, H2D));
	map2map_do<<<dim3(4, 4, 1), dim3(WRAP_TX, 4), 0, 0>>>
		(wrap_gpu, p, p+1, 0, 1, 1, alpha, 0, trans);
	DO(cudaMemcpy(&wrap, wrap_gpu(), sizeof(map2map_t), D2H));
	if(wrap.reverse){
		cudaFree(wrap.reverse);
	}
}

/*Test ray tracing*/
void test_prop(){
	cudaSetDevice(0);
	map_t* mapin=mapnew(30, 30, 0.5, 0.5);
	map_t* mapout=mapnew(200, 200, 0.1, 0.1);
	rand_t rstat;
	seed_rand(&rstat, 1);
	drandn(mapin->dmat, 1, &rstat);
	mapin->iac=0.3;
	writebin(mapin, "prop_mapin_cpu");
	cumap_t cumapin; cumapin=mapin;
	cumap_t cumapout; cumapout=mapout;
	cp2gpu(cumapin.p, mapin->dmat);
	cuwrite(cumapin.p, 0, "prop_mapin");
	Real alpha=1;
	Real dispx=0;
	Real dispy=0;
	curmat cc=iac2cc(mapin->iac);
	gpu_map2map(cumapout, cumapin, dispx, dispy, alpha, cc, 'n');
	cuwrite(cumapout.p, 0, "prop_mapout");
	propdata_t propdata={.mapin=mapin, .mapout=mapout, .alpha=alpha, .shiftx=dispx, .shifty=dispy};
	prop(&propdata);
	writebin(mapout, "prop_mapout_cpu");
}

int main(int argc, char** argv){
	const char* cmd=0;
	if(argc>1){
		cmd=argv[1];
	}
	if(!cmd){
		test_prop();
	} else if(!strcmp(cmd, "sum")){
		test_sum();
	} else if(!strcmp(cmd, "prop")){
		test_prop();
	}
}

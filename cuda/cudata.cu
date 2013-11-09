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
#include "cudata.h"
int nstream=0;

int gpu_recon;/**<GPU for reconstruction*/
int NGPU=0;
int* GPUS=NULL;
cudata_t *cudata_all=NULL;/*for all GPU. */

#ifdef __APPLE__
pthread_key_t cudata_key;
static __attribute((constructor)) void init(){
    pthread_key_create(&cudata_key, NULL);
}
#else
__thread cudata_t *cudata=NULL;/*for current thread and current GPU */
#endif

/**
   Get GPU info.
*/
void gpu_info(){
    struct cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    info("name=%s\n"
	 "TotalGlobalMem=%d\n"
	 "SharedMemPerBlock=%d\n"
	 "regsPerBlock=%d\n"
	 "warpSize=%d",
	 prop.name,
	 (int)prop.totalGlobalMem,
	 (int)prop.sharedMemPerBlock,
	 prop.regsPerBlock,
	 prop.warpSize);
}
/**
   Print memory consumption.
*/
void gpu_print_mem(const char *msg){
    size_t fr, tot;
    cudaDeviceSynchronize();
    DO(cudaMemGetInfo(&fr, &tot));
    info2("GPU (%d) mem used %ld MB (%s)\n",(int)(cudata-cudata_all),(long)(tot-fr)/1024/1024, msg);
}
/**
   Get available memory.
*/
long gpu_get_mem(void){
    size_t fr, tot;
    DO(cudaMemGetInfo(&fr, &tot));
    return (long)fr;
}
/**
   Get available memory.
*/
static long gpu_get_idle_mem(void){
    size_t fr, tot;
    DO(cudaMemGetInfo(&fr, &tot));
    if(tot-fr>500000000){//GPU used by some other process. do not use it.
	return 0;
    }else{
	return (long)fr;
    }
}
static int cmp_gpu_info(const long *a, const long *b){
    return b[1]>a[1]?1:0;
}
/**
   Initialize GPU. Return 1 if success.
   if gpus is not null, it is of length ngpu. gpus specifies gpu index to use.
   if gpus is null, ngpu specifies number of gpus to use. all if 0.

   when mix Tesla and GTX cards, the ordering the GPUs may be different in CUDA
   and NVML, causing the selection of GPUs to fail. Do not use NVML 
*/
int gpu_init(int *gpus, int ngpu){
    int ans, ngpu_tot=0;//total number of GPUs.
    if((ans=cudaGetDeviceCount(&ngpu_tot)) || ngpu_tot==0){//no GPUs available.
	info2("No GPUs available. ans=%d\n", ans);
	return 0;
    }
    //sem_lock("/maos.gpu");
    /*do not use sem_lock. if program terminators before unlock. subsequent programs cannot proceed.*/
    char fnlock[PATH_MAX];
    snprintf(fnlock, PATH_MAX, "%s/gpu.lock", TEMP);
    int fdlock=lock_file(fnlock, 1, 0);
    NGPU=0;
    /*
      User specified exact GPUs to use. We check every entry. 
      If <0 is found, do not use any GPU.
      If >=ngpu_tot is found, skip the GPU.
      If duplicates are found, use only once.
     */
    if(gpus && ngpu>0){
	if(!GPUS) GPUS=(int*)malloc(ngpu*sizeof(int));
	for(int ig=0; ig<ngpu; ig++){
	    if(gpus[ig]<0){
		info2("CUDA is disabled by user.\n");
		free(GPUS); GPUS=NULL; 
		NGPU=0;
		goto end;
	    }else{
		if(gpus[ig]>=ngpu_tot){
		    warning2("Skip GPU %d: not exist\n", gpus[ig]);
		}else{
		    GPUS[NGPU++]=gpus[ig];
		    /* Enable the following to disallow use GPUs in multiple threads
		      int j;
		    for(j=0; j<NGPU; j++){
			if(GPUS[j]==gpus[ig]){
			    warning2("Skip GPU %d: duplicated\n", gpus[ig]);
			    break;
			}
		    }
		    if(j==NGPU){
			GPUS[NGPU++]=gpus[ig];
			}*/
		}
	    }
	}
    }else{
	int repeat=0;
	if(ngpu<=0){
	    repeat=0;
	    ngpu=ngpu_tot;
	}
	GPUS=(int*)calloc(ngpu, sizeof(int));
	register_deinit(NULL, GPUS);
	/*For each GPU, query the available memory.*/
	long (*gpu_info)[2]=(long(*)[2])calloc(2*ngpu_tot, sizeof(long));
	int gpu_valid_count;
	do{
	    for(int ig=0; ig<ngpu_tot; ig++){
		gpu_info[ig][0]=ig;
		cudaSetDevice(ig);//this allocates context.
		gpu_info[ig][1]=gpu_get_idle_mem();
	    }
	    /*sort so that gpus with higest memory is in the front.*/
	    qsort(gpu_info, ngpu_tot, sizeof(long)*2, (int(*)(const void*, const void *))cmp_gpu_info);
	    gpu_valid_count=0;
	    for(int igpu=0; igpu<ngpu_tot; igpu++){
		if(gpu_info[igpu][1]>=500000000){
		    gpu_valid_count++;
		}
		info2("GPU %d has mem %.1f GB\n", (int)gpu_info[igpu][0], gpu_info[igpu][1]/1024/1024/1024.);
	    }
	}while(gpu_valid_count<ngpu && gpu_valid_count<ngpu_tot && sleep(60));

	for(int i=0, igpu=0; i<ngpu; i++, igpu++){
	    if(igpu==ngpu_tot || gpu_info[igpu][1]<500000000){
		if(repeat){
		    igpu=0; //reset to beginning.
		}else{
		    break; //stop
		}
	    }
	    GPUS[NGPU++]=(int)gpu_info[igpu][0];
	}
	free(gpu_info);
    }
    if(NGPU) {
	gpu_recon=0;/*first gpu in GPUS*/
	cudata_all=new cudata_t[NGPU];
	register_deinit(NULL, cudata_all);
	info2("Using GPU");
	for(int i=0; GPUS && i<NGPU; i++){
	    info2(" %d", GPUS[i]);
	    gpu_set(i);
	    //Reserve memory in GPU so the next maos will not pick this GPU.
	    DO(cudaMalloc(&cudata->reserve, 500000000));
	}
	info2("\n");
    }
 end:
    //sem_unlock("maos.gpu");
    close(fdlock);
    return NGPU;
}

/**
   Clean up device.
*/
void gpu_cleanup(void){
    for(int ig=0; ig<NGPU; ig++){
	cudaSetDevice(GPUS[ig]);
	cudaDeviceReset();
    }
}

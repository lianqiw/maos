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
#define MEM_RESERVE 200000000
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
int cudata_t::recongpu=0;
int *cudata_t::evlgpu=0;
int *cudata_t::wfsgpu=0;
cuwfs_t *cudata_t::wfs=0;
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
    info2("GPU (%d) mem used %ld MB (%s)\n",cudata->igpu,(long)(tot-fr)/1024/1024, msg);
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
static long gpu_get_free_mem(int igpu){
    size_t fr=0, tot=0;
    int ans;
    if((ans=cudaMemGetInfo(&fr, &tot))){
	warning2("cudaMemGetInfo failed with error %d\n", ans);
    }
    info2("GPU%2d has %.1fGB free, %.1fGB total device memory.\n", 
	  igpu, fr*9.3e-10, tot*9.3e-10);
    return (long)fr;
}
static int cmp_long2_descend(const long *a, const long *b){
    if(b[1]>a[1]){
	return 1;
    }else if(a[1]>b[1]){
	return -1;
    }else{
	return 0;
    }
}
static int cmp_long2_ascend(const long *a, const long *b){
    if(b[1]<a[1]){
	return 1;
    }else if(a[1]<b[1]){
	return -1;
    }else{
	return 0;
    }
}
struct task_t{
    float timing;/*based on empirical data*/
    int *dest;
};
static int task_cmp(const task_t *a, const task_t *b){
    if(b->timing > a->timing){
	return 1;
    }else if(b->timing < a->timing){
	return -1;
    }else{
	return 0;
    }
}
/**
   Initialize GPU. Return 1 if success.
   if gpus is not null, it is of length ngpu. gpus specifies gpu index to use.
   if gpus is null, ngpu specifies number of gpus to use. all if 0.

   when mix Tesla and GTX cards, the ordering the GPUs may be different in CUDA
   and NVML, causing the selection of GPUs to fail. Do not use NVML 
*/
int gpu_init(int *gpus, int ngpu, const PARMS_T *parms){
    if(gpus && ngpu>0){
	for(int ig=0; ig<ngpu; ig++){
	    if(gpus[ig]<0){
		info2("CUDA is disabled by user.\n");
		return 0;
	    }
	}
    }
    int ans, ngpu_tot=0;//total number of GPUs.
    if((ans=cudaGetDeviceCount(&ngpu_tot)) || ngpu_tot==0){//no GPUs available.
	info2("No GPUs available. ans=%d\n", ans);
	return 0;
    }
    long mem_minimum=0;
    if(parms){
	if(parms->gpu.evl || parms->gpu.wfs){
	    const int nps=parms->atm.nps;
	    long nxn=parms->atm.nxn;
	    long nyn=parms->atm.nyn;
	    mem_minimum+=sizeof(Real)*nps*nxn*nyn;
	}
	if(parms->gpu.evl){
	    mem_minimum+=sizeof(Real)*parms->evl.nevl*pow(parms->aper.d/parms->evl.dx, 2);
	}
	if(parms->gpu.wfs){
	    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		mem_minimum+=sizeof(Real)*parms->powfs[ipowfs].nwfs*pow(parms->aper.d/parms->powfs[ipowfs].dx, 2)*4;
	    }
	}
	if(parms->gpu.tomo || parms->gpu.fit){
	    mem_minimum+=sizeof(Real)*parms->atmr.nps*pow(parms->aper.d*parms->tomo.pos/parms->atmr.dx, 2)*4;
	}
	if(mem_minimum==0){//gpu is disabled
	    return 0;
	}else{
	    info2("CUDA: minimum memory requirement is %.1fGB\n", mem_minimum/(double)(1024*1024*1024));
	}
    }
    char fnlock[PATH_MAX];
    snprintf(fnlock, PATH_MAX, "%s/gpu.lock", TEMP);
    int fdlock=lock_file(fnlock, 1, 0);
    /*
      Create a mapping between CUDA device ordering and NVML ordering (nvidia-smi).
    */
    long gmap[ngpu_tot][2];
    for(int ig=0; ig<ngpu_tot; ig++){
	gmap[ig][0]=ig;
	cudaDeviceProp properties;
	if(!cudaSetDevice(ig) && !cudaGetDeviceProperties(&properties, ig)){
	    gmap[ig][1]=properties.pciBusID;
	}else{
	    error("Error getting information for GPU %d\n", ig);
	}
    }
    qsort(gmap, ngpu_tot, sizeof(long)*2, (int(*)(const void*, const void *))cmp_long2_ascend);
    //now gmap[igpu] is the cuda index of the nvml index igpu.
    NGPU=0;
    /*
      User specified exact GPUs to use. We check every entry. 
      If <0 is found, do not use any GPU.
      If >=ngpu_tot is found, skip the GPU and print warning.
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
		    warning("GPU %d: not exist\n", gpus[ig]);
		}else{
		    GPUS[NGPU++]=gmap[gpus[ig]][0];
		}
	    }
	}
    }else{
	int repeat=0;
	if(ngpu<=0){
	    repeat=0;
	    ngpu=ngpu_tot;
	}
	GPUS=(int*)calloc(ngpu, sizeof(int));//stores CUDA index
	register_deinit(NULL, GPUS);
	/*For each GPU, query the available memory.*/
	long (*gpu_info)[2]=(long(*)[2])calloc(2*ngpu_tot, sizeof(long));
	int gpu_valid_count;
	do{
	    gpu_valid_count=0;
	    for(int jg=0; jg<ngpu_tot; jg++){//jg: nvml index. ig: cuda index
		int ig=gmap[jg][0];
		gpu_info[ig][0]=ig;
		if(!cudaSetDevice(ig)){
		    //this allocates context and create a CPU thread for this GPU.
		    gpu_info[ig][1]=gpu_get_free_mem(ig);
		    if(gpu_info[ig][1]>=mem_minimum){
			gpu_valid_count++;
		    }
		}
	    }
	}while(0);
	if(gpu_valid_count){
	    /*sort so that gpus with higest memory is in the front.*/
	    qsort(gpu_info, ngpu_tot, sizeof(long)*2, (int(*)(const void*, const void *))cmp_long2_descend);
	    for(int i=0, ig=0; i<ngpu; i++, ig++){//ig: cuda index
		if(ig<ngpu_tot && gpu_info[ig][1]>=mem_minimum){
		    GPUS[NGPU++]=(int)gpu_info[ig][0];
		}else if(ig==ngpu_tot || gpu_info[ig][1]<mem_minimum){
		    if(NGPU && repeat){
			ig=0; //reset to beginning.
		    }else{
			break; //stop
		    }
		}
	    }
	}
	free(gpu_info);
    }
    if(NGPU) {
	cudata_all=new cudata_t[NGPU];
	register_deinit(NULL, cudata_all);
	info2("Using GPU");
	for(int i=0; GPUS && i<NGPU; i++){
	    gpu_set(i);
	    for(int j=0; j<NGPU; j++){
		if(j!=i){
		    cudaDeviceEnablePeerAccess(j, 0);
		}
	    }
	    for(int j=0; j<ngpu_tot; j++){
		if(GPUS[i]==gmap[j][0]){
		    cudata->igpu=j;
		    break;
		}
	    }
	    info2(" %d", cudata->igpu);
	    //Reserve memory in GPU so the next maos will not pick this GPU.
	    DO(cudaMalloc(&cudata->reserve, MEM_RESERVE));
	}
	info2("\n");
	if(parms){
	    /*Assign task to gpu evenly based on empirical data to maximum GPU
	     * usage. We first gather together all the tasks and assign a timing
	     * in ms to each. Sort all the tasks in descend order and then
	     * iteratively assign each task to the minimally used GPU*/
	    cudata_t::evlgpu=(int*)calloc(parms->evl.nevl, sizeof(int));
	    cudata_t::wfsgpu=(int*)calloc(parms->nwfs, sizeof(int));
	    int ntask=0;
	    if(parms->gpu.tomo || parms->gpu.fit) ntask++;
	    if(parms->gpu.evl) ntask+=parms->evl.nevl;
	    if(parms->gpu.wfs) ntask+=parms->nwfs;
	    if(ntask==0){
		delete [] cudata_all;
		free(GPUS);
		NGPU=0;GPUS=0;
		return NGPU;
	    }
	    struct task_t *tasks=(task_t*)calloc(ntask, sizeof(task_t));
	    //recon
	    int count=0;
	    if(parms->gpu.tomo || parms->gpu.fit){
		tasks[count].timing=20;//ms
		tasks[count].dest=&cudata_t::recongpu;
		count++;
	    }
	    //evl
	    for(int ievl=0; parms->gpu.evl && ievl<parms->evl.nevl; ievl++){
		tasks[count].timing=4.7;//ms
		tasks[count].dest=cudata_t::evlgpu+ievl;
		count++;
	    }
	    //wfs
	    for(int iwfs=0; parms->gpu.wfs && iwfs<parms->nwfs; iwfs++){
		const int ipowfs=parms->wfs[iwfs].powfs;
		tasks[count].timing=parms->powfs[ipowfs].usephy?17:1.5;
		tasks[count].dest=cudata_t::wfsgpu+iwfs;
		count++;
	    }
	    qsort(tasks, count, sizeof(task_t), (int(*)(const void*, const void *))task_cmp);
	    float timtot[NGPU];
	    for(int igpu=0; igpu<NGPU; igpu++){
		timtot[igpu]=0;
	    }
	    for(int it=0; it<count; it++){
		int min_gpu=0; /*current gpu with minimum task*/
		float min_val=INFINITY;
		for(int igpu=0; igpu<NGPU; igpu++){
		    if(timtot[igpu]<min_val){
			min_val=timtot[igpu];
			min_gpu=igpu;
		    }
		}
		*(tasks[it].dest)=min_gpu;
		timtot[min_gpu]+=tasks[it].timing;
	    }
	    if(parms->sim.nthread>NGPU && (parms->gpu.tomo || parms->gpu.fit) && parms->gpu.evl && parms->gpu.wfs){
		((PARMS_T*)parms)->sim.nthread=NGPU+1;
		NTHREAD=NGPU+1;
		info2("Reset nthread to %d\n", NTHREAD);
	    }
	    free(tasks);
	}
    }
 end:
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

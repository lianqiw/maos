/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#include <unistd.h>
#include "cudata.h"
int nstream=0;
#define MEM_RESERVE 50000000 
int NGPU=0;
int MAXGPU=0;
Array<int> GPUS;
cudata_t** cudata_all=NULL;/*for all GPU. */
cuglobal_t* cuglobal=NULL; /*global data*/
#ifdef __APPLE__
pthread_key_t cudata_key;
#else
__thread cudata_t* cudata=NULL;/*for current thread and current GPU */
#endif
static __attribute((constructor)) void init(){
#ifdef __APPLE__
	pthread_key_create(&cudata_key, NULL);
#endif
}

/**
   Get GPU info.
*/
void gpu_dbg(){
	struct cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, 0);
	dbg("name=%s\n"
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
int current_gpu(){
	int igpu;
	if(cudata){
		igpu=cudata->igpu;
	}else{
		cudaGetDevice(&igpu);
	}
	return igpu;
}
/**
   Print memory consumption.
*/
void gpu_print_mem(const char* msg){
	size_t fr, tot;
	cudaDeviceSynchronize();
	DO(cudaMemGetInfo(&fr, &tot));
	info("GPU %d: mem used %ld MB (%s)\n", current_gpu(), (long)(tot-fr)/1024/1024, msg);
}
/**
   Get available memory.
*/
long gpu_get_free_mem(void){
	size_t fr, tot;
	DO(cudaMemGetInfo(&fr, &tot));
	return (long)fr;
}
/**
   Get available memory/maximum memory ratio. Returns 0-100.
*/
static long gpu_get_usage_percentage(int igpu, long minimum){
	size_t fr=0, tot=0;
	cudaDeviceProp prop;
	int ans;
	if((ans=cudaSetDevice(igpu))){
		info("GPU%2d cudaSetDevice failed with error %d: %s\n", igpu, ans, cudaGetErrorString((cudaError_t)ans));
		return 0;
	}else if((ans=cudaMemGetInfo(&fr, &tot))){
		info("GPU%2d cudaMemGetInfo failed with error %d: %s\n", igpu, ans, cudaGetErrorString((cudaError_t)ans));
		return 0;
	}else if((ans=cudaGetDeviceProperties(&prop, igpu))){
		info("GPU%2d cudaGetDeviceProperties failed with error %d: %s\n", igpu, ans, cudaGetErrorString((cudaError_t)ans));
		return 0;
	}else{
		info("GPU%2d is %s with arch %d.%d, %.1fGB free, %.1fGB total device memory.\n", igpu, prop.name, prop.major, prop.minor, fr*9.3e-10, tot*9.3e-10);
		if((long)fr>minimum){
			return (long)((tot-fr)*100./tot);
		} else{
			return 0;
		}
	}
}
/*static int cmp_long_ascend(const long *a, const long *b){
	if(*b<*a){
		return 1;
	} else if(*a<*b){
		return -1;
	} else{
		return 0;
	}
}*/
static int cmp_long2_ascend(const long* a, const long* b){
	if(b[1]<a[1]){
		return 1;
	} else if(a[1]<b[1]){
		return -1;
	} else{
		return 0;
	}
}
static int cmp_long3_ascend(const long *a, const long *b){
	if(b[2]<a[2]){
		return 1;
	} else if(a[2]<b[2]){
		return -1;
	} else{
		return 0;
	}
}
struct task_t{
	float timing;/*based on empirical data*/
	int* dest;
	char name[64];
};
static int task_cmp(const task_t* a, const task_t* b){
	if(b->timing>a->timing){
		return 1;
	} else if(b->timing<a->timing){
		return -1;
	} else{
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
int gpu_init(const parms_t* parms, int* gpus, int ngpu){
	//Set the CUDA CACHE size big enough if user didn't set it.
	setenv("CUDA_CACHE_MAXSIZE", "1000000000", 0);

	int ans;//total number of GPUs.
	if((ans=cudaGetDeviceCount(&MAXGPU))||MAXGPU==0){//no GPUs available.
		info2("No GPUs available. ans=%d (%s)\n", ans, cudaGetErrorString((cudaError_t)ans));
		return 0;
	}
	if(gpus && !ngpu){
		warning("gpus is set, but ngpu=0. set gpus to NULL.\n");
		gpus=NULL;
	}
		
	if(gpus){//user prescribed list of GPUs
		int ngpu2=0;
		for(int ig=0; ig<ngpu; ig++){
			//info("gpus[%d]=%d, ngpu=%d, ngpu2=%d\n", ig, gpus[ig], ngpu, ngpu2);
			if(gpus[ig]<0){
				ngpu2=0;//-1 disable all previous entry
			}else{
				int found=0;//check for duplicate entry. Skip if found.
				for(int jg=0; jg<ngpu2; jg++){
					if(gpus[jg]==gpus[ig]){
						found=1;
						break;
					}
				}
				if(!found){
					gpus[ngpu2]=gpus[ig];
					ngpu2++;
				}
			}
		}
		if(ngpu2==0){
			ngpu=-1;
		}else{
			ngpu=ngpu2;
		}
	}
	if(ngpu<0){
		info("CUDA is disabled by user.\n");
		NGPU=0;
		return 0;
	}
	long mem_minimum=0;
	if(parms){
		if(parms->gpu.evl||parms->gpu.wfs){
			const int nps=parms->atm.nps;
			long nxn=parms->atm.nxnmax;
			mem_minimum+=sizeof(Real)*nps*nxn*nxn;
		}
		if(parms->gpu.evl){
			mem_minimum+=sizeof(Real)*parms->evl.nevl*(long)pow(parms->aper.d/parms->evl.dx, 2);
		}
		if(parms->gpu.wfs&&!parms->sim.idealtomo){
			for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
				mem_minimum+=sizeof(Real)*parms->powfs[ipowfs].nwfs*(long)pow(parms->aper.d/parms->powfs[ipowfs].dx, 2)*2;
			}
		}
		if(parms->gpu.tomo||parms->gpu.fit){
			mem_minimum+=sizeof(Real)*parms->atmr.nps*(long)pow(parms->aper.d*parms->tomo.pos/parms->atmr.dx, 2)*2;
		}
		//mem_minimum*=2;
		if(mem_minimum==0){//gpu is disabled
			return 0;
		} else{
			info("CUDA: minimum memory requirement is %.1fGB\n", mem_minimum/(real)(1024*1024*1024));
		}
	} else{
		mem_minimum=1000000000;//1GB.
	}
	char fnlock[PATH_MAX];
	snprintf(fnlock, PATH_MAX, "%s/gpu.lock", TEMP);
	int fdlock=lock_file(fnlock, 1);
	/*
	  Create a mapping between CUDA device ordering and NVML ordering (nvidia-smi).
	  [0] is CUDA index
	  [1] is PCI bus id. = INT_MAX if disabled
	  [2] is memory usage ratio percentage (100 is full)
	*/
	
	long gpu_info[MAXGPU][4];//cuda index, and nvidia-smi index, and memory usage
	int gpu_valid_count=0;
	for(int ic=0; ic<MAXGPU; ic++){
		cudaDeviceProp properties;
		gpu_info[ic][0]=ic;//cuda-index
		gpu_info[ic][1]=0;//nvidia-smi reported index (order of pci bus id)
		gpu_info[ic][2]=100;//preset to 100%
		gpu_info[ic][3]=0;//flag 0: unknown or not usable or disabled, 1: usable. 2: being used.
		if(!cudaSetDevice(ic) &&!cudaGetDeviceProperties(&properties, ic)){
			gpu_info[ic][1]=properties.pciBusID;//nvidia-smi order (pci index)
			gpu_info[ic][3]=1;//mark usable
			gpu_valid_count++;
		}else{
			warning("Error getting information for CUDA GPU %d, skip.\n", ic);
		}
	}
	if(!gpu_valid_count){
		goto end;
	}
	//sort to obtain [1] as nvidia-smi index
	qsort(gpu_info, MAXGPU, sizeof(gpu_info[0]), (int(*)(const void *, const void *))cmp_long2_ascend);
	//check whether GPU is disabled by environment variable
	for(int ig=0; ig<MAXGPU; ig++){
		gpu_info[ig][1]=ig;//PCI-id ordering. (nvidia-smi order)
		if(gpu_info[ig][3]){
			char tmp[64];
			snprintf(tmp, sizeof(tmp), "MAOS_GPU_DISABLE_%d", ig);
			const char *ctmp=NULL;
			if((ctmp=getenv(tmp))&&ctmp[0]!='0'){
				info("GPU %d is disabled by environment variable %s=%s\n", ig, tmp, ctmp);
				gpu_info[ig][3]=0;
			}
		}
	}
	//gmap[igpu] is sorted by PCI bus id (nvidia-smi index)
	NGPU=0;
	if(ngpu==0||ngpu>MAXGPU){
		ngpu=MAXGPU;
	}
	GPUS.init(ngpu, 2); //stores CUDA index
	if(gpus){
		/*
		  	User specified exact GPUs using PCI bus index to use. We check every
	  		entry. If <0 is found, do not use any GPU. If >=MAXGPU is found,
	  		skip the GPU and print warning. If duplicates are found, use only
	  		once.              
	 	*/
		//sort by nvidia-smi (pcie) order
		
		for(int ii=0; ii<ngpu; ii++){
			int ig=gpus[ii];//user uses nvidia-smi ordering
			if(ig>=MAXGPU){
				warning("GPU %d does not exist\n", ig);
			}else if(gpu_info[ig][3]==1){
				info("GPU %d is not usable.\n", ig);
			} else{
				GPUS(NGPU, 0)=gpu_info[ig][0];
				GPUS(NGPU, 1)=gpu_info[ig][1];
				gpu_info[ig][3]=2;//mark used.
				NGPU++;
			}
		}
	} else{//Sort and use least busy GPUs.
		for(int ig=0; ig<MAXGPU; ig++){
			int ic=gpu_info[ig][0];//cuda index
			if(gpu_info[ig][3]){
				gpu_info[ig][2]=gpu_get_usage_percentage(ic, mem_minimum);
			}
		}
		/*sort so that gpus with more available memory is in the front.*/
		qsort(gpu_info, MAXGPU, sizeof(gpu_info[0]), (int(*)(const void*, const void*))cmp_long3_ascend);
		for(int ig=0; ig<ngpu; ig++){
			if(gpu_info[ig][3] && gpu_info[ig][2]<100){
				GPUS(NGPU,0)=(int)gpu_info[ig][0];//cuda index
				GPUS(NGPU,1)=(int)gpu_info[ig][1];//nvidia-smi index
				gpu_info[ig][3]=2;//mark used
				NGPU++;
			}
		}
	}
	for(int ig=0; ig<MAXGPU; ig++){
		if(gpu_info[ig][3]==1){
			dbg("GPU %ld is not used; reset its context.\n", gpu_info[ig][1]);
			cudaSetDevice(gpu_info[ig][0]);
			cudaDeviceReset();//releases context, frees memory.
		}
	}
	if(NGPU){
		cuglobal=new cuglobal_t;
		cudata_all=new cudata_t*[NGPU];

		info2("Using GPU");
		for(int i=0; i<NGPU; i++){
			//cudaInitDevice(GPUS(i, 0), 0, 0);//Does not seem to be necessary. cudaSetDevice does the same thing and more.
			cudaSetDevice(GPUS(i,0));
			cudata_all[i]=new cudata_t;//make sure allocation on the right gpu.
			gpu_set(i);//set cudata
			cudata->igpu=GPUS(i,1);//nvidia-smi index
			info2(" %d", GPUS(i,1));//cudata->igpu);
			//Reserve memory in GPU so the next maos will not pick this GPU.
			cudata->reserve.init(MEM_RESERVE, 1);
		}
		info2("\n");
		if(parms){
			/*Assign task to gpu evenly based on empirical data to maximum GPU
			 * usage. We first gather together all the tasks and assign a timing
			 * in ms to each. Sort all the tasks in descend order and then
			 * iteratively assign each task to the minimally used GPU*/
			cuglobal->evlgpu.init(parms->evl.nevl, 1);
			cuglobal->wfsgpu.init(parms->nwfs, 1);
			int ntask=0;
			if(parms->gpu.tomo||parms->gpu.fit||parms->gpu.lsr) ntask++;
			if(parms->gpu.evl) ntask+=parms->evl.nevl;
			if(parms->gpu.wfs) ntask+=parms->nwfs;
			if(ntask==0){
				info("GPU is not needed\n");
				delete[] cudata_all;
				NGPU=0;
				GPUS.deinit();
				goto end;
			}
			struct task_t* tasks=(task_t*)calloc(ntask, sizeof(task_t));
			//recon
			int count=0;
			if(parms->gpu.tomo||parms->gpu.fit||parms->gpu.lsr){
				tasks[count].timing=20;//ms
				tasks[count].dest=&cuglobal->recongpu;
				snprintf(tasks[count].name, 64, "RECON");
				count++;
			}
			//evl
			for(int ievl=0; parms->gpu.evl&&ievl<parms->evl.nevl; ievl++){
				if(parms->evl.psfmean==1&&parms->evl.psf->p[ievl]){
					tasks[count].timing=100;//ms; more time for PSF
				} else{
					tasks[count].timing=10;//ms
				}
				tasks[count].dest=cuglobal->evlgpu()+ievl;
				snprintf(tasks[count].name, 64, "EVL %d", ievl);
				count++;
			}
			//wfs
			for(int iwfs=0; parms->gpu.wfs&&iwfs<parms->nwfs; iwfs++){
				const int ipowfs=parms->wfs[iwfs].powfs;
				if(parms->powfs[ipowfs].type==WFS_PY){//pwfs is slower
					tasks[count].timing=50;
				} else if(parms->powfs[ipowfs].usephy){
					tasks[count].timing=20;
				} else{
					tasks[count].timing=10;
				}
				tasks[count].dest=cuglobal->wfsgpu()+iwfs;
				snprintf(tasks[count].name, 64, "WFS %d", iwfs);
				count++;
			}
			qsort(tasks, count, sizeof(task_t), (int(*)(const void*, const void*))task_cmp);
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
				dbg3("%s --> GPU %d\n", tasks[it].name, GPUS(*tasks[it].dest,1));
			}
			free(tasks);

			{
				info("Use %d GPUs for%s%s%s%s%s\n", NGPU,
					parms->gpu.wfs?" WFS":"",
					parms->gpu.lsr?" LSR":"",
					parms->gpu.tomo?" Tomo":"",
					parms->gpu.fit?" Fit":"",
					parms->gpu.evl?" Evl":""
				);
			}
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
	if(NGPU){
		for(int ig=0; ig<NGPU; ig++){
			gpu_set(ig);
			delete cudata;
		}
		delete[] cudata_all;
		delete cuglobal;
		GPUS.deinit();
		NGPU=0;
	}
}

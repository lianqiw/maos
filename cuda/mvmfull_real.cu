/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include <errno.h>
#include "utils.h"
#include "curmat.h"
#include "cucmat.h"
/**
\file mvmfull_real.cu

2013-11-12

Involved from mvmfull_iwfs.cu. Take into account real WFS geometry.
Accuracy Verified.
Test MVM for a single WFS, with networking. ethtest is the server.

two ways: 
1) broadcast gradients to both GPUs with each GPU handle part of actuators. (not tested yet)
2) partition gradients to GPUs with GPU handle all actuators. (selected)

Use schedtool -a 0x1 PID to let the exe only run one CPU 0. This prevents jitter.

For GTX590, use NOVER=2 and NGRID=30 gives best performance.
*/

static unsigned int event_flag=cudaEventDisableTiming;
typedef struct{
    curmat cumvm;//active mvm control matrix
    curmat cumvm_next;//inactive mvm control matrix.
    curmat cumvm1;
    curmat cumvm2;
    curcell mtch;
    Array<short,Gpu> pix;//pixels. Each sa has 15x6=90 pixels.
    Array<short,Gpu> pixbias;
    curmat grad;
    curmat gdm;/*add to grad*/
    curmat act;
    curcell actelse;
    Real FSMdelta; /*FSM actual angle difference from command.*/
    curcell im0;
    int mtch_isa;
    cuimat saind;
    stream_t stream_p;//pixels and other transportation across PCI-E.
    stream_t stream_g;//grads
    Array<stream_t> stream_a;//act
    stream_t stream_b;//background process
    int ism;//index of stream for mvm
    int count;
    int gpu;//Which GPU this data is for
    int istep;//Which time step we are in
    int copy_mvm;//1: need to copy mvm.
    int ic;//the column that we are copying.
    cudaEvent_t *event_p;
    cudaEvent_t *event_g;
    cudaEvent_t event_pall;
    event_t *event_w;
}GPU_DATA_T;
/*Does matched filter
  threadIdx.x is 32, for 1/3 of pixels
  threadIdx.y is for a few subapertures
  no need to sync threads as every 32 is in the same warp.
*/
static void __global__ mtch_do(const Real *mtch, const short *pix, const short *pixbias, Real *grad,
			       int *saind, int nsa){
    extern __shared__ Real cum[];//for cumulation and reduction
    for(int isa=threadIdx.y+blockDim.y*blockIdx.x; isa<nsa; isa+=blockDim.y*gridDim.x){
	const int npix=saind[isa+1]-saind[isa];
	Real *cumx=cum+threadIdx.y*blockDim.x*2;
	Real *cumy=cumx+blockDim.x;
	const short *pixi=pix+saind[isa];
	const short *pixbiasi=pixbias+saind[isa];
	const Real *mtchx=mtch+saind[isa]*2;
	const Real *mtchy=mtchx+npix;
	const int npix3=npix/3;
	//sum 3 times for max 90 pixels.
	int ipix=threadIdx.x;
	cumx[ipix]=0;
	cumy[ipix]=0;
	if(ipix<npix3){
	    cumx[ipix]=mtchx[ipix]*(pixi[ipix]-pixbiasi[ipix])
		+mtchx[ipix+npix3]*(pixi[ipix+npix3]-pixbiasi[ipix+npix3])
		+mtchx[ipix+npix3*2]*(pixi[ipix+npix3*2]-pixbiasi[ipix+npix3*2]);
	    cumy[ipix]=mtchy[ipix]*(pixi[ipix]-pixbiasi[ipix])
		+mtchy[ipix+npix3]*(pixi[ipix+npix3]-pixbiasi[ipix+npix3])
		+mtchy[ipix+npix3*2]*(pixi[ipix+npix3*2]-pixbiasi[ipix+npix3*2]);
	}
	//reduction
	for(int step=16;step>0;step>>=1){
	    if(ipix<step){
		cumx[ipix]+=cumx[ipix+step];
		cumy[ipix]+=cumy[ipix+step];
	    }
	}
	if(ipix==0){
	    grad[isa*2]=cumx[0];
	    grad[isa*2+1]=cumy[0];
	}
    }
}
/*
  Accumulate statistics
*/
static void __global__ dither_acc_do(const short *pix,
				     Real *im0, Real *imx, Real *imy, 
				     Real cd, Real sd, int totpix){
    for(int ipix=threadIdx.x+blockIdx.x*blockDim.x; ipix<totpix; ipix+=blockDim.x*gridDim.x){
	short ii=pix[ipix];//-pixbias[ipix];
	im0[ipix]+=ii;
	imx[ipix]+=ii*cd;
	imy[ipix]+=ii*sd;
    }
}

    // First calibrate out pixels
    // im0/=imc; imx*=(2/a2m*imc); imy*=2/(a2m*imc)
/**
   A standalone routine that testes applying MVM for a single WFS and update mvm.
   The orderig of gradients are like xyxyxy instead of normal xxxyyy.

   Important: 
   1) Only page locked host memory can do async memcpy that overallps with computation
   2) Has to be Portable for multiple GPUs to do async memcpy concurrently.
   
*/
void mvmfull_real(int *gpus, int ngpu, int nstep){
    dbg("Using %d gpus. nstep=%d\n", ngpu, nstep);
    int nstep0=nstep>1?1000:0;//for warm up
    //Load subaperture actual pixel numbers along radial direction and offset of pixel of each subaperture.
    dmat *d_saind=dread("NFIRAOS_saind");
    const int nsa=d_saind->nx-1;
    int *saind=(int*)malloc(sizeof(int)*(1+nsa));
    for(int i=0; i<nsa+1; i++){
	saind[i]=(int)d_saind->p[i];
    }
    dfree(d_saind);
    const int totpix=saind[nsa];
    const int nact=6981;//active subapertures.
    int ng=nsa*2;
    X(mat) *mvm1, *mvm2, *pix1, *pix2, *mtch, *ptt, *pixbias;
    X(mat) *im0;
    if(zfexist("mvm2.bin")){
	mvm1=X(read)("mvm1");
	mvm2=X(read)("mvm2");
	pix1=X(read)("pix1");
	pix2=X(read)("pix2");
	mtch=X(read)("mtch");
	ptt=X(read)("ptt");
	pixbias=X(read)("pixbias");
    }else{
	mvm1=X(new)(nact, ng);
	mvm2=X(new)(nact, ng);
	pix1=X(new)(totpix,1);
	pix2=X(new)(totpix,1);
	mtch=X(new)(totpix*2,1);
	ptt=X(new)(ng, 2);
	pixbias=X(new)(totpix, 1);
	rand_t srand;
	seed_rand(&srand, 1);
	X(randu)(mvm1,1e-7,&srand);
	X(randu)(mvm2,1e-7,&srand);
	X(randu)(mtch, 1, &srand);
	X(randu)(pix1,50, &srand);
	memcpy(pix2->p, pix1->p, sizeof(Real)*totpix);
	X(randu)(ptt, 1, &srand);
	//srandn(pixbias, 1, &srand);
    }
    X(mat) *mvm=mvm1;
    X(mat) *pix=pix2;
    //To receive statistics from GPU
    im0=X(new)(totpix,3);
    if(nstep==1){//Verify accuracy
	//We use half of the array as short.
	writearr("pix", 1, sizeof(short), M_INT16, NULL, pix->p, totpix, 1);
	writearr("pixbias", 1, sizeof(short), M_INT16, NULL, pixbias->p, totpix, 1);
	writebin(mvm1, "mvm1");
	writebin(mtch, "mtch");
    }
    X(cell) *dmres=X(cellnew)(ngpu, 1);
    X(pagelock)(im0, pix1, pix2, mvm1, mvm2, mtch, dmres, NULL);

    int port=20000;
    int sock=-1;
    int ready=1;
    int mtch_ngrid=50;//30;//can change to utilize GPU fully. 16 is good for cassiopeia
    const int mtch_dimx=32;//must launch 32 threads so that they belong to single wrap.
    const int mtch_dimy=12;//number of subapertures
    const int naeach=128;//Each block handle this many subapertures
    int nover=9;//determining how many blocks to launch.
    int nsm=2;//number of streams
    {
	char *MVM_NSM=getenv("MVM_NSM");
	if(MVM_NSM){
	    nsm=strtol(MVM_NSM, NULL, 10);
	    info("nsm is set to %d\n", nsm);
	}
	char *MVM_NGRID=getenv("MVM_NGRID");
	if(MVM_NGRID){
	    mtch_ngrid=strtol(MVM_NGRID, NULL, 10);
	    info("mtch_ngrid is set to %d\n", mtch_ngrid);
	}
	char *MVM_NOVER=getenv("MVM_NOVER");
	if(MVM_NOVER){
	    nover=strtol(MVM_NOVER, NULL, 10);
	}
    }
    const int sastep=mtch_dimy*mtch_ngrid;
    {
	char *MVM_CLIENT=getenv("MVM_CLIENT");
	if(MVM_CLIENT){
	    char *MVM_PORT=getenv("MVM_PORT");
	    if(MVM_PORT){
		port=strtol(MVM_PORT, NULL, 10);
	    }
	    info("Connecting to server %s\n", MVM_CLIENT);
	    sock=connect_port(MVM_CLIENT, port, 0 ,1);
	    if(sock!=-1) {
		info("Connected");
		int cmd[7];
		cmd[0]=nact;
		cmd[1]=nsa;
		cmd[2]=sastep;
		cmd[3]=totpix;
		cmd[4]=nstep;
		cmd[5]=nstep0;
		cmd[6]=2;
		if(stwriteintarr(sock, cmd, 7) 
		   || stwriteintarr(sock, saind, nsa+1)
		   || stwrite(sock, pix->p, sizeof(short)*totpix)){
		    close(sock); sock=-1;
		    warning("Failed: %s\n", strerror(errno));
		}
	    } else {
		info("Failed\n");
	    }
	}
    }

    const int nbuf=2;//two buffers for im0.
    //int dither_nsa=10;//each time step compute this many subapertures for matched filter
    //int comp_mtch_done[ngpu];
    //Real imc, a2m;//PLL results
    int nc=10;//each time copy nc column of mvm.
    GPU_DATA_T **data=new GPU_DATA_T*[ngpu];
    const int sect_gpu=(nsa+sastep*ngpu-1)/(sastep*ngpu);
    for(int igpu=0; igpu<ngpu; igpu++){
	cudaSetDevice(gpus[igpu]);
	data[igpu]=new GPU_DATA_T;
	data[igpu]->cumvm1=curmat(nact, ng);
	data[igpu]->cumvm2=curmat(nact, ng);
	data[igpu]->cumvm=data[igpu]->cumvm1;
	data[igpu]->cumvm_next=data[igpu]->cumvm2;
	cp2gpu(data[igpu]->cumvm1, mvm);
	data[igpu]->pix=Array<short,Gpu>(totpix, 1);
	data[igpu]->pixbias=Array<short,Gpu>(totpix, 1);
	cp2gpu(data[igpu]->pixbias(), (short*)pixbias->p, totpix*sizeof(short), cudaMemcpyHostToDevice);
	data[igpu]->mtch=curcell(nbuf, 1, totpix*2,1);
	cp2gpu(data[igpu]->mtch[0], mtch);
	data[igpu]->grad=curmat(ng, 1);
	data[igpu]->act=curmat(nact, 1);
	data[igpu]->actelse=curcell(ngpu-1, 1, nact, 1);
	data[igpu]->im0=curcell(nbuf,1, totpix, 3);//two buffers
	data[igpu]->FSMdelta=0.2;
	data[igpu]->stream_a=Array<stream_t>(nsm, 1);
	data[igpu]->event_w=new event_t[nsm];
	data[igpu]->gpu=gpus[igpu];
	data[igpu]->event_g=new cudaEvent_t[sect_gpu];
	data[igpu]->event_p=new cudaEvent_t[sect_gpu];
	for(int i=0; i<sect_gpu; i++){
	    cudaEventCreateWithFlags(&data[igpu]->event_g[i],event_flag);
	    cudaEventCreateWithFlags(&data[igpu]->event_p[i],event_flag);
	}
	cudaEventCreateWithFlags(&data[igpu]->event_pall,event_flag);
	dmres->p[igpu]=X(new)(nact, 1);
	X(pagelock)(dmres->p[igpu], NULL);
	data[igpu]->saind=cuimat(nsa+1,1);
	cp2gpu(data[igpu]->saind(), saind, nsa+1, 1, 0);
    }
    X(mat) *timing=X(new)(nstep, 1);
    X(mat) *timing_tot=X(new)(nstep, 1);
    X(mat) *timing_sock=X(new)(nstep, 1);
    cudaProfilerStop();
    cudaProfilerStart();
    TIC;
    if(sock!=-1 && stwriteint(sock, ready)){
	warning("error send ready signal: %s\n", strerror(errno));
	close(sock); sock=-1;
    }
    int nblock;
    info("Ready\n");
    int ibuf=0;//buffer to use for statistics
    int ibuf_mtch=0;//buffer for matched filter
    int ibuf_stat=0;//buffer for matched filter computation
    int mtch_down=0;
    int nset=(nsa+sastep-1)/sastep;
    char *copied_mtch=(char*)calloc(nset*3, sizeof(char));
    Real tim_tot=0, tim_min=INFINITY, tim_max=0;
    for(int jstep=-nstep0; jstep<nstep; jstep++){
	//run 20 frames to warm up before timing.
	int istep=jstep<0?0:jstep;
	if(sock!=-1){//start signal
	    timing_sock->p[istep]=0;
	}
	tic;
	if(nover>0){
	    nblock=(nact*nover+naeach-1)/naeach;
	}else{
	    nblock=(nact*(1+istep/50)+naeach-1)/naeach;
	}
	if((1+istep)%8000==0){//need to update MVM
	    if(mvm==mvm1){//switch mvm on host.
		mvm=mvm2;
	    }else{
		mvm=mvm1;
	    }
	    for(int igpu=0; igpu<ngpu; igpu++){
		data[igpu]->copy_mvm=1;
		if(data[igpu]->ic!=0){
		    warning("Sync error, skip update request at step %d\n", istep);
		}
	    }
	}
	for(int igpu=0; igpu<ngpu; igpu++){
	    data[igpu]->ism=-1;
	    data[igpu]->count=0;
	    data[igpu]->istep=istep;
	}
	if(sock==-1){
	    if(pix==pix1){
		pix=pix2;
	    }else{
		pix=pix1;
	    }
	}

	for(int isa=0, igpu=0; isa<nsa; isa+=sastep, igpu=((igpu+1)%ngpu)){
	    cudaSetDevice(gpus[igpu]); 
	    GPU_DATA_T *datai=data[igpu];
	    int npixleft;
	    int nsaleft;
	    if(nsa<isa+sastep){//terminate
		npixleft=totpix-saind[isa];
		nsaleft=nsa-isa;
	    }else{
		npixleft=saind[isa+sastep]-saind[isa];
		nsaleft=sastep;
	    }
	    //One stream handling the memcpy
	    short *pcur=(short*)(pix->p)+saind[isa];
	    if(sock!=-1){
		real tmp0=myclockd();
		if(stread(sock, pcur, sizeof(short)*npixleft)){
		    warning("failed: %s\n", strerror(errno));
		    close(sock); sock=-1;
		    _Exit(1);
		}
		timing_sock->p[istep]+=myclockd()-tmp0;
	    }
	 
	    DO(cudaMemcpyAsync(datai->pix()+saind[isa], pcur, sizeof(short)*npixleft,
			       cudaMemcpyHostToDevice, datai->stream_p));
	    //Recored the event when the memcpy is finished
	    DO(cudaEventRecord(datai->event_p[datai->count], datai->stream_p));
	    //Start matched filter when pixel transfer is done.
	    DO(cudaStreamWaitEvent(datai->stream_g, datai->event_p[datai->count], 0));

	    mtch_do<<<mtch_ngrid, dim3(mtch_dimx, mtch_dimy), 
		mtch_dimx*mtch_dimy*sizeof(Real)*2, datai->stream_g>>>
		(datai->mtch[ibuf_mtch], datai->pix, datai->pixbias, 
		 datai->grad()+isa*2, datai->saind()+isa, nsaleft);
	    //Record the event when matched filter is done
	    DO(cudaEventRecord(datai->event_g[datai->count], datai->stream_g));
	    //Another stream does the matrix vector multiplication. Wait for the event before executing.
	    //The stream stream will wait only for the completion of the most recent host call to cudaEventRecord() on event
	    datai->ism=(datai->ism+1)%nsm;
	    /*
	      Todo: 
	      *) Project PTT modes. Generate LLT FSM signal.
	      *) Minus LPF focus mode
	      *) Add gradient from DM to form PSOL grads.
	    */
	    cudaStreamWaitEvent(datai->stream_a[datai->ism], datai->event_g[datai->count], 0);

#if 0
	    Real one=1;
	    DO(CUBL(gemv)(datai->stream_a[datai->ism], CUBLAS_OP_N, nact, nsaleft*2, &one, datai->cumvm->p+nact*isa*2, nact, datai->grad->p+isa*2, 1, &one, datai->act->p, 1));
#else
	    multimv_do<<<nblock, naeach, sizeof(Real)*naeach, datai->stream_a[datai->ism]>>>
		(datai->cumvm()+nact*isa*2, datai->act, datai->grad()+isa*2, 
		 nact, nsaleft*2);
#endif
	    DO(cudaEventRecord(datai->event_w[datai->ism], datai->stream_a[datai->ism]));
	    datai->count++;
	}
	for(int igpu=0; igpu<ngpu; igpu++){
	    GPU_DATA_T *datai=data[igpu];
	    //Record an event when pixel tranporting is over. So we can start transporting mvm matrix.
	    DO(cudaEventRecord(datai->event_pall, datai->stream_p));
	}
	/*Accumulate statistics in each cycle. The pixels are present in different GPUs*/
	for(int isa=0, igpu=0; isa<nsa; isa+=sastep, igpu=((igpu+1)%ngpu)){
	    cudaSetDevice(gpus[igpu]); 
	    GPU_DATA_T *datai=data[igpu];
	    int npixleft;
	    if(nsa<isa+sastep){//terminate
		npixleft=totpix-saind[isa];
	    }else{
		npixleft=saind[isa+sastep]-saind[isa];
	    }
	    real theta=M_PI*0.5*istep+datai->FSMdelta;
	    Real cd=cos(theta);
	    Real sd=cos(theta);
	    //Do not start before pixels are transported
#if 0
	    for(int ism=1; ism<nsm; ism++){//wait for MVM
		DO(cudaStreamWaitEvent(datai->stream_b, datai->event_w[ism], 0));
	    }
#else
	    //Wait for pixel transfer
	    cudaStreamWaitEvent(datai->stream_b, datai->event_pall, 0);
#endif
	    dither_acc_do<<<DIM(npixleft, 256), 0, datai->stream_b>>>
		(datai->pix()+saind[isa], datai->im0[ibuf]()+saind[isa], 
		 datai->im0[ibuf]()+totpix+saind[isa],datai->im0[ibuf]()+totpix*2+saind[isa],
		 cd, sd, npixleft);
	}
	//Download statistics to CPU for matched filter building.
	if(mtch_down || nstep==1){
	    cudaStream_t stream;
	    int iset=0;
	    for(int icol=0; icol<3; icol++){
		for(int isa=0, igpu=0; isa<nsa; isa+=sastep, igpu=((igpu+1)%ngpu), iset++){
		    if(copied_mtch[iset]) continue;
		    cudaSetDevice(gpus[igpu]); 
		    if(nstep==1){
			stream=data[igpu]->stream_b;
		    }else{
			stream=data[igpu]->stream_p;
		    }
		    int npixleft;
		    if(nsa<isa+sastep){//terminate
			npixleft=totpix-saind[isa];
		    }else{
			npixleft=saind[isa+sastep]-saind[isa];
		    }
		    DO(cudaMemcpyAsync(im0->p+saind[isa]+icol*totpix, 
				       data[igpu]->im0[ibuf_stat]()+saind[isa]+icol*totpix,
				       sizeof(Real)*npixleft,
				       cudaMemcpyDeviceToHost, stream));
		    copied_mtch[iset]=1;
		    if(nstep!=1) goto endhere;
		}
	    }
	    mtch_down=0;//completed
	  endhere:;
	}
	//Queue copying MVM matrix to second slot.
	for(int igpu=0; igpu<ngpu; igpu++){
	    GPU_DATA_T *datai=data[igpu];
	    if(datai->copy_mvm){
		int done=0, nsaleft;
		if(mvm->ny-datai->ic < nc){
		    done=1;
		    nsaleft=mvm->ny-datai->ic;
		}else{
		    nsaleft=nc;
		}
		if(datai->ic==0){
		    info("step %d: gpu %d uploading mvm\n", istep, igpu);
		}
		DO(cudaMemcpyAsync(datai->cumvm_next()+datai->ic*mvm->nx, 
				   mvm->p+datai->ic*mvm->nx, sizeof(Real)*mvm->nx*nsaleft, 
				   cudaMemcpyHostToDevice, datai->stream_p));
		
		datai->ic+=nsaleft;
		if(done){
		    datai->ic=0;
		    datai->copy_mvm=0;
		    curmat tmp=datai->cumvm;
		    datai->cumvm=datai->cumvm_next;
		    datai->cumvm_next=tmp;
		    info("gpu %d switched over at step %d\n", datai->gpu, datai->istep);
		}
	    }
	}
#if 1	//CPU sum
	//Copy DM commands back to CPU
	for(int igpu=0; igpu<ngpu; igpu++){
	    GPU_DATA_T *datai=data[igpu];
	    cudaSetDevice(gpus[igpu]); 
	    for(int ism=1; ism<nsm; ism++){
		DO(cudaStreamWaitEvent(datai->stream_a[0], datai->event_w[ism], 0));
	    }
	    cudaMemcpyAsync(dmres->p[igpu]->p, datai->act, nact*sizeof(Real), 
			    cudaMemcpyDeviceToHost, datai->stream_a[0]);
	    cuzero(datai->act, datai->stream_a[0]);
	}
	//CPU sums them together. sync first gpu
	data[0]->stream_a[0].sync();
	//sum other GPUs
	for(int igpu=1; igpu<ngpu; igpu++){
	    cudaSetDevice(gpus[igpu]); 
	    data[igpu]->stream_a[0].sync();
	    for(int iact=0; iact<nact; iact++){
		dmres->p[0]->p[iact]+=dmres->p[igpu]->p[iact];
	    }
	}
#else //GPU sum
	//First copy second GPU to first GPU.
	for(int igpu=1; igpu<ngpu; igpu++){
	    GPU_DATA_T *datai=data[igpu];
	    cudaSetDevice(gpus[igpu]); 
	    for(int ism=1; ism<nsm; ism++){
		DO(cudaStreamWaitEvent(datai->stream_a[0], datai->event_w[ism], 0));
	    }
	    cudaMemcpyAsync(data[0].actelse->p[igpu-1]->p, datai->act->p, nact*sizeof(Real), 
			    cudaMemcpyDeviceToDevice, datai->stream_a[0]);
	}
	if(ngpu>1){
	    if(ngpu==2){
		int igpu=0;
		cudaSetDevice(gpus[igpu]);
		GPU_DATA_T *datai=data[igpu];
		for(int ism=1; ism<nsm; ism++){
		    DO(cudaStreamWaitEvent(datai->stream_a[0], datai->event_w[ism], 0));
		}
		add_do<<<DIM(nact, 256), 0, datai->stream_a[0]>>>
		    (datai->act->p, datai->actelse->p[0]->p, (Real*)0, 1, nact);
		cudaMemcpyAsync(dmres->p[0]->p, datai->act->p, nact*sizeof(Real), 
				cudaMemcpyDeviceToHost, datai->stream_a[0]);
		datai->stream_a[0].sync();
	    }else{
		error("Please implement\n");
	    }
	}
#endif
	/*
	  Save resutls for debugging.
	 */
	//usleep(50);//yield
	if(nstep==1){//save result for verifying accuracy
	    writebin(dmres->p[0], "dmres");
	    for(int igpu=0; igpu<ngpu; igpu++){
		cudaSetDevice(gpus[igpu]); 
		cudaMemcpy(pix->p, data[igpu]->pix, sizeof(short)*totpix, cudaMemcpyDefault);
		char fn[PATH_MAX];
		snprintf(fn, PATH_MAX, "pix_gpu%d", igpu);
		writearr(fn, 1, sizeof(short), M_INT16, NULL, pix->p, totpix, 1);
		cuwrite(data[igpu]->grad, "grad_gpu%d", igpu);
	    }
	}
	/*
	  ToDO: the following background process.
	  *) Update PLL each time step
	  *) Output PLL results every 240 steps
	 */

	if(istep>0 && (istep+1)%2400==0){
	    ibuf_stat=ibuf;
	    info("Download statistics to CPU at step %d\n", istep);
	    ibuf=(ibuf+1)%nbuf;
	    mtch_down=1;
	    memset(copied_mtch, 0, sizeof(char)*nset*3);
	}
	for(int igpu=0; igpu<ngpu; igpu++){
	    GPU_DATA_T *datai=data[igpu];
	    cudaSetDevice(datai->gpu);
	    //no need to zero gradients.
	    datai->stream_b.sync();
	    datai->stream_p.sync();
	}
	if(sock!=-1){
	    real tmp0=myclockd();
	    if(stwrite(sock, dmres->p[0]->p, sizeof(Real)*nact)){
		warning("error write dmres: %s\n", strerror(errno));
		close(sock); sock=-1;
		_Exit(1);
	    }
	    if(streadint(sock, &ready)){//acknowledgement.
		warning("error read ack failed: %s\n", strerror(errno));
		close(sock), sock=-1;
		_Exit(1);
	    }
	    timing_sock->p[istep]+=myclockd()-tmp0;
	    timing->p[istep]=ready*1.e-6;
	}else{
	    timing->p[istep]=toc3;//do not tic.
	}
	timing_tot->p[istep]=toc3;
	if(istep>0){
	    tim_tot+=timing->p[istep];
	    if(tim_min>timing->p[istep]) tim_min=timing->p[istep];
	    if(tim_max<timing->p[istep]) tim_max=timing->p[istep];
	    if(istep%1000==0){
		info("Step %d, mean=%g, min=%g, max=%g ms\n", istep, tim_tot/(istep)*1e3, tim_min*1e3, tim_max*1e3);
	    }
	}
    }
    cudaProfilerStop();
    //writebin(dmres->p[0], "dmres");
    
    writebin(timing, "timing_%s_%dgpu", HOST, ngpu);
    writebin(timing_tot, "timing_tot_%s_%dgpu", HOST, ngpu);
    writebin(timing_sock, "timing_sock_%s_%dgpu", HOST, ngpu);
    X(pageunlock)(pix1, pix2, mvm1, mvm2, NULL);
    X(free)(mvm1);
    X(free)(mvm2);
    X(free)(pix1);
    X(free)(pix2);
    X(free)(mtch);
    X(cellfree)(dmres);
    X(free)(timing);
    X(free)(timing_tot);
    X(free)(timing_sock);
    for(int igpu=0; igpu<ngpu; igpu++){
	cudaSetDevice(gpus[igpu]);
	delete[] data[igpu]->event_w;
	delete[] data[igpu]->event_g;
	delete[] data[igpu]->event_p;
	cudaDeviceReset();
    }
    free(data);
  
}

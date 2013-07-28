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
#include "accphi.h"
#include "cucmat.h"
#include "curmat.h"
#include "cudata.h"
#include "tomo.h"

/*Prepare data for GP operation in GPU*/
void prep_GP(cumat<int>**GPp, float *GPscale, cusp **GPf,
	     const dsp *GP, const loc_t *saloc, const loc_t *ploc){
    if(!GP){ 
	error("GP is required\n");
    }
    int pos=(int)round(saloc->dx/ploc->dx);
    if(pos==1 || pos==2){//normally true
	dsp *GPt=sptrans(GP);
	const spint *pp=GPt->p;
	const spint *pi=GPt->i;
	const double *px=GPt->x;
	//convert the max float to max 2 byte integer
	const double pxscale=floor(32767./maxabs(px, GPt->nzmax));
	const int np1=pos+1;
	const int np=np1*np1;
	const int zmax=pos;
	int nsa=saloc->nloc;
	short2 *partxy=(short2*)calloc(sizeof(short2),np*nsa);//need to zero memory
	double dx1=1./ploc->dx;
	double dy1=1./ploc->dy;
	for(int ic=0; ic<GPt->n; ic++){
	    int isa=(ic<nsa)?ic:(ic-nsa);
	    for(spint ir=pp[ic]; ir<pp[ic+1]; ir++){
		int ix=pi[ir];
		double lx=ploc->locx[ix];
		double ly=ploc->locy[ix];
		double sx=saloc->locx[isa];
		double sy=saloc->locy[isa];
		int zx=(int)round((lx-sx)*dx1);
		int zy=(int)round((ly-sy)*dy1);
		/**
		   When the points used to generate GP align well
		   with the subaperture edge, the coupled points are
		   confined within the subaperture.
		*/
		if(zx<0 || zx>zmax || zy<0 || zy>zmax){
		    warning("isa=%d, zxy=%d %d\n", isa, zx, zy);
		}
		if(zx<0) zx=0;
		if(zx>zmax) zx=zmax;
		if(zy<0) zy=0;
		if(zy>zmax) zy=zmax;
		if(ic<nsa){
		    partxy[np*isa+zx+zy*np1].x+=(short)round(px[ir]*pxscale);
		}else{
		    partxy[np*isa+zx+zy*np1].y+=(short)round(px[ir]*pxscale);
		}
	    }
	}
	spfree(GPt);
	*GPp=new cumat<int>(np, nsa);
	cudaMemcpy((*GPp)->p, partxy, sizeof(int)*np*nsa, cudaMemcpyHostToDevice);
	*GPscale=1./pxscale;
	free(partxy);
    }else{/*use sparse */
	cp2gpu(GPf, GP, 1);
    }
}
static cumat<int> *
prep_saptr(loc_t *saloc, map_t *pmap){
    /*saloc mapped onto pmap*/
    int nsa=saloc->nloc;
    int (*saptr)[2]=new int[nsa][2];
    const float dx1=1./pmap->dx;
    const float dy1=1./pmap->dy;
    const float ox=pmap->ox;
    const float oy=pmap->oy;
    const double *salocx=saloc->locx;
    const double *salocy=saloc->locy;
    for(int isa=0; isa<nsa; isa++){
	saptr[isa][0]=(int)roundf((salocx[isa]-ox)*dx1);
	saptr[isa][1]=(int)roundf((salocy[isa]-oy)*dy1);
    }
    cumat<int> *saptr_gpu=new cumat<int>(2, nsa);
    DO(cudaMemcpy(saptr_gpu->p, saptr, nsa*2*sizeof(int), cudaMemcpyHostToDevice));
    delete [] saptr;
    return saptr_gpu;
}
static curmat *get_neai(dsp *nea){
    spint *pp=nea->p;
    spint *pi=nea->i;
    double *px=nea->x;
    int nsa=nea->ny/2;
    float (*neai)[3]=(float(*)[3])calloc(3*nsa, sizeof(float));
    for(int ic=0; ic<nea->n; ic++){
	for(spint ir=pp[ic]; ir<pp[ic+1]; ir++){
	    int ix=pi[ir];
	    int isa=ic<nsa?ic:ic-nsa;
	    float val=(float)px[ir]*TOMOSCALE;
	    if(ix==ic){/*diagonal part. */
		if(ic==isa){/*x */
		    neai[isa][0]=val;
		}else{/*y */
		    neai[isa][1]=val;
		}
	    }else if(ix==ic-nsa || ix==ic+nsa){/*cross part. symmetric. */
		neai[isa][2]=val;
	    }else{
		error("saneai has invalid format\n");
	    }
	}
    }
    curmat *neai_gpu=curnew(3, nsa);
    DO(cudaMemcpy(neai_gpu->p, neai, 3*nsa*sizeof(float), cudaMemcpyHostToDevice));
    free(neai);
    return neai_gpu;
}
#define TIMING 0
void cutomo_grid::init_hxdata(const PARMS_T *parms, const RECON_T *recon){
    PROP_WRAP_T *HXDATA=new PROP_WRAP_T[nwfs*recon->npsr];
    if(!hxdata){
	DO(cudaMalloc(&hxdata, sizeof(PROP_WRAP_T)*nwfs*recon->npsr));
    }
    for(int ips=0; ips<recon->npsr; ips++){ 
	const float ht=recon->ht->p[ips]; 
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    const int ipowfs = parms->wfsr[iwfs].powfs;
	    if(!parms->powfs[ipowfs].skip){
		const float hs = parms->powfs[ipowfs].hs; 
		const float scale = 1.f - ht/hs; 
		float dispx=parms->wfsr[iwfs].thetax*ht; 
		float dispy=parms->wfsr[iwfs].thetay*ht; 
		if(parms->tomo.predict && grid->vx){
		    dispx+=grid->vx[ips]*grid->dt*grid->delay;
		    dispy+=grid->vy[ips]*grid->dt*grid->delay;
		}
		cugrid_t pmapscale=grid->pmap.scale(scale);
		gpu_prop_grid_prep(HXDATA+iwfs+ips*nwfs, pmapscale, grid->xmap[ips],
				   dispx, dispy, NULL); 
		{
		    float tmp=laplacian_coef(recon->r0, recon->wt->p[ips], recon->xmap[ips]->dx)*0.25f;
		    HXDATA[iwfs+ips*nwfs].l2c=tmp*tmp*TOMOSCALE;
		    if(parms->tomo.piston_cr){
			HXDATA[iwfs+ips*nwfs].zzi=loccenter(recon->xloc[ips]);
			HXDATA[iwfs+ips*nwfs].zzv=tmp*tmp*TOMOSCALE*1e-6;
		    }else{
			HXDATA[iwfs+ips*nwfs].zzi=-1;
		    }
		}
	    }
	    HXDATA[iwfs+ips*nwfs].togpu(&hxdata[iwfs+ips*nwfs]);
	}
    }
    delete [] HXDATA;
}
void cutomo_grid::init_fdpcg(FDPCG_T *fdpcg, curecon_geom *grid){
    precond=new cufdpcg_t(fdpcg, grid);
}
void cutomo_grid::init(const PARMS_T *parms, const RECON_T *recon, const POWFS_T *powfs){
    if(!parms) return;
    nwfs=parms->nwfsr;

    if(recon->PTT && !PTT){//for t/t proj in 1)uplink t/t 2) recon
	cp2gpu(&PTT, recon->PTT);
    }
    ptt=!parms->recon.split || parms->tomo.splitlrt; 
    {
	PDF=curcellnew(parms->npowfs, 1);
	PDFTT=curcellnew(parms->npowfs, 1);
	dcell *pdftt=NULL;
	dcellmm(&pdftt, recon->PDF, recon->TT, "nn", 1);
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    if(parms->powfs[ipowfs].dfrs){//deprecated.	
		/*We only use the first diagonal block for each powfs. The
		  off diagonal is simply -1/(nlgswfs-1) times the diagonal block*/
		int iwfs1=parms->powfs[ipowfs].wfs[1];//second wfs
		cp2gpu(&PDF->p[ipowfs], recon->PDF->p[iwfs1*nwfs+iwfs1]);
		if(parms->powfs[ipowfs].trs){
		    /*coupling between TT and DF modes. 
		      We desire (I-DF*PDF)(I-TT*PTT)g=(I-TT*PTT-DF*PDF+DF*PDF*TT*PTT)g
		      So we first compute tt=PTT*g; df=PDF*g; then
		      g2=(I-TT*tt-DF*(df-(PDF*TT)*tt))
		      Here we record the values of PDF*TT
		    */
		    cp2gpu(&PDFTT->p[ipowfs], pdftt->p[iwfs1*nwfs+iwfs1]);
		}
	    }
	}
	dcellfree(pdftt);
    }
    {
	const int npowfs=parms->npowfs;
	GPp=new cucell<int>(npowfs, 1);
	GP=new cuspcell(npowfs, 1);
	GPscale=new float[npowfs];
	saptr=new cucell<int>(npowfs, 1);
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    if(parms->powfs[ipowfs].skip) continue;
	    prep_GP(&GPp->p[ipowfs], &GPscale[ipowfs], &GP->p[ipowfs],
		    recon->GP->p[ipowfs], powfs[ipowfs].saloc, recon->ploc);
	    saptr->p[ipowfs]=prep_saptr(powfs[ipowfs].saloc, recon->pmap);
	}
    }

    {
	neai=curcellnew(parms->nwfsr, 1);
	/*convert recon->saneai to gpu. */
	for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	    int ipowfs=parms->wfsr[iwfs].powfs;
	    int iwfs0=parms->recon.glao?iwfs:parms->powfs[ipowfs].wfs[0];/*first wfs in this group. */
	    if(iwfs!=iwfs0 && recon->saneai->p[iwfs+iwfs*parms->nwfsr]->p
	       ==recon->saneai->p[iwfs0+iwfs0*parms->nwfsr]->p){
		neai->p[iwfs]=neai->p[iwfs0]->ref();
	    }else{
		dsp *nea=recon->saneai->p[iwfs+iwfs*parms->nwfsr];
		neai->p[iwfs]=get_neai(nea);
	    }
	}/*for iwfs */
    }

    {
	init_hxdata(parms, recon);
    }
    {
	GPU_GP_T *GPDATA=new GPU_GP_T[nwfs];
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    const int ipowfs = parms->wfsr[iwfs].powfs;
	    if(parms->powfs[ipowfs].skip) continue;
	    GPDATA[iwfs].ipowfs=ipowfs;
	    GPDATA[iwfs].nwfs=parms->powfs[ipowfs].nwfsr;
	    GPDATA[iwfs].jwfs=parms->powfs[ipowfs].wfsind[iwfs];//wfs index in this group
	    GPDATA[iwfs].dsa=powfs[ipowfs].pts->dsa;
	    GPDATA[iwfs].pos=parms->tomo.pos;
	    GPDATA[iwfs].saptr=(int(*)[2])saptr->p[ipowfs]->p;
	    GPDATA[iwfs].GPp=(short2*)GPp->p[ipowfs]->p;
	    GPDATA[iwfs].GPscale=GPscale[ipowfs];
	    GPDATA[iwfs].PTT=PTT->p[iwfs+iwfs*nwfs]->p;
	    if(parms->powfs[ipowfs].dfrs){
		GPDATA[iwfs].PDF=PDF->p[ipowfs]->p;
		GPDATA[iwfs].PDFTT=PDFTT->p[ipowfs]->p;
	    }

	    GPDATA[iwfs].neai=(const float(*)[3])neai->p[iwfs]->p;
	    GPDATA[iwfs].nsa=powfs[ipowfs].pts->nsa;
	    GPDATA[iwfs].nxp=recon->pmap->nx;
	    GPDATA[iwfs].dxp=recon->pmap->dx;
	    GPDATA[iwfs].dyp=recon->pmap->dy;
	    GPDATA[iwfs].oxp=recon->pmap->ox;
	    GPDATA[iwfs].oyp=recon->pmap->oy;
	}
	DO(cudaMalloc(&gpdata, sizeof(GPU_GP_T)*nwfs));
	DO(cudaMemcpy(gpdata, GPDATA, sizeof(GPU_GP_T)*nwfs, cudaMemcpyHostToDevice));
	delete [] GPDATA;
    }

    if(parms->tomo.precond==1){
	init_fdpcg(recon->fdpcg, grid);
    }

    {/**Initialize temporary data*/
	int nxp=recon->pmap->nx;
	int nyp=recon->pmap->ny;
	int nxpw[nwfs], nypw[nwfs], ngw[nwfs];
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    const int ipowfs = parms->wfsr[iwfs].powfs;
	    if(parms->powfs[ipowfs].skip){
		nxpw[iwfs]=0;
		nypw[iwfs]=0;
		ngw[iwfs]=0;
	    }else{
		nxpw[iwfs]=nxp;
		nypw[iwfs]=nyp;
		ngw[iwfs]=powfs[ipowfs].pts->nsa*2;
	    }
	}

	opdwfs=curcellnew(nwfs, 1, nxpw, nypw);
	grad=curcellnew(nwfs, 1, ngw, (int*)NULL);
	ttf=curnew(3*nwfs, 1);
    }
}

/*
  If merge the operation in to gpu_prop_grid_do, need to do atomic
operation because previous retracing starts at offo, not 0.  */
__global__ void gpu_laplacian_do(PROP_WRAP_T *datai, float **outall, float **inall, int nwfs, float alpha){
    float *restrict in, *restrict out;
    
    const int ips=blockIdx.z;
    in=inall[ips];
    out=outall[ips];
    datai+=nwfs*ips;
    
    const int nx=datai->nxps;
    const int ny=datai->nyps;
    const float alpha2=datai->l2c*alpha;
    const int stepx=blockDim.x*gridDim.x;
    const int stepy=blockDim.y*gridDim.y;
    const int nx1=nx-1;
    const int ny1=ny-1;
    const int ix0=blockIdx.x*blockDim.x+threadIdx.x;
    const int iy0=blockIdx.y*blockDim.y+threadIdx.y;
    for(int iy=iy0; iy<ny; iy+=stepy){
	int iy1 =iy+1; if(iy1>ny1) iy1-=ny;
	int iy2 =iy+2; if(iy2>ny1) iy2-=ny;
	int iy1p=iy-1; if(iy1p<0) iy1p+=ny;
	int iy2p=iy-2; if(iy2p<0) iy2p+=ny;
	for(int ix=ix0; ix<nx; ix+=stepx){
	    int ix1 =ix+1; if(ix1>nx1) ix1-=nx;
	    int ix2 =ix+2; if(ix2>nx1) ix2-=nx;
	    int ix1p=ix-1; if(ix1p<0) ix1p+=nx;
	    int ix2p=ix-2; if(ix2p<0) ix2p+=nx;
	    out[ix+iy*nx]+=alpha2*(20.f*in[ix+iy*nx]
				   -8.f*(in[ix1p+iy*nx]+in[ix1+iy*nx]+in[ix+iy1p*nx]+in[ix+iy1*nx])
				   +2.f*(in[ix1p+iy1p*nx]+in[ix1+iy1p*nx]+in[ix1p+iy1*nx]+in[ix1+iy1*nx])
				   +(in[ix+iy2p*nx]+in[ix2p+iy*nx]+in[ix2+iy*nx]+in[ix+iy2*nx]));
	}
    }
    if(datai->zzi>-1){/*piston constaint*/
	if(threadIdx.x==0 && threadIdx.y==0 && blockIdx.x==0 && blockIdx.y==0){
	    out[datai->zzi]+=in[datai->zzi]*datai->zzv*alpha;
	}
    }
}

/*
  The third grid dimension tells the wfs to handle. 
*/
#define DIM_GP 128
__global__ static void gpu_gp_do(GPU_GP_T *data, float **gout, float *ttout, float *dfout, float **wfsopd, int ptt){
    __shared__ float gx[DIM_GP];
    __shared__ float gy[DIM_GP];
    __shared__ float gdf[DIM_GP];
    const int iwfs=blockIdx.z;
    const int nwfs=gridDim.z;
    GPU_GP_T *datai=data+iwfs;
    const int pos=datai->pos;
    if(!pos) return;
    const int nsa=datai->nsa;
    const int step=blockDim.x * gridDim.x;
    int (*restrict saptr)[2]=datai->saptr;
    float *restrict g=gout[iwfs];
    float GPscale=datai->GPscale;
    if(wfsopd){
	const float *restrict map=wfsopd[iwfs];
	const short2 *restrict pxy=datai->GPp;
	int nx=datai->nxp;
	/*GP operation.*/
	if(pos==1){
	    for(int isa=blockIdx.x * blockDim.x + threadIdx.x; isa<nsa; isa+=step){
		int ix=saptr[isa][0];
		int iy=saptr[isa][1];
	
		const short2 *restrict pxy2=pxy+isa*4;
		g[isa]=GPscale*(
		    +map[iy*nx+ix  ]*pxy2[0].x
		    +map[iy*nx+ix+1]*pxy2[1].x
		    +map[(iy+1)*nx+ix ] *pxy2[2].x
		    +map[(iy+1)*nx+ix+1]*pxy2[3].x);

		g[isa+nsa]=GPscale*(
		    +map[iy*nx+ix  ]*pxy2[0].y
		    +map[iy*nx+ix+1]*pxy2[1].y
		    +map[(iy+1)*nx+ix ] *pxy2[2].y
		    +map[(iy+1)*nx+ix+1]*pxy2[3].y);
	    }/*for isa */
	}else{
	    for(int isa=blockIdx.x * blockDim.x + threadIdx.x; isa<nsa; isa+=step){
		int ix=saptr[isa][0];
		int iy=saptr[isa][1];
		const short2 *restrict pxy2=pxy+isa*9;
		g[isa]=GPscale*(
		    +map[iy*nx+ix  ]*pxy2[0].x
		    +map[iy*nx+ix+1]*pxy2[1].x
		    +map[iy*nx+ix+2]*pxy2[2].x
		    +map[(iy+1)*nx+ix ] *pxy2[3].x
		    +map[(iy+1)*nx+ix+1]*pxy2[4].x
		    +map[(iy+1)*nx+ix+2]*pxy2[5].x
		    +map[(iy+2)*nx+ix ] *pxy2[6].x
		    +map[(iy+2)*nx+ix+1]*pxy2[7].x
		    +map[(iy+2)*nx+ix+2]*pxy2[8].x);
		g[isa+nsa]=GPscale*(
		    +map[iy*nx+ix  ]*pxy2[0].y
		    +map[iy*nx+ix+1]*pxy2[1].y
		    +map[iy*nx+ix+2]*pxy2[2].y
		    +map[(iy+1)*nx+ix ] *pxy2[3].y
		    +map[(iy+1)*nx+ix+1]*pxy2[4].y
		    +map[(iy+1)*nx+ix+2]*pxy2[5].y
		    +map[(iy+2)*nx+ix ] *pxy2[6].y
		    +map[(iy+2)*nx+ix+1]*pxy2[7].y
		    +map[(iy+2)*nx+ix+2]*pxy2[8].y);
	    }/*for isa */
	}
    }
    /* Global TT, Diff-Focus projection. Modifed from previous kernel so that
       each thread handle the same subaperture as previous gradient operation to
       avoid synchronization */
    if(datai->PTT && ptt){ //temp
	float (*restrict PTT)[2]=(float(*)[2])datai->PTT;
	gx[threadIdx.x]=0;
	gy[threadIdx.x]=0;
	for(int isa=blockIdx.x * blockDim.x + threadIdx.x; isa<nsa; isa+=step){/*ng is nsa*2. */
	    float tmpx=PTT[isa][0]*g[isa];
	    float tmpy=PTT[isa][1]*g[isa];
	    gx[threadIdx.x]+=tmpx+PTT[isa+nsa][0]*g[isa+nsa];
	    gy[threadIdx.x]+=tmpy+PTT[isa+nsa][1]*g[isa+nsa];
	}
	for(int step=(DIM_GP>>1); step>0; step>>=1){
	    __syncthreads();
	    if(threadIdx.x<step){
		gx[threadIdx.x]+=gx[threadIdx.x+step];
		gy[threadIdx.x]+=gy[threadIdx.x+step];
	    }
	}
	if(threadIdx.x==0){
	    atomicAdd(&ttout[iwfs*2], -gx[0]);
	    atomicAdd(&ttout[iwfs*2+1], -gy[0]);
	}
    }
    if(datai->PDF && ptt){
	float *restrict PDF=datai->PDF;//the diagonal block
	const int ipowfs=datai->ipowfs;
	float scale=-1./(datai->nwfs-1);//PDF*scale gives the off diagonal block
	for(int irow=0; irow<nwfs; irow++){//skip first row
	    if(data[irow].ipowfs!=ipowfs || data[irow].jwfs==0) continue;//different group or first wfs.
	    gdf[threadIdx.x]=0;
	    for(int isa=blockIdx.x * blockDim.x + threadIdx.x; isa<nsa; isa+=step){/*ng is nsa*2. */
		gdf[threadIdx.x]+=PDF[isa]*g[isa]+PDF[isa+nsa]*g[isa+nsa];
	    }
	    for(int step=(DIM_GP>>1); step>0; step>>=1){
		__syncthreads();
		if(threadIdx.x<step){
		    gdf[threadIdx.x]+=gdf[threadIdx.x+step];
		}
	    }
	    if(threadIdx.x==0){
		if(datai->PTT){//adjust by TT coupling
		    gdf[0]-=datai->PDFTT[0]*gx[0]+datai->PDFTT[1]*gy[0];
		}
		if(irow!=iwfs){
		    gdf[0]*=scale;
		}
		atomicAdd(&dfout[irow], -gdf[0]);
	    }
	}
    }
}
__global__ static void gpu_gpt_do(GPU_GP_T *data, float **wfsopd, float *ttin, float *dfin, float **gin, int ptt){
    const int iwfs=blockIdx.z;
    const int nwfs=gridDim.z;
    GPU_GP_T *datai=data+iwfs;
    const int pos=datai->pos;
    if(!pos) return;
    const int step=blockDim.x * gridDim.x;
    const int nsa=datai->nsa;
    int (*saptr)[2]=datai->saptr;
    const float (*restrict neai)[3]=datai->neai;
    float dxp=datai->dxp;
    float oxp=datai->oxp;
    float oyp=datai->oyp;
    float focus=0;
    if(datai->PDF && ptt){
	if(iwfs==0){
	    for(int id=1; id<nwfs; id++){
		focus+=dfin[id];
	    }
	}else{
	    focus=-dfin[iwfs];
	}
    }
    const float *restrict g=gin[iwfs];
    float *restrict map=wfsopd[iwfs];
    float GPscale=datai->GPscale;
    const short2 *restrict pxy=datai->GPp;
    float ttx=0, tty=0;
    if(datai->PTT && ptt){
	ttx=ttin[iwfs*2+0];
	tty=ttin[iwfs*2+1];
    }
    const int nx=datai->nxp;
    if(pos==1){
	for(int isa=blockIdx.x * blockDim.x + threadIdx.x; isa<nsa; isa+=step){
	    int ix=saptr[isa][0];
	    int iy=saptr[isa][1];
	    float cx=neai[isa][0];
	    float cy=neai[isa][1];
	    float cxy=neai[isa][2];
	    float gx=g[isa    ]+ttx+focus*(ix*dxp+oxp);
	    float gy=g[isa+nsa]+tty+focus*(iy*dxp+oyp);
	    float tmp=cxy*gx;//save data
	    gx=GPscale*(cx*gx+cxy*gy);
	    gy=GPscale*(tmp+cy*gy);
	    const short2 *restrict pxy2=pxy+isa*4;
	    atomicAdd(&map[iy    *nx+ix],   gx*pxy2[0].x + gy*pxy2[0].y);
	    atomicAdd(&map[iy    *nx+ix+1], gx*pxy2[1].x + gy*pxy2[1].y);
	    atomicAdd(&map[(iy+1)*nx+ix],   gx*pxy2[2].x + gy*pxy2[2].y);
	    atomicAdd(&map[(iy+1)*nx+ix+1], gx*pxy2[3].x + gy*pxy2[3].y);
	}
    }else if(pos==2){
	for(int isa=blockIdx.x * blockDim.x + threadIdx.x; isa<nsa; isa+=step){
	    int ix=saptr[isa][0];
	    int iy=saptr[isa][1];
	    float cx=neai[isa][0];
	    float cy=neai[isa][1];
	    float cxy=neai[isa][2];
	    float gx=g[isa    ]+ttx+focus*(ix*dxp+oxp);
	    float gy=g[isa+nsa]+tty+focus*(iy*dxp+oyp);
	    float tmp=cxy*gx;
	    gx=GPscale*(cx*gx+cxy*gy);
	    gy=GPscale*(tmp+cy*gy);
	    const short2 *restrict pxy2=pxy+isa*9;
	    atomicAdd(&map[iy    *nx+ix],   gx*pxy2[0].x + gy*pxy2[0].y);
	    atomicAdd(&map[iy    *nx+ix+1], gx*pxy2[1].x + gy*pxy2[1].y);
	    atomicAdd(&map[iy    *nx+ix+2], gx*pxy2[2].x + gy*pxy2[2].y);
	    atomicAdd(&map[(iy+1)*nx+ix],   gx*pxy2[3].x + gy*pxy2[3].y);
	    atomicAdd(&map[(iy+1)*nx+ix+1], gx*pxy2[4].x + gy*pxy2[4].y);
	    atomicAdd(&map[(iy+1)*nx+ix+2], gx*pxy2[5].x + gy*pxy2[5].y);
	    atomicAdd(&map[(iy+2)*nx+ix],   gx*pxy2[6].x + gy*pxy2[6].y);
	    atomicAdd(&map[(iy+2)*nx+ix+1], gx*pxy2[7].x + gy*pxy2[7].y);
	    atomicAdd(&map[(iy+2)*nx+ix+2], gx*pxy2[8].x + gy*pxy2[8].y);
	}
    }
}
/*Only do tt, NEA. do not to GP'*/
__global__ static void gpu_nea_do(GPU_GP_T *data, float *ttin, float *dfin, float **gin){
    const int iwfs=blockIdx.z;
    const int nwfs=gridDim.z;
    GPU_GP_T *datai=data+iwfs;
    const int pos=datai->pos;
    if(!pos) return;
    const int step=blockDim.x * gridDim.x;
    const int nsa=datai->nsa;
    int (*saptr)[2]=datai->saptr;
    const float (*restrict neai)[3]=datai->neai;
    float dxp=datai->dxp;
    float oxp=datai->oxp;
    float oyp=datai->oyp;
    float focus=0;
    if(datai->PDF){
	if(iwfs==0){
	    for(int id=1; id<nwfs; id++){
		focus+=dfin[id];
	    }
	}else{
	    focus=-dfin[iwfs];
	}
    }
    float *restrict g=gin[iwfs];
    float ttx=ttin[iwfs*2+0];
    float tty=ttin[iwfs*2+1];
    for(int isa=blockIdx.x * blockDim.x + threadIdx.x; isa<nsa; isa+=step){
	int ix=saptr[isa][0];
	int iy=saptr[isa][1];
	float cx=neai[isa][0];
	float cy=neai[isa][1];
	float cxy=neai[isa][2];
	float gx=g[isa    ]+ttx+focus*(ix*dxp+oxp);
	float gy=g[isa+nsa]+tty+focus*(iy*dxp+oyp);
	g[isa]=cx*gx+cxy*gy;
	g[isa+nsa]=cxy*gx+cy*gy;
    }
}

/*
  Tomography right hand size matrix. Computes xout = xout *beta + alpha * Hx' G' C * xin.
  xout is zeroed out before accumulation.
*/
void cutomo_grid::R(curcell **xout, float beta, const curcell *grad, float alpha, stream_t &stream){
    if(!*xout){
	*xout=curcellnew(grid->npsr, 1, grid->xnx, grid->xny);
    }else{
	curcellscale(*xout, beta, stream);
    }
    curcell *opdx=*xout;
    curzero(opdwfs->m, stream);
    curzero(ttf, stream);
    //just does low rank terms
    gpu_gp_do<<<dim3(24,1,nwfs), dim3(DIM_GP,1), 0, stream>>>
	(gpdata, grad->pm, ttf->p, ttf->p+nwfs*2, NULL, 1);
    gpu_gpt_do<<<dim3(24,1,nwfs), dim3(DIM_GP,1), 0, stream>>>
	(gpdata, opdwfs->pm, ttf->p, ttf->p+nwfs*2, grad->pm, 1);
    gpu_prop_grid_do<<<dim3(3,3,grid->npsr), dim3(16,16), 0, stream>>>
    	(hxdata, opdwfs->pm, opdx->pm, nwfs, grid->npsr, alpha, NULL, 't');
}

void cutomo_grid::Rt(curcell **gout, float beta, const curcell *xin, float alpha, stream_t &stream){
    if(!*gout){
	*gout=curcellnew(nwfs, 1, grid->ngrad, (long*)NULL);
    }else{
	curcellscale(*gout, beta, stream);
    }
    curcell *grad=*gout;
    curzero(opdwfs->m, stream);
    gpu_prop_grid_do<<<dim3(3,3, nwfs), dim3(16,16), 0, stream>>>
	(hxdata, opdwfs->pm, xin->pm, nwfs, grid->npsr, alpha, NULL, 'n');
    curzero(ttf, stream);
    gpu_gp_do<<<dim3(24,1,nwfs), dim3(DIM_GP,1), 0, stream>>>
	(gpdata, grad->pm, ttf->p, ttf->p+nwfs*2, opdwfs->pm, 1);
    gpu_nea_do<<<dim3(24,1,nwfs), dim3(DIM_GP,1), 0, stream>>>
	(gpdata, ttf->p, ttf->p+nwfs*2, grad->pm);
}

/*
  Tomography left hand size matrix. Computes xout = beta*xout + alpha * Hx' G' C Gp Hx * xin.
  xout is zeroed out before accumulation.

  Be Careful about big kernels 
  1) When a kernel thread reads locations written by other threads, synchronization may be required unless gauranteed in the same wrap.
  2) When a kernel thread writes to locations that may be written by others, atomic operations are required unless gauranteed in the same wrap.
*/
void cutomo_grid::L(curcell **xout, float beta, const curcell *xin, float alpha, stream_t &stream){
    const int npsr=grid->npsr;
    if(!*xout){
	*xout=curcellnew(grid->npsr, 1, grid->xnx, grid->xny);
    }else{
	curscale((*xout)->m, beta, stream);
    }
#if TIMING==2
    EVENT_INIT(6);
#define RECORD(i) EVENT_TIC(i)
#else
#define RECORD(i)
#endif
    RECORD(0);
    curcell *opdx=*xout;
    curzero(opdwfs->m, stream);
    //xin to opdwfs
    gpu_prop_grid_do<<<dim3(3,3, nwfs), dim3(16,16), 0, stream>>>
    	(hxdata, opdwfs->pm, xin->pm, nwfs, npsr, 1, NULL, 'n');
    RECORD(1);
    curzero(ttf, stream);
    //opdwfs to grad to ttf
    gpu_gp_do<<<dim3(24,1,nwfs), dim3(DIM_GP,1), 0, stream>>>
	(gpdata, grad->pm, ttf->p, ttf->p+nwfs*2, opdwfs->pm, ptt);
    RECORD(2);
    curzero(opdwfs->m, stream);
    //grad and ttf to opdwfs
    gpu_gpt_do<<<dim3(24,1,nwfs), dim3(DIM_GP,1), 0, stream>>>
	(gpdata, opdwfs->pm, ttf->p, ttf->p+nwfs*2, grad->pm, ptt);
    RECORD(3);
    //opdwfs to opdx
    gpu_prop_grid_do<<<dim3(3,3,grid->npsr), dim3(16,16), 0, stream>>>
    	(hxdata, opdwfs->pm, opdx->pm, nwfs, grid->npsr, alpha, NULL,'t');
    RECORD(4);
    /*This could be in parallel to above ones. */
    gpu_laplacian_do<<<dim3(3,3,grid->npsr),dim3(16,16), 0, stream>>>
	(hxdata, opdx->pm, xin->pm, nwfs, alpha);
    RECORD(5);
#if TIMING==2
    EVENT_TOC;
    info2("TomoL: Hx %.0f, Gp %.0f, Gpt %.0f, Hxt %.0f, L2 %.0f\n", 
	  times[1], times[2], times[3], times[4], times[5]);
    EVENT_DEINIT;
#endif
    //overhead of TomoL 27 micro-seconds (timing without synchornization).
}


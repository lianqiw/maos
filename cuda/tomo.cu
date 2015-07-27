/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "utils.h"
#include "accphi.h"
#include "cucmat.h"
#include "curmat.h"
#include "cudata.h"
#include "tomo.h"
namespace cuda_recon{
/*
  Prepare data for GP operation in GPU
  We convert GP from sparse matrix to coefficients of pos==1 or 2. By doing so, we can also convert the type to integer to save memory transfer. Short2 is used without performance penalty.
  For partially illiminated subapertures, this conversion brings in approximation as we limite the influence of gradient from each subaperture to (pos*2-1)x(pos*2-1) points. The system performance is largely not affected.
*/
void prep_GP(cumat<short2> &GPp, Real *GPscale, cusp &GPf,
	     const dsp *GP, const loc_t *saloc, const loc_t *ploc){
    if(!GP){ 
	error("GP is required\n");
    }
    int pos=(int)round(saloc->dx/ploc->dx);
    if(pos==1 || pos==2){//normally true
	dsp *GPt=dsptrans(GP);
	const spint *pp=GPt->p;
	const spint *pi=GPt->i;
	const double *px=GPt->x;
	//convert the max Real to max 2 byte integer
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
		if((zx<0 || zx>zmax || zy<0 || zy>zmax) & fabs(px[ir])>1e-7){
		    warning("isa=%d, G(%d, %d)=%g\n", isa, zx, zy, px[ir]);
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
	dspfree(GPt);
	GPp=cumat<short2>(np, nsa);
	cudaMemcpy(GPp.P(), partxy, sizeof(short2)*np*nsa, cudaMemcpyHostToDevice);
	*GPscale=1./pxscale;
	free(partxy);
    }else{/*use sparse */
	GPf=cusp(GP, 1);
    }
}
static cumat<int> 
prep_saptr(loc_t *saloc, map_t *pmap){
    /*saloc mapped onto pmap*/
    int nsa=saloc->nloc;
    int (*saptr)[2]=new int[nsa][2];
    const Real dx1=1./pmap->dx;
    const Real dy1=1./pmap->dy;
    const Real ox=pmap->ox;
    const Real oy=pmap->oy;
    const double *salocx=saloc->locx;
    const double *salocy=saloc->locy;
    for(int isa=0; isa<nsa; isa++){
	saptr[isa][0]=(int)roundf((salocx[isa]-ox)*dx1);
	saptr[isa][1]=(int)roundf((salocy[isa]-oy)*dy1);
    }
    cumat<int> saptr_gpu=cumat<int>(2, nsa);
    DO(cudaMemcpy(saptr_gpu.P(), saptr, nsa*2*sizeof(int), cudaMemcpyHostToDevice));
    delete [] saptr;
    return saptr_gpu;
}
static curmat convert_neai(dsp *nea){
    spint *pp=nea->p;
    spint *pi=nea->i;
    double *px=nea->x;
    int nsa=nea->ny/2;
    Real (*neai)[3]=(Real(*)[3])calloc(3*nsa, sizeof(Real));
    for(int ic=0; ic<nea->n; ic++){
	for(spint ir=pp[ic]; ir<pp[ic+1]; ir++){
	    int ix=pi[ir];
	    int isa=ic<nsa?ic:ic-nsa;
	    Real val=(Real)px[ir];
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
    curmat neai_gpu=curmat(3, nsa);
    DO(cudaMemcpy(neai_gpu.P(), neai, 3*nsa*sizeof(Real), cudaMemcpyHostToDevice));
    free(neai);
    return neai_gpu;
}
#define TIMING 0
void cutomo_grid::init_hx(const PARMS_T *parms, const RECON_T *recon){
    dir_t dir[nwfs];
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	const int ipowfs = parms->wfsr[iwfs].powfs;
	dir[iwfs].skip=parms->powfs[ipowfs].skip;
	dir[iwfs].hs=parms->wfsr[iwfs].hs;
	dir[iwfs].thetax=parms->wfsr[iwfs].thetax;
	dir[iwfs].thetay=parms->wfsr[iwfs].thetay;
    }
    Real dt=parms->tomo.predict?parms->sim.dt*2:0;
    if(dt>0){
	info("dt=%g. vx[0]=%g\n", dt, grid->xmap[0].vx);
    }
    hx.Init_l2d(grid->pmap, dir, nwfs, grid->xmap, dt);

    LAP_T lapc[recon->npsr];
    for(int ips=0; ips<recon->npsr; ips++){ 
	Real tmp=laplacian_coef(recon->r0, recon->wt->p[ips], recon->xmap->p[ips]->dx)*0.25f;
	lapc[ips].nxps=recon->xmap->p[ips]->nx;
	lapc[ips].nyps=recon->xmap->p[ips]->ny;
	lapc[ips].l2c=tmp*tmp*parms->tomo.cxxscale*TOMOSCALE;
	if(parms->tomo.piston_cr){
	    lapc[ips].zzi=loccenter(recon->xloc->p[ips]);
	    lapc[ips].zzv=lapc[ips].l2c*1e-6;
	}else{
	    lapc[ips].zzi=-1;
	}
    }
    lap=cumat<LAP_T>(recon->npsr, 1);
    cudaMemcpy(lap.P(), lapc, sizeof(LAP_T)*recon->npsr, cudaMemcpyHostToDevice);
}

cutomo_grid::cutomo_grid(const PARMS_T *parms, const RECON_T *recon,
			 const POWFS_T *powfs, curecon_geom *_grid)
    :cucg_t(parms?parms->tomo.maxit:0, parms?parms->recon.warm_restart:0),
     grid(_grid), GPscale(0),ptt(0),nwfs(0){
    nwfs=parms->nwfsr;

    if(recon->PTT && !PTT){//for t/t proj in 1)uplink t/t 2) recon
	cp2gpu(PTT, recon->PTT);
    }
    ptt=!parms->recon.split || (parms->tomo.splitlrt && parms->recon.mvm!=2); 
    {
	PDF=curcell(parms->npowfs, 1);
	PDFTT=curcell(parms->npowfs, 1);
	dcell *pdftt=NULL;
	dcellmm(&pdftt, recon->PDF, recon->TT, "nn", 1);
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    if(parms->powfs[ipowfs].dfrs){//deprecated.	
		/*We only use the first diagonal block for each powfs. The
		  off diagonal is simply -1/(nlgswfs-1) times the diagonal block*/
		int iwfs1=parms->powfs[ipowfs].wfs->p[1];//second wfs
		cp2gpu(PDF[ipowfs], recon->PDF->p[iwfs1*nwfs+iwfs1]);
		if(parms->powfs[ipowfs].trs){
		    /*coupling between TT and DF modes. 
		      We desire (I-DF*PDF)(I-TT*PTT)g=(I-TT*PTT-DF*PDF+DF*PDF*TT*PTT)g
		      So we first compute tt=PTT*g; df=PDF*g; then
		      g2=(I-TT*tt-DF*(df-(PDF*TT)*tt))
		      Here we record the values of PDF*TT
		    */
		    cp2gpu(PDFTT[ipowfs], pdftt->p[iwfs1*nwfs+iwfs1]);
		}
	    }
	}
	dcellfree(pdftt);
    }
    {
	const int npowfs=parms->npowfs;
	GPp=cucell<short2>(nwfs, 1);
	GP=cuspcell(nwfs, 1);
	GPscale=new Real[npowfs];
	saptr=cucell<int>(npowfs, 1);
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    if(parms->powfs[ipowfs].skip) continue;
	    cusp GPf;
	    cumat<short2> GPpi;
	    prep_GP(GPpi, &GPscale[ipowfs], GPf,
		    recon->GP->p[ipowfs], powfs[ipowfs].saloc, recon->ploc);
	    saptr[ipowfs]=prep_saptr(powfs[ipowfs].saloc, recon->pmap);
	    for(int iwfs=0; iwfs<nwfs; iwfs++){
		if(ipowfs==parms->wfsr[iwfs].powfs){
		    if(GPpi){
			GPp[iwfs]=GPpi;
		    }
		    if(GPf){
			GP[iwfs]=GPf;
		    }
		}
	    }
	}
    }

    {
	neai=curcell(parms->nwfsr, 1);
	/*convert recon->saneai to gpu. */
	for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	    int ipowfs=parms->wfsr[iwfs].powfs;
	    int iwfs0=parms->recon.glao?iwfs:parms->powfs[ipowfs].wfs->p[0];/*first wfs in this group. */
	    if(iwfs!=iwfs0 && recon->saneai->p[iwfs+iwfs*parms->nwfsr]->p
	       ==recon->saneai->p[iwfs0+iwfs0*parms->nwfsr]->p){
		neai[iwfs]=neai[iwfs0];
	    }else{
		dsp *nea=recon->saneai->p[iwfs+iwfs*parms->nwfsr];
		neai[iwfs]=convert_neai(nea);
	    }
	}/*for iwfs */
    }

    
    init_hx(parms, recon);
    
    {
	GPU_GP_T *GPDATA=new GPU_GP_T[nwfs];
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    const int ipowfs = parms->wfsr[iwfs].powfs;
	    if(parms->powfs[ipowfs].skip) continue;
	    GPDATA[iwfs].ipowfs=ipowfs;
	    GPDATA[iwfs].nwfs=parms->powfs[ipowfs].nwfsr;
	    GPDATA[iwfs].jwfs=parms->powfs[ipowfs].wfsind->p[iwfs];//wfs index in this group
	    GPDATA[iwfs].dsa=powfs[ipowfs].pts->dsa;
	    GPDATA[iwfs].pos=parms->tomo.pos;
	    GPDATA[iwfs].saptr=(int(*)[2])saptr[ipowfs].P();
	    GPDATA[iwfs].GPp=(short2*)GPp[iwfs].P();
	    GPDATA[iwfs].GPscale=GPscale[ipowfs];
	    if(parms->powfs[ipowfs].trs){
		GPDATA[iwfs].PTT=PTT(iwfs,iwfs).P();
	    }
	    if(parms->powfs[ipowfs].dfrs){
		GPDATA[iwfs].PDF=PDF[ipowfs].P();
		GPDATA[iwfs].PDFTT=PDFTT[ipowfs].P();
	    }

	    GPDATA[iwfs].neai=(const Real(*)[3])neai[iwfs].P();
	    GPDATA[iwfs].nsa=powfs[ipowfs].pts->nsa;
	    GPDATA[iwfs].nxp=recon->pmap->nx;
	    GPDATA[iwfs].dxp=recon->pmap->dx;
	    GPDATA[iwfs].dyp=recon->pmap->dy;
	    GPDATA[iwfs].oxp=recon->pmap->ox;
	    GPDATA[iwfs].oyp=recon->pmap->oy;
	}
	gpdata=cumat<GPU_GP_T>(nwfs,1);
	DO(cudaMemcpy(gpdata.P(), GPDATA, sizeof(GPU_GP_T)*nwfs, cudaMemcpyHostToDevice));
	delete [] GPDATA;
    }

    if(parms->tomo.precond==1){
	precond=new cufdpcg_t(recon->fdpcg, grid);
    }

    {/**Initialize run time data*/
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

	opdwfs=curcell(nwfs, 1, nxpw, nypw);
	grad=curcell(nwfs, 1, ngw, (int*)NULL);
	ttf=curmat(3*nwfs, 1);
    }
}

/*
  If merge the operation in to gpu_prop_grid_do, need to do atomic
  operation because previous retracing starts at offo, not 0.  */
__global__ void gpu_laplacian_do(LAP_T *datai, Real **outall, Real **inall, int nwfs, Real alpha){
    Real *restrict in, *restrict out;
    
    const int ips=blockIdx.z;
    in=inall[ips];
    out=outall[ips];
    datai+=ips;
    
    const int nx=datai->nxps;
    const int ny=datai->nyps;
    const Real alpha2=datai->l2c*alpha;
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

/**
   The third grid dimension tells the wfs to handle. 
   Handles TTDF*GP
*/
#define DIM_GP 128
__global__ static void gpu_gp_do(GPU_GP_T *data, Real **gout, Real *ttout, Real *dfout, Real **wfsopd, int ptt){
    __shared__ Real gx[DIM_GP];
    __shared__ Real gy[DIM_GP];
    __shared__ Real gdf[DIM_GP];
    const int iwfs=blockIdx.z;
    const int nwfs=gridDim.z;
    GPU_GP_T *datai=data+iwfs;
    const int pos=datai->pos;
    if(!pos) return;
    const int nsa=datai->nsa;
    const int step=blockDim.x * gridDim.x;
    int (*restrict saptr)[2]=datai->saptr;
    Real *restrict g=gout[iwfs];
    Real GPscale=datai->GPscale;
    if(datai->GPp && wfsopd){
	const Real *restrict map=wfsopd[iwfs];
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
    if(datai->PTT && ptt){ 
	Real (*restrict PTT)[2]=(Real(*)[2])datai->PTT;
	gx[threadIdx.x]=0;
	gy[threadIdx.x]=0;
	for(int isa=blockIdx.x * blockDim.x + threadIdx.x; isa<nsa; isa+=step){/*ng is nsa*2. */
	    Real tmpx=PTT[isa][0]*g[isa];
	    Real tmpy=PTT[isa][1]*g[isa];
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
	Real *restrict PDF=datai->PDF;//the diagonal block
	const int ipowfs=datai->ipowfs;
	Real scale=-1./(datai->nwfs-1);//PDF*scale gives the off diagonal block
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
/**
   Todo: Improve by removing atomic operation.
   Handles GP'*nea*(1-TTDF)
   Be carefulll about the ptt flag. It is always 1 for Right hand side, but may be zero for Left hand side.
*/
__global__ static void gpu_gpt_do(GPU_GP_T *data, Real **wfsopd, Real *ttin, Real *dfin, Real **gin, int ptt){
    const int iwfs=blockIdx.z;
    const int nwfs=gridDim.z;
    GPU_GP_T *datai=data+iwfs;
    const int pos=datai->pos;
    if(!pos) return;
    const int step=blockDim.x * gridDim.x;
    const int nsa=datai->nsa;
    int (*const saptr)[2]=datai->saptr;
    const Real (*restrict neai)[3]=datai->neai;
    const Real dxp=datai->dxp;
    const Real oxp=datai->oxp;
    const Real oyp=datai->oyp;
    Real focus=0;
    if(datai->PDF && ptt){
	if(iwfs==0){
	    for(int id=1; id<nwfs; id++){
		focus+=dfin[id];
	    }
	}else{
	    focus=-dfin[iwfs];
	}
    }
    Real *restrict g=gin[iwfs];
    Real *restrict map=wfsopd?wfsopd[iwfs]:0;
    Real GPscale=datai->GPscale;
    Real ttx=0, tty=0;
    if(datai->PTT && ptt){
	ttx=ttin[iwfs*2+0];
	tty=ttin[iwfs*2+1];
    }
    const int nx=datai->nxp;
    const short2 *restrict pxy=datai->GPp;

    if(!pxy || !map){
	for(int isa=blockIdx.x * blockDim.x + threadIdx.x; isa<nsa; isa+=step){
	    int ix=saptr[isa][0];
	    int iy=saptr[isa][1];
	    Real cx=neai[isa][0];
	    Real cy=neai[isa][1];
	    Real cxy=neai[isa][2];
	    Real gx=g[isa    ]+ttx+focus*(ix*dxp+oxp);
	    Real gy=g[isa+nsa]+tty+focus*(iy*dxp+oyp);
	    g[isa]=cx*gx+cxy*gy;
	    g[isa+nsa]=cxy*gx+cy*gy;
	}
    }else if(pos==1){
	for(int isa=blockIdx.x * blockDim.x + threadIdx.x; isa<nsa; isa+=step){
	    int ix=saptr[isa][0];
	    int iy=saptr[isa][1];
	    Real cx=neai[isa][0];
	    Real cy=neai[isa][1];
	    Real cxy=neai[isa][2];
	    Real gx=g[isa    ]+ttx+focus*(ix*dxp+oxp);
	    Real gy=g[isa+nsa]+tty+focus*(iy*dxp+oyp);
	    Real tmp=cxy*gx;//save data
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
	    Real cx=neai[isa][0];
	    Real cy=neai[isa][1];
	    Real cxy=neai[isa][2];
	    Real gx=g[isa    ]+ttx+focus*(ix*dxp+oxp);
	    Real gy=g[isa+nsa]+tty+focus*(iy*dxp+oyp);
	    Real tmp=cxy*gx;
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

void cutomo_grid::do_gp(curcell &_grad, const curcell &_opdwfs, int ptt2, stream_t &stream){
    if(_opdwfs){
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    if(GPp[iwfs] || !GP[iwfs]) continue;
	    cuzero(_grad[iwfs], stream);
	    cuspmul(_grad[iwfs].P(), GP[iwfs], _opdwfs[iwfs].P(), 1, 'n', 1, stream);
	}
    }
    cuzero(ttf, stream);
    gpu_gp_do<<<dim3(24,1,nwfs), dim3(DIM_GP,1), 0, stream>>>
	(gpdata.P(), _grad.pm, ttf.P(), ttf.P()+nwfs*2, _opdwfs?_opdwfs.pm:NULL, ptt2);
}
void cutomo_grid::do_gpt(curcell &_opdwfs, curcell &_grad, int ptt2, stream_t &stream){
    if(_opdwfs){
	cuzero(_opdwfs.M(), stream);
    }
    //Does  GP'*NEA*(1-TTDF) if _opdwfs!=0 and GPp!=0 or NEA*(1-TTDF)
    gpu_gpt_do<<<dim3(24,1,nwfs), dim3(DIM_GP,1), 0, stream>>>
	(gpdata.P(), _opdwfs?_opdwfs.pm:0, ttf.P(), ttf.P()+nwfs*2, _grad.pm, ptt2);
    
    if(_opdwfs){//Does GP' for GP with sparse
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    if(GPp[iwfs] || !GP[iwfs]) continue;
	    cuzero(_opdwfs[iwfs], stream);
	    cuspmul(_opdwfs[iwfs].P(), GP[iwfs], _grad[iwfs].P(), 1, 't', 1, stream);
	}
    }
}
/*
  Tomography right hand side matrix. Computes xout = xout *beta + alpha * Hx' G' C * xin.
  xout is zeroed out before accumulation.
*/
void cutomo_grid::R(curcell &xout, Real beta,  curcell &_grad, Real alpha, stream_t &stream){
    if(!xout){
	xout=curcell(grid->npsr, 1, grid->xnx, grid->xny);
    }else{
	curcellscale(xout, beta, stream);
    }
    do_gp(_grad, curcell(), 1, stream);
    do_gpt(opdwfs, _grad, 1, stream);
    hx.backward(opdwfs.pm, xout.pm, alpha, NULL, stream);
}

void cutomo_grid::Rt(curcell &gout, Real beta,  curcell &xin, Real alpha, stream_t &stream){
    if(!gout){
	gout=curcell(nwfs, 1, grid->ngrad, (long*)NULL);
    }else{
	curcellscale(gout, beta, stream);
    }
    opdwfs.M().zero(stream);
    hx.forward(opdwfs.pm, xin.pm, alpha, NULL, stream);
    do_gp(gout, opdwfs, 1, stream);
    curcell dummy;
    do_gpt(dummy, gout, 1, stream);
}

/*
  Tomography left hand size matrix. Computes xout = beta*xout + alpha * Hx' G' C Gp Hx * xin.
  xout is zeroed out before accumulation.

  Be Careful about big kernels 
  1) When a kernel thread reads locations written by other threads, synchronization may be required unless gauranteed in the same wrap.
  2) When a kernel thread writes to locations that may be written by others, atomic operations are required unless gauranteed in the same wrap.
*/
void cutomo_grid::L(curcell &xout, Real beta, const curcell &xin, Real alpha, stream_t &stream){
    if(!xout){
	xout=curcell(grid->npsr, 1, grid->xnx, grid->xny);
    }else{
	curscale(xout.M(), beta, stream);
    }
#if TIMING==2
    EVENT_INIT(6);
#define RECORD(i) EVENT_TIC(i)
#else
#define RECORD(i)
#endif
    RECORD(0);
    cuzero(opdwfs.M(), stream);
    //xin to opdwfs
    hx.forward(opdwfs.pm, xin.pm, 1, NULL, stream);
    RECORD(1);
    //opdwfs to grad to ttf
    do_gp(grad, opdwfs, ptt, stream);
    RECORD(2);
    //grad and ttf to opdwfs
    do_gpt(opdwfs, grad, ptt, stream);
    RECORD(3);
    //opdwfs to xout
    hx.backward(opdwfs.pm, xout.pm, alpha, NULL, stream);
    RECORD(4);
    /*This could be in parallel to hx->forward, do_gp, do_gpt*/
    gpu_laplacian_do<<<dim3(3,3,grid->npsr),dim3(16,16), 0, stream>>>
	(lap.P(), xout.pm, xin.pm, nwfs, alpha);
    RECORD(5);
#if TIMING==2
    EVENT_TOC;
    info2("TomoL: Hx %.0f, Gp %.0f, Gpt %.0f, Hxt %.0f, L2 %.0f\n", 
	  times[1], times[2], times[3], times[4], times[5]);
    EVENT_DEINIT;
#endif
    //overhead of TomoL is 27 micro-seconds (timing without synchornization).

}
}//namespace

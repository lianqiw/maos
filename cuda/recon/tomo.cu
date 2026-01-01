/*
  Copyright 2009-2026 Lianqi Wang
  
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
#define TIMING 0
#define SAVE_TOMO 0
#include "../sim/accphi.h"
#include "../sim/cudata.h"
#include "tomo.h"

/*
  Prepare data for GP operation in GPU
  We convert GP from sparse matrix to coefficients of pos==1 or 2. By doing so, we can also convert the type to integer to save memory transfer. Short2 is used without performance penalty.
  For partially illiminated subapertures, this conversion brings in approximation as we limite the influence of gradient from each subaperture to (pos*2-1)x(pos*2-1) points. The system performance is largely not affected.
*/
void prep_GP(Array<short2, Gpu>& GPp, Real* GPscale, cusp& GPf,
	const dsp* GP, const loc_t* saloc, const loc_t* ploc){
	if(!GP){
		error("GP is required\n");
	} else if(GP->nx!=saloc->nloc*2||GP->ny!=ploc->nloc){
		error("GP has wrong dimension: %ldx%ld. (%ldx%ld is expected).\n",
			GP->nx, GP->ny, saloc->nloc*2, ploc->nloc);
	}
	real pos=saloc->dx/ploc->dx;
	real xdiff=(ploc->locx[0]-saloc->locx[0])/ploc->dx;
	real ydiff=(ploc->locy[0]-saloc->locy[0])/ploc->dy;
	if((fabs(pos-1)<EPS||fabs(pos-2)<EPS) 
		&& fabs(xdiff-round(xdiff))<EPS && fabs(ydiff-round(ydiff))<EPS){//These well aligned cases are accelerated with matrix-free approach.
		dbg("GP uses matrix-free approach: pos=%g, xdiff=%g, ydiff=%g.\n", pos, xdiff, ydiff);
		dsp* GPt=dsptrans(GP);
		const spint* pp=GPt->pp;
		const spint* pi=GPt->pi;
		const real* px=GPt->px;
		//convert the max Real to max 2 byte integer
		const real pxscale=floor(32767./dvecmaxabs(px, GPt->nzmax));
		const int zmax=(int)round(pos);
		const int np1=zmax+1;
		const int np=np1*np1;
		int nsa=saloc->nloc;
		short2* partxy=mycalloc(np*nsa, short2);//need to zero memory
		const real dx1=1./ploc->dx;
		const real dy1=1./ploc->dy;
		for(int ic=0; ic<GPt->ny; ic++){
			int isa=(ic<nsa)?ic:(ic-nsa);
			for(spint ir=pp[ic]; ir<pp[ic+1]; ir++){
				const int ix=pi[ir];
				const real lx=ploc->locx[ix];
				const real ly=ploc->locy[ix];
				const real sx=saloc->locx[isa];
				const real sy=saloc->locy[isa];
				int zx=(int)round((lx-sx)*dx1);
				int zy=(int)round((ly-sy)*dy1);
				/**
				   When the points used to generate GP align well
				   with the subaperture edge, the coupled points are
				   confined within the subaperture.
				*/
				if((zx<0||zx>zmax||zy<0||zy>zmax)&&fabs(px[ir])>EPS){
					warning("isa=%d, G(%d, %d)=%g is not zero. ploc is at (%g,%g), saloc is at (%g, %g)\n", isa, zx, zy, px[ir], lx, ly, sx, sy);
				}
				if(zx<0) zx=0;
				if(zx>zmax) zx=zmax;
				if(zy<0) zy=0;
				if(zy>zmax) zy=zmax;
				if(ic<nsa){
					partxy[np*isa+zx+zy*np1].x+=(short)round(px[ir]*pxscale);
				} else{
					partxy[np*isa+zx+zy*np1].y+=(short)round(px[ir]*pxscale);
				}
			}
		}
		dspfree(GPt);
		GPp.init(np, nsa);
		DO(cudaMemcpy(GPp(), partxy, sizeof(short2)*np*nsa, H2D));
		*GPscale=1./pxscale;
		free(partxy);
	} else{/*use sparse */
		dbg("GP uses sparse matrix: pos=%g, xdiff=%g, ydiff=%g.\n", pos, xdiff, ydiff);
		GPf=cusp(GP);
	}
}
static void
prep_saptr(cuimat& saptr_gpu, loc_t* saloc, map_t* pmap){
	/*saloc mapped onto pmap*/
	int nsa=saloc->nloc;
	int(*saptr)[2]=new int[nsa][2];
	const Real dx1=1./pmap->dx;
	const Real dy1=1./pmap->dy;
	const Real ox=pmap->ox;
	const Real oy=pmap->oy;
	const real* salocx=saloc->locx;
	const real* salocy=saloc->locy;
	for(int isa=0; isa<nsa; isa++){
		saptr[isa][0]=(int)roundf((salocx[isa]-ox)*dx1);
		saptr[isa][1]=(int)roundf((salocy[isa]-oy)*dy1);
	}
	saptr_gpu.init(2, nsa);
	DO(cudaMemcpy(saptr_gpu, saptr, nsa*2*sizeof(int), H2D));
	delete[] saptr;
}
static curmat convert_neai(dsp* nea){
	spint* pp=nea->pp;
	spint* pi=nea->pi;
	real* px=nea->px;
	int nsa=nea->ny/2;
	Real(*neai)[3]=(Real(*)[3])calloc(3*nsa, sizeof(Real));
	for(int ic=0; ic<nea->ny; ic++){
		for(spint ir=pp[ic]; ir<pp[ic+1]; ir++){
			int ix=pi[ir];
			int isa=ic<nsa?ic:ic-nsa;
			Real val=(Real)px[ir];
			if(ix==ic){/*diagonal part. */
				if(ic==isa){/*x */
					neai[isa][0]=val;
				} else{/*y */
					neai[isa][1]=val;
				}
			} else if(ix==ic-nsa||ix==ic+nsa){/*cross part. symmetric. */
				neai[isa][2]=val;
			} else{
				error("saneai has invalid format\n");
			}
		}
	}
	curmat neai_gpu=curmat(3, nsa);
	DO(cudaMemcpy(neai_gpu(), neai, 3*nsa*sizeof(Real), H2D));
	free(neai);
	return neai_gpu;
}

void cutomo_grid::init_hx(const parms_t* parms, const recon_t* recon){
	//dbg("create raytracing operator.\n");
	dir_t* dir=new dir_t[nwfs];
	for(int iwfs=0; iwfs<nwfs; iwfs++){
		const int ipowfs=parms->wfsr[iwfs].powfs;
		dir[iwfs].skip=parms->powfs[ipowfs].skip;
		dir[iwfs].hs=parms->wfsr[iwfs].hs;
		dir[iwfs].thetax=parms->wfsr[iwfs].thetax;
		dir[iwfs].thetay=parms->wfsr[iwfs].thetay;
		if(parms->powfs[ipowfs].type==WFS_SH){
			dir[iwfs].misregx=parms->wfsr[iwfs].misregx;
			dir[iwfs].misregy=parms->wfsr[iwfs].misregy;
		}
		if(parms->tomo.predict){
			dir[iwfs].delay=parms->sim.dt*(parms->powfs[ipowfs].dtrat+1+parms->sim.alhi);
		}
	}
	hx.init_l2d(grid->pmap, dir, nwfs, grid->xmap);
	delete[] dir;
	Array<lap_t>lapc(recon->npsr);
	for(int ips=0; ips<recon->npsr; ips++){
		Real tmp=laplacian_coef(recon->r0, recon->wt->p[ips], recon->xmap->p[ips]->dx)*0.25f;
		lapc[ips].nxps=recon->xmap->p[ips]->nx;
		lapc[ips].nyps=recon->xmap->p[ips]->ny;
		lapc[ips].l2c=tmp*tmp*parms->tomo.cxxscale*TOMOSCALE;
		if(parms->tomo.piston_cr){
			lapc[ips].zzi=loccenter(recon->xloc->p[ips]);
			lapc[ips].zzv=lapc[ips].l2c*1e-6;
		} else{
			lapc[ips].zzi=-1;
		}
	}
	lap.init(recon->npsr, 1);
	DO(cudaMemcpy(lap(), lapc, sizeof(lap_t)*recon->npsr, H2D));
}

cutomo_grid::cutomo_grid(const parms_t* parms, const recon_t* recon, const curecon_geom* _grid)
	:cusolve_cg(parms?parms->tomo.maxit:0, parms?parms->tomo.cgwarm:0),
	grid(_grid), rhs_nttf(0), lhs_nttf(0), lhs_skip0(0), nwfs(0){
	nwfs=parms->nwfsr;
	if(!parms->tomo.square){
		error("cutomo_grid requires parms->tomo.square=1.\n");
	}
	if(recon->PTTF&&!PTTF){//for t/t proj in 1)uplink t/t 2) recon
		cp2gpu(PTTF, recon->PTTF);
	}
	if(recon->PFF && !PFF){
		cp2gpu(PFF, recon->PFF);
	}
	{
		//dbg("Copying GP\n");
		GPp=NumCell<short2, Gpu>(nwfs, 1);
		GP=cuspcell(nwfs, 1);
		GPscale.init(nwfs, 1);
		saptr=NumCell<int, Gpu>(nwfs, 1);
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			const int ipowfs=parms->wfsr[iwfs].powfs;
			const int iwfs0=parms->powfs[ipowfs].wfsr->p[0];
			if(parms->powfs[ipowfs].skip) continue;
			if(iwfs==iwfs0||recon->GP->p[iwfs]->pp!=recon->GP->p[iwfs0]->pp){
				prep_GP(GPp[iwfs], &GPscale[iwfs], GP[iwfs], recon->GP->p[iwfs], recon->saloc->p[ipowfs], recon->ploc);
			} else{
				GPp[iwfs]=GPp[iwfs0];
				GPscale[iwfs]=GPscale[iwfs0];
				GP[iwfs]=GP[iwfs0];
			}
			if(iwfs==iwfs0){
				prep_saptr(saptr[iwfs], recon->saloc->p[ipowfs], recon->pmap);
			} else{
				saptr[iwfs]=saptr[iwfs0];
			}
		}
	}

	{
		//dbg("Copying neai\n");
		neai=curcell(parms->nwfsr, 1);
		/*convert recon->saneai to gpu. */
		for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
			int ipowfs=parms->wfsr[iwfs].powfs;
			if(parms->powfs[ipowfs].skip) continue;
			int iwfs0=parms->recon.glao?iwfs:parms->powfs[ipowfs].wfs->p[0];/*first wfs in this group. */
			if(iwfs!=iwfs0&&P(recon->saneai,iwfs,iwfs)->pp==P(recon->saneai,iwfs0,iwfs0)->pp){
				//dbg("reference neai from %d to %d\n", iwfs0, iwfs);
				neai[iwfs]=neai[iwfs0];
			} else{
				//dbg("Copy neai from cpu for %d\n", iwfs);
				neai[iwfs]=convert_neai(P(recon->saneai, iwfs, iwfs));
			}
		}/*for iwfs */
	}

	init_hx(parms, recon);

	{
		Array<gpu_gp_t> GPDATA(nwfs, 1);
		//gpu_gp_t *GPDATA=new gpu_gp_t[nwfs];
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			const int ipowfs=parms->wfsr[iwfs].powfs;
			if(parms->powfs[ipowfs].skip) continue;
			GPDATA[iwfs].ipowfs=ipowfs;
			GPDATA[iwfs].nwfs=parms->powfs[ipowfs].nwfsr;
			GPDATA[iwfs].jwfs=parms->powfs[ipowfs].wfsind->p[iwfs];//wfs index in this group
			GPDATA[iwfs].dsa=recon->saloc->p[ipowfs]->dx;
			GPDATA[iwfs].pos=parms->tomo.pos;
			GPDATA[iwfs].saptr=(int(*)[2])saptr[iwfs]();
			GPDATA[iwfs].GPp=(short2*)GPp[iwfs]();
			GPDATA[iwfs].GPscale=GPscale[iwfs];
			if(parms->powfs[ipowfs].trs || parms->powfs[ipowfs].frs){
				GPDATA[iwfs].PTTF=PTTF(iwfs, iwfs)();
				if(!rhs_nttf) rhs_nttf=PTTF(iwfs, iwfs).Nx();
				else if(rhs_nttf!=PTTF(iwfs, iwfs).Nx()){
					error("PTTF for different WFS does not match: %d vs %ld\n", rhs_nttf, PTTF(iwfs, iwfs).Nx());
				}
				
				if(parms->recon.split==1&&parms->tomo.splitlrt==1&&parms->powfs[ipowfs].frs){
					GPDATA[iwfs].PFF=PFF(iwfs, iwfs)();
					if(!lhs_nttf) lhs_nttf=PFF(iwfs, iwfs).Nx();
					else if(lhs_nttf!=PFF(iwfs, iwfs).Nx()){
						error("PTTF for different WFS does not match: %d vs %ld\n", rhs_nttf, PTTF(iwfs, iwfs).Nx());
					}
				}else if(!parms->recon.split || parms->tomo.splitlrt==2){
					lhs_nttf=rhs_nttf;
				}
				if(parms->recon.split && parms->tomo.splitlrt){
					if(!lhs_skip0) lhs_skip0=1+iwfs;
				}
			}
			GPDATA[iwfs].neai=(const Real(*)[3])neai[iwfs]();
			//dbg("GPU: nea[%d]=%p\n", iwfs, GPDATA[iwfs].neai);
			GPDATA[iwfs].nsa=recon->saloc->p[ipowfs]->nloc;
			GPDATA[iwfs].nxp=recon->pmap->nx;
			GPDATA[iwfs].dxp=recon->pmap->dx;
			GPDATA[iwfs].dyp=recon->pmap->dy;
			GPDATA[iwfs].oxp=recon->pmap->ox;
			GPDATA[iwfs].oyp=recon->pmap->oy;
		}
		info("GPU tomography: lhs_skip0=%d, lhs_nttf=%d, rhs_nttf=%d\n", lhs_skip0, lhs_nttf, rhs_nttf);
		gpdata.init(nwfs, 1);
		DO(cudaMemcpy(gpdata(), GPDATA(), sizeof(gpu_gp_t)*nwfs, H2D));
	}

	{/**Initialize run time data*/
		//dbg("init runtime temporary data\n");
		int nxp=recon->pmap->nx;
		int nyp=recon->pmap->ny;
		NumArray<int> nxpw(nwfs), nypw(nwfs), ngw(nwfs);

		smat* wfsrot2=0;
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			if(parms->wfsr[iwfs].misregc){
				if(!wfsrot2){
					wfsrot2=snew(2, nwfs);
				}
				P(wfsrot2, 0, iwfs)=cos(parms->wfsr[iwfs].misregc);
				P(wfsrot2, 1, iwfs)=sin(parms->wfsr[iwfs].misregc);
			}
			const int ipowfs=parms->wfsr[iwfs].powfs;
			if(parms->powfs[ipowfs].skip){
				nxpw[iwfs]=0;
				nypw[iwfs]=0;
				ngw[iwfs]=0;
			} else{
				nxpw[iwfs]=nxp;
				nypw[iwfs]=nyp;
				ngw[iwfs]=recon->saloc->p[ipowfs]->nloc*2;
			}
		}
		if(wfsrot2){
			opdwfs2=curcell(nwfs, 1, nxpw(), nypw());
			cp2gpu(wfsrot, wfsrot2);
			sfree(wfsrot2);
		}
		opdwfs=curcell(nwfs, 1, nxpw(), nypw());

		grad=curcell(nwfs, 1, ngw(), (int*)NULL);
		ttf=curmat(rhs_nttf,nwfs);
	}
}

/*
  If merge the operation in to gpu_prop_grid_do, need to do atomic
  operation because previous retracing starts at offo, not 0.  */
__global__ void gpu_laplacian_do(lap_t* datai, Real* const* outall, const Real* const* inall, int nwfs, Real alpha){
	const int ips=blockIdx.z;
	const Real* in=inall[ips];
	Real* out=outall[ips];
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
		int iy1=iy+1; if(iy1>ny1) iy1-=ny;
		int iy2=iy+2; if(iy2>ny1) iy2-=ny;
		int iy1p=iy-1; if(iy1p<0) iy1p+=ny;
		int iy2p=iy-2; if(iy2p<0) iy2p+=ny;
		for(int ix=ix0; ix<nx; ix+=stepx){
			int ix1=ix+1; if(ix1>nx1) ix1-=nx;
			int ix2=ix+2; if(ix2>nx1) ix2-=nx;
			int ix1p=ix-1; if(ix1p<0) ix1p+=nx;
			int ix2p=ix-2; if(ix2p<0) ix2p+=nx;
			out[ix+iy*nx]+=alpha2*(20.f*in[ix+iy*nx]
				-8.f*(in[ix1p+iy*nx]+in[ix1+iy*nx]+in[ix+iy1p*nx]+in[ix+iy1*nx])
				+2.f*(in[ix1p+iy1p*nx]+in[ix1+iy1p*nx]+in[ix1p+iy1*nx]+in[ix1+iy1*nx])
				+(in[ix+iy2p*nx]+in[ix2p+iy*nx]+in[ix2+iy*nx]+in[ix+iy2*nx]));
		}
	}
	if(datai->zzi>-1){/*piston constaint*/
		if(threadIdx.x==0&&threadIdx.y==0&&blockIdx.x==0&&blockIdx.y==0){
			out[datai->zzi]+=in[datai->zzi]*datai->zzv*alpha;
		}
	}
}

/**
   The third grid dimension tells the wfs to handle.
   gout = GP*wfsopd    #if wfsopd is set
   ttf = - PTTF * gout #if nttf is set
*/
#define DIM_GP 128
__global__ static void gpu_gp_do(gpu_gp_t* data, Real* const* gout, Real* ttfout, Real* const* wfsopd, int nttf, int skip0){
	__shared__ Real ggf[3][DIM_GP];
	//__shared__ Real gy[DIM_GP];
	//__shared__ Real gf[DIM_GP];
	const int iwfs=blockIdx.z;
	//const int nwfs=gridDim.z;
	gpu_gp_t* datai=data+iwfs;
	const int pos=datai->pos;
	if(!pos) return;
	const int nsa=datai->nsa;
	const int step=blockDim.x*gridDim.x;
	int(*restrict saptr)[2]=datai->saptr;
	Real* restrict g=gout[iwfs];
	Real GPscale=datai->GPscale;
	if(datai->GPp&&wfsopd){
		const Real* restrict map=wfsopd[iwfs];
		const short2* restrict pxy=datai->GPp;
		int nx=datai->nxp;
		/*GP operation.*/
		if(pos==1){
			for(int isa=blockIdx.x*blockDim.x+threadIdx.x; isa<nsa; isa+=step){
				int ix=saptr[isa][0];
				int iy=saptr[isa][1];

				const short2* restrict pxy2=pxy+isa*4;
				g[isa]=GPscale*(
					+map[iy*nx+ix]*pxy2[0].x
					+map[iy*nx+ix+1]*pxy2[1].x
					+map[(iy+1)*nx+ix]*pxy2[2].x
					+map[(iy+1)*nx+ix+1]*pxy2[3].x);

				g[isa+nsa]=GPscale*(
					+map[iy*nx+ix]*pxy2[0].y
					+map[iy*nx+ix+1]*pxy2[1].y
					+map[(iy+1)*nx+ix]*pxy2[2].y
					+map[(iy+1)*nx+ix+1]*pxy2[3].y);
			}/*for isa */
		} else{
			for(int isa=blockIdx.x*blockDim.x+threadIdx.x; isa<nsa; isa+=step){
				int ix=saptr[isa][0];
				int iy=saptr[isa][1];
				const short2* restrict pxy2=pxy+isa*9;
				g[isa]=GPscale*(
					+map[iy*nx+ix]*pxy2[0].x
					+map[iy*nx+ix+1]*pxy2[1].x
					+map[iy*nx+ix+2]*pxy2[2].x
					+map[(iy+1)*nx+ix]*pxy2[3].x
					+map[(iy+1)*nx+ix+1]*pxy2[4].x
					+map[(iy+1)*nx+ix+2]*pxy2[5].x
					+map[(iy+2)*nx+ix]*pxy2[6].x
					+map[(iy+2)*nx+ix+1]*pxy2[7].x
					+map[(iy+2)*nx+ix+2]*pxy2[8].x);
				g[isa+nsa]=GPscale*(
					+map[iy*nx+ix]*pxy2[0].y
					+map[iy*nx+ix+1]*pxy2[1].y
					+map[iy*nx+ix+2]*pxy2[2].y
					+map[(iy+1)*nx+ix]*pxy2[3].y
					+map[(iy+1)*nx+ix+1]*pxy2[4].y
					+map[(iy+1)*nx+ix+2]*pxy2[5].y
					+map[(iy+2)*nx+ix]*pxy2[6].y
					+map[(iy+2)*nx+ix+1]*pxy2[7].y
					+map[(iy+2)*nx+ix+2]*pxy2[8].y);
			}/*for isa */
		}
	}
	/* Global TT, Focus projection. Make sure  each thread handle the same
	subaperture as previous gradient operation to avoid synchronization issue.*/
	if(datai->PTTF&&nttf&&(!skip0 || iwfs+1!=skip0)){//npttf has number of modes
		Real *PTTF=(nttf==1)?datai->PFF:datai->PTTF;
		for(int im=0; im<nttf; im++){
			ggf[im][threadIdx.x]=0;
		}
		for(int isa=blockIdx.x*blockDim.x+threadIdx.x; isa<nsa; isa+=step){/*ng is nsa*2. */
			for (int im=0; im<nttf; im++){
				ggf[im][threadIdx.x]+=PTTF[isa*nttf+im]*g[isa]+PTTF[(isa+nsa)*nttf+im]*g[isa+nsa];
			}
		}
		for(int nstep=(DIM_GP>>1); nstep>0; nstep>>=1){
			__syncthreads();
			if(threadIdx.x<nstep){
				for(int im=0; im<nttf; im++){
					ggf[im][threadIdx.x]+=ggf[im][threadIdx.x+nstep];
				}
			}
		}
		if(threadIdx.x==0){
			for(int im=0; im<nttf; im++){
				atomicAdd(&ttfout[iwfs*nttf+im], -ggf[im][0]);
			}
		}
	}
}
/**
   Todo: Improve by removing atomic operation.
   This functions two tasks:
   1) gin=nea*(gin + ttf)
   2) wfsopd = GP' * gin #if GP' is available.
   Be carefull about the ptt flag. It is always 1 for Right hand side, but may be zero for Left hand side.
*/
__global__ static void gpu_gpt_do(gpu_gp_t* data, Real* const* wfsopd, const Real* ttfin, Real* const* gin, int nttf){
	const int iwfs=blockIdx.z;
	//const int nwfs=gridDim.z;
	gpu_gp_t* datai=data+iwfs;
	const int pos=datai->pos;
	if(!pos) return;
	const int step=blockDim.x*gridDim.x;
	const int nsa=datai->nsa;
	int(* const saptr)[2]=datai->saptr;
	const Real(*restrict neai)[3]=datai->neai;
	const Real dxp=datai->dxp;
	const Real oxp=datai->oxp;
	const Real oyp=datai->oyp;

	Real* restrict g=gin[iwfs];
	Real* restrict map=wfsopd?wfsopd[iwfs]:0;
	Real GPscale=datai->GPscale;
	Real ttx=0, tty=0, focus=0;
	if(datai->PTTF&&nttf){
		if(nttf>1){
			ttx=ttfin[iwfs*nttf+0];
			tty=ttfin[iwfs*nttf+1];
			if(nttf>2) focus=ttfin[iwfs*nttf+2];
		} else if(nttf==1){//only focus
			focus=ttfin[iwfs*nttf];
		}
	}
	const int nx=datai->nxp;
	const short2* restrict pxy=datai->GPp;

	if(!pxy||!map){//No GP' operation
		for(int isa=blockIdx.x*blockDim.x+threadIdx.x; isa<nsa; isa+=step){
			int ix=saptr[isa][0];
			int iy=saptr[isa][1];
			Real cx=neai[isa][0];
			Real cy=neai[isa][1];
			Real cxy=neai[isa][2];
			Real gx=g[isa]+ttx+focus*(ix*dxp+oxp);
			Real gy=g[isa+nsa]+tty+focus*(iy*dxp+oyp);
			g[isa]=cx*gx+cxy*gy;
			g[isa+nsa]=cxy*gx+cy*gy;
		}
	} else if(pos==1){//Both GP' and NEA operation
		for(int isa=blockIdx.x*blockDim.x+threadIdx.x; isa<nsa; isa+=step){
			int ix=saptr[isa][0];
			int iy=saptr[isa][1];
			Real cx=neai[isa][0];
			Real cy=neai[isa][1];
			Real cxy=neai[isa][2];
			Real gx=g[isa]+ttx+focus*(ix*dxp+oxp);
			Real gy=g[isa+nsa]+tty+focus*(iy*dxp+oyp);
			Real tmp=cxy*gx;//save data
			gx=GPscale*(cx*gx+cxy*gy);
			gy=GPscale*(tmp+cy*gy);
			const short2* restrict pxy2=pxy+isa*4;
			atomicAdd(&map[iy*nx+ix], gx*pxy2[0].x+gy*pxy2[0].y);
			atomicAdd(&map[iy*nx+ix+1], gx*pxy2[1].x+gy*pxy2[1].y);
			atomicAdd(&map[(iy+1)*nx+ix], gx*pxy2[2].x+gy*pxy2[2].y);
			atomicAdd(&map[(iy+1)*nx+ix+1], gx*pxy2[3].x+gy*pxy2[3].y);
		}
	} else if(pos==2){
		for(int isa=blockIdx.x*blockDim.x+threadIdx.x; isa<nsa; isa+=step){
			int ix=saptr[isa][0];
			int iy=saptr[isa][1];
			Real cx=neai[isa][0];
			Real cy=neai[isa][1];
			Real cxy=neai[isa][2];
			Real gx=g[isa]+ttx+focus*(ix*dxp+oxp);
			Real gy=g[isa+nsa]+tty+focus*(iy*dxp+oyp);
			Real tmp=cxy*gx;
			gx=GPscale*(cx*gx+cxy*gy);
			gy=GPscale*(tmp+cy*gy);
			const short2* restrict pxy2=pxy+isa*9;
			atomicAdd(&map[iy*nx+ix], gx*pxy2[0].x+gy*pxy2[0].y);
			atomicAdd(&map[iy*nx+ix+1], gx*pxy2[1].x+gy*pxy2[1].y);
			atomicAdd(&map[iy*nx+ix+2], gx*pxy2[2].x+gy*pxy2[2].y);
			atomicAdd(&map[(iy+1)*nx+ix], gx*pxy2[3].x+gy*pxy2[3].y);
			atomicAdd(&map[(iy+1)*nx+ix+1], gx*pxy2[4].x+gy*pxy2[4].y);
			atomicAdd(&map[(iy+1)*nx+ix+2], gx*pxy2[5].x+gy*pxy2[5].y);
			atomicAdd(&map[(iy+2)*nx+ix], gx*pxy2[6].x+gy*pxy2[6].y);
			atomicAdd(&map[(iy+2)*nx+ix+1], gx*pxy2[7].x+gy*pxy2[7].y);
			atomicAdd(&map[(iy+2)*nx+ix+2], gx*pxy2[8].x+gy*pxy2[8].y);
		}
	}
}

void cutomo_grid::do_gp(curcell& _grad, const curcell& _opdwfs, int ptt2, int skip0, stream_t& stream){
	if(_opdwfs){//Sparse based is only needed when _opdwfs is set.
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			if(GPp[iwfs]||!GP[iwfs]) continue;
			_grad[iwfs].Zero(stream);//sparse based, need to zero it first.
			cuspmul(_grad[iwfs], GP[iwfs], _opdwfs[iwfs], 1, 'n', 1., stream);
		}
	}
	cuzero(ttf, stream);
	gpu_gp_do<<<dim3(24, 1, nwfs), dim3(DIM_GP, 1), 0, stream>>> //_opdwfs is initialized in the kernel.
		(gpdata(), _grad.pm(), ttf(), _opdwfs?_opdwfs.pm():NULL, ptt2, skip0);
}
void cutomo_grid::do_gpt(curcell& _opdwfs, curcell& _grad, int ptt2, stream_t& stream){
	if(_opdwfs){
		_opdwfs.Zero(stream);
	}
	//Does  GP'*NEA*(1-TTDF) if _opdwfs!=0 and GPp!=0 or NEA*(1-TTDF)
	gpu_gpt_do<<<dim3(24, 1, nwfs), dim3(DIM_GP, 1), 0, stream>>>
		(gpdata(), _opdwfs?_opdwfs.pm():0, ttf(), _grad.pm(), ptt2);

	if(_opdwfs){//Does GP' for GP with sparse
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			if(!GPp[iwfs] && GP[iwfs]){
				//info("wfs %d: opdwfs=%p, grad=%p\n", iwfs, _opdwfs[iwfs](), _grad[iwfs]());
				cuspmul(_opdwfs[iwfs], GP[iwfs], _grad[iwfs], 1, 't', 1., stream);
			}
		}
	}
}
/**
   apply HX from xin to opdwfs.
*/
void cutomo_grid::HX(const curcell& xin, Real alpha, stream_t& stream){
	if(wfsrot){
		opdwfs2.M().Zero(stream);
		hx.forward(opdwfs2.pm, xin.pm, alpha, NULL, stream);
		map_rot(opdwfs, opdwfs2, wfsrot, -1, stream);
	} else{
		opdwfs.M().Zero(stream);
		hx.forward(opdwfs.pm, xin.pm, alpha, NULL, stream);
	}
}
/**
   Apply HX' from opdwfs to xout
*/
void cutomo_grid::HXT(curcell& xout, Real alpha, stream_t& stream){
	if(wfsrot){
		map_rot(opdwfs2, opdwfs, wfsrot, 1, stream);//map_rot zeros the output array
		hx.backward(opdwfs2.pm, xout.pm, alpha, NULL, stream);
	} else{
		hx.backward(opdwfs.pm, xout.pm, alpha, NULL, stream);
	}
}
/*
  Tomography right hand side matrix. Computes xout = xout *beta + alpha * Hx' G' C * xin.
  xout is zeroed out before accumulation.
*/
void cutomo_grid::R(curcell& xout, Real beta, curcell& _grad, Real alpha, stream_t& stream){
	if(!xout){
		xout=curcell(grid->npsr, 1, grid->xnx, grid->xny);
	} else{
		curcellscale(xout, beta, stream);
	}
	do_gp(_grad, curcell(), rhs_nttf, 0, stream);
	do_gpt(opdwfs, _grad, rhs_nttf, stream);
	HXT(xout, alpha, stream);
#if SAVE_TOMO
	static int count=0; count++;
	cuwrite(_grad, stream, "tomo_R_grad_%d", count);
	cuwrite(ttf, stream, "tomo_R_gp_ttf_%d", count);
	cuwrite(opdwfs, stream, "tomo_R_opdwfs_%d", count);
	cuwrite(xout, stream, "tomo_R_xout_%d", count);
#endif	
}

void cutomo_grid::Rt(curcell& gout, Real beta, const curcell& xin, Real alpha, stream_t& stream){
	if(!gout){
		gout=curcell(nwfs, 1, grid->ngrad, (long*)NULL);
	} else{
		curcellscale(gout, beta, stream);
	}
	HX(xin, alpha, stream);
	do_gp(gout, opdwfs, rhs_nttf, 0, stream);
	curcell dummy;
	do_gpt(dummy, gout, rhs_nttf, stream);
}

/*
  Tomography left hand size matrix. Computes xout = beta*xout + alpha * Hx' G' C Gp Hx * xin.
  xout is zeroed out before accumulation.

  Be Careful about big kernels
  1) When a kernel thread reads locations written by other threads, synchronization may be required unless gauranteed in the same wrap.
  2) When a kernel thread writes to locations that may be written by others, atomic operations are required unless gauranteed in the same wrap.
*/
void cutomo_grid::L(curcell& xout, Real beta, const curcell& xin, Real alpha, stream_t& stream){
	if(!xout){
		xout=curcell(grid->npsr, 1, grid->xnx, grid->xny);
	} else{
		Scale(xout.M(), beta, stream);
	}
	ctoc_init(10);
	//xin to opdwfs
	HX(xin, 1, stream);
	ctoc("Hx");
#if SAVE_TOMO
	static int count=0; count++;
	cuwrite(xin, stream, "tomo_L_xin_%d", count);
	cuwrite(opdwfs, stream, "tomo_L_opdwfs_%d", count);
#endif
	//opdwfs to grad to ttf
	do_gp(grad, opdwfs, lhs_nttf, lhs_skip0, stream);
#if SAVE_TOMO
	cuwrite(grad, stream, "tomo_L_grad_%d", count);
	cuwrite(ttf, stream, "tomo_L_ttf_%d", count);
#endif	
	ctoc("Gp");
	//grad and ttf to opdwfs
	do_gpt(opdwfs, grad, lhs_nttf, stream);
	ctoc("Gpt");
	//opdwfs to xout
	HXT(xout, alpha, stream);
#if SAVE_TOMO
	cuwrite(opdwfs, stream, "tomo_L_opdwfs2_%d", count);
	cuwrite(xout, stream, "tomo_L_xout_%d", count);
#endif
	ctoc("HxT");
	/*This could be in parallel to hx->forward, do_gp, do_gpt*/
	gpu_laplacian_do<<<dim3(3, 3, grid->npsr), dim3(16, 16), 0, stream>>>
		(lap(), xout.pm, xin.pm, nwfs, alpha);
#if SAVE_TOMO
	cuwrite(xout, stream, "tomo_L_xout2_%d", count);
#endif
	ctoc("L2");
	ctoc_final("TomoL");
	//overhead of TomoL is 27 micro-seconds (timing without synchornization).
}


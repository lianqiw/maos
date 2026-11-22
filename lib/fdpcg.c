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
/*
   Changelog
   2010-05-27:
   Change to follow  LAOS convention.
   1) use normalized FFT. changed prop accordingly by *0.5
   2) use Laplacian to generate Cxx instead of PSD to follow forward reg.
   3) Fixed a bug in dispx,dispy. should not include cone effect since dk doesn't

   2011-10-19:
   1) No need to do fftshift after the first fft. this is done by the incorporting
   the fftshift into the permutation vector.

   2012-07-11:
   Do not need cfft2s any more. Minv is properly scaled to account for FFT scaling.
 */

#include "fdpcg.h"
#if USE_CUDA
#define PRE_PERMUT 0
#else
#define PRE_PERMUT 0 /*1: apply permutation to the inverse of the sparse matrix.
			   0: apply permutation to the vectors: faster*/
#endif
#undef TIMING
#define TIMING 0
#if !TIMING
#define TIC_tm
#define tic_tm
#define toc_tm(A)
#else
#define TIC_tm TIC
#define tic_tm tic
#define toc_tm(A) toc2(A);tic
#endif
/**
   Create aperture selection function that selects the gradients for valid
   subapertures from ground layer xloc (or ploc).

   2012-04-06: Found that by clipping xsel in FD, the strength of xsel in the
   spatial domain (SD) is reduced because of reduced strength of xsel in FD. Need
   to renormalize xsel in FD to have the same total, which corresponds to the center
   value in SD.

   clipping xsel in FD essentially removed the boundary conditiona, making
   infinite subapertures. The sum of xsel in FD has to be 1.

   @param[in] nx Grid size along x.
   @param[in] ny Grid size along y.
   @param[in] dx Grid sampling (m).
   @param[in] saloc Subaperture lower left corner coordinate.
   @param[in] saa Normalized subaperture area. Subapertures with area above 0.9 are selected.
   @return Newly allocated Fourier domain convolution operator of the selection function. The caller frees it.
*/
csp*
fdpcg_saselect(long nx, long ny, real dx, loc_t* saloc, real* saa){
	const long threas=1;
	cmat* xsel=cnew(nx, ny);
	//cfft2plan(xsel,-1);
	cmat* pxsel=xsel/*PCMAT*/;
	real dx1=1./dx;
	long offx=nx/2;
	long offy=ny/2;
	for(long isa=0; isa<saloc->nloc; isa++){
		if(saa[isa]>0.9){
			long ix=(saloc->locx[isa])*dx1+offx;/*subaperture lower left corner. */
			long iy=(saloc->locy[isa])*dx1+offy;
			P(pxsel, ix, iy)=1./(real)(nx*ny);/*cancel FFT scaling.*/
		}
	}
	cfftshift(xsel);
	cfft2(xsel, -1);
	cfftshift(xsel);
	real xselc=creal(P(pxsel, nx/2, ny/2))*threas;/*Fourier center */

	for(long ix=0; ix<nx; ix++){
		for(long iy=0; iy<ny; iy++){
			if(cabs(P(pxsel, ix, iy))<xselc){
				P(pxsel, ix, iy)=0;
			}
		}
	}
	/*scale xsel to have the same sum, to preserve its value in spatial domain. (2012-04-06)*/
	//cscale(xsel, 1./(real)csum(xsel)); //degrades performance.
	csp* sel=cspconvolvop(xsel);
	cfree(xsel);
	return sel;
}

/**
   Create a permulation acting on xhat (fft of x on xloc) so that points of the
   same spatial frequency in all layers are grouped together.

   In this implementation, it implies that all the screens have the same spatial
   sampling frequency, or they have the same tot spatial size.

   The ray tracing couples the screens but only on points with same or harmonic
   spatial frequencies. Group them together.

   shift=1: apply a fftshift with permutation vector.
   half=1: only using (half+1) of the inner most dimension (FFT is hermitian for real matrix)

   @param[in] nx Grid size along x of each layer.
   @param[in] ny Grid size along y of each layer.
   @param[in] os Over sampling ratio of each layer relative to the subaperture.
   @param[in] bs Size of each diagonal block.
   @param[in] nps Number of layers.
   @param[in] shift Whether to apply fftshift in the permutation vector.
   @param[in] half Whether to only use about half of the frequencies.
   @return Newly allocated permutation vector. Negative entries indicate a conjugate is needed. The caller frees it.
*/
lmat*
fdpcg_perm(const long* nx, const long* ny, const real* os, int bs, int nps, int shift, int half){
	long nx2[nps], ny2[nps];
	long noff[nps];
	long xloctot=0;
	for(long ips=0; ips<nps; ips++){
		nx2[ips]=nx[ips]/2;
		ny2[ips]=ny[ips]/2;
		noff[ips]=xloctot;
		xloctot+=nx[ips]*ny[ips];
	}
	if(xloctot>(1<<29)-1){
	//we reserved 2 bit for load and store flag.
		error("Current storage format overflows\n");
	}
	const long adimx=nx[0]/os[0];//subaperture dimension.
	const long adimy=ny[0]/os[0];
	const long adimx2=adimx/2;
	if(half){
		xloctot=(adimx2+1)*adimy*bs;//about half of xloctot.
	}
	lmat* perm=lnew(xloctot, 1);
	long count=0;
	/*We select the positive x frequencies first and negative x frequencies
	  later so that in GPU, when at os0 or os6, we only need to do postive
	  frequency operations */
	int nxh=half?1:2;
	long ix0, ix1;
	for(int ixh=0; ixh<nxh; ixh++){
		if(ixh==0){
			ix0=-adimx2;
			ix1=1;
		} else{
			ix0=1;
			ix1=adimx2;
		}
		for(long iy=-adimy/2; iy<adimy/2; iy++){
			for(long ix=ix0; ix<ix1; ix++){
			/*(ix,iy) is the frequency in SALOC grid. We are going to group all
			  points couple to this subaperture together. First loop over all
			  points on PLOC and find all coupled points in XLOC*/
				for(long ips=0; ips<nps; ips++){
					for(long jos=0; jos<os[ips]; jos++){
						long jy=(iy+adimy*jos);
						if(shift){
							if(jy<0) jy+=ny[ips];
						} else{
							jy=(jy+ny2[ips]);
							if(jy>=ny[ips]) jy-=ny[ips];
						}
						for(long ios=0; ios<os[ips]; ios++){
							long jx=(ix+adimx*ios);
							int ratio=1;/*no scale by default*/
							//jx may be negative already here.
							if(shift){
								if(jx<0) jx+=nx[ips];/*we embeded a fftshift here.*/
								if(half&&jx>nx2[ips]){
									ratio=-1; /*reverse the perm value to indicate a conjugate is needed*/
									jx=nx[ips]-jx;/*FFT of jx>n/2 equals to conj(n-jx)*/
								}
							} else{
								jx=(jx+nx2[ips]);
								if(jx>=nx[ips]) jx-=nx[ips];
							}
							P(perm, count)=ratio*(noff[ips]+jx+(ratio<0?(jy==0?0:ny[ips]-jy):jy)*nx[ips]);
							count++;
						}
					}
				}
			}
		}
	}
	if(count!=xloctot){
		error("count=%ld, xloctot=%ld\n", count, xloctot);
	}
	return perm;
}

/**
   Compute gradient operator in Fourier domain. Subapertures are denoted with lower left corner, so no shift is required.

   @param[out] gx Gradient operator along x. Newly allocated.
   @param[out] gy Gradient operator along y. Newly allocated.
   @param[in] nx Grid size along x.
   @param[in] ny Grid size along y.
   @param[in] dx Grid sampling (m).
   @param[in] dsa Subaperture size (m). Must be multiple of dx.
*/
void
fdpcg_g(cmat** gx, cmat** gy, long nx, long ny, real dx, real dsa){
	long os=(long)round(dsa/dx);
	if(fabs(dsa-dx*os)>1.e-10){
		error("dsa must be multiple of dx. dsa=%g, dx=%g, diff=%g\n", dsa, dx, fabs(dsa-dx*os));
	}

	real* wt=mycalloc(os+1, real);
	real* st=mycalloc(os+1, real);
	/*Trapzoidal weights for averaging. */
	wt[os]=wt[0]=0.5/(real)os/dsa;
	for(long ios=1; ios<os; ios++){
		wt[ios]=1./(real)os/dsa;
	}
	for(long ios=0; ios<os+1; ios++){
		st[ios]=ios*dx;
	}
	long ny2=ny/2;
	long nx2=nx/2;
	comp cf=COMPLEX(0, 2*M_PI);
	real dfy=1/(ny*dx);
	real dfx=1/(nx*dx);
	*gx=cnew(nx*ny, 1);
	*gy=cnew(nx*ny, 1);
	comp* pgx=P(*gx);
	comp* pgy=P(*gy);
	for(long iy=0; iy<ny; iy++){
		real fy=(real)(iy-ny2)*dfy;
		for(long ix=0; ix<nx; ix++){
			real fx=(real)(ix-nx2)*dfx;
			comp tx=0;
			comp ty=0;
			for(int ios=0; ios<os+1; ios++){
				tx+=wt[ios]*(cexp(cf*(fx*dsa+fy*st[ios]))-cexp(cf*(fy*st[ios])));
				ty+=wt[ios]*(cexp(cf*(fy*dsa+fx*st[ios]))-cexp(cf*(fx*st[ios])));
			}
			pgx[ix+iy*nx]=tx;
			pgy[ix+iy*nx]=ty;
		}
	}
	free(wt);
	free(st);
	/*Remove global tt is same as remove zero frequency value in fourier space.*/
	/*if(ttr){
	warning("fdpcg: Remove global t/t (%g %g) (%g, %g)",
		creal(pgx[nx/2+ny/2*nx]), cimag(pgx[nx/2+ny/2*nx]),
		creal(pgy[nx/2+ny/2*nx]), cimag(pgy[nx/2+ny/2*nx]));
	pgx[nx/2+ny/2*nx]=0;
	pgy[nx/2+ny/2*nx]=0;
	}*/
}

/**
   Propagate operator for nlayer screens to ground of size
   nxp*nxp, sampling dx, with displacement of dispx, dispy.

   @param[in] nps Number of layers.
   @param[in] nxp Ground grid size along x.
   @param[in] nyp Ground grid size along y.
   @param[in] nx Grid size along x of each layer.
   @param[in] ny Grid size along y of each layer.
   @param[in] dx Ground grid sampling (m).
   @param[in] dispx Displacement along x of each layer (m).
   @param[in] dispy Displacement along y of each layer (m).
   @return Newly allocated transpose of the Fourier domain propagation operator. The caller frees it.
*/
csp*
fdpcg_prop_trans(long nps, long nxp, long nyp, long* nx, long* ny, real dx, real* dispx, real* dispy){
	long nx2[nps], nx3[nps];
	long ny2[nps], ny3[nps];
	long noff[nps];
	long nxtot=0;
	real dkx=1./(nxp*dx);
	real dky=1./(nyp*dx);
	for(long ips=0; ips<nps; ips++){
		nx2[ips]=nx[ips]/2;
		nx3[ips]=nx2[ips]+nx[ips];
		ny2[ips]=ny[ips]/2;
		ny3[ips]=ny2[ips]+ny[ips];
		noff[ips]=nxtot;
		nxtot+=nx[ips]*ny[ips];
	}
	long nxp2=nxp/2;
	long nyp2=nyp/2;
	/*We build the transposed matrix.*/
	csp* propt=cspnew(nxtot, nxp*nyp, nxp*nyp*nps);
	spint* pp=propt->pp;
	spint* pi=propt->pi;
	comp* px=propt->px;
	long count=0;
	for(long iy=0; iy<nyp; iy++){
		real fyg=(iy-nyp2)*dky;
		for(long ix=0; ix<nxp; ix++){
			real fxg=(ix-nxp2)*dkx;/*spatial frequency in pupil. */
			long icol=ix+iy*nxp;
			pp[icol]=count;
			for(long ips=0; ips<nps; ips++){
				long jx=((ix-nxp2)+nx3[ips])%nx[ips];/*map to layer ips with wrapping because FFT wraps around. */
				long jy=((iy-nyp2)+ny3[ips])%ny[ips];/*wrapping improves performance (periodic boundary condition)*/
				real fx=(jx-nx2[ips])*dkx;/*spatial frequency in plane ips. */
				real fy=(jy-ny2[ips])*dky;
				pi[count]=jx+jy*nx[ips]+noff[ips];
				comp shiftc=EXPI(-TWOPI*(fx*dispx[ips]+fy*dispy[ips]));//conjugate of shift
				switch(nxp/nx[ips]){
				case 0:
				case 1:
					px[count]=shiftc;
					break;
				case 2:
				{
				/*
				  2010-05-28:
				  We are using hermitian FFT in that FFT'=IFFT so that
				  we can apply G'G. We need to multiply 0.5 to account
				  for this scaling effect.
				*/
				/*
				  comp shiftx=cexp(cf*(fxg*dx));
				  comp shifty=cexp(cf*(fyg*dx));
				  px[count]=conj(shift*(1+0.5*shiftx+0.5*conj(shiftx))
				  *(1+0.5*shifty+0.5*conj(shifty)));
				  //the following is equivalent.
				  */
					px[count]=shiftc*((1.+cos(TWOPI*(fxg*dx)))*(1.+cos(TWOPI*(fyg*dx)))*0.5);
				}
				break;
				default:
					error("Invalid: %ld\n", nxp/nx[ips]);
				}
				count++;
			}
		}
	}
	pp[nxp*nyp]=count;
	/*we put conj above because csptrans applies a conjugation. */
	//csp* propf=csptrans(propt);
	//cspfree(propt);
	return propt;
}

/**
   Assemble the Fourier domain preconditioner matrix: the inverse PSD of each
   layer plus, for each WFS, the propagation of the layers to the ground,
   followed by the gradient operator and subaperture selection, weighted by the
   inverse of the NEA squared:
   Mhat=diag(invpsd)+sum_wfs prop'*(gx'*sel'*sel*gx+gy'*sel'*sel*gy)*prop/nea^2.

   @param[in] nps Number of layers.
   @param[in] nx Grid size along x of each layer.
   @param[in] ny Grid size along y of each layer.
   @param[in] nxp Ground (ploc) grid size along x.
   @param[in] nyp Ground (ploc) grid size along y.
   @param[in] dxp Ground grid sampling (m).
   @param[in] saloc Subaperture lower left corner coordinate.
   @param[in] saa Normalized subaperture area. Subapertures with area above 0.9 are selected.
   @param[in] invpsd Inverse PSD of all layers concatenated (nxtot x 1), with the zero frequency at the center.
   @param[in] nwfs Number of WFS.
   @param[in] dispx Displacement along x of each layer for each WFS (m), nps x nwfs.
   @param[in] dispy Displacement along y of each layer for each WFS (m), nps x nwfs.
   @param[in] neai NEA of each WFS (rad).
   @param[in] iwfs WFS number of each WFS, used in the file names when save is set. May be NULL to use the index.
   @param[in] save Save the intermediate matrices to files fdpcg_*.
   @return Newly allocated matrix. The caller frees it.
*/
csp* fdpcg_mhat(long nps, long* nx, long* ny, long nxp, long nyp, real dxp, loc_t* saloc, real* saa, const cmat* invpsd, long nwfs, const dmat* dispx, const dmat* dispy, const real* neai, const long* iwfs, int save){
	long nxtot=0;
	for(long ips=0; ips<nps; ips++){
		nxtot+=nx[ips]*ny[ips];
	}
	/*Subaperture selection operator */
	csp* sel=fdpcg_saselect(nxp, nyp, dxp, saloc, saa);
	/*Gradient operator. */
	cmat* gx=0, * gy=0;
	fdpcg_g(&gx, &gy, nxp, nyp, dxp, saloc->dx);/*tested ok. */
	/*make it sparse diagonal operator */
	csp* Mhat=cspnewdiag(nxtot, P(invpsd), 1);

	csp* Mmid=NULL;
	/*Compute gx'*sel'*sel*gx+gy'*sel'*sel*gy as Mmid */
	for(int i=0; i<2; i++){
		cmat* g;
		if(i==0){
			g=gx;
		} else{
			g=gy;
		}
		csp* tmp=cspdup(sel);
		cspmuldiag(tmp, P(g), 1);
		csp* tmp2=cspmulsp(tmp, tmp, "cn");
		cspadd(&Mmid, 1, tmp2, 1);
		cspfree(tmp);
		cspfree(tmp2);
	}
	/*Ray tracing operator for each WFS */
	if(save){
		writebin(sel, "fdpcg_sel");
		writebin(gx, "fdpcg_gx");
		writebin(gy, "fdpcg_gy");
		writebin(Mhat, "fdpcg_invpsd");
		writebin(Mmid, "fdpcg_Mmid");
	}
	cfree(gx);
	cfree(gy);
	cspfree(sel);
	/* Mhat = Mhat + propx' * Mmid * propx */
OMP_FOR(nwfs)
	for(int jwfs=0; jwfs<nwfs; jwfs++){
		csp* propxt=fdpcg_prop_trans(nps, nxp, nyp, nx, ny, dxp, (real*)PCOL(dispx, jwfs), (real*)PCOL(dispy, jwfs));
		if(save){
			writebin(propxt, "fdpcg_prop_trans_wfs%ld", iwfs?iwfs[jwfs]:(long)jwfs);
		}
		/*need to test this in spatial domain. */
		//cspscale(propxt, 1./neai); //absorbed below
		/*prop is not real for off axis wfs. */
		/*Compute propx'*Mmid*propx and add to Mhat; */
		csp* tmp=cspmulsp(Mmid, propxt, "nc");
		csp* tmp2=cspmulsp(propxt, tmp, "nn");
		cspfree(tmp);
		cspfree(propxt);
OMP_CRITICAL
		cspadd(&Mhat, 1, tmp2, 1./(neai[jwfs]*neai[jwfs]));
		cspfree(tmp2);
	}
	cspfree(Mmid);
	//cspdroptol(Mhat, EPS); 
	//cspsym(&Mhat);//uncommet if not good (2024-07-23)
	return Mhat;
}
/**
   Create the preconditioner from the matrix assembled by fdpcg_mhat(): the
   matrix is permuted to block diagonal form and each block is inverted with
   truncated SVD.

   @param[in] Mhat The matrix from fdpcg_mhat(). It is freed.
   @param[in] nps Number of layers.
   @param[in] nx Grid size along x of each layer.
   @param[in] ny Grid size along y of each layer.
   @param[in] os Over sampling ratio of each layer relative to the subaperture.
   @param[in] half Also create the permutation for real to complex FFT (permhf, used by GPU).
   @param[in] xloc Grid of each layer. Referenced, not copied.
   @param[in] square Whether xloc is square (the full grid).
   @param[in] save Save the intermediate matrices to files fdpcg_*.
   @return Newly allocated preconditioner data. Free with fdpcg_free().
*/
fdpcg_t* fdpcg_new(csp* Mhat, long nps, const long* nx, const long* ny, const real* os, int half, loccell* xloc, int square, int save){
	long nxtot=0;
	int os0=os[0];
	int needscale=0;
	for(long ips=0; ips<nps; ips++){
		if(os[ips]!=os0){
			needscale=1;
		}
		nxtot+=nx[ips]*ny[ips];
	}
	if(!needscale){  /*scale Mhat to avoid scaling FFT. */
		cmat* sc=cnew(nxtot, 1);
		comp* psc=P(sc);
		for(int ips=0; ips<nps; ips++){
			real scale=(real)(nx[ips]*ny[ips]);
			for(long ix=0; ix<nx[ips]*ny[ips]; ix++){
				*(psc++)=scale;
			}
		}
		cspmuldiag(Mhat, P(sc), 1);
		cfree(sc);
	}

	if(save){
		writebin(Mhat, "fdpcg_Mhat");
	}
	fdpcg_t* fdpcg=mycalloc(1, fdpcg_t);
	/*Now invert each block. */
	/*bs: blocksize. */
	int bs=0;
	for(long ips=0; ips<nps; ips++){
		bs+=os[ips]*os[ips];
	}
	fdpcg->bs=bs;
	long nb=(nx[0]/os[0])*(ny[0]/os[0]);
	dbg("fdpcg: Block size is %d, there are %ld blocks\n", bs, nb);
	/*Permutation vector */
	fdpcg->scale=needscale;
	lmat* perm=fdpcg_perm(nx, ny, os, bs, (int)nps, 0, 0);
	csp* Mhatp=cspperm(Mhat, 0, P(perm), P(perm));/*forward permutation. */
	cspfree(Mhat);
	lfree(perm);

	perm=fdpcg_perm(nx, ny, os, bs, (int)nps, 1, 0); /*contains fft shift information*/
	if(save){
		writebin(perm, "fdpcg_perm");
	}
#if PRE_PERMUT == 1//Absorbes the permutation operation in the sparse matrix.
	csp* Minvp=cspinvbdiag(Mhatp, bs);
	csp* Minv=cspperm(Minvp, 1, P(perm), P(perm));/*revert permutation */
	lfree(perm); perm=NULL;
	cspfree(Minvp);
	fdpcg->Minv=Minv;
	if(save){
		writebin(Minv, "fdpcg_Minv");
	}
#else//use block diagonal matrix
	fdpcg->perm=perm; perm=NULL;
	if(half){
		fdpcg->permhf=fdpcg_perm(nx, ny, os, bs, (int)nps, 1, 1);
		if(save){
			writebin(fdpcg->permhf, "fdpcg_permhf");
		}
	}

	fdpcg->Mbinv=cspblockextract(Mhatp, bs);
	if(save){
		writebin(fdpcg->Mbinv, "fdpcg_Mhatb");
	}
OMP_FOR(NTHREAD)
	for(long ib=0; ib<NX(fdpcg->Mbinv); ib++){
	/*2012-04-07: was using inv_inplace that calls gesv that does not truncate svd. In
	  one of the cells the conditional is more than 1e8. This creates
	  problem in GPU code causing a lot of pistion to accumulate and
	  diverges the pcg. Use svd to truncate smaller eigen values.
	*/
		csvd_pow2(P(fdpcg->Mbinv,ib), -1, 1e-7,0);
	}
	if(save){
		writebin(fdpcg->Mbinv, "fdpcg_Minvb");
	}
#endif
	cspfree(Mhatp);
	fdpcg->xloc=xloc;
	fdpcg->square=square;
	fdpcg->nxtot=nxtot;
	fdpcg->nx=lnew(nps, 1);
	fdpcg->ny=lnew(nps, 1);
	for(long ips=0; ips<nps; ips++){
		P(fdpcg->nx, ips)=nx[ips];
		P(fdpcg->ny, ips)=ny[ips];
	}
	return fdpcg;
}
typedef struct{
	cmat* xhat;
	cmat* xhat2;
	ccell* xhati;
	ccell* xhat2i;
	const fdpcg_t* fdpcg;
	const dcell* xin;
	dcell* xout;
	unsigned int ips;
	unsigned int nps;
	unsigned int ib;
	unsigned int nb;
}fdpcg_info_t;
/**
   Copy x vector and do FFT on each layer

   @param[in,out] pdata Shared thread data. data->xhati is updated. Layers are distributed among threads with data->ips.
   @return NULL.
*/
static void* fdpcg_fft(void *pdata){
	fdpcg_info_t *data=(fdpcg_info_t *)pdata;
	const fdpcg_t* fdpcg=data->fdpcg;
	ccell* xhati=data->xhati;
	const dcell* xin=data->xin;
	unsigned int ips;
	while((ips=atomic_fetch_add(&data->ips, 1))<data->nps){
		if(fdpcg->square){
			for(long i=0; i<P(xhati,ips)->nx*P(xhati,ips)->ny; i++){
				P(P(xhati,ips),i)=P(P(xin,ips),i);
			}
		} else{
			/*
			  czero takes 0.000074
			  cembed_locstat takes 0.000037
			*/
			czero(P(xhati,ips));
			cembed_locstat(&P(xhati,ips), 0, P(fdpcg->xloc,ips), P(P(xin,ips)), 1, 0);
		}
		if(fdpcg->scale){
			cfft2s(P(xhati,ips), -1);
		} else{
			cfft2(P(xhati,ips), -1);
		}
	}
	return NULL;
}

/**
   Multiply each block in parallel

   @param[in,out] pdata Shared thread data. data->xhat is accumulated. Blocks are distributed among threads with
   data->ib.
   @return NULL.
*/
static void* fdpcg_mulblock(void *pdata){
	fdpcg_info_t *data=(fdpcg_info_t *)pdata;
	const fdpcg_t* fdpcg=data->fdpcg;
	long bs=P(fdpcg->Mbinv,0)->nx;//size of block
	unsigned int ib;
	while((ib=atomic_fetch_add(&data->ib, 1))<data->nb){
		cmulvec(&P(data->xhat,ib*bs), P(fdpcg->Mbinv,ib), &P(data->xhat2,ib*bs), 1);
	}
	return NULL;
}

/**
   Inverse FFT for each block. Put result in xout, replace content, do not accumulate.

   @param[in,out] pdata Shared thread data. data->xout is overwritten. Layers are distributed among threads with
   data->ips.
   @return NULL.
*/
static void* fdpcg_ifft(void *pdata){
	fdpcg_info_t *data=(fdpcg_info_t *)pdata;
	const fdpcg_t* fdpcg=data->fdpcg;
	ccell* xhat2i=data->xhat2i;
	dcell* xout=data->xout;
	unsigned int ips;
	while((ips=atomic_fetch_add(&data->ips, 1))<data->nps){
		if(fdpcg->scale){
			cfft2s(P(xhat2i,ips), 1);
		} else{
			cfft2(P(xhat2i,ips), 1);
		}
		if(fdpcg->square){
			for(long i=0; i<P(xhat2i,ips)->nx*P(xhat2i,ips)->ny; i++){
				P(P(xout,ips),i)=creal(P(P(xhat2i,ips),i));
			}
		} else{
			dzero(P(xout,ips));
			cembed_locstat(&P(xhat2i,ips), 1, P(fdpcg->xloc,ips), P(P(xout,ips)), 0, 1);
		}
	}
	return NULL;
}

/**
   Apply the Fourier domain preconditioner. Preconditioner solves xout=M^-1*xin
   where M is the preconditioner matrix. Here M^-1 is applied in Fourier Domain
   via xout=IFFT(M*FFT(xin));

   @param[out] xout Output. Allocated if NULL. The content is replaced.
   @param[in] fdpcg The preconditioner from fdpcg_new().
   @param[in] xin Input on the grid of each layer.
   @param[in] nthread Number of threads for multiplying the blocks.
*/
void fdpcg_apply(dcell** xout, const fdpcg_t* fdpcg, const dcell* xin, int nthread){
	TIC_tm; tic_tm;
	if(NY(xin)!=1){
		error("Invalid\n");
	}
	const unsigned int nps=NX(fdpcg->nx);
	long bs=fdpcg->bs;
	unsigned int nb=(unsigned int)(fdpcg->nxtot/bs);//nxtot<2^29 is checked in fdpcg_perm
	ccell* xhati=ccellnew3(nps, 1, P(fdpcg->nx), P(fdpcg->ny));
	ccell* xhat2i=ccellnew3(nps, 1, P(fdpcg->nx), P(fdpcg->ny));
	cmat* xhat=cref(xhati->m);
	cmat* xhat2=cref(xhat2i->m);
	if(!*xout){
		*xout=dcellnew2(xin);
	}
	fdpcg_info_t data={xhat, xhat2, xhati, xhat2i, fdpcg, xin, *xout, 0, nps, 0, nb};
	/*apply forward FFT */
#define DBG_FD 0
#if DBG_FD
	writebin(xin, "fdc_xin");
#endif
	data.ips=0;
	CALL(fdpcg_fft, &data, 1, 1);//fast
#if DBG_FD
	writebin(xhati, "fdc_fft");
#endif
#if PRE_PERMUT
	czero(xhat2);
	cspmv(xhat2, fdpcg->Minv, xhat, 'n', 1);//double check. 
#else/*permute vectors and apply block diagonal matrix */
	/*permute xhat and put into xhat2 */
	cvecperm(xhat2, xhat, P(fdpcg->perm));
#if DBG_FD
	writebin(xhat2i, "fdc_perm");
#endif
	czero(xhat);
	data.ib=0;
	CALL(fdpcg_mulblock, &data, nthread, 1);//slower than FFT. **Wrong result in taurus if use nthread=4**
#if DBG_FD
	writebin(xhati, "fdc_mul");
#endif
	cvecpermi(xhat2, xhat, P(fdpcg->perm));
#if DBG_FD
	writebin(xhat2i, "fdc_permi");
#endif
#endif
	/*Apply inverse FFT */
	data.ips=0;
	CALL(fdpcg_ifft, &data, 1, 1);
#if DBG_FD
	writebin(*xout, "fdc_xout");
	exit(1);
#endif
	ccellfree(xhati);
	ccellfree(xhat2i);
	cfree(xhat);
	cfree(xhat2);
	toc_tm("fdpcg")
}

/**
   Free fdpcg related data structs.

   @param[in] fdpcg The struct to free. May be NULL.
*/
void fdpcg_free(fdpcg_t* fdpcg){
	if(!fdpcg) return;
	cspfree(fdpcg->Minv);
	lfree(fdpcg->perm);
	lfree(fdpcg->permhf);
	lfree(fdpcg->nx);
	lfree(fdpcg->ny);
	ccellfree(fdpcg->Mbinv);
	free(fdpcg);
}

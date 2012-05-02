/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

/**
   \file fdpcg.c
   Fourier Domain Preconditioner for Tomography step.
   cfft2s(A,-1) and cfft2s(A,1) are hermitian ffts because they have the same scale.
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

 */

#include "maos.h"
#include "fdpcg.h"
#if USE_CUDA
#define PRE_PERMUT 0
#else
#define PRE_PERMUT 0 /*1: apply permutation to the inverse of the sparse matrix.
		       0: apply permutation to the vectors: faster*/
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
*/
csp* fdpcg_saselect(long nx, long ny, double dx,loc_t *saloc, double *saa){
    const long threas=1;
    cmat *xsel=cnew(nx,ny);
    cfft2plan(xsel,-1);
    PCMAT(xsel,pxsel);
    double dx1=1./dx;
    long offx=nx/2;
    long offy=ny/2;
    for(long isa=0; isa<saloc->nloc; isa++){
	if(saa[isa]>0.9){
	    long ix=(saloc->locx[isa])*dx1+offx;/*subaperture lower left corner. */
	    long iy=(saloc->locy[isa])*dx1+offy;
	    pxsel[iy][ix]=1./(double)(nx*ny);/*cancel FFT scaling.*/
	}
    }
    cfftshift(xsel);
    cfft2(xsel,-1);
    cfftshift(xsel);
    double xselc=creal(pxsel[ny/2][nx/2])*threas;/*Fourier center */
   
    for(long ix=0; ix<nx; ix++){
	for(long iy=0; iy<ny; iy++){
	    if(cabs(pxsel[iy][ix])<xselc){
		pxsel[iy][ix]=0;
	    }
	}
    }
    /*scale xsel to have the same sum, to preserve its value in spatial domain. (2012-04-06)*/
    cscale(xsel, 1./(double)csum(xsel));
    csp *sel=cspconvolvop(xsel);
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
*/
long *fdpcg_perm(const long *nx, const long *ny, const int *os, int nps, int shift, int half){
    long nx2[nps],ny2[nps];
    long noff[nps];
    long xloctot=0;
    for(long ips=0; ips<nps; ips++){
	nx2[ips]=nx[ips]/2;
	ny2[ips]=ny[ips]/2;
	noff[ips]=xloctot;
	if(half){
	    xloctot+=(nx[ips]/2+1)*ny[ips];
	}else{
	    xloctot+=nx[ips]*ny[ips];
	}
    }
    long *perm=calloc(xloctot, sizeof(long));
    // long use_os=pos;
    long adim=nx[0]/os[0];//subaperture dimension.
    //long osx=nx[0]/2;
    long count=0;

    for(long iy=-adim/2; iy<adim/2; iy++){
	for(long ix=-adim/2; ix<adim/2; ix++){
	    /*(ix,iy) is the frequency in SALOC grid. We are going to group all
	      points couple to this subaperture together. First loop over all
	      points on PLOC and find all coupled points in XLOC*/
	    for(long ips=0; ips<nps; ips++){
		for(long jos=0; jos<os[ips]; jos++){
		    long jy=(iy+adim*jos); 
		    if(jy>=ny2[ips]) jy-=ny[ips];
		    for(long ios=0; ios<os[ips]; ios++){
			long jx=(ix+adim*ios); 
			if(jx>=nx2[ips]) jx-=nx[ips];
			if(shift){ /*we embeded a fftshift here.*/
			    perm[count]=noff[ips]+(jx<0?jx+nx[ips]:jx)+(jy<0?jy+ny[ips]:jy)*nx[ips];
			}else{
			    perm[count]=noff[ips]+(jx+nx2[ips])+(jy+ny2[ips])*nx[ips];
			}
			count++;
		    }
		}
	    }
	}
    }
    return perm;
}

/**
   Compute gradient operator in Fourier domain. Subapertures are denoted with lower left corner, so no shift is required.
*/
void fdpcg_g(cmat **gx, cmat **gy, long nx, long ny, double dx, double dsa, int ttr){
    long os=(long)(dsa/dx);
    if(fabs(dsa-dx*os)>1.e-10){
	error("dsa must be multiple of dx");
    }
 
    double *wt=alloca(sizeof(double)*(os+1));
    double *st=alloca(sizeof(double)*(os+1));
    /*Trapzoidal weights for averaging. */
    wt[os]=wt[0]=0.5/(double)os/dsa;
    for(long ios=1; ios<os; ios++){
	wt[ios]=1./(double)os/dsa;
    }
    for(long ios=0; ios<os+1; ios++){
	st[ios]=ios*dx;
    }
    long ny2=ny/2;
    long nx2=nx/2;
    dcomplex cf=2*M_PI*I;
    double dfy=1/(ny*dx);
    double dfx=1/(nx*dx);
    *gx=cnew(nx*ny,1);
    *gy=cnew(nx*ny,1);
    dcomplex *pgx=(*gx)->p;
    dcomplex *pgy=(*gy)->p;
    for(long iy=0; iy<ny; iy++){
	double fy=(double)(iy-ny2)*dfy;
	for(long ix=0; ix<nx; ix++){
	    double fx=(double)(ix-nx2)*dfx;
	    dcomplex tx=0;
	    dcomplex ty=0;
	    for(int ios=0; ios<os+1; ios++){
		tx+=wt[ios]*(cexp(cf*(fx*dsa+fy*st[ios]))-cexp(cf*(fy*st[ios])));
		ty+=wt[ios]*(cexp(cf*(fy*dsa+fx*st[ios]))-cexp(cf*(fx*st[ios])));
	    }
	    pgx[ix+iy*nx]=tx;
	    pgy[ix+iy*nx]=ty;
	}
    }
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
*/
csp *fdpcg_prop(long nps, long pos, long nxp, long nyp, long *nx, long *ny, double dx, double *dispx, double *dispy){
    long nx2[nps],nx3[nps];
    long ny2[nps],ny3[nps];
    long noff[nps];
    long nxtot=0;
    double dkx=1./(nxp*dx);
    double dky=1./(nyp*dx);
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
    csp *propt=cspnew(nxtot,nxp*nyp,nxp*nyp*nps);
    spint *pp=propt->p;
    spint *pi=propt->i;
    dcomplex *px=propt->x;
    long count=0;
    dcomplex cf=2*M_PI*I;
    double cfr=2*M_PI;
    for(long iy=0; iy<nyp; iy++){
	double fyg=(iy-nyp2)*dky;
	for(long ix=0; ix<nxp; ix++){
	    double fxg=(ix-nxp2)*dkx;/*spatial frequency in pupil. */
	    long icol=ix+iy*nxp;
	    pp[icol]=count;
	    for(long ips=0; ips<nps; ips++){
		long jx=((ix-nxp2)+nx3[ips])%nx[ips];/*map to layer ips. */
		long jy=((iy-nyp2)+ny3[ips])%ny[ips];/*map to layer ips. */
		double fx=(jx-nx2[ips])*dkx;/*spatial frequency in plane ips. */
		double fy=(jy-ny2[ips])*dky;
		pi[count]=jx+jy*nx[ips]+noff[ips];
		dcomplex shift=cexp(cf*(fx*dispx[ips]+fy*dispy[ips]));
		switch(nxp/nx[ips]){
		case 0:
		case 1:
		    px[count]=conj(shift);
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
			  dcomplex shiftx=cexp(cf*(fxg*dx));
			  dcomplex shifty=cexp(cf*(fyg*dx));
			  px[count]=conj(shift*(1+0.5*shiftx+0.5*conj(shiftx))
			  *(1+0.5*shifty+0.5*conj(shifty)));
			  //the following is equivalent. 
			  */
			px[count]=conj(shift*(1+cos(cfr*(fxg*dx)))*(1+cos(cfr*(fyg*dx))))*0.5;
		    }
		    break;
		default:
		    error("Invalid: %ld\n", nxp/nx[ips]);
		}
		count++;
	    }
	}
    }
    pp[nxp*nxp]=count;
    /*we put conj above because csptrans applies a conjugation. */
    csp *propf=csptrans(propt);
    cspfree(propt);
    return propf;
}

/**
  Prepare data for Tomography Fourier Domain Preconditioner
*/
FDPCG_T *fdpcg_prepare(const PARMS_T *parms, const RECON_T *recon, const POWFS_T *powfs){
    int hipowfs=-1;
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(!parms->powfs[ipowfs].lo){
	    if(hipowfs==-1){
		hipowfs=ipowfs;
	    }else{
		error("Multiple High order powfs found\n");
	    }
	}
    }
    loc_t *saloc=powfs[hipowfs].saloc;
    const long nps=recon->npsr;
    long pos=parms->tomo.pos;
    const int* os=parms->atmr.os;
    if(pos!=os[0]){
	warning("pupil does not equal to ground layer over sampling. Debug required.\n");
    }
    long* nx=recon->xnx;
    long* ny=recon->xny;
    const long nxp=nx[0]/os[0]*parms->tomo.pos;
    const long nyp=ny[0]/os[0]*parms->tomo.pos;
    const double dxp=recon->ploc->dx;
    const double *ht=parms->atmr.ht;
    long nxtot=0;
    for(long ips=0; ips<nps; ips++){
	nxtot+=nx[ips]*ny[ips];
	if(os[ips]>pos){
	    warning("Layer %ld os is greater than ground\n", ips);
	}
    }
    /*Subaperture selection operator */
    csp *sel=fdpcg_saselect(nxp,nyp,dxp, saloc, powfs[hipowfs].saa->p);
    /*Gradient operator. */
    cmat *gx, *gy;
    int ttr=parms->recon.split?1:0;
    fdpcg_g(&gx,&gy,nxp,nyp,dxp,saloc->dx,ttr);/*tested ok. */
    /*Concatenate invpsd; */
    dcomplex *invpsd=calloc(nxtot, sizeof(dcomplex));
    long offset=0;
    switch(parms->tomo.cxx){
    case 0:/*forward matrix uses biharmonic approx. We use here also. */
	for(long ips=0; ips<nps; ips++){
	    cmat *psd=cnew(nx[ips],ny[ips]);
	    cfft2plan(psd,-1);
	    dsp *L2;
	    if(parms->tomo.square){
		L2=spref(recon->L2->p[ips+nps*ips]);
	    }else{/*L2 is for non square xloc. need to build L2 for square xloc. */
		L2=mklaplacian_map(nx[ips], ny[ips],
				   recon->xloc[ips]->dx, recon->r0,
				   recon->wt->p[ips]);
	    }
	    dsp *tmp=sptmulsp(L2, L2);
	    spfree(L2);
	    for(long irow=tmp->p[0]; irow<tmp->p[1]; irow++){/*first column of tmp to psf. */
		psd->p[tmp->i[irow]]=tmp->x[irow];
	    }
	    spfree(tmp);
	    cfft2(psd,-1);
	    cfftshift(psd);
	    /*look for a way to obtain this automatically. */
	    const double eps=2.220446049250313e-16;
	    double max;
	    maxmincmp(psd->p,psd->nx*psd->ny,&max,NULL,NULL);
	    max=max*sqrt(eps);
	    for(long i=0; i<nx[ips]*ny[ips]; i++){
		invpsd[offset+i]=creal(psd->p[i])+max;
	    }
	    offset+=nx[ips]*ny[ips];
	    cfree(psd);
	}
	break;
    case 1:
    case 2:
	/*forward matrix uses inverse PSD or fractal. we use PSF here. */
	for(long ips=0; ips<nps; ips++){
	    dmat *tmp=ddup(recon->invpsd->invpsd->p[ips]);
	    dfftshift(tmp);
	    /*cancel the scaling applied in invpsd routine. */
	    dscale(tmp,(double)(nx[ips]*ny[ips]));
	    for(long i=0; i<nx[ips]*ny[ips]; i++){
		invpsd[offset+i]=tmp->p[i];
	    }
	    offset+=nx[ips]*ny[ips];
	    dfree(tmp);
	}
	break;
    }

    /*make it sparse diagonal operator */
    csp *Mhat=cspnewdiag(nxtot,invpsd,TOMOSCALE);/*parms->gpu.tomo?(TOMOSCALE):(1.));*/
    free(invpsd);

    csp *Mmid=NULL;
    /*Compute gx'*sel'*sel*gx+gy'*sel'*sel*gy as Mmid */
    for(int i=0; i<2; i++){
	cmat *g;
	if(i==0){
	    g=gx;
	}else{
	    g=gy;
	}
	csp *tmp=cspdup(sel);
	cspmuldiag(tmp,g->p,1);
	csp *tmp2=csptmulsp(tmp,tmp);
	cspadd(&Mmid, tmp2);
	cspfree(tmp);
	cspfree(tmp2);
    }
    /*Ray tracing operator for each WFS */
    if(parms->save.setup){
	cspwrite(sel,"%s/fdpcg_sel",dirsetup);
	cwrite(gx,"%s/fdpcg_gx",dirsetup);
	cwrite(gy,"%s/fdpcg_gy",dirsetup);
	cspwrite(Mhat,"%s/fdpcg_invpsd",dirsetup);
	cspwrite(Mmid,"%s/fdpcg_Mmid",dirsetup);
    }
    cfree(gx);
    cfree(gy);
    cspfree(sel);
    double dispx[nps];
    double dispy[nps];
    for(int jwfs=0; jwfs<parms->powfs[hipowfs].nwfs; jwfs++){
	int iwfs=parms->powfs[hipowfs].wfs[jwfs];
	double neai=recon->neam->p[iwfs];
	info2("fdpcg: mean sanea used for wfs %d is %g mas\n",iwfs, 206265000*neai);
	for(long ips=0; ips<nps; ips++){
	    /*
	      2010-05-28: The cone effect cancels with the cone
	      coordinate. so we are like doing parallel beam
	      propagation. Removed scaling by 1/(1-ht[ips]/hs);
	    */
	    dispx[ips]=ht[ips]*parms->wfs[iwfs].thetax;
	    dispy[ips]=ht[ips]*parms->wfs[iwfs].thetay;
	}
	csp *propx=fdpcg_prop(nps,pos,nxp,nyp,nx,ny,dxp,dispx,dispy);
	if(parms->save.setup){
	    cspwrite(propx,"%s/fdpcg_prop_wfs%d",dirsetup,iwfs);
	}
	/*need to test this in spatial domain. */
	cspscale(propx,sqrt(TOMOSCALE)/neai);/*prop is not real for off axis wfs. */
	/*Compute propx'*Mmid*propx and add to Mhat; */
	csp *tmp=cspmulsp(Mmid,propx);
	csp *tmp2=csptmulsp(propx,tmp);
	cspfree(tmp);
	cspfree(propx);
	cspadd(&Mhat,tmp2);
	cspfree(tmp2);
    }
  
    cspfree(Mmid);
    cspdroptol(Mhat,EPS);
    cspsym(Mhat);
    if(parms->save.setup){
	cspwrite(Mhat,"%s/fdpcg_Mhat",dirsetup);
    }
    /*Now invert each block. */
    /*bs: blocksize. */
    long bs=0;
    for(long ips=0; ips<nps; ips++){
	bs+=os[ips]*os[ips];
    }
    long nb=Mhat->m/bs;
    info("Block size is %ld, there are %ld blocks\n",bs,nb);
    /*Permutation vector */
    FDPCG_T *fdpcg=calloc(1, sizeof(FDPCG_T));
    long *perm=fdpcg_perm(nx,ny, os, nps,0,0);
    csp *Mhatp=cspperm(Mhat,0,perm,perm);/*forward permutation. */
    cspfree(Mhat);
    free(perm);
    perm=fdpcg_perm(nx,ny, os, nps, 1,0); /*contains fft shift information*/
    if(parms->save.setup){
	writelong(perm, nxtot, 1, "%s/fdpcg_perm", dirsetup);
    }
#if PRE_PERMUT == 1
    csp *Minvp=cspinvbdiag(Mhatp,bs);
    csp *Minv=cspperm(Minvp,1,perm, perm);/*revert permutation */
    free(perm);
    cspfree(Minvp);
    fdpcg->Minv=Minv;
    if(parms->save.setup){
	cspwrite(Minv,"%s/fdpcg_Minv",dirsetup);
    }
#else
    fdpcg->perm=perm;
    fdpcg->Mbinv=cspblockextract(Mhatp,bs);
    if(parms->save.setup){
	ccellwrite(fdpcg->Mbinv,"%s/fdpcg_Mhatb",dirsetup);
    }
    for(long ib=0; ib<nb; ib++){
	/*2012-04-07: was using inv_inplace that calls gesv that does not truncate svd. In
	  one of the cells the conditional is more than 1e8. This creates
	  problem in GPU code causing a lot of pistion to accumulate and
	  diverges the pcg.

	  2012-04-26: Even with 1e-7 threshold, some blocks has bad conditions
	  that blows the error if we are assemble MVM in GPU.
	*/
	csvd_pow(fdpcg->Mbinv->p[ib], -1, 0, -0.01);
    }
    if(parms->save.setup){
	ccellwrite(fdpcg->Mbinv,"%s/fdpcg_Minvb",dirsetup);
    }
#endif
    cspfree(Mhatp);
    fdpcg->xloc=recon->xloc;
    fdpcg->square=parms->tomo.square;
    fdpcg->scale=calloc(nps,sizeof(double));
    for(int ips=0; ips<nps; ips++){
	fdpcg->scale[ips]=1./(double)(nx[ips]*ny[ips]);
    }

    fdpcg->nxtot=nxtot;
    return fdpcg;
}
typedef struct{
    cmat *xhat;
    cmat *xhat2;
    ccell *xhati;
    ccell *xhat2i;
    FDPCG_T *fdpcg;
    const dcell *xin;
    dcell *xout;
}fdpcg_info_t;
/*typedef struct thread_info{
    int ips;
    int ib;
    pthread_mutex_t lock;
    FDPCG_T *fdpcg;
    const dcell *xin;
    dcell *xout;
}thread_info;
*/
/**
   Copy x vector and do FFT on each layer
*/
static void fdpcg_fft(thread_t *info){
    fdpcg_info_t *data=info->data;
    FDPCG_T *fdpcg=data->fdpcg;
    ccell *xhati=data->xhati;
    const dcell *xin=data->xin;
    for(int ips=info->start; ips<info->end; ips++){
	if(fdpcg->square){
	    for(long i=0; i<xhati->p[ips]->nx*xhati->p[ips]->ny; i++){
		xhati->p[ips]->p[i]=xin->p[ips]->p[i];
	    }
	}else{
	    /*
	      czero takes 0.000074 
	      cembed_locstat takes 0.000037
	    */
	    czero(xhati->p[ips]);
	    cembed_locstat(&xhati->p[ips], 0, fdpcg->xloc[ips],  xin->p[ips]->p, 1, 0);
	}
	cfft2s(xhati->p[ips],-1);
    }
}

/**
   Multiply each block in pthreads
 */
static void fdpcg_mulblock(thread_t *info){
    fdpcg_info_t *data=info->data;
    FDPCG_T *fdpcg=data->fdpcg;
    long bs=fdpcg->Mbinv->p[0]->nx;
    for(int ib=info->start; ib<info->end; ib++){
	cmulvec(data->xhat->p+ib*bs, fdpcg->Mbinv->p[ib], data->xhat2->p+ib*bs,1);
    }
}

/**
   Inverse FFT for each block. Put result in xout, replace content, do not accumulate.
 */
static void fdpcg_ifft(thread_t *info){
    fdpcg_info_t *data=info->data;
    FDPCG_T *fdpcg=data->fdpcg;
    ccell *xhat2i=data->xhat2i;
    dcell *xout=data->xout;
    for(int ips=info->start; ips<info->end; ips++){
	cfft2s(xhat2i->p[ips],1);
	if(fdpcg->square){
	    for(long i=0; i<xhat2i->p[ips]->nx*xhat2i->p[ips]->ny; i++){
		xout->p[ips]->p[i]=creal(xhat2i->p[ips]->p[i]);
	    }
	}else{
	    dzero(xout->p[ips]);
	    cembed_locstat(&xhat2i->p[ips], 1, fdpcg->xloc[ips], xout->p[ips]->p, 0, 1);
	}
    }
}

/**
   Apply the Fourier domain preconditioner. Preconditioner solves xout=M^-1*xin
   where M is the preconditioner matrix. Here M^-1 is applied in Fourier Domain
   via xout=IFFT(M*FFT(xin));
 */
void fdpcg_precond(dcell **xout, const void *A, const dcell *xin){
    //static int count=-1; count++; //temp
    //dcellwrite(xin, "fdpcg_xin_%d", count);
    //dcellwrite(*xout, "fdpcg_xout1_%d", count);
    const RECON_T *recon=(RECON_T*)A;
    if(xin->ny!=1){
	error("Invalid\n");
    }
    const long nps=recon->npsr;
    FDPCG_T *fdpcg=recon->fdpcg;
    long nxtot=fdpcg->nxtot;
    cmat *xhat=cnew(nxtot,1);
    cmat *xhat2=cnew(nxtot,1);
    ccell *xhati=ccellnew(nps,1);/*references the data in xhat. */
    ccell *xhat2i=ccellnew(nps,1);
    long* nx=recon->xnx;
    long* ny=recon->xny;
    long offset=0;
    for(int ips=0; ips<nps; ips++){
	xhati->p[ips]=cnew_ref(nx[ips],ny[ips],xhat->p+offset);
	xhat2i->p[ips]=cnew_ref(nx[ips],ny[ips],xhat2->p+offset);
	cfft2plan(xhati->p[ips],-1);
	cfft2plan(xhat2i->p[ips],1);
	offset+=nx[ips]*ny[ips];
    }
    /*
    thread_info info;
    info.ips=0;
    info.ib=0;
    PINIT(info.lock);
    info.fdpcg=recon->fdpcg;
    info.xin=xin;
    */
    if(!*xout){
	*xout=dcellnew2(xin);
    }
    //info.xout=*xout;
    fdpcg_info_t data={xhat, xhat2, xhati, xhat2i, recon->fdpcg, xin, *xout};
    int NTH=recon->nthread;
    thread_t info_fft[NTH],info_ifft[NTH],info_mulblock[NTH];
    thread_prep(info_fft, 0, nps, NTH, fdpcg_fft, &data);
    thread_prep(info_ifft, 0, nps, NTH, fdpcg_ifft, &data);
    long bs=recon->fdpcg->Mbinv->p[0]->nx;
    long nb=recon->fdpcg->nxtot/bs;
    thread_prep(info_mulblock, 0, nb, NTH, fdpcg_mulblock, &data);
    /*apply forward FFT */
    //CALL(fdpcg_fft,&info,NTH,1);
    CALL_THREAD(info_fft, NTH, 1);
    //ccellwrite(fdpcg->xhati, "fdc_fft_%d", count);
#if PRE_PERMUT
    czero(xhat2);
    cspmulvec(xhat2->p, recon->fdpcg->Minv, xhat->p, 1);
#else/*permute vectors and apply block diagonal matrix */
    /*permute xhat and put into xhat2 */
    cvecperm(xhat2->p,xhat->p,recon->fdpcg->perm,nxtot);
    //ccellwrite(fdpcg->xhat2i, "fdc_perm_%d", count);
    czero(xhat);
    //CALL(fdpcg_mulblock,&info,NTH,1);
    CALL_THREAD(info_mulblock, NTH, 1);
    //ccellwrite(fdpcg->xhati, "fdc_mul_%d", count);
    cvecpermi(xhat2->p,xhat->p,fdpcg->perm,nxtot);
    //ccellwrite(fdpcg->xhat2i, "fdc_iperm_%d", count);
#endif
    //info.ips=0;
    /*Apply inverse FFT */
    //CALL(fdpcg_ifft,&info,NTH,1);
    CALL_THREAD(info_ifft, NTH, 1);
    //ccellwrite(fdpcg->xhat2i, "fdc_ifft_%d", count);
    //dcellwrite(*xout, "fdpcg_xout2_%d", count);
    ccellfree(xhati);
    ccellfree(xhat2i);
    cfree(xhat);
    cfree(xhat2);
}

/**
   Free fdpcg related data structs.
 */
void fdpcg_free(FDPCG_T *fdpcg){
    if(!fdpcg) return;
    cspfree(fdpcg->Minv);
    if(fdpcg->Mbinv){
	free(fdpcg->perm);
	ccellfree(fdpcg->Mbinv);
    }
    free(fdpcg->scale);
    free(fdpcg);
}

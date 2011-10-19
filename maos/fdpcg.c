/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
   Changed to use frequencies from 0 to N-1 to avoid fftshift.
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
   subapertures from ground layer xloc (or ploc).  */
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
	    long ix=(saloc->locx[isa])*dx1+offx;/*subaperture center. */
	    long iy=(saloc->locy[isa])*dx1+offy;
	    pxsel[iy][ix]=1;
	}
    }
    cfftshift(xsel);
    cfft2(xsel,-1);
    //cfftshift(xsel);
    cscale(xsel,1./(double)(nx*ny));/*cancel FFT effect. */
    double xselc=creal(pxsel[ny/2][nx/2])*threas;/*Fourier center */
   
    for(long ix=0; ix<nx; ix++){
	for(long iy=0; iy<ny; iy++){
	    if(cabs(pxsel[iy][ix])<xselc){
		pxsel[iy][ix]=0;
	    }
	}
    }
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

*/
long *fdpcg_perm(const long *nx, const long *ny, long pos, int nps){
    //long nx2[nps],ny2[nps];
    long noff[nps];
    long xloctot=0;
    for(long ips=0; ips<nps; ips++){
	//nx2[ips]=nx[ips]/2;
	//ny2[ips]=ny[ips]/2;
	noff[ips]=xloctot;
	xloctot+=nx[ips]*ny[ips];
    }
    long *perm=calloc(xloctot, sizeof(long));
    long use_os=pos;
    long adim=nx[0]/use_os;
    //long osx=nx[0];
    long count=0;
    //for(long iy=-adim/2; iy<adim/2; iy++){
    //	for(long ix=-adim/2; ix<adim/2; ix++){
    for(long iy=0; iy<adim; iy++){
	for(long ix=0; ix<adim; ix++){
	    /*(ix,iy) is the frequency in SALOC grid. */
	    /*loop over oversampling.*/
	    for(long juse_os=0; juse_os<use_os; juse_os++){
		long jy=(iy+adim*juse_os);
		while(jy>=ny[0]) jy-=ny[0];
		for(long iuse_os=0; iuse_os<use_os; iuse_os++){
		    long jx=(ix+adim*iuse_os);
		    while(jx>=nx[0]) jx-=nx[0];
		    /*jx, jy is the frequency in XLOC grid. */
		    for(long ips=0; ips<nps; ips++){
			/*if this layer has such freq. */
			//if(jy>=-ny2[ips] && jy<ny2[ips] && jx>=-nx2[ips] && jx<nx2[ips]){
			if(jy>=0 && jy<ny[ips] && jx>=0 && jx<nx[ips]){
			    perm[count]=noff[ips]+(jx)+(jy)*nx[ips];
			    count++;
			}
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
void fdpcg_g(cmat **gx, cmat **gy, long nx, long ny, double dx, double dsa){
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
    double dsa2=dsa*0.5;
    for(long iy=0; iy<ny; iy++){
	//double fy=(double)(iy-ny2)*dfy;
	double fy=iy*dfy;
	for(long ix=0; ix<nx; ix++){
	    //double fx=(double)(ix-nx2)*dfx;
	    double fx=ix*dfx;
	    dcomplex tx=0;
	    dcomplex ty=0;
	    dcomplex offset=1;
	    for(int ios=0; ios<os+1; ios++){
		tx+=wt[ios]*(cexp(cf*(fx*dsa+fy*st[ios]))-cexp(cf*(fy*st[ios])));
		ty+=wt[ios]*(cexp(cf*(fy*dsa+fx*st[ios]))-cexp(cf*(fx*st[ios])));
	    }
	    pgx[ix+iy*nx]=offset*tx;
	    pgy[ix+iy*nx]=offset*ty;
	}
    }
}

/**
   Propagate operator for nlayer screens to ground of size
   nxg*nxg, sampling dx, with displacement of dispx, dispy.
*/
csp *fdpcg_prop(long nps, long pos, const int *os, long nxg, double dx, double *dispx, double *dispy){
    long nxi[nps],nxi2[nps],nxi3[nps],noff[nps];
    long nxtot=0;
    double dk=1./(nxg*dx);
    for(long ips=0; ips<nps; ips++){
	nxi[ips]=nxg/pos*os[ips];
	nxi2[ips]=nxi[ips]/2;
	nxi3[ips]=nxi[ips]/2+nxi[ips];
	noff[ips]=nxtot;
	nxtot+=nxi[ips]*nxi[ips];
    }
    long nxg2=nxg/2;
    /*We build the transposed matrix.*/
    csp *propt=cspnew(nxtot,nxg*nxg,nxg*nxg*nps);
    spint *pp=propt->p;
    spint *pi=propt->i;
    dcomplex *px=propt->x;
    long count=0;
    dcomplex cf=2*M_PI*I;
    double cfr=2*M_PI;
    for(long iy=0; iy<nxg; iy++){
	for(long ix=0; ix<nxg; ix++){
	    //double fxg=(ix-nxg2)*dk;/*spatial frequency in ground layer. */
	    //double fyg=(iy-nxg2)*dk;
	    double fxg=ix*dk;
	    double fyg=iy*dk;
	    long icol=ix+iy*nxg;
	    pp[icol]=count;
	    for(long ips=0; ips<nps; ips++){
		//long jx=((ix-nxg2)+nxi3[ips])%nxi[ips];/*map to layer ips. */
		//long jy=((iy-nxg2)+nxi3[ips])%nxi[ips];/*map to layer ips. */
		long jx=ix%nxi[ips];
		long jy=iy%nxi[ips];
		//double fx=(jx-nxi2[ips])*dk;/*spatial frequency in plane ips. */
		//double fy=(jy-nxi2[ips])*dk;
		double fx=jx*dk;
		double fy=jy*dk;
		pi[count]=jx+jy*nxi[ips]+noff[ips];
		dcomplex shift=cexp(cf*(fx*dispx[ips]+fy*dispy[ips]));
		switch(pos/os[ips]){
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
		    error("Invalid\n");
		}
		count++;
	    }
	}
    }
    pp[nxg*nxg]=count;
    /*we put conj above because csptrans applies a conjugation. */
    csp *propf=csptrans(propt);
    cspfree(propt);
    return propf;
}

/**
  Prepare data for Tomography Fourier Domain Preconditioner
*/
FDPCG_T *fdpcg_prepare(const PARMS_T *parms, const RECON_T *recon, const POWFS_T *powfs){
    loc_t **xloc=recon->xloc;
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
    const double hs=parms->powfs[hipowfs].hs;
    const long nps=recon->npsr;
    long pos=parms->tomo.pos;
    const int* os=parms->atmr.os;
    if(pos!=os[0]){
	error("pupil does not equal to ground layer over sampling. Please change.\n");
    }
    long* nx=recon->xnx;
    long* ny=recon->xny;
    const double *ht=parms->atmr.ht;
    long nxtot=0;
    for(long ips=0; ips<nps; ips++){
	nxtot+=nx[ips]*ny[ips];
	if(os[ips]>pos){
	    error("Layer %ld os is greater than ground\n", ips);
	}
    }
    /*Subaperture selection operator */
    csp *sel=fdpcg_saselect(nx[0],ny[0],xloc[0]->dx, 
			    saloc, powfs[hipowfs].saa->p);
    if(parms->save.setup){
	cspwrite(sel,"%s/fdpcg_sel",dirsetup);
    }
    /*Gradient operator. */
    cmat *gx, *gy;
    fdpcg_g(&gx,&gy,nx[0],ny[0],xloc[0]->dx,saloc->dx);/*tested ok. */
    if(parms->save.setup){
	cwrite(gx,"%s/fdpcg_gx",dirsetup);
	cwrite(gy,"%s/fdpcg_gy",dirsetup);
    }
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
    csp *Mhat=cspnewdiag(nxtot,invpsd,1);
    free(invpsd);
    if(parms->save.setup){
	cspwrite(Mhat,"%s/fdpcg_invpsd",dirsetup);
    }

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
    cfree(gx);
    cfree(gy);
    cspfree(sel);
    double dispx[nps];
    double dispy[nps];
    /*Ray tracing operator for each WFS */
    if(parms->save.setup){
	cspwrite(Mmid,"%s/fdpcg_Mmid",dirsetup);
    }
    for(int jwfs=0; jwfs<parms->powfs[hipowfs].nwfs; jwfs++){
	int iwfs=parms->powfs[hipowfs].wfs[jwfs];
	double neai=recon->neam->p[iwfs];
	info("sanea used for wfs %d is %g mas\n",iwfs, 206265000*neai);
	for(long ips=0; ips<nps; ips++){
	    /*
	      2010-05-28: The cone effect cancels with the cone
	      coordinate. so we are like doing parallel beam
	      propagation. Removed scaling by 1/(1-ht[ips]/hs);
	    */
	    dispx[ips]=ht[ips]*parms->wfs[iwfs].thetax;
	    dispy[ips]=ht[ips]*parms->wfs[iwfs].thetay;
	}
	csp *propx=fdpcg_prop(nps,pos,os,nx[0],xloc[0]->dx,dispx,dispy);
	if(parms->save.setup){
	    cspwrite(propx,"%s/fdpcg_prop_wfs%d",dirsetup,iwfs);
	}
	/*need to test this in spatial domain. */
	cspscale(propx,1./neai);/*prop is not real for off axis wfs. */
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
    long *perm=fdpcg_perm(nx,ny, pos, nps);
    csp *Mhatp=cspperm(Mhat,0,perm,perm);/*forward permutation. */
    cspfree(Mhat);
#if PRE_PERMUT == 1
    csp *Minvp=cspinvbdiag(Mhatp,bs);
    csp *Minv=cspperm(Minvp,1,perm, perm);/*revert permutation */
    cspfree(Minvp);
    free(perm);
    fdpcg->Minv=Minv;
    if(parms->save.setup){
	cspwrite(Minv,"%s/fdpcg_Minv",dirsetup);
    }
#else
    fdpcg->perm=perm;
    fdpcg->Mbinv=cspblockextract(Mhatp,bs);
    for(long ib=0; ib<nb; ib++){
	cinv_inplace(fdpcg->Mbinv->p[ib]);
    }
#endif
    cspfree(Mhatp);

    fdpcg->xhat=cnew(nxtot,1);
    fdpcg->xhat2=cnew(nxtot,1);
    fdpcg->xhati=ccellnew(nps,1);/*references the data in xhat. */
    fdpcg->xhat2i=ccellnew(nps,1);
    fdpcg->xloc=recon->xloc;
    fdpcg->square=parms->tomo.square;
    fdpcg->scale=calloc(nps,sizeof(double));
    offset=0;
    for(int ips=0; ips<nps; ips++){
	fdpcg->xhati->p[ips]=cnew_ref(nx[ips],ny[ips],fdpcg->xhat->p+offset);
	fdpcg->xhat2i->p[ips]=cnew_ref(nx[ips],ny[ips],fdpcg->xhat2->p+offset);
	cfft2plan(fdpcg->xhati->p[ips],-1);
	cfft2plan(fdpcg->xhat2i->p[ips],1);
	fdpcg->scale[ips]=1./(double)(nx[ips]*ny[ips]);
	offset+=nx[ips]*ny[ips];
    }

    fdpcg->nxtot=nxtot;
    return fdpcg;
}
typedef struct thread_info{
    int ips;
    int ib;
    pthread_mutex_t lock;
    FDPCG_T *fdpcg;
    const dcell *xin;
    dcell *xout;
}thread_info;

/**
   Copy x vector and do FFT on each layer
*/
static void fdpcg_fft(void *data){
    thread_info *info=data;
    int ips;
    FDPCG_T *fdpcg=info->fdpcg;
    int nps=fdpcg->xhati->nx;
    ccell *xhati=fdpcg->xhati;
    const dcell *xin=info->xin;
    while(LOCKADD(ips, info->ips, 1)<nps){
	if(fdpcg->square){
	    for(long i=0; i<xhati->p[ips]->nx*xhati->p[ips]->ny; i++){
		xhati->p[ips]->p[i]=xin->p[ips]->p[i];
	    }
	}else{
	    /*czero takes 0.000074 */
	    /*cembed_locstat takes 0.000037 */
	    /*embedc_in takes 0.000062 */
	    czero(xhati->p[ips]);
	    cembed_locstat(&xhati->p[ips], 0, fdpcg->xloc[ips],  xin->p[ips]->p, 1, 0);
	}
	/*Apply FFT. first fftshift is not necessary. */
	/*cfftshift(xhati->p[ips]);//enable this needs enable the one in fdpcg_ifft. */
	/*cfft2(xhati->p[ips],-1); */
	cfft2s(xhati->p[ips],-1);
	//cfftshift(xhati->p[ips]);
    }
}

/**
   Multiply each block in pthreads
 */
static void fdpcg_mulblock(void *p){
    thread_info *info=p;
    int ib;
    FDPCG_T *fdpcg=info->fdpcg;
    long bs=fdpcg->Mbinv->p[0]->nx;
    long nb=fdpcg->nxtot/bs;
    cmat *xhat=fdpcg->xhat;
    cmat *xhat2=fdpcg->xhat2;
    while(LOCKADD(ib, info->ib, 1) < nb){
	cmulvec(xhat->p+ib*bs, fdpcg->Mbinv->p[ib], xhat2->p+ib*bs,1);
    }
}

/**
   Inverse FFT for each block. Put result in xout, replace content, do not accumulate.
 */
static void fdpcg_ifft(void *p){
    thread_info *info=p;
    int ips;
    FDPCG_T *fdpcg=info->fdpcg;
    int nps=fdpcg->xhati->nx;
    ccell *xhat2i=fdpcg->xhat2i;
    dcell *xout=info->xout;
    /*const dcell *xin=info->xin; */
    while(LOCKADD(ips, info->ips, 1)<nps){
	//cfftshift(xhat2i->p[ips]);
	cfft2s(xhat2i->p[ips],1);
	/*cfftshift(xhat2i->p[ips]);//enable this needs enable the one in fdpcg_fft. */
	if(fdpcg->square){
	    for(long i=0; i<xhat2i->p[ips]->nx*xhat2i->p[ips]->ny; i++){
		xout->p[ips]->p[i]=creal(xhat2i->p[ips]->p[i]);
	    }
	}else{
	    dzero(xout->p[ips]);
	    cembed_locstat(&xhat2i->p[ips], 1, fdpcg->xloc[ips], xout->p[ips]->p, 0, 1);
	    /*embedc_out(xhat2i->p[ips]->p, xout->p[ips]->p, xout->p[ips]->nx, fdpcg->xembed[ips]); */
	}
    }
}

/**
   Apply the Fourier domain preconditioner. Preconditioner solves xout=M^-1*xin
   where M is the preconditioner matrix. Here M^-1 is applied in Fourier Domain
   via xout=IFFT(M*FFT(xin));
 */
void fdpcg_precond(dcell **xout, const void *A, const dcell *xin){
    const RECON_T *recon=(RECON_T*)A;
    if(xin->ny!=1){
	error("Invalid\n");
    }
    FDPCG_T *fdpcg=recon->fdpcg;
    long nxtot=fdpcg->nxtot;
    cmat *xhat=fdpcg->xhat;
    cmat *xhat2=fdpcg->xhat2;
    thread_info info;
    info.ips=0;
    info.ib=0;
    PINIT(info.lock);
    info.fdpcg=recon->fdpcg;
    info.xin=xin;
    if(!*xout){
	*xout=dcellnew2(xin);
    }
    info.xout=*xout;
    /*apply forward FFT */
    CALL(fdpcg_fft,&info,recon->nthread,1);
    ccellwrite(recon->fdpcg->xhati, "fdpcg_fft");
    if(recon->fdpcg->Minv){/*use sparse matrix */
	czero(xhat2);
	cspmulvec(xhat2->p, recon->fdpcg->Minv, xhat->p, 1);
    }else{/*permute vectors and apply block diagonal matrix */
	/*permute xhat and put into xhat2 */
	cvecperm(xhat2->p,xhat->p,recon->fdpcg->perm,nxtot);
	cwrite(xhat2, "fdpcg_perm");
	czero(xhat);
	CALL(fdpcg_mulblock,&info,recon->nthread,1);
	cwrite(xhat, "fdpcg_mul");
	/*permute back to have natural order. */
	cvecpermi(xhat2->p,xhat->p,fdpcg->perm,nxtot);
	ccellwrite(fdpcg->xhat2i, "fdpcg_perm_i");
    }
    info.ips=0;
    /*Apply inverse FFT */
    CALL(fdpcg_ifft,&info,recon->nthread,1);
    ccellwrite(fdpcg->xhat2i, "fdpcg_ifft");exit(0);
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
    ccellfree(fdpcg->xhati);
    ccellfree(fdpcg->xhat2i);
    cfree(fdpcg->xhat);
    cfree(fdpcg->xhat2);
    free(fdpcg->scale);
    free(fdpcg);
}

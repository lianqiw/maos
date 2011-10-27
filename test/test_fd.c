/**
   Try to implement the FDPCG. Test the subroutines here before committing.
   A few requirements for the implementation:
   The spatial size of all screens must be identical. They may have different sampling though.
   
*/
#include "../lib/aos.h"

static csp* fdpcg_sa(loc_t *xloc, loc_t *saloc, double *saa){
    /**
       Create aperture selection function selects the gradients
       for valid subapertures from ground layer xloc (or ploc).
     */
    const long threas=1;

    long nxloc=xloc->nloc;
    long os=(long)round(saloc->dx/xloc->dx);
    if(fabs(os*xloc->dx-saloc->dx)>1.e-10){
	warning("saloc->dx=%g is not multiple of xloc->dx=%g\n",
		saloc->dx,xloc->dx);
    }
    long nx=(long)sqrt((double)nxloc);
    long ny=nx;
    if(nxloc!=nx*ny){
	warning("xloc is not square. The code may not work\n");
    }
    cmat *xsel=cnew(nx,ny);
    cfft2plan(xsel,-1);
    PCMAT(xsel,pxsel);
    double dx1=1./xloc->dx;
    long offx=-xloc->locx[0]*dx1+saloc->dx*0.5*dx1;
    long offy=-xloc->locy[0]*dx1+saloc->dx*0.5*dx1;
    for(long isa=0; isa<saloc->nloc; isa++){
	if(saa[isa]>0.9){
	    long ix=(saloc->locx[isa])*dx1+offx;/*subaperture center. */
	    long iy=(saloc->locy[isa])*dx1+offy;
	    pxsel[iy][ix]=1;
	}
    }
    cwrite(xsel,"xsel0");
    cfftshift(xsel);
    cfft2(xsel,-1);
    cfftshift(xsel);
    /*temporary. cancel fft effect and compare with laos. */
    cscale(xsel,1./(double)(nx*ny));
    double xselc=creal(pxsel[ny/2][nx/2])*threas;/*Fourier center */
   
    for(long ix=0; ix<nx; ix++){
	for(long iy=0; iy<ny; iy++){
	    if(cabs(pxsel[iy][ix])<xselc){
		pxsel[iy][ix]=0;
	    }
	}
    }
    cwrite(xsel,"xsel1");
    csp *sel=cspconvolvop(xsel);
    cfree(xsel);
    return sel;
}
static long *fdpcg_perm(long *nperm, loc_t **xloc, int nps, loc_t *saloc){
    /**
       Create a permulation acting on xhat (fft of x on xloc) so that points of
       the same spatial frequency in all layers are grouped together.

       In this implementation, it implies that all the screens
       have the same spatial sampling frequency, or they have the
       same tot spatial size.

       The ray tracing couples the screens but only on points
       with same or harmonic spatial frequencies. Group them together.

    */
    long *xdim=calloc(nps, sizeof(long));
    long *os=calloc(nps, sizeof(long));
    long *fxlim=calloc(nps, sizeof(long));
    long *noff=calloc(nps, sizeof(long));
    long xloctot=0;
    for(long ips=0; ips<nps; ips++){
	xdim[ips]=(long)sqrt((double)xloc[ips]->nloc);
	fxlim[ips]=xdim[ips]/2;
	noff[ips]=xloctot;
	xloctot+=xloc[ips]->nloc;
	if(xloc[ips]->nloc!=xdim[ips]*xdim[ips]){
	    error("xloc must be square\n");
	}
	os[ips]=(long)(saloc->dx/xloc[ips]->dx);
	if(fabs(os[ips]*xloc[ips]->dx-saloc->dx)>1.e-10){
	    warning("saloc->dx=%g is not multiple of xloc->dx=%g\n",saloc->dx,xloc[ips]->dx);
	}
	if(os[ips]>os[0]){
	    error("Layer %ld oversampling ratio is greater than ground layer\n",ips);
	}
    }
    long *perm=calloc(xloctot, sizeof(long));
    long use_os=os[0];
    long adim=xdim[0]/use_os;
    long osx=xdim[0]/2;
    long count=0;
    for(long iy=-adim/2; iy<adim/2; iy++){
	for(long ix=-adim/2; ix<adim/2; ix++){
	    /*(ix,iy) is the frequency in SALOC grid. */
	    for(long juse_os=0; juse_os<use_os; juse_os++){
		long jy=(iy+adim*juse_os);
		if(jy>=osx) jy-=xdim[0];
		for(long iuse_os=0; iuse_os<use_os; iuse_os++){
		    long jx=(ix+adim*iuse_os);
		    if(jx>=osx) jx-=xdim[0];

		    /*jx, jy is the frequency in XLOC grid. */
		    for(long ips=0; ips<nps; ips++){
			if(jy>=-fxlim[ips] && jy<fxlim[ips] /*this layer has such freq. */
			   && jx>=-fxlim[ips] && jx<fxlim[ips]){
			    perm[count]=noff[ips]+(jx+fxlim[ips])+(jy+fxlim[ips])*xdim[ips];
			    count++;
			}
		    }
		}
	    }
	}
    }
    free(fxlim);
    free(noff);
    free(xdim);
    free(os);
    *nperm=xloctot;
    return perm;
}
static void fdpcg_g(cmat **gx, cmat **gy, long nx, long ny, double dx, double dsa){
    /**
       Compute gradient operator in Fourier domain
     */
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
	double fy=(double)(iy-ny2)*dfy;
	for(long ix=0; ix<nx; ix++){
	    double fx=(double)(ix-nx2)*dfx;
	    dcomplex tx=0;
	    dcomplex ty=0;
	    dcomplex offset=0;
	    if(os>1){
		offset=cexp(-cf*(fx+fy)*dsa2);/*shift by half a subaperture */
	    }
	    for(int ios=0; ios<os+1; ios++){
		tx+=wt[ios]*(cexp(cf*(fx*dsa+fy*st[ios]))-cexp(cf*(fy*st[ios])));
		ty+=wt[ios]*(cexp(cf*(fy*dsa+fx*st[ios]))-cexp(cf*(fx*st[ios])));
	    }
	    pgx[ix+iy*nx]=offset*tx;
	    pgy[ix+iy*nx]=offset*ty;
	}
    }
}
static csp *fdpcg_prop(long nps, const long *os, long nxg, double dx, 
		    double *dispx, double *dispy){
    /**
       Propagate operator for nlayer screens to ground of size
       nxg*nxg, sampling dx, with displacement of dispx, dispy.
    */
    long nxi[nps],nxi2[nps],nxi3[nps],noff[nps];
    long nxtot=0;
    double dk=1./(nxg*dx);
    for(long ips=0; ips<nps; ips++){
	nxi[ips]=nxg/os[0]*os[ips];
	nxi2[ips]=nxi[ips]/2;
	nxi3[ips]=nxi[ips]/2+nxi[ips];
	noff[ips]=nxtot;
	nxtot+=nxi[ips]*nxi[ips];
    }
    long nxg2=nxg/2;
    csp *propt=cspnew(nxtot,nxg*nxg,nxg*nxg*nps);
    spint *pp=propt->p;
    spint *pi=propt->i;
    dcomplex *px=propt->x;
    long count=0;
    dcomplex cf=2*M_PI*I;
    for(long iy=0; iy<nxg; iy++){
	for(long ix=0; ix<nxg; ix++){
	    long icol=ix+iy*nxg;
	    pp[icol]=count;
	    for(long ips=0; ips<nps; ips++){
		long jx=((ix-nxg2)+nxi3[ips])%nxi[ips];/*map to layer ips. */
		long jy=((iy-nxg2)+nxi3[ips])%nxi[ips];/*map to layer ips. */
		double fx=(jx-nxi2[ips])*dk;
		double fy=(jy-nxi2[ips])*dk;
		pi[count]=jx+jy*nxi[ips]+noff[ips];
		dcomplex shift=cexp(cf*(fx*dispx[ips]+fy*dispy[ips]));
		switch(os[0]/os[ips]){
		case 1:
		    px[count]=conj(shift);
		    break;
		case 2:
		    {
			dcomplex shiftx=cexp(cf*(fx*dx));
			dcomplex shifty=cexp(cf*(fy*dx));
			px[count]=conj(shift*(1+0.25*
					      (shiftx*shifty
					       +conj(shiftx)*shifty
					       +shiftx*conj(shifty)
					       +conj(shiftx)*conj(shifty))));
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
    csp *prop0=csptrans(propt);
    cspfree(propt);
    return prop0;
}
int main(){
    loc_t *xloc=mksqloc_auto(256,256,0.25);
    loc_t *saloc=locread("saloc.bin");
    saloc->dx=0.5;
    long nps=6;
    loc_t **xlocs=calloc(nps, sizeof(loc_t*));
    long os[nps];
    for(long ips=0; ips<nps; ips++){
	xlocs[ips]=xloc;
	os[ips]=2;
    }
    double ht[]={0,1000,2000,4000,8000,16000};
    double hs=90000;
    double thetax[]={0,0,-0.0001613797784,-9.973820086e-05,9.973820086e-05,0.0001613797784};
    double thetay[]={0,0.000169684629,5.243545924e-05,-0.0001372777737,-0.0001372777737,5.243545924e-05};
    double NEA[]={1.712908313e-07,1.712908313e-07,1.712908313e-07,1.712908313e-07,1.712908313e-07,1.712908313e-07};
    double wt[]={0.383615641256048,0.161455623630726,0.065928670909859,0.140502287357995,0.121599369179309,0.126898407666063};
    dmat *saa=dread("SAA.bin");
    csp *sel=fdpcg_sa(xloc,saloc,saa->p);/*Tested. fully agree with laos mkapfd */
    cspwrite(sel,"sel.bin");

    long nperm;
    long *perm=fdpcg_perm(&nperm,xlocs,nps,saloc);/*tested ok. */
    writelong(perm, nperm,1,"perm");
    cmat *gx, *gy;
    fdpcg_g(&gx,&gy,256,256,xloc->dx,saloc->dx);/*tested ok. */
    cwrite(gx,"gx");
    cwrite(gy,"gy");
    double dispx[nps],dispy[nps];

    double r0=0.198749305619780;
    dcomplex *invpsd=calloc(256*256*6, sizeof(dcomplex));
    long offset=0;
    for(long ips=0; ips<nps; ips++){
	dmat *tmp=turbpsd(256,256,0.25*(1-ht[ips]/hs), r0*pow(wt[ips],-3./5.),30,-1);
	dfftshift(tmp);
	dscale(tmp,256*256);
	for(int i=0; i<256*256; i++){
	    invpsd[offset+i]=tmp->p[i];
	}
	offset+=256*256;
	dwrite(tmp,"invpsd_%ld",ips);
	dfree(tmp);
    }
    int load=0;
    
    csp *Mhatp;
    if(load){
	Mhatp=cspread("Mhatperm");
    }else{
	csp *Mpsd=cspnewdiag(256*256*6,invpsd,1);
	cspwrite(Mpsd,"Mpsd.bin");

	/*csp *Mhat=cspnewdiag(256*256*6,invpsd,1); */
	/*cspwrite(Mhat,"cinvpsd"); */
	/*Compute gx'*sel'*sel*gx+gy'*sel'*sel*gy */
	csp *Mmid=NULL;
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
	cspwrite(Mmid,"Mmid");
	TIC;
	tic;
	csp *Mhat=NULL;
	for(int iwfs=0; iwfs<6; iwfs++){
	    for(long ips=0; ips<nps; ips++){
		dispx[ips]=ht[ips]*thetax[iwfs]/(1-ht[ips]/hs);
		dispy[ips]=ht[ips]*thetay[iwfs]/(1-ht[ips]/hs);
	    }
	    /*tested ok if all layers share the same sampling. */
	    csp *prop1=fdpcg_prop(nps,os,256,0.25,dispx,dispy);
	    cspwrite(prop1,"prop_%d",iwfs);

	    /*need to test this in spatial domain. */
	    toc("prop");
	    cspscale(prop1,1./NEA[iwfs]);
	    /*Compute prop'*Mmid*prop and add to Mhat; */
	    csp *tmp=cspmulsp(Mmid,prop1);
	    toc("Mul1");
	    csp *tmp2=csptmulsp(prop1,tmp);
	    toc("Mul2");
	    cspwrite(tmp2,"Mhat_%d.bin",iwfs);
	    cspadd(&Mhat,tmp2);
	    cspfree(tmp);
	    cspfree(tmp2);
	    toc("done");
	    cspfree(prop1);
	    toc("done");
	}
	cspwrite(Mhat,"Mhat_prop.bin");

	cspadd(&Mhat,Mpsd);
	cspfree(Mpsd);
	cspwrite(Mhat,"Mhat.bin");/*Verified with manual computation in MATLAB so far. */
	cspdropeps(Mhat);
	cspsym(Mhat);
	Mhatp=cspperm(Mhat,0,perm,perm);
	cspfree(Mhat);
	cspwrite(Mhatp,"Mhatperm.bin");/*Verified so far. */
    }
    /*Now invert each block. */
    /*First blocksize. */
    long bs=0;
    for(long ips=0; ips<nps; ips++){
	if(os[0]==2){
	    if(os[ips]==2){
		bs+=4;
	    }else if(os[ips]==1){
		bs+=1;
	    }else{
		error("Invalid");
	    }
	}else{
	    bs+=1;
	}
    }
    long nb=Mhatp->m/bs;
    info("Block size is %ld, there are %ld blocks\n",bs,nb);
    csp *Minv=cspinvbdiag(Mhatp,bs);
   
    cspwrite(Minv,"Minv.bin");
    free(perm);
    locfree(xloc);
    locfree(saloc);
    free(xlocs);
    cspfree(Mhatp);
}

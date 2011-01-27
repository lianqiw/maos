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

#define USE_COMPLEX
#include "mat.c"
#include "cell.c"
#include "matbin.c"


/**
   \file cmat.c
   The following are functions that are only useful for
   cmat. notice that fma and modf are slower than native
   functions.
 */

/**
   compute maximum of abs of the numbers
*/
double cmaxabs(const cmat *A){
    double max,min,sum;
    maxmincmp(A->p,A->nx*A->ny,&max,&min,&sum);
    return max;
}
/**
   compute minimum of abs of the numbers
*/
double cminabs(const cmat *A){
    double max,min,sum;
    maxmincmp(A->p,A->nx*A->ny,&max,&min,&sum);
    return min;
}
/**
   compute sum of abs of the numbers
*/
double csumabs(const cmat *A){
    double max,min,sum;
    maxmincmp(A->p,A->nx*A->ny,&max,&min,&sum);
    return sum;
}

/**
   component-wise multiply of three matrices.
   A=A.*B.*C
*/
void ccwm3(cmat *restrict A, const cmat *restrict B, const cmat *restrict C){
    if(!B){
	ccwm(A,C);
    }else if (!C){
	ccwm(A,B);
    }else{
	assert(A && B && C);
	assert(A->nx==B->nx && A->nx==C->nx&&A->ny==B->ny && A->ny==C->ny);
	//component-wise multiply A=A.*B
	const size_t ntot=A->nx*A->ny;
	for(size_t i=0; i<ntot; i++){
	    A->p[i]=A->p[i]*B->p[i]*C->p[i];
	}
    }
}

/**
   component-wise multiply A=A.*B*alpha
*/
void ccwm2(cmat *restrict A, const cmat *restrict B, const double alpha){
    assert(A && A->p);
    const size_t ntot=A->nx*A->ny;
    if(B){
	if(fabs(alpha-1)>1.e-15){
	    for(size_t i=0; i<ntot; i++){
		A->p[i]*=B->p[i]*alpha;
	    }
	}else{
	    for(size_t i=0; i<ntot; i++){
		A->p[i]*=B->p[i];
	    }
	}
    }else{
	if(fabs(alpha-1)>1.e-15){	
	    for(size_t i=0; i<ntot; i++){
		A->p[i]*=alpha;
	    }
	}
    }
}
/**
   Component-wise multiply each column of A with B
   A(:,i)=A(:,i).*B;
 */
void ccwmcol(cmat *restrict A, const cmat *restrict B){
    if (!B) return;
    assert(A->nx==B->nx && B->ny==1);
    dcomplex (*As)[A->nx]=(dcomplex(*)[A->nx])A->p;
    dcomplex *B1=B->p;
    for(unsigned long iy=0; iy<A->ny; iy++){
	for(unsigned long ix=0; ix<A->nx; ix++){
	    As[iy][ix]*=B1[ix];
	}
    }
}

/**
   component wise multiply of 2d complex matrix A,W and 1d vector B.
   A(:,i)=A(:,i).*W(:,i).*B;
*/
void ccwm3col(cmat *restrict A,const cmat *restrict W,const cmat *restrict B){

    if(!W){
	ccwmcol(A,B);
    }else{
	assert(A->nx==B->nx&& A->nx==W->nx&&A->ny==W->ny&&B->ny==1);
	dcomplex (*As)[A->nx]=(dcomplex(*)[A->nx])A->p;
	dcomplex (*Ws)[W->nx]=(dcomplex(*)[W->nx])W->p;
	dcomplex *B1=B->p;

	for(unsigned long iy=0; iy<A->ny; iy++){
	    for(unsigned long ix=0; ix<A->nx; ix++){
		As[iy][ix]=As[iy][ix]*Ws[iy][ix]*B1[ix];
	    }
	}
    }
}
/**
   Component wise multiply each row of A with B.
   A(i,:)=A(i,:)*B
*/
void ccwmrow(cmat *restrict A, const cmat *restrict B){
    if(!A || !B) return;
    dcomplex (*As)[A->nx]=(dcomplex(*)[A->nx])A->p;
    dcomplex *B1=B->p;
    assert(A->ny==B->nx && B->ny==1);
    for(int iy=0; iy<A->ny; iy++){
	dcomplex junk=B1[iy];
	for(int ix=0; ix<A->nx; ix++){
	    As[iy][ix]*=junk;
	}
    }
}
/**
   component-wise multiply of columns of A with combination of B1 and B2:
   A(:,i)=A(:,i)*(B1*wt1+B2*wt2);
 */
void ccwmcol2(cmat *restrict A, 
	      const dcomplex *restrict B1, const double wt1,
	      const dcomplex *restrict B2, const double wt2){
    assert(A && A->p); 
    assert(B1);
    dcomplex (*As)[A->nx]=(dcomplex(*)[A->nx])A->p;
    if(B2){
	for(int ix=0; ix<A->nx; ix++){
	    dcomplex junk=B1[ix]*wt1+B2[ix]*wt2;
	    for(int iy=0; iy<A->ny; iy++){
		As[iy][ix]*=junk;
	    }
	}
    }else{
	for(int ix=0; ix<A->nx; ix++){
	    dcomplex junk=B1[ix]*wt1;
	    for(int iy=0; iy<A->ny; iy++){
		As[iy][ix]*=junk;
	    }
	}
    }
}
/**
   component-wise multiply of rows of A with combination of B1 and B2:
   A(i,:)=A(i,:).*(B1*wt1+B2*wt2);
 */
void ccwmrow2(cmat *restrict A, 
	     const dcomplex *restrict B1, const double wt1,
	     const dcomplex *restrict B2, const double wt2){
    assert(A && A->p); 
    assert(B1);
    dcomplex (*As)[A->nx]=(dcomplex(*)[A->nx])A->p;
    if(B2){
	for(int iy=0; iy<A->ny; iy++){
	    dcomplex junk=B1[iy]*wt1+B2[iy]*wt2;
	    for(int ix=0; ix<A->nx; ix++){
		As[iy][ix]*=junk;
	    }
	}
    }else{
	for(int iy=0; iy<A->ny; iy++){
	    dcomplex junk=B1[iy]*wt1;
	    for(int ix=0; ix<A->nx; ix++){
		As[iy][ix]*=junk;
	    }
	}
    }
}
/**
   component-wise multiply A with conjugate of B: A=A.*conj(B)*alpha;
 */
void ccwmc(cmat *restrict A, const cmat *restrict B, const double alpha){
    assert(A && A->p);
    const size_t ntot=A->nx*A->ny;
    if(B){
	if(fabs(alpha-1)>1.e-15){
	    for(size_t i=0; i<ntot; i++){
		A->p[i]*=conj(B->p[i])*alpha;
	    }
	}else{
	    for(size_t i=0; i<ntot; i++){
		A->p[i]*=conj(B->p[i]);
	    }
	}
    }else{
	if(fabs(alpha-1)>1.e-15){	
	    for(size_t i=0; i<ntot; i++){
		A->p[i]*=alpha;
	    }
	}
    }
}
/**
   component multiply dmat A with cmat B. A=A.*B;
 */
void ccwmd(cmat *restrict A, const dmat *restrict B, const double alpha){
    assert(A && A->p);
    const size_t ntot=A->nx*A->ny;
    if(B){
	if(fabs(alpha-1)>1.e-15){
	    for(size_t i=0; i<ntot; i++){
		A->p[i]*=B->p[i]*alpha;
	    }
	}else{
	    for(size_t i=0; i<ntot; i++){
		A->p[i]*=B->p[i];
	    }
	}
    }else{
	if(fabs(alpha-1)>1.e-15){	
	    for(size_t i=0; i<ntot; i++){
		A->p[i]*=alpha;
	    }
	}
    }
}

/**
   reshape amp.*exp(-2*pi/lambda*opd) into square and embed into center of A,
   with optional rotation of theta CW.

   benchmark: embed 32x32 into 64x64,
   without rotation takes 0.161 ms.
   with rotation takes 0.330 ms.
*/
void cembed_wvf(cmat *restrict A, const double *opd, const double *amp,
	       const int nopdx, const int nopdy, 
	       const double wvl, const double theta){
    dcomplex *psf=A->p;
    const int npsfx=A->nx;
    const int npsfy=A->ny;
 
    dcomplex wvk=2.*M_PI/wvl*I;
    memset(psf, 0, sizeof(dcomplex)*npsfx*npsfy);
    if(fabs(theta)<1.e-10){//no rotation.
	const int skipx=(npsfx-nopdx)/2;
	const int skipy=(npsfy-nopdy)/2;
	assert(skipx>=0 && skipy>=0);
	dcomplex *psf0=psf+skipy*npsfx+skipx;
	for(int iy=0; iy<nopdy; iy++){
	    dcomplex *psfi=psf0+iy*npsfx;
	    const double *opdi=opd+iy*nopdx;
	    const double *ampi=amp+iy*nopdx;
	    for(int ix=0; ix<nopdx; ix++){
		//sqrt(1-cos^2)!=sin due to sign problem. 
		psfi[ix]=ampi[ix]*cexp(wvk*opdi[ix]);
		//most time is spend in cexp.
	    }
	}
    }else{
	//rotated for LGS.
	//rotate image CCW to -theta. coordinate rotate in reverse way.
	//xnew=x*cos(theta)-y*sin(theta)
	//ynew=x*sin(theta)+y*cos(theta);
	//rotation is around the fft center: (nx/2, ny/2);
	/*
	  The original method of adding in to out is not right.
	*/
	dcomplex (*psfs)[npsfx]=(dcomplex(*)[npsfx])psf;
	double (*amps)[nopdx]=(double(*)[nopdx])amp;
	double (*opds)[nopdx]=(double(*)[nopdx])opd;
	const double ctheta=cos(theta);
	const double stheta=sin(theta);
	double x2,y2;
	double x,y;
	int nopdx2=nopdx/2;
	int npsfx2=npsfx/2;
	int nopdy2=nopdy/2;
	int npsfy2=npsfy/2;
	int ix2, iy2;
	double iopd,iamp;
	const int maxr=iceil(sqrt(nopdx*nopdx+nopdy*nopdy));
	int xskip=npsfx>maxr?(npsfx-maxr)/2:0;
	int yskip=npsfy>maxr?(npsfy-maxr)/2:0;
	for(int iy=yskip; iy<npsfy-yskip; iy++){
	    y=iy-npsfy2;
	    for(int ix=xskip; ix<npsfx-xskip; ix++){
		x=ix-npsfx2;
		x2=x*ctheta+y*stheta+nopdx2;
		y2=-x*stheta+y*ctheta+nopdy2;
		if(x2>0 && x2<nopdx-1 && y2>0 && y2<nopdy-1){
		    ix2=ifloor(x2);
		    iy2=ifloor(y2);
		    x2=x2-ix2;
		    y2=y2-iy2;
		    iopd=opds[iy2][ix2]*(1.-x2)*(1.-y2)
			+opds[iy2][ix2+1]*(x2*(1.-y2))
			+opds[iy2+1][ix2]*((1-x2)*y2)
			+opds[iy2+1][ix2+1]*(x2*y2);
		    iamp=amps[iy2][ix2]*(1.-x2)*(1.-y2)
			+amps[iy2][ix2+1]*(x2*(1.-y2))
			+amps[iy2+1][ix2]*((1-x2)*y2)
			+amps[iy2+1][ix2+1]*(x2*y2);
		    psfs[iy][ix]=iamp*cexp(wvk*iopd);
		}
	    }
	}
    }
}
#define cmpcpy(A,B,S) memcpy(A,B,S*sizeof(dcomplex)); /**<macro to do copying*/
/**
   do a modules square and copy to real part of output
*/
inline static void sq2cpy(dcomplex *out, const dcomplex *in, const size_t length){
    for(size_t i=0; i<length; i++){
	out[i]=cabs2(in[i]);
    }
}
/**
   Copy the real part.
*/
inline static void realcpy(dcomplex *out, const dcomplex *in, const size_t length){
    for(size_t i=0; i<length; i++){
	out[i]=creal(in[i]);
    }
}
/**
   Copy the absolute value to output.
 */
inline static void abscpy(dcomplex *out, const dcomplex *in, const size_t length){
    for(size_t i=0; i<length; i++){
	out[i]=cabs(in[i]);
    }
}
#define RA2XY(A) (creal(A)*(cos(cimag(A))+I*sin(cimag(A)))) /**<macro to convert r/a to x/y*/
#define XY2RA(A) (cabs(A)+I*atan2(cimag(A),creal(A))) /**<macro to convert x/y to r/a*/


/**
   Embed array B into A with rotation theta CW.  Current version, preferred */
void cembed(cmat *restrict A, const cmat *restrict B, const double theta, CEMBED flag)
{
    const dcomplex *restrict in=B->p;
    const int ninx=B->nx;
    const int niny=B->ny;
    dcomplex *restrict out=A->p;
    const int noutx=A->nx;
    const int nouty=A->ny;
    /*
      rotate (around fft center: (nx/2,ny/2)) CCW theta and embed in into A. 
      flag==0: pure copy
      flag==1: copy abs2;
      flag==2: copy real only.
    */
    memset(out, 0, sizeof(dcomplex)*noutx*nouty);
    if(fabs(theta)<1.e-10){//no rotation.
	const int skipx=(noutx-ninx)/2;
	const int skipy=(nouty-niny)/2;
	int ixstart=0, ixend=ninx;
	int iystart=0, iyend=niny;
	if(skipx<0){
	    ixstart=-skipx;
	    ixend=ninx+skipx;
	}
	if(skipy<0){
	    iystart=-skipy;
	    iyend=niny+skipy;
	}
	dcomplex *restrict out2=out+skipy*noutx+skipx;
	for(int iy=iystart; iy<iyend; iy++){
	    dcomplex *restrict outi=out2+iy*noutx;
	    const dcomplex *restrict ini=in+iy*ninx;
	    switch(flag){//this switch does not affect speed.
	    case C_FULL:
		cmpcpy(outi+ixstart,ini+ixstart,(ixend-ixstart));
		break;
	    case C_ABS2:
		sq2cpy(outi+ixstart,ini+ixstart,(ixend-ixstart));
		break;
	    case C_REAL:
		realcpy(outi+ixstart,ini+ixstart,(ixend-ixstart));
		break;
	    case C_ABS:
		abscpy(outi+ixstart,ini+ixstart,(ixend-ixstart));
		break;
	    case C_LITERAL://same as 0.
		cmpcpy(outi+ixstart,ini+ixstart,(ixend-ixstart));
		break;
	    default:
		error("Invalid flag\n");
	    }
	}
    }else{
	/*
	  rotated for LGS.
	  rotate image CCW to -theta. coordinate rotate in reverse way.
	  xnew=x*cos(theta)-y*sin(theta)
	  ynew=x*sin(theta)+y*cos(theta);
	  
	  The original method of adding in to out is not right.
	 */
	dcomplex (*restrict outs)[noutx]=(dcomplex(*)[noutx])out;
	dcomplex (*restrict ins)[ninx]=(dcomplex(*)[ninx])in;
	const double ctheta=cos(theta);
	const double stheta=sin(theta);
	const double negstheta=-stheta;
	//const double negctheta=-ctheta;
	//use long double to reduce accumulation error. but slow. error is on order of 1.e-14 
	double x2,y2;
	double x4,y4;
	double x3,y3,x31;
	double noutx2=noutx>>1;
	double nouty2=nouty>>1;
	int ix2, iy2;
	
#define DO_LOOP(AFTER,CMD)						\
	x4=(ninx>>1)-noutx2*ctheta-nouty2*stheta;			\
	y4=(niny>>1)+noutx2*stheta-nouty2*ctheta;			\
	for(int iy=0; iy<nouty; iy++){					\
	    double xbd1=-x4/ctheta; double xbd2=(ninx-1-x4)/ctheta;	\
	    double ybd1=-y4/negstheta; double ybd2=(ninx-1-y4)/negstheta; \
	    if(xbd1>xbd2){double tmp=xbd1; xbd1=xbd2; xbd2=tmp;}	\
	    if(ybd1>ybd2){double tmp=ybd1; ybd1=ybd2; ybd2=tmp;}	\
	    int sx=iceil(fmax(xbd1,ybd1));				\
	    int mx=1+ifloor(fmin(xbd2,ybd2));				\
	    sx=sx>0?sx:0; mx=mx<noutx?mx:noutx;				\
	    x2=x4+ctheta*sx; y2=y4+negstheta*sx;			\
	    for(int ix=sx; ix<mx; ix++){				\
		ix2=ifloor(x2); x3=x2-ix2;x31=1.-x3;			\
		iy2=ifloor(y2); y3=y2-iy2;				\
		outs[iy][ix] =						\
		    AFTER((CMD(ins[iy2][ix2])*(x31)			\
			   +CMD(ins[iy2][ix2+1])*x3)*(1.-y3)		\
			  +(CMD(ins[iy2+1][ix2])*(x31)			\
			    +CMD(ins[iy2+1][ix2+1])*x3)*y3);		\
		x2+=ctheta;						\
		y2+=negstheta;						\
	    }								\
	    x4+=stheta;							\
	    y4+=ctheta;							\
	} 
	//it is not good to embed flag in the inner most loop.
	switch(flag){
	case C_FULL:
	    DO_LOOP(,);
	    break;
	case C_ABS2:
	    DO_LOOP(,cabs2);
	    break;
	case C_REAL:
	    DO_LOOP(,creal);
	    break;
	case C_ABS:
	    DO_LOOP(,cabs);
	    break;
	case C_LITERAL:
	    DO_LOOP(RA2XY,XY2RA);
	    break;
	default:
	    error("Invalid flag\n");
	}
#undef DO_LOOP
    }
}
void cembedd(cmat *restrict A, dmat *restrict B, const double theta){
    double *restrict in=B->p;
    long ninx=B->nx;
    long niny=B->ny;
    dcomplex *out=A->p;
    const long noutx=A->nx;
    const long nouty=A->ny;
    memset(out, 0, sizeof(dcomplex)*noutx*nouty);
    if(fabs(theta)<1.e-10){//no rotation.
	const long skipx=(noutx-ninx)/2;
	const long skipy=(nouty-niny)/2;
	long ixstart=0, ixend=ninx;
	long iystart=0, iyend=niny;
	if(skipx<0){
	    ixstart=-skipx;
	    ixend=ninx+skipx;
	}
	if(skipy<0){
	    iystart=-skipy;
	    iyend=niny+skipy;
	}
	dcomplex *out2=out+skipy*noutx+skipx;
	for(long iy=iystart; iy<iyend; iy++){
	    dcomplex *outi=out2+iy*noutx;
	    const double *ini=in+iy*ninx;
	    for(long ix=ixstart; ix<ixend; ix++){
		outi[ix]=ini[ix];
	    }
	}
    }else{
	dcomplex (*outs)[noutx]=(void*)out;
	double (*ins)[ninx]=(void*)in;
	const double ctheta=cos(theta);
	const double stheta=sin(theta);
	double x2,y2;
	double x,y;
	long ninx2=ninx/2;
	long noutx2=noutx/2;
	long niny2=niny/2;
	long nouty2=nouty/2;
	long ix2, iy2;
	for(long iy=0; iy<nouty; iy++){ 
	    y=(double)(iy-nouty2); 
	    for(long ix=0; ix<noutx; ix++){ 
		x=(double)(ix-noutx2); 
		x2=x*ctheta+y*stheta+ninx2; 
		y2=-x*stheta+y*ctheta+niny2; 
		if(x2>0 && x2<ninx-1 && y2>0 && y2<niny-1){ 
		    ix2=ifloor(x2); 
		    iy2=ifloor(y2); 
		    x2=x2-ix2; 
		    y2=y2-iy2; 
		    outs[iy][ix] =
			ins[iy2][ix2]*((1.-x2)*(1.-y2))
			+ins[iy2][ix2+1]*(x2*(1.-y2))
			+ins[iy2+1][ix2]*((1-x2)*y2)
			+ins[iy2+1][ix2+1]*(x2*y2); 
		} 
	    } 
	} 
    }
}

/**
   rotate (around fft center: (nx/2,ny/2)) CCW theta and embed in into A. 
*/
void cembedscaleout(cmat *restrict A, const cmat *B, 
		    double xoutscale,double youtscale,
		    const double theta, CEMBED flag){
    if(fabs(xoutscale-1)<1.e-10 && fabs(youtscale-1)<1.e-10){
	cembed(A,B,theta,flag);
	return;
    }
    const dcomplex *in=B->p;
    const int ninx=B->nx;
    const int niny=B->ny;
    dcomplex *out=A->p;
    const int noutx=A->nx;
    const int nouty=A->ny;
  

    dcomplex (*outs)[noutx]=(dcomplex(*)[noutx])out;
    dcomplex (*ins)[ninx]=(dcomplex(*)[ninx])in;
    const double ctheta=cos(theta);
    const double stheta=sin(theta);
    double x2,y2;
    double x,y;
    int ninx2=ninx/2;
    int noutx2=noutx/2;
    int niny2=niny/2;
    int nouty2=nouty/2;
    int ix2, iy2;
#define DO_LOOP(AFTER,CMD)						\
    for(int iy=0; iy<nouty; iy++){					\
	y=(iy-nouty2)*youtscale;					\
	for(int ix=0; ix<noutx; ix++){					\
	    x=(ix-noutx2)*xoutscale;					\
	    x2=(x*ctheta+y*stheta+ninx2);				\
	    y2=(-x*stheta+y*ctheta+niny2);				\
	    if(x2>0 && x2<ninx-1 && y2>0 && y2<niny-1){			\
		ix2=ifloor(x2);						\
		iy2=ifloor(y2);						\
		x2=x2-ix2;						\
		y2=y2-iy2;						\
		outs[iy][ix] =AFTER(+CMD(ins[iy2][ix2])*((1.-x2)*(1.-y2)) \
				    +CMD(ins[iy2][ix2+1])*(x2*(1.-y2))	\
				    +CMD(ins[iy2+1][ix2])*((1-x2)*y2)	\
				    +CMD(ins[iy2+1][ix2+1])*(x2*y2));	\
	    }else outs[iy][ix]=0;					\
	}								\
    }
    //it is not good to embed flag in the inner most loop.
    switch(flag){
    case C_FULL:
	DO_LOOP(,);
	break;
    case C_ABS2:
	DO_LOOP(,cabs2);
	break;
    case C_REAL:
	DO_LOOP(,creal);
	break;
    case C_ABS:
	DO_LOOP(,cabs);
	break;
    case C_LITERAL:
	DO_LOOP(RA2XY,XY2RA);
	break;
    default:
	error("Invalid flag\n");
    }
#undef DO_LOOP
}


/**
   copy and embed/crop psfin into psfout
   \verbatim
   4 * * 3
   * * * *
   * * * *
   2 * * 1
   \endverbatim
   to
   \verbatim
   4 * * * 3
   * * * * *
   * * * * *
   * * * * *
   2 * * * 1 
   \endverbatim
   }
*/
void ccpcorner(cmat *A, const cmat *restrict B, CEMBED flag){
    const size_t nx=A->nx;
    const size_t ny=A->ny;
    const size_t ninx=B->nx;
    const size_t niny=B->ny;
    dcomplex *psfout=A->p;
    const dcomplex *restrict psfin=B->p;
    assert((nx&1)==0 && (ny&1)==0 && (ninx&1)==0 && (niny&1)==0);
    memset(psfout, 0, sizeof(dcomplex)*nx*ny);
    const int ny2=(ny<niny)?ny/2:niny/2;
    const int nx2=(nx<ninx)?nx/2:ninx/2;

#define DO_COPY2(CCMD)			\
	CCMD(psfout+(ny-ny2+i)*nx+(nx-nx2),		\
	     psfin+(niny-ny2+i)*ninx+(ninx-nx2),nx2);	\
	CCMD(psfout+(ny-ny2+i)*nx,			\
	     psfin+(niny-ny2+i)*ninx, nx2);		\
	CCMD(psfout+i*nx+(nx-nx2),			\
	     psfin+i*ninx+(ninx-nx2), nx2);		\
	CCMD(psfout+i*nx,				\
	     psfin+i*ninx, nx2); 
    for(int i=0; i<ny2; i++){
	switch(flag){
	case C_FULL:
	    DO_COPY2(cmpcpy);
	    break;
	case C_ABS2:
	    DO_COPY2(sq2cpy);
	    break;
	case C_REAL:
	    DO_COPY2(realcpy);
	    break;
	case C_ABS:
	    DO_COPY2(abscpy);
	    break;
	default:
	    error("Invalid flag\n");
	}	
    }
}
/**
   Put abs square of each element into its realpart.
   A=abs(A).^2;
*/
void cabs2toreal(cmat *A){
    for(int i=0; i<A->nx*A->ny; i++){
	A->p[i]=cabs2(A->p[i]);
    }
}
/**
   Put abs of each elemnt into its realpart. A=abs(A);
*/
void cabstoreal(cmat *A){
    //put abs to real
    for(int i=0; i<A->nx*A->ny; i++){
	A->p[i]=cabs(A->p[i]);
    }
}
/**
   Copy a dmat into real part of cmat.
 */
void ccpd(cmat**restrict A0, const dmat *restrict B){
    cmat *restrict A=*A0;
    if(!A){
	*A0=A=cnew(B->nx, B->ny);
    }else{
	assert((A->nx==B->nx && A->ny==B->ny));
    }
    for(int i=0; i<B->nx*B->ny; i++){
	A->p[i]=B->p[i];
    }
}
/**
   Copy real part of a cmat to dmat with optional scaling:
   A0=A0.*alpha+real(B)*beta
*/
void creal2d(dmat**restrict A0, double alpha,
	     const cmat *restrict B, double beta){
    dmat *restrict A=*A0;
    if(!A){
	*A0=A=dnew(B->nx, B->ny);
    }else{
	assert(A->nx==B->nx && A->ny==B->ny);
    }
    if(fabs(alpha)<EPS){
	memset(A->p, 0,sizeof(double)*B->nx*B->ny);
	for(int i=0; i<B->nx*B->ny; i++){
	    A->p[i]=creal(B->p[i])*beta;
	}
    }else{
	for(int i=0; i<B->nx*B->ny; i++){
	    A->p[i]=A->p[i]*alpha+creal(B->p[i])*beta;
	}
    }
}
/**
   Copy abs squared of a cmat to dmat with optional scaling:
   A0=A0*alpha+abs(B).^2*beta;
 */
void cabs22d(dmat**restrict A0, double alpha,
	     const cmat *restrict B, double beta){
    dmat *restrict A=*A0;
    if(!A){
	*A0=A=dnew(B->nx, B->ny);
    }else{
	assert(A->nx==B->nx && A->ny==B->ny);
    }
    if(fabs(alpha)<1.e-60){
	memset(A->p, 0,sizeof(double)*B->nx*B->ny);
	for(int i=0; i<B->nx*B->ny; i++){
	    A->p[i]=cabs2(B->p[i])*beta;
	}
    }else{
	for(int i=0; i<B->nx*B->ny; i++){
	    A->p[i]=A->p[i]*alpha+cabs2(B->p[i])*beta;
	}
    }
}

/**
   Tilt the otf to make the image shift. 
   apply exp(-2*pi*sx/nx)*exp(-2*pi*sy/ny) to otf to move the image
   sx, sy are shifts in terms of pixel.
   sx/nx is equivalent to sx*dtheta*du.
   pinct=1: peak is in center
   pinct=0: peak is in corner
*/
void ctilt(cmat *otf, double sx, double sy, int pinct){
    int nx=otf->nx;
    int ny=otf->ny;
    double dux=1./(double)nx;
    double duy=1./(double)ny;
    dcomplex ux[nx], uy[ny];
    dcomplex cx=cexp(-2*M_PI*I*dux*sx);
    dcomplex cy=cexp(-2*M_PI*I*duy*sy);
    if(pinct==1){//peak in center
	ux[0]=cexp(-2*M_PI*I*dux*sx*(-nx/2));
	for(int i=1; i<nx; i++){
	    ux[i]=ux[i-1]*cx;
	}
	uy[0]=cexp(-2*M_PI*I*duy*sy*(-ny/2));
	for(int i=1; i<ny; i++){
	    uy[i]=uy[i-1]*cy;
	}
    }else{
	ux[0]=1;
	for(int i=1; i<nx/2; i++){
	    ux[i]=ux[i-1]*cx;
	}
	ux[nx/2]=cexp(-2*M_PI*I*dux*sx*(-nx/2));
	for(int i=nx/2+1; i<nx; i++){
	    ux[i]=ux[i-1]*cx;
	}
	uy[0]=1;
	for(int i=1; i<ny/2; i++){
	    uy[i]=uy[i-1]*cy;
	}
	uy[ny/2]=cexp(-2*M_PI*I*duy*sy*(-ny/2));
	for(int i=ny/2+1; i<ny; i++){
	    uy[i]=uy[i-1]*cy;
	}
    }
    for(int iy=0; iy<ny; iy++){
	dcomplex *p=otf->p+iy*nx;
	for(int ix=0; ix<nx; ix++){
	    p[ix]*=ux[ix]*uy[iy];
	}
    }
}

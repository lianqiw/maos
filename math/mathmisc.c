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


#include <stdint.h>
#include "mathmisc.h"
#include "blas.h"

/**
   Compute the factorial from n1 to n2. Overflow LONG if n2>20, so we use double as output.*/
double factorial(long n1, long n2){
    double fact=1;
    while(n2>=n1){
	fact*=n2--;
    }
    if(!isfinite(fact)){
	error("Factorial overflows\n");
    }
    return fact;
}
/**
   normalize vector to sum to norm;*/
void normalize_sum(double *p, long nloc, double norm){
    if(!nloc) return;
    double ss=norm/dotdbl(p,NULL,NULL,nloc);
    for(int i=0; i<nloc; i++){
	p[i]*=ss;
    }
}
/**
   normalize vector to max to max;*/
void normalize_max(double *p, long nloc, double max){
    if(!nloc) return;
    double ss=max/maxabs(p,nloc);
    for(int i=0; i<nloc; i++){
	p[i]*=ss;
    }
}
/**
   compute sum(p1.*p2.*p3)*/

double dotdbl(const double *restrict p1, const double *restrict p2, 
	   const double *restrict p3, long n){
    double ans=0;
    if(p1 && p2 && p3){
	for(long i=0; i<n; i++) ans+=p1[i]*p2[i]*p3[i];
    }else if(!p1 && p2 && p3){
	for(long i=0; i<n; i++) ans+=p2[i]*p3[i];
    }else if(p1 && !p2 && p3){
	for(long i=0; i<n; i++) ans+=p1[i]*p3[i];
    }else if(!p1 && !p2 && p3){
	for(long i=0; i<n; i++) ans+=p3[i];
    }else if(p1 && p2 && !p3){
	for(long i=0; i<n; i++) ans+=p1[i]*p2[i];
    }else if(!p1 && p2 && !p3){
	for(long i=0; i<n; i++) ans+=p2[i];
    }else if(p1 && !p2 && !p3){
	for(long i=0; i<n; i++) ans+=p1[i];
    }else if(!p1 && !p2 && !p3){
	ans=(double)n;/*assume all ones. */
    }
    return  ans;
}

/**
   compute sum(p1.*p2.*p3)*/

float dotflt(const float *restrict p1, const float *restrict p2, 
	   const float *restrict p3, long n){
    float ans=0;
    if(p1 && p2 && p3){
	for(long i=0; i<n; i++) ans+=p1[i]*p2[i]*p3[i];
    }else if(!p1 && p2 && p3){
	for(long i=0; i<n; i++) ans+=p2[i]*p3[i];
    }else if(p1 && !p2 && p3){
	for(long i=0; i<n; i++) ans+=p1[i]*p3[i];
    }else if(!p1 && !p2 && p3){
	for(long i=0; i<n; i++) ans+=p3[i];
    }else if(p1 && p2 && !p3){
	for(long i=0; i<n; i++) ans+=p1[i]*p2[i];
    }else if(!p1 && p2 && !p3){
	for(long i=0; i<n; i++) ans+=p2[i];
    }else if(p1 && !p2 && !p3){
	for(long i=0; i<n; i++) ans+=p1[i];
    }else if(!p1 && !p2 && !p3){
	ans=(float)n;/*assume all ones. */
    }
    return  ans;
}
/**
   compute sum(p1.*p2.*p3) for complex p1, p2 and double p3*/
dcomplex dotcmp(const dcomplex *restrict p1, const dcomplex *restrict p2, 
		const double *restrict p3, long n){
    dcomplex ans=0;
    if(p1 && p2 && p3){
	for(long i=0; i<n; i++) ans+=p1[i]*p2[i]*p3[i];
    }else if(!p1 && p2 && p3){
	for(long i=0; i<n; i++) ans+=p2[i]*p3[i];
    }else if(p1 && !p2 && p3){
	for(long i=0; i<n; i++) ans+=p1[i]*p3[i];
    }else if(!p1 && !p2 && p3){
	for(long i=0; i<n; i++) ans+=p3[i];
    }else if(p1 && p2 && !p3){
	for(long i=0; i<n; i++) ans+=p1[i]*p2[i];
    }else if(!p1 && p2 && !p3){
	for(long i=0; i<n; i++) ans+=p2[i];
    }else if(p1 && !p2 && !p3){
	for(long i=0; i<n; i++) ans+=p1[i];
    }else if(!p1 && !p2 && !p3){
	ans=(double)n;/*assume all ones. */
    }
    return  ans;
}

/**
   compute sum(p1.*p2.*p3) for complex p1, p2 and float p3*/
fcomplex dotzmp(const fcomplex *restrict p1, const fcomplex *restrict p2, 
		const float *restrict p3, long n){
    fcomplex ans=0;
    if(p1 && p2 && p3){
	for(long i=0; i<n; i++) ans+=p1[i]*p2[i]*p3[i];
    }else if(!p1 && p2 && p3){
	for(long i=0; i<n; i++) ans+=p2[i]*p3[i];
    }else if(p1 && !p2 && p3){
	for(long i=0; i<n; i++) ans+=p1[i]*p3[i];
    }else if(!p1 && !p2 && p3){
	for(long i=0; i<n; i++) ans+=p3[i];
    }else if(p1 && p2 && !p3){
	for(long i=0; i<n; i++) ans+=p1[i]*p2[i];
    }else if(!p1 && p2 && !p3){
	for(long i=0; i<n; i++) ans+=p2[i];
    }else if(p1 && !p2 && !p3){
	for(long i=0; i<n; i++) ans+=p1[i];
    }else if(!p1 && !p2 && !p3){
	ans=(float)n;/*assume all ones. */
    }
    return  ans;
}
/**
   compute the maximum of double vector */
double maxdbl(const double *p, long n){
    if(!p) return 0;
    double a=p[0];
    for(long i=1; i<n; i++){
	if(p[i]>a) a=p[i];
    }
    return a;
}
/**
   compute the maximum of the abs of double vector*/
double maxabs(const double *p, long n){
    if(!p) return 0;
    double a=fabs(p[0]);
    for(long i=1; i<n; i++){
	if(fabs(p[i])>a) a=fabs(p[i]);
    }
    return a;
}
/**
   compute the maximum of the abs of float vector*/
float maxabsf(const float *p, long n){
    if(!p) return 0;
    float a=fabs(p[0]);
    for(long i=1; i<n; i++){
	if(fabsf(p[i])>a) a=fabsf(p[i]);
    }
    return a;
}
/**
   Compute the sum of double vector*/
double dblsum(const double *p, long nx){
    double sum=0;
    for(long ix=0; ix<nx; ix++){
	sum+=p[ix];
    }
    return sum;
}

/**
   Compute the sum of double vector*/
float fltsum(const float *p, long nx){
    float sum=0;
    for(long ix=0; ix<nx; ix++){
	sum+=p[ix];
    }
    return sum;
}

/**
   add two double vectors: out=out*alpha+in*beta+theta
*/
void adddbl(double *restrict out, double alpha, 
	    const double *in, int N, double beta, double theta){

    if (fabs(alpha)<EPS){
	memset(out,0,sizeof(double)*N);
	if(in){
	    for(int i=0; i<N; i++){
		out[i]=in[i]*beta+theta;
	    }
	}else{
	    for(int i=0; i<N; i++){
		out[i]=theta;
	    }
	}
    }else{
	if(in){
	    for(int i=0; i<N; i++){
		out[i]=out[i]*alpha+in[i]*beta+theta;
	    }
	}else{
	    for(int i=0; i<N; i++){
		out[i]=out[i]*alpha+theta;
	    }
	}
    }
}


/**
   Reverse an permution vector.
   A=B(p) <==> A(pi)=B
*/ 
long *invperm(long *p, long np){
    long *restrict pi=malloc(sizeof(long)*np);/*inverse of p */
    for(long irow=0; irow<np; irow++){
	pi[p[irow]]=irow;
    }
    return pi;
}
/**
   Permute the vector so that
   out(:)=in(p);
*/
void cvecperm(dcomplex *restrict out, const dcomplex *restrict in, const long *perm, long nx){
    for(long i=0; i<nx; i++){
	if(perm[i]>0){
	    out[i]=in[perm[i]];
	}else{
	    out[i]=conj(in[-perm[i]]);
	}
    }
}
/**
   Reverse permute the vector so that
   out(p)=in(:);
*/
void cvecpermi(dcomplex *restrict out, const dcomplex *restrict in, const long *perm, long nx){
    for(long i=0; i<nx; i++){
	if(perm[i]>0){
	    out[perm[i]]=in[i];
	}else{
	    out[-perm[i]]=conj(in[i]);
	}
    }
}
/**
 Compute max, min and sum of a long vector*/
void maxminlong(const long *restrict p, long N,
		long *restrict max, long *restrict min){
    if(!p) return;
    long a=p[0];
    long b=p[0];
    long i=0;
    for(i=0; i<N; i++){
	if(p[i]>a) a=p[i];
	if(p[i]<b) b=p[i];
    }
    if(max)*max=a; 
    if(min)*min=b; 
}
/**
   Remove piston from p of length n
*/
void remove_piston(double *p, long n){
    double piston=0;
    for(long i=0; i<n; i++){
	piston+=p[i];
    }
    piston/=-n;
    for(long i=0; i<n; i++){
	p[i]+=piston;
    }
}
/**
   Find the next number that is power of 2.
*/
long nextpow2(long n){
    n--;
    for(long i=1; i<sizeof(long)*8; i<<=1){
	n = n | (n >> i);
    }
    return n+1;
}
/**
   Find the next number suitable for FFT. Use radix of 2, 3, 5, 7
*/
long nextfftsize(long n){
    if(n==0) n=1;
    const int nradix=4;
    const int radixs[]={2,3,5,7};
    int selected[4];
    long n2;/*division*/
    long n3;/*accumulated value so far*/
    do{
	n2=n;
	n3=1;
	memset(selected, 0, nradix*sizeof(int));
	/*We only allow mix of two radices. More will slow down fft*/
	for(int irad=0; irad<nradix; ){
	    int radix=radixs[irad];
	    int ratio=n2/radix;
	    if(ratio*radix==n2){/*no remainder*/
		selected[irad]=1;
		n2=ratio;
		n3*=radix;
	    }else{
		irad++;
	    }
	}
	/*If there is remainder not divisble by 2, 3, 5, or 7. Increase n2 by 1 and test again.*/
	int count=0;
	for(int i=0; i<nradix; i++){
	    count+=selected[i];
	}
	if(count>2){
	    n2=n;
	}
	n++; 
    }while(n2>1 || (n3&1)==1);
    return n3;
}
unsigned long mylog2(unsigned long n){/*find m so that pow(2,m)==n. */
    assert((n & (n-1))==0);
    unsigned long m=-1;
    for(;n;n>>=1){
	m++;
    }
    return m;
}
double golden_section_search(golden_section_fun f, void *param, 
			     double x1, double x4, double tau){
    static double resphi= 0.381966011250105;/*2-0.5*(1+sqrt(5)); */
    double x2=(x4-x1)*resphi+x1;
    double f2=f(param, x2);
    double x3, f3;
    /*stop searching. */
#undef DBG_GF
    do{
	x3=(x4-x2)*resphi+x2;
	f3=f(param, x3);
#ifdef DBG_GF 
	double f1=f(param,x1);
	double f4=f(param,x4);
	info2("x=%.4g, %.4g, %.4g, %.4g\n", x1, x2, x3, x4);
	info2("f=%.4g, %.4g, %.4g, %.4g\n", f1, f2, f3, f4);
#endif
	if(f3<f2){
	    x1=x2;
	    x2=x3;
#ifdef  DBG_GF
	    f1=f2;
#endif
	    f2=f3;
	}else{
	    x4=x1;
	    x1=x3;
#ifdef  DBG_GF
	    f4=f1;
	    f1=f3;
#endif
	}
    }while(fabs(x4-x1) > tau * (fabs(x1)+fabs(x4)));
    /*do not return x1 or x4 which may be in unstable region*/
    return f2<f3?x2:x3;
}


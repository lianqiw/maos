/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#include <math.h>
#include <stdint.h>
#include "mathmisc.h"
#include "blas.h"
/**
   \file mathmisc.c
   A few math routines
*/
/**
   normalize vector to sum to norm;*/
void normalize(double *p, long nloc, double norm){
    double ss=norm/dotdbl(p,NULL,NULL,nloc);
    for(int i=0; i<nloc; i++){
	p[i]*=ss;
    }
}
/**
   normalize vector to max to max;*/
void normalize_max(double *p, long nloc, double max){
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
	ans=(double)n;//assume all ones.
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
	ans=(double)n;//assume all ones.
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
   Compute the sum of double vector*/
double dblsum(double *p, long nx){
    double sum=0;
    for(long ix=0; ix<nx; ix++){
	sum+=p[ix];
    }
    return sum;
}


/**
   inplace inversion of square SPD matrix. A=A^-1*/
void invsq(long n, double *restrict A){
    int N=n;
    int info;
    char uplo='U';
    /*B is identity matrix*/
    double *B=calloc(N*N, sizeof(double));
    for(int i=0;i<N;i++)
	B[i+i*N]=1;
    
    dposv_(&uplo,&N,&N,A,&N,B,&N,&info);
    memcpy(A,B,N*N);
    free(B);
}

/**
   add two double vectors: out=out*alpha+in*beta
*/
void adddbl(double *restrict out, double alpha, 
	    const double *in, int N, double beta){

    if (fabs(alpha)<1.e-60){
	memset(out,0,sizeof(double)*N);
	for(int i=0; i<N; i++){
	    out[i]=in[i]*beta;
	}
    }else{
	for(int i=0; i<N; i++){
	    out[i]=out[i]*alpha+in[i]*beta;
	}
    }
}


/**
   Reverse an permution vector.
   A=B(p) <==> A(pi)=B
*/ 
long *invperm(long *p, long np){
    long *restrict pi=malloc(sizeof(long)*np);//inverse of p
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
	out[i]=in[perm[i]];
    }
}
/**
   Reverse permute the vector so that
   out(p)=in(:);
*/
void cvecpermi(dcomplex *restrict out, const dcomplex *restrict in, const long *perm, long nx){
 
    for(long i=0; i<nx; i++){
	out[perm[i]]=in[i];
    }
}
/**
 Compute max, min and sum of a double vector*/
void maxmindbl(const double *restrict p, long N, 
	       double *restrict max, double *restrict min){
    if(N==0){
	*max=0;
	*min=0;
	return;
    }
    double a,b;
    long i;
    a=-INFINITY;
    b=INFINITY;
    for(i=0; i<N; i++){
	if(p[i]>a) a=p[i];
	if(p[i]<b) b=p[i];
    }
    if(max) *max=a;
    if(min) *min=b;
}
/**
 Compute max, min and sum of a long vector*/
void maxminlong(const long *restrict p, long N,
		long *restrict max, long *restrict min){
    long a,b;
    long i;
    a=INTMAX_MIN;
    b=INTMAX_MAX;
    for(i=0; i<N; i++){
	if(p[i]>a) a=p[i];
	if(p[i]<b) b=p[i];
    }
    if(max)*max=a; 
    if(min)*min=b; 
}
/**
 Compute max, min and sum of a complex vector*/
void maxmincmp(const dcomplex *restrict p, long N,
	       double *restrict max, double *restrict min, double *restrict sum){
    double a,b,s;
    long i;
    a=cabs(p[0]); 
    b=cabs(p[0]);
    s=0;
    for(i=0; i<N; i++){
	double tmp=cabs(p[i]);
	s+=tmp;
	if(tmp>a) a=tmp;
	if(tmp<b) b=tmp;
    }
    if(max)*max=a; 
    if(min)*min=b; 
    if(sum)*sum=s;
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
unsigned long mylog2(unsigned long n){//find m so that pow(2,m)==n.
    assert((n & (n-1))==0);
    unsigned long m=-1;
    for(;n;n>>=1){
	m++;
    }
    return m;
}

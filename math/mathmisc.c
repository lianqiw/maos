/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "numtype.h"
#include "mathmisc.h"
#include "blas.h"
const real AS2RAD=4.848136811095360e-06; //arcsec in unit of radian
const real MAS2RAD=4.848136811095360e-09; //arcsec in unit of radian
const real RAD2AS=206264.8062470964; //radian in unit of arcsec
const real RAD2MAS=206264806.2470964; //radian to milli-arcsecond
#define check_vec(p,N) (p?(N?1:0):(N?(error("p is null but N is not 0\n"),0):0))
/**
   Compute the factorial from n1 to n2. Overflow LONG if n2>20, so we use double as output.*/
long double factorial(long n1, long n2){
	long double fact=1;
	while(n2>=n1){
		fact*=n2--;
	}
	if(!isfinite(fact)){
		error("Factorial overflows\n");
	}
	return fact;
}

/**
   Reverse an permution vector.
   A=B(p) <==> A(pi)=B
*/
long* invperm(long* p, long np){
	if(!check_vec(p, np)) return NULL;
	long* restrict pi=mymalloc(np, long);/*inverse of p */
	for(long irow=0; irow<np; irow++){
		pi[p[irow]]=irow;
	}
	return pi;
}

/**
 Compute max, min and sum of a long vector*/
void maxminlong(const long* restrict p, long N,
	long* restrict max, long* restrict min){
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
   Find the next number that is power of 2.
*/
long nextpow2(long n){
	n--;
	for(size_t i=1; i<sizeof(long)*8; i<<=1){
		n=n|(n>>i);
	}
	return n+1;
}
/**
   Find the next number suitable for FFT. Use radix of 2, 3, 5, 7, 11
*/
long nextfftsize(long n){
	if(n==0) n=1;
	n=(n+1)/2*2;//make it a even number first
#define nradix 4
	const int radixs[nradix]={2,3,5,7};
	int selected[nradix];
	long n2;/*division*/
	long n3;/*accumulated value so far*/
	do{
		n2=n;
		n3=1;
		memset(selected, 0, nradix*sizeof(int));

		for(int irad=0; irad<nradix; ){
			int radix=radixs[irad];
			int ratio=n2/radix;
			if(ratio*radix==n2){/*no remainder*/
				selected[irad]=1;
				n2=ratio;
				n3*=radix;
			} else{
				irad++;
			}
		}
		/*
		  //Limit to mix of two radices. More will slow down fft.
		  int count=0;
		  for(int i=0; i<nradix; i++){
		  count+=selected[i];
		  }
		  if(count>2){
		  n2=n;
		  }*/
		/*If there is remainder not divisble by the radix. Increase n2 by 2 and test again.*/
		n+=2;
	} while(n2>1||(n3&1)==1);
	return n3;
}
unsigned long mylog2(unsigned long n){/*find m so that pow(2,m)==n. */
	assert((n&(n-1))==0);
	unsigned long m=-1;
	for(;n;n>>=1){
		m++;
	}
	return m;
}
real golden_section_search(golden_section_fun f, void* param,
	real x1, real x4, real tau){
	static real resphi=0.381966011250105;/*2-0.5*(1+sqrt(5)); */
	real x2=(x4-x1)*resphi+x1;
	real f2=f(param, x2);
	real x3, f3;
	/*stop searching. */
#undef DBG_GF
	do{
		x3=(x4-x2)*resphi+x2;
		f3=f(param, x3);
#ifdef DBG_GF 
		real f1=f(param, x1);
		real f4=f(param, x4);
		info("x=%.4g, %.4g, %.4g, %.4g\n", x1, x2, x3, x4);
		info("f=%.4g, %.4g, %.4g, %.4g\n", f1, f2, f3, f4);
#endif
		if(f3<f2){
			x1=x2;
			x2=x3;
#ifdef  DBG_GF
			f1=f2;
#endif
			f2=f3;
		} else{
			x4=x1;
			x1=x3;
#ifdef  DBG_GF
			f4=f1;
			f1=f3;
#endif
		}
	} while(fabs(x4-x1)>tau*(fabs(x1)+fabs(x4)));
	/*do not return x1 or x4 which may be in unstable region*/
	return f2<f3?x2:x3;
}

/**
   Read spint array of size nx*ny from file and do optional data conversion.
*/
spint* readspint(file_t* fp, long* nx, long* ny){
	header_t header;
	read_header(&header, fp);
	free(header.str);
	spint* out=NULL;
	if(nx!=0&&ny!=0){
		*nx=(long)header.nx;
		*ny=(long)header.ny;
		size_t len=(*nx)*(*ny);
		out=mymalloc(len, spint);
		if(readvec(out, M_SPINT, header.magic, sizeof(spint), len, fp)){
			free(out);
			return NULL;
		}
	} else{
		*nx=0;
		*ny=0;
	}
	return out;
}
/**
   Read a vector from file and perform conversion if type mismatches.
*/
int readvec(void* p, uint32_t magic_p, uint32_t magic_file, size_t size, size_t nmemb, file_t* fp){
	if(!p||!fp) return -1;
	if(nmemb==0) return 0;
	int ans=0;
	magic_p&=0xFFFF;
	magic_file&=0xFFFF;
#define TEST_CONVERT(M1, T1, M2, T2)		\
    else if(magic_p==M1 && magic_file==M2){	\
		T2 *p2=mymalloc(nmemb,T2);			\
		if(!(ans=zfread(p2, sizeof(T2), nmemb, fp))){\
		for(size_t i=0; i<nmemb; i++){		\
			((T1*)p)[i]=p2[i];				\
		}}							        \
		free(p2);					        \
    }

	if(magic_p==magic_file){
		ans=zfread(p, size, nmemb, fp);
	}
		TEST_CONVERT(M_DBL, double, M_FLT, float)
		TEST_CONVERT(M_DBL, double, M_CMP, dcomplex)
		TEST_CONVERT(M_DBL, double, M_ZMP, fcomplex)
		TEST_CONVERT(M_DBL, double, M_INT32, int)
		TEST_CONVERT(M_DBL, double, M_INT64, int64_t)
		TEST_CONVERT(M_DBL, double, M_SPINT, spint)
		TEST_CONVERT(M_FLT, float, M_DBL, double)
		TEST_CONVERT(M_FLT, float, M_CMP, dcomplex)
		TEST_CONVERT(M_FLT, float, M_ZMP, fcomplex)
		TEST_CONVERT(M_FLT, float, M_INT32, int)
		TEST_CONVERT(M_FLT, float, M_INT64, int64_t)
		TEST_CONVERT(M_FLT, float, M_SPINT, spint)
		TEST_CONVERT(M_CMP, dcomplex, M_DBL, double)
		TEST_CONVERT(M_CMP, dcomplex, M_FLT, float)
		TEST_CONVERT(M_CMP, dcomplex, M_ZMP, fcomplex)
		TEST_CONVERT(M_CMP, dcomplex, M_INT32, int)
		TEST_CONVERT(M_CMP, dcomplex, M_INT64, int64_t)
		TEST_CONVERT(M_CMP, dcomplex, M_SPINT, spint)
		TEST_CONVERT(M_ZMP, fcomplex, M_DBL, double)
		TEST_CONVERT(M_ZMP, fcomplex, M_FLT, float)
		TEST_CONVERT(M_ZMP, fcomplex, M_CMP, dcomplex)
		TEST_CONVERT(M_ZMP, fcomplex, M_INT32, int)
		TEST_CONVERT(M_ZMP, fcomplex, M_INT64, int64_t)
		TEST_CONVERT(M_ZMP, fcomplex, M_SPINT, spint)
		TEST_CONVERT(M_INT32, int, M_DBL, double)
		TEST_CONVERT(M_INT32, int, M_FLT, float)
		TEST_CONVERT(M_INT32, int, M_CMP, dcomplex)
		TEST_CONVERT(M_INT32, int, M_ZMP, fcomplex)
		TEST_CONVERT(M_INT32, int, M_INT64, int64_t)
		TEST_CONVERT(M_INT32, int, M_SPINT, spint)
		TEST_CONVERT(M_INT64, int64_t, M_DBL, double)
		TEST_CONVERT(M_INT64, int64_t, M_FLT, float)
		TEST_CONVERT(M_INT64, int64_t, M_CMP, dcomplex)
		TEST_CONVERT(M_INT64, int64_t, M_ZMP, fcomplex)
		TEST_CONVERT(M_INT64, int64_t, M_INT32, int)
		TEST_CONVERT(M_INT64, int64_t, M_SPINT, spint)
		TEST_CONVERT(M_SPINT, spint, M_DBL, double)
		TEST_CONVERT(M_SPINT, spint, M_FLT, float)
		TEST_CONVERT(M_SPINT, spint, M_CMP, dcomplex)
		TEST_CONVERT(M_SPINT, spint, M_ZMP, fcomplex)
		TEST_CONVERT(M_SPINT, spint, M_INT32, int)
		TEST_CONVERT(M_SPINT, spint, M_INT64, int64_t)
	else{
		ans=-1;
		warning("Please implement conversion from %x to %x\n", magic_file, magic_p);
	}
	return ans;
}


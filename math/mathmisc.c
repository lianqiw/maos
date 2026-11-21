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
#include <stdint.h>
#include "numtype.h"
#include "mathmisc.h"
#include "blas.h"
int fft_disable_threads=0;
#define check_vec(p,N) (p?(N?1:0):(N?(error("p is null but N is not 0\n"),0):0))
/**
   Compute the factorial from n1 to n2, i.e., the product n1*(n1+1)*...*n2. Overflow LONG if n2>20, so we use long
   double as output.

   @param[in] n1 First factor.
   @param[in] n2 Last factor.
   @return The product n1*(n1+1)*...*n2, or 1 if n2<n1.
*/
long double factorial(long n1, long n2){
	long double fact=1;
	while(n2>=n1){
		fact*=n2--;
	}
	if(isinf(fact)){
		error("Factorial overflows\n");
	}
	return fact;
}

/**
   Compute the greatest common divisor.

   @param[in] a One value.
   @param[in] b Another value. The order of a and b does not matter.
   @return The greatest common divisor of a and b (non-negative).
*/
long calc_gcd(long a, long b){
	if(a<0) a=-a;
	if(b<0) b=-b;
	if(a==b || b==0){
		return a;
	}else if(a>b){
		return calc_gcd(b, a%b);
	}else{
		return calc_gcd(b, a);
	}
}
/**
   Compute the least common multiple.

   @param[in] a One value.
   @param[in] b Another value.
   @return The least common multiple of a and b (non-negative), or 0 if both are 0.
*/
long calc_lcm(long a, long b){
	if(a<0) a=-a;
	if(b<0) b=-b;
	if(a==0&&b==0){
		return 0;
	}else{
		return a*(b/calc_gcd(a,b));
	}
}
/**
   Invert a permutation vector.
   A=B(p) <==> A(pi)=B

   @param[in] p The permutation vector.
   @param[in] np Length of p.
   @return The inverse permutation vector, or NULL if p is NULL. Newly allocated; the caller frees it.
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
   Compute max and min of a long vector.

   @param[in] p The vector.
   @param[in] N Length of p.
   @param[out] max The maximum value. Can be NULL. Not modified if p is NULL or empty.
   @param[out] min The minimum value. Can be NULL. Not modified if p is NULL or empty.
*/
void maxminlong(const long* restrict p, long N, long* restrict max, long* restrict min){
	if(!p||N<=0) return;
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

   @param[in] n The input number.
   @return The smallest power of 2 that is no less than n.
*/
long nextpow2(long n){
	n--;
	for(size_t i=1; i<sizeof(long)*8; i<<=1){
		n=n|(n>>i);
	}
	return n+1;
}
/**
   Find the next even number suitable for FFT. Use radix of 2, 3, 5, 7.

   @param[in] n The input number.
   @return The smallest even number no less than n that only has factors of 2, 3, 5, and 7.
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
			long ratio=n2/radix;
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
/**
   Compute integer log2 of n, which must be a power of 2.

   @param[in] n The input number. Must be a power of 2.
   @return m so that pow(2,m)==n.
*/
unsigned long mylog2(unsigned long n){/*find m so that pow(2,m)==n. */
	assert((n&(n-1))==0);
	unsigned long m=-1;
	for(;n;n>>=1){
		m++;
	}
	return m;
}
/**
   Find the minimum of f(param, x) in the interval [x1, x4] using golden section
   search. Stops when the interval is smaller than tau relative to |x1|+|x4|.

   @param[in] f The function to minimize.
   @param[in] param Parameter passed to f.
   @param[in] x1 Lower end of the interval.
   @param[in] x4 Upper end of the interval.
   @param[in] tau Relative tolerance.
   @return The abscissa of the minimum.
*/
real golden_section_search(golden_section_fun f, void* param, real x1, real x4, real tau){
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

   @param[in,out] fp The file to read from.
   @param[out] nx Number of rows read. Can be NULL.
   @param[out] ny Number of columns read. Can be NULL.
   @return The array, or NULL if reading failed. Newly allocated; the caller frees it.
*/
spint* readspint(file_t* fp, long* nx, long* ny){
	header_t header;
	read_header(&header, fp);
	free(header.str);
	spint* out=NULL;
	size_t len=(size_t)header.nx*header.ny;
	if(nx) *nx=len?(long)header.nx:0;
	if(ny) *ny=len?(long)header.ny:0;
	if(len){
		out=mymalloc(len, spint);
		if(readvec(out, M_SPINT, header.id, sizeof(spint), len, fp)){
			free(out);
			return NULL;
		}
	}
	return out;
}
/**
   Read nmemb elements of type magic from file. Complex numbers are read as
   pairs of real numbers so that the byte swapping for fits files is done on
   each real component.

   @param[out] p The output array.
   @param[in] magic Magic number (data type) of the data in the file and p.
   @param[in] size Size in bytes of each element.
   @param[in] nmemb Number of elements to read.
   @param[in,out] fp The file to read from.
   @return 0 if succeeded, nonzero if failed.
*/
static int readvec_do(void* p, uint32_t magic, size_t size, size_t nmemb, file_t* fp){
	if(magic_iscomplex(magic)){
		return zfread(p, size/2, nmemb*2, fp);
	}else{
		return zfread(p, size, nmemb, fp);
	}
}
/**
   Read a vector from file and perform conversion if type mismatches.

   @param[out] p The output array.
   @param[in] magic_p Magic number (data type) of p.
   @param[in] magic_file Magic number (data type) of data in the file.
   @param[in] size Size in bytes of each element of p. Only used if no conversion is needed.
   @param[in] nmemb Number of elements to read.
   @param[in,out] fp The file to read from.
   @return 0 if succeeded, nonzero if failed.
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
		if(!(ans=readvec_do(p2, M2, sizeof(T2), nmemb, fp))){\
		for(size_t i=0; i<nmemb; i++){		\
			((T1*)p)[i]=(T1)(p2[i]);		\
		}}							        \
		free(p2);					        \
    }

	if(magic_p==magic_file){
		ans=readvec_do(p, magic_file, size, nmemb, fp);
	}
		TEST_CONVERT(M_DBL, double, M_FLT, float)
		//TEST_CONVERT(M_DBL, double, M_CMP, dcomplex)
		//TEST_CONVERT(M_DBL, double, M_ZMP, fcomplex)
		TEST_CONVERT(M_DBL, double, M_INT32, int)
		TEST_CONVERT(M_DBL, double, M_INT64, int64_t)
		TEST_CONVERT(M_DBL, double, M_SPINT, spint)
		TEST_CONVERT(M_FLT, float, M_DBL, double)
		//TEST_CONVERT(M_FLT, float, M_CMP, dcomplex)
		//TEST_CONVERT(M_FLT, float, M_ZMP, fcomplex)
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
		//TEST_CONVERT(M_INT32, int, M_CMP, dcomplex)
		//TEST_CONVERT(M_INT32, int, M_ZMP, fcomplex)
		TEST_CONVERT(M_INT32, int, M_INT64, int64_t)
		TEST_CONVERT(M_INT32, int, M_SPINT, spint)
		TEST_CONVERT(M_INT64, int64_t, M_DBL, double)
		TEST_CONVERT(M_INT64, int64_t, M_FLT, float)
		//TEST_CONVERT(M_INT64, int64_t, M_CMP, dcomplex)
		//TEST_CONVERT(M_INT64, int64_t, M_ZMP, fcomplex)
		TEST_CONVERT(M_INT64, int64_t, M_INT32, int)
		TEST_CONVERT(M_INT64, int64_t, M_SPINT, spint)
		TEST_CONVERT(M_SPINT, spint, M_DBL, double)
		TEST_CONVERT(M_SPINT, spint, M_FLT, float)
		//TEST_CONVERT(M_SPINT, spint, M_CMP, dcomplex)
		//TEST_CONVERT(M_SPINT, spint, M_ZMP, fcomplex)
		TEST_CONVERT(M_SPINT, spint, M_INT32, int)
		TEST_CONVERT(M_SPINT, spint, M_INT64, int64_t)
	else{
		ans=-1;
		warning("Please implement conversion from %x to %x\n", magic_file, magic_p);
	}
	return ans;
}
/**
   Force MKL to treat CPU as Intel
   https://danieldk.eu/Posts/2020-08-31-MKL-Zen.html
   For dsvd_inv of 2714x2714 array, it reduced the timing from 2 to 1.4 s

   @return Always 1.
*/
int mkl_serv_intel_cpu_true(){
	//dbg("mkl_serv_intel_cpu_true is called\n");
	return 1;
}

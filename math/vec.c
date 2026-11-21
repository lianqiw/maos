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
#include "mathdef.h"
#include "defs.h"/*Defines T, X, etc */

#define check_vec(p,N) (p?(N?1:0):(N?(error("p is null but N is not 0\n"),0):0))

/**
   sum of a vector using pointer and length

   @param[in] p Pointer to the vector.
   @param[in] nx Number of elements.
   @return Sum of the elements.
*/
T X(vecsum)(const T *restrict p, long nx){
    TD sum=0;
    if(check_vec(p, nx)){
        for(long ix=0; ix<nx; ix++){
            sum+=p[ix];
        }
    }
    return (T)sum;
}
/**
   normalize vector so that the sum is norm.

   @param[in,out] p Pointer to the vector, scaled in place.
   @param[in] nx Number of elements.
   @param[in] norm Target sum.
*/
void X(vecnormalize_sum)(T *restrict p, long nx, T norm){
    if(!check_vec(p, nx)) return;
    T ss=norm/X(vecdot)(p,NULL,NULL,nx);
    for(int i=0; i<nx; i++){
		p[i]*=ss;
    }
}
/**
   normalize vector so that the sum of abs is norm.

   @param[in,out] p Pointer to the vector, scaled in place.
   @param[in] nx Number of elements.
   @param[in] norm Target sum of absolute values.
*/
void X(vecnormalize_sumabs)(T *restrict p, long nx, T norm){
    if(!check_vec(p, nx)) return;
    TD ss=0;
    for (long i=0; i<nx; i++){
		ss+=ABS(p[i]);
    }
    ss=norm/(T)ss;
    for(int i=0; i<nx; i++){
		p[i]*=ss;
    }
}
/**
   Compute max and min. NAN values are skipped except for the first element.
   Complex values are converted into magnitude during comparison.

   @param[in] p Pointer to the vector.
   @param[in] N Number of elements.
   @param[out] max Maximum value. May be NULL.
   @param[out] min Minimum value. May be NULL.
*/
void X(vecmaxmin)(const T *restrict p, long N, R *max, R *min){
    if(!check_vec(p, N)) return; 
    R a,b;
    long i;
#ifdef COMP_COMPLEX
    a=b=ABS(p[0]);
#else
    a=b=p[0];
#endif    
    for(i=1; i<N; i++){
#ifdef COMP_COMPLEX
	    R tmp=ABS(p[i]);
#else
	    R tmp=p[i];
#endif
#ifndef COMP_LONG	
	    if(!isnan(tmp))
#endif
	    {
	        if(tmp>a) a=tmp;
	        if(tmp<b) b=tmp;
	    }
    }
    if(max)*max=a; 
    if(min)*min=b; 
}
/**
   compute sum(p1.*p2.*p3). A NULL vector is treated as all ones.

   @param[in] p1 First vector. May be NULL.
   @param[in] p2 Second vector. May be NULL.
   @param[in] p3 Third (real) vector. May be NULL.
   @param[in] n Number of elements.
   @return The weighted dot product.
*/
T X(vecdot)(const T *restrict p1, const T *restrict p2, const R *restrict p3, long n){
    TD ans=0;
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
	ans=(T)n;/*assume all ones. */
    }
    return  (T)ans;
}

/**
   find the maximum of abs of all elements using pointer and count.

   @param[in] p Pointer to the vector.
   @param[in] n Number of elements.
   @return Maximum absolute value. 0 if empty.
*/
R X(vecmaxabs)(const T*restrict p, long n){
    if(!check_vec(p, n)) return 0;
    R max,min;
    X(vecmaxmin)(p, n, &max, &min);
    max=FABS(max);
    min=FABS(min);
    return max>min?max:min;
}

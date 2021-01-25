/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
   normalize vector to sum to norm;*/
void X(normalize_sum)(T *restrict p, long nx, T norm){
    if(!check_vec(p, nx)) return;
    T ss=norm/X(vecdot)(p,NULL,NULL,nx);
    for(int i=0; i<nx; i++){
		p[i]*=ss;
    }
}
/**
   normalize vector to sum of abs to norm;*/
void X(normalize_sumabs)(T *restrict p, long nx, T norm){
    if(!check_vec(p, nx)) return;
    TD ss=0;
    for (long i=0; i<nx; i++){
		ss+=fabs(p[i]);
    }
    ss=norm/ss;
    for(int i=0; i<nx; i++){
		p[i]*=ss;
    }
}
/**
   normalize vector to max to max;*/
void X(normalize_max)(T *restrict p, long nx, T max){
    if(!check_vec(p,nx)) return;
    T ss=max/X(vecmaxabs)(p, nx);
    for(int i=0; i<nx; i++){
		p[i]*=ss;
    }
}
/**
   Compute max, min and sum. Has to handle NAN nicely. Complex values are
   converted into magnitude during comparison. */
void X(maxmin)(const T *restrict p, long N, R *max, R *min){
    if(!check_vec(p, N)) return; 
    R a,b;
    long i;
#ifdef COMP_LONG
    a=-INT_MAX;
    b=INT_MAX;
#else
    a=-INFINITY;
    b=INFINITY;
#endif
    for(i=0; i<N; i++){
#ifdef COMP_COMPLEX
	R tmp=fabs(p[i]);
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
   compute sum(p1.*p2.*p3)*/

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
 */
R X(vecmaxabs)(const T*restrict p, long n){
    if(!check_vec(p, n)) return 0;
    R max,min;
    X(maxmin)(p, n, &max, &min);
    max=fabs(max);
    min=fabs(min);
    return max>min?max:min;
}

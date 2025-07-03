/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#include "comb.h"
/**
   Computes the number of possibilities of selection k items from n items \f$C_n^k\f$:
   factorial(n)/(factorial(k)*factorial(n-k));
*/
static long comb_tot(long n, long k){
	return (long)round(factorial(n-k+1, n)/factorial(1, k));
}
/**
   initialize an initial combination composed a vector of non-negative numbers 0,1,2,...
*/
static int* comb_init(long k){
	int* comb=mycalloc(k, int);
	for(int i=0; i<k; i++){
		comb[i]=i;
	}
	return comb;
}
/**
   Find the next combination of selecting k from n. 
 */
static int comb_next(int* comb, long n, long k){
	if(n<1||k<1){
		return 0;
	}
	int i=k-1;
	comb[i]++;/*increment to next */
	while(comb[i]+k>=n+i+1&&i>0){/*out of range, increment previous one */
		i--;
		comb[i]++;
	}
	if(comb[0]+k>n){
		return 0;/*no more */
	}
	for(i=i+1; i<k; i++){
		comb[i]=comb[i-1]+1;
	}
	return 1;
}
/**
 * @brief Determine all valid combinations of stars for all wfs. 
	Assume that the stars are sorted in descending S/N order. 
	valid is set to -1 for stars that need to be skipped.
	flag=0: keep order. nstar should be equal to nwfstot. do not reorder stars. Only 1 combination.
	flag=1: regular.
	flag=2: make wfs[0] brighter than wfs[1] (stars are sorted from bright to dim)
	only implemented for npowfs=1 or npowfs=2
	@param nwfsmax 		npowfs*1 array for maximum number of wfs
	@param starvalid 	nstar*npowfs array for validity of star at powfs mode. -1 means invalid
	@param flag			see above

*/
lmat* comb_stars(lmat* nwfsmax, lmat *starvalid, int flag){
	lmat *res=NULL;
	const int npowfs=NX(nwfsmax);
	const int nstar=NX(starvalid);
	if(NY(starvalid)<npowfs){
		error("starvalid should have at least %d columns\n", npowfs);
		return NULL;
	}
	if(npowfs<1)return NULL;
	if(npowfs>2){
		error("To implement\n");
		return NULL;
	}
	int nleft=nstar;
	int nwfstot=0;
	int ncomb=1;
	int nwfs[npowfs];//actual nwfs for each powfs
	for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
		if(nleft>=P(nwfsmax,ipowfs)){
			nwfs[ipowfs]=P(nwfsmax,ipowfs);
		} else{
			nwfs[ipowfs]=nleft;
		}
		nwfstot+=nwfs[ipowfs];
		if(flag!=0){
			ncomb*=comb_tot(nleft, nwfs[ipowfs]);
		}
		nleft-=nwfs[ipowfs];
	}
	if(flag==2){
		ncomb/=comb_tot(nwfstot, nwfs[0]);
	}
retry:
	res=lnew(nwfstot, ncomb);
	int *comb=comb_init(nwfstot);//select all stars
	int count=0;
	do{
		if(count>=NY(res)){//this may happen if there are invalid stars and a retry
			dbg("count should be smaller than %ld.\n", NY(res));
			break;
		}
		if(npowfs==1){
			int ipowfs=0;
			int skip=0;
			for(int iwfs=0; iwfs<nwfs[ipowfs]; iwfs++){
				int istar=P(res, iwfs, count)=comb[iwfs];
				if(P(starvalid, istar, ipowfs)==-1){
					skip=1;
				}
			}
			if(!skip) count++;
		}else if(npowfs==2){
			int mask[nwfstot];//length is the same as comb[]
			int *comb2=comb_init(nwfs[0]);//select stars for powfs0
			do{
				if(count==NY(res)){//this may happen if there are invalid stars and a retry
					dbg("count should be smaller than %ld.\n", NY(res));
					break;
				}
				int skip=0;
				memset(mask, 0, sizeof(int)*nwfstot);
				int ipowfs=0;
				for(int iwfs=0; iwfs<nwfs[ipowfs]; iwfs++){
					int istar=P(res, iwfs, count)=comb[comb2[iwfs]];
					if(P(starvalid, istar, ipowfs)==-1){ 
						skip=1;
					}
					mask[comb2[iwfs]]=1;//mark already used
				}
				int jstar=0;
				ipowfs=1;
				for(int iwfs=0; iwfs<nwfs[ipowfs]; iwfs++){
					while(mask[jstar]) jstar++;
					int istar=P(res, iwfs+nwfs[0], count)=comb[jstar];
					if(P(starvalid, istar, ipowfs)==-1){
						skip=1;
					}
					mask[jstar]=1;
				}
				if(!skip) count++;
			}while(!(count&&flag!=1)&&comb_next(comb2, nwfstot, nwfs[0]));
			free(comb2);
		}else{
			error("Please implement for npowfs=%d\n", npowfs);
		}
	}while(!(count&&flag==0) && comb_next(comb, nstar, nwfstot));
	free(comb);
	if(count==0){
		lfree(res);
		if(nwfstot>1){//reduce number of wfs and retry
			nwfstot--;
			if(nwfs[1]>0){
				nwfs[1]--;
			}else if(nwfs[0]>0){
				nwfs[0]--;
			}
			goto retry;
		}
	}else if(count<NY(res)){
		lresize(res, NX(res), count);
	}else if(count>NY(res)){
		error("Memory overflows\n");//should never happen due to test above
	}
	return res;
}

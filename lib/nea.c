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
#include "nea.h"
/**
   Convert column vector nsax2 or nsax3 to sparse matrix. The third column of nea is coupling.

   @param[in] nea NEA (or its inverse or Cholesky factor) of size nsa x ng or nsa x 3.
   @param[in] ll Whether to fill the lower left (yx) coupling term.
   @param[in] ur Whether to fill the upper right (xy) coupling term.
   @param[in] ng Number of gradients per subaperture. If ng>2, treat as diagonal matrix (for raw PWFS).
   @return The (nsa*ng) x (nsa*ng) sparse matrix. Newly allocated; the caller frees it.
*/
dsp* nea2sp(dmat* nea, int ll, int ur, int ng){
	if(NY(nea) != ng && NY(nea)!=3){
		error("nea has wrong format:%ldx%ld. Expect 2 or 3 columns.\n", NX(nea), NY(nea));
	}
	if(NY(nea)==ng){//off-diagonal not available
		ll=0;
		ur=0;
	}
	const long nsa=NX(nea);
	dsp* sanea=dspnew(nsa*ng, nsa*ng, (ng+(ll?1:0)+(ur?1:0))*nsa);
	spint* pp=sanea->pp;
	spint* pi=sanea->pi;
	real* px=sanea->px;
	long count=0;
	for(long ig=0; ig<ng; ig++){
		for(long isa=0; isa<nsa; isa++){
			pp[isa+nsa*ig]=count;
			if(ig==1&&ur&&P(nea, isa, 2)){/*Cxy */
				pi[count]=isa;
				px[count]=P(nea, isa, 2);
				count++;
			}
			{/*Cxx */
				pi[count]=isa+nsa*ig;
				px[count]=P(nea, isa, ig);
				count++;
			}
			if(ig==0&&ll&&P(nea, isa, 2)){/*Cyx */
				pi[count]=isa+nsa;
				px[count]=P(nea, isa, 2);
				count++;
			}
		}
	}
	pp[nsa*ng]=count;
	if(count>sanea->nzmax){
		error("memory overflow\n");
	}else if(count<sanea->nzmax){
		dspsetnzmax(sanea, count);
	}
	return sanea;
}
/**
   Check the dimension of NEA and fix if possible

   @param[in,out] nea The NEA. A single column is reshaped to nsa columns. Errors out if NULL or the format is wrong.
   @param[in] nsa Number of subapertures.
   @param[in] ng Number of gradients per subaperture.
*/
void nea_check(dmat *nea, long nsa, int ng){
	if(!nea){
		error("nea is not defined.\n");
	} else if(NY(nea)==1){
		reshape(nea, nsa, NX(nea)/nsa);
	}
	if(NX(nea)!=nsa||(NY(nea)!=ng&&NY(nea)!=3)){
		error("nea has wrong format (%ldx%ld), should be (%ldx%d).\n", NX(nea), NY(nea), nsa, 3);
	}
}
/**
   Common routine for nea_chol, nea_mm, and nea_inv to prepare the output.

   @param[in,out] pout The output. Allocated if NULL; resized if it needs the cross term column.
   @param[in] in The input, nsa x ng or nsa x 3.
   @param[in] ng Number of gradients per subaperture.
   @return 1 if the cross term (third column) is present and ng is 2, 0 otherwise.
*/
static int nea_prep(dmat **pout, const dmat *in, const int ng){
	if(NY(in)!=ng&&NY(in)!=3){
		error("in has wrong format: %ldx%ld\n", NX(in), NY(in));
	}
	int isxy=(NY(in)>ng&&ng==2)?1:0;//cross term
	int ncol=ng==2?3:ng;
	if(!*pout){
		*pout=dnew(NX(in), ncol);
	}else if(isxy&&NY((*pout))!=ncol){
		dresize(*pout, NX(in), ncol);
	}
	return isxy;
}
/**
   Apply cholesky in a 2x2 symmetric matrix packed in [a[0,0],a[1,1],a[0,1]] as row vector.
   Input and output may be the same.

   @param[in,out] pout The lower triangular Cholesky factor in the same packing. Allocated if NULL.
   @param[in] in The symmetric matrix of each subaperture, nsa x 2 or nsa x 3.
   @param[in] ng Number of gradients per subaperture. If ng>2, treat as diagonal matrix.
*/
void nea_chol(dmat **pout, const dmat *in, const int ng){
	const int isxy=nea_prep(pout, in, ng);
	dmat *out=*pout;
	if(ng==2){
		for(int isa=0; isa<NX(in); isa++){
		//Use temporary variable to handle the case that out and in is the same.
			real a=sqrt(P(in, isa, 0));
			real b=(isxy && a)?(P(in, isa, 2)/a):0;
			real c=sqrt(P(in, isa, 1)-b*b);
			P(out, isa, 0)=a;
			P(out, isa, 1)=c;
			if(isxy) P(out, isa, 2)=b;
		}
	} else{
		for(int ig=0; ig<NX(in)*ng; ig++){
			P(out, ig)=sqrt(P(in, ig));
		}
	}
}
/**
   Apply LL' to lower diagonal matrix packed in [a[0,0],a[1,1],a[0,1]] as row vector.
   Input and output may be the same.

   @param[in,out] pout The product LL' in the same packing. Allocated if NULL.
   @param[in] in The lower triangular matrix L of each subaperture, nsa x 2 or nsa x 3.
   @param[in] ng Number of gradients per subaperture. If ng>2, treat as diagonal matrix.
*/
void nea_mm(dmat **pout, const dmat *in, const int ng){
	const int isxy=nea_prep(pout, in, ng);
	dmat *out=*pout;
	if(ng==2){
		for(int isa=0; isa<NX(in); isa++){
		//Use temporary variable to handle the case that out and in is the same.
			real a=P(in, isa, 0);
			real b=isxy?(P(in, isa, 2)):0;
			real c=P(in, isa, 1);
			P(out, isa, 0)=a*a;
			P(out, isa, 1)=c*c+b*b;
			if(isxy) P(out, isa, 2)=a*b;
		}
	} else{
		for(int ig=0; ig<NX(in)*ng; ig++){
			P(out, ig)=P(in, ig)*P(in, ig);
		}
	}
}
/**
   Apply matrix inversion in a 2x2 symmetric matrix packed in row vector [a[0,0],a[1,1],a[0,1]].
   Input and output may be the same. Entries that are singular are not updated.

   @param[in,out] pout The scaled inverse in the same packing. Allocated if NULL.
   @param[in] in The symmetric matrix of each subaperture, nsa x 2 or nsa x 3.
   @param[in] ng Number of gradients per subaperture. If ng>2, treat as diagonal matrix.
   @param[in] scale Scaling factor applied to the inverse.
*/
void nea_inv(dmat **pout, const dmat *in, int ng, real scale){
	const int isxy=nea_prep(pout, in, ng);
	dmat *out=*pout;
	if(ng==2){
		for(int isa=0; isa<NX(in); isa++){
			//use double to prevent overflow of intermediate value
			double xx=P(in, isa, 0);
			double yy=P(in, isa, 1);
			double xy=isxy?P(in, isa, 2):0;
			double invdet=scale/(xx*yy-xy*xy);
			if(isfinite(invdet)){
				P(out, isa, 0)=invdet*yy;
				P(out, isa, 1)=invdet*xx;
				if(xy) P(out, isa, 2)=-invdet*xy;
			}
		}
	} else{
		for(int ig=0; ig<NX(in)*ng; ig++){
			if(P(in, ig)){
				P(out, ig)=scale/P(in, ig);
			}
		}
	}
}

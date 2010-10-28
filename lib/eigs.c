/*
  Copyright 2009,2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#include "eigs.h"
/**
   \file eigs.c
   Compute a few eigenvalues of a sparse matrix.
   Wrapper of the arpack library to provide convenient interface.   
   eigtype:
   "LA" - compute the NEV largest (algebraic) eigenvalues.
   "SA" - compute the NEV smallest (algebraic) eigenvalues.
   "LM" - compute the NEV largest (in magnitude) eigenvalues.
   "SM" - compute the NEV smallest (in magnitude) eigenvalues. 
   "BE" - compute NEV eigenvalues, half from each end of the

  the returned eig is in ascending order.
*/
/**
   compute the eigenvallues. Full interface, use other want interfaces for
   simple calling. This is a wrapper to the Arpack library.  
   
   \bug The Arpack does not work in 32 bit machines.*/
void eigs(dmat **eig,      /**<[out] computed eigen values*/
	  dmat **eigv,     /**<[out] computed eigen vectors. NULL to disable*/
	  const void *A,   /**<[in] the matrix, dsp or dspcell*/
	  ATYPE Atype,     /**<[in] type of A, A_SPARSE or A_SPCELL*/
	  const int neig,  /**<[in] number of eigen values wanted*/
	  char *eigtype   /**<[in] type of eigen values. See the file description.*/
	  ){
    if(!(*eig)){
	*eig=dnew(neig,1);
    }
    
    //make sure eig is neig long.
    int maxiter=3000;

    //input to dsaupd
    int ido=0;
    char* bmat="I";
    int N;
    if(Atype==A_SPARSE){
	N=((dsp*)A)->m;//dimension of the eigen problem.
    }else if(Atype==A_SPCELL){
	N=0;
	for(int ix=0; ix<((spcell*)A)->nx; ix++){
	    N+=((spcell*)A)->p[ix]->m;
	}
    }else{
	error("Invalid Atype");
    }

    char* which=eigtype;
    int nev=neig;//# of eig val
    double tol=0;
    double *resid=malloc(sizeof(double)*N);
    int ncv=20;//# if columns used.
    if(ncv<neig){
	if(neig>N) error("Invalid neig\n");
	warning("Too many eigen values requested\n");
	ncv=neig+20;
	if(ncv>N) ncv=N;
    }
    double *V=malloc(sizeof(double)*N*ncv);
    int ldv=N;
    int iparam[11]; 
    iparam[0]=1;
    iparam[2]=maxiter;
    iparam[6]=1;
    int ipntr[11];//output
    double *workd=malloc(sizeof(double)*N*3);
    int lworkl=ncv*(ncv+8);
    double *workl=malloc(sizeof(double)*lworkl);
    int info=0;

    //input to dseupd
    int rvec=(eigv!=NULL);//whether want eigen vectors.
    char *howmany="A";
    int select0[nev];
    double *D=(*eig)->p;//output ritz values.(eig val)
    //Z,lvz are the save as above V, ldv
    double sigma;//not really used.
    int ierr=0;
    do{
	dsaupd_(&ido,bmat,&N,which,&nev,&tol,resid,&ncv,V,&ldv,
		iparam,ipntr,workd,workl,&lworkl,&info);
	//set to 0 because spmulvec accumulates.
	memset(workd+ipntr[1]-1,0,sizeof(double)*N);
	if(Atype==A_SPARSE){
	    spmulvec(workd+ipntr[1]-1,(const dsp*)A,workd+ipntr[0]-1,1);
	}else if(Atype==A_SPCELL){
	    spcellmulvec(workd+ipntr[1]-1,(const spcell*)A,workd+ipntr[0]-1,1);
	}else{
	    error("Invalid A type\n");
	}
    }while(ido == -1 || ido == 1);
    if(info<0){
	error("Error with saupd with info is %d\n",info);
    }else{
	switch (info){
	case 1:
	    warning("Maximum number oof iterations reached.\n"); break;
	case 3:
	    warning("No shifted could be applied during implicit Arnoldi update"
		    ", please try to increase ncv"); break;
	}
	
	dseupd_(&rvec,howmany,select0,D,V,&ldv,&sigma,
		bmat,&N,which,&nev,&tol,resid,&ncv,V,&ldv,
		iparam,ipntr,workd,workl,&lworkl,&ierr);
	if(ierr!=0){
	    error("Error with seupd, info is %d\n",ierr);
	}else{
	    int nconv=iparam[4];//# is converged eigs
	    //info("nconv=%d\n",nconv);
	    if(nconv<nev){
		error("Want %d eigs, only %d converged\n",nev,nconv);
	    }
	    if(eigv){
		*eigv=dnew(N,neig);
		memcpy((*eigv)->p,V,sizeof(double)*N*nev);
	    }
	}
    }
    free(workl);
    free(workd);
    free(V);
    free(resid);
}
/**
   compute the maximum eigenvalue.
*/
double spmaxeig(const dsp *A){
    dmat *eig=NULL;
    eigs(&eig,NULL,(const void*)A,A_SPARSE,1,"LM");
    double eig0=eig->p[0];
    dfree(eig);
    return eig0;
}
/**
   compute the minimum eigenvalue.
*/
double spmineig(const dsp *A){
    dmat *eig=NULL;
    eigs(&eig,NULL,(const void*)A,A_SPARSE,1,"SM");
    double eig0=eig->p[0];
    dfree(eig);
    return eig0;
}
/**
   compute the maximum eigenvalue.*/
double spcellmaxeig(const spcell *A){
    dmat *eig=NULL;
    eigs(&eig,NULL,(const void*)A,A_SPCELL,1,"LM");
    double eig0=eig->p[0];
    dfree(eig);
    return eig0;
}
/**
   compute the minimum eigenvalue.*/
double spcellmineig(const spcell *A){
    dmat *eig=NULL;
    eigs(&eig,NULL,(const void*)A,A_SPCELL,1,"SM");
    double eig0=eig->p[0];
    dfree(eig);
    return eig0;
}
/**
   Apply tikhonov regulation with threshold of thres*max(eig(A)).  The
   maximum eigen value of A is on the same order of magnitude as the maximum
   number of abs(A).*/
void sptikcr(dsp *A, double thres){
    double eig=spmaxeig(A);
    dsp *Reg=spnewdiag(A->m,NULL,thres*eig);
    spadd(&A,Reg);
    spfree(Reg);
}

/**
   Apply tikhonov regulation with threshold of thres*max(eig(A)).  The
   maximum eigen value of A is on the same order of magnitude as the maximum
   number of abs(A).*/
void spcelltikcr(spcell *A, double thres){
    double eig=spcellmaxeig(A);
    info("max eigenvalue is %g\n",eig);
    for(int ii=0; ii<A->ny; ii++){
	//Add to diagonals only.
	dsp *iA=A->p[ii+ii*A->nx];
	if(iA->m!=iA->n)
	    error("diagonal block is not diagonal matrix\n");
	dsp *Reg=spnewdiag(iA->m,NULL,thres*eig);
	spadd(&iA,Reg);
	spfree(Reg);
    }
}

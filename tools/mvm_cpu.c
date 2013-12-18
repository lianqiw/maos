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

/**
   \file mvm_cpu.c
   Testing mvm on CPU
*/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 
#endif
#include <errno.h>
#include "../lib/aos.h"

static int use_trans=1;
static void mvmt_do(const float *mvmt, float *g, float *dm, int nact, int ngleft, int ngtot){
 
#pragma omp parallel for
    for(int ia=0; ia<nact; ia++){
	register float tmp=dm[ia];
#ifdef __ICC
#pragma unroll
	__assume_aligned(mvmt,128);
	__assume_aligned(g,128);
	__assume_aligned(dm,128);
#endif
#ifdef __ICC
#pragma vector aligned
#pragma ivdep
#pragma simd vectorlength(16) assert
#endif
	for(int ig=0; ig<ngleft; ig++){
	    tmp+=mvmt[ig+ia*ngtot]*g[ig];
	}
	dm[ia]=tmp;
    }
}
static void mvm_do(const float *mvm, float *g, float *dm, int nact, int ngleft){

    for(int ig=0; ig<ngleft; ig++){
#pragma omp parallel for
	for(int ia=0; ia<nact; ia++){
	    dm[ia]+=mvm[ia+ig*nact]*g[ig];
	}
    }
}
static void mtch_do(const float *mtch, const short *pix, const short *pixbias, 
		    float *grad, float *im0, float *imx, float *imy,
		    const int *saind, int nsa, float ct, float st){
#pragma omp parallel for
    for(int isa=0; isa<nsa; isa++){
	const int npix=saind[isa+1]-saind[isa];
	const float *mtchx=mtch+saind[isa]*2;
	const float *mtchy=mtchx+npix;
	const short *pixi=pix+saind[isa];
	const short *bias=pixbias+saind[isa];
	float *im0i=im0+saind[isa];
	float *imxi=imx+saind[isa];
	float *imyi=imy+saind[isa];
	float gx=0,gy=0;
	for(int ipix=0; ipix<npix; ipix++){
	    short et=(pixi[ipix]-bias[ipix]);
	    gx+=mtchx[ipix]*et;
	    gy+=mtchy[ipix]*et;
	    im0i[ipix]+=et;
	    imxi[ipix]+=et*ct;
	    imyi[ipix]+=et*st;
	}
	grad[isa*2]=gx;
	grad[isa*2+1]=gy;
    }	
}
int main(int argc, char *argv[]){
    enum{
	P_EXE,
	P_FRAC,
	P_NSTEP,
	P_TOT,
    };
    if(argc!=P_TOT){
	info2("Usage: \n\tenv MVM_CLIENT=hostname MVM_PORT=port MVM_SASTEP=sastep ./mvm_cpu fraction nstep\n");
	_Exit(0);
    }
    int fraction=strtol(argv[P_FRAC], NULL, 10);
    int nstep=strtol(argv[P_NSTEP], NULL, 10);
    int nstep0=nstep>1?20:0;//warm up
    dmat *d_saind=dread("NFIRAOS_saind");
    const int nsa=(d_saind->nx-1)/fraction;
    int *saind=(int*)malloc(sizeof(int)*(1+nsa));
    for(int i=0; i<nsa+1; i++){
	saind[i]=(int)d_saind->p[i];
    }
    dfree(d_saind);
    const int totpix=saind[nsa];
    const int nact=6981;//active subapertures.
    int ng=nsa*2;
    float FSMdelta=-0.2;
    smat *dm=snew(nact,1);
    smat *mvm=snew(nact, ng);
    smat *mtch=snew(totpix*2,1);
    smat *grad=snew(ng,1);
    smat *im0=snew(totpix,3);
    short *pix=malloc(totpix*sizeof(short));
    short *pixbias=malloc(totpix*sizeof(short));
    {
	rand_t rseed;
	seed_rand(&rseed, 1);
	srandu(mvm, 1e-7, &rseed);
	srandu(mtch, 1, &rseed);
	for(int i=0; i<totpix; i++){
	    pix[i]=(short)(randu(&rseed)*25565);
	    pixbias[i]=(short)(randu(&rseed)*1000);
	}
    }
    smat *mvmt=strans(mvm);
    int sastep=200;//how many subapertures each time
    int nrep=1;
    if(getenv("MVM_NREP")){
	nrep=strtol(getenv("MVM_NREP"), NULL, 10);
    }
    if(getenv("MVM_SECT")){
	sastep=nsa/strtol(getenv("MVM_SECT"), NULL, 10);
    }
    if(getenv("MVM_TRANS")){
	use_trans=strtol(getenv("MVM_TRANS"), NULL, 10);
    }
    if(getenv("MVM_SASTEP")){
	sastep=strtol(getenv("MVM_SASTEP"), NULL, 10);
    }
    info2("use_trans=%d, nrep=%d, sastep=%d\n", use_trans, nrep, sastep);
    int sock=-1;
    char* MVM_CLIENT=getenv("MVM_CLIENT");
    if(MVM_CLIENT){
	short port=(short)strtol(getenv("MVM_PORT"), NULL, 10);
	sock=connect_port(MVM_CLIENT, port, 0 ,1);
	if(sock!=-1) {
	    info2("Connected\n");
	    int cmd[7];
	    cmd[0]=nact;
	    cmd[1]=nsa;
	    cmd[2]=sastep;
	    cmd[3]=totpix;
	    cmd[4]=nstep;
	    cmd[5]=nstep0;
	    cmd[6]=2;
	    if(stwriteintarr(sock, cmd, 7) 
	       || stwriteintarr(sock, saind, nsa+1)
	       || stwrite(sock, pix, sizeof(short)*totpix)){
		close(sock); sock=-1;
		warning("Failed: %s\n", strerror(errno));
	    }
	}
    }
    int ready=0;
    if(sock!=-1 && stwriteint(sock, ready)){
	warning("error send ready signal: %s\n", strerror(errno));
	close(sock); sock=-1;
    }
    smat *timing=snew(nstep, 1);
    TIC;
    float timtot=0, timmax=0, timmin=INFINITY;
    set_realtime(-1, -20);
    for(int jstep=-nstep0; jstep<nstep; jstep++){
	int istep=jstep<0?0:jstep;
	tic;
	double theta=M_PI*0.5*istep+FSMdelta;
	float cd=cos(theta);
	float sd=cos(theta);
	szero(dm);
	for(int isa=0; isa<nsa; isa+=sastep){
	    int npixleft;
	    int nsaleft;
	    if(nsa<isa+sastep){//terminate
		npixleft=totpix-saind[isa];
		nsaleft=nsa-isa;
	    }else{
		npixleft=saind[isa+sastep]-saind[isa];
		nsaleft=sastep;
	    }

	    short *pcur=pix+saind[isa];
	    if(sock!=-1){
		if(stread(sock, pcur, sizeof(short)*npixleft)){
		    warning("failed: %s\n", strerror(errno));
		    close(sock); sock=-1;
		    _Exit(1);
		}
		if(isa==0) tic;
	    }
	    //Matched filter
	    mtch_do(mtch->p, pix, pixbias, 
		    grad->p+isa*2, im0->p, im0->p+totpix, im0->p+totpix*2,
		    saind+isa, nsaleft, cd, sd);
	    //MVM
	    for(int irep=0; irep<nrep; irep++){
		if(use_trans){
		    mvmt_do(mvmt->p+isa*2, grad->p+isa*2,dm->p, nact, nsaleft*2, ng);
		}else{
		    mvm_do(mvm->p+isa*2*nact, grad->p+isa*2, dm->p, nact, nsaleft*2);
		}
	    }
	}//for isa
	if(sock!=-1){
	    if(stwrite(sock, dm->p, sizeof(float)*nact)){
		warning("error write dmres: %s\n", strerror(errno));
		close(sock); sock=-1;
		_Exit(1);
	    }
	    if(streadint(sock, &ready)){//acknowledgement.
		warning("error read ack failed: %s\n", strerror(errno));
		close(sock), sock=-1;
		_Exit(1);
	    }
	    timing->p[istep]=ready*1.e-6;
	}else{
	    timing->p[istep]=toc3;//do not tic.
	}
	if(jstep==istep){
	    timtot+=timing->p[istep];
	    if(timmax<timing->p[istep]){
		timmax=timing->p[istep];
	    }
	    if(timmin>timing->p[istep]){
		timmin=timing->p[istep];
	    }
	}
    }//for istep
    float timmean=timtot/nstep;
    info2("Timing is mean %.3f, max %.3f min %.3f. BW is %.1f of 51.2GB/s\n",
	  timmean*1e3, timmax*1e3, timmin*1e3, nrep*(nact*ng+nact+ng)*sizeof(float)/timmean/(1024*1024*1024));
    swrite(timing, "cpu_timing_%s", myhostname());
    if(nstep==1){
	do_write("cpu_pix", 1, sizeof(short), M_INT16, NULL, pix, totpix, 1);
	do_write("cpu_pixbias", 1, sizeof(short), M_INT16, NULL, pixbias, totpix, 1);
	swrite(dm, "cpu_dm");
	swrite(grad, "cpu_grad");
	swrite(mvm, "cpu_mvm");
	swrite(mtch, "cpu_mtch");
    }
   
}

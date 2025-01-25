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
#include "../lib/aos.h"
/**
  Test and find good FFT numbers.


 
int main(){
    int nlen=200;
    cmat *data[nlen];
    rand_t stat;
    seed_rand(&stat,1);
    int i;
    for(i=2; i<nlen; i++){
	data[i]=cnew(nlen,nlen);
	//cfft2plan(data[i],-1);
	//cfft2plan(data[i],1);
	crandn(data[i],20,&stat);
    }
    char fn[64];
    FILE *fp=fopen("fft.txt","w");
    real tk;
    for(i=2; i<nlen; i++){
	sprintf(fn,"%d",i);
	for(int j=0; j<10;j++){
	    cfft2(data[i],-1);
	}
	tk=myclockd();
	for(int j=0; j<100;j++){
	    cfft2(data[i],-1);
	    cfft2(data[i],1);
	    cscale(data[i],1./(j*j));
	}
	tk=myclockd()-tk;
	fprintf(stderr,"%d %g seconds\n",i, tk);
	fprintf(fp,"%d %g\n",i, tk);
    }
    fflush(fp);
    fclose(fp);
}
*/
/*
static int test_fft_speed_small(){
    int nis=128;
    int *is=mycalloc(nis,int);
    dmat *tim=dnew(nis,1);
    for(int ii=0; ii<nis; ii++){
	is[ii]=ii+1;
    }
    ccell *ac=cellnew(nis,1);
    rand_t stat;
    seed_rand(&stat,1);
    for(int ii=0; ii<nis; ii++){
	P(ac,ii)=cnew(is[ii],is[ii]);
	//cfft2plan(P(ac,ii),-1);
	crandn(P(ac,ii),20,&stat);
    }
    TIC;
    for(int ii=0; ii<nis; ii++){
	info("size %4d: ",is[ii]);
	tic;
	for(int i=0; i<1000; i++){
	    cfft2(P(ac,ii),-1);
	}
	toc("fft");
	P(tim,ii)=toc3;
    }
    writebin(tim,"fft_timing");
}

static void test_fft_speed(){
    int nis=2048;
    int *is=mycalloc(nis,int);
    dmat *tim=dnew(nis,1);
    for(int ii=0; ii<nis; ii++){
	is[ii]=(ii+1)*2;
    }
    ccell *ac=cellnew(nis,1);
    rand_t stat;
    seed_rand(&stat,1);
    TIC;
    for(int ii=0; ii<nis; ii++){
	info("size %4d: ",is[ii]);
	tic;
	P(ac,ii)=cnew(is[ii],is[ii]);
	//cfft2plan(P(ac,ii),-1);
	crandn(P(ac,ii),20,&stat);
	toc("plan");
    }
    toc("plan");
    for(int ii=0; ii<nis; ii++){
	info("size %4d: ",is[ii]);
	tic;
	int nrepeat;
	if(is[ii]<300)
	    nrepeat=100;
	else if(is[ii]<1000)
	    nrepeat=10;
	else
	    nrepeat=1;

	for(int i=0; i<nrepeat; i++){
	    cfft2(P(ac,ii),-1);
	}
	toc("fft");
	P(tim,ii)=toc3/nrepeat;
    }
    writebin(tim,"fft_timing");
    }*/
static void test_fft_special(){
    int nis=2;
    int *is=mycalloc(nis,int);
    dmat *tim=dnew(nis,1);
    is[0]=3824;
    is[1]=4096;
    ccell *ac=ccellnew(nis,1);
    rand_t rstat;
    seed_rand(&rstat,1);
    TIC;
    for(int ii=0; ii<nis; ii++){
	info("size %4d: ",is[ii]);
	tic;
	P(ac,ii)=cnew(is[ii],is[ii]);
	//cfft2plan(P(ac,ii),-1);
	//cfft2partialplan(P(ac,ii),512,-1);
	crandn(P(ac,ii),20,&rstat);
	toc("plan");
    }

    for(int ii=0; ii<nis; ii++){
	info("size %4d: ",is[ii]);
	tic;
	int nrepeat;
	if(is[ii]<300)
	    nrepeat=100;
	else if(is[ii]<1000)
	    nrepeat=10;
	else
	    nrepeat=4;

	for(int i=0; i<nrepeat; i++){
	    cfft2(P(ac,ii),-1);
	}
	toc("fft");
	for(int i=0; i<nrepeat; i++){
	    cfft2partial(P(ac,ii),512,-1);
	}
	toc("fft2partial");
	P(tim,ii)=toc3/nrepeat;
    }
    writebin(tim,"fft_timing");

}
int main(){
    test_fft_special();
}

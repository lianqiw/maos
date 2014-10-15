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
    double tk;
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
    int *is=calloc(nis,sizeof(int));
    dmat *tim=dnew(nis,1);
    for(int ii=0; ii<nis; ii++){
	is[ii]=ii+1;
    }
    ccell *ac=cellnew(nis,1);
    rand_t stat;
    seed_rand(&stat,1);
    for(int ii=0; ii<nis; ii++){
	ac->p[ii]=cnew(is[ii],is[ii]);
	//cfft2plan(ac->p[ii],-1);
	crandn(ac->p[ii],20,&stat);
    }
    TIC;
    for(int ii=0; ii<nis; ii++){
	info2("size %4d: ",is[ii]);
	tic;
	for(int i=0; i<1000; i++){
	    cfft2(ac->p[ii],-1);
	}
	toc2("fft");
	tim->p[ii]=toc3;
    }
    writebin(tim,"fft_timing");
}

static void test_fft_speed(){
    int nis=2048;
    int *is=calloc(nis,sizeof(int));
    dmat *tim=dnew(nis,1);
    for(int ii=0; ii<nis; ii++){
	is[ii]=(ii+1)*2;
    }
    ccell *ac=cellnew(nis,1);
    rand_t stat;
    seed_rand(&stat,1);
    TIC;
    for(int ii=0; ii<nis; ii++){
	info2("size %4d: ",is[ii]);
	tic;
	ac->p[ii]=cnew(is[ii],is[ii]);
	//cfft2plan(ac->p[ii],-1);
	crandn(ac->p[ii],20,&stat);
	toc("plan");
    }
    toc("plan");
    for(int ii=0; ii<nis; ii++){
	info2("size %4d: ",is[ii]);
	tic;
	int nrepeat;
	if(is[ii]<300)
	    nrepeat=100;
	else if(is[ii]<1000)
	    nrepeat=10;
	else
	    nrepeat=1;

	for(int i=0; i<nrepeat; i++){
	    cfft2(ac->p[ii],-1);
	}
	toc2("fft");
	tim->p[ii]=toc3/nrepeat;
    }
    writebin(tim,"fft_timing");
    }*/
static void test_fft_special(){
    int nis=2;
    int *is=calloc(nis,sizeof(int));
    dmat *tim=dnew(nis,1);
    is[0]=3824;
    is[1]=4096;
    ccell *ac=cellnew(nis,1);
    rand_t rstat;
    seed_rand(&rstat,1);
    TIC;
    for(int ii=0; ii<nis; ii++){
	info2("size %4d: ",is[ii]);
	tic;
	ac->p[ii]=cnew(is[ii],is[ii]);
	//cfft2plan(ac->p[ii],-1);
	//cfft2partialplan(ac->p[ii],512,-1);
	crandn(ac->p[ii],20,&rstat);
	toc("plan");
    }

    for(int ii=0; ii<nis; ii++){
	info2("size %4d: ",is[ii]);
	tic;
	int nrepeat;
	if(is[ii]<300)
	    nrepeat=100;
	else if(is[ii]<1000)
	    nrepeat=10;
	else
	    nrepeat=4;

	for(int i=0; i<nrepeat; i++){
	    cfft2(ac->p[ii],-1);
	}
	toc2("fft");
	for(int i=0; i<nrepeat; i++){
	    cfft2partial(ac->p[ii],512,-1);
	}
	toc2("fft2partial");
	tim->p[ii]=toc3/nrepeat;
    }
    writebin(tim,"fft_timing");

}
int main(){
    test_fft_special();
}

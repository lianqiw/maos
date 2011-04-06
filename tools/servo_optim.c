/**
   standalone routine that does servo filtering.
*/
#include "aos.h"
int main(int argc, char **argv){
    if(argc<6){
	info2("Usage: %s psd.bin dtrat dt sigma gainout.bin\n", argv[0]);
	exit(1);
    }
    dmat *psd=dread("%s",argv[1]);
    long dtrat=strtol(argv[2],NULL,10);
    double dt=strtod(argv[3],NULL);
    double sigma=strtod(argv[4],NULL);//m^2
    dmat *sigma2=dnew(1,1); sigma2->p[0]=sigma;
    dcell *gain=servo_typeII_optim(psd,dtrat,dt,sigma2);
    dwrite(gain->p[0],"%s",argv[5]);
    dfree(sigma2);
    dcellfree(gain);
    dfree(psd);
}


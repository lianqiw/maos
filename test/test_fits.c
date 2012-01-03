#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/mman.h>
#include "../lib/aos.h"
int main(){
    int N=4;
    dmat *A=dnew(N,N);
    dcell *B=dcellnew(2,2);
    B->p[0]=dnew(N,N);
    B->p[2]=dnew(N,N);
    rand_t rstat;
    seed_rand(&rstat, 1);
    drandn(A, 1, &rstat);
    drandn(B->p[0], 1, &rstat);
    drandn(B->p[2], 1, &rstat);
    A->header=strdup("dx=1/64\ndy=1/64\nLong=kaflasdjfl ksflksjlfjasdflkj asldfj als fals fl asfjalsdf jalsf djasldf jal fsalfkjasdflkj asldfkj asldf \nadf\nwvf=1.e-4");
    dwrite(A, "A.bin");
    dwrite(A, "A.fits");
    dwrite(A, "A.fits.gz");
    dcellwrite(B, "B.bin");
    dcellwrite(B, "B.fits");
    dmat *C=dread("A.fits");
    dwrite(C, "C.fits");
    dcell *D=dcellread("B.fits");
    dcellwrite(D, "D.fits");
    dcell *E=dcellread("A.fits");
    dcellwrite(E, "E.fits");
    dcellwrite(E, "E.fits.gz");
}

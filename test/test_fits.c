



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
    writebin(A, "A.bin");
    writebin(A, "A.fits");
    writebin(A, "A.fits.gz");
    writebin(B, "B.bin");
    writebin(B, "B.fits");
    dmat *C=dread("A.fits");
    writebin(C, "C.fits");
    dcell *D=dcellread("B.fits");
    writebin(D, "D.fits");
    dcell *E=dcellread("A.fits");
    writebin(E, "E.fits");
    writebin(E, "E.fits.gz");
}

/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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





#include <sys/mman.h>
#include "../lib/aos.h"
int main(){
    int N=4;
    dmat *A=dnew(N,N);
    dcell *B=dcellnew(2,2);
    P(B,0)=dnew(N,N);
    P(B,2)=dnew(N,N);
    rand_t rstat;
    seed_rand(&rstat, 1);
    drandn(A, 1, &rstat);
    drandn(P(B,0), 1, &rstat);
    drandn(P(B,2), 1, &rstat);
    A->keywords=strdup("dx=1/64\ndy=1/64\nLong=kaflasdjfl ksflksjlfjasdflkj asldfj als fals fl asfjalsdf jalsf djasldf jal fsalfkjasdflkj asldfkj asldf \nadf\nwvf=1.e-4");
    writebin(A, "A.bin");
    writebin(A, "A.fits");
    writebin(A, "Azip.fits.gz");
    writebin(B, "B.bin");
    writebin(B, "B.fits");
    dfree(A);
    dcellfree(B);
    dmat *C=dread("A.fits");
    dfree(C);
    C=dread("Azip");
    dfree(C);
    C=dread("A.bin");
    dfree(C);

    dcell *D=dcellread("B.fits");
    dcellfree(D);
    D=dcellread("B.bin");
    dcellfree(D);
    
    //read dmat as dcell
    D=dcellread("A.fits");
    dcellfree(D);
}

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

#include "../lib/aos.h"
/*TIC; */
/*
void benchmark()
    tic;
    dspcell *G0=dspcellread("G0.bin");
    toc("");tic;
    dspcell *G0T=dspcelltrans(G0);
    toc("");tic;
    dspcell *A=dspcellmulspcell(G0T,G0,1);//this is slow.
    toc("");
    dcell *x=dcellread("rhs.bin");
    tic;
    real res=dspcellwddot(x,A,x);//this takes 0.15s in debugging 0.04 in optim
    toc("");
    printf("%g\n",res);
    tic;
    dcell *y=NULL;
    dspcellmulmat(&y, A,x,1);
    real res2=dcellinn(x,y);//this takes 0.15s in debugging 0.04 in optim
    toc("");
    printf("%g\n",res2);
    //dcell *saneai=dcellread("saneai.bin");
    tic;
    dcell *z=NULL;
    dspcellmulmat(&z, G0,x,1);toc("");
    dcell *z2=NULL;
    dspcellmulmat(&z2, G0T,z,1);toc("");
    real res3=dcellinn(x,z2);
    //this three steps takes 0.02 s, 0.01 s in optim
    toc("");
    printf("%g\n",res3);
}*/
/*
static void test_spmul(){
    dmat *w=dnew(10,1);
    dset(w,2);
    dshow(w);
    dsp *spw=dspnewdiag(w->nx,P(w),1);
    dspdisp(spw);
    rand_t stat;
    seed_rand(&stat,1);
    dsp *A=dspnewrandu(10,10,10,0.2,&stat);
    dsp *Aw=NULL;
    dsp *Aw2=NULL;
    Aw2=dspmulsp(A,spw);

    dspdisp(Aw);
    dspdisp(Aw2);
    dsp *diff=dspadd2(Aw,1,Aw2,-1);
    dspdisp(diff);
    }*/
/*
static void test_spsum(){//passed
    dspcell *FLMc=dspcellread("FLM.bin");
    dsp *FLM=P(FLMc,0);
    dmat *sum1=dspsum(FLM,1);
    writebin(sum1,"sum_1");   dfree(sum1);
    sum1=dspsum(FLM,2);
    writebin(sum1,"sum_2");   dfree(sum1);
    sum1=dspsumabs(FLM,1);
    writebin(sum1,"sumabs_1"); dfree(sum1);
    sum1=dspsumabs(FLM,2);
    writebin(sum1,"sumabs_2"); dfree(sum1);
    dspcellfree(FLMc);
    }*/
/*static void test_L2(){
    dsp*L2=dspread("L2tmp.bin");
    dspcheck(L2);
    }*/
static void test_spmul(){
    TIC;
    long nx=2048;
    long ny=1024;
    rand_t rstat;
    seed_rand(&rstat,1);
    dsp *A=dspnewrandu(nx,ny,1,0.01,&rstat);
    dspaddI(A,1);
    dmat *x=dnew(A->ny,1);
    drandn(x,1,&rstat);
    dmat *y=dnew(A->nx,1);
    dbg("P(x)=%p\n",x);
    tic;
    dspmulvec(P(y),A,P(x),'n',1);
    toc("mul");
    writebin(x,"x");
    writebin(y,"y");
    dzero(y);
    tic;
    dspmulvec(P(y),A,P(x),'n',1);
    toc("mul");;
    dzero(y);
    /*
    tic;
    spmulvec_mkl(P(y),A,P(x),1);
    toc("mul_mkl");
    writebin(y,"y_mkl");
    */
    dzero(y);
    tic;
    dspmulvec(P(y),A,P(x),'n',1);
    toc("mul");
    dzero(y);
    tic;
    dspmulvec(P(y),A,P(x),'n',2);
    toc("mul 2");
    writebin(y,"y2");
    tic;
    dspmulvec(P(y),A,P(x),'n',2);
    toc("mul 2");
    tic;
    dspmulvec(P(y),A,P(x),'n',3);
    toc("mul 3");
    tic;
    dspmulvec(P(y),A,P(x),'n',4);
    toc("mul 4");
    dzero(x);
    tic;
    dspmulvec(P(x),A,P(y),'t',1);
    toc("sptmul");
    writebin(x,"x");
    tic;
    dspmulvec(P(x),A,P(y),'t',1);
    toc("sptmul 1");
    tic;
    dspmulvec(P(x),A,P(y),'t',1);
    toc("sptmul 2");
    writebin(x,"x2");
    dzero(x);
    tic;
    dspmulvec(P(x),A,P(y),'t',1);
    toc("sptmul 3");
    tic;
    dspmulvec(P(x),A,P(y),'t',1);
    toc("sptmul 4");
    tic;
    dspmulvec(P(x),A,P(y),'t',1);
    toc("sptmul 5");
    exit(0);
}

int main(){
    THREAD_POOL_INIT(2);
    test_spmul();
    /*test_L2(); */
    /*test_spsum(); */
    /*test_spmul();    */
}

#include "../lib/aos.h"
//TIC;
/*
void benchmark()
    tic;
    spcell *G0=spcellread("G0.bin.gz");
    toc("");tic;
    spcell *G0T=spcelltrans(G0);
    toc("");tic;
    spcell *A=spcellmulspcell(G0T,G0,1);//this is slow.
    toc("");
    dcell *x=dcellread("rhs.bin.gz");
    tic;
    double res=spcellwddot(x,A,x);//this takes 0.15s in debugging 0.04 in optim
    toc("");
    printf("%g\n",res);
    tic;
    dcell *y=NULL;
    spcellmulmat(&y, A,x,1);
    double res2=dcellinn(x,y);//this takes 0.15s in debugging 0.04 in optim
    toc("");
    printf("%g\n",res2);
    //dcell *saneai=dcellread("saneai.bin.gz");
    tic;
    dcell *z=NULL;
    spcellmulmat(&z, G0,x,1);toc("");
    dcell *z2=NULL;
    spcellmulmat(&z2, G0T,z,1);toc("");
    double res3=dcellinn(x,z2);
    //this three steps takes 0.02 s, 0.01 s in optim
    toc("");
    printf("%g\n",res3);
}*/
/*
static void test_spmul(){
    dmat *w=dnew(10,1);
    dset(w,2);
    dshow(w);
    dsp *spw=spnewdiag(w->nx,w->p,1);
    spdisp(spw);
    rand_t stat;
    seed_rand(&stat,1);
    dsp *A=spnewrandu(10,10,10,0.2,&stat);
    dsp *Aw=NULL;
    dsp *Aw2=NULL;
    Aw2=spmulsp(A,spw);
    
    spdisp(Aw);
    spdisp(Aw2);
    dsp *diff=spadd2(Aw,Aw2,1,-1);
    spdisp(diff);
    }*/
/*
static void test_spsum(){//passed
    spcell *FLMc=spcellread("FLM.bin.gz");
    dsp *FLM=FLMc->p[0];
    dmat *sum1=spsum(FLM,1);
    dwrite(sum1,"sum_1");   dfree(sum1);
    sum1=spsum(FLM,2);
    dwrite(sum1,"sum_2");   dfree(sum1);
    sum1=spsumabs(FLM,1);
    dwrite(sum1,"sumabs_1"); dfree(sum1);
    sum1=spsumabs(FLM,2);
    dwrite(sum1,"sumabs_2"); dfree(sum1);
    spcellfree(FLMc);
    }*/
/*static void test_L2(){
    dsp*L2=spread("L2tmp.bin.gz");
    spcheck(L2);
    }*/
static void test_spmul(){
    TIC;
    spcell *Ac=spcellread("RLM.bin");
    dsp *A=Ac->p[0];
    rand_t rstat;
    seed_rand(&rstat,1);
    dmat *x=dnew(A->n,1);
    drandn(x,1,&rstat);
    dmat *y=dnew(A->m,1);
    info("x->p=%p\n",x);
    tic;
    spmulvec(y->p,A,x->p,1);
    toc("mul");
    dwrite(x,"x"); 
    dwrite(y,"y");
    dzero(y);
    tic;
    spmulvec(y->p,A,x->p,1);
    toc("mul");;
    dzero(y);
    /*
    tic;
    spmulvec_mkl(y->p,A,x->p,1);
    toc("mul_mkl");
    dwrite(y,"y_mkl");
    */
    dzero(y);
    tic;
    spmulvec_thread(y->p, A, x->p, 1, 1);
    toc("mul_thread");
    dzero(y);
    tic;
    spmulvec_thread(y->p, A, x->p, 1, 2);
    toc("mul_thread 2");
    dwrite(y,"y2");
    tic;
    spmulvec_thread(y->p, A, x->p, 1, 2);
    toc("mul_thread 2");
    tic;
    spmulvec_thread(y->p, A, x->p, 1, 3);
    toc("mul_thread 3");
    tic;
    spmulvec_thread(y->p, A, x->p, 1, 4);
    toc("mul_thread 4");
    dzero(x);
    tic;
    sptmulvec(x->p, A, y->p, 1);
    toc("sptmul");
    dwrite(x,"x");
    tic;
    sptmulvec_thread(x->p, A, y->p, 1, 1);
    toc("sptmul_thread 1");
    tic;
    sptmulvec_thread(x->p, A, y->p, 1, 2);
    toc("sptmul_thread 2");
    dwrite(x,"x2");
    dzero(x);
    tic;
    sptmulvec_thread(x->p, A, y->p, 1, 3);
    toc("sptmul_thread 3");
    tic;
    sptmulvec_thread(x->p, A, y->p, 1, 4);
    toc("sptmul_thread 4");
    tic;
    sptmulvec_thread(x->p, A, y->p, 1, 5);
    toc("sptmul_thread 5");
    
}
int main(){
    THREAD_POOL_INIT(2);
    test_spmul();
    //test_L2();
    //test_spsum();
    //test_spmul();   
}

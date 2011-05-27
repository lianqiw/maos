#include "../lib/aos.h"
/*
int main(){
    spcell *FLM=spcellread("FLM.bin.gz");
    dsp *F=FLM->p[0];
    spchol *Fchol=chol_factorize(F);
 
    rand_t rrand;
    seed_rand(&rrand,1);
    dmat *y=dnew(F->m,10);
    drandu(y,1,&rrand);
    dmat *x=NULL;
    chol_solve(&x,Fchol,y);
    dwrite(y,"y");
    dwrite(x,"x");
    chol_save(Fchol,"Chol");
    chol_free(Fchol);
    dfree(x);
    dfree(y);
    spcellfree(FLM);
}
*/
TIC;
int main(int argc, char* argv[]){
    //dsp *RLMc1=spread("RLMc_old.bin");
    if(argc!=2){
	error("Need 1 argument\n");
    }
    spcell *RLM=spcellread("%s",argv[1]);
    dsp *RLMc=spcell2sp(RLM);
    tic;info2("chol ...");
    spchol *R1=chol_factorize(RLMc);
    toc("done");
    rand_t rstat;
    seed_rand(&rstat,1);
    dmat *y=dnew(RLMc->m, 1);
    drandn(y, 1, &rstat);
    dmat *x=NULL, *x2=NULL, *x3=NULL;
    dsp *S1=NULL;
    long *perm=NULL;
    chol_convert(R1, 1);
    tic;
    chol_solve(&x, R1, y);
    toc("cholmod");tic;
    chol_solve(&x, R1, y);
    toc("cholmod");tic;
    chol_solve_upper(&x3, R1, y);
    toc("upper");tic;
    chol_solve_upper(&x3, R1, y);
    toc("upper");tic;
    chol_solve_lower(&x2, R1,y);
    toc("lower");tic;
    chol_solve_lower(&x2, R1,y);
    toc("lower");tic;
    chol_solve(&x, R1, y);
    toc("cholmod");tic;
    chol_solve(&x, R1, y);
    toc("cholmod");tic;
    dwrite(y,"y");
    dwrite(x,"x");
    dwrite(x2,"x2");
    dwrite(x3,"x3");
    chol_free(R1);
    spfree(RLMc);
    spcellfree(RLM);
}

#include "../lib/aos.h"

int main(){
    spcell *FLM=spcellread("FLMc.bin.gz");
    //dsp *F=FLM->p[0];
    dsp *FF=spcell2sp(FLM);
    spwrite(FF,"FF.bin.gz");
    dmat *eig=NULL, *eigv=NULL;
    int neig=1;
    dmat *y=dnew(FF->n,1);
    struct_rand rrand;
    seed_rand(&rrand,1);
    drandu(y,1,&rrand);
    dmat *x=dnew(FF->m,1);
    dmat *x2=dnew(FF->m,1);
    spmulvec(x->p,FF,y->p,1);
    spcellmulvec(x2->p,FLM,y->p,1);
    dwrite(y,"eig_y.bin.gz");
    dwrite(x,"eig_x.bin.gz");
    dwrite(x2,"eig_x2.bin.gz");
    eigs(&eig,&eigv,FF,A_SPARSE,neig,"LA");
    info("Max eig is %g\n",eig->p[neig-1]);
    info("Max eig 2 is %g\n",spmaxeig(FF));
    info("Max eigcell is %g\n",spcellmaxeig(FLM));
    dwrite(eig,"eig");
    dwrite(eigv,"eigv");
    spcellfree(FLM);
    spfree(FF);
}





#include <sys/mman.h>
#include "../lib/aos.h"
static void test_dpinv(){
    rand_t rstat;
    seed_rand(&rstat,1);
    dmat *A=dnew(10,4);
    drandn(A,1,&rstat);
    dmat *w=dnew(10,1);
    dset(w,2);
    dmat *Ap=dpinv(A,w);
    writebin(A,"A");
    writebin(Ap,"Ap");
    dmat *ApA=NULL;
    dmm(&ApA,0,Ap,A,"nn",1);
    dsp *spw=dspnewdiag(w->nx, w->p, 1);
    dmat *Ap2=dpinv(A,spw);
    writebin(w,"w");
    dspwrite(spw,"spw");
    writebin(Ap2,"Ap2");
    dfree(A); dfree(w); dfree(Ap); dfree(ApA); dspfree(spw); dfree(Ap2);
}
TIC;
static void test_dcp(){
    dmat *A=dnew(10240,1);
    dmat *B=dnew(10240,1);
    tic;
    for(int i=0; i<100;i++){
	dcp(&A,B);
    }
    toc("dcp");
    tic;
    for(int i=0; i<100; i++){
	memcpy(A->p,B->p,A->nx*A->ny*sizeof(double));
    }
    toc("cpy");
    tic;
    for(int i=0; i<100;i++){
	dcp(&A,B);
    }
    toc("dcp");
    dfree(A);
    dfree(B);
}
static void test_dref(){
    dmat *A=dnew(10,10);
    dmat *B=dref(A);
    dmat *C=dref(A);
    dfree(A);
    dfree(B);
    dfree(C);
}
static void test_dcircle(){
    dmat *A=dnew(100,100);
    double r=45;
    dcircle(A,50,50,1,1,r,1);
    writebin(A,"dcircle_linear");
    dzero(A);
    writebin(A,"dcircle");
    double r2=sqrt(dsum(A)/M_PI);
    info("r=%g, r2=%g\n",r,r2);
}
static void test_d2cell(){
    dcell *A=dcellread("bcell.bin");
    dmat *A2=dcell2m(A);
    writebin(A2,"bcell2m.bin");
    dcell *B=NULL;
    d2cell(&B,A2,A);
    writebin(B,"bcell2m2cell.bin");
    dcellfree(A);
    dfree(A2);

    A=dcellread("ccell.bin");
    long *dims=calloc(A->nx,sizeof(long));
    for(int ix=0; ix<A->nx; ix++){
	dims[ix]=A->p[ix]->nx;
    }
    A2=dcell2m(A);
    dcell *B2=d2cellref(A2,dims,A->nx);
    writebin(B2,"ccell2m2cell2.bin");
    
}
static void test_dshift2center(){

    dcell *pis=dcellread("pistat_seed1_wfs6.bin");
    cmat *B=NULL;
    for(int ip=0; ip<pis->nx; ip++){
	ccpd(&B,pis->p[ip]);
	cshift2center(B,0.5,0.5);
	dshift2center(pis->p[ip],0.5,0.5);
    }
}
static void test_clip(){
    int N=64;
    dmat *A=dnew(N,N);
    rand_t rstat;
    seed_rand(&rstat,1);
    drandn(A,1,&rstat);
    writebin(A,"A1");
    dclip(A,-0.5, 0.5);
    writebin(A,"A2");
    /*exit(0); */
}
static void test_hist(){
    int N=64;
    dmat *A=dnew(N,N);
    dmat *B=NULL;
    rand_t rstat;
    seed_rand(&rstat,1);
    for(int i=0; i<100; i++){
	drandn(A, 1, &rstat);
	dhistfill(&B, A, 0, 0.1, 100);
    }
    writebin(B,"hist");
    /*exit(0); */
}
static void test_dcellpinv(){
    dcell *TT=dcellread("TT.bin");
    dspcell *saneai=dspcellread("saneai.bin");
    dcell *PTT=dcellpinv(TT,saneai);
    writebin(PTT,"TT_pinv.bin");
    dcell *PTT2=dcellnew2(TT);
    for(int i=0; i<PTT->ny; i++){
	PTT2->p[i+i*PTT->nx]=dpinv(TT->p[i+i*TT->nx],saneai->p[i+i*TT->nx]);
    }
    writebin(PTT2,"TT2_pinv.bin");
    /*exit(0); */
}
static void test_dcellcat(int argc, char**argv){
    if(argc!=2){
	error("Usage: %s FILE.bin", argv[0]);
    }
    dcell *A=dcellread("%s",argv[1]);
    dcell *B=dcellreduce(A,1);
    dcell *C=dcellreduce(A,2);
    writebin(B,"B");
    writebin(C,"C");
    exit(0);
}
static void test_save(void){/*passed */
    /*info("page size is %ld\n",sysconf(_SC_PAGE_SIZE)); */
    dmat *a=dnew_mmap(20,20,NULL,"a");
    rand_t rstat;
    seed_rand(&rstat,1);
    drandn(a,1,&rstat);
    writebin(a,"a2.bin");
    long nnx[6]={2,3,4,5,6,0};
    long nny[6]={3,4,2,5,6,3};
    dcell *b=dcellnew_mmap(2,3, nnx, nny, NULL,NULL,"ac.bin");
    for(int ix=0; ix<b->nx*b->ny; ix++){
	drandn(b->p[ix], 1, &rstat);
    }
    writebin(b, "ac2.bin");
    dcellfree(b);
    dfree(a);
    exit(0);
}
static void test_spline(void){
    long nx=20;
    dmat *x=dnew(nx,1);
    dmat *y=dnew(nx,1);
    dmat *x2=dnew(nx*100,1);
    for(int i=0; i<nx; i++){
	x->p[i]=i;
	y->p[i]=sin(i*M_PI/2.);
    }
    for(int i=0; i<nx*100; i++){
	x2->p[i]=0.01*i;
    }
    dmat *coeff=dspline_prep(x,y);
    dmat *y2=dspline_eval(coeff,x,x2);
    dmat *y3=dspline(x,y,x2);
    writebin(x,"x");
    writebin(y,"y");
    writebin(x2,"x2");
    writebin(y2,"y2");
    writebin(y3,"y3");
    writebin(coeff,"coeff");
    exit(0);
}
static void test_spline_2d(void){
    long nx=20;
    long ny=30;
    dmat *x=dnew(nx,1);
    dmat *y=dnew(ny,1);
    dmat *z=dnew(nx,ny);
    for(int i=0; i<ny; i++){
	y->p[i]=i*1.2+10;
    }
    for(int i=0; i<nx; i++){
	x->p[i]=i*2+10;
    }
    PDMAT(z,pz);
    for(int iy=0; iy<ny; iy++){
	for(int ix=0; ix<nx; ix++){
	    pz[iy][ix]=sin(ix*M_PI/14.2)*sin(iy*M_PI/12.3);
	}
    }
    long nxnew=nx*11;
    long nynew=ny*11;
    dmat *xnew=dnew(nxnew,nynew);
    dmat *ynew=dnew(nxnew,nynew);
    PDMAT(xnew,pxnew);
    PDMAT(ynew,pynew);
    for(int iy=0; iy<nynew; iy++){
	for(int ix=0; ix<nxnew; ix++){
	    pxnew[iy][ix]=ix*0.2+6;
	    pynew[iy][ix]=iy*0.12+6;
	}
    }
    dcell *coeff=dbspline_prep(x,y,z);
    dmat *z22=dbspline_eval(coeff,x,y,xnew,ynew);

    writebin(x,"x");
    writebin(y,"y");
    writebin(z,"z");

    writebin(xnew,"xnew");
    writebin(ynew,"ynew");
    writebin(z22,"z22");
    writebin(coeff,"coeff");
    exit(0);
}
static void test_svd(void){
    dspcell *a=dspcellread("SVD");
    dmat *A=NULL;
    dspfull(&A, a->p[0], 1);
    if(0){
	dmat *U, *S, *VT;
	tic;
	dsvd(&U, &S, &VT, A);
	toc("dsvd");
	tic;
	writebin(U,"U.bin");
	writebin(S, "S.bin");
	writebin(VT,"VT.bin");
    }else{
	writebin(A,"A.bin");
	tic;
	dsvd_pow(A, -1, 1e-15);
	toc("dsvd_pow");
	writebin(A,"AI.bin");
	tic;
	dsvd_pow(A, -1, 1e-15);
	toc("dsvd_pow svd");
	writebin(A,"A2.bin");
    }
    exit(0);
}
static void test_psd1d(){
    dmat *tmp=dnew(100,1);
    for(int i=0; i<100; i++){
	tmp->p[i]=sin(i/10.);
    }
    dmat *psd=psd1d(tmp, 10);
    writebin(tmp, "tmp");
    writebin(psd, "psd");
    exit(0);
}
static void test_svd2(void){
    cmat *A=cread("SVD.bin");
    csvd_pow(A, -1, 1.e-7);
    writebin(A, "SVDI");
}
static void test_kalman(){
    //dmat *psd=dread("MODE_TT");
    info("sde_fit\n");
    //dmat *coeff0=dread("coeff0");
    //dmat *coeff=sde_fit(psd, coeff0, 0.1, 0, 1e5, 0);
    dmat *coeff=dread("coeff");
    //dmat *Gwfs=dnew(1,1);daddI(Gwfs,1);
    //dmat *Rwfs=dnew(1,1);Rwfs->p[0]=1;
    dcell *Gwfs=dcellread("Gwfs");
    dcell *Rwfs=dcellread("Rwfs");
    dmat *dtrat_wfs=dread("dtrat_wfs");
    dmat *proj=dread("Proj");
    kalman_t *k=sde_kalman(coeff, 1./800, dtrat_wfs, Gwfs, Rwfs, proj);
    writebin(k->Ad, "mex_Ad");
    writebin(k->Cd, "mex_Cd");
    writebin(k->M, "mex_M");
    writebin(k->P, "mex_P");
    /*rand_t rstat; seed_rand(&rstat, 1);
    //dmat *ts=psd2time(psd, &rstat, 1./800, 5000);
    dmat *ts=dread("ts");
    for(int i=0; i<100; i++){
	dmat *res=kalman_test(k, ts);
	writebin(res, "kalman_res");
	dfree(res);
	}*/
    exit(0);
}
/*static void test_reccati(){
    dmat *A=dread("test_A");
    dmat *Qn=dread("test_Qn");
    dmat *C=dread("test_C");
    dmat *Rn=dread("test_Rn");
    dmat *M=0, *P=0;
    M=reccati(&P, A, Qn, C, Rn);
    writebin(M, "test_M");
    writebin(P, "test_P");
    exit(0);
}*/
static void test_expm(){
    dmat *A=dread("expm_in");
    dmat *B=0;
    dexpm(&B, 1, A, 1);
    writebin(B, "expm_out");
    exit(0);
}
static void test_mm(){
    dcell *A=cellnew(1,1);
    dcell *B=cellnew(1,1);
    A->p[0]=dnew(3,4);
    dset(A->p[0],1);
    B->p[0]=dnew(4,5);
    dset(B->p[0],1);
    dcell *C=dcellmm2(A,B,"nn");
    writebin(A, "A");
    writebin(B, "B");
    writebin(C, "C");
    exit(0);
}
int main(int argc, char **argv){
    test_mm();
    test_expm();
    //test_reccati();
    test_kalman();
    test_svd2();
    test_svd();
    test_psd1d();
    test_spline_2d();
    test_spline();
    test_dpinv();

    test_save();
    test_dcellcat(argc, argv);
    test_dcellpinv();
    test_hist();
    test_clip();
    test_dshift2center();
    test_dcircle();
    test_d2cell();
    test_dref();

    test_dcp();
}

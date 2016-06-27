void TomoR(dcell **xout, const void *A, 
	   const dcell *xin, const double alpha){
    
    const RECON_T *recon=(const RECON_T *)A;
    dcell *x2=NULL;
    dcellcp(&x2, xin);
    TTFR(x2, recon->TTF, recon->PTTF);
    dcell *x3=NULL;
    dspcellmulmat_thread(&x3,recon->saneai, x2, 1,recon->nthread);
    sptcellmulmat_thread(xout,recon->G0tomo,x3,alpha,recon->nthread);
    dcellfree(x2);
    dcellfree(x3);
}
void TomoL(dcell **xout, const void *A, 
	   const dcell *xin, const double alpha){

    const RECON_T *recon=(const RECON_T *)A;
    dcell *gg=NULL;
    dspcellmulmat_thread(&gg, recon->G0tomo, xin, 1.,recon->nthread);
    TTFR(gg, recon->TTF, recon->PTTF);
    dcell *gg2=NULL;
    dspcellmulmat_thread(&gg2, recon->saneai, gg,1,recon->nthread);
    sptcellmulmat_thread(xout,recon->G0tomo,gg2,alpha,recon->nthread);
    dcellfree(gg);
    dcellfree(gg2);
    if(recon->L2){
	apply_L2(xout, recon->L2, xin, alpha, recon->nthread);
    }else{
	apply_invpsd(xout, recon->xloc_embed, recon->invpsd,
		     recon->fftxopd, xin, alpha);
    }
    if(recon->ZZT){/*single point piston constraint */
	sptcellmulmat_thread(xout, recon->ZZT, xin, alpha,
			     recon->nthread);
    }
    /*Tikhonov regularization is not added because it is not
      necessary in CG mode.*/
}

void MUV(dcell **xout, const void *A, 
	 const dcell *xin, const double alpha){
    /**
       Apply the sparse plug low rand compuation to xin
       without scaling of alpha:
       xout=(muv.M-muv.U*muv.V')*xin*alpha; U,V are low
       rank.
    */
    const MUV_T *muv=(const MUV_T *)A;

    dspcellmulmat(xout, muv->M, xin, alpha);
    dcell *tmp=NULL;
    dcellmm(&tmp,muv->V, xin, "tn", -1.);
    dcellmm(xout,muv->U, tmp, "nn", alpha);
    dcellfree(tmp);
    if(muv->invpsd){
	apply_invpsd(xout, muv->xembed, muv->invpsd, muv->fftxopd, xin, alpha);
    }
}

TIC;
static int test_tomo(){
    info("Tomo\n");
    dcell *grad=dcellread("grad.bin");
    MUV_T RR;
    RR.M=dspcellread("RRM.bin"); 
    RR.U=dcellread("RRU.bin");
    RR.V=dcellread("RRV.bin");
    tic;
    dcell *rhs=NULL;
    MUV(&rhs, &RR, grad, 1);
    toc("");
    writebin(rhs,"rhs.bin");
    MUV_T RL;
    RL.M=dspcellread("RLM.bin");
    RL.U=dcellread("RLU.bin");
    RL.V=dcellread("RLV.bin");
    dcell *junk=NULL;
    tic; 
    for(int i=0; i<1; i++){
	MUV(&junk, &RL, rhs, 1);
    } 
    toc("");
    writebin(junk,"MUV1.bin");
    RECON_T *recon=mycalloc(1,RECON_T);
    recon->G0=dspcellread("G0.bin");
    recon->TT=dcellread("TT.bin");
    recon->PTT=dcellread("PTT.bin");
    recon->L2=dspcellread("L2.bin");
    recon->ZZT=dspcellread("ZZT.bin");
    recon->saneai=dspcellread("saneai.bin");
    tic;
    dcell *rhs2=NULL;
    TomoR(&rhs2, recon, grad, 1);
    toc("");
    writebin(rhs2,"rhs2.bin");
    dcell *junk2=NULL;
    tic;
    TomoL(&junk2, recon, rhs, 1);
    toc("");
    writebin(junk2,"TomoL1.bin");
    info("Diff between rhs is %g\n", dcelldiff(rhs, rhs2));
    info("Diff between lhs is %g\n", dcelldiff(junk,junk2));

    return 0;
}
static int test_fit(){
    info("Fit\n");
    dcell *opdr=dcellread("opdr.bin");
    MUV_T FR;
    FR.M=dspcellread("FRM.bin");
    FR.U=dcellread("FRU.bin");
    FR.V=dcellread("FRV.bin");
    MUV_T FL;
    FL.M=dspcellread("FLM.bin");
    FL.U=dcellread("FLU.bin");
    FL.V=dcellread("FLV.bin");
    dcell *rhs=NULL;
    tic;
    for(int i=0; i<10; i++)
	MUV(&rhs, &FR, opdr, 1);
    toc("");
    dcell *MUV_f=NULL;
    tic;
    for(int i=0; i<10; i++)
	MUV(&MUV_f, &FL, rhs, 1);
    toc("");
    writebin(rhs,"fit_rhs1.bin");
    writebin(MUV_f,"MUV_f.bin");
    RECON_T *recon=mycalloc(1,RECON_T);
    recon->HX=dspcellread("HX.bin");
    recon->HA=dspcellread("HA.bin");
    recon->W1=dread("W1.bin");
    recon->W0=dspread("W0.bin");
    recon->NW=dcellread("NW.bin");
    recon->fitwt=dread("fitwt.bin");
    dcell *rhs2=NULL;
    tic;
    for(int i=0; i<10; i++)
	FitR(&rhs2, recon, opdr, 1);
    toc("");
    writebin(rhs2,"fit_rhs2.bin");
    tic;
    dcell *FitL_f=NULL;
    for(int i=0; i<10; i++)
	FitL(&FitL_f, recon, rhs2, 1);
    toc("");
    writebin(FitL_f,"FitL_f.bin");
    info("Diff between rhs is %g\n", dcelldiff(rhs, rhs2));
    info("Diff between lhs is %g\n", dcelldiff(MUV_f, FitL_f));
    dcellfree(rhs);
    dcellfree(MUV_f);
    dcellfree(rhs2);
    dcellfree(FitL_f);
    return 0;
}
int main(){
    test_tomo();
    test_fit();
    /*tested ok.  */
    /*For tomography, MUV is slower than assembly on the fly. */
    /*For fit, MUV is faster than assembly on the fly primarily due to many fit directions. */
}
#endif

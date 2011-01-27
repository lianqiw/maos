/*
  This file contains the code of an unsuccess attempt to unwrap the wavefront for interpolation. It failed because for partially compensated turbulence, the OPD between neighboring pixels (in such coarsely sampled pupil function) may be very large.
*/


static void do_unwrap(cmat *phi, cmat *wvf, dmat *unwrap, dmat *diff, dmat *phirecon){
    /*
      Do the actual unwrapping. We do not check the
      dimension.s The Caller must make sure the matrices
      agree
     */
    int npsf=wvf->nx;
    PCMAT(wvf, pwvf);
    PDMAT(diff,pdiff);
    //TIC;tic;
    for(int ix=1; ix<npsf; ix++){
	pdiff[0][ix]=carg(pwvf[0][ix]*conj(pwvf[0][ix-1]));
	pdiff[ix][npsf]=carg(pwvf[ix][0]*conj(pwvf[ix-1][0]));
	for(int iy=1; iy<npsf; iy++){
	    pdiff[iy][ix]=carg(pwvf[iy][ix]*conj(pwvf[iy][ix-1]));
	    pdiff[iy][ix+npsf]=carg(pwvf[iy][ix]*conj(pwvf[iy-1][ix]));
	}
    }
    //toc("assemble");tic;
    dzero(phirecon);
    //dwrite(diff,"diff");
    dmulvec(phirecon->p, unwrap, diff->p, 1);
    //toc("mul");tic;
    //assert(phi->nx==npsf && phi->ny==npsf && npsf*npsf==unwrap->nx);
    for(int ix=0; ix<npsf*npsf; ix++){
	phi->p[ix]=phirecon->p[ix]*I+log(cabs(wvf->p[ix]));//real part saves amplitude.
    }
    //toc("assign");
}

static void convert_wvf(GENPISTAT_S *data){
    const PARMS_S *parms=data->parms;
    //POWFS_S *powfs=data->powfs;
    long icase=0;
    //Do not run this function.
    return;
    while((icase=lockadd(&data->icase, 1))<data->ncase){
    TIC;tic;
    double thetax=data->ngsgrid*data->cases[icase][0];
    double thetay=data->ngsgrid*data->cases[icase][1];
    long ipowfs=data->cases[icase][2];
    long ncomp=parms->maos.ncomp[ipowfs];
    long seed=data->cases[icase][3];
    long msa=parms->maos.msa[ipowfs];//in 1-d
    char fnwvf[PATH_MAX],fnphase[PATH_MAX];
    mymkdir("%s/phase",dirstart);
    snprintf(fnwvf,PATH_MAX,"%s/wvfout/wvfout_seed%ld_sa%ld_x%g_y%g.bin",
	     dirstart,seed, msa, thetax, thetay);
    snprintf(fnphase,PATH_MAX,"%s/phase/phase_seed%ld_sa%ld_x%g_y%g.bin",
	     dirstart,seed, msa, thetax, thetay);
    if(!exist(fnwvf) || exist(fnphase)){
	continue;
    }
    info("processing %s\n", fnwvf);
    file_t *fp_wvf=zfopen(fnwvf,"rb");
    uint32_t magic;
    //zfread(&magic, sizeof(uint32_t),1,fp_wvf);
    magic=read_magic(fp_wvf, NULL);
    if(!iscell(magic)){
	error("expected data type: %u, got %u\n",(uint32_t)MCC_ANY,magic);
    }
    long nstep,junk;
    zfread(&nstep,sizeof(uint64_t),1,fp_wvf);
    zfread(&junk,sizeof(uint64_t),1,fp_wvf);
    cellarr *phase=cellarr_init(nstep,"%s",fnphase);
    const int nsa=msa*msa;
    const int nwvl=parms->maos.nwvl;
    ccell *phi=ccellnew(nsa,nwvl);
    const long npsf=ncomp/2;
    dmat *phirecon=dnew(npsf,npsf);
    dmat *diff=dnew(npsf*2,npsf);
    for(long ic=0; ic<nsa*nwvl; ic++){
	phi->p[ic]=cnew(npsf,npsf);
    }
    for(long istep=0; istep<nstep; istep++){
	LOCK(data->mutex_read);
	ccell *wvfi=ccellreaddata(fp_wvf, 0);
	UNLOCK(data->mutex_read);
	if(wvfi){
	    for(long ic=0; ic<nsa*nwvl; ic++){
		do_unwrap(phi->p[ic], wvfi->p[ic], data->unwrap->p[ipowfs], diff, phirecon);
	    }
	}
	cellarr_ccell(phase,phi);
	ccellfree(wvfi);
    }
    ccellfree(phi);
    dfree(diff);
    dfree(phirecon);
    toc2("Processing %s:", fnphase);
    }
}
  //convert wvf of a+bi to log(a+bi) for interpolation.
    /*
    data->unwrap=dcellnew(parms->maos.npowfs,1);
    for(int ipowfs=0; ipowfs<parms->maos.npowfs; ipowfs++){
	long ncomp=parms->maos.ncomp[ipowfs];
	data->unwrap->p[ipowfs]=gen_unwrap(ncomp/2,ncomp/2);
	}*/
    /*  
	//test phase unwrap.
	{
	rand_t rstat;
	seed_rand(&rstat,1);
	double wt=1;
	const int npsf=16;
	double pi2l=M_PI/2.2e-6;
	MAP_S **screen=vonkarman_screen(&rstat, 16, 16, 1/2., .2, 30, &wt, 1, 1);
	cmat *wvf=cnew(16,16);
	dmat *opd=dnew(16,16);
	for(int ix=0; ix<16*16; ix++){
	    opd->p[ix]=screen[0]->p[ix]*pi2l*2.;
	    wvf->p[ix]=cexp(opd->p[ix]*I);
	}
	cmat *phi=cnew(16,16);
	dmat *diff=dnew(npsf*2,npsf);
	dmat *phirecon=dnew(npsf,npsf);
	do_unwrap(phi,wvf,data->unwrap->p[0],diff,phirecon);
	cwrite(wvf,"wvf");
	cwrite(phi,"phi");
	dwrite(opd,"opd");
	//	exit(0);
	}
    data->icase=0;
    CALL(convert_wvf, data, parms->skyc.nthread);
    */


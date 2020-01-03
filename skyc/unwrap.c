/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
/*
  This file contains the code of an unsuccess attempt to unwrap the wavefront for interpolation. It failed because for partially compensated turbulence, the OPD between neighboring pixels (in such coarsely sampled pupil function) may be very large.
*/


static dmat* gen_unwrap(long nx, long ny){
    /*Create phase unwrap operator */
    char fn[PATH_MAX];
    snprintf(fn,PATH_MAX,"%s/.aos/unwrap_%ld_%ld.bin",HOME,nx,ny);
    /*if(exist(fn)){ */
    /*	cmat *out=cread(fn); */
    /*	return out; */
    /*} */
    dsp *Ht=dspnew(nx*ny,nx*ny*2, nx*ny*2*2);/*forward */
    long count=0;
    spint *pp=Ht->p;
    spint *pi=Ht->i;
    real *px=Ht->x;
    long col=0;
    for(long iy=0; iy<ny; iy++){
	long offx=iy*nx;
	/*X difference */
	for(long ix=0; ix<nx; ix++){
	    pp[col]=count;
	    if(ix==0 && iy==0){
		px[count]=1;/*piston constraint. */
		pi[count]=0;
		count++;
	    }
	    col++;
	    if(ix>0){
		pi[count]=ix-1+offx;
		px[count]=-1;
		count++;
		
		pi[count]=ix+offx;
		px[count]=1;
		count++;
		/*}else{
		pi[count]=ix+1+offx;
		px[count]=1;
		count++;
		
		pi[count]=ix+offx;
		px[count]=-1;
		count++;*/
	    }
	}
	/*Y difference */
	for(long ix=0; ix<nx; ix++){
	    pp[col]=count;
	    col++;
	    if(iy>0){
		pi[count]=ix+offx-nx;
		px[count]=-1;
		count++;
		
		pi[count]=ix+offx;
		px[count]=1;
		count++;
		/*}else{
		pi[count]=ix+offx+nx;
		px[count]=1;
		count++;
		
		pi[count]=ix+offx;
		px[count]=-1;
		count++;*/
	    }
	    
	}
    }
    pp[col]=count;
    
    dbg("col=%ld,count=%ld\n",col,count);
    dspsetnzmax(Ht,count);
    /*writebin(Ht,"Ht"); */
    dsp *H=dsptrans(Ht);
    /*writebin(H,"H"); */
    dsp *HtH=dspmulsp(Ht,H);
    dmat *cHtH=NULL;
    dspfull(&cHtH,HtH,'n',1);
    /*writebin(cHtH,"cHtH"); */
    /*caddI(cHtH,1e-9); */
    dmat *IcHtH=dinv(cHtH);
    /*writebin(IcHtH,"IcHtH"); */
    dmat *out=NULL;
    dmulsp(&out, IcHtH, Ht, 1);
    /*writebin(out,"HI"); */
    dfree(IcHtH);
    dfree(cHtH);
    dspfree(HtH);
    dspfree(H);
    dspfree(Ht);
    writebin(out,"%s",fn);
    return out;
}

static void do_unwrap(cmat *phi, cmat *wvf, dmat *unwrap, dmat *diff, dmat *phirecon){
    /*
      Do the actual unwrapping. We do not check the
      dimension.s The Caller must make sure the matrices
      agree
     */
    int npsf=wvf->nx;
    cmat*  pwvf=wvf;
    dmat* pdiff=diff;
    /*TIC;tic; */
    for(int ix=1; ix<npsf; ix++){
	P(pdiff,ix,0)=carg(P(pwvf,ix,0)*conj(P(pwvf,ix-1,0)));
	P(pdiff,npsf,ix)=carg(P(pwvf,0,ix)*conj(P(pwvf,0,ix-1)));
	for(int iy=1; iy<npsf; iy++){
	    P(pdiff,ix,iy)=carg(P(pwvf,ix,iy)*conj(P(pwvf,ix-1,iy)));
	    P(pdiff,ix+npsf,iy)=carg(P(pwvf,ix,iy)*conj(P(pwvf,ix,iy-1)));
	}
    }
    /*toc("assemble");tic; */
    dzero(phirecon);
    /*writebin(diff,"diff"); */
    dmulvec(phirecon->p, unwrap, diff->p, 1);
    /*toc("mul");tic; */
    /*assert(phi->nx==npsf && phi->ny==npsf && npsf*npsf==unwrap->nx); */
    for(int ix=0; ix<npsf*npsf; ix++){
	phi->p[ix]=COMPLEX(log(cabs(wvf->p[ix])), phirecon->p[ix]);/*real part saves amplitude. */
    }
    /*toc("assign"); */
}

static void convert_wvf(GENPISTAT_S *data){
    const PARMS_S *parms=data->parms;
    /*POWFS_S *powfs=data->powfs; */
    long icase=0;
    /*Do not run this function. */
    return;
    while((LOCKADD(icase,data->icase, 1))<data->ncase){
    TIC;tic;
    real thetax=data->ngsgrid*data->P(cases,0,icase);
    real thetay=data->ngsgrid*data->P(cases,1,icase);
    long ipowfs=data->P(cases,2,icase);
    long ncomp=parms->maos.ncomp[ipowfs];
    long seed=data->P(cases,3,icase);
    long msa=parms->maos.msa[ipowfs];/*in 1-d */
    char fnwvf[PATH_MAX],fnphase[PATH_MAX];
    mymkdir("%s/phase",dirstart);
    snprintf(fnwvf,PATH_MAX,"%s/wvfout/wvfout_seed%ld_sa%ld_x%g_y%g",
	     dirstart,seed, msa, thetax, thetay);
    snprintf(fnphase,PATH_MAX,"%s/phase/phase_seed%ld_sa%ld_x%g_y%g",
	     dirstart,seed, msa, thetax, thetay);
    if(!zfexist(fnwvf) || zfexist(fnphase)){
	continue;
    }
    dbg("processing %s\n", fnwvf);
    file_t *fp_wvf=zfopen(fnwvf,"rb");
    uint32_t magic;
    /*zfread(&magic, sizeof(uint32_t),1,fp_wvf); */
    magic=read_magic(fp_wvf, NULL);
    if(!iscell(&magic)){
	error("expected data type: %u, got %u\n",(uint32_t)MCC_ANY,magic);
    }
    long nstep,junk;
    zfread(&nstep,sizeof(uint64_t),1,fp_wvf);
    zfread(&junk,sizeof(uint64_t),1,fp_wvf);
    zfarr *phase=zfarr_init(nstep,1,"%s",fnphase);
    const int nsa=msa*msa;
    const int nwvl=parms->maos.nwvl;
    ccell *phi=cellnew(nsa,nwvl);
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
	zfarr_push(phase,phi);
	ccellfree(wvfi);
    }
    ccellfree(phi);
    dfree(diff);
    dfree(phirecon);
    toc("Processing %s:", fnphase);
    }
}
  /*convert wvf of a+bi to log(a+bi) for interpolation. */
    /*
    data->unwrap=cellnew(parms->maos.npowfs,1);
    for(int ipowfs=0; ipowfs<parms->maos.npowfs; ipowfs++){
	long ncomp=parms->maos.ncomp[ipowfs];
	data->unwrap->p[ipowfs]=gen_unwrap(ncomp/2,ncomp/2);
	}*/
    /*  
	/*test phase unwrap. */
	{
	rand_t rstat;
	seed_rand(&rstat,1);
	real wt=1;
	const int npsf=16;
	real pi2l=M_PI/2.2e-6;
	MAP_S **screen=vonkarman_screen(&rstat, 16, 16, 1/2., .2, 30, &wt, 1, 1);
	cmat *wvf=cnew(16,16);
	dmat *opd=dnew(16,16);
	for(int ix=0; ix<16*16; ix++){
	    opd->p[ix]=screen[0]->p[ix]*pi2l*2.;
	    wvf->p[ix]=cexp(COMPLEX(0, opd->p[ix]));
	}
	cmat *phi=cnew(16,16);
	dmat *diff=dnew(npsf*2,npsf);
	dmat *phirecon=dnew(npsf,npsf);
	do_unwrap(phi,wvf,data->unwrap->p[0],diff,phirecon);
	writebin(wvf,"wvf");
	writebin(phi,"phi");
	writebin(opd,"opd");
	/*	exit(0); */
	}
    data->icase=0;
    CALL(convert_wvf, data, parms->skyc.nthread);
    */


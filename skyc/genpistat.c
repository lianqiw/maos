/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
  calculates pistat using PSF and ztilt.
*/
#include <sys/stat.h>
#include <sys/types.h>
#include "skyc.h"
#include "parms.h"
#include "types.h"
#include "skysim_utils.h"
#include "genpistat.h"

typedef struct GENPISTAT_S{
    long ncase;
    long(*cases)[4];
    long icase;
#if USE_PTHREAD > 0
    pthread_mutex_t mutex_icase;
    pthread_mutex_t mutex_read;//don't let them read in the same time.
#endif
    const PARMS_S *parms;
    POWFS_S *powfs;
    double ngsgrid;
    dcell *unwrap;//matrix that does unwraping.
}GENPISTAT_S;

static void* calc_pistat(GENPISTAT_S *data){
    const PARMS_S *parms=data->parms;
    POWFS_S *powfs=data->powfs;
    long icase=0;
 retry:
    LOCK(data->mutex_icase);
    icase=data->icase++;
    UNLOCK(data->mutex_icase);
    if(icase>=data->ncase) return NULL;
    double thetax=data->ngsgrid*data->cases[icase][0];
    double thetay=data->ngsgrid*data->cases[icase][1];
    long ipowfs=data->cases[icase][2];
    long ncomp=parms->maos.ncomp[ipowfs];
    long seed=data->cases[icase][3];
    long msa=parms->maos.msa[ipowfs];//in 1-d
    const long phystart=parms->skyc.phystart;
    char fnwvf[PATH_MAX],fnztilt[PATH_MAX], fnpistat[PATH_MAX], fngstat[PATH_MAX];
    snprintf(fnwvf,PATH_MAX,"%s/wvfout/wvfout_seed%ld_sa%ld_x%g_y%g.bin",
	     dirstart,seed, msa, thetax, thetay);
    snprintf(fnztilt,PATH_MAX,"%s/ztiltout/ztiltout_seed%ld_sa%ld_x%g_y%g.bin",
	     dirstart,seed, msa, thetax, thetay);
    snprintf(fnpistat,PATH_MAX,"%s/pistat/pistat_seed%ld_sa%ld_x%g_y%g.bin.gz",
	     dirstart,seed, msa, thetax, thetay);
    snprintf(fngstat,PATH_MAX,"%s/pistat/gstat_seed%ld_sa%ld_x%g_y%g.bin.gz",
	     dirstart,seed, msa, thetax, thetay);
    if(!exist(fnwvf)){
	goto retry;
    }
    if(exist(fnpistat) && exist(fngstat)){
	//warning("%s already exist, skip\n", fnpistat);
    }else{
	dmat *mapply=dnew(2,1);
	char dirstat[PATH_MAX];
	snprintf(dirstat, PATH_MAX, "%s/pistat", dirstart);
	mymkdir(dirstat,0700);
	TIC;tic;
	if(!exist(fnztilt)){
	    error("%s exist, but %s doesn't exist\n", fnwvf, fnztilt);
	}
	file_t *fp_wvf=zfopen(fnwvf,"rb");
	uint32_t magic;
	zfread(&magic, sizeof(uint32_t),1,fp_wvf);
	if(magic!=MCC_CMP && magic !=MCC_ANY){
	    error("expected data type: %u, got %u\n",(uint32_t)MCC_CMP,magic);
	}
	long nstep,junk;
	zfread(&nstep,sizeof(uint64_t),1,fp_wvf);
	zfread(&junk,sizeof(uint64_t),1,fp_wvf);
		
	file_t *fp_ztilt=zfopen(fnztilt,"rb");
	zfread(&magic, sizeof(uint32_t),1,fp_ztilt);
	if(magic!=MC_DBL && magic !=MCC_ANY){
	    error("expected data type: %u, got %u\n",(uint32_t)MC_DBL,magic);
	} 
	long nstep2;
	zfread(&nstep2,sizeof(uint64_t),1,fp_ztilt);
	zfread(&junk,sizeof(uint64_t),1,fp_ztilt);
		
	if(nstep!=nstep2){
	    error("length of wvfout doesn't equal to length of ztiltout\n");
	}
	const int nsa=msa*msa;
	const int nwvl=parms->maos.nwvl;
	dcell *pistat=dcellnew(nsa, nwvl);//pixel intensity mean(I)
	dmat *gstat=dnew(nsa*2, nstep);//original gradient at each time step.
	cmat *wvf=cnew(ncomp,ncomp);
	cmat *wvfc=cnew(ncomp/2,ncomp/2);
	cfft2plan(wvf,-1);
	PDCELL(pistat, ppistat);
	PDMAT(gstat, pgstat);
	cmat *otf=cnew(ncomp,ncomp);
	cfft2plan(otf,1);
	cfft2plan(otf,-1);
	dmat *psf=NULL;
	double nwvli=1./(nwvl);
	for(long istep=0; istep<nstep; istep++){
	    LOCK(data->mutex_read);
	    ccell *wvfi=ccellreaddata(fp_wvf, 0);
	    dmat *ztilti=dreaddata(fp_ztilt, 0);
	    UNLOCK(data->mutex_read);
	    PCCELL(wvfi,wvfout);
	    if(istep>=phystart){
		for(long iwvl=0; iwvl<nwvl; iwvl++){
		    double wvl=parms->maos.wvl[iwvl];
		    for(long isa=0; isa<nsa; isa++){
			mapply->p[0]=-ztilti->p[isa];
			mapply->p[1]=-ztilti->p[isa+nsa];
			ccp(&wvfc, wvfout[iwvl][isa]);
			ngsmod2wvf(wvfc, wvl, mapply,powfs->cloc[isa],
				   powfs->fpc[isa], 
				   thetax, thetay, parms);
			cembed(wvf,wvfc,0,C_FULL);
			cfft2(wvf,-1);
			cabs22d(&ppistat[iwvl][isa],1, wvf, 1);
		    }
		}
		for(long iwvl=0; iwvl<nwvl; iwvl++){//don't apply ztilt.
		    for(long isa=0; isa<nsa; isa++){
			double grad[2]={0,0};
			ccp(&wvfc, wvfout[iwvl][isa]);
			cembed(wvf,wvfc,0,C_FULL);
			cfft2(wvf,-1);
			dzero(psf);
			cabs22d(&psf, 1, wvf, 1);
			//Compute gradients.
			dfftshift(psf);//shift peak to center.
			double pmax=dmax(psf);
			dcog(grad,psf,0.5, 0.5, 0.1*pmax, 0.1*pmax);
			pgstat[istep][isa]+=grad[0]*parms->skyc.pixtheta[ipowfs]*nwvli;
			pgstat[istep][isa+nsa]+=grad[1]*parms->skyc.pixtheta[ipowfs]*nwvli;
			dfree(psf);
		    }
		}
	    }
	    ccellfree(wvfi);
	    dfree(ztilti);  
	}//for istep
	zfeof(fp_wvf);
	zfeof(fp_ztilt);
	zfclose(fp_ztilt);
	zfclose(fp_wvf);
	cfree(wvf);
	cfree(wvfc);
	double pgrad[2];

	//Shift PSF so that the images are centered on FFT center.
	for(int i=0; i<nwvl*nsa; i++){
	    psf=pistat->p[i];
	    dfftshift(psf);//shift peak to center.
	    double pmax=dmax(psf);
	    dcog(pgrad,psf,0.5,0.5,0.1*pmax,0.1*pmax);
	    dfftshift(psf); //Shift peak to corner.
	    ccpd(&otf, psf);
	    cfft2(otf,-1);
	    ctilt(otf,-pgrad[0],-pgrad[1],0);
	    cifft2(otf,1);
	    creal2d(&psf,0,otf,1);
	}
	/*
	  Saved pistat should have peak on the corner
	 */
	cfree(otf);
	dcellscale(pistat, 1./(nstep-phystart));
	dcellwrite(pistat, "%s",fnpistat);
	dwrite(gstat, "%s",fngstat);
	dcellfree(pistat);
	dfree(gstat);
	dfree(mapply);
	toc2("Processing %s:", fnwvf);
    }//if exist
    goto retry;
}

static dmat* gen_unwrap(long nx, long ny){
    //Create phase unwrap operator
    char fn[PATH_MAX];
    snprintf(fn,PATH_MAX,"%s/.aos/unwrap_%ld_%ld.bin",HOME,nx,ny);
    //if(exist(fn)){
    //	cmat *out=cread(fn);
    //	return out;
    //}
    dsp *Ht=spnew(nx*ny,nx*ny*2, nx*ny*2*2);//forward
    long count=0;
    spint *pp=Ht->p;
    spint *pi=Ht->i;
    double *px=Ht->x;
    long col=0;
    for(long iy=0; iy<ny; iy++){
	long offx=iy*nx;
	//X difference
	for(long ix=0; ix<nx; ix++){
	    pp[col]=count;
	    if(ix==0 && iy==0){
		px[count]=1;//piston constraint.
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
	//Y difference
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
    
    info("col=%ld,count=%ld\n",col,count);
    spsetnzmax(Ht,count);
    //spwrite(Ht,"Ht");
    dsp *H=sptrans(Ht);
    //spwrite(H,"H");
    dsp *HtH=spmulsp(Ht,H);
    dmat *cHtH=NULL;
    spfull(&cHtH,HtH,1);
    //dwrite(cHtH,"cHtH");
    //caddI(cHtH,1e-9);
    dmat *IcHtH=dinv(cHtH);
    //dwrite(IcHtH,"IcHtH");
    dmat *out=NULL;
    dmulsp(&out, IcHtH, Ht, 1);
    //dwrite(out,"HI");
    dfree(IcHtH);
    dfree(cHtH);
    spfree(HtH);
    spfree(H);
    spfree(Ht);
    dwrite(out,"%s",fn);
    return out;
}

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

static void *convert_wvf(GENPISTAT_S *data){
    const PARMS_S *parms=data->parms;
    //POWFS_S *powfs=data->powfs;
    long icase=0;
    return NULL;
 retry:
    LOCK(data->mutex_icase);
    icase=data->icase++;
    UNLOCK(data->mutex_icase);
    if(icase>=data->ncase) return NULL;
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
	goto retry;
    }
    info("processing %s\n", fnwvf);
    file_t *fp_wvf=zfopen(fnwvf,"rb");
    uint32_t magic;
    zfread(&magic, sizeof(uint32_t),1,fp_wvf);
    if(magic!=MCC_CMP && magic !=MCC_ANY){
	error("expected data type: %u, got %u\n",(uint32_t)MCC_CMP,magic);
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
    //  exit(0);
    //(void)phase;(void)ncomp;(void)powfs;
    goto retry;
}
void genpistat(const PARMS_S *parms, POWFS_S *powfs){
    double patfov=parms->skyc.patfov;
    double ngsgrid=parms->maos.ngsgrid;
    long ng=ceil(patfov/2/ngsgrid);
    info("Genpistat\n");
    GENPISTAT_S *data=calloc(1, sizeof(GENPISTAT_S));
    data->parms=parms;
    data->powfs=powfs;
    data->ncase=parms->maos.nseed*(2*ng+1)*(2*ng+1)*parms->maos.npowfs;
    data->cases=calloc(4*data->ncase,sizeof(long));
    data->ngsgrid=ngsgrid;
    PINIT(data->mutex_icase);
    long count=0;
    for(int iseed=0; iseed<parms->maos.nseed; iseed++){
	int seed=parms->maos.seeds[iseed];//loop over seed
	for(long gy=-ng; gy<=ng; gy++){
	    for(long gx=-ng; gx<=ng; gx++){
		//double thetax=gx*ngsgrid;
		//double thetay=gy*ngsgrid;
		for(int ipowfs=0; ipowfs<parms->maos.npowfs; ipowfs++){//for ipowfs
		    data->cases[count][0]=gx;
		    data->cases[count][1]=gy;
		    data->cases[count][2]=ipowfs;
		    data->cases[count][3]=seed;
		    count++;
		    if(count>data->ncase){
			warning("Initial guess of %ld is too low: %ld\n",data->ncase,count);
			data->ncase=data->ncase*2;
			data->cases=realloc(data->cases,sizeof(long)*data->ncase*4);
		    }
		}
	    }//for gx
	}//for gy
    }//for iseed
    //info("count=%ld, data->ncase=%ld\n",count,data->ncase);
    data->ncase=count;
    data->icase=0;
    data->cases=realloc(data->cases, sizeof(long)*4*data->ncase);
    CALL(calc_pistat, data, parms->skyc.nthread);
    //convert wvf of a+bi to log(a+bi) for interpolation.
    data->icase=0;
    data->unwrap=dcellnew(parms->maos.npowfs,1);
    for(int ipowfs=0; ipowfs<parms->maos.npowfs; ipowfs++){
	long ncomp=parms->maos.ncomp[ipowfs];
	data->unwrap->p[ipowfs]=gen_unwrap(ncomp/2,ncomp/2);
    }
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
	}*/
    CALL(convert_wvf, data, parms->skyc.nthread);
    dcellfree(data->unwrap);
    free(data->cases);
    free(data);
}
/**
   Read in the average pixel intensities on the ngs grid
*/
void prep_bspstrehl(SIM_S *simu){
    const PARMS_S *parms = simu->parms;
    double patfov  = parms->skyc.patfov;
    double ngsgrid = parms->maos.ngsgrid;
    long ng  = ceil(patfov/2/ngsgrid);
    long ng2 = ng*2+1;
    long seed = simu->seed_maos;
    dmat *gg=dnew(ng*2+1,1);
    for(long gx=-ng; gx<=ng; gx++){
	gg->p[gx+ng]=(double)gx;
    }
    simu->bspstrehlxy=gg;

    simu->bspstrehl=calloc(parms->maos.npowfs, sizeof(dcell**));

    long ngnew=ng*10;
    dmat *xnew=dnew(ngnew*2+1, ngnew*2+1);
    dmat *ynew=dnew(ngnew*2+1, ngnew*2+1);
    for(int iy=0; iy<xnew->ny; iy++){
	for(int ix=0; ix<xnew->nx; ix++){
	    xnew->p[ix+iy*xnew->nx]=(ix-ngnew)*0.1;
	    ynew->p[ix+iy*xnew->nx]=(iy-ngnew)*0.1;
	}
    }
    dwrite(xnew,"xnew");
    dwrite(ynew,"ynew");
    for(int ipowfs=0; ipowfs<parms->maos.npowfs; ipowfs++){
	long msa=parms->maos.msa[ipowfs];
	long nsa=parms->maos.nsa[ipowfs];
	long nwvl=parms->maos.nwvl;
	dcell *strehlgrid = dcellnew(nsa,nwvl);
	for(long ic=0; ic<nsa*nwvl; ic++){
	    strehlgrid->p[ic]=dnew(ng2,ng2);
	}
	for(long gy=-ng; gy<=ng; gy++){
	    double thetay = gy * ngsgrid;
	    for(long gx=-ng; gx<=ng; gx++){
		double thetax = gx * ngsgrid;
		char fnpistat[PATH_MAX];
		snprintf(fnpistat,PATH_MAX,
			 "%s/pistat/pistat_seed%ld_sa%ld_x%g_y%g.bin.gz",
			 dirstart, seed, msa, thetax, thetay);

		if(exist(fnpistat)){
		    dcell *tmp=dcellread("%s",fnpistat);
		    for(long ic=0; ic<nsa*nwvl; ic++){
			//peak is in the corner
			strehlgrid->p[ic]->p[gx+ng+(gy+ng)*ng2]=tmp->p[ic]->p[0];
		    }
		}
	    }
	}

	simu->bspstrehl[ipowfs]=calloc(nsa*nwvl, sizeof(dcell*));
	dcellwrite(strehlgrid, "strehlgrid_%d",ipowfs);
	for(long ic=0; ic<nsa*nwvl; ic++){
	    simu->bspstrehl[ipowfs][ic]=dbspline_prep(gg,gg,strehlgrid->p[ic]);
	    dmat *strehlnew=dbspline_eval(simu->bspstrehl[ipowfs][ic],gg,gg,xnew,ynew);
	    dwrite(strehlnew, "strehlnew_%d_%ld", ipowfs, ic);
	    dfree(strehlnew);
	}
	dcellfree(strehlgrid);
    }
    dfree(xnew);
    dfree(ynew);
}

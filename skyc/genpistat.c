/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
    int ncase;
    long(*cases)[4];
    int icase;
#if USE_PTHREAD > 0
    pthread_mutex_t mutex_read;/*don't let them read in the same time. */
#endif
    const PARMS_S *parms;
    POWFS_S *powfs;
    double ngsgrid;
    dcell *unwrap;/*matrix that does unwraping. */
}GENPISTAT_S;

static void calc_pistat(GENPISTAT_S *data){
    const PARMS_S *parms=data->parms;
    POWFS_S *powfs=data->powfs;
    int icase=0;
    char fnlock[PATH_MAX];
    snprintf(fnlock, PATH_MAX, "%s/wvfout/wvfout.lock", dirstart);
    /*Obtain exclusive lock before proceeding. */
    //int fd=lock_file(fnlock, 1, 0);
    while(LOCKADD(icase, data->icase, 1)<data->ncase){
	double thetax=data->ngsgrid*data->cases[icase][0];
	double thetay=data->ngsgrid*data->cases[icase][1];
	long ipowfs=data->cases[icase][2];
	long ncomp=parms->maos.ncomp[ipowfs];
	long seed=data->cases[icase][3];
	long msa=parms->maos.msa[ipowfs];/*in 1-d */
	const long phystart=parms->skyc.phystart;
	char fnwvf[PATH_MAX],fnztilt[PATH_MAX], fnpistat[PATH_MAX], fngstat[PATH_MAX];
	snprintf(fnwvf,PATH_MAX,"%s/wvfout/wvfout_seed%ld_sa%ld_x%g_y%g",
		 dirstart,seed, msa, thetax, thetay);
	snprintf(fnztilt,PATH_MAX,"%s/ztiltout/ztiltout_seed%ld_sa%ld_x%g_y%g",
		 dirstart,seed, msa, thetax, thetay);
	snprintf(fnpistat,PATH_MAX,"%s/pistat/pistat_seed%ld_sa%ld_x%g_y%g",
		 dirstart,seed, msa, thetax, thetay);
	snprintf(fngstat,PATH_MAX,"%s/pistat/gstat_seed%ld_sa%ld_x%g_y%g",
		 dirstart,seed, msa, thetax, thetay);
	if(!zfexist(fnwvf)){
	    continue;
	}
	if(zfexist(fnpistat) && zfexist(fngstat)){
	    /*warning("%s already exist, skip\n", fnpistat); */
	}else{
	    dmat *mapply=dnew(2,1);
	    char dirstat[PATH_MAX];
	    snprintf(dirstat, PATH_MAX, "%s/pistat", dirstart);
	    mymkdir(dirstat,0700);
	    TIC;tic;
	    if(!zfexist(fnztilt)){
		error("%s exist, but %s doesn't exist\n", fnwvf, fnztilt);
	    }
	    file_t *fp_wvf=zfopen(fnwvf,"rb");
	    header_t header;
	    read_header(&header, fp_wvf);
	    long nstep=header.nx;
	    if(!iscell(header.magic)){
		error("expected data type: %u, got %u\n", (uint32_t)MCC_ANY, header.magic);
	    }
	    free(header.str); header.str=NULL;
	    file_t *fp_ztilt=zfopen(fnztilt,"rb");
	    read_header(&header, fp_ztilt);
	    if(!iscell(header.magic)){
		error("expected data type: %u, got %u\n",(uint32_t)MCC_ANY, header.magic);
	    } 
	    if(nstep!=header.nx){
		error("length of wvfout doesn't equal to length of ztiltout\n");
	    }
	    free(header.str); header.str=NULL;
	    const int nsa=msa*msa;
	    const int nwvl=parms->maos.nwvl;
	    dcell *pistat=dcellnew(nsa, nwvl);/*pixel intensity mean(I) */
	    dmat *gstat=dnew(nsa*2, nstep);/*original gradient at each time step. */
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
				       powfs->fpc[isa], thetax, thetay, parms);
			    cembed(wvf,wvfc,0,C_FULL);
			    cfft2(wvf,-1);
			    cabs22d(&ppistat[iwvl][isa],1, wvf, 1);
			}
		    }
		    for(long iwvl=0; iwvl<nwvl; iwvl++){/*don't apply ztilt. */
			for(long isa=0; isa<nsa; isa++){
			    double grad[2]={0,0};
			    ccp(&wvfc, wvfout[iwvl][isa]);
			    cembed(wvf,wvfc,0,C_FULL);
			    cfft2(wvf,-1);
			    dzero(psf);
			    cabs22d(&psf, 1, wvf, 1);
			    /*Compute gradients. */
			    dfftshift(psf);/*shift peak to center. */
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
	    }/*for istep */
	    zfeof(fp_wvf);
	    zfeof(fp_ztilt);
	    zfclose(fp_ztilt);
	    zfclose(fp_wvf);
	    cfree(wvf);
	    cfree(wvfc);
	    double pgrad[2];

	    /*Shift PSF so that the images are centered on FFT center. */
	    for(int i=0; i<nwvl*nsa; i++){
		psf=pistat->p[i];
		dfftshift(psf);/*shift peak to center. */
		double pmax=dmax(psf);
		dcog(pgrad,psf,0.5,0.5,0.1*pmax,0.1*pmax);
		dfftshift(psf); /*Shift peak to corner. */
		ccpd(&otf, psf);
		cfft2(otf,-1);
		ctilt(otf,-pgrad[0],-pgrad[1],0);
		cfft2i(otf,1);
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
	}/*if exist */
    }
    //close(fd);
}

void genpistat(const PARMS_S *parms, POWFS_S *powfs){
    double patfov=parms->skyc.patfov;
    double ngsgrid=parms->maos.ngsgrid;
    long ng=ceil(patfov/2/ngsgrid);
    info2("Genpistat..");
    GENPISTAT_S *data=calloc(1, sizeof(GENPISTAT_S));
    data->parms=parms;
    data->powfs=powfs;
    data->ncase=parms->maos.nseed*(2*ng+1)*(2*ng+1)*parms->maos.npowfs;
    data->cases=calloc(4*data->ncase,sizeof(long));
    data->ngsgrid=ngsgrid;
    long count=0;
    for(int iseed=0; iseed<parms->maos.nseed; iseed++){
	int seed=parms->maos.seeds[iseed];/*loop over seed */
	for(long gy=-ng; gy<=ng; gy++){
	    for(long gx=-ng; gx<=ng; gx++){
		/*double thetax=gx*ngsgrid; */
		/*double thetay=gy*ngsgrid; */
		for(int ipowfs=0; ipowfs<parms->maos.npowfs; ipowfs++){/*for ipowfs */
		    data->cases[count][0]=gx;
		    data->cases[count][1]=gy;
		    data->cases[count][2]=ipowfs;
		    data->cases[count][3]=seed;
		    count++;
		    if(count>data->ncase){
			data->ncase=data->ncase*2;
			data->cases=realloc(data->cases,sizeof(long)*data->ncase*4);
		    }
		}
	    }/*for gx */
	}/*for gy */
    }/*for iseed */
    /*info("count=%ld, data->ncase=%ld\n",count,data->ncase); */
    data->ncase=count;
    data->icase=0;
    data->cases=realloc(data->cases, sizeof(long)*4*data->ncase);
    CALL(calc_pistat, data, parms->skyc.nthread,0);
    info2("done\n");
  
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
    simu->bspstrehlxy=dref(gg);

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
			 "%s/pistat/pistat_seed%ld_sa%ld_x%g_y%g",
			 dirstart, seed, msa, thetax, thetay);
		
		if(zfexist(fnpistat)){
		    dcell *tmp=dcellread("%s",fnpistat);
		    for(long ic=0; ic<nsa*nwvl; ic++){
			/*peak is in the corner */
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
    dfree(gg);
}

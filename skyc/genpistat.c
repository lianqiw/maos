/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
    pthread_mutex_t mutex_read;/*don't let them read in the same time. */
    const PARMS_S *parms;
    POWFS_S *powfs;
    double ngsgrid;
    dcell *unwrap;/*matrix that does unwraping. */
}GENPISTAT_S;

static void calc_pistat(GENPISTAT_S *data){
    const PARMS_S *parms=data->parms;
    POWFS_S *powfs=data->powfs;
    int icase=0;
    const int ndtrat=9;
    dmat *dtrats=dnew(ndtrat,1);
    for(int i=0; i<ndtrat; i++){
	dtrats->p[i]=(1<<i);
    }
    mymkdir("%s/pistat", dirstart);
    mymkdir("%s/neaspec", dirstart);
    mymkdir("%s/phygrad", dirstart);
    while(LOCKADD(icase, data->icase, 1)<data->ncase)
#if _OPENMP>=200805
#pragma omp task default(shared) firstprivate(icase)
#endif
    {
	if(icase==0){
	    dwrite(dtrats, "%s/neaspec/neaspec_dtrats", dirstart);
	}
	double thetax=data->ngsgrid*data->cases[icase][0];
	double thetay=data->ngsgrid*data->cases[icase][1];
	long ipowfs=data->cases[icase][2];
	long ncomp=parms->maos.ncomp[ipowfs];
	long seed=data->cases[icase][3];
	long msa=parms->maos.msa[ipowfs];/*in 1-d */
	const long avgstart=100;
	char fnwvf[PATH_MAX], fnpistat[PATH_MAX], fnphygrad[PATH_MAX], fnneaspec[PATH_MAX];
	snprintf(fnwvf,PATH_MAX,"%s/wvfout/wvfout_seed%ld_sa%ld_x%g_y%g",
		 dirstart,seed, msa, thetax, thetay);
	snprintf(fnpistat,PATH_MAX,"%s/pistat/pistat_seed%ld_sa%ld_x%g_y%g",
		 dirstart,seed, msa, thetax, thetay);
	snprintf(fnphygrad,PATH_MAX,"%s/phygrad/phygrad_seed%ld_sa%ld_x%g_y%g",
		 dirstart,seed, msa, thetax, thetay);
	snprintf(fnneaspec,PATH_MAX,"%s/neaspec/neaspec_seed%ld_sa%ld_x%g_y%g",
		 dirstart,seed, msa, thetax, thetay);

	if(zfexist(fnwvf) && (!zfexist(fnpistat) || !zfexist(fnphygrad) || !zfexist(fnneaspec))){
	    dmat *mapply=dnew(2,1);
	    TIC;tic;
	    file_t *fp_wvf=zfopen(fnwvf,"rb");
	    header_t header;
	    read_header(&header, fp_wvf);
	    long nstep=header.nx;
	    if(!iscell(header.magic)){
		error("expected data type: %u, got %u\n", (uint32_t)MCC_ANY, header.magic);
	    }
	    free(header.str); header.str=NULL;
	    const int nsa=msa*msa;
	    const int nwvl=parms->maos.nwvl;
	    dcell *pistat=dcellnew(nsa, nwvl);/*pixel intensity mean(I) */
	    dcell *neaspec=dcellnew(nsa*2, nwvl);
	    dcell **avgpsf=calloc(nwvl, sizeof(dcell*));
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		avgpsf[iwvl]=dcellnew(ndtrat, nsa);
	    }
	    for(long ig=0; ig<2*nsa*nwvl; ig++){
		neaspec->p[ig]=dnew(ndtrat, 1);
	    }
	    dmat *phygrad=dnew(nsa*2, nstep);/*original gradient at each time step. */
	    cmat *wvf=cnew(ncomp,ncomp);
	    cmat *wvfc=NULL;
	    cfft2plan(wvf,-1);
	    PDCELL(pistat, ppistat);
	    PDMAT(phygrad, pphygrad);
	    cmat *otf=cnew(ncomp,ncomp);
	    cfft2plan(otf,1);
	    cfft2plan(otf,-1);
	    dmat *psf=NULL;
	    double nwvli=1./(nwvl);
	    double dtheta[nwvl];
	    const int embfac=parms->maos.embfac[ipowfs];
	    const double dxsa=parms->maos.dxsa[ipowfs];
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		const double wvl=parms->maos.wvl[iwvl];
		dtheta[iwvl]=wvl/(dxsa*embfac);
	    }
	    dmat *gmean=dnew(2,nsa);
	    PDMAT(gmean, pgmean);
	    for(long istep=0; istep<nstep; istep++){
		LOCK(data->mutex_read);
		ccell *wvfi=ccellreaddata(fp_wvf, 0);
		UNLOCK(data->mutex_read);
		PCCELL(wvfi,wvfout);
		for(long iwvl=0; iwvl<nwvl; iwvl++){
		    double wvl=parms->maos.wvl[iwvl];
		    for(long isa=0; isa<nsa; isa++){
			//first compute PSF from WVF and compute CoG
			ccp(&wvfc, wvfout[iwvl][isa]);
			cembed(wvf,wvfc,0,C_FULL);
			cfft2(wvf,-1);
			cabs22d(&psf, 0, wvf, 1);//peak in corner.
			dfftshift(psf);//peak in center
			    
			double grad[2]={0,0};
			double pmax=dmax(psf);
			dcog(grad,psf,0.5, 0.5, 0.*pmax, 0.2*pmax);
			grad[0]*=dtheta[iwvl];//convert to Radian
			grad[1]*=dtheta[iwvl];
			pphygrad[istep][isa]+=grad[0]*nwvli;//record the value
			pphygrad[istep][isa+nsa]+=grad[1]*nwvli;

			if(istep>=avgstart){
			    pgmean[isa][0]+=grad[0];//record the average
			    pgmean[isa][1]+=grad[1];
			    //Then remove the CoG from the WVF and accumulate PSF.
			    mapply->p[0]=-grad[0];
			    mapply->p[1]=-grad[1];
			    int nframe=istep-avgstart+1;
			    for(int idtrat=0; idtrat<ndtrat; idtrat++){
				dmat **pavgpsf=&avgpsf[iwvl]->p[idtrat+isa*ndtrat];
				dadd(pavgpsf, 1, psf, 1);
				if(nframe%(int)dtrats->p[idtrat]==0){
				    grad[0]=grad[1]=0;
				    pmax=dmax(*pavgpsf);
				    dcog(grad,*pavgpsf,0.5, 0.5, 0.*pmax, 0.2*pmax);
				    dzero(*pavgpsf);
				    neaspec->p[isa+iwvl*nsa*2]->p[idtrat]+=pow(grad[0]*dtheta[iwvl],2);
				    neaspec->p[isa+nsa+iwvl*nsa*2]->p[idtrat]+=pow(grad[1]*dtheta[iwvl],2);
				}
			    }
		
			    ngsmod2wvf(wvfc, wvl, mapply, powfs+ipowfs, isa, thetax, thetay, parms);
			    cembed(wvf,wvfc,0,C_FULL);
			    cfft2(wvf,-1);
			    cabs22d(&ppistat[iwvl][isa], 1, wvf, 1);
			}
		    }
		    
		}
		ccellfree(wvfi);
	    }/*for istep */
	    for(int idtrat=0; idtrat<ndtrat; idtrat++){
		int navg=(nstep-avgstart)/dtrats->p[idtrat];
		for(long ig=0; ig<2*nsa*nwvl; ig++){
		    neaspec->p[ig]->p[idtrat]=sqrt(neaspec->p[ig]->p[idtrat]/navg);
		}
	    }
	    zfeof(fp_wvf);
	    zfclose(fp_wvf);
	    cfree(wvf);
	    cfree(wvfc);
	    dfree(mapply);
	    dfree(psf);
	    dscale(gmean, 1./(nwvl*(nstep-avgstart)));
	    dcellscale(pistat, 1./(nstep-avgstart));	

	    //Put back the average gradient to PSF.
	    for(int isa=0; isa<nsa; isa++){
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    int i=isa+nsa*iwvl;
		    psf=pistat->p[i];//peak in corner
		    ccpd(&otf, psf);
		    cfft2(otf,-1);//turn to otf. peak in corner
		    ctilt(otf,pgmean[isa][0]/dtheta[iwvl],pgmean[isa][1]/dtheta[iwvl],0);
		    cfft2i(otf,1);//turn to psf, peak in corner
		    creal2d(&psf,0,otf,1);
		}
	    }
	    dfree(gmean);
	    psf=NULL;
	    /* Saved pistat should have peak on the corner */
	    cfree(otf);
	    dcellwrite(pistat, "%s",fnpistat);
	    dcellwrite(neaspec, "%s",fnneaspec);
	    dwrite(phygrad, "%s",fnphygrad);
	    dcellfree(pistat);
	    dcellfree(neaspec);
	    dcellfreearr(avgpsf, nwvl);
	    dfree(phygrad);
	    toc2("Processing %s:", fnwvf);
	}/*if exist */
    }
#if _OPENMP >= 200805
#pragma omp taskwait
#endif
    dfree(dtrats);
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
		    dcellfree(tmp);
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
#include "mtch.h"
/**Determine WFS nonlinearity.*/
dcell** wfs_nonlinearity(const PARMS_S *parms, POWFS_S *powfs, long seed){
    const int npowfs=parms->maos.npowfs;
    const int nwvl=parms->maos.nwvl;
    double patfov=parms->skyc.patfov;
    const double ngsgrid=parms->maos.ngsgrid;
    dmat *siglevs=dlogspace(1, 2, 20);//dlinspace(10, 5, 20);
    int nstep=1000;
    rand_t rstat;
    seed_rand(&rstat, 1);
    long ng=round(patfov/2/ngsgrid)+1;
    dcell **nonlin=calloc(npowfs, sizeof(dcell*));
    for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
	char fnnonlin[PATH_MAX];
	snprintf(fnnonlin, PATH_MAX, "%s/powfs%d_nonlin", dirstart, ipowfs);
	if(zfexist(fnnonlin)){
	    nonlin[ipowfs]=dcellread("%s", fnnonlin);
	}else{
	    dcell *avgpi=0;
	    const long msa=parms->maos.msa[ipowfs];
	    const long nsa=parms->maos.nsa[ipowfs];
	    const long pixpsa=parms->skyc.pixpsa[ipowfs];
	    const double pixtheta=parms->skyc.pixtheta[ipowfs];
	    long pixpsas[nsa];
	    for(int isa=0; isa<nsa; isa++){
		pixpsas[isa]=pixpsa;
	    }
	    dcell *i0=dcellnew3(nsa, 1, pixpsas, pixpsas);
	    dcell *gx=dcellnew3(nsa, 1, pixpsas, pixpsas);
	    dcell *gy=dcellnew3(nsa, 1, pixpsas, pixpsas);
	    dcell *i0s=dcellnew3(nsa,1, pixpsas, pixpsas);
	    dcell *gxs=dcellnew3(nsa,1, pixpsas, pixpsas);
	    dcell *gys=dcellnew3(nsa,1, pixpsas, pixpsas);
	    dcell *is=dcellnew3(nsa,1, pixpsas, pixpsas);
	    dcell *mtche=0;
	    dmat *sanea=0;
	    ccell *otf1=ccellnew(nsa, nwvl);
	    ccell *otf2=ccellnew(nsa, nwvl);
	    long ncomp=parms->maos.ncomp[ipowfs];
	    for(int isa=0; isa<nsa; isa++){
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    otf1->p[isa+iwvl*nsa]=cnew(ncomp, ncomp);
		    otf2->p[isa+iwvl*nsa]=cnew(ncomp, ncomp);
		    cfft2plan(otf1->p[isa+iwvl*nsa], -1);
		    cfft2plan(otf2->p[isa+iwvl*nsa], 1);
		}
	    }
	    double dtheta[nwvl];
	    const int embfac=parms->maos.embfac[ipowfs];
	    const double dxsa=parms->maos.dxsa[ipowfs];
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		const double wvl=parms->maos.wvl[iwvl];
		dtheta[iwvl]=wvl/(dxsa*embfac);
	    }
	    dcell *nonxy=dcellnew(ng,1);
	    for(int ig=0; ig<ng; ig++){
		nonxy->p[ig]=dnew(siglevs->nx, 2);
		char fnpistat[PATH_MAX];
		int count=0;
		dcellzero(avgpi);
		for(int rx=-1; rx<=1; rx++){
		    for(int ry=-1; ry<=1; ry++){
			if((ig==0 && (abs(rx)+abs(ry))!=0) || (abs(rx)+abs(ry))!=1) continue;
			snprintf(fnpistat, PATH_MAX, "%s/pistat/pistat_seed%ld_sa%ld_x%g_y%g",
				 dirstart,seed, msa, ig*ngsgrid*rx, ig*ngsgrid*ry);
			if(!zfexist(fnpistat)){
			    error("%s doesn't exist\n", fnpistat);
			}else{
			    info2("reading %s\n", fnpistat);
			    dcell *pistat=dcellread("%s", fnpistat);
			    dcelladd(&avgpi, 1, pistat, 1);
			    dcellfree(pistat);
			    count++;
			}
		    }
		}
		if(count>0){
		    PDCELL(avgpi, pavgpi);
		    dcellscale(avgpi, 1./count);
		    dcellzero(i0); dcellzero(gx); dcellzero(gy);
		    for(int isa=0; isa<nsa; isa++){
			/*Assume each WVL has same weighting*/
			for(long iwvl=0; iwvl<nwvl; iwvl++){
			    dcellwrite(avgpi, "avgpi");
			    psf2i0gxgy(i0->p[isa], gx->p[isa], gy->p[isa], pavgpi[iwvl][isa], powfs[ipowfs].dtf+iwvl);
			    ccpd(&otf1->p[isa+iwvl*nsa], pavgpi[iwvl][isa]);
			    cfft2(otf1->p[isa+iwvl*nsa], -1);//turn to otf, peak in corner
			}
		    }
		    /*Build matched filter for different siglevs and test linearity*/
		    for(int isig=0; isig<siglevs->nx; isig++){
			double sig=siglevs->p[isig];
			dcelladd(&i0s, 0, i0, sig);
			dcelladd(&gxs, 0, gx, sig);
			dcelladd(&gys, 0, gy, sig);
			mtch(&mtche, &sanea, i0s, gxs, gys, pixtheta, 3, 0, 0);
			/*dcellwrite(mtche, "mtche_%.0f", sig);
			  dwrite(sanea, "sanea_%.0f", sig);
			  dcellwrite(i0s, "i0s_%.0f", sig);
			  dcellwrite(gxs, "gxs_%.0f", sig);
			  dcellwrite(gys, "gys_%.0f", sig);*/
#define SCAN 0
			double sxe=0,sye=0,neam=0;
			for(int isa=0; isa<nsa; isa++){
			    neam+=sanea->p[isa*2]+sanea->p[isa*2+1];
#if SCAN
			    dmat *resp=dnew(nstep, 4);
#endif
			    for(int istep=0; istep<nstep; istep++){
				dzero(is->p[isa]);
#if SCAN
				double sx=pixtheta*istep/nstep;
				double sy=0;
#else
				double sx=randn(&rstat)*sanea->p[isa*2];
				double sy=randn(&rstat)*sanea->p[isa*2+1];
#endif
				double sout[2]={0,0};
				for(long iwvl=0; iwvl<nwvl; iwvl++){
				    ccp(&otf2->p[isa+iwvl*nsa], otf1->p[isa+iwvl*nsa]);
				    ctilt(otf2->p[isa+iwvl*nsa], sx/dtheta[iwvl], sy/dtheta[iwvl],0);
				    ccwm(otf2->p[isa+iwvl*nsa], powfs[ipowfs].dtf[iwvl].nominal);
				    cfft2i(otf2->p[isa+iwvl*nsa], 1);//turn to psd space
				    spmulcreal(is->p[isa]->p, powfs[ipowfs].dtf[iwvl].si, 
					       otf2->p[isa+iwvl*nsa]->p, sig);
				}
				dmulvec(sout, mtche->p[isa], is->p[isa]->p, 1);
#if SCAN
				resp->p[istep]=sx;
				resp->p[istep+nstep]=sy;
				resp->p[istep+nstep*2]=sout[0];
				resp->p[istep+nstep*3]=sout[1];
#else
				sxe+=(sx-sout[0])*(sx-sout[0]);
				sye+=(sy-sout[1])*(sy-sout[1]);
#endif
			    }
#if SCAN
			    dwrite(resp, "powfs%d_sa%d_sig%g_response", ipowfs, isa, sig);
#endif
			}//for isa
			nonxy->p[ig]->p[isig]=(neam/(nsa*2));
			nonxy->p[ig]->p[isig+siglevs->nx]=sqrt((sxe+sye)/(nstep*2*nsa));
		    }
		}
	    }//for ng
	    dcellwrite(nonxy, "powfs%d_nonxy", ipowfs);
	    dcellfree(avgpi);
	    dcellfree(i0);
	    dcellfree(gx);
	    dcellfree(gy);

	    dcellfree(i0s);
	    dcellfree(gxs);
	    dcellfree(gys);

	    dcellfree(is);
	    dcellfree(mtche);
	    dfree(sanea);
	    ccellfree(otf1);
	    ccellfree(otf2);
	    nonlin[ipowfs]=nonxy;
	    dcellwrite(nonxy, "%s", fnnonlin);
	}
    }//for ipowfs
    dfree(siglevs);
    return nonlin;
}

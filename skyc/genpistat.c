/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
typedef long long4[4];
typedef struct GENPISTAT_S{
	int ncase;
	long4* cases;
	_Atomic(int) icase;
	pthread_mutex_t mutex_read;/*don't let them read in the same time. */
	const PARMS_S* parms;
	POWFS_S* powfs;
	real ngsgrid;
	dcell* unwrap;/*matrix that does unwraping. */
}GENPISTAT_S;

static void calc_pistat(GENPISTAT_S* data){
	const PARMS_S* parms=data->parms;
	POWFS_S* powfs=data->powfs;
	int icase=0;
	const int ndtrat=9;
	dmat* dtrats=dnew(ndtrat, 1);
	for(int i=0; i<ndtrat; i++){
		P(dtrats,i)=(1<<i);
	}
	mymkdir("%s/pistat", dirstart);
	mymkdir("%s/neaspec", dirstart);
	mymkdir("%s/phygrad", dirstart);
	while((icase=atomic_fetch_add(&data->icase, 1))<data->ncase)
#if _OPENMP>=200805
#pragma omp task default(shared) firstprivate(icase)
#endif
	{
		if(icase==0){
			writebin(dtrats, "%s/neaspec/neaspec_dtrats", dirstart);
		}
		real thetax=data->ngsgrid*data->cases[icase][0];
		real thetay=data->ngsgrid*data->cases[icase][1];
		long ipowfs=data->cases[icase][2];
		long ncomp=parms->maos.ncomp[ipowfs];
		long seed=data->cases[icase][3];
		long msa=parms->maos.msa[ipowfs];/*in 1-d */
		const long avgstart=100;
		char fnwvf[PATH_MAX], fnpistat[PATH_MAX], fnphygrad[PATH_MAX], fnneaspec[PATH_MAX];
		snprintf(fnwvf, PATH_MAX, "%s/wvfout/wvfout_seed%ld_sa%ld_x%g_y%g",
			dirstart, seed, msa, thetax, thetay);
		snprintf(fnpistat, PATH_MAX, "%s/pistat/pistat_seed%ld_sa%ld_x%g_y%g",
			dirstart, seed, msa, thetax, thetay);
		snprintf(fnphygrad, PATH_MAX, "%s/phygrad/phygrad_seed%ld_sa%ld_x%g_y%g",
			dirstart, seed, msa, thetax, thetay);
		snprintf(fnneaspec, PATH_MAX, "%s/neaspec/neaspec_seed%ld_sa%ld_x%g_y%g",
			dirstart, seed, msa, thetax, thetay);

		if(zfexist("%s",fnwvf)&&(!zfexist("%s",fnpistat)||!zfexist("%s",fnphygrad)||!zfexist("%s",fnneaspec))){
			dmat* mapply=dnew(2, 1);
			TIC;tic;
			file_t* fp_wvf=zfopen(fnwvf, "rb");
			header_t header={0,0,0,0};
			read_header(&header, fp_wvf);
			long nstep=header.nx;
			if(!iscell(&header.magic)){
				error("expected data type: %u, got %u\n", (uint32_t)MCC_ANY, header.magic);
			}
			free(header.str); header.str=NULL;
			const int nsa=msa*msa;
			const int nwvl=parms->maos.nwvl;
			dcell* pistat=dcellnew(nsa, nwvl);/*pixel intensity mean(I) */
			dcell* neaspec=dcellnew(nsa*2, nwvl);
			dcell** avgpsf=mycalloc(nwvl, dcell*);
			for(int iwvl=0; iwvl<nwvl; iwvl++){
				avgpsf[iwvl]=dcellnew(ndtrat, nsa);
			}
			for(long ig=0; ig<2*nsa*nwvl; ig++){
				P(neaspec,ig)=dnew(ndtrat, 1);
			}
			dmat* phygrad=dnew(nsa*2, nstep);/*original gradient at each time step. */
			cmat* wvf=cnew(ncomp, ncomp);
			cmat* wvfc=NULL;
			//cfft2plan(wvf,-1);
			dcell* ppistat=pistat;
			dmat* pphygrad=phygrad;
			cmat* otf=cnew(ncomp, ncomp);
			//cfft2plan(otf,1);
			//cfft2plan(otf,-1);
			dmat* psf=NULL;
			real nwvli=1./(nwvl);
			real dtheta[nwvl];
			const int embfac=parms->maos.embfac[ipowfs];
			const real dxsa=parms->maos.dxsa[ipowfs];
			for(int iwvl=0; iwvl<nwvl; iwvl++){
				const real wvl=parms->maos.wvl[iwvl];
				dtheta[iwvl]=wvl/(dxsa*embfac);
			}
			dmat* gmean=dnew(2, nsa);
			dmat* pgmean=gmean;
			for(long istep=0; istep<nstep; istep++){
				LOCK(data->mutex_read);
				ccell* wvfi=ccellreaddata(fp_wvf, 0);
				UNLOCK(data->mutex_read);
				if(!wvfi||wvfi->nx==0||wvfi->ny==0) continue;
				for(long iwvl=0; iwvl<nwvl; iwvl++){
					real wvl=parms->maos.wvl[iwvl];
					for(long isa=0; isa<nsa; isa++){
					//first compute PSF from WVF and compute CoG
						ccp(&wvfc, P(wvfi, isa, iwvl));
						cembedc(wvf, wvfc, 0, C_FULL);
						cfft2(wvf, -1);
						cabs22d(&psf, 0, wvf, 1);//peak in corner.
						dfftshift(psf);//peak in center

						real grad[2]={0,0};
						real pmax=dmax(psf);
						dcog(grad, psf, 0.5, 0.5, 0.*pmax, 0.2*pmax, 0);
						grad[0]*=dtheta[iwvl];//convert to Radian
						grad[1]*=dtheta[iwvl];
						P(pphygrad, isa, istep)+=grad[0]*nwvli;//record the value
						P(pphygrad, isa+nsa, istep)+=grad[1]*nwvli;

						if(istep>=avgstart){
							P(pgmean, 0, isa)+=grad[0];//record the average
							P(pgmean, 1, isa)+=grad[1];
							//Then remove the CoG from the WVF and accumulate PSF.
							P(mapply,0)=-grad[0];
							P(mapply,1)=-grad[1];
							int nframe=istep-avgstart+1;
							for(int idtrat=0; idtrat<ndtrat; idtrat++){
								dmat** pavgpsf=&P(avgpsf[iwvl],idtrat,isa);
								dadd(pavgpsf, 1, psf, 1);
								if(nframe%(int)P(dtrats,idtrat)==0){
									grad[0]=grad[1]=0;
									pmax=dmax(*pavgpsf);
									dcog(grad, *pavgpsf, 0.5, 0.5, 0.*pmax, 0.2*pmax, 0);
									dzero(*pavgpsf);
									P(P(neaspec,isa,iwvl),idtrat)+=pow(grad[0]*dtheta[iwvl], 2);
									P(P(neaspec,isa+nsa,iwvl),idtrat)+=pow(grad[1]*dtheta[iwvl], 2);
								}
							}

							ngsmod2wvf(wvfc, wvl, mapply, powfs+ipowfs, isa, thetax, thetay, parms);
							cembedc(wvf, wvfc, 0, C_FULL);
							cfft2(wvf, -1);
							cabs22d(&P(ppistat, isa, iwvl), 1, wvf, 1);
						}
					}

				}
				ccellfree(wvfi);
			}/*for istep */
			for(int idtrat=0; idtrat<ndtrat; idtrat++){
				int navg=(nstep-avgstart)/P(dtrats,idtrat);
				for(long ig=0; ig<2*nsa*nwvl; ig++){
					P(P(neaspec,ig),idtrat)=sqrt(P(P(neaspec,ig),idtrat)/navg);
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
					psf=P(pistat,i);//peak in corner
					ccpd(&otf, psf);
					cfft2(otf, -1);//turn to otf. peak in corner
					ctilt(otf, P(pgmean, 0, isa)/dtheta[iwvl], P(pgmean, 1, isa)/dtheta[iwvl], 0);
					cfft2i(otf, 1);//turn to psf, peak in corner
					creal2d(&psf, 0, otf, 1);
				}
			}
			dfree(gmean);
			psf=NULL;
			/* Saved pistat should have peak on the corner */
			cfree(otf);
			writebin(pistat, "%s", fnpistat);
			writebin(neaspec, "%s", fnneaspec);
			writebin(phygrad, "%s", fnphygrad);
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

void genpistat(const PARMS_S* parms, POWFS_S* powfs){
	real patfov=parms->skyc.patfov;
	real ngsgrid=parms->maos.ngsgrid;
	long ng=ceil(patfov/2/ngsgrid);
	info("Genpistat..");
	GENPISTAT_S* data=mycalloc(1, GENPISTAT_S);
	atomic_init(&data->icase, 0);
	data->parms=parms;
	data->powfs=powfs;
	data->ncase=parms->maos.nseed*(2*ng+1)*(2*ng+1)*parms->maos.npowfs;
	data->cases=mycalloc(data->ncase, long4);
	data->ngsgrid=ngsgrid;
	long count=0;
	for(int iseed=0; iseed<parms->maos.nseed; iseed++){
		int seed=parms->maos.seeds[iseed];/*loop over seed */
		for(long gy=-ng; gy<=ng; gy++){
			for(long gx=-ng; gx<=ng; gx++){
			/*real thetax=gx*ngsgrid; */
			/*real thetay=gy*ngsgrid; */
				for(int ipowfs=0; ipowfs<parms->maos.npowfs; ipowfs++){/*for ipowfs */
					data->cases[count][0]=gx;
					data->cases[count][1]=gy;
					data->cases[count][2]=ipowfs;
					data->cases[count][3]=seed;
					count++;
					if(count>data->ncase){
						data->ncase=data->ncase*2;
						data->cases=myrealloc(data->cases, data->ncase, long4);
					}
				}
			}/*for gx */
		}/*for gy */
	}/*for iseed */
	/*dbg("count=%ld, data->ncase=%ld\n",count,data->ncase); */
	data->ncase=count;
	data->cases=myrealloc(data->cases, data->ncase, long4);
	CALL((thread_wrapfun)calc_pistat, (void*)data, parms->skyc.nthread, 0);
	info("done\n");

	dcellfree(data->unwrap);
	free(data->cases);
	free(data);
}
/**
   Read in the average pixel intensities on the ngs grid
*/
void prep_bspstrehl(SIM_S* simu){
	const PARMS_S* parms=simu->parms;
	real patfov=parms->skyc.patfov;
	real ngsgrid=parms->maos.ngsgrid;
	long ng=ceil(patfov/2/ngsgrid);
	long ng2=ng*2+1;
	long seed=simu->seed_maos;
	dmat* gg=dnew(ng*2+1, 1);
	for(long gx=-ng; gx<=ng; gx++){
		P(gg,gx+ng)=(real)gx;
	}
	simu->bspstrehlxy=dref(gg);

	simu->bspstrehl=mycalloc(parms->maos.npowfs, dcell**);

	long ngnew=ng*10;
	dmat* xnew=dnew(ngnew*2+1, ngnew*2+1);
	dmat* ynew=dnew(ngnew*2+1, ngnew*2+1);
	for(int iy=0; iy<xnew->ny; iy++){
		for(int ix=0; ix<xnew->nx; ix++){
			P(xnew,ix,iy)=(ix-ngnew)*0.1;
			P(ynew,ix,iy)=(iy-ngnew)*0.1;
		}
	}
	for(int ipowfs=0; ipowfs<parms->maos.npowfs; ipowfs++){
		long msa=parms->maos.msa[ipowfs];
		long nsa=parms->maos.nsa[ipowfs];
		long nwvl=parms->maos.nwvl;
		dcell* strehlgrid=dcellnew(nsa, nwvl);
		for(long ic=0; ic<nsa*nwvl; ic++){
			P(strehlgrid,ic)=dnew(ng2, ng2);
		}
		for(long gy=-ng; gy<=ng; gy++){
			real thetay=gy*ngsgrid;
			for(long gx=-ng; gx<=ng; gx++){
				real thetax=gx*ngsgrid;
				char fnpistat[PATH_MAX];
				snprintf(fnpistat, PATH_MAX,
					"%s/pistat/pistat_seed%ld_sa%ld_x%g_y%g",
					dirstart, seed, msa, thetax, thetay);

				if(zfexist("%s",fnpistat)){
					dcell* tmp=dcellread("%s", fnpistat);
					for(long ic=0; ic<nsa*nwvl; ic++){
					/*peak is in the corner */
						P(P(strehlgrid,ic),gx+ng,gy+ng)=P(P(tmp,ic),0);
					}
					dcellfree(tmp);
				}
			}
		}

		simu->bspstrehl[ipowfs]=mycalloc(nsa*nwvl, dcell*);
		writebin(strehlgrid, "strehlgrid_%d", ipowfs);
		for(long ic=0; ic<nsa*nwvl; ic++){
			simu->bspstrehl[ipowfs][ic]=dbspline_prep(gg, gg, P(strehlgrid,ic));
			dmat* strehlnew=dbspline_eval(simu->bspstrehl[ipowfs][ic], gg, gg, xnew, ynew);
			writebin(strehlnew, "strehlnew_%d_%ld", ipowfs, ic);
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
dcell** wfs_nonlinearity(const PARMS_S* parms, POWFS_S* powfs, long seed){
	const int npowfs=parms->maos.npowfs;
	const int nwvl=parms->maos.nwvl;
	real patfov=parms->skyc.patfov;
	const real ngsgrid=parms->maos.ngsgrid;
	dmat* siglevs=dlogspace(1, 2, 20);//dlinspace(10, 5, 20);
	int nstep=1000;
	rand_t rstat;
	seed_rand(&rstat, 1);
	long ng=round(patfov/2/ngsgrid)+1;
	dcell** nonlin=mycalloc(npowfs, dcell*);
	for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
		char fnnonlin[PATH_MAX];
		snprintf(fnnonlin, PATH_MAX, "%s/powfs%d_nonlin", dirstart, ipowfs);
		if(zfexist("%s",fnnonlin)){
			nonlin[ipowfs]=dcellread("%s", fnnonlin);
		} else{
			dcell* avgpi=0;
			const long msa=parms->maos.msa[ipowfs];
			const long nsa=parms->maos.nsa[ipowfs];
			const long pixpsa=parms->skyc.pixpsa[ipowfs];
			const real pixtheta=parms->skyc.pixtheta[ipowfs];
			long pixpsas[nsa];
			for(int isa=0; isa<nsa; isa++){
				pixpsas[isa]=pixpsa;
			}
			dcell* i0=dcellnew3(nsa, 1, pixpsas, pixpsas);
			dcell* gx=dcellnew3(nsa, 1, pixpsas, pixpsas);
			dcell* gy=dcellnew3(nsa, 1, pixpsas, pixpsas);
			dcell* i0s=dcellnew3(nsa, 1, pixpsas, pixpsas);
			dcell* gxs=dcellnew3(nsa, 1, pixpsas, pixpsas);
			dcell* gys=dcellnew3(nsa, 1, pixpsas, pixpsas);
			dcell* is=dcellnew3(nsa, 1, pixpsas, pixpsas);
			dcell* mtche=0;
			dmat* sanea=0;
			ccell* otf1=ccellnew(nsa, nwvl);
			ccell* otf2=ccellnew(nsa, nwvl);
			long ncomp=parms->maos.ncomp[ipowfs];
			for(int isa=0; isa<nsa; isa++){
				for(int iwvl=0; iwvl<nwvl; iwvl++){
					P(otf1,isa,iwvl)=cnew(ncomp, ncomp);
					P(otf2,isa,iwvl)=cnew(ncomp, ncomp);
				}
			}
			real dtheta[nwvl];
			const int embfac=parms->maos.embfac[ipowfs];
			const real dxsa=parms->maos.dxsa[ipowfs];
			for(int iwvl=0; iwvl<nwvl; iwvl++){
				const real wvl=parms->maos.wvl[iwvl];
				dtheta[iwvl]=wvl/(dxsa*embfac);
			}
			dcell* nonxy=dcellnew(ng, 1);
			for(int ig=0; ig<ng; ig++){
				P(nonxy,ig)=dnew(siglevs->nx, 2);
				char fnpistat[PATH_MAX];
				int count=0;
				dcellzero(avgpi);
				for(int rx=-1; rx<=1; rx++){
					for(int ry=-1; ry<=1; ry++){
						if((ig==0&&(abs(rx)+abs(ry))!=0)||(abs(rx)+abs(ry))!=1) continue;
						snprintf(fnpistat, PATH_MAX, "%s/pistat/pistat_seed%ld_sa%ld_x%g_y%g",
							dirstart, seed, msa, ig*ngsgrid*rx, ig*ngsgrid*ry);
						if(!zfexist("%s",fnpistat)){
							error("%s doesn't exist\n", fnpistat);
						} else{
							info("reading %s\n", fnpistat);
							dcell* pistat=dcellread("%s", fnpistat);
							dcelladd(&avgpi, 1, pistat, 1);
							dcellfree(pistat);
							count++;
						}
					}
				}
				if(count>0){
					dcell* pavgpi=avgpi;
					dcellscale(avgpi, 1./count);
					dcellzero(i0); dcellzero(gx); dcellzero(gy);
					for(int isa=0; isa<nsa; isa++){
					/*Assume each WVL has same weighting*/
						for(long iwvl=0; iwvl<nwvl; iwvl++){
							writebin(avgpi, "avgpi");
							psf2i0gxgy(P(i0,isa), P(gx,isa), P(gy,isa), P(pavgpi, isa, iwvl), powfs[ipowfs].dtf+iwvl, 1);
							ccpd(&P(otf1,isa,iwvl), P(pavgpi, isa, iwvl));
							cfft2(P(otf1,isa,iwvl), -1);//turn to otf, peak in corner
						}
					}
					/*Build matched filter for different siglevs and test linearity*/
					for(int isig=0; isig<siglevs->nx; isig++){
						real sig=P(siglevs,isig);
						dcelladd(&i0s, 0, i0, sig);
						dcelladd(&gxs, 0, gx, sig);
						dcelladd(&gys, 0, gy, sig);
						genmtch(&mtche, &sanea, i0s, gxs, gys, pixtheta, 3, 0, 0);
						/*writebin(mtche, "mtche_%.0f", sig);
						  writebin(sanea, "sanea_%.0f", sig);
						  writebin(i0s, "i0s_%.0f", sig);
						  writebin(gxs, "gxs_%.0f", sig);
						  writebin(gys, "gys_%.0f", sig);*/
#define SCAN 0
						real sxe=0, sye=0, neam=0;
						for(int isa=0; isa<nsa; isa++){
							neam+=P(sanea,isa*2)+P(sanea,isa*2+1);
#if SCAN
							dmat* resp=dnew(nstep, 4);
#endif
							for(int istep=0; istep<nstep; istep++){
								dzero(P(is,isa));
#if SCAN
								real sx=pixtheta*istep/nstep;
								real sy=0;
#else
								real sx=randn(&rstat)*P(sanea,isa*2);
								real sy=randn(&rstat)*P(sanea,isa*2+1);
#endif
								real sout[2]={0,0};
								for(long iwvl=0; iwvl<nwvl; iwvl++){
									ccp(&P(otf2,isa,iwvl), P(otf1,isa,iwvl));
									ctilt(P(otf2,isa,iwvl), sx/dtheta[iwvl], sy/dtheta[iwvl], 0);
									ccwm(P(otf2,isa,iwvl), powfs[ipowfs].dtf[iwvl].nominal);
									cfft2i(P(otf2,isa,iwvl), 1);//turn to psd space
									dspmulcreal(P(P(is,isa)), powfs[ipowfs].dtf[iwvl].si,
										P(P(otf2,isa,iwvl)), sig);
								}
								dmulvec(sout, P(mtche,isa), P(P(is,isa)), 1);
#if SCAN
								P(resp,istep, 0)=sx;
								P(resp,istep, 1)=sy;
								P(resp,istep, 2)=sout[0];
								P(resp,istep, 3)=sout[1];
#else
								sxe+=(sx-sout[0])*(sx-sout[0]);
								sye+=(sy-sout[1])*(sy-sout[1]);
#endif
							}
#if SCAN
							writebin(resp, "powfs%d_sa%d_sig%g_response", ipowfs, isa, sig);
#endif
						}//for isa
						P(P(nonxy,ig),isig,0)=(neam/(nsa*2));
						P(P(nonxy,ig),isig,1)=sqrt((sxe+sye)/(nstep*2*nsa));
					}
				}
			}//for ng
			writebin(nonxy, "powfs%d_nonxy", ipowfs);
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
			writebin(nonxy, "%s", fnnonlin);
		}
	}//for ipowfs
	dfree(siglevs);
	return nonlin;
}

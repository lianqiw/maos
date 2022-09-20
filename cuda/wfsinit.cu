/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#include "utils.h"
#include "accphi.h"
#include <curand_kernel.h>
#include "cucmat.h"
#include "wfs.h"
#include "cudata.h"
#include "../maos/pywfs.h"

static void etf2gpu(cucmat& cuetf, etf_t* etf, int icol){
	cmat* etfm=ccell_col(etf->etf, icol);
	cp2gpu(cuetf, etfm);
	cfree(etfm);
}
/**
   Initialize or update etf.
*/
void gpu_wfsgrad_update_etf(const parms_t* parms, const powfs_t* powfs, int ipowfs){
	const int* wfsgpu=cuglobal->wfsgpu();
	TIC;tic;
	if(parms->powfs[ipowfs].usephy&&parms->powfs[ipowfs].llt){
		Array<cuwfs_t>& cuwfs=cuglobal->wfs;
		const int nwvl=parms->powfs[ipowfs].nwvl;
		const int iwfs0=parms->powfs[ipowfs].wfs->p[0];
		for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
			int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
			gpu_set(wfsgpu[iwfs]);/*Only initialize WFS in assigned GPU. */
			for(int iwvl=0; iwvl<nwvl; iwvl++){
				int skip_1=0;
				if(powfs[ipowfs].etfsim2){
					if(fabs(cuwfs[iwfs].dtf[iwvl].etf[1].hs-powfs[ipowfs].etfsim[iwvl].hs)<1
					&&cuwfs[iwfs].dtf[iwvl].etf[1].icol==powfs[ipowfs].etfsim[iwvl].icol){
						//dbg("Reuse etf2 as etf for next step for wfs %d\n", iwfs);
						skip_1=1;
						cucmat tmp=cuwfs[iwfs].dtf[iwvl].etf[0].etf;
						cuwfs[iwfs].dtf[iwvl].etf[0].etf=cuwfs[iwfs].dtf[iwvl].etf[1].etf;
						cuwfs[iwfs].dtf[iwvl].etf[1].etf=tmp;
					}
				}
				if(jwfs==0||parms->powfs[ipowfs].llt->n>1){
					int kwfs=parms->powfs[ipowfs].llt->n>1?jwfs:0;
					if(powfs[ipowfs].etfsim && !skip_1){
						//dbg("Copy etf from wfs %d from CPU to GPU\n", iwfs);
						etf2gpu(cuwfs[iwfs].dtf[iwvl].etf[0].etf, &powfs[ipowfs].etfsim[iwvl], kwfs);
					}
					if(powfs[ipowfs].etfsim2){
						//dbg("Copy etf2 from wfs %d from CPU to GPU\n", iwfs);
						etf2gpu(cuwfs[iwfs].dtf[iwvl].etf[1].etf, &powfs[ipowfs].etfsim2[iwvl], kwfs);
					}
				} else {//copy across GPU
					int iwfs2;
					for(iwfs2=iwfs0; iwfs2<iwfs; iwfs2++){
						if(wfsgpu[iwfs2]==wfsgpu[iwfs]){
							break;
						}
					}
					if(iwfs2==iwfs){
						if(powfs[ipowfs].etfsim&&!skip_1){
							//dbg("Copy etf from wfs %d to %d within GPU\n", iwfs0, iwfs);
							cuwfs[iwfs].dtf[iwvl].etf[0].etf.Copy(cuwfs[iwfs0].dtf[iwvl].etf[0].etf);
						}
						if(powfs[ipowfs].etfsim2){
							//dbg("Copy etf2 from wfs %d to %d within GPU\n", iwfs0, iwfs);
							cuwfs[iwfs].dtf[iwvl].etf[1].etf.Copy(cuwfs[iwfs0].dtf[iwvl].etf[1].etf);
						}
					}else{//reference
						if(powfs[ipowfs].etfsim&&!skip_1){
							//dbg("Reference etf from wfs %d to %d within GPU\n", iwfs2, iwfs);
							cuwfs[iwfs].dtf[iwvl].etf[0].etf=cuwfs[iwfs2].dtf[iwvl].etf[0].etf;
						}
						if(powfs[ipowfs].etfsim2){
							//dbg("Reference etf2 from wfs %d to %d within GPU\n", iwfs2, iwfs);
							cuwfs[iwfs].dtf[iwvl].etf[1].etf=cuwfs[iwfs2].dtf[iwvl].etf[1].etf;
						}
					}
				}
				if(powfs[ipowfs].etfsim){
					cuwfs[iwfs].dtf[iwvl].etf[0].icol=powfs[ipowfs].etfsim[iwvl].icol;
					cuwfs[iwfs].dtf[iwvl].etf[0].hs=powfs[ipowfs].etfsim[iwvl].hs;
				}
				if(powfs[ipowfs].etfsim2){
					cuwfs[iwfs].dtf[iwvl].etf[1].icol=powfs[ipowfs].etfsim2[iwvl].icol;
					cuwfs[iwfs].dtf[iwvl].etf[1].hs=powfs[ipowfs].etfsim2[iwvl].hs;
				}
			}
		}
	}
	toc2("gpu_wfsgrad_update_etf");
}

void gpu_wfsgrad_update_mtche(const parms_t* parms, const powfs_t* powfs, int ipowfs){
	const int* wfsgpu=cuglobal->wfsgpu();
	Array<cuwfs_t>& cuwfs=cuglobal->wfs;
	if(parms->powfs[ipowfs].usephy&&powfs[ipowfs].intstat){
		dbg("powfs%d: updating matched filter in GPU\n", ipowfs);
		const int iwfs0=parms->powfs[ipowfs].wfs->p[0];
		const int nsa=powfs[ipowfs].saloc->nloc;
		const int multi_mf=parms->powfs[ipowfs].phytype_sim==PTYPE_MF&&powfs[ipowfs].intstat->mtche->ny>1;
		for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
			int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
			gpu_set(wfsgpu[iwfs]);/*Only initialize WFS in assigned GPU. */
			int iwfs2;
			for(iwfs2=iwfs0; iwfs2<iwfs; iwfs2++){
				if(wfsgpu[iwfs2]==wfsgpu[iwfs]){
					break;
				}
			}
			if(multi_mf||iwfs2==iwfs){
				if(parms->powfs[ipowfs].phytype_sim==PTYPE_MF){//matched filter
					int icol=multi_mf?jwfs:0;
					dmat* mtche=dcell_col(powfs[ipowfs].intstat->mtche, icol);
					//Delete old referenced values.
					cuwfs[iwfs].mtche.deinit();
					cp2gpu(cuwfs[iwfs].mtche, mtche);
					dfree(mtche);
				}
				if(powfs[ipowfs].intstat->i0sum){
					//Delete old referenced values.
					cuwfs[iwfs].i0sum.deinit();
					cp2gpu(cuwfs[iwfs].i0sum, &PR(powfs[ipowfs].intstat->i0sum, 0, jwfs), nsa, 1);
					cuwfs[iwfs].i0sumsum=PR(powfs[ipowfs].intstat->i0sumsum, jwfs, 0);
				}
			} else{
				cuwfs[iwfs].mtche=cuwfs[iwfs2].mtche;
				cuwfs[iwfs].i0sum=cuwfs[iwfs2].i0sum;
				cuwfs[iwfs].i0sumsum=cuwfs[iwfs2].i0sumsum;
			}
		}
	}
}
/**
   Initialize or update mtched filter
*/

/**
   Initialize other arrays
*/
void gpu_wfsgrad_init(const parms_t* parms, const powfs_t* powfs){
	const int* wfsgpu=cuglobal->wfsgpu();
	cuglobal->wfs.init(parms->nwfs, 1);
	for(int im=0; im<NGPU; im++){
		gpu_set(im);
		cudata->powfs.init(parms->npowfs, 1);
		Array<cupowfs_t>& cupowfs=cudata->powfs;

		/* Setup information that are same for wfs in each powfs*/
		for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
			if(parms->powfs[ipowfs].nwfs==0) continue;
			loc_t* loc=powfs[ipowfs].loc;
			cupowfs[ipowfs].loc=culoc_t(loc);
			if(parms->powfs[ipowfs].type==WFS_PY){//only for pywfs
				cupowfs[ipowfs].saloc=culoc_t(powfs[ipowfs].saloc);
			}
			pts_t* pts=powfs[ipowfs].pts;
			if(pts) cupowfs[ipowfs].pts=cupts_t(pts);

			if(powfs[ipowfs].llt&&parms->powfs[ipowfs].trs){
				pts=powfs[ipowfs].llt->pts;
				loc=powfs[ipowfs].llt->loc;
				cupowfs[ipowfs].llt.pts=cupts_t(pts);
				cupowfs[ipowfs].llt.loc=culoc_t(loc);
			}
			/*cupowfs[ipowfs].skip=parms->powfs[ipowfs].skip; */
			locfft_t* locfft=0;
			if(powfs[ipowfs].pywfs){
				locfft=powfs[ipowfs].pywfs->locfft;
			} else if(powfs[ipowfs].fieldstop){
				locfft=powfs[ipowfs].fieldstop;
			}
			if(locfft){
				const int nwvl=parms->powfs[ipowfs].nwvl;
				cupowfs[ipowfs].embed.init(nwvl, 1);
				cupowfs[ipowfs].nembed.init(nwvl, 1);
				if(locfft->fieldmask){
					cupowfs[ipowfs].fieldstop=curcell(nwvl, 1);
				}
				for(int iwvl=0; iwvl<nwvl; iwvl++){
					cp2gpu(cupowfs[ipowfs].embed[iwvl], P(locfft->embed->p[iwvl]), powfs[ipowfs].loc->nloc, 1);
					cupowfs[ipowfs].nembed[iwvl]=locfft->nembed->p[iwvl];
					if(locfft->fieldmask){
						cp2gpu(cupowfs[ipowfs].fieldstop[iwvl], locfft->fieldmask->p[iwvl]);
					}
				}
			}
			cp2gpu(cupowfs[ipowfs].saa, powfs[ipowfs].saa);
			cp2gpu(cupowfs[ipowfs].pixoffx, powfs[ipowfs].pixoffx);
			cp2gpu(cupowfs[ipowfs].pixoffy, powfs[ipowfs].pixoffy);
			if(powfs[ipowfs].pywfs){
				powfs[ipowfs].pywfs->gpu=1;
				cupowfs[ipowfs].pywfs=powfs[ipowfs].pywfs;
				cp2gpu(cupowfs[ipowfs].pyramid, powfs[ipowfs].pywfs->pyramid);
				cp2gpu(cupowfs[ipowfs].pynominal, powfs[ipowfs].pywfs->nominal);
				cp2gpu(cupowfs[ipowfs].pyoff, powfs[ipowfs].pywfs->gradoff);
				if(powfs[ipowfs].pywfs->msaloc){
					cupowfs[ipowfs].msaloc.init(powfs[ipowfs].pywfs->msaloc->nx, 1);
					for(int i=0; i<powfs[ipowfs].pywfs->msaloc->nx;i++){
						cupowfs[ipowfs].msaloc[i]=culoc_t(powfs[ipowfs].pywfs->msaloc->p[i]);
					}
				}
			}
		}
	}

	/* setup information that maybe different for wfs in same powfs due to
	   misregistration or NCPA.*/
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		const int ipowfs=parms->wfs[iwfs].powfs;
		gpu_set(wfsgpu[iwfs]);/*Only initialize WFS in assigned GPU. */
		Array<cupowfs_t>& cupowfs=cudata->powfs;
		Array<cuwfs_t>& cuwfs=cuglobal->wfs;
		cuwfs[iwfs].stream.reset();//Recreate streams in current GPU.
		cuwfs[iwfs].powfs=&cupowfs[ipowfs];
		const int nsa=powfs[ipowfs].saloc->nloc;
		const int nwvl=parms->powfs[ipowfs].nwvl;
		const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
		const int iwfs0=parms->powfs[ipowfs].wfs->p[0];
		const int nwfsp=parms->powfs[ipowfs].nwfs;
		const int ndm=parms->ndm;
		/*imcc for ztilt. */
		cuwfs[iwfs].loc_dm.init(ndm, 1);
		for(int idm=0; idm<ndm; idm++){
			if(powfs[ipowfs].loc_dm){
				cuwfs[iwfs].loc_dm[idm]=culoc_t(powfs[ipowfs].loc_dm->p[wfsind+idm*nwfsp]);
			} else{
				cuwfs[iwfs].loc_dm[idm]=culoc_t(powfs[ipowfs].loc);
			}
		}
		if(powfs[ipowfs].loc_tel){
			cuwfs[iwfs].loc_tel=culoc_t(powfs[ipowfs].loc_tel->p[wfsind]);
		} else{
			cuwfs[iwfs].loc_tel=culoc_t(powfs[ipowfs].loc);
		}
		cuwfs[iwfs].phiout=curmat(powfs[ipowfs].loc->nloc, 1);
		if(cupowfs[ipowfs].nembed){
			DO(cufftPlan2d(&cuwfs[iwfs].plan_fs, cupowfs[ipowfs].nembed[0], cupowfs[ipowfs].nembed[0], FFT_T_C2C));
			DO(cufftSetStream(cuwfs[iwfs].plan_fs, cuwfs[iwfs].stream));
		}
		if(powfs[ipowfs].saimcc){
			if(powfs[ipowfs].saimcc->nx>1||wfsind==0||wfsgpu[iwfs]!=wfsgpu[iwfs0]){
				int icol=powfs[ipowfs].saimcc->nx>1?wfsind:0;
				dmat* imcc=dcell_col(powfs[ipowfs].saimcc->p[icol], 0);
				cp2gpu(cuwfs[iwfs].imcc, imcc);
				dfree(imcc);
			} else{
				cuwfs[iwfs].imcc=cuwfs[iwfs0].imcc;
			}
		}
		/*GS0 for gtilt. */
		if(powfs[ipowfs].GS0){
			dsp* t=powfs[ipowfs].GS0->p[powfs[ipowfs].GS0->nx>1?wfsind:0];
			cuwfs[iwfs].GS0=cusp(t, 1);
		}
		/*wfs amplitude map on loc */
		cp2gpu(cuwfs[iwfs].amp, powfs[ipowfs].realamp->p[wfsind]);
		if(powfs[ipowfs].neasim){//neasim is the LL' decomposition of nea
			dmat* neasim=PR(powfs[ipowfs].neasim, wfsind, 0);
			if(neasim) cp2gpu(cuwfs[iwfs].neasim, neasim);
		}
		/* * Now start physical optics setup * */
		if(parms->powfs[ipowfs].type==WFS_PY){//Pyramid WFS
			cuwfs[iwfs].pywvf=cuccell(nwvl, 1);
			for(int iwvl=0; iwvl<nwvl; iwvl++){
				cuwfs[iwfs].pywvf[iwvl]=cucmat(cupowfs[ipowfs].nembed[iwvl], cupowfs[ipowfs].nembed[iwvl]);
			}
			const int nxotf=cupowfs[ipowfs].pyramid[0].Nx();;
			const int nyotf=cupowfs[ipowfs].pyramid[0].Ny();
			cuwfs[iwfs].pyotf=cucmat(nxotf, nyotf);
			cuwfs[iwfs].pypsf=curmat(nxotf, nyotf);
			DO(cufftPlan2d(&cuwfs[iwfs].plan_py, nyotf, nxotf, FFT_T_C2C));
			cufftSetStream(cuwfs[iwfs].plan_py, cuwfs[iwfs].stream);
			cuwfs[iwfs].isum=curmat(1, 1);
		} else if(parms->powfs[ipowfs].usephy||parms->powfs[ipowfs].psfout||parms->powfs[ipowfs].pistatout){
			/*If there is llt. */
			if(powfs[ipowfs].llt&&parms->powfs[ipowfs].trs){
				if(powfs[ipowfs].llt->ncpa){
					if(powfs[ipowfs].llt->ncpa->nx>1||wfsind==0||wfsgpu[iwfs]!=wfsgpu[iwfs0]){
						cp2gpu(cuwfs[iwfs].lltncpa,
							powfs[ipowfs].llt->ncpa->p[powfs[ipowfs].llt->ncpa->nx>1?wfsind:0]);
					} else{
						cuwfs[iwfs].lltncpa=cuwfs[iwfs0].lltncpa;
					}
				}
				if(wfsind==0||wfsgpu[iwfs]!=wfsgpu[iwfs0]){
					cp2gpu(cuwfs[iwfs].lltimcc, powfs[ipowfs].llt->imcc->p[0]);
					cp2gpu(cuwfs[iwfs].lltamp, powfs[ipowfs].llt->amp);
				} else{
					cuwfs[iwfs].lltimcc=cuwfs[iwfs0].lltimcc;
					cuwfs[iwfs].lltamp=cuwfs[iwfs0].lltamp;
				}
			}
			/*CUFFTW is row major. */
			int nwvf=powfs[ipowfs].pts->nxsa*parms->powfs[ipowfs].embfac;/*size of fft */
			int nwvf2[2]={nwvf, nwvf};
			const int notfx=powfs[ipowfs].notfx;
			const int notfy=powfs[ipowfs].notfy;
			int notf2[2]={notfy, notfx};
			/*limit the number of subapertures in each batch to less than 1024
			  to save memory. The speed is actually a tiny bit faster for NFIRAOS.*/
			cuwfs[iwfs].msa=(nsa>1024)?((int)ceil((Real)nsa/(Real)(nsa/800))):nsa;
			if(cufftPlanMany(&cuwfs[iwfs].plan_wvf, 2, nwvf2, NULL, 1, 0, NULL, 1, 0,
							 FFT_T_C2C, cuwfs[iwfs].msa)){
				error("CUFFT plan failed\n");
			}
			cufftSetStream(cuwfs[iwfs].plan_wvf, cuwfs[iwfs].stream);

			if(notfx==nwvf&&notfx==notfy){
				cuwfs[iwfs].plan_psf=cuwfs[iwfs].plan_wvf;
			} else{
				if(cufftPlanMany(&cuwfs[iwfs].plan_psf, 2, notf2, NULL, 1, 0, NULL, 1, 0, FFT_T_C2C, cuwfs[iwfs].msa)){
					error("CUFFT plan failed\n");
				}
				cufftSetStream(cuwfs[iwfs].plan_psf, cuwfs[iwfs].stream);
			}
			if(parms->powfs[ipowfs].llt){
				int nlwvf=powfs[ipowfs].llt->pts->nxsa*parms->powfs[ipowfs].embfac;
				int nlwvf2[2]={nlwvf, nlwvf};
				cuwfs[iwfs].lltwvf=cucmat(nlwvf, nlwvf);
				if(cufftPlanMany(&cuwfs[iwfs].lltplan_wvf, 2, nlwvf2, NULL, 1, 0, NULL, 1, 0, FFT_T_C2C, 1)){
					error("CUFFT plan failed\n");
				}
				cufftSetStream(cuwfs[iwfs].lltplan_wvf, cuwfs[iwfs].stream);
				if(nlwvf!=nwvf){
					cuwfs[iwfs].lltotfc2=cucmat(nwvf, nwvf);
					if(cufftPlanMany(&cuwfs[iwfs].lltplan_lotfc2, 2, nwvf2, NULL, 1, 0, NULL, 1, 0, FFT_T_C2C, 1)){
						error("CUFFT plan failed\n");
					}
					cufftSetStream(cuwfs[iwfs].lltplan_lotfc2, cuwfs[iwfs].stream);
				}else{
					cuwfs[iwfs].lltplan_lotfc2=cuwfs[iwfs].lltplan_wvf;
				}
				if(nwvf!=notfx||notfx!=notfy){
					cuwfs[iwfs].lltotfc=cucmat(notfx, notfy);
					if(cufftPlanMany(&cuwfs[iwfs].lltplan_lotfc, 2, notf2, NULL, 1, 0, NULL, 1, 0, FFT_T_C2C, 1)){
						error("CUFFT plan failed\n");
					}
					cufftSetStream(cuwfs[iwfs].lltplan_lotfc, cuwfs[iwfs].stream);
				}else{
					cuwfs[iwfs].lltplan_lotfc=cuwfs[iwfs].lltplan_lotfc2;
				}
			}
			/*DTF. */
			if(parms->powfs[ipowfs].usephy){
				if((parms->powfs[ipowfs].llt&&parms->powfs[ipowfs].llt->n>1)||wfsind==0||wfsgpu[iwfs]!=wfsgpu[iwfs0]){
					 /*Need one per wfs in this powfs, or the first wfs. */
					cuwfs[iwfs].dtf.init(nwvl, 1);
					for(int iwvl=0; iwvl<nwvl; iwvl++){
						if(!parms->powfs[ipowfs].llt){
							int icol=powfs[ipowfs].dtf[iwvl].nominal->ny>1?wfsind:0;
							cmat* nominal=ccell_col(powfs[ipowfs].dtf[iwvl].nominal, icol);
							cp2gpu(cuwfs[iwfs].dtf[iwvl].nominal, nominal);
							cfree(nominal);
						}
						//ETF moved to gpu_wfsgrad_update_etf();
					}/*for iwvl. */
					if(parms->powfs[ipowfs].llt){
						cp2gpu(cuwfs[iwfs].srot, PR(powfs[ipowfs].srot, wfsind, 0));
					}
				} else{
					cuwfs[iwfs].dtf=cuwfs[iwfs0].dtf;
					cuwfs[iwfs].srot=cuwfs[iwfs0].srot;
				}
				if(wfsind==0||wfsgpu[iwfs]!=wfsgpu[iwfs0]){
					cp2gpu(cuwfs[iwfs].qe, parms->powfs[ipowfs].qe);
				} else{
					cuwfs[iwfs].qe=cuwfs[iwfs0].qe;
				}
				//Matched filter is separated with gpu_wfsgrad_upate_mtche();
			
				if(powfs[ipowfs].bkgrnd){
					if(powfs[ipowfs].bkgrnd->ny==1||wfsind==0||wfsgpu[iwfs]!=wfsgpu[iwfs0]){
						int icol=(powfs[ipowfs].bkgrnd->ny==1?wfsind:0);
						dmat* bkgrnd=dcell_col(powfs[ipowfs].bkgrnd, icol);
						cp2gpu(cuwfs[iwfs].bkgrnd2, bkgrnd);
						dfree(bkgrnd);
					} else{
						cuwfs[iwfs].bkgrnd2=cuwfs[iwfs0].bkgrnd2;
					}
				}
				if(powfs[ipowfs].bkgrndc){
					if(powfs[ipowfs].bkgrndc->ny==1||wfsind==0||wfsgpu[iwfs]!=wfsgpu[iwfs0]){
						int icol=(powfs[ipowfs].bkgrndc->ny==1?wfsind:0);
						dmat* bkgrnd=dcell_col(powfs[ipowfs].bkgrndc, icol);
						cp2gpu(cuwfs[iwfs].bkgrnd2c, bkgrnd);
						dfree(bkgrnd);

					} else{
						cuwfs[iwfs].bkgrnd2c=cuwfs[iwfs0].bkgrnd2c;
					}
				}
				if(parms->powfs[ipowfs].dither==1){
					cuwfs[iwfs].dither=Dither_t(nsa, powfs[ipowfs].pixpsax, powfs[ipowfs].pixpsay, 1);
				} else if(parms->powfs[ipowfs].dither==-1){
					cuwfs[iwfs].dither=Dither_t(nsa, powfs[ipowfs].pixpsax, powfs[ipowfs].pixpsay, 0);
				}
			}
			const int msa=cuwfs[iwfs].msa;
			cuwfs[iwfs].wvf=cucmat(nwvf*nwvf, msa);
			if(nwvf!=notfx||notfx!=notfy){
				cuwfs[iwfs].psf=cucmat(notfx*notfy, msa);
			}
			/*if(notfx!=notf||notfy!=notf){
				cuwfs[iwfs].otf=cucmat(notfx*notfy, msa);//no longer used
			}*/
			if(parms->powfs[ipowfs].psfout){
				cuwfs[iwfs].wvfout=cuccell(nsa, nwvl, notfx/2+2, notfy/2+2);
				cuwfs[iwfs].psfout=cucmat(notfx*notfy, msa);
			}
			if(parms->powfs[ipowfs].pistatout){
				cuwfs[iwfs].psfstat=cucmat(notfx*notfy, msa);
			}
			if(powfs[ipowfs].llt&&parms->powfs[ipowfs].trs){
				int nlx=powfs[ipowfs].llt->pts->nxsa;
				int nlwvf=nlx*parms->powfs[ipowfs].embfac;
				cuwfs[iwfs].lltopd=curmat(nlx, nlx);
				if(parms->powfs[ipowfs].pistatout||parms->sim.idealfsm){
					cuwfs[iwfs].lltg.init(2, 1);
				}
				if(nlwvf!=nwvf){
					cuwfs[iwfs].lltotfc2=cucmat(nwvf, nwvf);
				}
				if(nwvf!=notfx||notfx!=notfy){
					cuwfs[iwfs].lltotfc=cucmat(notfx, notfy);
				}
			}
		}/*if phy */
		CUDA_SYNC_DEVICE;
	}/*for iwfs */
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		gpu_wfsgrad_update_etf(parms, powfs, ipowfs);
		gpu_wfsgrad_update_mtche(parms, powfs, ipowfs);
	}
	gpu_print_mem("wfs init");
}
void gpu_wfs_init_sim(const parms_t* parms, powfs_t* powfs){
	int* wfsgpu=cuglobal->wfsgpu();
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		gpu_set(wfsgpu[iwfs]);/*Only initialize WFS in assigned GPU. */
		Array<cuwfs_t>& cuwfs=cuglobal->wfs;
		int ipowfs=parms->wfs[iwfs].powfs;
		int nsa=powfs[ipowfs].saloc->nloc;
		//gradacc is used for accumulation in geom mode and for output in phy mode
		int ng=parms->powfs[ipowfs].ng;
		initzero(cuwfs[iwfs].gradacc, nsa*ng, 1);
		initzero(cuwfs[iwfs].gradcalc, nsa*ng, 1);
		if(parms->powfs[ipowfs].usephy||parms->powfs[ipowfs].dither){
			if(!cuwfs[iwfs].ints){
				if(parms->powfs[ipowfs].type==WFS_PY){//PYWFS
					cuwfs[iwfs].ints=curcell(1, 1, nsa, powfs[ipowfs].pywfs->nside);
				} else{
					cuwfs[iwfs].ints=curcell(nsa, 1, powfs[ipowfs].pixpsax, powfs[ipowfs].pixpsay);
				}
			} else{
				cuzero(cuwfs[iwfs].ints);
			}
		}
		if(parms->powfs[ipowfs].pistatout){
			if(parms->powfs[ipowfs].pistatstc){
				error("pistatstc is not supported yet.\n");
			}
			if(!cuwfs[iwfs].pistatout){
				const int notfx=powfs[ipowfs].notfx;/*necessary size to build detector image. */
				const int notfy=powfs[ipowfs].notfy;
				const int npsf=MAX(notfx, notfy);
				cuwfs[iwfs].pistatout=curcell(nsa, parms->powfs[ipowfs].nwvl, npsf, npsf);
			} else{
				cuzero(cuwfs[iwfs].pistatout);
			}
		}
		if(parms->powfs[ipowfs].i0save){
			cuzero(cuwfs[iwfs].intsout);
		}
		CUDA_SYNC_DEVICE;
	}
}
void gpu_wfssurf2gpu(const parms_t* parms, powfs_t* powfs){
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		gpu_set(cuglobal->wfsgpu[iwfs]);
		Array<cuwfs_t>& cuwfs=cuglobal->wfs;
		int ipowfs=parms->wfs[iwfs].powfs;
		int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
		if(powfs[ipowfs].opdadd&&powfs[ipowfs].opdadd->p[wfsind]){
			cp2gpu(cuwfs[iwfs].opdadd, powfs[ipowfs].opdadd->p[wfsind]);
		}
	}
}
__global__ static void setup_rand(curandState* rstat, int seed){
	int id=threadIdx.x+blockIdx.x*blockDim.x;
	curand_init(seed, id, 0, &rstat[id]);
}
/**
   Seed the random number genrator
*/
void gpu_wfsgrad_seeding(const parms_t* parms, const powfs_t* powfs, rand_t* rstat){
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		gpu_set(cuglobal->wfsgpu[iwfs]);
		Array<cuwfs_t>& cuwfs=cuglobal->wfs;
		int seed=lrand(rstat);/*don't put this after continue. */
		int ipowfs=parms->wfs[iwfs].powfs;
		if(!parms->powfs[ipowfs].noisy) continue;
		int nsa=powfs[ipowfs].saloc->nloc*2;
		if(nsa<RAND_THREAD){
			cuwfs[iwfs].custatt=nsa;//number of threads
			cuwfs[iwfs].custatb=1;//number of blocks
		} else if(nsa<RAND_THREAD*RAND_BLOCK){
			cuwfs[iwfs].custatt=RAND_THREAD;
			cuwfs[iwfs].custatb=nsa/RAND_THREAD+(nsa%RAND_THREAD)?1:0;
		} else{
			cuwfs[iwfs].custatt=RAND_THREAD;
			cuwfs[iwfs].custatb=RAND_BLOCK;
		}
		DO(cudaMalloc(&cuwfs[iwfs].custat, (cuwfs[iwfs].custatt*cuwfs[iwfs].custatb)*sizeof(curandState)));
		setup_rand<<<cuwfs[iwfs].custatb, cuwfs[iwfs].custatt>>>(cuwfs[iwfs].custat, seed);
	}
	CUDA_SYNC_DEVICE;
	gpu_print_mem("wfs seeding");
}

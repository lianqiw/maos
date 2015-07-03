/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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




#include "common.h"
static TIC;
#include "setup_aper.h"
/**
   \file setup_aper.c
   Contains routines that setup the telescope aperture.
*/


/**
   Setting up aperture cordinate grid aper_locs, and amplitude map for
performance evaluation. */
APER_T * setup_aper(const PARMS_T *const parms){
    APER_T *aper = calloc(1, sizeof(APER_T));
    tic;
    if(parms->aper.fnamp){
	info2("Reading aperture amplitude map from %s\n", parms->aper.fnamp);
	aper->ampground=mapread("%s",parms->aper.fnamp);
	if(fabs(aper->ampground->h)>1.e-14){
	    warning("ampground is not on ground, this is not tested\n");
	}else{
	    double amp_d, amp_din;
	    map_d_din(aper->ampground, &amp_d, &amp_din);
	    if(fabs(parms->aper.d - amp_d) > 1 ||
	       fabs(parms->aper.din - amp_din) > 0.5){
		if(!parms->aper.fnampuser){
		    warning2("Amplitude map does not match aperture diameter: amp.d=(%g, %g) aper.d=(%g, %g). Disabled\n", 
			     amp_d, amp_din, parms->aper.d, parms->aper.din);
		    mapfree(aper->ampground);
		    free(((PARMS_T*)parms)->aper.fnamp);
		    ((PARMS_T*)parms)->aper.fnamp=0;
		}else{
		    error("Use supplied amplitude map does not match aperture diameter: amp.d=(%g, %g) aper.d=(%g, %g).\n", 
			  amp_d, amp_din, parms->aper.d, parms->aper.din);
		}
	    }
	}
	if(fabs(parms->aper.rotdeg)>1.e-12){
	    warning("Pupil is rotated by %g deg\n",parms->aper.rotdeg);
	    const long nx=aper->ampground->nx;
	    const long ny=aper->ampground->ny;
	    dmat *B=dnew_data(nx, ny, aper->ampground->p);
	    aper->ampground->p=calloc(nx*ny, sizeof(double));
	    dembed((dmat*)aper->ampground,B,parms->aper.rotdeg/180.*M_PI);
	    dfree(B);
	}
    }
    if(parms->load.locs){
	char *fn=parms->load.locs;
	if(!zfexist(fn)) error("%s doesn't exist\n",fn);
	warning("Loading plocs from %s\n",fn);
	aper->locs=locread("%s",fn);
    }else{/* locs act as a pupil mask. no points outside will be evaluated. */
	map_t *smap=0;
	create_metapupil(&smap,0,0,parms->dirs,parms->aper.d,0,parms->evl.dx,parms->evl.dx,0,0.5,0,0,0,0);
	aper->locs=map2loc(smap, 0);
	mapfree(smap);
    }
    if(!aper->amp){
	/*Use negative misreg.pupil because it is the misreg of telescope entrace, not our pupil.*/
	aper->amp=mkamp(aper->locs, aper->ampground, 
			-parms->misreg.pupil->p[0], -parms->misreg.pupil->p[1], 
			parms->aper.d, parms->aper.din);
    }
    if(parms->aper.pupmask){
	map_t *mask=mapread("%s",parms->aper.pupmask);
	if(fabs(parms->aper.rotdeg)>1.e-12){
	    warning("Pupil mask is rotated by %g deg\n",parms->aper.rotdeg);
	    dmat *B=dnew_data(mask->nx, mask->ny, mask->p);
	    mask->p=calloc(mask->nx*mask->ny, sizeof(double));
	    dembed((dmat*)mask, B, parms->aper.rotdeg/180.*M_PI);
	    dfree(B);
	}
	dmat *ampmask=dnew(aper->locs->nloc, 1);
	prop_grid_stat(mask, aper->locs->stat, ampmask->p, 1, 0, 0, 1, 0, 0, 0);
	dcwm(aper->amp, ampmask);
	dfree(ampmask);
	mapfree(mask);
    }else{/*apply an annular mask */
	//locannularmask(aper->amp->p, aper->locs, 0,0, parms->aper.d*0.5, parms->aper.din*0.5);
    }
    if(!parms->load.locs){
	loc_reduce(aper->locs, aper->amp, EPS, 1, NULL);
    }
    loc_create_stat(aper->locs);
    if(parms->misreg.dm2sci){
	int isset=0;
	int nevl=parms->evl.nevl;
	aper->locs_dm=cellnew(parms->ndm, nevl);
	for(int idm=0; idm<parms->ndm; idm++){
	    for(int ievl=0; ievl<nevl; ievl++){
		if(parms->misreg.dm2sci[ievl+idm*nevl])
#pragma omp task shared(isset)
		{
		    aper->locs_dm->p[ievl+idm*nevl]
			=loctransform(aper->locs, parms->misreg.dm2sci[ievl+idm*nevl]);
		    isset=1;
		}
	    }
	}
#pragma omp taskwait
	if(!isset){
	    cellfree(aper->locs_dm);
	}else if(parms->save.setup){
	    writebin(aper->locs_dm, "locs_dm");
	}
    }
    /*Set the amp for plotting. */
    aper->amp1=ddup(aper->amp);
    /*normalize amp to sum to 1. */
    normalize_sum(aper->amp->p, aper->locs->nloc, 1);
    aper->sumamp2=dsumsq(aper->amp);
    aper->mcc=loc_mcc_ptt(aper->locs, aper->amp->p);
    aper->ipcc=1./aper->mcc->p[0];/*piston inverse. should be 1 since amp is normlaized. */
    aper->imcc=dinvspd(aper->mcc);/*pttr inverse */
    /*piston term correction in focus mode */
    aper->fcp=(aper->mcc->p[aper->mcc->nx+1]+aper->mcc->p[2*(aper->mcc->nx+1)])*aper->ipcc;
    if(parms->evl.rmax!=1){
	aper->mod=zernike(aper->locs, parms->aper.d, 0, parms->evl.rmax, 0);
	if(parms->save.setup){
	    writebin(aper->mod,"aper_mode");
	}
	dgramschmidt(aper->mod, aper->amp->p);
	if(parms->save.setup){
	    writebin(aper->mod,"aper_mode_gramschmidt");
	}
    }
    if(parms->save.setup){
	locwrite(aper->locs, "aper_locs");
	writebin(aper->amp, "aper_amp");
	writebin(aper->mcc, "aper_mcc");
    }

    if(parms->evl.psfmean || parms->evl.psfhist){
	aper->embed=locfft_init(aper->locs, aper->amp, parms->evl.wvl, 
				parms->evl.psfgridsize, 2, 0);
	for(int iwvl=0; iwvl<parms->evl.nwvl; iwvl++){
	    long nembed=aper->embed->nembed->p[iwvl];
	    if(parms->evl.psfsize->p[iwvl]<1 || parms->evl.psfsize->p[iwvl] > nembed){
		parms->evl.psfsize->p[iwvl] = nembed;
	    }
	    info2("iwvl %d: Science PSF is using grid size of %ld. The PSF will sum to %.15g\n",
		  iwvl, nembed, aper->sumamp2*nembed*nembed);
	}
    }
    toc2("setup_aper");
    return aper;
}
/**
   Free the aper structure after simulation*/
void free_aper(const PARMS_T *parms, APER_T *aper){
    /*aper->ampground is freed on setup_recon*/
    locfree(aper->locs);
    dfree(aper->amp);
    dfree(aper->amp1);
    dfree(aper->imcc);
    dfree(aper->mcc);
    dfree(aper->mod);
    locfft_free(aper->embed);
    dcellfree(aper->opdadd);
    mapfree(aper->ampground);
    free(aper);
}

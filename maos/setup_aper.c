/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "maos.h"
TIC;
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
	aper->ampground=mapread("%s",parms->aper.fnamp);
	if(fabs(aper->ampground->h)>1.e-14){
	    warning("ampground is not on ground, this is not tested\n");
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
	if(!exist(fn)) error("%s doesn't exist\n",fn);
	warning("Loading plocs from %s\n",fn);
	aper->locs=locread("%s",fn);
    }else{//locs act as a pupil mask. no points outside will be evaluated.
	/*
	  Commented out on 2011-02-14
	  if(aper->ampground && fabs(aper->ampground->dx-dx)<1.e-6 && !parms->evl.ismisreg){
	  //LOCSTAT records the starting of each row to speed up accphi
	  locstat_t *locstat=calloc(1, sizeof(locstat_t));
	  info2("Using amplitude map to generate aper grid\n");
	  aper->locs=mkcirloc_amp(&(aper->amp), locstat, aper->ampground, 
	  d, dx,parms->aper.cropamp);
	  aper->locs->stat=locstat;
	  }else{
	
	  map_t *pmap=create_metapupil_wrap(parms,0,dx,0,0,0,T_PLOC,0,0);
	  aper->locs=map2loc(pmap);
	  mapfree(pmap);
	  warning("No amplitude map defined or matched to aperture dx.\n");
	  }  */
	//We choose to make circular map so that it acts as pupil stop in case of misregistration.
	aper->locs=mkcirloc(parms->aper.d, parms->aper.dx);
	if(parms->save.setup){
	    locwrite(aper->locs, "%s/aper_locs.bin.gz",dirsetup);
	}
    }
    loc_create_stat(aper->locs);
    if(!aper->amp){
	aper->amp=dnew(aper->locs->nloc, 1);
	if(aper->ampground){
	    prop_grid_stat(aper->ampground, aper->locs->stat, aper->amp->p, 1,
			   parms->evl.misreg[0]-parms->aper.misreg[0],
			   parms->evl.misreg[1]-parms->aper.misreg[1],
			   1, 0, 0, 0);
	}else{
	    warning("Using locannular to create a gray pixel aperture\n");
	    locannular(aper->amp->p, aper->locs,
		       parms->aper.misreg[0]-parms->evl.misreg[0],
		       parms->aper.misreg[1]-parms->evl.misreg[1],
		       parms->aper.d*0.5,parms->aper.din*0.5,1);
	}
    }
    if(parms->aper.pupmask){
	map_t *mask=mapread("%s",parms->aper.pupmask);
	dmat *ampmask=dnew(aper->locs->nloc, 1);
	prop_grid_stat(mask, aper->locs->stat, ampmask->p, 1, 0, 0, 1, 0, 0, 0);
	dcwm(aper->amp, ampmask);
	dfree(ampmask);
	mapfree(mask);
    }else{//apply an annular mask
	locannularmask(aper->amp->p, aper->locs, 0,0, parms->aper.d*0.5, parms->aper.din*0.5);
    }
    //Set the amp for plotting.
    if(parms->plot.setup ||parms->plot.run){
	aper->amp1=ddup(aper->amp);
    }
    //normalize amp to sum to 1.
    normalize(aper->amp->p, aper->locs->nloc, 1);
    if(parms->plot.setup){
	drawopd("amp",aper->locs,aper->amp1->p,NULL,"Aperture Amplitude Map",
		"x (m)","y (m)","aper");
    }
    aper->mcc=loc_mcc_ptt(aper->locs, aper->amp->p);
    aper->ipcc=1./aper->mcc->p[0];//piston inverse. should be 1 since amp is normlaized.
    aper->imcc=dinvspd(aper->mcc);//pttr inverse
    //piston term correction in focus mode
    aper->fcp=(aper->mcc->p[aper->mcc->nx+1]+aper->mcc->p[2*(aper->mcc->nx+1)])*aper->ipcc;
    if(parms->evl.rmax!=1){
	aper->mod=loc_zernike(aper->locs, parms->aper.d/2, parms->evl.rmax);
	if(parms->save.setup){
	    dwrite(aper->mod,"%s/aper_mode.bin.gz",dirsetup);
	}
	dgramschmidt(aper->mod, aper->amp->p);
	if(parms->save.setup){
	    dwrite(aper->mod,"%s/aper_mode_gramschmidt.bin.gz",dirsetup);
	}
    }
    if(parms->save.setup){
	dwrite(aper->amp, "%s/aper_amp.bin.gz",dirsetup);
	dwrite(aper->mcc, "%s/aper_mcc.bin.gz",dirsetup);
    }
    aper->sumamp2=0;
    if(aper->amp){
	for(long iloc=0; iloc<aper->locs->nloc; iloc++){
	    aper->sumamp2+=pow(aper->amp->p[iloc],2);
	}
    }
    if(parms->evl.psfmean || parms->evl.psfhist){
	const int nwvl=parms->evl.nwvl;
	aper->nembed=calloc(nwvl, sizeof(int));
	aper->embed=calloc(nwvl, sizeof(int*));
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    aper->nembed[iwvl]=parms->evl.psfgridsize[iwvl];
	    aper->embed[iwvl]=loc_create_embed(&(aper->nembed[iwvl]), aper->locs);
	    if(parms->evl.psfsize[iwvl]<1 || parms->evl.psfsize[iwvl] > aper->nembed[iwvl]){
		parms->evl.psfsize[iwvl] = aper->nembed[iwvl];
	    }
	    info2("iwvl %d: Science PSF is using grid size of %d. The PSF will sum to %.15g\n",
		  iwvl, aper->nembed[iwvl], aper->sumamp2*aper->nembed[iwvl]*aper->nembed[iwvl]);
	}
    }
    toc2("setup_aper");
    return aper;
}
/**
   Free the aper structure after simulation*/
void free_aper(APER_T *aper, const PARMS_T *parms){
    /*aper->ampground is freed on setup_recon*/
    locfree(aper->locs);
    dfree(aper->amp);
    dfree(aper->amp1);
    dfree(aper->imcc);
    dfree(aper->mcc);
    dfree(aper->mod);
    if(aper->embed){
	for(int iwvl=0; iwvl<parms->evl.nwvl; iwvl++){
	    free(aper->embed[iwvl]);
	}
	free(aper->embed);
	free(aper->nembed);
    }
    free(aper);
}

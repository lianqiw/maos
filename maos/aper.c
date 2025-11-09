/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
/**
   \file aper.h

   Contains routines that setup the telescope aperture for performance evaluation.
*/


#include "common.h"
static TIC;
#include "aper.h"

/**
   Setting up aperture cordinate grid aper_locs, and amplitude map for
performance evaluation. */
aper_t* setup_aper(const parms_t* const parms){
	aper_t* aper=mycalloc(1, aper_t);
	tic;

	if(parms->load.locs){
		char* fn=parms->load.locs;
		if(!zfexist("%s",fn)) error("%s doesn't exist\n", fn);
		if(!(aper->locs=locread("%s", fn))){
			error("Failed to Load locs from %s\n", fn);
		}else{
			warning("Loaded locs from %s\n", fn);
		}
		if(fabs(aper->locs->dx-parms->evl.dx)>1e-7*parms->evl.dx){
			warning("Loaded locs has unexpected sampling of %g, should be %g\n",
				aper->locs->dx, parms->evl.dx);
		}
	} else{/* locs act as a pupil mask. no points outside will be evaluated. */
		if(parms->aper.amp && fabs(parms->evl.dx-parms->aper.amp->dx)<EPS){
			aper->locs=loc_from_map(parms->aper.amp, 0);
		}else{
			map_t* smap=0;
			create_metapupil(&smap, 0, 0, parms->dirs, parms->aper.d, 0, parms->evl.dx, parms->evl.dx, 0, parms->evl.dx, 0, 0, 0, 0);
			aper->locs=loc_from_map(smap, 0);
			mapfree(smap);
		}
	}
	if(!aper->amp){
		/*Use negative misreg.pupil because it is the misreg of telescope entrace, not our pupil.*/
		aper->amp=mkamp(aper->locs, parms->aper.amp, 
			-P(parms->aper.misreg, 0), -P(parms->aper.misreg, 1),
			parms->aper.d, parms->aper.din);
	}
	if(parms->aper.pupmask){
		map_t* mask=mapread("%s", parms->aper.pupmask);
		if(!mask){
			error("Failed to load pupil mask from %s\n", parms->aper.pupmask);
		}else{
			//aper.miserg is between the telescope pupil and the pupil mask in the instrument.
			info("Loaded pupil mask from %s\n", parms->aper.pupmask);
		 	if(parms->aper.rot){
				dbg("Pupil mask is rotated by %g deg\n", parms->aper.rot*180./M_PI);
				dmaprot(mask, parms->aper.rot);
			}
			dmat* ampmask=dnew(aper->locs->nloc, 1);
			prop(&(propdata_t){.mapin=mask, .locout=aper->locs, .phiout=P(ampmask)}, 0, 0);
			dcwm(aper->amp, ampmask);
			dfree(ampmask);
			mapfree(mask);
		}
	}
	if(!parms->load.locs){
		loc_reduce(aper->locs, aper->amp, EPS, 1, NULL);
	}
	loc_create_stat(aper->locs);
	if(parms->distortion.dm2sci){
		int isset=0;
		int nevl=parms->evl.nevl;
		aper->locs_dm=loccellnew(nevl, parms->ndm);
		for(int idm=0; idm<parms->ndm; idm++){
			for(int ievl=0; ievl<nevl; ievl++){
				if(parms->distortion.dm2sci[ievl+idm*nevl])
				{
					P(aper->locs_dm, ievl, idm)
						=loctransform(aper->locs, parms->distortion.dm2sci[ievl+idm*nevl]);
					isset=1;
				}
			}
		}
		if(!isset){
			cellfree(aper->locs_dm);
		} else if(parms->save.setup){
			writebin(aper->locs_dm, "aper_locs_dm");
		}
	}
	/*Set the amp for plotting. */
	aper->amp1=ddup(aper->amp);
	/*normalize amp to sum to 1. */
	dnormalize_sum(aper->amp, 1);
	aper->sumamp2=dsumsq(aper->amp);
	aper->mcc=loc_mcc_ptt(aper->locs, P(aper->amp));
	aper->ipcc=1./P(aper->mcc,0);/*piston inverse. should be 1 since amp is normlaized. */
	aper->imcc=dinvspd(aper->mcc);/*pttr inverse */
	/*piston term correction in focus mode */
	aper->fcp=(P(aper->mcc,1,1)+P(aper->mcc,2,2))*aper->ipcc;
	if(parms->evl.rmax!=1){
		aper->mod=zernike(aper->locs, parms->aper.d, 0, parms->evl.rmax, 0);
		if(parms->save.setup){
			writebin(aper->mod, "aper_mode");
		}
		dgramschmidt(aper->mod, P(aper->amp));
		if(parms->save.setup){
			writebin(aper->mod, "aper_mode_gramschmidt");
		}
	}
	if(parms->save.setup){
		locwrite(aper->locs, "aper_locs");
		writebin(aper->amp, "aper_amp");
		writebin(aper->mcc, "aper_mcc");
	}

	if(parms->evl.psfmean||parms->evl.psfhist){
		aper->locfft=locfft_init(aper->locs, aper->amp, parms->evl.wvl,
			parms->evl.psfgridsize, 2, 0);
		long nembed_old=0;
		for(int iwvl=0; iwvl<parms->evl.nwvl; iwvl++){
			long nembed=P(aper->locfft->nembed,iwvl);
			if(P(parms->evl.psfsize,iwvl)<1||P(parms->evl.psfsize,iwvl) > nembed){
				P(parms->evl.psfsize,iwvl)=nembed;
			}
			if(nembed!=nembed_old){
				nembed_old=nembed;
				info2("iwvl %d: Science PSF is using grid size of %ld. The PSF will sum to %.15g\n",
					iwvl, nembed, aper->sumamp2*nembed*nembed);
			}
		}
	}
	toc2("setup_aper");
	return aper;
}
/**
   Free the aper structure after simulation*/
void free_aper(aper_t* aper){
	/*parms->aper.amp is freed on setup_recon*/
	locfree(aper->locs);
	cellfree(aper->locs_dm);
	dfree(aper->amp);
	dfree(aper->amp1);
	dfree(aper->imcc);
	dfree(aper->mcc);
	dfree(aper->mod);
	locfft_free(aper->locfft);
	dcellfree(aper->opdadd);
	free(aper);
}

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
#include "../math/mathdef.h"
#include "../math/mathdef.h"
#include "slaving.h"

/**
   Compute the actuator coupling coefficient to be used to identify non-coupled
   actuators. W1 is optional weighting function. The max value is 1.
*/
dcell *genactcpl(const spcell *HA, const dmat *W1){
    int ndm=HA->ny;
    dcell *actcplc=dcellnew(ndm, 1);
    for(int idm=0; idm<ndm; idm++){
	for(int ifit=0; ifit<HA->nx; ifit++){
	    dsp *ha=HA->p[ifit+idm*HA->nx];
	    if(!ha) continue;
	    if(W1){
		sptmulmat(&actcplc->p[idm], ha, W1, 1);
	    }else{
		dmat *tmp=spsumabs(ha, 1);
		tmp->nx=tmp->ny; tmp->ny=1;
		dadd(&actcplc->p[idm], 1, tmp, 1);
		dfree(tmp);
	    }
	}
	/* If a point is fully coupled to any direction, it is treated as fully
	   coupled. */
	double thres=dmax(actcplc->p[idm])/HA->nx;
	double scale=1./thres;
	double *p=actcplc->p[idm]->p;
	for(long i=0; i<actcplc->p[idm]->nx; i++){
	    if(p[i]>thres){
		p[i]=1;
	    }else{
		p[i]*=scale;
	    }
	}
    }
    return actcplc;
}
/**
   Compute slaving actuator regularization. HA (or use GA) is used to compute
   active actuators. If NW is non NULL, orthogonalize it with the slaving
   regularization.  When the actuators are in the NULL space of HA, we want to
   contraint their values to be close to the ones that are active. We put an
   additional term in the fitting matrix to force this. Be careful with it when
   using tip/tilt constraint and cholesky back solve.  */
spcell *slaving(loc_t **aloc,  /**<[in]The actuator grid*/
		const dcell *actcplc,/**<[in]Actuator coupling coefficiency*/
		const dcell *NW,     /**<[in]The low rank terms that need to be orthogonal to the output (optional)*/
		const icell *actstuck,/**<[in]List of stuck actuators that will not be slaved, but have value constrained.*/
		const icell *actfloat,/**<[in]List of stuck actuators that will be slaved, but not have value constrained.*/
		const double thres,  /**<[in]The threshold that an actuator is deemed slave*/
		const double sclsq   /**<[in]The square of scaling of the overall strength*/
		){
    if(!actcplc && !actfloat) {
	error("Both actcplc and actfloat are not supplied\n");
    }
    double scl=sqrt(sclsq);
    if(scl<EPS){
	error("scl=%g is too small\n", scl);
    }
    int ndm=actcplc?actcplc->nx:actfloat->nx;
    spcell *actslavec=spcellnew(ndm, ndm);/*block diagonal. */
    PDSPCELL(actslavec, actslave);
    int nslavetot=0;
    /*Next process stuck and floating actuators. Adjust actcplc and compute slaving matrix.*/
    for(int idm=0; idm<ndm; idm++){
	int nact=aloc[idm]->nloc;
	int *stuck=calloc(nact, sizeof(int));
	int *floated=calloc(nact, sizeof(int));
	int nstuck=0;
	int nfloat=0;
	if(actstuck && actstuck->p[idm]){
	    nstuck=actstuck->p[idm]->nx;
	    for(int jact=0; jact<nstuck; jact++){
		int iact=actstuck->p[idm]->p[jact];
		stuck[iact]=1;
		actcplc->p[idm]->p[iact] = 1;/*Skip the stuck actuators. */
		warning2("Skip actuator %d in slaving\n", iact);
	    }
	}
	if(actfloat && actfloat->p[idm]){
	    nfloat=actfloat->p[idm]->nx;
	    for(int jact=0; jact<nfloat; jact++){
		int iact=actfloat->p[idm]->p[jact];
		floated[iact]=1;
		if(actcplc && actcplc->p[idm]->p[iact]>0.1){
		    error("actuator %d is floating, but actcpl is not zero\n", iact);
		}else{
		    actcplc->p[idm]->p[iact]=0;/*Slave this actuator. */
		}
	    }
	}
	
  	double *actcpl= actcplc->p[idm]->p;
	double *actcpl0 = actcpl-1;
	int nslave = 0;

	for(int iact=0; iact<nact; iact++){
	    if(actcpl[iact]<thres){
		nslave++;
	    }else if(stuck[iact]){
		nslave++;
	    }
	}
	nslavetot+=nslave;
	info2("dm %d: there are %d slave actuators\n", idm, nslave);
	if(nslave==0) {
	    continue;
	}
	loc_create_map(aloc[idm]);
	map_t *map=aloc[idm]->map;
	const double ox=map->ox;
	const double oy=map->oy;
	const double dx1=1./aloc[idm]->dx;
	const double dy1=1./aloc[idm]->dy;
	//dsp *slavet=spnew(nact,nslave,nslave*5);
	dsp *slavet=spnew(nact,nact,nslave*5);
	spint *pp=slavet->p;
	spint *pi=slavet->i;
	double *px=slavet->x;
	const double *locx=aloc[idm]->locx;
	const double *locy=aloc[idm]->locy;
	long count=0;
	for(int iact=0; iact<nact; iact++){
	    pp[iact]=count;
	    if(stuck && stuck[iact]){/*limit the strength of stuck actuators. */
		pi[count]=iact;
		px[count]=scl;
		count++;
	    }else if(actcpl[iact]<thres){/*slave actuators */
		long mapx=(long)round((locx[iact]-ox)*dx1);
		long mapy=(long)round((locy[iact]-oy)*dy1);
		int near_active=0;//neighbor is more coupled
		int near_exist=0;//neighbor exist
		double thres2=MAX(0.1, actcpl[iact]);//was 0.1
		for(int iy=-1; iy<2; iy++){
		    for(int ix=-1; ix<2; ix++){
			if(abs(ix+iy)!=1){
			    continue;/*skip self and corner */
			}
			int kact1=loc_map_get(map, mapx+ix, mapy+iy);
			if(kact1>0 && !stuck[kact1-1]){
			    near_exist++;
			    if(actcpl0[kact1]>thres2){//better than this one.
				near_active++;
			    }
			}
		
		    }
		}
		if(!near_exist){
		    error("This is an isolated actuator\n");
		}
		/*
		  neighbors are defined as the four pixels to the left, right,
		  top and bottom.  If some of the neighbors are active, use the
		  average of them for my value, otherwise, use the average of
		  all neighbors.
		*/
		double value=0;
		if(!near_active) value=-scl*0.1/near_exist;
		double valsum=0;
		/*part 1*/
		for(int iy=-1; iy<2; iy++){
		    for(int ix=-1; ix<2; ix++){
			if(abs(ix+iy)!=1){
			    continue;//skip self and corner
			}
			int kact1=loc_map_get(map, mapx+ix, mapy+iy);
			if(kact1>0 && !stuck[kact1-1]){
			    if(!near_active){
				valsum+=(px[count]=value);
				pi[count]=kact1-1;
				count++;
			    }else if(actcpl0[kact1]>actcpl[iact]){
				valsum+=(px[count]=-scl*MAX(0.1,actcpl0[kact1]));
				pi[count]=kact1-1;
				count++;
			    }
			}
		    }
		}
		
		/*part 2, matches negative sum of part 1*/
		pi[count]=iact;
		px[count]=-valsum;
		count++;
	    }/*if */
	}/*for iact */
	free(stuck);
	free(floated);
	pp[nact]=count;
	spsetnzmax(slavet, count);
	dsp *slave=sptrans(slavet);
	/*if(!disable_save){
	  spwrite(slave, "slave");
	  }*/
	actslave[idm][idm]=spmulsp(slavet, slave);
	if(NW && 0){
	    /*Now we need to make sure NW is in the NULL
	      space of the slaving regularization, especially
	      the tip/tilt constraints. NW=NW-slavet*inv(slavet)*NW */
	    if(NW->p[idm]){
		dmat *H=NULL;
		spfull(&H, slavet, 1);
		dmat *Hinv=dpinv(H,NULL,NULL);
		dmat *mod=NULL;
		dmm(&mod, 0, Hinv, NW->p[idm], "nn", 1);
		dmm(&NW->p[idm], 1, H, mod,"nn", -1);
		dfree(H);
		dfree(Hinv);
		dfree(mod);
		if(nfloat || nstuck){/*Remove corresponding rows in NW */
		    PDMAT(NW->p[idm], pNW);
		    for(int iy=0; iy<NW->p[idm]->ny; iy++){
			for(int jact=0; jact<nstuck; jact++){
			    int iact=actstuck->p[idm]->p[jact];
			    pNW[iy][iact]=0;
			}
			for(int jact=0; jact<nfloat; jact++){
			    int iact=actfloat->p[idm]->p[jact];
			    pNW[iy][iact]=0;
			}
		    }
		}
	    }
	}
	spfree(slave);
	spfree(slavet);
    }/*idm */
    if(nslavetot==0){
	spcellfree(actslavec);
	actslavec=NULL;
    }
    return actslavec;
}
/**
   When some actuators are stuck, remove the corresponding column in HA and/or HB
*/
void act_stuck(loc_t **aloc, const spcell *HA, const dcell *HB, const icell *stuck){
    if(!stuck || (!HA && !HB)) return; 
    int ndm=stuck->nx;
    for(int idm=0; idm<ndm; idm++){
	if(!stuck->p[idm]){
	    continue;
	}
	if(HA){
	    for(int ifit=0; ifit<HA->nx; ifit++){
		spint *pp=HA->p[idm*HA->nx+ifit]->p;
		double *px=HA->p[idm*HA->nx+ifit]->x;
		assert(HA->p[idm*HA->nx+ifit]->n==aloc[idm]->nloc);
		for(int jact=0; jact<stuck->p[idm]->nx; jact++){
		    int iact=stuck->p[idm]->p[jact];
		    for(int ic=pp[iact]; ic<pp[iact+1]; ic++){
			px[ic]=0;
		    }
		}
	    }
	}else if(HB){
	    for(int ifit=0; ifit<HB->nx; ifit++){
		dmat *hb=HB->p[idm*HB->nx+ifit];
		assert(hb->ny==aloc[idm]->nloc);
		for(int jact=0; jact<stuck->p[idm]->nx; jact++){
		    int iact=stuck->p[idm]->p[jact];
		    memset(hb->p+hb->nx*iact, 0, sizeof(double)*hb->nx);
		}
	    }
	}
    }
}
/**
   Zero out rows of dead actuators in mode vector.
*/
void act_zero(loc_t **aloc, const dcell *HB, const icell *dead){
    if(!dead || !HB) return;
    for(int idm=0; idm<dead->nx; idm++){
	if(!dead->p[idm]){
	    continue;
	}
	for(int imod=0; imod<HB->ny; imod++){
	    dmat *hb=HB->p[idm+imod*HB->nx];
	    assert(hb->nx==aloc[idm]->nloc);
	    PDMAT(hb, phb);
	    for(int jact=0; jact<dead->p[idm]->nx; jact++){
		int iact=dead->p[idm]->p[jact];
		for(int iy=0; iy<hb->ny; iy++){
		    phb[iy][iact]=0;
		}
	    }
	}
    }
}

/**
   When some actuators are float, remove the corresponding column in HA and/or HB,
   and add to neigh boring actuators. This is implemented using a second matrix and
   do matrix addition.*/
void act_float(loc_t **aloc, spcell **HA, const dcell *HB, const icell *actfloat){
    if(!actfloat || ((!HA || !*HA) && !HB)) return;
    int ndm=actfloat->nx;
    spcell *dHA=NULL;
    if(HA && *HA){
	int nfit=(*HA)->nx;
	dHA=spcellnew(nfit,ndm);
    }
    for(int idm=0; idm<ndm; idm++){
	if(!actfloat->p[idm]) continue;
	loc_create_map(aloc[idm]);
	map_t *map=aloc[idm]->map;
	double ox=map->ox;
	double oy=map->oy;
	double dx1=1./aloc[idm]->dx;
	double dy1=1./aloc[idm]->dy;
	const double *locx=aloc[idm]->locx;
	const double *locy=aloc[idm]->locy;
	long nact=aloc[idm]->nloc;
	long ndead=actfloat->p[idm]->nx;
	/*actuator that is floating */
	long *isfloat=calloc(nact, sizeof(long));
	/*which floating actuator to assign for this one. */
	long (*indfloat)[4]=calloc(nact, sizeof(long)*4);
	/*number of floating actuators that is assigned for this one. */
	long *nindfloat=calloc(nact, sizeof(long));
	/*active neighbors of each dead act. */
	long *neighbor=calloc(nact, sizeof(long));
	long nzmax=0;
	
	for(long jact=0; jact<ndead; jact++){
	    int iact=actfloat->p[idm]->p[jact];
	    isfloat[iact]=1;
	}
	/*loop through the floating actuators */
	for(long jact=0; jact<ndead; jact++){
	    int iact=actfloat->p[idm]->p[jact];
	    long mapx=(long)round((locx[iact]-ox)*dx1);
	    long mapy=(long)round((locy[iact]-oy)*dy1);
	    /*find all its neighbors. */
	    for(int iy=-1; iy<2; iy++){
		for(int ix=-1; ix<2; ix++){
		    if((ix!=0 && iy!=0) || (ix==0 && iy==0)){
			continue;/*skip center and corner */
		    }
		    int kact=loc_map_get(map, mapx+ix, mapy+iy)-1;
		    if(kact>-1){
			if(!isfloat[kact]){
			    indfloat[kact][nindfloat[kact]]=iact;
			    nindfloat[kact]++;
			    neighbor[iact]++;
			    nzmax+=4;/*assume 1 point couples to 4 points max. */
			}
		    }
		}
	    }
	}
	if(HA && *HA){
	    PDSPCELL(*HA, pHA);
	    PDSPCELL(dHA, pdHA);
	    /*Create dHA to assign weights of floating actuators to neighbors. */
	    for(int ifit=0; ifit<(*HA)->nx; ifit++){
		spint *pp=pHA[idm][ifit]->p;
		spint *pi=pHA[idm][ifit]->i;
		double *px=pHA[idm][ifit]->x;

		pdHA[idm][ifit]=spnew(pHA[idm][ifit]->m, pHA[idm][ifit]->n, nzmax);
		spint *pp2=pdHA[idm][ifit]->p;
		spint *pi2=pdHA[idm][ifit]->i;
		double *px2=pdHA[idm][ifit]->x;
		long count=0;
		for(long iact=0; iact<nact; iact++){
		    pp2[iact]=count;
		    if(nindfloat[iact]){
			for(long in=0; in<nindfloat[iact]; in++){
			    long jact=indfloat[iact][in];/*the floating act. */
			    double scale=1./neighbor[jact];
			    for(long ie=pp[jact]; ie<pp[jact+1]; ie++){
				if(count>=nzmax){
				    nzmax*=2;
				    spsetnzmax(pdHA[idm][ifit], nzmax);
				    pp2=pdHA[idm][ifit]->p;
				    pi2=pdHA[idm][ifit]->i;
				    px2=pdHA[idm][ifit]->x;
				}	
				pi2[count]=pi[ie];
				px2[count]=px[ie]*scale;
				count++;
			    }
			}
		    }
		}
		pp2[nact]=count;
		spsetnzmax(pdHA[idm][ifit], count);
		/*Remove weights of floating actuatorsf from HA. */
		for(int jact=0; jact<actfloat->p[idm]->nx; jact++){
		    int iact=actfloat->p[idm]->p[jact];
		    for(int ic=pp[iact]; ic<pp[iact+1]; ic++){
			px[ic]=0;
		    }
		}
	    }/*ifit */
	}else{/*Do dense matrix */
	    PDCELL(HB, hb);
	    for(long iact=0; iact<nact; iact++){
		if(nindfloat[iact]){
		    for(long in=0; in<nindfloat[iact]; in++){
			long jact=indfloat[iact][in];/*the floating act. */
			double scale=1./neighbor[jact];
			for(int ifit=0; ifit<HB->nx; ifit++){
			    dmat *hbi=hb[idm][ifit];
			    PDMAT(hbi, phbi);
			    for(long ix=0; ix<hbi->nx; ix++){
				phbi[iact][ix]+=phbi[jact][ix]*scale;
			    }
			}
		    }
		}
	    }
	    for(int jact=0; jact<actfloat->p[idm]->nx; jact++){
		int iact=actfloat->p[idm]->p[jact];
		for(int ifit=0; ifit<HB->nx; ifit++){
		    dmat *hbi=hb[idm][ifit];
		    PDMAT(hbi, phbi);
		    memset(phbi[iact], 0, hbi->nx*sizeof(double));
		}
	    }
	}
	free(nindfloat);
	free(indfloat);
	free(isfloat);
	free(neighbor);	
    }/*idm */
    if(HA){
	spcelladd(HA, dHA);
	spcellfree(dHA);
    }
}
/**
   Make DM actuator commands zero at stuck actuator locations.
*/
void act_stuck_cmd(loc_t **aloc, /**<[in] Actuator grid array*/
		   const dcell *adm,   /**<[in,out] Actuator command to process*/
		   const icell *stuck  /**<[in] List of stuck actuators*/
		   ){
    if(!adm || !stuck) return;
    PDCELL(adm, pa);
    int ndm=adm->nx;
    for(int idm=0; idm<ndm; idm++){
	if(!stuck->p[idm]) continue;
	for(int iy=0; iy<adm->ny; iy++){
	    dmat *pai=pa[iy][idm];
	    PDMAT(pai, px);
	    assert(pai->nx==aloc[idm]->nloc);
	    for(int icol=0; icol<pai->ny; icol++){
		for(int jact=0; jact<stuck->p[idm]->nx; jact++){
		    int iact=stuck->p[idm]->p[jact];
		    px[icol][iact]=0;
		}
	    }
	}
    }
}
/**
   Create an interpreter that make floating actuators equal to their neighbors.
*/
spcell* act_float_interp(loc_t **aloc,  /**<[in] Actuator grid array*/
			 const icell *actfloat/**<[in] List of floating actuators*/
			 ){
    if(!actfloat) return NULL;
    int ndm=actfloat->nx;
    spcell *out=spcellnew(ndm, ndm);
    for(int idm=0; idm<ndm; idm++){
	loc_create_map(aloc[idm]);
	map_t *map=aloc[idm]->map;
	double ox=map->ox;
	double oy=map->oy;
	double dx1=1./aloc[idm]->dx;	
	double dy1=1./aloc[idm]->dy;
	const double *locx=aloc[idm]->locx;
	const double *locy=aloc[idm]->locy;
	long nact=aloc[idm]->nloc;
	long ndead=actfloat->p[idm]?actfloat->p[idm]->nx:0;
	/*actuator that is floating */
	long *isfloat=calloc(nact, sizeof(long));
	for(long jact=0; jact<ndead; jact++){
	    int iact=actfloat->p[idm]->p[jact];
	    isfloat[iact]=1;
	}
	dsp *outit=spnew(nact, nact, nact*4);
	double *px=outit->x;
	spint *pp=outit->p;
	spint *pi=outit->i;
	long count=0;
	for(long iact=0; iact<nact; iact++){
	    pp[iact]=count;
	    if(isfloat[iact]){
		long mapx=(long)round((locx[iact]-ox)*dx1);
		long mapy=(long)round((locy[iact]-oy)*dy1);
		int count2=count;
		for(int idy=-1; idy<2; idy++){
		    for(int idx=-1; idx<2; idx++){
			if(abs(idx+idy)!=1){
			    continue;/*skip center and corner */
			}
			int kact=loc_map_get(map, mapx+idx, mapy+idy)-1;
			if(kact>-1){
			    pi[count]=kact;
			    px[count]=1.;
			    count++;
			}
		    }
		}
		if(count>count2){
		    double scl=1./(count-count2);
		    for(; count2<count; count2++){
			px[count2]=scl;
		    }
		}
	    }else{/*just copy over data */
		pi[count]=iact;
		px[count]=1;
		count++;
	    }
	}
	pp[nact]=count;
	out->p[idm+ndm*idm]=sptrans(outit);
	spfree(outit);
	free(isfloat);
    }
    return out;
}

/**
   Create an interpreter that make inactive actuators equal avearage of active
   neighbors if exist or all other neighbors.
*/
dsp* act_extrap_do(loc_t *aloc,        /**<[in] Actuator grid array*/
		  const dmat *actcplc,/**<[in] Actuator coupling coefficiency*/
		  const double thres  /**<[in] Threshold of coupling to turn on interpolation*/
		  ){
    dsp *out=0;
    const double *cpl=actcplc->p;
    const double *cpl0 = cpl-1;
    loc_create_map(aloc);
    map_t *map=aloc->map;
    double ox=map->ox;
    double oy=map->oy;
    double dx1=1./aloc->dx;
    double dy1=1./aloc->dy;
    const double *locx=aloc->locx;
    const double *locy=aloc->locy;
    long nact=aloc->nloc;
    dsp *outit=spnew(nact, nact, nact*4);
    double *px=outit->x;
    spint *pp=outit->p;
    spint *pi=outit->i;
    long count=0;
    for(long iact=0; iact<nact; iact++){
	pp[iact]=count;
	if(cpl[iact]<thres){
	    long mapx=(long)round((locx[iact]-ox)*dx1);
	    long mapy=(long)round((locy[iact]-oy)*dy1);
	    int count2=count;
	    double sum=0;
	    /*first, interpolate from neighbors of higher cpl*/
	    int near_active=0;
	    double thres2=0.1;
	    for(int iy=-1; iy<2; iy++){
		for(int ix=-1; ix<2; ix++){
		    if(abs(ix)+abs(iy)>1){
			continue;/*include self and neighbor*/
		    }
		    int kact1=loc_map_get(map, mapx+ix, mapy+iy);
		    if(kact1>0 && cpl0[kact1]>=thres2){
			near_active++;
		    }
		}
	    }
	    for(int iy=-1; iy<2; iy++){
		for(int ix=-1; ix<2; ix++){
		    if(abs(ix)+abs(iy)>1){
			continue;/*skip corner only*/
		    }
		    int kact1=loc_map_get(map, mapx+ix, mapy+iy);
		    if(kact1>0){
			if(!near_active){
			    pi[count]=kact1-1;
			    sum+=(px[count]=1);
			    count++;
			}else if(cpl0[kact1] >=thres2){
			    pi[count]=kact1-1;
			    sum+=(px[count]=cpl0[kact1]);
			    count++;
			}
		    }
		}
	    }
	    if(count>count2){
		double scl=1./sum;
		for(;count2<count; count2++){
		    px[count2]*=scl;
		}
	    }
	}else{
	    /*just copy over data */
	    pi[count]=iact;
	    px[count]=1;
	    count++;
	}
    }
    pp[nact]=count;
    out=sptrans(outit);
    spfree(outit);
	
    /*The above interpolation only propagate the value one step. Multiple the
      interpolator a few times to propagate longer.*/
    for(int i=0; i<5; i++){
	dsp *tmp=spmulsp(out,out);
	spfree(out);
	out=tmp;
    }
    return out;
}
spcell* act_extrap(loc_t **aloc,     /**<[in] Actuator grid array*/
		   const dcell *actcplc,/**<[in] Actuator coupling coefficiency*/
		   const double thres /**<[in] Threshold of coupling to turn on interpolation*/
		   ){
    int ndm=actcplc->nx;
    spcell *out=spcellnew(ndm, ndm);
    for(int idm=0; idm<ndm; idm++){
	out->p[idm+ndm*idm]=act_extrap_do(aloc[idm], actcplc->p[idm], thres);
    }
    return out;
}

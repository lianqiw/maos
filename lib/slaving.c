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
#include "../math/mathdef.h"
#include "slaving.h"

/**
   Compute the actuator coupling coefficient to be used to identify non-coupled
   actuators. W1 is optional weighting function. The max value is 1.
*/
dcell* genactcpl(const dspcell* HA, const dmat* W1){
	int ndm=HA->ny;
	dcell* actcplc=dcellnew(ndm, 1);
	for(int idm=0; idm<ndm; idm++){
		for(int ifit=0; ifit<HA->nx; ifit++){
			dsp* ha=HA->p[ifit+idm*HA->nx];
			if(!ha) continue;
			if(W1){
				dspmm(&actcplc->p[idm], ha, W1, "tn", 1);
			} else{
				dmat* tmp=dspsumabs(ha, 1);
				tmp->nx=tmp->ny; tmp->ny=1;
				dadd(&actcplc->p[idm], 1, tmp, 1);
				dfree(tmp);
			}
		}
		/* If a point is fully coupled to any direction, it is treated as fully
		   coupled. */
		real thres=dmax(actcplc->p[idm])/HA->nx;
		real scale=1./thres;
		real* p=actcplc->p[idm]->p;
		for(long i=0; i<actcplc->p[idm]->nx; i++){
			if(p[i]>thres||thres==0){
				p[i]=1;
			} else{
				p[i]*=scale;
			}
		}
	}
	return actcplc;
}
/**
   Compute slaving actuator regularization term.


   Inactive actuators (coupling coefficiency below the threshold) are slaved to neighbors that are active.

   The result is S=H'*H; Let b=H*a, we have
   if mode==1: b(isa)=alpha*(a(isa)-sum(a(jsa).*weight(jsa))), where isa is inactive and jsa are neighoring active actuators.
   if mode==2: b(igroup)=alpha*sum(a(jsa)), where jsa are fully active actuators belonging to igroup
   if mode==3: b(isa)=alpha*sum(a(jsa)-a(ksa)), where isa is inactive, and jsa and ksa are neighrboring and opposite active actuator

*/
dspcell* slaving(loccell* aloc,        /**<[in]The actuator grid*/
	const dcell* actcplc, /**<[in]Actuator coupling coefficiency*/
	const lcell* actstuck,/**<[in]mask for stuck actuators that will not be slaved, but have value constrained.*/
	const lcell* actfloat,/**<[in]mask for float actuators that will be slaved, but not have value constrained.*/
	const real thres,   /**<[in]The threshold that an actuator is deemed slave*/
	const real sclsq,   /**<[in] Expected norm of the slaving matrix.*/
	const int mode        /**<[in] Mode of operation. (1 or 2)*/
){
	if(!actcplc&&!actfloat){
		error("Both actcplc and actfloat are not supplied\n");
	}
	int ndm=aloc->nx;
	dspcell* actslavec=(dspcell*)cellnew(ndm, ndm);/*block diagonal. */
	dspcell* actslave=actslavec;
	int nslavetot=0;
	/*Next process stuck and floating actuators. Adjust actcplc and compute slaving matrix.*/
	for(int idm=0; idm<ndm; idm++){
		int nact=aloc->p[idm]->nloc;
		int nslave=0;
		real* actcpl=actcplc->p[idm]->p;
		real* actcpl0=actcpl-1;
		const long* isstuck=(actstuck&&actstuck->p[idm])?actstuck->p[idm]->p:0;
		const long* isfloat=(actfloat&&actfloat->p[idm])?actfloat->p[idm]->p:0;
		for(int iact=0; iact<nact; iact++){
			if(isstuck&&isstuck[iact]){
				actcpl[iact]=1;/*always Skip the stuck actuators. */
				nslave++;
			}
			if(isfloat&&isfloat[iact]){
				actcpl[iact]=0;/*always include the float actuators */
				nslave++;
			}
			if(actcpl[iact]<thres){
				nslave++;
			}
		}

		nslavetot+=nslave;
		info("dm %d, mode %d: there are %d slave actuators\n", idm, mode, nslave);

		loc_create_map(aloc->p[idm]);
		map_t* map=aloc->p[idm]->map;
		const real ox=map->ox;
		const real oy=map->oy;
		const real dx1=1./aloc->p[idm]->dx;
		const real dy1=1./aloc->p[idm]->dy;

		const real* locx=aloc->p[idm]->locx;
		const real* locy=aloc->p[idm]->locy;
		dsp* slavet=0;
		long count=0;
		if(mode==2){//falls back to tikhonov when ngroup==0
			lmat* group=lnew(nact, 1);
			lmat* groupc=lnew(nact, 1);
			long ngroup=1;//number of isolated islands
			long ngsub=0;//over counted number of ngroup.
			for(int iact=0; iact<nact; iact++){
				if(isstuck&&isstuck[iact]) continue;
				if(actcpl[iact]>thres){/*active actuators */
					long mapx=(long)round((locx[iact]-ox)*dx1);
					long mapy=(long)round((locy[iact]-oy)*dy1);
					long found=0;
					//Check if any of its neighbor has group assigned
					for(int iy=-1; iy<2; iy++){
						for(int ix=-1; ix<2; ix++){
							int kact1=loc_map_get(map, mapx+ix, mapy+iy);
							if(actcpl0[kact1]>thres){
								if(P(group, kact1-1)){
									if(found){
										if(found!=P(group, kact1-1)){
											int toreplace=P(group, kact1-1);
											for(int jact=0; jact<nact; jact++){
												if(P(group, jact)==toreplace){
													P(group, jact)=found;
												}
											}
											ngsub++;
											P(groupc, toreplace)=0;
										}
									} else{
										found=P(group, kact1-1);
									}
								}
							}
						}
					}
					//Assign a value if no neighbor is assigned
					if(!found){
						found=ngroup;
						P(groupc, found)=1;
						ngroup++;
					}
					//Assign all neighboring ones and itself
					for(int iy=-1; iy<2; iy++){
						for(int ix=-1; ix<2; ix++){
							int kact1=loc_map_get(map, mapx+ix, mapy+iy);
							if(actcpl0[kact1]>thres){
								P(group, kact1-1)=found;
							}
						}
					}
				}
			}
			lmat* groupu=lnew(ngroup-ngsub-1, 1);//list of uniq group value
			int ic=0;
			for(int ig=1; ig<ngroup; ig++){
				if(P(groupc, ig)){
					P(groupu, ic)=ig;
					ic++;
				}
			}
			lfree(groupc);
			ngroup=ic;
			/*if(ngroup<2){
			warning("Only one discontinuous region found\n");
			}*/

			slavet=dspnew(nact, ngroup, nact);
			spint* pp=slavet->p;
			spint* pi=slavet->i;
			real* px=slavet->x;
			for(int igroup=0; igroup<ngroup; igroup++){
				pp[igroup]=count;
				int count2=count;
				for(int iact=0; iact<nact; iact++){
					if(P(group, iact)==P(groupu, igroup)){
						pi[count]=iact;
						count++;
					}
				}
				real scale=1./(count-count2);
				for(; count2<count; count2++){
					px[count2]=scale;
				}
			}
			pp[slavet->ny]=count;

			lfree(groupu);
			lfree(group);
		} else if(nslave>0){
			slavet=dspnew(nact, nact, nslave*9);
			lmat* actct=lnew(nact, 1);
			spint* pp=slavet->p;
			spint* pi=slavet->i;
			real* px=slavet->x;

			for(int iact=0; iact<nact; iact++){
				pp[iact]=count;
				if(isstuck&&isstuck[iact]){/*limit the strength of stuck actuators. */
					pi[count]=iact;
					px[count]=1;
					count++;
				} else if(actcpl[iact]<thres){/*slave actuators */
					long mapx=(long)round((locx[iact]-ox)*dx1);
					long mapy=(long)round((locy[iact]-oy)*dy1);
					int near_activer=0;//neighbor is more coupled
					int near_exist=0;//neighbor exist
					int near_inactive=0;//inactive neighbor
					real thres2=MAX(0.1, actcpl[iact]);//was 0.1
					for(int iy=-1; iy<2; iy++){
						for(int ix=-1; ix<2; ix++){
							if(ix==0&&iy==0) continue;//skip self
							int kact1=loc_map_get(map, mapx+ix, mapy+iy);
							if(kact1>0&&(!isstuck||!isstuck[kact1-1])){
								if(abs(ix+iy)==1){
									near_exist++;
									if(actcpl0[kact1]>thres2){//better than this one.
										near_activer++;
									}

									if(actcpl0[kact1]<thres){
										near_inactive++;
									}
								}
							}

						}
					}
					if(!near_exist){
						warning("This is an isolated actuator\n");
					} else if(mode==1){//Determining inactive actuators using neighbors.
					/*
					  neighbors are defined as the four pixels to the left, right,
					  top and bottom.  If some of the neighbors are active, use the
					  average of them for my value, otherwise, use the average of
					  all neighbors.
					*/
					/*part 1*/
						real valsum=0;
						for(int iy=-1; iy<2; iy++){
							for(int ix=-1; ix<2; ix++){
								if(abs(ix+iy)!=1){
									continue;//skip self and corner
								}
								int kact1=loc_map_get(map, mapx+ix, mapy+iy);
								if(kact1>0&&(!isstuck||!isstuck[kact1-1])){
									if(!near_activer){//simply average neighoring ones
										valsum+=(px[count]=-0.1/near_exist);
										pi[count]=kact1-1;
										count++;
									} else if(actcpl0[kact1]>actcpl[iact]){//weighted average neighrboring ones.
										valsum+=(px[count]=-MAX(0.1, actcpl0[kact1]));
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
					} else if(mode==3&&near_inactive>1&&near_inactive<7){

					/*Limit difference of neighbors actuators of inactive
					 * actuators. To minimize island effect*/

					//This one works in dual DM case.
						for(int iy=-1; iy<2; iy++){
							for(int ix=-1; ix<=MIN(iy, 0); ix++){
							//if(abs(ix+iy)!=1) continue;//skip self and corner
								if(ix==0&&iy==0) continue; //skip self
								int kact1=loc_map_get(map, mapx+ix, mapy+iy);
								int kact2=loc_map_get(map, mapx-ix, mapy-iy);
								if(kact1>0&&(!isstuck||!isstuck[kact1-1])
									&&kact2>0&&(!isstuck||!isstuck[kact2-1])
									&&actcpl0[kact1]>thres&&actcpl0[kact2]>thres
									&&!P(actct, kact1-1)&&!P(actct, kact2-1)
									){
									px[count]=1;
									pi[count]=kact1-1;
									P(actct, kact1-1)++;
									count++;

									px[count]=-1;
									pi[count]=kact2-1;
									P(actct, kact2-1)++;
									count++;
								}
							}
						}
					} else if(mode==4){
					//this one does no work
						int sign=aloc->p[idm]->locy[iact]>0?1:-1;
						int count2=0;
						int offx=0, offy=0;
						for(int id=0; id<4; id++){
							switch(id){
							case 0:
								offx=0; offy=0; break;
							case 1:
								offx=1; offy=0; break;
							case 2:
								offx=0; offy=1; break;
							case 3:
								offx=1; offy=1; break;
							default:
								continue;
							}

							int last_ic=-1;
							for(int ic=0; ic<4; ic++){
								int ix, iy;
								switch(ic){
								case 0:
									ix=-1; iy=-1; break;
								case 1:
									ix=-1; iy=0; break;
								case 2:
									ix=-1; iy=1; break;
								case 3:
									ix=0; iy=1; break;
								default:
									continue;
								}
								int kact1=loc_map_get(map, mapx+ix*(1+offx), mapy+iy*(1+offy));
								int kact2=loc_map_get(map, mapx-ix, mapy-iy);
								//if actuators at opposite edge/corners are active
								if(kact1>0&&(!isstuck||!isstuck[kact1-1])&&
									kact2>0&&(!isstuck||!isstuck[kact2-1])&&
									actcpl0[kact1]>thres&&actcpl0[kact2]>thres
									&&!P(actct, kact1-1)&&!P(actct, kact2-1)
									){
									if(ic>last_ic+1) sign=-sign;
									last_ic=ic;
									px[count]=sign;
									pi[count]=kact1-1;
									P(actct, kact1-1)++;
									count++;
									px[count]=-sign;
									pi[count]=kact2-1;
									P(actct, kact2-1)++;
									count++;
									count2++;
								}
							}
							if(count2>0) break;
						}
						/*if(count2==1){//false positive
							count-=2;
							}*/
					}
				}/*if */
			}/*for iact */
			pp[nact]=count;
			lfree(actct);
		}
		dspsetnzmax(slavet, count);
		//writebin(slavet, "slave_%d_%d", mode, idm);
		P(actslave, idm, idm)=dspmulsp(slavet, slavet, "nt");
		dspscale(P(actslave, idm, idm), sclsq);
		dspfree(slavet);
	}/*idm */

	return actslavec;
}
/**
   When some actuators are stuck, zero the corresponding column in HA
*/
void act_stuck(loccell* aloc, void* HA_, const lcell* stuck){
	if(!stuck||!HA_) return;
	cell* HA=(cell*)HA_;
	int ndm=aloc->nx;
	int nfit=0;
	if(HA->ny==ndm){
		nfit=HA->nx;
	} else if(HA->ny==1&&HA->nx==ndm){
		nfit=1;
	} else{
		error("HA: Invalid format %ldx%ld\n", HA->nx, HA->ny);
	}
	for(int idm=0; idm<ndm; idm++){
		if(!stuck->p[idm]){
			continue;
		}
		const int nact=aloc->p[idm]->nloc;
		for(int ifit=0; ifit<nfit; ifit++){
			cell* HAi=HA->p[idm*nfit+ifit];
			if(HAi->id==M_REAL){//dense
				dmat* hb=(dmat*)HAi;
				if(hb->nx>1&&hb->ny==aloc->p[idm]->nloc){
					//modifying interaction matrix
					for(int iact=0; iact<nact; iact++){
						if(stuck->p[idm]->p[iact]){
							memset(hb->p+hb->nx*iact, 0, sizeof(real)*hb->nx);
						}
					}
				} else if(hb->ny==1&&hb->nx==aloc->p[idm]->nloc){
					//modifying coupling vector
					for(int iact=0; iact<nact; iact++){
						if(stuck->p[idm]->p[iact]){
							hb->p[iact]=0;
						}
					}
				} else{
					error("Invalid input: hb is %ldx%ld, nloc=%ld\n", hb->nx, hb->ny, aloc->p[idm]->nloc);
				}
			} else if(HAi->id==M_DSP){//sparse
				dsp* ha=(dsp*)HAi;
				spint* pp=ha->p;
				real* px=ha->x;
				assert(ha->ny==aloc->p[idm]->nloc);
				for(int iact=0; iact<nact; iact++){
					if(stuck->p[idm]->p[iact]){
						for(int ic=pp[iact]; ic<pp[iact+1]; ic++){
							px[ic]=0;
						}
					}
				}
			} else{
				error("Invalid parameter: HA.id=%u\n", HAi->id);
			}
		}
	}
}
/**
   Zero out rows of dead actuators in mode vector.
*/
void act_zero(loccell* aloc, const dcell* HB, const lcell* dead){
	if(!dead||!HB) return;
	for(int idm=0; idm<dead->nx; idm++){
		if(!dead->p[idm]){
			continue;
		}
		const int nact=aloc->p[idm]->nloc;
		if(HB->nx!=aloc->nx){
			error("HB: Invalid format\n");
		}
		for(int imod=0; imod<HB->ny; imod++){
			dmat* hb=HB->p[idm+imod*HB->nx];
			if(hb->nx!=aloc->p[idm]->nloc){
				error("hb: Invalid format\n");
			}
			dmat* phb=hb;
			for(int iact=0; iact<nact; iact++){
				if(dead->p[idm]->p[iact]){
					for(int iy=0; iy<hb->ny; iy++){
						P(phb, iact, iy)=0;
					}
				}
			}
		}
	}
}

/**
   When some actuators are float, remove the corresponding column in HA and/or HB,
   and add to neigh boring actuators. This is implemented using a second matrix and
   then add to the original matrix.*/
void act_float(loccell* aloc, dspcell** HA, const dcell* HB, const lcell* actfloat){
	if(!actfloat||((!HA||!*HA)&&!HB)) return;
	int ndm=actfloat->nx;
	dspcell* dHA=NULL;
	if(HA&&*HA){
		int nfit=(*HA)->nx;
		dHA=dspcellnew(nfit, ndm);
	}
	for(int idm=0; idm<ndm; idm++){
		if(!actfloat->p[idm]) continue;
		loc_create_map(aloc->p[idm]);
		map_t* map=aloc->p[idm]->map;
		real ox=map->ox;
		real oy=map->oy;
		real dx1=1./aloc->p[idm]->dx;
		real dy1=1./aloc->p[idm]->dy;
		const real* locx=aloc->p[idm]->locx;
		const real* locy=aloc->p[idm]->locy;
		long nact=aloc->p[idm]->nloc;
		/*which floating actuator to assign for this one. */
		lmat* indfloat=lnew(4, nact);
		/*number of floating actuators that is assigned for this one. */
		long* nindfloat=mycalloc(nact, long);
		/*active neighbors of each dead act. */
		long* neighbor=mycalloc(nact, long);
		long nzmax=0;
		long* isfloat=actfloat->p[idm]->p;
		/*loop through the floating actuators */
		for(int iact=0; iact<nact; iact++){
			if(!actfloat->p[idm]->p[iact]) continue;
			long mapx=(long)round((locx[iact]-ox)*dx1);
			long mapy=(long)round((locy[iact]-oy)*dy1);
			/*find all its neighbors. */
			for(int iy=-1; iy<2; iy++){
				for(int ix=-1; ix<2; ix++){
					if((ix!=0&&iy!=0)||(ix==0&&iy==0)){
						continue;/*skip center and corner */
					}
					int kact=loc_map_get(map, mapx+ix, mapy+iy)-1;
					if(kact>-1){
						if(!isfloat[kact]){
							P(indfloat, nindfloat[kact], kact)=iact;
							nindfloat[kact]++;
							neighbor[iact]++;
							nzmax+=4;/*assume 1 point couples to 4 points max. */
						}
					}
				}
			}
		}
		if(HA&&*HA){
			dspcell* pHA=*HA;
			dspcell* pdHA=dHA;
			/*Create dHA to assign weights of floating actuators to neighbors. */
			for(int ifit=0; ifit<(*HA)->nx; ifit++){
				spint* pp=P(pHA, ifit, idm)->p;
				spint* pi=P(pHA, ifit, idm)->i;
				real* px=P(pHA, ifit, idm)->x;

				P(pdHA, ifit, idm)=dspnew(P(pHA, ifit, idm)->nx, P(pHA, ifit, idm)->ny, nzmax);
				spint* pp2=P(pdHA, ifit, idm)->p;
				spint* pi2=P(pdHA, ifit, idm)->i;
				real* px2=P(pdHA, ifit, idm)->x;
				long count=0;
				for(long iact=0; iact<nact; iact++){
					pp2[iact]=count;
					if(nindfloat[iact]){
						for(long in=0; in<nindfloat[iact]; in++){
							long jact=P(indfloat, in, iact);/*the floating act. */
							real scale=1./neighbor[jact];
							for(long ie=pp[jact]; ie<pp[jact+1]; ie++){
								if(count>=nzmax){
									nzmax*=2;
									dspsetnzmax(P(pdHA, ifit, idm), nzmax);
									pp2=P(pdHA, ifit, idm)->p;
									pi2=P(pdHA, ifit, idm)->i;
									px2=P(pdHA, ifit, idm)->x;
								}
								pi2[count]=pi[ie];
								px2[count]=px[ie]*scale;
								count++;
							}
						}
					}
				}
				pp2[nact]=count;
				dspsetnzmax(P(pdHA, ifit, idm), count);
				/*Remove weights of floating actuatorsf from HA. */
				for(int jact=0; jact<actfloat->p[idm]->nx; jact++){
					int iact=actfloat->p[idm]->p[jact];
					for(int ic=pp[iact]; ic<pp[iact+1]; ic++){
						px[ic]=0;
					}
				}
			}/*ifit */
		} else{/*Do dense matrix */
			for(long iact=0; iact<nact; iact++){
				if(nindfloat[iact]){
					for(long in=0; in<nindfloat[iact]; in++){
						long jact=P(indfloat, in, iact);/*the floating act. */
						real scale=1./neighbor[jact];
						for(int ifit=0; ifit<HB->nx; ifit++){
							dmat* hbi=P(HB, ifit, idm);
							dmat* phbi=hbi;
							for(long ix=0; ix<hbi->nx; ix++){
								P(phbi, ix, iact)+=P(phbi, ix, jact)*scale;
							}
						}
					}
				}
			}
			for(int jact=0; jact<actfloat->p[idm]->nx; jact++){
				int iact=actfloat->p[idm]->p[jact];
				for(int ifit=0; ifit<HB->nx; ifit++){
					dmat* hbi=P(HB, ifit, idm);
					memset(PCOL(hbi, iact), 0, hbi->nx*sizeof(real));
				}
			}
		}
		free(nindfloat);
		lfree(indfloat);
		free(neighbor);
	}/*idm */
	if(HA){
		dcelladd(HA, 1, dHA, 1);
		dspcellfree(dHA);
	}
}
/**
   Make DM actuator commands zero at stuck actuator locations.
*/
void act_stuck_cmd(loccell* aloc, /**<[in] Actuator grid array*/
	const dcell* adm,   /**<[in,out] Actuator command to process*/
	const lcell* stuck  /**<[in] List of stuck actuators*/
){
	if(!adm||!stuck) return;
	const int ndm=aloc->nx;
	for(int idm=0; idm<ndm; idm++){
		if(!stuck->p[idm]) continue;
		const int nact=aloc->p[idm]->nloc;
		for(int iy=0; iy<adm->ny; iy++){
			dmat* pai=P(adm, idm, iy);
			dmat* px=pai;
			assert(pai->nx==aloc->p[idm]->nloc);
			for(int icol=0; icol<pai->ny; icol++){
				for(int iact=0; iact<nact; iact++){
					if(stuck->p[idm]->p[iact]){
						P(px, iact, icol)=0;
					}
				}
			}
		}
	}
}
/**
   Create an interpreter that make floating actuators equal to their neighbors.
*/
dspcell* act_float_interp(loccell* aloc,  /**<[in] Actuator grid array*/
	const lcell* actfloat/**<[in] List of floating actuators*/
){
	int ndm=aloc->nx;
	dspcell* out=dspcellnew(ndm, ndm);
	for(int idm=0; idm<ndm; idm++){
		loc_create_map(aloc->p[idm]);
		map_t* map=aloc->p[idm]->map;
		real ox=map->ox;
		real oy=map->oy;
		real dx1=1./aloc->p[idm]->dx;
		real dy1=1./aloc->p[idm]->dy;
		const real* locx=aloc->p[idm]->locx;
		const real* locy=aloc->p[idm]->locy;
		long nact=aloc->p[idm]->nloc;
		/*actuator that is floating */
		long* isfloat=(actfloat&&actfloat->p[idm])?actfloat->p[idm]->p:0;
		dsp* outit=dspnew(nact, nact, nact*4);
		real* px=outit->x;
		spint* pp=outit->p;
		spint* pi=outit->i;
		long count=0;
		for(long iact=0; iact<nact; iact++){
			pp[iact]=count;
			if(isfloat&&isfloat[iact]){
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
					real scl=1./(count-count2);
					for(; count2<count; count2++){
						px[count2]=scl;
					}
				}
			} else{/*just copy over data */
				pi[count]=iact;
				px[count]=1;
				count++;
			}
		}
		pp[nact]=count;
		out->p[idm+ndm*idm]=dsptrans(outit);
		dspfree(outit);
	}
	return out;
}

/**
   Create an interpreter that make inactive actuators equal avearage of active
   neighbors if exist or all other neighbors.
*/
static dsp* act_extrap_do(loc_t* aloc,        /**<[in] Actuator grid array*/
	const dmat* actcplc,/**<[in] Actuator coupling coefficiency*/
	const real thres  /**<[in] Threshold of coupling to turn on interpolation*/
){
	TIC;tic;

	dsp* out=0;
	const real* cpl=actcplc->p;
	const real* cpl0=cpl-1;
	loc_create_map(aloc);
	map_t* map=aloc->map;
	real ox=map->ox;
	real oy=map->oy;
	real dx1=1./aloc->dx;
	real dy1=1./aloc->dy;
	const real* locx=aloc->locx;
	const real* locy=aloc->locy;
	long nact=aloc->nloc;
	dsp* outit=dspnew(nact, nact, nact*4);
	real* px=outit->x;
	spint* pp=outit->p;
	spint* pi=outit->i;
	long count=0;
	for(long iact=0; iact<nact; iact++){
		pp[iact]=count;
		if(cpl[iact]<thres){
			long mapx=(long)round((locx[iact]-ox)*dx1);
			long mapy=(long)round((locy[iact]-oy)*dy1);
			int count2=count;
			real sum=0;
			/*first, interpolate from neighbors of higher cpl*/
			int near_active=0;
			//real thres2=0.1;
			real thres2=MAX(0.1, cpl[iact]);
			//Find active neighbors
			for(int iy=-1; iy<2; iy++){
				for(int ix=-1; ix<2; ix++){
					if(abs(ix)+abs(iy)==2){
						continue;//skip corner
					}
					//include self and neighbor
					int kact1=loc_map_get(map, mapx+ix, mapy+iy);
					if(kact1>0&&cpl0[kact1]>=thres2){
						near_active++;
					}
				}
			}
			for(int iy=-1; iy<2; iy++){
				for(int ix=-1; ix<2; ix++){
					if(abs(ix)+abs(iy)==2){
						continue;
					}
					int kact1=loc_map_get(map, mapx+ix, mapy+iy);
					if(kact1>0){
						if(!near_active){
							//there are no active neighbors, use the average
							pi[count]=kact1-1;
							sum+=(px[count]=1);
							count++;
						} else if(cpl0[kact1]>=thres2||(ix==0&&iy==0)){
							//there are a few active neighbors
							pi[count]=kact1-1;
							sum+=(px[count]=cpl0[kact1]);
							count++;
						}
					}
				}
			}
			if(count>count2){
				real scl=1./sum;
				for(;count2<count; count2++){
					px[count2]*=scl;
				}
			}
		} else{
			/*just copy over data */
			pi[count]=iact;
			px[count]=1;
			count++;
		}
	}
	pp[nact]=count;
	out=dsptrans(outit);
	dspfree(outit);

	//New method. Test convergence
	dmat* x=dnew(out->nx, 1); dset(x, 1);
	dmat* y1=dnew(out->nx, 1);
	dmat* y2=dnew(out->nx, 1);
	dcellmm(&y1, out, x, "tn", 1);
	real diff;
	count=0;
	do{
		dsp* tmp=dspmulsp(out, out, "nn");
		dspfree(out);
		out=tmp;
		dzero(y2);
		dcellmm(&y2, out, x, "tn", 1);
		diff=ddiff(y1, y2);
		dcp(&y1, y2);
		count++;
	} while(diff>EPS);
	dfree(x);
	dfree(y1);
	dfree(y2);
	dspdroptol(out, 1e-3);
	toc("act_extrap: %ld iterations", count);

	return out;
}
/**
   Create an interpreter that make inactive actuators equal avearage of active
   neighbors if exist or all other neighbors.
*/
dspcell* act_extrap(loccell* aloc,     /**<[in] Actuator grid array*/
	const dcell* actcplc,/**<[in] Actuator coupling coefficiency*/
	const real thres /**<[in] Threshold of coupling to turn on interpolation*/
){
	int ndm=actcplc->nx;
	dspcell* out=dspcellnew(ndm, ndm);
	for(int idm=0; idm<ndm; idm++){
		out->p[idm+ndm*idm]=act_extrap_do(aloc->p[idm], actcplc->p[idm], thres);
	}
	return out;
}

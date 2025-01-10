/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "zernike.h"
/**
   Compute the actuator coupling coefficient to be used to identify non-coupled
   actuators. W1 is optional weighting function. The max value is 1.
*/
dcell* genactcpl(const_anyarray HA_, const dmat* W1){
	const cell *HA=HA_.c;
	int ndm=NY(HA);
	dcell* actcplc=dcellnew(ndm, 1);
	for(int idm=0; idm<ndm; idm++){
		for(int ifit=0; ifit<NX(HA); ifit++){
			cell* ha=P(HA, ifit, idm);
			if(!ha) continue;
			if(W1){
				dcellmm(&P(actcplc,idm), ha, W1, "tn", 1);
			} else{
				dmat* tmp=dspsumabs(dsp_cast(ha), 1);
				reshape(tmp, NY(tmp), 1);
				dadd(&P(actcplc,idm), 1, tmp, 1);
				dfree(tmp);
			}
		}
		/* If a point is fully coupled to any direction, it is treated as fully
		   coupled. */
		real thres=dmax(P(actcplc,idm))/NX(HA);
		real scale=1./thres;
		real* p=P(P(actcplc,idm));
		for(long i=0; i<P(actcplc,idm)->nx; i++){
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
	int ndm=NX(aloc);
	dspcell* actslavec=(dspcell*)cellnew(ndm, ndm);/*block diagonal. */
	dspcell* actslave=actslavec;
	//int nslavetot=0;
	/*Next process stuck and floating actuators. Adjust actcplc and compute slaving matrix.*/
	for(int idm=0; idm<ndm; idm++){
		int nact=P(aloc,idm)->nloc;
		int nslave=0;
		real* actcpl=P(P(actcplc,idm));
		real* actcpl0=actcpl-1;
		const long* isstuck=(actstuck&&P(actstuck,idm))?P(P(actstuck,idm)):0;
		const long* isfloat=(actfloat&&P(actfloat,idm))?P(P(actfloat,idm)):0;
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

		//nslavetot+=nslave;
		dbg("dm %d, mode %d: there are %d slave actuators\n", idm, mode, nslave);

		loc_create_map(P(aloc,idm));
		map_t* map=P(aloc,idm)->map;
		const real ox=map->ox;
		const real oy=map->oy;
		const real dx1=1./P(aloc,idm)->dx;
		const real dy1=1./P(aloc,idm)->dy;

		const real* locx=P(aloc,idm)->locx;
		const real* locy=P(aloc,idm)->locy;
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
			spint* pp=slavet->pp;
			spint* pi=slavet->pi;
			real* px=slavet->px;
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
			pp[NY(slavet)]=count;

			lfree(groupu);
			lfree(group);
		} else if(nslave>0){
			slavet=dspnew(nact, nact, nslave*9);
			lmat* actct=lnew(nact, 1);
			spint* pp=slavet->pp;
			spint* pi=slavet->pi;
			real* px=slavet->px;

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
						int sign=P(aloc,idm)->locy[iact]>0?1:-1;
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
void act_stuck(loccell* aloc, anyarray HA_, const lcell* stuck){
	cell *HA=HA_.c;
	if(!stuck||!HA) return;
	int ndm=NX(aloc);
	int nfit=0;
	if(NY(HA)==ndm){
		nfit=NX(HA);
	} else if(NY(HA)==1&&NX(HA)==ndm){
		nfit=1;
	} else{
		error("HA: Invalid format %ldx%ld\n", NX(HA), NY(HA));
	}
	for(int idm=0; idm<ndm; idm++){
		if(!P(stuck,idm)){
			continue;
		}
		const int nact=P(aloc,idm)->nloc;
		for(int ifit=0; ifit<nfit; ifit++){
			cell* HAi=nfit==1?P(HA,idm):P(HA,ifit,idm);
			if(HAi->id==M_REAL){//dense
				dmat* hb=(dmat*)HAi;
				if(NX(hb)>1&&NY(hb)==P(aloc,idm)->nloc){
					//modifying interaction matrix
					for(int iact=0; iact<nact; iact++){
						if(P(P(stuck,idm),iact)){
							memset(PCOL(hb, iact), 0, sizeof(real)*NX(hb));
						}
					}
				} else if(NY(hb)==1&&NX(hb)==P(aloc,idm)->nloc){
					//modifying coupling vector
					for(int iact=0; iact<nact; iact++){
						if(P(P(stuck,idm),iact)){
							P(hb,iact)=0;
						}
					}
				} else{
					error("Invalid input: hb is %ldx%ld, nloc=%ld\n", NX(hb), NY(hb), P(aloc,idm)->nloc);
				}
			} else if(HAi->id==M_DSP){//sparse
				dsp* ha=(dsp*)HAi;
				spint* pp=ha->pp;
				real* px=ha->px;
				assert(NY(ha)==P(aloc,idm)->nloc);
				for(int iact=0; iact<nact; iact++){
					if(P(P(stuck,idm),iact)){
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
	for(int idm=0; idm<NX(dead); idm++){
		if(!P(dead,idm)){
			continue;
		}
		const int nact=P(aloc,idm)->nloc;
		if(NX(HB)!=NX(aloc)){
			error("HB: Invalid format\n");
		}
		for(int imod=0; imod<NY(HB); imod++){
			dmat* hb=P(HB, idm, imod);
			if(NX(hb)!=P(aloc,idm)->nloc){
				error("hb: Invalid format\n");
			}
			dmat* phb=hb;
			for(int iact=0; iact<nact; iact++){
				if(P(P(dead,idm),iact)){
					for(int iy=0; iy<NY(hb); iy++){
						P(phb, iact, iy)=0;
					}
				}
			}
		}
	}
}

/**
   When some actuators are float, remove the corresponding column in HA and/or HB,
   and add to neighboring actuators. This is implemented using a second matrix and
   then add to the original matrix.*/
void act_float(loccell* aloc, 	///<coordinate of actuators
	dspcell** HA, 				///<sparse matrix to modify
	const dcell* HB, 			///<dense matrix to modify
	const lcell* actfloat		///<floating actuator mask
	){
	if(!actfloat||((!HA||!*HA)&&!HB)){
		warning("Nothing to do\n");
		return;
	}
	int ndm=NX(actfloat);
	dspcell* dHA=NULL;
	if(HA&&*HA){
		int nfit=(*HA)->nx;
		dHA=dspcellnew(nfit, ndm);
	}
	for(int idm=0; idm<ndm; idm++){
		if(!P(actfloat,idm)) continue;
		loc_create_map(P(aloc,idm));
		map_t* map=P(aloc,idm)->map;
		real ox=map->ox;
		real oy=map->oy;
		real dx1=1./P(aloc,idm)->dx;
		real dy1=1./P(aloc,idm)->dy;
		const real* locx=P(aloc,idm)->locx;
		const real* locy=P(aloc,idm)->locy;
		long nact=P(aloc,idm)->nloc;
		/*which floating actuator to assign for this one. */
		lmat* indfloat=lnew(4, nact);
		/*number of floating actuators that is assigned for this one. */
		long* nindfloat=mycalloc(nact, long);
		/*active neighbors of each dead act. */
		long* neighbor=mycalloc(nact, long);
		long nzmax=0;
		long* isfloat=P(P(actfloat,idm));
		/*loop through the floating actuators */
		for(int iact=0; iact<nact; iact++){
			if(!isfloat[iact]) continue;
			//warning("DM %d act %d is floating\n", idm, iact);
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
				spint* pp=P(pHA, ifit, idm)->pp;
				spint* pi=P(pHA, ifit, idm)->pi;
				real* px=P(pHA, ifit, idm)->px;

				P(pdHA, ifit, idm)=dspnew(P(pHA, ifit, idm)->nx, P(pHA, ifit, idm)->ny, nzmax);
				spint* pp2=P(pdHA, ifit, idm)->pp;
				spint* pi2=P(pdHA, ifit, idm)->pi;
				real* px2=P(pdHA, ifit, idm)->px;
				long count=0;
				for(long iact=0; iact<nact; iact++){
					pp2[iact]=count;
					if(nindfloat[iact]){//active actuator with a float neighbor
						for(long in=0; in<nindfloat[iact]; in++){
							long jact=P(indfloat, in, iact);/*the floating act. */
							real scale=1./neighbor[jact];
							for(long ie=pp[jact]; ie<pp[jact+1]; ie++){
								if(count>=nzmax){
									nzmax*=2;
									dspsetnzmax(P(pdHA, ifit, idm), nzmax);
									pp2=P(pdHA, ifit, idm)->pp;
									pi2=P(pdHA, ifit, idm)->pi;
									px2=P(pdHA, ifit, idm)->px;
								}
								pi2[count]=pi[ie];
								px2[count]=px[ie]*scale;//move weight from floating to active actuator
								count++;
							}
						}
					}
				}
				pp2[nact]=count;
				dspsetnzmax(P(pdHA, ifit, idm), count);
				/*Remove weights of floating actuatorsf from HA. */
				for(int iact=0; iact<nact; iact++){
					if(isfloat[iact]){
						for(int ic=pp[iact]; ic<pp[iact+1]; ic++){
							px[ic]=0;
						}
					}
				}
			}/*ifit */
		} else{/*Do dense matrix */
			for(long iact=0; iact<nact; iact++){
				if(nindfloat[iact]){//active actuator with a float neighbor
					for(long in=0; in<nindfloat[iact]; in++){
						long jact=P(indfloat, in, iact);/*the floating act. */
						real scale=1./neighbor[jact];
						for(int ifit=0; ifit<NX(HB); ifit++){
							dmat* hbi=P(HB, ifit, idm);
							for(long ix=0; ix<NX(hbi); ix++){//move weight from floating to active actuator
								P(hbi, ix, iact)+=P(hbi, ix, jact)*scale;
							}
						}
					}
				}
			}
			for(long iact=0; iact<nact; iact++){
				if(isfloat[iact]){
					for(int ifit=0; ifit<NX(HB); ifit++){
						dmat* hbi=P(HB, ifit, idm);
						memset(PCOL(hbi, iact), 0, NX(hbi)*sizeof(real));
					}
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
	const int ndm=NX(aloc);
	for(int idm=0; idm<ndm; idm++){
		if(!P(stuck,idm)) continue;
		const int nact=P(aloc,idm)->nloc;
		for(int iy=0; iy<NY(adm); iy++){
			dmat* pai=P(adm, idm, iy);
			dmat* px=pai;
			assert(NX(pai)==P(aloc,idm)->nloc);
			for(int icol=0; icol<NY(pai); icol++){
				for(int iact=0; iact<nact; iact++){
					if(P(P(stuck,idm),iact)){
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
	int ndm=NX(aloc);
	dspcell* out=dspcellnew(ndm, ndm);
	for(int idm=0; idm<ndm; idm++){
		loc_create_map(P(aloc,idm));
		map_t* map=P(aloc,idm)->map;
		real ox=map->ox;
		real oy=map->oy;
		real dx1=1./P(aloc,idm)->dx;
		real dy1=1./P(aloc,idm)->dy;
		const real* locx=P(aloc,idm)->locx;
		const real* locy=P(aloc,idm)->locy;
		long nact=P(aloc,idm)->nloc;
		/*actuator that is floating */
		long* isfloat=(actfloat&&P(actfloat,idm))?P(P(actfloat,idm)):0;
		dsp* outit=dspnew(nact, nact, nact*4);
		real* px=outit->px;
		spint* pp=outit->pp;
		spint* pi=outit->pi;
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
		P(out,idm,idm)=dsptrans(outit);
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
	//TIC;tic;
	loc_create_map(aloc);
	map_t *map=aloc->map;
	const real ox=map->ox;
	const real oy=map->oy;
	const real dx1=1./aloc->dx;
	const real dy1=1./aloc->dy;
	const real *locx=aloc->locx;
	const real *locy=aloc->locy;
	long nact=aloc->nloc;
	/*dmat *actcpl2=NULL;
	if(dmin(actcplc)==0){//there are acts with no active neighbors. need to process actcplc
		dmat *actcpl1=ddup(actcplc);
		actcpl2=ddup(actcplc);
		cpl=P(actcpl1); cpl0=cpl-1;
		while(dmin(actcpl2)==0){
			dcp(&actcpl1, actcpl2);
			//warning("dmin=%g, repeat\n", dmin(actcpl2));
			for(long iact=0; iact<nact; iact++){
				for(int iy=-1; iy<2; iy++){
					for(int ix=-1; ix<2; ix++){
						if(abs(ix)+abs(iy)==2){
							continue;//skip corner
						}
						//include self and neighbor
						long mapx=(long)round((locx[iact]-ox)*dx1);
						long mapy=(long)round((locy[iact]-oy)*dy1);
						int kact1=loc_map_get(map, mapx+ix, mapy+iy);
						if(kact1>0){
							P(actcpl2,iact)+=cpl0[kact1];
						}
					}
				}
			}
		}
		//writebin(actcpl2, "actcpl2");
		dscale(actcpl2, 1./dmax(actcpl2));
		cpl=P(actcpl2);
		cpl0=cpl-1;
		dfree(actcpl1);
	}*/
	const real *cpl=P(actcplc);
	const real *cpl0=cpl-1;
	dsp *outt=NULL;
	long nmissing=0;
	int nrepeat=0;
	do{
		long count=0;
		nmissing=0;
		dmat *actcpl2=0;
		if(outt){
			dspmm(&actcpl2, outt, actcplc, "tn", 1);
			cpl=P(actcpl2);
			cpl0=cpl-1;
		}
		dsp *outit=dspnew(nact, nact, nact*5);
		real *px=outit->px;
		spint *pp=outit->pp;
		spint *pi=outit->pi;
		for(long iact=0; iact<nact; iact++){
			pp[iact]=count;
			if(cpl[iact]<thres){
				long mapx=(long)round((locx[iact]-ox)*dx1);
				long mapy=(long)round((locy[iact]-oy)*dy1);
				int count2=count;
				real sum=0;
				for(int iy=-1; iy<2; iy++){
					for(int ix=-1; ix<2; ix++){
						if(abs(ix)+abs(iy)==2){
							continue;
						}
						int kact1=loc_map_get(map, mapx+ix, mapy+iy);
						if(kact1>0 &&(cpl0[kact1]>=thres||(ix==0&&iy==0))){
							//there are a few active neighbors
							pi[count]=kact1-1;
							sum+=(px[count]=cpl0[kact1]);
							count++;
						}
					}
				}
				//info("iact=%ld, sum=%g, count=%ld\n", iact, sum, count-pp[iact]);
				if(count==count2+1){
					//warning("actuator %ld receives no data\n", iact);
					nmissing++;
				}
				if(sum>0){
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
		if(!outt){
			outt=dspref(outit);
		}else{
			dsp *tmp=dspmulsp(outt, outit, "nn");
			dspfree(outt);
			outt=tmp;
		}
		dspfree(outit);
		dfree(actcpl2);
		if(nmissing) {
			if(nrepeat==0){
				dbg("act_extrap: nmissing=%ld", nmissing);
			}else{
				dbg(" %ld", nmissing);
			}
			nrepeat++;
		}
	}while(nmissing>0 && nrepeat<20);
	if(nrepeat>0){
		dbg("\n");
	}
	dsp* out=dsptrans(outt);
	dspfree(outt);
	/*
	//New method. Test convergence
	dmat* x=dnew(NX(out), 1); dset(x, 1);
	dmat* y1=dnew(NX(out), 1);
	dmat* y2=dnew(NX(out), 1);
	dspmm(&y1, out, x, "tn", 1);
	real diff;
	count=0;
	do{
		dsp* tmp=dspmulsp(out, out, "nn");
		dspfree(out);
		out=tmp;
		dzero(y2);
		dspmm(&y2, out, x, "tn", 1);
		diff=ddiff(y1, y2);
		dcp(&y1, y2);
		count++;
	} while(diff>EPS);
	dfree(x);
	dfree(y1);
	dfree(y2);*/
	
	dspdroptol(out, 1e-12);
	//toc2("act_extrap: %ld iterations", count);
	return out;
}
/**
 * Create an extrapolator that does the following
 * 1. Reduce the information to active actuators
 * 2. Project low order modes out from the active actuators
 * 3. Project back to full vector
 * 4. Extrapolate to inactive actuators
 * 5. Add back low order modes to the full vector
 * extrap = H_extrap*H_reduce*(I-M*M^+)+M*M^
 * */

dsp *act_extrap_each(loc_t *aloc,
	const dmat *actcplc,
	const real thres
){
	//step1: reduce to active actuators
	dmat *mask=dnew(NX(actcplc), NY(actcplc));
	for(long i=0; i<PN(actcplc); i++){
		if(P(actcplc,i)>thres){
			P(mask,i)=1;
		}
	}
	dsp *Hr=dspnewdiag(PN(mask), P(mask), 1); 
	//writebin(Hr, "Hr");
	//step2: project 
	dmat *Mz=zernike(aloc, loc_diam(aloc), 1, 3, 0);//Modes
	dmat *Md=dpinv(Mz, mask);//projector
	//writebin(mask, "mask");
	dfree(mask);
	dmat *MMp=NULL;//MMp=1-Mz*Md. //removal projector
	dmm(&MMp, 0, Mz, Md, "nn", -1);
	//writebin(MMp, "MMp");
	//writebin(Md, "Md");
	//writebin(Mz, "Mz");
	dfree(Mz);
	dfree(Md);
	daddI(MMp, 1); 
	dmat *MMr=NULL;//MMr=Hr*MMp
	dspmm(&MMr, Hr, MMp, "nn", 1);
	//writebin(MMr, "MMr");
	cellfree(Hr);
	dsp *He=act_extrap_do(aloc, actcplc, thres);
	//writebin(He, "He");
	dmat *MMe=NULL;//MMe=He*MMr+MMp
	dspmm(&MMe, He, MMr, "nn", 1);
	dspfree(He);
	//writebin(MMe, "MMe");
	daddI(MMp, -1);
	dadd(&MMe, 1, MMp, -1); //+Mz*Md
	//writebin(MMp, "MMp2");
	//writebin(MMe, "MMe2");
	cellfree(MMp);
	cellfree(MMr);
	cellfree(Hr);
	dsp *res=d2sp(MMe, 1e-12);
	dfree(MMe);
	//writebin(res, "MMe3");
	return res;
}

/**
   Create an extrapolator that make inactive actuators equal avearage of active
   neighbors if exist or all other neighbors.
*/
dspcell* act_extrap(loccell* aloc,     /**<[in] Actuator grid array*/
	const dcell* actcplc,/**<[in] Actuator coupling coefficiency*/
	real thres, /**<[in] Threshold of coupling to turn on interpolation*/
	int lor /**<[in] Low order mode removal before extrapolation*/ 
){
	int ndm=NX(actcplc);
	dspcell* out=dspcellnew(ndm, ndm);
	for(int idm=0; idm<ndm; idm++){
		if(lor){
			//when enabled, the resulting matrix is much less sparse.
			P(out,idm,idm)=act_extrap_each(P(aloc,idm), P(actcplc,idm), thres);
		}else{
			P(out, idm, idm)=act_extrap_do(P(aloc, idm), P(actcplc, idm), thres);
		}
	}
	return out;
}

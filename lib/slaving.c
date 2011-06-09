#include "bin.h"
#include "loc.h"
#include "dmat.h"
#include "dsp.h"
#include "mathmisc.h"
#include "slaving.h"
#include "imat.h"
/**
   Compute slaving actuator regularization. HA (or use GA) is used to compute
   active actuators. If NW is non NULL, orthogonalize it with the slaving
   regularization.  When the actuators are in the NULL space of HA, we want to
   contraint their values to be close to the ones that are active. We put an
   additional term in the fitting matrix to force this. Be careful with it when
   using tip/tilt constraint and cholesky back solve.  */
spcell *slaving(loc_t **aloc,  /**<[in]The actuator grid*/
		spcell *HA,    /**<[in]The influence function from actuator to destination*/
		dmat *W1,      /**<[in]The weighting function on the destination grid. (optional)*/
		dcell *NW,     /**<[in]The low rank terms that need to be orthogonal to the output (optional)*/
		icell *actstuck,/**<[in]List of stuck actuators that will not be slaved, but have value constrained.*/
		icell *actfloat,/**<[in]List of stuck actuators that will be slaved, but not have value constrained.*/
		double thres,  /**<[in]The threshold that an actuator is deemed slave*/
		double scl     /**<[in]The scaling of the overal value*/
		){
    if(!HA) {
	error("HA is not supplied\n");
    }
    if(scl<EPS){
	error("scl=%g is too small\n", scl);
    }
    PSPCELL(HA,pHA);
    int ndm=HA->ny;
    dcell *actcplc=dcellnew(ndm, 1);
    spcell *actslavec=spcellnew(ndm, ndm);//block diagonal.
    PSPCELL(actslavec, actslave);
    int nslavetot=0;
    /*
      We first sum the absolute value of the weights for each actuator to figure
      out how much it is coupled to destination.
    */
    for(int idm=0; idm<ndm; idm++){
	int nact=aloc[idm]->nloc;
	for(int ifit=0; ifit<HA->nx; ifit++){
	    if(!pHA[idm][ifit]) continue;
	    if(W1){
		sptmulmat(&actcplc->p[idm], pHA[idm][ifit], W1, 1);
	    }else{
		dmat *tmp=spsumabs(pHA[idm][ifit], 1);
		dadd(&actcplc->p[idm], 1, tmp, 1);
		dfree(tmp);
	    }
	}
	if(actcplc->p[idm]->nx*actcplc->p[idm]->ny!=nact){
	    error("Invalid actcplc\n");
	}
	normalize_max(actcplc->p[idm]->p, nact, 1);//bring max to 1;
	int *stuck=calloc(nact, sizeof(int));
	int *floated=calloc(nact, sizeof(int));
	int nstuck=0;
	int nfloat=0;
	if(actstuck && actstuck->p[idm]){
	    nstuck=actstuck->p[idm]->nx;
	    for(int jact=0; jact<nstuck; jact++){
		int iact=actstuck->p[idm]->p[jact];
		stuck[iact]=1;
		actcplc->p[idm]->p[iact] = 1;//Skip the stuck actuators.
		warning2("Skip actuator %d in slaving\n", iact);
	    }
	}
	if(actfloat && actfloat->p[idm]){
	    nfloat=actfloat->p[idm]->nx;
	    for(int jact=0; jact<nfloat; jact++){
		int iact=actfloat->p[idm]->p[jact];
		floated[iact]=1;
		if(actcplc->p[idm]->p[iact]>0.1){
		    error("actuator %d is floating, but actcpl is not zero\n", iact);
		}
		warning2("Don't constrain value of actuator %d in slaving\n", iact);
	    }
	}
  	double *actcpl= actcplc->p[idm]->p;
	double *actcpl0 = actcpl-1;
	int  nslave   = 0;
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
	loc_create_map_npad(aloc[idm],1);
	long (*map)[aloc[idm]->map->nx]=(void*)aloc[idm]->map->p;
	double ox=aloc[idm]->map->ox;
	double oy=aloc[idm]->map->oy;
	double dx1=1./aloc[idm]->dx;
	dsp *slavet=spnew(nact,nslave,nslave*5);
	spint *pp=slavet->p;
	spint *pi=slavet->i;
	double *px=slavet->x;
	const double *locx=aloc[idm]->locx;
	const double *locy=aloc[idm]->locy;
	long count=0;
	long icol=0;
	for(int iact=0; iact<nact; iact++){
	    if(stuck && stuck[iact]){//limit the strength of stuck actuators.
		pp[icol]=count;
		pi[count]=iact;
		px[count]=scl;
		count++;
		icol++;
	    }else if(actcpl[iact]<thres){//slave actuators
		pp[icol]=count;
		long mapx=(long)round((locx[iact]-ox)*dx1);
		long mapy=(long)round((locy[iact]-oy)*dx1);
		assert(map[mapy][mapx]-1==iact);
		int near_active=0;
		int near_exist=0;
		for(int idy=-1; idy<2; idy++){
		    for(int idx=-1; idx<2; idx++){
			if((idx!=0 && idy!=0) || (idx==0 && idy==0)){
			    continue;//skip center and corner
			}
			int kact1=map[mapy+idy][mapx+idx];
			if(kact1 && !stuck[kact1-1]){
			    near_exist++;
			    if(actcpl0[kact1]>0.1){
				near_active++;
			    }
			}
		
		    }
		}
		if(!near_exist){
		    error("This is an isolated actuator\n");
		}
		/*10x bigger constraint make the slaving work more strictly for floating actuators.*/
		pi[count]=iact;
		px[count]=scl;
		/*
		  We limits the strength of these slaving actuators.
		*/
		
		if(actcpl[iact]<0.1 && !floated[iact]){
		    px[count]+=scl*1;
		}else if(actcpl[iact]<0.1 && floated[iact]){
		    warning2("Actuator %d is floating, don't limit its strength\n", iact);
		}
		count++;

		if(near_active>0){
		    /*
		      neighbors are defined as the four pixels to the left, right, top and bottom.
		      If some of the neighbors are active, use the average of them for my value.
		     */
		    double value=-scl/near_active;
		    for(int idy=-1; idy<2; idy++){
			for(int idx=-1; idx<2; idx++){
			    if((idx!=0 && idy!=0) || (idx==0 && idy==0)){
				continue;//skip center and corner
			    }
			    int kact1=map[mapy+idy][mapx+idx];
			    if(kact1 && actcpl0[kact1]>0.1 && !stuck[kact1-1]){
				pi[count]=kact1-1;
				px[count]=value;
				count++;
			    }
			}
		    }
		}else{
		    /*
		      If none of the neighbors are active, use the average of
		      all neighbors for my value.
		     */
		    double value=-scl/near_exist;
		    for(int idy=-1; idy<2; idy++){
			for(int idx=-1; idx<2; idx++){
			    if((idx!=0 && idy!=0) || (idx==0 && idy==0)){
				continue;//skip center and corner
			    }
			    int kact1=map[mapy+idy][mapx+idx];
			    if(kact1 && !stuck[kact1-1]){
				pi[count]=kact1-1;
				px[count]=value;
				count++;
			    }
			}
		    }
	
		}
		icol++;
	    }//if
	}//for iact
	free(stuck);
	free(floated);
	pp[icol]=count;
	assert(icol==nslave);
	spsetnzmax(slavet, count);

	dsp *slave=sptrans(slavet);
	spwrite(slave, "slave_%d", idm);
	actslave[idm][idm]=spmulsp(slavet, slave);

	if(NW){
	    /*Now we need to make sure NW is in the NULL
	      space of the slaving regularization, especially
	      the tip/tilt constraints.*/
	    if(NW->p[idm]){
		dmat *H=NULL;
		spfull(&H, slavet, 1);
		dmat *Hinv=dpinv(H,NULL,NULL);
		dmat *mod=NULL;
		dmm(&mod, Hinv, NW->p[idm], "nn", 1);
		dmm(&NW->p[idm], H, mod,"nn", -1);
		dfree(H);
		dfree(Hinv);
		dfree(mod);
		if(nfloat || nstuck){//Remove corresponding rows in NW
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
    }//idm
    dcellfree(actcplc);
    if(nslavetot==0){
	spcellfree(actslavec);
	actslavec=NULL;
    }
    return actslavec;
}
/**
   When some actuators are stuck, remove the corresponding column in HA or GA.
*/
void act_stuck(loc_t **aloc,
	       spcell *HA,
	       icell *stuck){
    if(!stuck) return;
    PSPCELL(HA,pHA);
    int ndm=HA->ny;
    for(int idm=0; idm<ndm; idm++){
	if(!stuck->p[idm]){
	    continue;
	}
	for(int ifit=0; ifit<HA->nx; ifit++){
	    spint *pp=pHA[idm][ifit]->p;
	    double *px=pHA[idm][ifit]->x;
	    for(int jact=0; jact<stuck->p[idm]->nx; jact++){
		int iact=stuck->p[idm]->p[jact];
		for(int ic=pp[iact]; ic<pp[iact+1]; ic++){
		    px[ic]=0;
		}
	    }
	}
    }
}
/**
   When some actuators are float, remove the corresponding column in HA or GA,
and add to neigh boring actuators. This is implemented using a second matrix and
do matrix addition.*/
void act_float(loc_t **aloc,
	       spcell **HA,
	       icell *actfloat){
    if(!actfloat || !HA || !*HA) return;
    int ndm=(*HA)->ny;
    int nfit=(*HA)->nx;
    PSPCELL(*HA,pHA);
    spcell *dHA=spcellnew(nfit,ndm);
    PSPCELL(dHA,pdHA);
    for(int idm=0; idm<ndm; idm++){
	if(!actfloat->p[idm]) continue;
	loc_create_map_npad(aloc[idm],1);
	long (*map)[aloc[idm]->map->nx]=(void*)aloc[idm]->map->p;
	double ox=aloc[idm]->map->ox;
	double oy=aloc[idm]->map->oy;
	double dx1=1./aloc[idm]->dx;
	const double *locx=aloc[idm]->locx;
	const double *locy=aloc[idm]->locy;
	long nact=aloc[idm]->nloc;
	long ndead=actfloat->p[idm]->nx;
	//actuator that is floating
	long *isfloat=calloc(nact, sizeof(long));
	//which floating actuator to assign for this one.
	long (*indfloat)[4]=calloc(nact, sizeof(long)*4);
	//number of floating actuators that is assigned for this one.
	long *nindfloat=calloc(nact, sizeof(long));
	//active neighbors of each dead act.
	long *neighbor=calloc(nact, sizeof(long));
	long nzmax=0;
	
	for(long jact=0; jact<ndead; jact++){
	    int iact=actfloat->p[idm]->p[jact];
	    isfloat[iact]=1;
	}
	//loop through the floating actuators
	for(long jact=0; jact<ndead; jact++){
	    int iact=actfloat->p[idm]->p[jact];
	    long mapx=(long)round((locx[iact]-ox)*dx1);
	    long mapy=(long)round((locy[iact]-oy)*dx1);
	    //find all its neighbors.
	    for(int idy=-1; idy<2; idy++){
		for(int idx=-1; idx<2; idx++){
		    if((idx!=0 && idy!=0) || (idx==0 && idy==0)){
			continue;//skip center and corner
		    }
		    if(map[mapy+idy][mapx+idx]){
			int kact=map[mapy+idy][mapx+idx]-1;
			if(!isfloat[kact]){
			    indfloat[kact][nindfloat[kact]]=iact;
			    nindfloat[kact]++;
			    neighbor[iact]++;
			    nzmax+=4;//assume 1 point couples to 4 points max.
			}
		    }
		}
	    }
	}
	//Create dHA to assign weights of floating actuators to neighbors.
	for(int ifit=0; ifit<nfit; ifit++){
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
			long jact=indfloat[iact][in];//the floating act.
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
	    //Remove weights of floating actuatorsf from HA.
	    for(int jact=0; jact<actfloat->p[idm]->nx; jact++){
		int iact=actfloat->p[idm]->p[jact];
		for(int ic=pp[iact]; ic<pp[iact+1]; ic++){
		    px[ic]=0;
		}
	    }
	}//ifit
	free(nindfloat);
	free(indfloat);
	free(isfloat);
	free(neighbor);	
    }//idm
    spcelladd(HA, dHA);
    spcellfree(dHA);
}

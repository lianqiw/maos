#include "bin.h"
#include "loc.h"
#include "dmat.h"
#include "dsp.h"
#include "mathmisc.h"
#include "slaving.h"
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
		double thres,  /**<[in]The threshold that an actuator is deemed slave*/
		double scl     /**<[in]The scaling of the overal value*/
		){
    if(!HA) {
	error("HA is not supplied\n");
    }
    PSPCELL(HA,pHA);
    int ndm=HA->ny;
    dcell *actcplc=dcellnew(ndm, 1);
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
    }
    int nslavetot=0;
    spcell *actslavec=spcellnew(ndm, ndm);//block diagonal.
    PSPCELL(actslavec, actslave);
    
    for(int idm=0; idm<ndm; idm++){
	int  nact     = aloc[idm]->nloc;
	double *actcpl= actcplc->p[idm]->p;
	double *actcpl0 = actcpl-1;
	int  nslave   = 0;
	for(int iact=0; iact<nact; iact++){
	    if(actcpl[iact]<thres){
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
	    if(actcpl[iact]<thres){//slave actuators
		pp[icol]=count;
		long mapx=(long)round((locx[iact]-ox)*dx1);
		long mapy=(long)round((locy[iact]-oy)*dx1);
		if(map[mapy][mapx]-1!=iact){
		    error("mapping is used incorrectly\n");
		}
		int near_active=0;
		int near_exist=0;
		for(int idy=-1; idy<2; idy++){
		    for(int idx=-1; idx<2; idx++){
			if((idx!=0 && idy!=0) || (idx==0 && idy==0)){
			    continue;//skip center and corner
			}
			if(map[mapy+idy][mapx+idx]){
			    near_exist++;
			    if(actcpl0[map[mapy+idy][mapx+idx]]>0.1){
				near_active++;
			    }
			}
		
		    }
		}
		if(!near_exist){
		    error("This is an isolated actuator\n");
		}
		pi[count]=iact;
		px[count]=scl;
		/*
		  We limits the strength of these slaving actuators.
		*/
		if(actcpl[iact]<0.1){
		    px[count]+=scl*1;
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
			    if(map[mapy+idy][mapx+idx] && actcpl0[map[mapy+idy][mapx+idx]]>0.1){
				pi[count]=map[mapy+idy][mapx+idx]-1;
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
			    if(map[mapy+idy][mapx+idx]){
				pi[count]=map[mapy+idy][mapx+idx]-1;
				px[count]=value;
				count++;
			    }
			}
		    }
	
		}
		icol++;
	    }
	}
	pp[icol]=count;
	if(icol!=nslave){
	    error("Doesnot match\n");
	}
	spsetnzmax(slavet, count);
	//spwrite(slavet,"slavet");
	loc_free_map(aloc[idm]);
	dsp *slave=sptrans(slavet);
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

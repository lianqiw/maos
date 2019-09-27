/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
  Included by accphi.c for prop_grid_map and prop_grid_map_tranpose
  derived from prop_grid_stat.c on 2011-10-05.
*/

#undef FUN_NAME_MAP
#undef FUN_NAME_PTS
#undef FUN_NAME_STAT
#undef FUN_NAME_BLOCK
#undef CONST_IN
#undef CONST_OUT

#if TRANSPOSE == 0
#define CONST_IN const
#define CONST_OUT 
#define FUN_NAME_MAP prop_grid_map
#define FUN_NAME_PTS prop_grid_pts
#define FUN_NAME_STAT prop_grid_stat
#define FUN_NAME_BLOCK prop_grid_block
#else
#define CONST_IN 
#define CONST_OUT const
#define FUN_NAME_MAP prop_grid_map_transpose
#define FUN_NAME_PTS prop_grid_pts_transpose
#define FUN_NAME_STAT prop_grid_stat_transpose
#define FUN_NAME_BLOCK prop_grid_block_transpose
#endif

/**
   Ray tracing from phiin with size of nxin*nyin to 
   phiout with size of nxout*nyout, normalized spacing of dxout, dyout, origin offset of oxout, oyout.
 */
static int 
FUN_NAME_BLOCK(CONST_IN double *phiin, long nxin, long nyin,
	       CONST_OUT double *phiout, long nxout, long nyout, 
	       double dxout, double dyout, double oxout, double oyout,
	       const double alpha, int wrap){
    
    const long wrapx  = nxin - 1;
    const long wrapy  = nyin - 1;
    int missing=0;
    if(PROP_GRID_MAP_OPTIM){
	if(fabs(dxout-1.)<EPS && fabs(dyout-1.)<EPS){
	    //loc_out and loc_in has the same sampling.
	    int irows=0;
	    double dplocx=oxout;
	    double dplocy=oyout;
	    int nplocx, nplocy0;
	    if(wrap){
		dplocx-=floor(dplocx/(double)nxin)*nxin;
	    }else if(dplocx<0){
		irows=iceil(-dplocx);
		missing+=irows*nyout;
	    }
	    SPLIT(dplocx,dplocx,nplocx);

	    const int do_edge=dplocx<EPS || wrap;
	    int rowdiv0=wrapx-nplocx;
	    if(rowdiv0>nxout) rowdiv0=nxout;
	    if(rowdiv0<0) rowdiv0=0;

	    int icols=0;
	    if(wrap){
		dplocy-=floor(dplocy/(double)nyin)*nyin;
	    }else{
		//limit column range to within [0, wrapy]
		if(dplocy<0){
		    icols=iceil(-dplocy);
		    missing+=icols*nxout;
		}
		int coldiv=ifloor(wrapy-oyout)+1;
		if(coldiv<nyout){
		    missing+=(nyout-coldiv)*nxout;
		    nyout=coldiv;
		}
	    }
	    SPLIT(dplocy,dplocy,nplocy0);
	    const double bl=dplocx*dplocy;
	    const double br=(1.-dplocx)*dplocy;
	    const double tl=dplocx*(1-dplocy);
	    const double tr=(1.-dplocx)*(1-dplocy);
	    phiin+=nplocx;
	    if(nyout==640){
		dbg("icols=%d, nyout=%ld\n", icols, nyout);
	    }
	    OMPTASK_FOR(icol, icols, nyout){
		CONST_IN double *restrict phicol, *restrict phicol2;
		CONST_OUT double *restrict phiout2=phiout+icol*nxout;//starting address of output
		long nplocy=nplocy0+icol;
		while(nplocy>=nyin){
		    nplocy-=nyin;
		}
		//The points on [wrapy nyout) is handled correctly
		phicol  = phiin+nplocy*nxin;
		phicol2 = phiin+(nplocy==wrapy?0:(nplocy+1))*nxin;
#if TRANSPOSE == 0				
#define GRID_ADD(irow, offset)			\
		phiout2[irow]+=alpha*		\
		    (bl*phicol2[irow+offset]	\
		     +br*phicol2[irow]		\
		     +tl*phicol[irow+offset]	\
		     +tr*phicol[irow]);
#else
#define GRID_ADD(irow,offset)			\
		double tmp=alpha*phiout2[irow];	\
		phicol2[irow+offset]+=tmp*bl;	\
		phicol2[irow]+=tmp*br;		\
		phicol[irow+offset]+=tmp*tl;	\
		phicol[irow]+=tmp*tr;
#endif
		int irow=irows;
		int rowdiv=rowdiv0;
		do{
		    //First handle points fall within [0, wrapx)
		    for(; irow<rowdiv; irow++){
			GRID_ADD(irow, 1);
		    }
		    //Check if point in [wrapx, nxin) if wrap or  [wrapx, wrapx] if not wrap
		    if(irow<nxout && do_edge){
			GRID_ADD(irow, -wrapx);
			irow++;
		    }
		    //wrap around the data
		    rowdiv=irow+wrapx;
		    phicol2-=nxin;
		    phicol-=nxin;
		    if(rowdiv>nxout) rowdiv=nxout;
		}while(wrap && irow<nxout);
#undef GRID_ADD
	    }/*end for icol*/
	    OMPTASK_END;
	}else{
	    //grid size of loc_in and loc_out doesn't agree
	    assert(dyout>0);
	    const double dxout1=1./dxout;
	    const double dyout1=1./dyout;
	    /*Figure out the array offset. icols is start of y. irows is start of x. */
	    long icols=0;
	    if(wrap){//wrap starting point to within the active region
		oyout-=floor(oyout/(double)nyin)*nyin;
	    }else{
		//limit column range to within [0, wrapy]
		if(oyout<0){
		    icols=iceil(-oyout*dyout1);
		}
		//ifloor()+1 gives [0, wrapy]. 
		//iceil() gives [0, wrapy)
		int coldiv=ifloor((wrapy-oyout)*dyout1)+1;
		if(coldiv<nyout){
		    nyout=coldiv;
		}
	
	    }
	    int irows=0;
	    if(wrap){
		oxout-=floor(oxout/(double)nxin)*nxin;
	    }else if(oxout<0){
		irows=iceil(-oxout*dxout1);
	    }  
	    //rowdiv is number of points falls within [0, wrapx)
	    int rowdiv0 = iceil((wrapx-oxout)*dxout1);
	    if(rowdiv0<=0){
		if(!wrap) {
		    missing=nxout*nyout;
		    return missing;
		}else{
		    rowdiv0=0;
		}
	    }else if(rowdiv0>nxout){
		rowdiv0=nxout;
	    }
	    const double dplocx00 = oxout+irows*dxout;
	    //allow points falls between [wrapx, wrapx+xover)
	    const double xover=wrap?1:EPS;
	    double *dplocx_arr=0;
	    int *nplocx_arr=0;
	    int *nplocx2_arr=0;
	    if(nyout>10 && nxout<2000){//cache SPLIT() results
		//alloca is not good in mac due to stack limitation.
		dplocx_arr=mymalloc(nxout,double);
		nplocx_arr=mymalloc(nxout,int);
		nplocx2_arr=mymalloc(nxout,int);

		double rowdiv=rowdiv0;
		double dplocx=dplocx00;
		int irow=irows;
		double dplocx0;
		int nplocx0;
		do{
		    //First handle points fall within [0, wrapx)
		    for(; irow<rowdiv; irow++, dplocx+=dxout){
			SPLIT(dplocx,dplocx0,nplocx0);
			dplocx_arr[irow]=dplocx0;
			nplocx_arr[irow]=nplocx0;
			nplocx2_arr[irow]=1;
		    }
		    //Then handle points fall within [wrapx, wrapx+xover)
		    //xover is EPS if wrap==0, 1 if wrap==1.
		    for(; irow < nxout && dplocx<wrapx+xover; irow++, dplocx+=dxout){//right on edge
			SPLIT(dplocx,dplocx0,nplocx0);
			dplocx_arr[irow]=dplocx0;
			nplocx_arr[irow]=nplocx0;
			nplocx2_arr[irow]=-wrapx;
		    }
		    //wrap the points around;
		    dplocx-=nxin;			
		    rowdiv=irow+iceil((wrapx-dplocx)*dxout1);
		    if(rowdiv>nxout) rowdiv=nxout;
		}while(wrap && irow<nxout);
	    }
	    for(long icol=icols; icol<nyout; icol++){
		CONST_OUT double *phiout2=phiout+icol*nxout;
		double dplocy=myfma(icol, dyout, oyout);

		CONST_IN double *phicol, *phicol2;
		double dplocx,dplocx0,dplocy0;
		int nplocx0,nplocy;
		int rowdiv;
		while(dplocy>=nyin){
		    //only happens when wrap=1. wrap to within [0, nyin)
		    dplocy-=nyin;
		    oyout-=nyin;
		}
		SPLIT(dplocy,dplocy0,nplocy);
		phicol  = phiin+nplocy*nxin;
                //The points on [wrapy nyout) is handled correctly
		phicol2 = phiin+(nplocy==wrapy?0:(nplocy+1))*nxin;
		const double dplocy1=1.-dplocy0;

#if TRANSPOSE == 0
#define GRID_ADD(irow,offset)						\
		phiout2[irow]+=alpha*					\
		    (+dplocy1*((1-dplocx0)*phicol[nplocx0]		\
			       +dplocx0*phicol[nplocx0+offset])		\
		     +dplocy0*((1-dplocx0)*phicol2[nplocx0]		\
			       +dplocx0*phicol2[nplocx0+offset]));
#else
#define GRID_ADD(irow,offset)					\
		double tmp=phiout2[irow]*alpha;			\
		phicol[nplocx0]+=tmp*dplocy1*(1-dplocx0);	\
		phicol[nplocx0+offset]+=tmp*dplocy1*dplocx0;	\
		phicol2[nplocx0]+=tmp*dplocy0*(1-dplocx0);	\
		phicol2[nplocx0+offset]+=tmp*dplocy0*dplocx0;
#endif

		rowdiv=rowdiv0;
		dplocx=dplocx00;

		if(dplocx_arr){//use cached results
		    //First handle points fall within [0, wrapx)
#pragma omp simd
		    for(int irow=irows; irow<rowdiv; irow++){
			dplocx0=dplocx_arr[irow];
			nplocx0=nplocx_arr[irow];
			GRID_ADD(irow, (nplocx2_arr[irow]));
		    }
		    //Then handle points fall within [wrapx, wrapx+xover)
		    //xover is EPS if wrap==0, 1 if wrap==1.
		    if(wrap){
#pragma omp simd
			for(int irow=rowdiv; irow < nxout ; irow++){//right on edge
			    dplocx0=dplocx_arr[irow];
			    nplocx0=nplocx_arr[irow];
			    GRID_ADD(irow, (nplocx2_arr[irow]));
			}
		    }
		}else{
		    int irow=irows;
		    do{
			//First handle points fall within [0, wrapx)
			for(; irow<rowdiv; irow++, dplocx+=dxout){
			    SPLIT(dplocx,dplocx0,nplocx0);
			    GRID_ADD(irow,1);
			}
			//Then handle points fall within [wrapx, wrapx+xover)
			//xover is EPS if wrap==0, 1 if wrap==1.
			for(; irow < nxout && dplocx<wrapx+xover; irow++, dplocx+=dxout){//right on edge
			    SPLIT(dplocx,dplocx0,nplocx0);
			    GRID_ADD(irow,-wrapx);
			}
			//wrap the points around;
			dplocx-=nxin;			
			rowdiv=irow+iceil((wrapx-dplocx)*dxout1);
			if(rowdiv>nxout) rowdiv=nxout;
		    }while(wrap && irow<nxout);
		}
#undef GRID_ADD
	    }/*for icol*/
	    if(dplocx_arr){
		free(dplocx_arr);
		free(nplocx_arr);
		free(nplocx2_arr);
	    }
	}/*fabs(dx_in2*dxout-1)<EPS*/
    }else{
	//warning_once("Using unoptmized prop_grid_map\n");
	/*non optimized case. slower, but hopefully accurate*/
	OMPTASK_FOR(icol, 0, nyout){
	    double dplocy1;
	    double dplocx,dplocx0;
	    int nplocx,nplocy,nplocy1, nplocx1;
	    CONST_OUT double *phiout2=phiout+icol*nxout;
	    double dplocy=oyout+icol*dyout;
	    if(wrap){
		dplocy=dplocy-floor(dplocy/(double)nyin)*nyin;
	    }
	    if(dplocy>=0 && dplocy <=wrapy){
		SPLIT(dplocy,dplocy,nplocy);
		dplocy1=1.-dplocy;
		nplocy1=(nplocy==wrapx?0:nplocy+1);
		dplocx0=oxout;

		for(int irow=0; irow<nxout; irow++){
		    dplocx=dplocx0+irow*dxout;
		    if(wrap){
			dplocx=dplocx-floor(dplocx/(double)nxin)*nxin;
		    }
		    if(dplocx>=0 && dplocx<=wrapx){
			SPLIT(dplocx,dplocx,nplocx);
			nplocx1=(nplocx==wrapx?0:nplocx+1);
#if TRANSPOSE == 0
			phiout2[irow]+= alpha*
			    (dplocx * (dplocy * phiin[(nplocx1) + (nplocy1)*nxin]
				       +dplocy1 * phiin[(nplocx1) + nplocy*nxin])
			     + (1-dplocx) * (dplocy * phiin[nplocx + (nplocy1)*nxin]
					     +dplocy1 * phiin[nplocx + nplocy*nxin]));
#else
			double tmp=alpha*phiout2[irow];
			phiin[(nplocx1) + (nplocy1)*nxin]+=tmp*dplocx*dplocy;
			phiin[(nplocx1) + nplocy*nxin]+=tmp*dplocx*dplocy1;
			phiin[nplocx + (nplocy1)*nxin]+=tmp*(1-dplocx)*dplocy;
			phiin[nplocx + nplocy*nxin]+=tmp*(1-dplocx)*dplocy1;
#endif		
		    }else{
			missing++;
		    }
		}/*for irow*/
	    }else{
		missing+=nxout;
	    }
	}//OMPTASK_FOR;
	OMPTASK_END;
    }
    return missing;
}/*function*/



/**
   Propagate OPD defines on grid mapin to grid mapout.  alpha is the scaling of
   data. displacex, displacy is the displacement of the center of the beam on
   the input grid. scale is the cone effect.
*/

/*
  Implement prop_grid_stat or prop_grid_stat_tranpose in a unified way. 
  Let mapin be x, phiout be y. H be propagator.
  prop_grid_stat does y+=H*x;
  prop_grid_stat_transpose does x+=H'*y;
 
*/
void FUN_NAME_MAP (CONST_IN map_t *mapin,   /**<[in] OPD defind on a square grid*/
		   CONST_OUT map_t *mapout, /**<[in,out] output OPD defined on a square grid*/
		   ARG_PROP,
		   int wrap,           /**<[in] wrap input OPD or not*/
		   long colstart,      /**<[in] First column to do ray tracing*/
		   long colend         /**<[in] Last column (exclusive) to do ray tracing*/
    ){
    if(mapin->iac){
#if TRANSPOSE == 0
	if(wrap) error("wrap=1 is invalid for cubic\n");
	prop_grid_map_cubic(mapin, mapout, ARG_PROP2, mapin->iac, colstart, colend);
	return;
#else
	error("transpose ray tracing is not available with iac\n");
#endif
    }
	
    CONST_OUT double *phiout=mapout->p;
    CONST_IN double *phiin  = mapin->p;
    /*With OpenMP compiler complained uninitialized value for the following
      because they are treated as firstprivate*/
    if(colend==0) colend = mapout->ny;
    const long nxin = mapin->nx;
    const long nyin = mapin->ny;
    /*
      convert to unitless variables
    */
    const double dx_in1 = 1./mapin->dx;
    const double dy_in1 = 1./mapin->dy;
    const double dxout  = mapout->dx*dx_in1*scale;
    const double dyout  = mapout->dy*dy_in1*scale;
    displacex = (displacex-mapin->ox)*dx_in1;
    displacey = (displacey-mapin->oy)*dy_in1;
    const double oxout=mapout->ox*dx_in1*scale+displacex;
    const double oyout=mapout->oy*dy_in1*scale+displacey;
    const long nxout=mapout->nx;
    int missing=FUN_NAME_BLOCK(phiin, nxin, nyin, 
			       phiout+colstart*nxout, nxout, colend-colstart,
			       dxout, dyout, oxout, oyout+dyout*colstart, 
			       alpha, wrap);
    WARN_MISSING;
}

void FUN_NAME_PTS(CONST_IN map_t *mapin, /**<[in] OPD defind on a square grid*/
		  const pts_t *pts,/**<[in] coordinate of destination grid*/
		  CONST_OUT double *phiout, /**<[in,out] OPD defined on locout*/
		  ARG_PROP,
		  int wrap,           /**<[in] wrap input OPD or not*/
		  long sastart,       /**<[in] The starting subaperture to trace ray*/
		  long saend          /**<[in] The last (exclusive) subaperture to trace ray*/
    ){
    if(mapin->iac){
#if TRANSPOSE == 0
	if(wrap) error("wrap=1 is invalid for cubic\n");
	prop_grid_pts_cubic(mapin, pts, phiout, ARG_PROP2, mapin->iac, sastart, saend);
	return;
#else
	error("transpose ray tracing is not available with iac\n");
#endif
    }
    const long nxin = mapin->nx;
    const long nyin = mapin->ny;
    /*
      convert to unitless variables
    */
    const double dx_in1 = 1./mapin->dx;
    const double dy_in1 = 1./mapin->dy;
    const double dxout  = pts->dx*dx_in1*scale;
    const double dyout  = pts->dy*dy_in1*scale;
    displacex = (displacex-mapin->ox)*dx_in1;
    displacey = (displacey-mapin->oy)*dy_in1;
    long nxout=pts->nx;
    long nyout=pts->ny?pts->ny:pts->nx;
    if(!saend) saend=pts->nsa;
    int missing=0;
    OMPTASK_FOR(isa, sastart, saend,shared(missing)){
	//for(long isa=sastart; isa<saend; isa++){
	const double oxout=pts->origx[isa]*dx_in1*scale+displacex;
	const double oyout=pts->origy[isa]*dy_in1*scale+displacey;
	missing+=FUN_NAME_BLOCK(mapin->p, nxin, nyin, 
				phiout+isa*nxout*nyout, nxout, nyout,
				dxout, dyout, oxout, oyout, 
				alpha,wrap);
    }
    OMPTASK_END;
    WARN_MISSING;
}

void FUN_NAME_STAT (CONST_IN map_t *mapin, /**<[in] OPD defind on a square grid*/
		    const locstat_t *ostat, /**<[in] information about each clumn of a output loc grid*/
		    CONST_OUT double *phiout,    /**<[in,out] OPD defined on ostat*/
		    ARG_PROP,
		    int wrap,           /**<[in] wrap input OPD or not*/
		    long colstart,      /**<[in] First column to do ray tracing*/
		    long colend         /**<[in] Last column (exclusive) to do ray tracing*/
    ){
    if(mapin->iac){
#if TRANSPOSE == 0
	prop_grid_stat_cubic(mapin, ostat, phiout, ARG_PROP2, mapin->iac, colstart, colend);
	return;
#else
	error("transpose ray tracing is not available with iac\n");
#endif
    }
    const long nxin = mapin->nx;
    const long nyin = mapin->ny;
    /*
      convert to unitless variables
    */
    const double dx_in1 = 1./mapin->dx;
    const double dy_in1 = 1./mapin->dy;
    const double dxout  = ostat->dx*dx_in1*scale;
    const double dyout  = ostat->dy*dy_in1*scale;
    displacex = (displacex-mapin->ox)*dx_in1;
    displacey = (displacey-mapin->oy)*dy_in1;
    const long nyout=1;
    int missing=0;
    if(colend==0) colend = ostat->ncol;
    OMPTASK_FOR(icol, colstart, colend,shared(missing)){
	const long offset=ostat->cols[icol].pos;
	const long nxout=ostat->cols[icol+1].pos-offset;
	const double oxout=ostat->cols[icol].xstart*dx_in1*scale+displacex;
	const double oyout=ostat->cols[icol].ystart*dy_in1*scale+displacey;
	missing+=FUN_NAME_BLOCK(mapin->p, nxin, nyin, 
				phiout+offset, nxout, nyout,
				dxout, dyout, oxout, oyout, 
				alpha, wrap);
    }
    OMPTASK_END;
    WARN_MISSING;
}

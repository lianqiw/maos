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
/*
  Included by accphi.c for prop_grid_stat and prop_grid_stat_tranpose
  Implemented and test on 2011-04-28
*/
#undef CONST_IN
#undef CONST_OUT
#undef FUN_NAME
#if TRANSPOSE == 0
#define CONST_IN const
#define CONST_OUT 
#define FUN_NAME prop_grid_stat
#else
#define CONST_IN 
#define CONST_OUT const
#define FUN_NAME prop_grid_stat_transpose
#endif
/*
  Implement prop_grid_stat or prop_grid_stat_tranpose in a unified way. 
  Let mapin be x, phiout be y. H be propagator.
  prop_grid_stat does y+=H*x;
  prop_grid_stat_transpose does x+=H'*y;
 
*/
void FUN_NAME (CONST_IN map_t *mapin, /**<[in] OPD defind on a square grid*/
	       const locstat_t *ostat, /**<[in] information about each clumn of a output loc grid*/
	       CONST_OUT double *phiout,    /**<[in,out] OPD defined on ostat*/
	       double alpha,       /**<[in] scaling of OPD*/
	       double displacex,   /**<[in] displacement of the ray */
	       double displacey,   /**<[in] displacement of the ray */
	       double scale,       /**<[in] scaling of the beam diameter (cone)*/
	       int wrap,           /**<[in] wrap input OPD or not*/
	       long colstart,      /**<[in] First column to do ray tracing*/
	       long colend         /**<[in] Last column (exclusive) to do ray tracing*/
	       ){
    CONST_OUT double *phiout2=0;
    double dplocx=0, dplocy=0;
    int nplocx=0, nplocy=0;
    int icol=0,irow=0;
    long offset=0;
    int collen=0;
    int missing=0;
    if(colend==0) colend = ostat->ncol;

    const int wrapx1 = mapin->nx;
    const int wrapy1 = mapin->ny;
    const int wrapx  = wrapx1 - 1;
    const int wrapy  = wrapy1 - 1;
    
    const double dx_in1 = 1./mapin->dx;
    const double dx_in2 = scale*dx_in1;
    const double dxout  = ostat->dx;
    const double dy_in1 = 1./mapin->dy;
    const double dy_in2 = scale*dy_in1;
    const double dyout  = ostat->dy;
    displacex = (displacex-mapin->ox)*dx_in1;
    displacey = (displacey-mapin->oy)*dy_in1;
    CONST_IN double *phiin  = mapin->p;
    const double xratio  = dxout*dx_in2;
    const double yratio  = dyout*dy_in2;
#if USE_OPTIM == 1
    if(fabs(xratio-1)<EPS && fabs(yratio-1)<EPS){
	/*loc_out and loc_in has the same grid sampling.*/
	for(icol=colstart; icol<colend; icol++)
#if TRANSPOSE == 0 && defined(__INTEL_COMPILER) && _OPENMP >= 200805
#pragma omp task
#endif
	    {
		/*grid size of loc_in and loc_out agree*/
		CONST_IN double *phicol, *phicol2;
		double bl,br,tl,tr;
		int rowdiv,rowdiv2;
		int irows;
		/*starting address of that col*/
		offset=ostat->cols[icol].pos;
		collen=ostat->cols[icol+1].pos-offset;/*exclusive*/

		phiout2=phiout+offset;
		dplocy=ostat->cols[icol].ystart*dy_in2+displacey;
		if(wrap){
		    dplocy=dplocy-floor(dplocy/(double)wrapy1)*wrapy1;
		}
		SPLIT(dplocy,dplocy,nplocy);
		dplocx=ostat->cols[icol].xstart*dx_in2+displacex;
		if(wrap){
		    dplocx-=wrapx1*floor(dplocx/(double)wrapx1);
		    irows=0;
		}else{
		    if(dplocx<0)
			irows=iceil(-dplocx);
		    else
			irows=0;
		}
		SPLIT(dplocx,dplocx,nplocx);

		if(wrap){
		    phicol  = phiin+nplocy*wrapx1+nplocx;
		    phicol2 = phiin+(nplocy==wrapy?0:nplocy+1)*wrapx1+nplocx;
		}else{
		    if(nplocy<wrapy && nplocy>=0){
			phicol=phiin+nplocy*wrapx1+nplocx;
			phicol2=phiin+(nplocy+1)*wrapx1+nplocx;
		    }else{
			if(nplocy!=wrapy){
			    missing+=collen;
			}
			goto end1;
		    }
		}
 
		bl=dplocx*dplocy;
		br=(1.-dplocx)*dplocy;
		tl=dplocx*(1-dplocy);
		tr=(1.-dplocx)*(1-dplocy);
		
		rowdiv=wrapx-nplocx;
		/*max number of rows possible before wraping*/
		if(rowdiv>collen) rowdiv=collen;
		if(rowdiv<0)rowdiv=0;
		/*reduce to number of needed*/

		for(irow=irows; irow<rowdiv; irow++){
#if TRANSPOSE == 0
		    phiout2[irow]+=alpha*
			(bl*phicol2[irow+1]
			 +br*phicol2[irow]
			 +tl*phicol[irow+1]
			 +tr*phicol[irow]);
#else
		    double tmp=alpha*phiout2[irow];
		    phicol2[irow+1]+=tmp*bl;
		    phicol2[irow]+=tmp*br;
		    phicol[irow+1]+=tmp*tl;
		    phicol[irow]+=tmp*tr;
#endif
		}

		if(wrap){
		    while(rowdiv < collen){/*can wrap several times*/ 
			irow=rowdiv;
#if TRANSPOSE == 0
			phiout2[irow]+=alpha*
			    (bl*phicol2[irow-wrapx]
			     +br*phicol2[irow]
			     +tl*phicol[irow-wrapx]
			     +tr*phicol[irow]);
#else
			double tmp=alpha*phiout2[irow];
			phicol2[irow-wrapx]+=tmp*bl;
			phicol2[irow]+=tmp*br;
			phicol[irow-wrapx]+=tmp*tl;
			phicol[irow]+=tmp*tr;
#endif
			rowdiv++;
			
			rowdiv2=rowdiv+wrapx;
			if(rowdiv2>collen) rowdiv2=collen;

			phicol-=wrapx1;/*wrap the second part*/
			phicol2-=wrapx1;

			for(irow=rowdiv; irow<rowdiv2; irow++){
#if TRANSPOSE == 0
			    phiout2[irow]+=alpha*
				(+bl*phicol2[irow+1]
				 +br*phicol2[irow]
				 +tl*phicol[irow+1]
				 +tr*phicol[irow]);
#else
			    double tmp2=alpha*phiout2[irow];
			    phicol2[irow+1]+=tmp2*bl;
			    phicol2[irow]+=tmp2*br;
			    phicol[irow+1]+=tmp2*tl;
			    phicol[irow]+=tmp2*tr;
#endif
			}
			rowdiv=rowdiv2;
		    }
		}
    end1:;
	    }/*end for icol*/
    }else{
	/*grid size of loc_in and loc_out doesn't agree*/
	for(icol=colstart; icol<colend; icol++)
#if TRANSPOSE == 0 && defined(__INTEL_COMPILER) && _OPENMP >= 200805
#pragma omp task
#endif
	    {
		CONST_IN double *phicol, *phicol2;
		double dplocx0;
		int nplocx0;
		int rowdiv,rowdiv2;
		int irows;
		/*starting address of that col*/
		offset=ostat->cols[icol].pos;
		collen=ostat->cols[icol+1].pos-offset;/*exclusive*/

		phiout2=phiout+offset;
		dplocy=ostat->cols[icol].ystart*dy_in2+displacey;
		if(wrap){
		    dplocy=dplocy-floor(dplocy/(double)wrapy1)*wrapy1;
		}
		SPLIT(dplocy,dplocy,nplocy);
		const double dplocy1=1.-dplocy;
		dplocx=ostat->cols[icol].xstart*dx_in2+displacex;
		if(wrap){
		    dplocx-=wrapx1*floor(dplocx/(double)wrapx1);
		    irows=0;
		}else{
		    if(dplocx<0){
			irows=iceil(-dplocx/xratio);
		    }else
			irows=0;
		}

		if(wrap){
		    phicol  = phiin+nplocy*wrapx1;
		    phicol2 = phiin+(nplocy==wrapy?0:(nplocy+1))*wrapx1;
		}else{
		    if(nplocy<wrapy && nplocy>=0){
			phicol=phiin+nplocy*wrapx1;
			phicol2=phiin+(nplocy+1)*wrapx1;
		    }else{
			if(nplocy!=wrapy){
			    missing+=collen;
			}
			goto end2;
		    }
		}
 	  
		/*last row to do if not wrap*/
		/*figure out which row we can go in this segment. */
		  
		rowdiv = iceil((wrapx-dplocx)/xratio);
		if(rowdiv<0) rowdiv=0;/*rowdiv may be -1 if dplocx=wrapx+0.* */
		if(rowdiv>collen) rowdiv=collen;
		dplocx = dplocx + xratio*irows;
		for(irow=irows; irow<rowdiv; irow++){/*no wrap*/
		    SPLIT(dplocx,dplocx0,nplocx0);
#if TRANSPOSE == 0
		    phiout2[irow]+=alpha*
			(+dplocy1*((1-dplocx0)*phicol[nplocx0]
				   +dplocx0*phicol[nplocx0+1])
			 +dplocy *((1-dplocx0)*phicol2[nplocx0]
				   +dplocx0*phicol2[nplocx0+1]));
#else
		    double tmp=phiout2[irow]*alpha;
		    phicol[nplocx0]+=tmp*dplocy1*(1-dplocx0);
		    phicol[nplocx0+1]+=tmp*dplocy1*dplocx0;
		    phicol2[nplocx0]+=tmp*dplocy*(1-dplocx0);
		    phicol2[nplocx0+1]+=tmp*dplocy*dplocx0;
#endif
		    dplocx+=xratio;
		}
		if(wrap){
		    while(rowdiv < collen){/*can wrap several times*/ 
			/*falls on the edge. */
			while(dplocx<wrapx1 && rowdiv<collen){
			    SPLIT(dplocx,dplocx0,nplocx0);
#if TRANSPOSE == 0
			    phiout2[rowdiv]+=alpha*
				(+dplocy1*((1-dplocx0)*phicol[nplocx0]
					   +dplocx0*phicol[nplocx0-wrapx])
				 +dplocy *((1-dplocx0)*phicol2[nplocx0]
					   +dplocx0*phicol2[nplocx0-wrapx]));
#else
			    double tmp=phiout2[rowdiv]*alpha;
			    phicol[nplocx0]+=tmp*dplocy1*(1-dplocx0);
			    phicol[nplocx0-wrapx]+=tmp*dplocy1*dplocx0;
			    phicol2[nplocx0]+=tmp*dplocy*(1-dplocx0);
			    phicol2[nplocx0-wrapx]+=tmp*dplocy*dplocx0;
#endif
			    dplocx+=xratio;
			    rowdiv++;
			}
			dplocx=dplocx-wrapx1;			
			rowdiv2=iceil((wrapx-dplocx)/xratio);
			if(rowdiv2>collen) rowdiv2=collen;
			for(irow=rowdiv; irow<rowdiv2; irow++){
			    SPLIT(dplocx,dplocx0,nplocx0);
#if TRANSPOSE == 0
			    phiout2[irow]+=alpha*
				(+dplocy1*((1-dplocx0)*phicol[nplocx0]
					   +dplocx0*phicol[nplocx0+1])
				 +dplocy *((1-dplocx0)*phicol2[nplocx0]
					   +dplocx0*phicol2[nplocx0+1]));
#else
			    double tmp=phiout2[irow]*alpha;
			    phicol[nplocx0]+=tmp*dplocy1*(1-dplocx0);
			    phicol[nplocx0+1]+=tmp*dplocy1*dplocx0;
			    phicol2[nplocx0]+=tmp*dplocy*(1-dplocx0);
			    phicol2[nplocx0+1]+=tmp*dplocy*dplocx0;
#endif
			    dplocx+=xratio;
			}
			rowdiv=rowdiv2;
		    }
		}
	    end2:;
	    }/*end for icol*/
    }/*fabs(dx_in2*dxout-1)<EPS*/
#else
    /*non optimized case. slower, but hopefully accurate*/
    for(icol=colstart; icol<colend; icol++)
#if TRANSPOSE == 0 && defined(__INTEL_COMPILER) && _OPENMP >= 200805
#pragma omp task
#endif
	{
	    double dplocx0;
	    int nplocy1, nplocx1;

	    offset=ostat->cols[icol].pos;
	    collen=ostat->cols[icol+1].pos-offset;
	    phiout2=phiout+offset;
	    dplocy=ostat->cols[icol].ystart*dy_in2+displacey;
	    SPLIT(dplocy,dplocy,nplocy);
	    const double dplocy1=1.-dplocy;
	    if(nplocy<0||nplocy>=wrapy){
		if(wrap){
		    while(nplocy<0)
			nplocy+=wrapy1;
		    while(nplocy>wrapy)
			nplocy-=wrapy1;
		    nplocy1=(nplocy==wrapx?0:nplocy+1);
		}else{
		    if(nplocy!=wrapy){
			missing+=collen;
		    }
		    continue;
		}
	    }else{
		nplocy1=nplocy+1;
	    }
	    dplocx0=(ostat->cols[icol].xstart)*dx_in2+displacex;
	    for(irow=0; irow<collen; irow++){
		SPLIT(dplocx0,dplocx,nplocx);
		dplocx0+=xratio;
		if(nplocx>=0 && nplocx<wrapx){
		    nplocx1=nplocx+1;
		}else{
		    if(wrap){
			while(nplocx<0)
			    nplocx+=wrapx1;
			while(nplocx>wrapx)
			    nplocx-=wrapx1;
			nplocx1=(nplocx==wrapx?0:nplocx+1);
		    }else{
			if(nplocx!=wrapx){
			    missing++;
			}
			continue;
		    }
		}
#if TRANSPOSE == 0
		phiout2[irow]+= alpha*
		    (dplocx * (dplocy * phiin[(nplocx1) + (nplocy1)*wrapx1]
			       +dplocy1 * phiin[(nplocx1) + nplocy*wrapx1])
		     + (1-dplocx) * (dplocy * phiin[nplocx + (nplocy1)*wrapx1]
				     +dplocy1 * phiin[nplocx + nplocy*wrapx1]));
#else
		double tmp=alpha*phiout2[irow];
		phiin[(nplocx1) + (nplocy1)*wrapx1]+=tmp*dplocx*dplocy;
		phiin[(nplocx1) + nplocy*wrapx1]+=tmp*dplocx*dplocy1;
		phiin[nplocx + (nplocy1)*wrapx1]+=tmp*(1-dplocx)*dplocy;
		phiin[nplocx + nplocy*wrapx1]+=tmp*(1-dplocx)*dplocy1;
#endif		
	    }/*for irow*/
	}/*for icol*/
#endif
    if(missing>0){
	warning(" %d points not covered by input screen\n", missing);
	print_backtrace();
    }
#if TRANSPOSE == 0 && defined(__INTEL_COMPILER) && _OPENMP >= 200805
#pragma omp taskwait
#endif
}/*function*/

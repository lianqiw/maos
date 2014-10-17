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
  Included by accphi.c for prop_grid_map and prop_grid_map_tranpose
  derived from prop_grid_stat.c on 2011-10-05.
*/
#undef CONST_IN
#undef CONST_OUT
#undef FUN_NAME
#if TRANSPOSE == 0
#define CONST_IN const
#define CONST_OUT 
#define FUN_NAME prop_grid_map
#else
#define CONST_IN 
#define CONST_OUT const
#define FUN_NAME prop_grid_map_transpose
#endif

/**
   Propagate OPD defines on grid mapin to grid mapout.  alpha is the scaling of
   data. displacex, displacy is the displacement of the center of the beam on
   the input grid. scale is the cone effect.*/

/*
  Implement prop_grid_stat or prop_grid_stat_tranpose in a unified way. 
  Let mapin be x, phiout be y. H be propagator.
  prop_grid_stat does y+=H*x;
  prop_grid_stat_transpose does x+=H'*y;
 
*/
void FUN_NAME (CONST_IN map_t *mapin,   /**<[in] OPD defind on a square grid*/
	       CONST_OUT map_t *mapout, /**<[in,out] output OPD defined on a square grid*/
	       ARG_PROP,
	       int wrap,           /**<[in] wrap input OPD or not*/
	       long colstart,      /**<[in] First column to do ray tracing*/
	       long colend         /**<[in] Last column (exclusive) to do ray tracing*/
    ){
    CONST_OUT double *phiout=mapout->p;
    /*With OpenMP compiler complained uninitialized value for the following
      because they are treated as firstprivate*/
    if(colend==0) colend = mapout->ny;
    const int wrapx1 = mapin->nx;
    const int wrapy1 = mapin->ny;
    const int wrapx  = wrapx1 - 1;
    const int wrapy  = wrapy1 - 1;
    
    const double dx_in1 = 1./mapin->dx;
    const double dx_in2 = scale*dx_in1;
    const double dy_in1 = 1./mapin->dy;
    const double dy_in2 = scale*dy_in1;
    const double dxout  = mapout->dx;
    const double dyout  = mapout->dy;
    displacex = (displacex-mapin->ox)*dx_in1;
    displacey = (displacey-mapin->oy)*dy_in1;
    CONST_IN double *phiin  = mapin->p;
    const double xratio  = dxout*dx_in2;
    const double yratio  = dyout*dy_in2;
    DEF_ENV_FLAG(PROP_GRID_MAP_OPTIM, 1);
    if(PROP_GRID_MAP_OPTIM){
	if(fabs(xratio-1.)<EPS && fabs(yratio-1.)<EPS){
	    int irows;
	    double dplocx=mapout->ox*dx_in2+displacex;
	    int nplocx;
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
	    /*loc_out and loc_in has the same grid sampling.*/
	    OMPTASK_FOR(icol, colstart, colend){
		int rowdiv,rowdiv2, nplocy;
		/*grid size of loc_in and loc_out agree*/
		double bl, br, tl, tr;
		CONST_IN double *phicol, *phicol2;
		/*starting address of that col*/
		int offset=icol*mapout->nx;
		int collen=mapout->nx;
		CONST_OUT double *phiout2=phiout+offset;
		double dplocy=(mapout->oy+icol*mapout->dy)*dy_in2+displacey;
		if(wrap){
		    dplocy=dplocy-floor(dplocy/(double)wrapy1)*wrapy1;
		}else if(dplocy<0 || dplocy>wrapy){
		    goto end1;
		}
		/*The points right on dplocy is handled correctly*/
		SPLIT(dplocy,dplocy,nplocy);
		phicol  = phiin+nplocy*wrapx1+nplocx;
		phicol2 = phiin+(nplocy==wrapy?0:(nplocy+1))*wrapx1+nplocx;
	    
		bl=dplocx*dplocy;
		br=(1.-dplocx)*dplocy;
		tl=dplocx*(1-dplocy);
		tr=(1.-dplocx)*(1-dplocy);
	    
		rowdiv=wrapx-nplocx;
		/*max number of rows possible before wraping*/
		if(rowdiv>collen) rowdiv=collen;
		if(rowdiv<0)rowdiv=0;
		/*reduce to number of needed*/
		for(int irow=irows; irow<rowdiv; irow++){
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
			int irow=rowdiv;
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
	    OMPTASK_END;
	}else{
	    OMPTASK_FOR(icol, colstart, colend){
		/*grid size of loc_in and loc_out doesn't agree*/
		CONST_IN double *phicol, *phicol2;
		double dplocx0;
		int nplocx0,nplocy;
		int rowdiv,rowdiv2;
		int irows;

		/*starting address of that col*/
		int offset=icol*mapout->nx;
		int collen=mapout->nx;
		CONST_OUT double *phiout2=phiout+offset;
		double dplocx=mapout->ox*dx_in2+displacex;
		if(wrap){
		    dplocx-=wrapx1*floor(dplocx/(double)wrapx1);
		    irows=0;
		}else{
		    if(dplocx<0)
			irows=iceil(-dplocx);
		    else
			irows=0;
		}
		double dplocy=(mapout->oy+icol*mapout->dy)*dy_in2+displacey;
		if(wrap){
		    dplocy-=floor(dplocy/(double)wrapy1)*wrapy1;
		}else if(dplocy<0 || dplocy>wrapy){
		    goto end2;
		}
		SPLIT(dplocy,dplocy,nplocy);
		phicol  = phiin+nplocy*wrapx1;
		phicol2 = phiin+(nplocy==wrapy?0:(nplocy+1))*wrapx1;
		const double dplocy1=1.-dplocy;
 	  
		/*last row to do if not wrap*/
		/*figure out which row we can go in this segment. */
		  
		rowdiv = iceil((wrapx-dplocx)/xratio);
		if(rowdiv<0) rowdiv=0;/*rowdiv may be -1 if dplocx=wrapx+0.* */
		if(rowdiv>collen) rowdiv=collen;
		dplocx += xratio*irows;
		for(int irow=irows; irow<rowdiv; irow++){/*no wrap*/
		    SPLIT(dplocx,dplocx0,nplocx0);
		    dplocx+=xratio;
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
		}
		if(wrap){
		    while(rowdiv < collen){/*can wrap several times*/ 
			/*falls on the edge. */
			while(dplocx<wrapx1 && rowdiv<collen){
			    SPLIT(dplocx,dplocx0,nplocx0);
			    dplocx+=xratio;
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
			    rowdiv++;
			}
			dplocx=dplocx-wrapx1;			
			rowdiv2=iceil((wrapx-dplocx)/xratio);
			if(rowdiv2>collen) rowdiv2=collen;
			for(int irow=rowdiv; irow<rowdiv2; irow++){
			    SPLIT(dplocx,dplocx0,nplocx0);
			    dplocx+=xratio;
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
			}
			rowdiv=rowdiv2;
		    }
		}/*if wrap*/
	      end2:;
	    }/*for icol*/
	    OMPTASK_END;
	}/*fabs(dx_in2*dxout-1)<EPS*/
    }else{
	/*non optimized case. slower, but hopefully accurate*/
	OMPTASK_FOR(icol, colstart, colend){
	    double dplocx,dplocx0;
	    int nplocx,nplocy,nplocy1, nplocx1;
	    int offset=icol*mapout->nx;
	    int collen=mapout->nx;
	    CONST_OUT double *phiout2=phiout+offset;
	    double dplocy=(mapout->oy+icol*mapout->dy)*dy_in2+displacey;
	    info2("dplocy=%g, map->dx=%g, map->dy=%g\n", dplocy, mapout->dx, mapout->dy);
	    if(wrap){
		dplocy=dplocy-floor(dplocy/(double)wrapy1)*wrapy1;
	    }else if (dplocy<0 || dplocy>wrapy){
		goto skip;
	    }
	    SPLIT(dplocy,dplocy,nplocy);
	    const double dplocy1=1.-dplocy;
	    nplocy1=(nplocy==wrapx?0:nplocy+1);
	    dplocx0=(mapout->ox)*dx_in2+displacex-xratio;
	    for(int irow=0; irow<collen; irow++){
		dplocx0+=xratio;
		info2("dplocx=%g\n", dplocx0);
		if(wrap){
		    dplocx0=dplocx0-floor(dplocx0/(double)wrapx1)*wrapx1;
		}else if(dplocx0<0 || dplocx0>wrapx){
		    continue;
		}
		SPLIT(dplocx0,dplocx,nplocx);
		nplocx1=(nplocx==wrapx?0:nplocx+1);
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
	  skip:;
	}/*for icol*/
	OMPTASK_END;
    }
}/*function*/

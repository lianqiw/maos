/**
   Propagate OPD defines on grid mapin to subapertures.  alpha is the scaling of
   data. displacex, displacy is the displacement of the center of the beam on
   the input grid.  scale is the cone effect.*/
void prop_grid_pts(ARGIN_GRID,
		   ARGOUT_PTS,
		   ARG_PROP,
		   int wrap,           /**<[in] wrap input OPD or not*/
		   long sastart,       /**<[in] The starting subaperture to trace ray*/
		   long saend          /**<[in] The last (exclusive) subaperture to trace ray*/
		   ){
    /*
       2010-01-02: Improved by saving x interpolations for each subaperture
       during non-matched case. original function renamed to prop_grid_pts_old
       time cut by almost 3 fold for each LGS WFS.  0.40s->0.16s
     */
    (void)ampout;
    int isa, ipix,jpix;
    int mx,my;
    int nplocx2;
    const int ninx   = mapin->nx;
    const int niny   = mapin->ny;
    const int nx     = pts->nx;
    const double dx_in1 = 1/mapin->dx;
    const double dxout  = pts->dx;
    const double dx_in2 = scale*dx_in1;
    displacex = (displacex-mapin->ox)*dx_in1;
    displacey = (displacey-mapin->oy)*dx_in1;
    if(!saend) saend=pts->nsa;
    if(USE_OPTIM && fabs(dx_in2*dxout-1)<EPS){
	double dplocx, dplocy;
	int nplocx, nplocy;
        for(isa=sastart; isa<saend; isa++){
	    double (*phioutsq)[nx]=(double(*)[nx])(phiout+nx*nx*isa);
	    double (*phiin)[ninx]=(double(*)[ninx])(mapin->p);
	    const double origx=pts->origx[isa];
	    const double origy=pts->origy[isa];
	    double w11,w10,w01,w00;
	    
	    dplocy = myfma(origy,dx_in2,displacey);
	    dplocx = myfma(origx,dx_in2,displacex);
	    SPLIT(dplocx,dplocx,nplocx);
	    SPLIT(dplocy,dplocy,nplocy);
	    if(!wrap){
		int sx, sy;
		if(nplocy<0){
		    sy=-nplocy;
		}else
		    sy=0;
		if(nplocx<0){
		    sx=-nplocx;
		}else
		    sx=0;
		
		my=niny-nplocy-1;//remaining possible points.
		if(my>nx){
		    my=nx;
		}
		if(my<=sy){
		    continue;
		}
		mx=ninx-nplocx-1;
		if(mx>nx){
		    mx=nx;
		}
		if(mx<=sx){
		    continue;
		}

		if((dplocx)<EPS && (dplocy)<EPS){
		    /*aligned perfectly.*/
		    for(jpix=sy; jpix<my; jpix++){
			//nplocy2=nplocy+jpix;
			double *restrict phiin2=phiin[nplocy+jpix];
			double *restrict phiout2=phioutsq[jpix];
			for(ipix=sx; ipix<mx; ipix++){
			    phiout2[ipix]+=alpha*phiin2[nplocx+ipix];
			}
		    }
		}else{	
		    w11=dplocx*dplocy;
		    w10=(1.-dplocx)*dplocy;
		    w01=dplocx*(1.-dplocy);
		    w00=(1.-dplocx)*(1.-dplocy);	
		    for(jpix=sy; jpix<my; jpix++){
			//nplocy2=nplocy+jpix;
			double *restrict phiin2=phiin[nplocy+jpix];
			double *restrict phiin3=phiin[nplocy+jpix+1];
			double *restrict phiout2=phioutsq[jpix];
			for(ipix=sx; ipix<mx; ipix++){
			    phiout2[ipix]+=
				alpha*(+phiin2[nplocx+ipix]*w00
				       +phiin2[nplocx+ipix+1]*w01
				       +phiin3[nplocx+ipix]*w10
				       +phiin3[nplocx+ipix+1]*w11);
			}
		    }
		}
	    }else{//wraping
		double *phiin_1, *phiin_2;

		if(ninx < nx || niny < nx){
		    error("Input map is too small. wraps more than once\n");
		}
		w11=dplocx*dplocy;
		w10=(1.-dplocx)*dplocy;
		w01=dplocx*(1.-dplocy);
		w00=(1.-dplocx)*(1.-dplocy);

		nplocy-=niny*(nplocy/niny);
		if(nplocy<0)
		    nplocy+=niny;
		nplocx-=ninx*(nplocx/ninx);
		if(nplocx<0)
		    nplocx+=ninx;
		my=niny-nplocy-1;//remaining possible points.
		mx=ninx-nplocx-1;
		if(my>nx) my=nx;
		if(mx>nx) mx=nx;

		for(jpix=0; jpix<nx; jpix++){
		    if(jpix<my){
			phiin_1=phiin[nplocy];
			phiin_2=phiin[nplocy+1];
		    }else if(jpix==my){
			phiin_1=phiin[nplocy];
			phiin_2=phiin[nplocy+1-niny];
		    }else{
			phiin_1=phiin[nplocy-niny];
			phiin_2=phiin[nplocy+1-niny];
		    }
		    nplocx2=nplocx;
		    for(ipix=0; ipix<mx; ipix++){
			phioutsq[jpix][ipix]+=alpha*
			    (phiin_1[nplocx2]*w00
			     +phiin_1[nplocx2+1]*w01
			     +phiin_2[nplocx2]*w10
			     +phiin_2[nplocx2+1]*w11);
			nplocx2++;
		    }
		    if(mx<nx){
			ipix=mx;
			phioutsq[jpix][ipix]+=alpha*
			    (phiin_1[nplocx2]*w00
			     +phiin_1[nplocx2+1-ninx]*w01
			     +phiin_2[nplocx2]*w10
			     +phiin_2[nplocx2+1-ninx]*w11);
			nplocx2++;

			nplocx2-=ninx;
			for(ipix=mx+1; ipix<nx; ipix++){
			    phioutsq[jpix][ipix]+=alpha*
				(phiin_1[nplocx2]*w00
				 +phiin_1[nplocx2+1]*w01
				 +phiin_2[nplocx2]*w10
				 +phiin_2[nplocx2+1]*w11);
			    nplocx2++;
			}
		    }
		    nplocy++;
		}
	    }	    
	}
    }else{ /*different spacing. or non optim*/
	double dplocx, dplocy;
	double ratio;
	const double (*phiin)[ninx]=(const double(*)[ninx])(mapin->p);
	int nplocx, nplocy;
	ratio = dxout*dx_in2;
	for(isa=sastart; isa<saend; isa++){
	    double (*phioutsq)[nx]=(double(*)[nx])(phiout+nx*nx*isa);
	    const double origx=pts->origx[isa];
	    const double origy=pts->origy[isa];
	    double dplocx0, dplocy0;
	    dplocy0 = myfma(origy,dx_in2,displacey);
	    dplocx0 = myfma(origx,dx_in2,displacex);
	    if(!wrap){
		int sx, sy;
		if(dplocx0<0){
		    sx=iceil(-dplocx0/ratio);
		}else
		    sx=0;
		if(dplocy0<0){
		    sy=iceil(-dplocy0/ratio);
		}else
		    sy=0;

		my=iceil((niny-1-dplocy0)/ratio);
		if(my>nx){
		    my=nx;
		}
		mx=iceil((ninx-1-dplocx0)/ratio);
		if(mx>nx){
		    mx=nx;
		}
		
		int nplocxs[mx];
		double dplocxs[mx];

		dplocy0  = myfma(origy+(double)sy*dxout,dx_in2,displacey);
		dplocx0  = myfma(origx+(double)sx*dxout,dx_in2,displacex);
		
		for(ipix=sx; ipix<mx; ipix++){
		    SPLIT(dplocx0,dplocx,nplocx);
		    nplocxs[ipix]=nplocx;
		    dplocxs[ipix]=dplocx;
		    dplocx0+=ratio;
		}

		for(jpix = sy; jpix < my; jpix++){
		    SPLIT(dplocy0,dplocy,nplocy);
		    const double dplocy1=1.-dplocy;
		    const double *phiin2=phiin[nplocy];
		    const double *phiin3=phiin[nplocy+1];
		    double *phiout2=phioutsq[jpix];
		    for(ipix=sx; ipix<mx; ipix++){
			nplocx=nplocxs[ipix];
			dplocx=dplocxs[ipix];
			phiout2[ipix]+=alpha*
			    (+(phiin2[nplocx]
			       +(phiin2[nplocx+1]-phiin2[nplocx])*dplocx)
			     *dplocy1
			     +(phiin3[nplocx]
			       +(phiin3[nplocx+1]-phiin3[nplocx])*dplocx)
			     *dplocy);
		    }
		    dplocy0+=ratio;
		}
	    }else{
		const double *phiin_1, *phiin_2;
		double dplocy1;
		int mx0;
		if(ninx < nx*ratio || niny < nx*ratio){
		    info("nx=%d, ratio=%g, ninx=%d\n",nx,ratio,ninx);
		    error("Input map is too small. wraps more than once\n");
		}
		dplocy0-=niny*ifloor(dplocy0/niny);
		dplocx0-=ninx*ifloor(dplocx0/ninx);
		mx0=iceil((ninx-1.-dplocx0)/ratio);
		if(mx0>nx) mx0=nx;
		if(mx0<0) mx0=0;
		
		int nplocxs[nx];
		int nplocxs2[nx];
		double dplocxs[nx];
		mx=mx0;
		for(ipix=0; ipix<mx; ipix++){
		    SPLIT(dplocx0,dplocx,nplocx);
		    nplocxs[ipix]=nplocx;
		    nplocxs2[ipix]=nplocx+1;
		    dplocxs[ipix]=dplocx;
		    dplocx0+=ratio;
		}
		if(mx<nx){
		    while(dplocx0<ninx && mx<nx){//falls on the edge
			SPLIT(dplocx0,dplocx,nplocx);
			nplocxs[mx]=nplocx;
			nplocxs2[mx]=0;
			dplocxs[mx]=dplocx;
			dplocx0+=ratio;
			mx++;
		    }
		    dplocx0-=ninx;
		    for(ipix=mx; ipix<nx; ipix++){
			SPLIT(dplocx0,dplocx,nplocx);
			nplocxs[ipix]=nplocx;
			nplocxs2[ipix]=nplocx+1;
			dplocxs[ipix]=dplocx;
			dplocx0+=ratio;
		    }
		}
		for(jpix=0; jpix<nx; jpix++){
		    mx=mx0;
		    SPLIT(dplocy0,dplocy,nplocy);
		    dplocy1=1.-dplocy;
		    if(nplocy<niny-1){
			phiin_1=phiin[nplocy];
			phiin_2=phiin[nplocy+1];
		    }else if(nplocy==niny-1){
			phiin_1=phiin[nplocy];
			phiin_2=phiin[0];
		    }else{
			phiin_1=phiin[nplocy-niny];
			phiin_2=phiin[nplocy+1-niny];
		    }
		    double *phiout2=phioutsq[jpix];
		    for(ipix=0; ipix<nx; ipix++){
			nplocx=nplocxs[ipix];
			nplocx2=nplocxs2[ipix];
			dplocx=dplocxs[ipix];
			phiout2[ipix]+=alpha*
			    (+(phiin_1[nplocx]
			       +(phiin_1[nplocx2]-phiin_1[nplocx])*dplocx)
			     *dplocy1
			     +(phiin_2[nplocx]
			       +(phiin_2[nplocx2]-phiin_2[nplocx])*dplocx)
			     *dplocy);
			
		    }
		    dplocy0+=ratio;
		}
	    }/*wrap*/
	}/*end isa*/
    }
}

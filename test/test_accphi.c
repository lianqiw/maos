#include <math.h>
#include <time.h>
#include <string.h>
#include "../lib/aos.h"
TIC;
//const double EPS=1e-12;
/*for propagation to WFS PTS*/
typedef struct PROP_GRD_PTS2_T{
    /*input*/
    map_t *mapin;
    /*output*/
    pts_t *pts;           /*out*/
    double *phiout;       /*stores out opd. needs to be preallocated.*/
    double alpha;
    double displace[2];
    double scale;
    int wrap;
    int optim;
    pthread_mutex_t mutex;
    int isa;
    int saend;
}PROP_GRD_PTS2_T;
static void prop_grid_pts_wrap2(PROP_GRD_PTS2_T *data){
    /**
       2010-01-02: Improved by saving x interpolations for each
       subaperture during non-matched case. original function
       renamed to prop_grid_pts_old time cut by almost 3 fold for
       each LGS WFS.  0.40s->0.16s
     */
    map_t *mapin=data->mapin;
    pts_t *pts=data->pts;
    double *phiout0=data->phiout;
    double alpha=data->alpha;
    double displacex=data->displace[0];
    double displacey=data->displace[1];
    double scale=data->scale;
    int wrap=data->wrap;
    int saend=data->saend;
    int optim=1;

    int isa, ipix,jpix;
    int mx,my;
    int nplocx2;
    //const int optim = 0;//for debugging
    const int ninx   = mapin->nx;
    const int niny   = mapin->ny;
    const int nx     = pts->nx;
    const double dx_in1 = 1/mapin->dx;
    const double dxout  = pts->dx;
    const double dx_in2 = scale*dx_in1;
    displacex = (displacex-mapin->ox)*dx_in1;
    displacey = (displacey-mapin->oy)*dx_in1;
    if(saend==0)
	saend=pts->nsa;
    if(optim && fabs(dx_in2*dxout-1)<EPS){
	double dplocx, dplocy;
	int nplocx, nplocy;
	while(1){
	    LOCK(data->mutex);
	    isa=data->isa++;
	    UNLOCK(data->mutex);
	    if(isa>=saend){
		break;
	    }else{
		double (*phiout)[nx]=(double(*)[nx])(phiout0+nx*nx*isa);
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
			    double *restrict phiout2=phiout[jpix];
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
			    double *restrict phiout2=phiout[jpix];
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
			    phiout[jpix][ipix]+=alpha*
				(phiin_1[nplocx2]*w00
				 +phiin_1[nplocx2+1]*w01
				 +phiin_2[nplocx2]*w10
				 +phiin_2[nplocx2+1]*w11);
			    nplocx2++;
			}
			if(mx<nx){
			    ipix=mx;
			    phiout[jpix][ipix]+=alpha*
				(phiin_1[nplocx2]*w00
				 +phiin_1[nplocx2+1-ninx]*w01
				 +phiin_2[nplocx2]*w10
				 +phiin_2[nplocx2+1-ninx]*w11);
			    nplocx2++;

			    nplocx2-=ninx;
			    for(ipix=mx+1; ipix<nx; ipix++){
				phiout[jpix][ipix]+=alpha*
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
	}
    }else{ /*different spacing. or non optim*/
	double dplocx, dplocy;
	double ratio;
	const double (*phiin)[ninx]=(const double(*)[ninx])(mapin->p);
	int nplocx, nplocy;
	ratio = dxout*dx_in2;
	
	while(LOCK(data->mutex), isa=data->isa++, UNLOCK(data->mutex), isa<saend){
	    double (*phiout)[nx]=(double(*)[nx])(phiout0+nx*nx*isa);
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
		    double *phiout2=phiout[jpix];
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
		    double *phiout2=phiout[jpix];
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
static void test_accuracy(void){
    double r0=0.15;
    double L0=30;
    rand_t rstat;
    //long m=2048, n=2048;
    double dx=1./2.;
    long m=(long)round(10/dx/16)*16;
    long n=(long)round(20/dx/16)*16;

    seed_rand(&rstat, 2);
    int nlayer=1;
    double wt[1]={1};
    map_t **screens;
    GENSCREEN_T data;
    memset(&data, 0, sizeof(GENSCREEN_T));
    data.rstat=&rstat;
    data.nx=m;
    data.ny=n;
    data.dx=dx;
    data.r0=r0;
    data.l0=L0;
    data.wt=wt;
    data.nlayer=nlayer;
    screens=vonkarman_screen(&data);
    //screens=vonkarman_screen(&rstat,m,n,dx,r0,L0,wt,nlayer,0,1);
    map_t *screen=screens[0];
    mapwrite(screen,"accphi_screen");
    //loc for the map
    loc_t *locin=mksqloc(m,n,dx,screen->ox,screen->oy);

    dmat *A=dnew(60,60);
    dcircle(A,30,30,30,1);
    dwrite(A,"accphi_A");
    pts_t *pts=calloc(1,sizeof(pts_t));
    map_t *maptmp=mapnew(A->nx, A->ny, 0.5, A->p);
    loc_t *tmp=map2loc(maptmp);
    memcpy(pts,tmp,sizeof(loc_t));
    free(tmp);
    pts->dx=1./64.;
    pts->nx=32;
    loc_t *loc=pts2loc(pts);
    
    loc_create_stat(loc);
    locstat_t *locstat=loc->stat;
    locwrite((loc_t*)pts,"accphi_pts");
    locwrite(loc,"accphi_loc");
    double *phi_pts, *phi_loc, *phi_stat, *phi_pts1, *phi_stat1;
    double *phi_loc2loc, *phi_h, *phi_cub, *phi_cubh;

    double displacex,displacey,scale;
    displacex=0;
    displacey=0;
    scale=1.2;//.414065;
    int ic;
    int wrap=0;
    int nc=1;
    double cubic=0.3;
    int ii;
	    

    for(ic=0; ic<nc; ic++){
	if(nc>1){
	    displacex=(randu(&rstat)-0.5)*80;
	    displacey=(randu(&rstat)-0.5)*80;
	    scale=randu(&rstat)*3+0.1;
	    if(randu(&rstat)>0.9)
		scale=1;
	}
	for(wrap=0; wrap<1; wrap++){
	    info("displacex=%g, displacey=%g, scale=%g, wrap=%d\n",
		 displacex, displacey,scale,wrap);
	    
	    double diff1, diff2,diff3,diff11,diff31,diff131,
		diff231,diff14,diff15,diff6;
	    diff1=0;
	    diff2=0;
	    diff3=0;
	    diff11=0;
	    diff31=0;
	    diff131=0;
	    diff231=0;
	    diff14=0;
	    diff15=0;
	    diff6=0;

	    phi_pts=calloc(loc->nloc, sizeof(double));
	    phi_loc=calloc(loc->nloc, sizeof(double));
	    phi_stat=calloc(loc->nloc, sizeof(double));
	    phi_pts1=calloc(loc->nloc, sizeof(double));
	    phi_stat1=calloc(loc->nloc, sizeof(double));
	    phi_loc2loc=calloc(loc->nloc, sizeof(double));

	    memset(phi_pts, 0, sizeof(double)*loc->nloc);
	    memset(phi_loc, 0, sizeof(double)*loc->nloc);
	    memset(phi_stat, 0, sizeof(double)*loc->nloc);
	    memset(phi_pts1,0, sizeof(double)*loc->nloc);
	    memset(phi_stat1,0, sizeof(double)*loc->nloc);
	    memset(phi_loc2loc, 0, sizeof(double)*loc->nloc);
	   

	    tic;
	    prop_grid_pts(screen, pts, phi_pts, -1, displacex, displacey,
			  scale, wrap, 0,0);
	    toc("pts optim\t\t");
	    tic;
	    prop_grid_pts(screen, pts, phi_pts1,-1, displacex, displacey,
			  scale, wrap, 0,0);
	    toc("pts nonoptim\t");
	    tic;
	    prop_grid(screen, loc, phi_loc, -1,displacex, displacey, scale, 
		      wrap, 0,0);
	    toc("loc\t\t\t");
	    tic;
	    prop_grid_stat(screen, locstat, phi_stat, -1,displacex, displacey, 
			   scale, wrap, 0,0);
	    toc("locstat optim\t");
	    tic;
	    prop_grid_stat(screen, locstat, phi_stat1,-1, displacex, displacey, 
			   scale, wrap, 0,0);
	    toc("locstat nonoptim\t");
	    tic;
	    prop_nongrid(locin, screen->p,loc, NULL,phi_loc2loc, -1,displacex, displacey, 
			 scale,0,0);
	    toc("nongrid\t\t");


	    phi_h=calloc(loc->nloc,sizeof(double));
	    phi_cub=calloc(loc->nloc,sizeof(double));
	    phi_cubh=calloc(loc->nloc, sizeof(double));
	    tic;
	    prop_nongrid_cubic(locin,screen->p,loc,NULL,phi_cub,-1,
			       displacex, displacey, 
			       scale, cubic,0,0);
	    toc("nongrid, cubic\t");
	    tic;
	    dsp *hfor=mkh(locin, loc, NULL,displacex, displacey, scale,0,0);
	    toc("mkh\t\t\t");
	    tic;
	    spmulvec(phi_h,hfor, screen->p, -1);
	    toc("mul h\t\t");
	    spfree(hfor);
	    
	    dsp *hforcubic;
	    tic;
	    hforcubic=mkh(locin, loc, NULL, displacex, displacey, scale, 1, cubic);
	    toc("mkh cubic \t\t");
	    tic;
	    spmulvec(phi_cubh, hforcubic,screen->p,-1);
	    toc("cubic mul h\t\t");
	    double diffc12=0,diff45=0;
	    for(ii=0; ii<loc->nloc; ii++){
		diff1+=fabs(phi_loc[ii]-phi_pts[ii]);
		diff2+=fabs(phi_stat[ii]-phi_loc[ii]);
		diff3+=fabs(phi_stat[ii]-phi_pts[ii]);
		diff11+=fabs(phi_pts[ii]-phi_pts1[ii]);
		diff31+=fabs(phi_stat[ii]-phi_stat1[ii]);
		diff131+=fabs(phi_pts[ii]-phi_stat1[ii]);
		diff231+=fabs(phi_loc[ii]-phi_stat1[ii]);
		diff14+=fabs(phi_loc2loc[ii]-phi_pts[ii]);
		diff15+=fabs(phi_h[ii]-phi_pts[ii]);
		diff45+=fabs(phi_loc2loc[ii]-phi_h[ii]);
		diffc12+=fabs(phi_cub[ii]-phi_cubh[ii]);
	    }
	    info2("(pts-loc)=\t%g\n(loc-stat)=\t%g\n(stat-pts)=\t%g\n"
		  "(pts-pts1)=\t%g\n(stat-stat1)=\t%g\n"
		  "(loc2loc-pts)=\t%g\n(h-pts)=\t%g\n"
		  "(loc2loc-h)=\t%g\n(cub-cub_h)=\t%g\n",
		  diff1, diff2,diff3,diff11,diff31,diff14,diff15,diff45,diffc12);

	    if(diff1>1.e-12||diff2>1.e-12||diff3>1.e-12||
	       diff11>1.e-12||diff31>1.e-12||diff14>1.e-12||diff15>1.e-12||diff6>1.e-12
	       ||diffc12>1.e-12 || diff45>0){
		writedbl(phi_pts1,loc->nloc,1,"accphi_pts1.bin");
		writedbl(phi_loc,loc->nloc,1,"accphi_loc.bin");
		writedbl(phi_stat,loc->nloc,1,"accphi_stat.bin");
		writedbl(phi_stat1,loc->nloc,1,"accphi_stat1.bin");
		writedbl(phi_loc2loc,loc->nloc,1,"accphi_loc2loc.bin");
		writedbl(phi_h,loc->nloc,1,"accphi_h.bin");
		writedbl(phi_cub,loc->nloc,1,"accphi_cub.bin");
		writedbl(phi_cubh,loc->nloc,1,"accphi_cubh.bin");
		info("saved\n");
		if(nc>1) getchar();
	    }
	    writedbl(phi_pts,loc->nloc,1,"accphi_pts.bin");
	    writedbl(phi_cub,loc->nloc,1,"accphi_cub.bin");
		
	    free(phi_loc); free(phi_stat);
	    free(phi_pts1);free(phi_stat1);
	    free(phi_loc2loc);

	    free(phi_pts);
	    free(phi_h); 
	    free(phi_cub); 
	    free(phi_cubh);
	}
    }
    free(screen->p);
    free(screens);
    locfree(locin);
}
/*
  matched grip, non wrap, test ok.
*/

static void test_speed(int nthread){
   double r0=0.15;
    double L0=30;
    rand_t rstat;
    //long m=2048, n=2048;
    double dx=1./2.;
    long m=(long)round(10/dx/16)*16;
    long n=(long)round(20/dx/16)*16;

    seed_rand(&rstat, 2);
    int nlayer=1;
    double wt[1]={1};
    map_t **screens;
    GENSCREEN_T data;
    memset(&data, 0,sizeof(GENSCREEN_T));
    data.rstat=&rstat;
    data.nx=m;
    data.ny=n;
    data.dx=dx;
    data.r0=r0;
    data.l0=L0;
    data.wt=wt;
    data.nlayer=nlayer;
    screens=vonkarman_screen(&data);
    //screens=vonkarman_screen(&rstat,m,n,dx,r0,L0,wt,nlayer,0,1);
    info("\n");
    map_t *screen=screens[0];
    double dsa=0.1;
    dmat *A=dnew(30/dsa,30/dsa);
    dcircle(A,30,30,15/dsa,1);
    pts_t *pts=calloc(1,sizeof(pts_t));
    map_t *maptmp=mapnew(A->nx, A->ny, dsa, A->p);
    loc_t *tmp=map2loc(maptmp);
    memcpy(pts,tmp,sizeof(loc_t));
    free(tmp);
    pts->dx=1./64.;
    pts->nx=32;
    loc_t *loc=pts2loc(pts);
    double *phi_pts=calloc(loc->nloc, sizeof(double));
    double displacex,displacey,scale;
    displacex=0;
    displacey=0;
    scale=1.2;//.414065;
    int wrap=0;
    
    PROP_GRD_PTS2_T *propdata= calloc(1, sizeof(PROP_GRD_PTS2_T));
    propdata->mapin=screen;
    propdata->pts=pts;
    propdata->phiout=phi_pts;
    propdata->alpha=-1;
    propdata->displace[0]=displacex;
    propdata->displace[1]=displacey;
    propdata->scale=scale;
    propdata->wrap=wrap;
    propdata->saend=pts->nsa;
    //prop_index(propdata);
    PINIT(propdata->mutex);
   
    read_self_cpu();
    tic;
    prop_grid_pts(screen, pts, phi_pts, -1, displacex, displacey,
		  scale, wrap, 0,0);
    
    prop_grid_pts(screen, pts, phi_pts, -1, displacex, displacey,
		  scale, wrap, 0,0);

    prop_grid_pts(screen, pts, phi_pts, -1, displacex, displacey,
		  scale, wrap, 0,0);

    prop_grid_pts(screen, pts, phi_pts, -1, displacex, displacey,
		  scale, wrap, 0,0);

    prop_grid_pts(screen, pts, phi_pts, -1, displacex, displacey,
		  scale, wrap, 0,0);

    prop_grid_pts(screen, pts, phi_pts, -1, displacex, displacey,
		  scale, wrap, 0,0);

    prop_grid_pts(screen, pts, phi_pts, -1, displacex, displacey,
		  scale, wrap, 0,0);

    prop_grid_pts(screen, pts, phi_pts, -1, displacex, displacey,
		  scale, wrap, 0,0);

    toc("prop_grid_pts");
    info2("cpu: %.2f\n", read_self_cpu());
    tic;
    propdata->isa=0;
    CALL(prop_grid_pts_wrap2, propdata, nthread, 0);

    propdata->isa=0;
    CALL(prop_grid_pts_wrap2, propdata, nthread, 0);

    propdata->isa=0;
    CALL(prop_grid_pts_wrap2, propdata, nthread, 0);
    
    propdata->isa=0;
    CALL(prop_grid_pts_wrap2, propdata, nthread, 0);

    propdata->isa=0;
    CALL(prop_grid_pts_wrap2, propdata, nthread, 0);

    propdata->isa=0;
    CALL(prop_grid_pts_wrap2, propdata, nthread, 0);

    propdata->isa=0;
    CALL(prop_grid_pts_wrap2, propdata, nthread, 0);
    
    propdata->isa=0;
    CALL(prop_grid_pts_wrap2, propdata, nthread, 0);
    toc("prop_grid_pts_wrap");
    info2("cpu: %.2f\n", read_self_cpu());
}
int main(int argc, char** argv){
    int nthread;
    if(argc>1){
	nthread=strtol(argv[1],NULL,10);
    }else{
	nthread=2;
    }
    THREAD_POOL_INIT(nthread);
    test_speed(nthread);
    
    test_accuracy();
}




#include "../lib/aos.h"
#ifdef __clang__
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wvariadic-macros"
#pragma clang diagnostic ignored "-Wc99-extensions"
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wpadded"
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#elif defined __GNUC__
#pragma GCC diagnostic ignored "-Wvariadic-macros"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wpadded"
#endif
TIC;

static void test_accuracy(int argc, char **argv){
    double displacex=0.01;
    double displacey=0.05;
    double scale=1.01;/*.414065; */
    int wrap=0;  

    double D=30;
    double D2=32;
    int save=1;
    if(argc>1){
	scale=strtod(argv[1], 0);
    }
    if(argc>2){
	displacex=strtod(argv[2], 0);
    }
    if(argc>3){
	displacey=strtod(argv[3], 0);
    }
    if(argc>4){
	wrap=strtol(argv[4], 0, 10);
    }
    if(argc>5){
	D=strtod(argv[5], 0);
    }
    if(argc>6){
	D2=strtod(argv[6], 0);
    }
    if(argc>7){
	save=strtol(argv[7], 0, 10);
    }

    double dx=1/64.; 
    double dsa=0.5;
    map_t *screen=mapnew(D2/dx, D2/dx, dx, dx, 0);
    dset((dmat*)screen, 1);
    long len=screen->nx*screen->ny;
    for(long iy=0; iy<screen->ny; iy++){
	for(long ix=0; ix<screen->nx; ix++){
	    screen->p[ix+iy*screen->nx]=sin((double)ix/screen->nx*2*M_PI)*sin((double)iy/screen->ny*2*M_PI);
	}
    }

    /*loc for the map */
    loc_t *locin=mksqloc(screen->nx,screen->ny,dx,dx,screen->ox,screen->oy);
    
    //pts_t *pts=realloc(mkcirloc(D, 1./2.), sizeof(pts_t));
    pts_t *pts=realloc(mksqloc_auto(D/dsa, D/dsa, dsa, dsa), sizeof(pts_t));
    pts->dx=dx;
    pts->dy=dx;
    pts->nx=dsa/dx;
    pts->ny=pts->nx;
    loc_t *loc=pts2loc(pts);

    loc_create_map(locin);
    loc_create_stat(loc);
    loc_create_map(loc);
    locstat_t *locstat=loc->stat;
  
    map_t *screen2=mapnew2(locin->map);
    dset((dmat*)screen2, NAN);
    loc_embed(screen2, locin, screen->p);


    double *phi_pts, *phi_loc, *phi_stat;
    double *phi_loc2loc, *phi_h, *phi_cub, *phi_cub2, *phi_cubh;

    double cubic=0.3;
    int ii;
	    
    info("displacex=%g, displacey=%g, scale=%g, wrap=%d\n",
	 displacex, displacey,scale,wrap);
	
    double diff1, diff2,diff3,diff14,diff15;
    diff1=0;
    diff2=0;
    diff3=0;
    diff14=0;
    diff15=0;
	
    phi_pts=calloc(loc->nloc, sizeof(double));
    phi_loc=calloc(loc->nloc, sizeof(double));
    phi_stat=calloc(loc->nloc, sizeof(double));
    phi_loc2loc=calloc(loc->nloc, sizeof(double));

    map_t *map1=mapnew2(loc->map);
    prop_grid_map(screen, map1, -2, displacex, displacey, scale, wrap, 0,0);
    tic;
    prop_grid_map(screen, map1, 1, displacex, displacey, scale, wrap, 0,0);
    toc("map\t\t");

    prop_grid_pts(screen, pts, NULL,phi_pts, -2, displacex, displacey, scale, wrap, 0,0);
    tic;
    prop_grid_pts(screen, pts, NULL,phi_pts,  1, displacex, displacey, scale, wrap, 0,0);
    toc("pts\t\t");

    prop_grid_stat(screen, locstat, phi_stat, -2,displacex, displacey, scale, wrap, 0,0);
    tic;
    prop_grid_stat(screen, locstat, phi_stat , 1, displacex, displacey, scale, wrap, 0,0);
    toc("stat\t");tic;

    prop_grid(screen, loc, NULL,phi_loc, -2,displacex, displacey, scale, wrap, 0,0);
    tic;
    prop_grid(screen, loc, NULL,phi_loc,  1,displacex, displacey, scale, wrap, 0,0);
    toc("loc\t\t");


    prop_nongrid(locin, screen->p,loc, NULL,phi_loc2loc, -2,displacex, displacey, scale,0,0);
    toc("nongrid\t"); tic;
    prop_nongrid(locin, screen->p,loc, NULL,phi_loc2loc, 1,displacex, displacey, scale,0,0);
    toc("nongrid\t");
	
    phi_h=calloc(loc->nloc,sizeof(double));
 
    tic;
    dsp *hfor=mkh(locin, loc, NULL,displacex, displacey, scale,0,0);
    toc("mkh\t\t\t");
    dspmulvec(phi_h,hfor, screen->p, 'n', -2);
    tic;
    dspmulvec(phi_h,hfor, screen->p, 'n', 1);
    toc("mul h\t\t");


    phi_cub=calloc(loc->nloc,sizeof(double));
    phi_cub2=calloc(loc->nloc,sizeof(double));
    double *phi_cub3=calloc(loc->nloc,sizeof(double));
    phi_cubh=calloc(loc->nloc, sizeof(double));

    prop_nongrid_cubic(locin,screen->p,loc,NULL,phi_cub,-2, displacex, displacey, scale, cubic,0,0);
    tic;
    prop_nongrid_cubic(locin,screen->p,loc,NULL,phi_cub,1, displacex, displacey, scale, cubic,0,0);
    toc("nongrid, cubic\t");
    prop_grid_cubic(screen, loc, NULL,phi_cub2, -2,displacex, displacey, scale,  cubic, 0,0);
    tic;
    prop_grid_cubic(screen, loc, NULL,phi_cub2, 1,displacex, displacey, scale,  cubic, 0,0);
    toc("grid,  cubic\t");
    prop_grid_cubic(screen2, loc, NULL,phi_cub3, -2,displacex, displacey, scale,  cubic, 0,0);
    tic;
    prop_grid_cubic(screen2, loc, NULL,phi_cub3, 1,displacex, displacey, scale,  cubic, 0,0);
    toc("grid2, cubic\t");
 
    dsp *hforcubic;
    tic;
    hforcubic=mkh(locin, loc, NULL, displacex, displacey, scale, 1, cubic);
    toc("mkh cubic \t\t");
    dspmulvec(phi_cubh, hforcubic,screen->p,'n',-2);
    tic;
    dspmulvec(phi_cubh, hforcubic,screen->p,'n',1);
    toc("cubic mul h\t\t");
    double diffc12=0,diff45=0,diff46=0,diff47=0;
    for(ii=0; ii<loc->nloc; ii++){
	diff1+=fabs(phi_loc[ii]-phi_pts[ii]);
	diff2+=fabs(phi_stat[ii]-phi_loc[ii]);
	diff3+=fabs(phi_stat[ii]-phi_pts[ii]);
	diff14+=fabs(phi_loc2loc[ii]-phi_pts[ii]);
	diff15+=fabs(phi_h[ii]-phi_pts[ii]);
	diff45+=fabs(phi_loc2loc[ii]-phi_h[ii]);
	diffc12+=fabs(phi_cub[ii]-phi_cubh[ii]);
	diff46+=fabs(phi_cub[ii]-phi_cub2[ii]);
	diff47+=fabs(phi_cub[ii]-phi_cub3[ii]);
    }
    info2("(pts-loc)=\t%g\n(loc-stat)=\t%g\n(stat-pts)=\t%g\n"
	  "(loc2loc-pts)=\t%g\n(h-pts)=\t%g\n"
	  "(loc2loc-h)=\t%g\n"
	  "(cub:h-loc2loc)=\t%g\n"
	  "(cub:map2loc-loc2loc)=\t%g\n"
	  "(cub:locmap2loc-loc2loc=\t%g\n"
	  ,diff1, diff2,diff3,diff14,diff15,diff45,diffc12, 
	  diff46, diff47);

//    exit(0);
    if(!save) return;
    mapwrite(screen2,"accphi_screen2");
    mapwrite(screen,"accphi_screen");
    locwrite((loc_t*)pts,"accphi_pts");
    locwrite(loc,"accphi_loc");
    locwrite(locin, "accphi_locin");
    loc_create_map_npad(locin, 0, 0, 0);
    mapwrite(locin->map, "accphi_locin_map");
    mapwrite(loc->map, "accphi_loc_map");
    mapwrite(map1, "accphi_map2map.bin");
    writedbl(phi_pts,loc->nloc,1,"accphi_pts1.bin");
    writedbl(phi_loc,loc->nloc,1,"accphi_loc0.bin");
    writedbl(phi_stat,loc->nloc,1,"accphi_stat.bin");
    writedbl(phi_loc2loc,loc->nloc,1,"accphi_loc2loc.bin");
    writedbl(phi_h,loc->nloc,1,"accphi_loc2h.bin");
    writedbl(phi_cub,loc->nloc,1,"accphi_cub_loc2loc.bin");
    writedbl(phi_cub2,loc->nloc,1,"accphi_cub_map2loc.bin");
    writedbl(phi_cub3,loc->nloc,1,"accphi_cub_locmap2loc.bin");
    writedbl(phi_cubh,loc->nloc,1,"accphi_cub_loc2h.bin");
    info("saved\n");

    writedbl(phi_pts,loc->nloc,1,"accphi_pts.bin");
    writedbl(phi_cub,loc->nloc,1,"accphi_cub.bin");

    writebin(hfor, "accphi_hfor");
    writebin(hforcubic, "accphi_cub_hfor");
    dspfree(hfor);
    dspfree(hforcubic);
    free(phi_loc); free(phi_stat);
    free(phi_loc2loc);

    free(phi_pts);
    free(phi_h); 
    free(phi_cub); 
    free(phi_cubh);
	
    cellfree(screen);
    locfree(locin);
}
/*
  matched grip, non wrap, test ok.
*/

int main(int argc, char** argv){
    test_accuracy(argc, argv);
}

/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

static void test_accuracy(int argc, char** argv){
    real displacex=0.01;
    real displacey=0.05;
    real scale=0.7;/*.414065; */
    int wrap=1;

    real D=10;
    real D2=12;
    int save=0;
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

    real dx=1/64.;
    real dsa=0.5;
    map_t* screen=mapnew(D2/dx, D2/dx, dx, dx);
    dset(DMAT(screen), 1);
    dmat* tmp=dnew(screen->nx, 1);
    for(long ix=0; ix<screen->nx; ix++){
        P(tmp, ix)=sin((real)ix/screen->nx*2*M_PI);
    }
    for(long iy=0; iy<screen->ny; iy++){
        for(long ix=0; ix<screen->nx; ix++){
            P(screen, ix, iy)=P(tmp, ix)*P(tmp, iy);
        }
    }
	const real ht=10000;
	const real hs=ht/(1-scale);
	screen->ht=ht;
	
    //pts_t *pts=realloc(mkannloc(D, 0, 1./2.,0), sizeof(pts_t));
    pts_t* pts=(pts_t*)realloc(mksqloc_auto(D/dsa, D/dsa, dsa, dsa), sizeof(pts_t));
    pts->dx=dx;
    pts->dy=dx;
    pts->nxsa=dsa/dx;
    pts->nysa=pts->nxsa;
    loc_t* locin=pts2loc(pts);
	locin->ht=ht;
    loc_create_map(locin);
	loc_t* locout=pts2loc(pts);
	locout->ht=0;
	loc_create_stat(locout);
	loc_create_map(locout);
	locstat_t* locstat=locout->stat;
	
    /*loc for the map */
    loc_t* locin2=mksqloc(screen->nx, screen->ny, dx, dx, screen->ox, screen->oy);
	locin2->ht=ht;
    loc_create_map(locin2);
    

    map_t* screen2=mapnew2(locin2->map);
    dset(DMAT(screen2), NAN);
    loc_embed(screen2, locin2, screen);
	screen2->ht=locin2->ht;
    dmat *phi_h=NULL, *phi_cubh=NULL;

    
    int ii;

    info("displacex=%g, displacey=%g, scale=%g, wrap=%d\n",
        displacex, displacey, scale, wrap);

    real diff1, diff2, diff3, diff14, diff15;
    diff1=0;
    diff2=0;
    diff3=0;
    diff14=0;
    diff15=0;

    dmat* phi_pts=dnew(locout->nloc, 1);
    dmat* phi_loc=dnew(locout->nloc, 1);
    dmat* phi_stat=dnew(locout->nloc, 1);
    dmat* phi_loc2loc=dnew(locout->nloc, 1);

    map_t* map1=mapnew2(locout->map);

    prop(&(propdata_t){.mapin=screen, .locout=locout, .phiout=P(phi_loc), .alpha=-2, .hs=hs, .shiftx=displacex, .shifty=displacey, .wrap=wrap});
    tic;
    prop(&(propdata_t){.mapin=screen, .locout=locout, .phiout=P(phi_loc), .alpha=1, .hs=hs, .shiftx=displacex, .shifty=displacey, .wrap=wrap});
    toc("loc\t");
	prop(&(propdata_t){.mapin=screen, .ptsout=pts, .phiout=P(phi_pts), .alpha=-2, .hs=hs, .shiftx=displacex, .shifty=displacey, .wrap=wrap});
	tic;
	prop(&(propdata_t){.mapin=screen, .ptsout=pts, .phiout=P(phi_pts), .alpha=1, .hs=hs, .shiftx=displacex, .shifty=displacey, .wrap=wrap});
    toc("pts\t");

	prop(&(propdata_t){.mapin=screen, .mapout=map1, .alpha=-2, .hs=hs, .shiftx=displacex, .shifty=displacey, .wrap=wrap});
    tic;
	prop(&(propdata_t){.mapin=screen, .mapout=map1, .alpha=1, .hs=hs, .shiftx=displacex, .shifty=displacey, .wrap=wrap});
    toc("map\t");

	prop(&(propdata_t){.mapin=screen, .ostat=locstat, .phiout=P(phi_stat), .alpha=-2, .hs=hs, .shiftx=displacex, .shifty=displacey, .wrap=wrap});
    tic;
	prop(&(propdata_t){.mapin=screen, .ostat=locstat, .phiout=P(phi_stat), .alpha=1, .hs=hs, .shiftx=displacex, .shifty=displacey, .wrap=wrap});
    toc("stat\t");tic;


	prop(&(propdata_t){.locin=locin, .phiin=P(screen->dmat), .locout=locout, .phiout=P(phi_loc2loc), .alpha=-2, .hs=hs, .shiftx=displacex, .shifty=displacey, .wrap=wrap});
    toc("nongrid\t"); tic;
	prop(&(propdata_t){.locin=locin, .phiin=P(screen->dmat), .locout=locout, .phiout=P(phi_loc2loc), .alpha=1, .hs=hs, .shiftx=displacex, .shifty=displacey, .wrap=wrap});
    toc("nongrid\t");


    phi_h=dnew(locout->nloc, 1);
    tic;
    dsp* hfor=mkh(locin2, locout, displacex, displacey, scale, 0);
    toc("mkh\t\t");
    dspmv(phi_h, hfor, dmat_cast(screen), 'n', -2);
    tic;
    dspmv(phi_h, hfor, dmat_cast(screen), 'n', 1);
    toc("mul h\t");
	//for cubic
	const real cubic=0.3;
	screen->iac=cubic;
	locin->iac=cubic;
	screen2->iac=cubic;
    dmat* phi_cub=dnew(locout->nloc, 1);
    dmat* phi_cub2=dnew(locout->nloc, 1);
    dmat* phi_cub3=dnew(locout->nloc, 1);
    dmat* phi_cub4=dnew(locout->nloc, 1);
    phi_cubh=dnew(locout->nloc, 1);
	prop(&(propdata_t){.locin=locin2, .phiin=P(screen->dmat), .locout=locout, .phiout=P(phi_cub), .alpha=-2, .hs=hs, .shiftx=displacex, .shifty=displacey});
    tic;
	prop(&(propdata_t){.locin=locin2, .phiin=P(screen->dmat), .locout=locout, .phiout=P(phi_cub), .alpha=1, .hs=hs, .shiftx=displacex, .shifty=displacey});
    toc("nongrid, cubic\t");
	prop(&(propdata_t){.mapin=screen, .locout=locout, .phiout=P(phi_cub2), .alpha=-2, .hs=hs, .shiftx=displacex, .shifty=displacey});
    tic;
	prop(&(propdata_t){.mapin=screen, .locout=locout, .phiout=P(phi_cub2), .alpha=1, .hs=hs, .shiftx=displacex, .shifty=displacey});
    toc("grid,  cubic\t");
	prop(&(propdata_t){.mapin=screen2, .locout=locout, .phiout=P(phi_cub3), .alpha=-2, .hs=hs, .shiftx=displacex, .shifty=displacey});
    tic;
	prop(&(propdata_t){.mapin=screen2, .locout=locout, .phiout=P(phi_cub3), .alpha=1, .hs=hs, .shiftx=displacex, .shifty=displacey});
    toc("grid2, cubic\t");
	prop(&(propdata_t){.mapin=screen, .ostat=locstat, .phiout=P(phi_cub4), .alpha=-2, .hs=hs, .shiftx=displacex, .shifty=displacey});
    tic;
	prop(&(propdata_t){.mapin=screen, .ostat=locstat, .phiout=P(phi_cub4), .alpha=1, .hs=hs, .shiftx=displacex, .shifty=displacey});
    toc("grid  2stat, cubic\t");
    dsp* hforcubic;
    tic;
    hforcubic=mkh_cubic(locin2, locout, displacex, displacey, scale, 0, cubic);
    toc("mkh cubic \t\t");
    dspmv(phi_cubh, hforcubic, dmat_cast(screen), 'n', -2);
    tic;
    dspmv(phi_cubh, hforcubic, dmat_cast(screen), 'n', 1);
    toc("cubic mul h\t\t");
    real diffc12=0, diff45=0, diff46=0, diff47=0;
    for(ii=0; ii<locout->nloc; ii++){
        diff1+=fabs(P(phi_loc,ii)-P(phi_pts,ii));
        diff2+=fabs(P(phi_stat,ii)-P(phi_loc,ii));
        diff3+=fabs(P(phi_stat,ii)-P(phi_pts,ii));
        diff14+=fabs(P(phi_loc2loc,ii)-P(phi_pts,ii));
        diff15+=fabs(P(phi_h,ii)-P(phi_pts,ii));
        diff45+=fabs(P(phi_loc2loc,ii)-P(phi_h,ii));
        diffc12+=fabs(P(phi_cub,ii)-P(phi_cubh,ii));
        diff46+=fabs(P(phi_cub,ii)-P(phi_cub2,ii));
        diff47+=fabs(P(phi_cub,ii)-P(phi_cub3,ii));
    }
    info("(pts-loc)=\t%g\n(stat-loc)=\t%g\n(stat-pts)=\t%g\n"
        "(loc2loc-pts)=\t%g\n(h-pts)=\t%g\n"
        "(loc2loc-h)=\t%g\n"
        "(cub:h-loc2loc)=\t%g\n"
        "(cub:map2loc-loc2loc)=\t%g\n"
        "(cub:locmap2loc-loc2loc=\t%g\n"
        , diff1, diff2, diff3, diff14, diff15, diff45, diffc12,
        diff46, diff47);

  //    exit(0);
    if(save){
        mapwrite(screen2, "accphi_screen2");
        mapwrite(screen, "accphi_screen");
        locwrite(pts, "accphi_pts");
        locwrite(locout, "accphi_loc");
        locwrite(locin2, "accphi_locin");
        loc_create_map_npad(locin2, 0, 0, 0);
        mapwrite(locin2->map, "accphi_locin_map");
        mapwrite(locout->map, "accphi_loc_map");
        mapwrite(map1, "accphi_map2map.bin");
        /*
        writedbl(phi_pts,locout->nloc,1,"accphi_pts1.bin");
        writedbl(phi_loc,locout->nloc,1,"accphi_loc0.bin");
        writedbl(phi_stat,locout->nloc,1,"accphi_stat.bin");
        writedbl(phi_loc2loc,locout->nloc,1,"accphi_loc2loc.bin");
        writedbl(phi_h,locout->nloc,1,"accphi_loc2h.bin");
        writedbl(phi_cub,locout->nloc,1,"accphi_cub_loc2loc.bin");
        writedbl(phi_cub2,locout->nloc,1,"accphi_cub_map2loc.bin");
        writedbl(phi_cub3,locout->nloc,1,"accphi_cub_locmap2loc.bin");
        writedbl(phi_cub4,locout->nloc,1,"accphi_cub_locmap2stat.bin");
        writedbl(phi_cubh,locout->nloc,1,"accphi_cub_loc2h.bin");
        info("saved\n");

        writedbl(phi_pts,locout->nloc,1,"accphi_pts.bin");
        writedbl(phi_cub,locout->nloc,1,"accphi_cub.bin");
        */
        writebin(hfor, "accphi_hfor");
        writebin(hforcubic, "accphi_cub_hfor");

    }
    dspfree(hfor);
    dspfree(hforcubic);
    free(phi_loc); free(phi_stat);
    free(phi_loc2loc);

    free(phi_pts);
    dfree(phi_h);
    free(phi_cub);
    dfree(phi_cubh);

    cellfree(screen);
    locfree(locin2);
}

int main(int argc, char** argv){
    OMPTASK_SINGLE
        test_accuracy(argc, argv);
}

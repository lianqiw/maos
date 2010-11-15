/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#include "locbin.h"
#include "matbin.h"
#include "cell.h"
/**
   \file locbin.c
   i/o functions for loc_t, map_t.
*/
/**
   Read in the loc_t data after magic and dimension has been read and verified.
*/
loc_t *locreaddata2(file_t *fp, long nx, long ny){
    if(ny!=2){
	error("This is not a LOC file\n");
    }
    loc_t *out=calloc(1, sizeof(loc_t));
    out->nloc=nx;
    out->locx=malloc(sizeof(double)*nx);
    zfread(out->locx, sizeof(double), nx, fp);
    out->locy=malloc(sizeof(double)*nx);
    zfread(out->locy, sizeof(double), nx, fp);
    if((out->locx[2]+out->locx[0]-2.*out->locx[1])<1.e-10){
	out->dx=out->locx[1]-out->locx[0];
    }else{
	out->dx=NAN;
	warning("Unable to determine dx from the loc file %s. Please set manually.\n", fp->fn);
    }
    out->map=NULL;
    return out;
}
/**
   Verify the magic, dimension and read in the loc_t by calling locreaddata2().
 */
loc_t *locreaddata(file_t *fp){
    uint32_t magic;
    zfread(&magic, sizeof(uint32_t), 1, fp);
    if(magic!=0x6402){
	error("This is not a LOC file\n");
    }
    uint64_t nx,ny;
    zfread(&nx, sizeof(uint64_t), 1, fp);
    zfread(&ny, sizeof(uint64_t), 1, fp);
    loc_t *out;
    if(nx==0 || ny==0){
	out=NULL;
    }else{
	out=locreaddata2(fp, nx, ny);
    }
    return out;
}
void locwritedata(file_t *fp, const loc_t *loc){
    uint32_t magic=M_DBL;
    zfwrite(&magic, sizeof(uint32_t),1,fp);
    if(loc){
	uint64_t nx=loc->nloc;
	uint64_t ny=2;
	zfwrite(&nx,sizeof(uint64_t),1,fp);
	zfwrite(&ny,sizeof(uint64_t),1,fp);
	zfwrite(loc->locx, sizeof(double),loc->nloc,fp);
	zfwrite(loc->locy, sizeof(double),loc->nloc,fp);
    }else{
	uint64_t nx=0;
	zfwrite(&nx, sizeof(uint64_t),1,fp);
	zfwrite(&nx, sizeof(uint64_t),1,fp);
    }
}
void locwrite(const loc_t *loc, const char *format,...){
    format2fn;
    file_t *fp=zfopen(fn,"wb");
    locwritedata(fp, loc);
    zfclose(fp);
}

void locarrwrite(loc_t ** loc, int nloc, const char *format,...){
    format2fn;
    file_t *fp=zfopen(fn,"wb");
    uint32_t magic=MC_DBL;
    uint64_t nx=nloc;
    uint64_t ny=1;
    zfwrite(&magic, sizeof(uint32_t),1, fp);
    zfwrite(&nx,sizeof(uint64_t),1, fp);
    zfwrite(&ny,sizeof(uint64_t),1, fp);
    for(unsigned long iloc=0; iloc<nloc; iloc++){
	locwritedata(fp, loc[iloc]);
    }
    zfclose(fp);
}
loc_t *locread(const char *format,...){
    format2fn;
    file_t *fp=zfopen(fn, "rb");
    loc_t *loc=locreaddata(fp);
    zfclose(fp);
    return loc;
}
loc_t ** locarrread(int *nloc, const char*format,...){
    format2fn;
    file_t *fp=zfopen(fn,"rb");
    uint32_t magic;
    zfread(&magic, sizeof(uint32_t), 1, fp);
    if(magic!=MC_DBL){
	error("This is not a locarr file");
    }
    uint64_t nx,ny;
    zfread(&nx, sizeof(uint64_t), 1, fp);
    zfread(&ny, sizeof(uint64_t), 1, fp);
    *nloc=nx*ny;
    loc_t **locarr=calloc(nx*ny, sizeof(loc_t*));
    for(long ix=0; ix<nx*ny; ix++){
	locarr[ix]=locreaddata(fp);
    }
    zfclose(fp);
    return locarr;
}
static void sqmapwritedata(file_t *fp, const map_t *map){
    uint32_t magic=M_DBL;
    zfwrite(&magic, sizeof(uint32_t),1,fp);
    if(map){
	uint64_t nx=map->nx;
	uint64_t ny=map->ny;
	zfwrite(&nx,sizeof(uint64_t),1,fp);
	zfwrite(&ny,sizeof(uint64_t),1,fp);
	zfwrite(map->p, sizeof(double),nx*ny,fp);
    }else{
	uint64_t nx=0;
	zfwrite(&nx, sizeof(uint64_t),1,fp);
	zfwrite(&nx, sizeof(uint64_t),1,fp);
    }
}

map_t *sqmapread(const char *format,...){
    format2fn;
    map_t *sqmap=calloc(1, sizeof(map_t));
    dcell *tmp=dcellread("%s",fn); 
    /*
      File should contain two cells. The first cell contains the information. 
      Cell 1: 5x1 vector, dx, dy, ox, oy, height.
      Cell 2: nxm array.
     */

    if(fabs(tmp->p[0]->p[0]-tmp->p[0]->p[1])>1.e-14){
	error("Map should be square\n");
    }
    sqmap->dx=tmp->p[0]->p[0];
    sqmap->ox=tmp->p[0]->p[2];
    sqmap->oy=tmp->p[0]->p[3];
    sqmap->h=tmp->p[0]->p[4];
    
    dmat *ampg=dref(tmp->p[1]);
    sqmap->p=ampg->p;
    sqmap->nx=ampg->nx;
    sqmap->ny=ampg->ny;
    dcellfree(tmp);
    dfree_keepdata(ampg);
    if(sqmap->ox/sqmap->dx*2+sqmap->nx > 2
       ||sqmap->oy/sqmap->dx*2+sqmap->ny > 2){
	warning("Ampground %s is not centered.\n",fn);
    }
    return sqmap;
}

rectmap_t *rectmapread(const char *format,...){
    format2fn;
    rectmap_t *rectmap=calloc(1, sizeof(map_t));
    dcell *tmp=dcellread("%s",fn); 
    rectmap->dx=tmp->p[0]->p[0];
    rectmap->dy=tmp->p[0]->p[1];
    rectmap->ox=tmp->p[0]->p[2];
    rectmap->oy=tmp->p[0]->p[3];

    rectmap->txdeg=tmp->p[0]->p[4];
    rectmap->tydeg=tmp->p[0]->p[5];
    rectmap->ftel =tmp->p[0]->p[6];
    rectmap->fexit=tmp->p[0]->p[7];
    rectmap->fsurf=tmp->p[0]->p[8];
 
    dmat *ampg=dref(tmp->p[1]);
    rectmap->p=ampg->p;
    rectmap->nx=ampg->nx;
    rectmap->ny=ampg->ny;
    dcellfree(tmp);
    dfree_keepdata(ampg);
    if(rectmap->ox/rectmap->dx*2+rectmap->nx > 2
       ||rectmap->oy/rectmap->dy*2+rectmap->ny > 2){
	warning("Ampground %s is not centered.\n",fn);
    }
    return rectmap;
}

void sqmapwrite(const map_t *map, const char *format,...){
    format2fn;
    file_t *fp=zfopen(fn,"wb");
    sqmapwritedata(fp, map);
    zfclose(fp);
}
void sqmaparrwrite(map_t ** map, int nmap, 
		   const char *format,...){
    format2fn;
    file_t *fp=zfopen(fn,"wb");
    uint32_t magic=MC_DBL;
    uint64_t nx=2;
    uint64_t ny=nmap;
    zfwrite(&magic, sizeof(uint32_t),1, fp);
    zfwrite(&nx,sizeof(uint64_t),1, fp);
    zfwrite(&ny,sizeof(uint64_t),1, fp);
    dmat *info=dnew(5,1);
    for(int imap=0; imap<nmap; imap++){
	info->p[0]=map[imap]->dx;
	info->p[1]=map[imap]->dx;
	info->p[2]=map[imap]->ox;
	info->p[3]=map[imap]->oy;
	info->p[4]=map[imap]->h;
	dwritedata(fp,info);
	sqmapwritedata(fp, map[imap]);
    }
    dfree(info);
    zfclose(fp);
}
map_t **sqmaparrread(int*nlayer, const char *format,...){
    format2fn;
    dcell *X=dcellread("%s",fn);
    map_t **screens;
    if(X->ny==1 && (X->nx!=2 || X->p[0]->nx*X->p[0]->ny!=5) ){
	warning("Deprecated format\n");
	*nlayer=X->nx;
	screens=calloc(X->nx,sizeof(map_t*));
        double dx=1./64.;//assume this is 1/64
	warning("dx is assumed to be %g\n",dx);
	for(int ilayer=0; ilayer<X->nx; ilayer++){
	    screens[ilayer]=calloc(1, sizeof(map_t));
	    screens[ilayer]->nx=X->p[ilayer]->nx;
	    screens[ilayer]->ny=X->p[ilayer]->ny;
	    screens[ilayer]->p=X->p[ilayer]->p;
	    screens[ilayer]->dx=dx;
	    screens[ilayer]->ox=X->p[ilayer]->nx*(-dx/2.);
	    screens[ilayer]->oy=X->p[ilayer]->ny*(-dx/2.);
	    dfree_keepdata(X->p[ilayer]);
	}
    }else if(X->nx==2){//New format. file contains 2*n
	*nlayer=X->ny;
    	screens=calloc(X->ny, sizeof(map_t*));
	PDCELL(X,pX);
	for(int ilayer=0; ilayer<X->ny; ilayer++){
	    screens[ilayer]=calloc(1, sizeof(map_t));
	    screens[ilayer]->nx=pX[ilayer][1]->nx;
	    screens[ilayer]->ny=pX[ilayer][1]->ny;
	    screens[ilayer]->p =pX[ilayer][1]->p;
	    screens[ilayer]->dx=pX[ilayer][0]->p[0];
	    screens[ilayer]->ox=pX[ilayer][0]->p[2];
	    screens[ilayer]->oy=pX[ilayer][0]->p[3];
	    screens[ilayer]->h=pX[ilayer][0]->p[4];
	    if(fabs(pX[ilayer][0]->p[0]-pX[ilayer][0]->p[1])>1.e-12){
		error("This map doesn't have equal sampling over x and y\n");
	    }
	    dfree_keepdata(X->p[ilayer]);
	}
    }else{
	screens=NULL;
	error("Invalid format\n");
    }
    free(X->p);
    free(X);
    return screens;
}

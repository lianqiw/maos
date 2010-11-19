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
   Verify the magic, dimension and read in the loc_t by calling locreaddata2().
 */
loc_t *locreaddata(file_t *fp, uint32_t magic, char *header0){
    char *header=NULL;
    int free_header;
    if(!magic){
	magic=read_magic(fp, &header);
	free_header=1;
    }else{
	header=header0;
	free_header=0;
    }
    if(magic!=M_DBL){
	error("magic=%x. Expect %x\n", magic, M_DBL);
    }
    double dx=search_header_num(header,"dx");
    if(free_header){
	free(header);
    }
    uint64_t nx,ny;
    zfread(&nx, sizeof(uint64_t), 1, fp);
    zfread(&ny, sizeof(uint64_t), 1, fp);
    loc_t *out;
    if(nx==0 || ny==0){
	out=NULL;
    }else{
	out=calloc(1, sizeof(loc_t));
	out->nloc=nx;
	out->locx=malloc(sizeof(double)*nx);
	zfread(out->locx, sizeof(double), nx, fp);
	out->locy=malloc(sizeof(double)*nx);
	zfread(out->locy, sizeof(double), nx, fp);
	if(fabs(dx)<1e-100){//dx is not available.
	    for(long i=0; i<out->nloc-1; i++){//we assume the rows are continuous.
		if(out->locy[i+1]>out->locy[i]){
		    dx=out->locy[i+1]-out->locy[i];
		}
	    }
	}else{
	    out->dx=dx;
	}
	out->map=NULL;
    }
    return out;
}
/**
   Read a loc_t from file. dx is contained in header
*/
loc_t *locread(const char *format,...){
    format2fn;
    file_t *fp=zfopen(fn, "rb");
    loc_t *loc=locreaddata(fp, 0, NULL);
    zfclose(fp);
    return loc;
}
/**
   Read an array of loc_t form file.
 */
loc_t ** locarrread(int *nloc, const char*format,...){
    format2fn;
    file_t *fp=zfopen(fn,"rb");
    uint32_t magic=read_magic(fp, NULL);
    if(magic!=MC_DBL && magic!=MCC_ANY){
	error("This is not a locarr file");
    }
    uint64_t nx,ny;
    zfread(&nx, sizeof(uint64_t), 1, fp);
    zfread(&ny, sizeof(uint64_t), 1, fp);
    *nloc=nx*ny;
    loc_t **locarr=calloc(nx*ny, sizeof(loc_t*));
    for(long ix=0; ix<nx*ny; ix++){
	locarr[ix]=locreaddata(fp, 0, NULL);
    }
    zfclose(fp);
    return locarr;
}
void locwritedata(file_t *fp, const loc_t *loc){
    char header[80];
    snprintf(header,80,"dx=%.15g\n",loc->dx);
    write_header(header,fp);
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
/**
   Read data of a map_t.
*/
map_t *sqmapreaddata(file_t *fp, uint32_t magic, char *header0){
    char *header=NULL;
    int free_header;
    if(!magic){
	magic=read_magic(fp, &header);
	free_header=1;
    }else{
	header=header0;
	free_header=0;
    }
    if(magic!=M_DBL){
	error("magic=%x. Expect %x\n", magic, M_DBL);
    }
    map_t *map=calloc(1, sizeof(map_t));
    if(header){//there is a header. so we don't need cell array.
	map->ox=search_header_num(header,"ox");
	map->oy=search_header_num(header,"oy");
	map->dx=search_header_num(header,"dx");
	map->h=search_header_num(header,"h");
	map->vx=search_header_num(header,"vx");
	map->vy=search_header_num(header,"vy");
	dmat *tmp=dreaddata(fp, magic);
	map->p=tmp->p;
	map->nx=tmp->nx;
	map->ny=tmp->ny;

	dfree_keepdata(tmp);
	if(free_header){
	    free(header);
	}
    }else{//there is no header. we require cell.
	dcell *tmp=dcellreaddata(fp, magic); 
	/*
	  File should contain two cells. The first cell contains the information. 
	  Cell 1: 5x1 vector, dx, dy, ox, oy, height.
	  Cell 2: nxm array.
	*/
	warning("Please update %s to newer format with headers\n", fp->fn);
	if(fabs(tmp->p[0]->p[0]-tmp->p[0]->p[1])>1.e-14){
	    error("Map should be square\n");
	}
	map->dx=tmp->p[0]->p[0];
	map->ox=tmp->p[0]->p[2];
	map->oy=tmp->p[0]->p[3];
	map->h=tmp->p[0]->p[4];
    
	dmat *ampg=dref(tmp->p[1]);
	map->p=ampg->p;
	map->nx=ampg->nx;
	map->ny=ampg->ny;
	dcellfree(tmp);
	dfree_keepdata(ampg);
    }
    if(map->ox/map->dx*2+map->nx > 2 ||map->oy/map->dx*2+map->ny > 2){
	warning("map_t %s is not centered.\n",fp->fn);
    }
    return map;
}
/**
   Read map_t from file.
*/
map_t *sqmapread(const char *format,...){
    format2fn;
    file_t *fp=zfopen(fn,"rb");
    map_t *map=sqmapreaddata(fp, 0, NULL);
    zfclose(fp);
    return map;
}
/**
   Read a rectmap data from file_t
*/
rectmap_t *rectmapreaddata(file_t *fp){
    rectmap_t *map=calloc(1, sizeof(map_t));
    char *header=NULL;
    uint32_t magic=read_magic(fp, &header);
    if(header){//there is a header. so we don't need cell array.
	if(magic!=M_DBL){
	    error("Invalid format %x for rectmap_t: %s\n", magic, fp->fn);
	}
	map->dx=search_header_num(header,"dx");
	map->dy=search_header_num(header,"dy");
	map->ox=search_header_num(header,"ox");
	map->oy=search_header_num(header,"oy");
	map->txdeg=search_header_num(header,"txdeg");
	map->tydeg=search_header_num(header,"tydeg");
	map->ftel =search_header_num(header,"ftel");
	map->fexit=search_header_num(header,"fexit");
	map->fsurf=search_header_num(header,"fsurf");
	dmat *tmp=dreaddata(fp, magic);
	map->p=tmp->p;
	map->nx=tmp->nx;
	map->ny=tmp->ny;
	dfree_keepdata(tmp);
	free(header);
    }else{
	if(magic!=MC_DBL && magic!=MCC_ANY){
	    error("Invalid format %x for rectmap_t: %s\n", magic, fp->fn);
	}
	warning("Please update %s to newer format with headers\n", fp->fn);
	dcell *tmp=dcellreaddata(fp, magic); 
	map->dx=tmp->p[0]->p[0];
	map->dy=tmp->p[0]->p[1];
	map->ox=tmp->p[0]->p[2];
	map->oy=tmp->p[0]->p[3];

	map->txdeg=tmp->p[0]->p[4];
	map->tydeg=tmp->p[0]->p[5];
	map->ftel =tmp->p[0]->p[6];
	map->fexit=tmp->p[0]->p[7];
	map->fsurf=tmp->p[0]->p[8];
 
	dmat *ampg=dref(tmp->p[1]);
	map->p=ampg->p;
	map->nx=ampg->nx;
	map->ny=ampg->ny;
	dcellfree(tmp);
	dfree_keepdata(ampg);
    }
    if(map->ox/map->dx*2+map->nx > 2
       ||map->oy/map->dy*2+map->ny > 2){
	warning("map_t %s is not centered.\n",fp->fn);
    }
    return map;
}
/**
   Read a rectmap_t from file.
*/
rectmap_t *rectmapread(const char *format,...){
    format2fn;
    file_t *fp=zfopen(fn,"rb");
    rectmap_t *map=rectmapreaddata(fp);
    zfclose(fp);
    return map;
}
static void sqmapwritedata(file_t *fp, const map_t *map){
    if(map){
	char header[1024];
	snprintf(header,1024,"ox=%.15g\noy=%.15g\ndx=%.15g\nh=%.15g\nvx=%.15g\nvy=%.15g\n",
		 map->ox, map->oy, map->dx, map->h, map->vx, map->vy);
	write_header(header,fp);
    }
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
void sqmapwrite(const map_t *map, const char *format,...){
    format2fn;
    file_t *fp=zfopen(fn,"wb");
    sqmapwritedata(fp, map);
    zfclose(fp);
}
void sqmaparrwrite(map_t ** map, int nmap, const char *format,...){
    format2fn;
    file_t *fp=zfopen(fn,"wb");
    uint32_t magic=MC_DBL;
    uint64_t nx=nmap;
    uint64_t ny=1;
    zfwrite(&magic, sizeof(uint32_t),1, fp);
    zfwrite(&nx,sizeof(uint64_t),1, fp);
    zfwrite(&ny,sizeof(uint64_t),1, fp);
    for(int imap=0; imap<nmap; imap++){
	sqmapwritedata(fp, map[imap]);
    }
    zfclose(fp);
}
map_t **sqmaparrread(int*nlayer, const char *format,...){
    format2fn;
    file_t *fp=zfopen(fn,"rb");
    int has_header=0;

    uint32_t magic=read_magic(fp, NULL);//magic of the cell. 
    if(magic!=MC_DBL && magic!=MCC_ANY){
	error("Invalid file. magic=%x\n",magic);
    }
    map_t **screens;
    uint64_t nx, ny;
    zfreadlarr(fp, 2, &nx, &ny);
    //test the first cell.
    char *header=NULL;
    uint32_t magic2=read_magic(fp, &header);
    if(magic2!=M_DBL){
	error("Wrong data type in the file %s\n", fp->fn);
    }
    if(header){
	has_header=1;
	free(header);
    }

    zfrewind(fp);//start from the begining.

    if(has_header){//there is header.
	//read the cell header.
	magic=read_magic(fp, NULL);
	zfreadlarr(fp, 2, &nx, &ny);
	assert(ny==1);
	screens=calloc(nx, sizeof(map_t*));
	for(int ilayer=0; ilayer<nx; ilayer++){
	    screens[ilayer]=sqmapreaddata(fp, 0, NULL);
	}
	*nlayer=nx;
    }else{//old format.
	dcell *X=dcellreaddata(fp,0);
	warning("Deprecated format: %s\n", fp->fn);
	if(X->ny==1 && (X->nx!=2 || X->p[0]->nx*X->p[0]->ny!=5) ){//oldest format.
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
	}else if(X->nx==2){
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
    }//old format
    zfclose(fp);
    return screens;
}

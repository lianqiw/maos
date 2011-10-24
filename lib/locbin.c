/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
	if(fabs(dx)<EPS || isnan(dx)){/*dx is not available. */
	    for(long i=0; i<out->nloc-1; i++){/*we assume the rows are continuous. */
		if(out->locy[i+1]>out->locy[i]){
		    dx=out->locy[i+1]-out->locy[i];
		    info("Guessing: dx=%g\n", dx);
		    break;
		}
	    }
	}
	out->dx=dx;
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
    if(!iscell(magic)){
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
    write_magic(M_DBL, fp);
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
    uint64_t nx=nloc;
    uint64_t ny=1;
    write_magic(MCC_ANY, fp);
    zfwrite(&nx,sizeof(uint64_t),1, fp);
    zfwrite(&ny,sizeof(uint64_t),1, fp);
    for(unsigned long iloc=0; iloc<nloc; iloc++){
	locwritedata(fp, loc[iloc]);
    }
    zfclose(fp);
}

static void mapwritedata(file_t *fp, map_t *map){
    if(map){
	if(!map->header){
	    char header[1024];
	    snprintf(header,1024,"ox=%.15g\noy=%.15g\ndx=%.15g\nh=%.15g\nvx=%.15g\nvy=%.15g\n",
		     map->ox, map->oy, map->dx, map->h, map->vx, map->vy);
	    map->header=strdup(header);
	}
    }
    dmat *in=(dmat*) map;
    dwritedata(fp, in);
}
void mapwrite(map_t *map, const char *format,...){
    format2fn;
    file_t *fp=zfopen(fn,"wb");
    mapwritedata(fp, map);
    zfclose(fp);
}
void maparrwrite(map_t ** map, int nmap, const char *format,...){
    format2fn;
    file_t *fp=zfopen(fn,"wb");
    uint64_t nx=nmap;
    uint64_t ny=1;
    write_magic(MCC_ANY, fp);
    zfwritelarr(fp, 2, &nx, &ny);
    for(int imap=0; imap<nmap; imap++){
	mapwritedata(fp, map[imap]);
    }
    zfclose(fp);
}

/**
   convert a dmat to map_t.
*/
map_t* d2map(dmat *in){
    map_t *map=realloc(dref(in), sizeof(map_t));
    char *header=in->header;
    if(header){
	map->ox=search_header_num(header,"ox");
	map->oy=search_header_num(header,"oy");
	map->dx=search_header_num(header,"dx");
	map->h =search_header_num(header,"h");
	map->vx=search_header_num(header,"vx");
	map->vy=search_header_num(header,"vy");
    }
    if(isnan(map->ox) || isnan(map->oy) || isnan(map->dx)){
	warning("Valid header is needed to convert dmat to map_t\n");
	map->dx=1;
	map->ox=map->nx/2*map->dx;
	map->oy=map->ny/2*map->dx;
	map->h=map->vx=map->vy=0;
    }
    if(isnan(map->ox) || isnan(map->oy) || isnan(map->dx)){
	warning("The map has no valid ox, oy and dx\n");
    }
    if(isnan(map->h)) map->h=0;
    if(isnan(map->vx)) map->vx=0;
    if(isnan(map->vy)) map->vy=0;
    return map;
}

/**
 * convert a mmap'ed dcell to map_t array
 */
map_t **dcell2map(int *nlayer, dcell *in){
    *nlayer=in->nx*in->ny;
    map_t **map=calloc(in->nx*in->ny, sizeof(map_t*));
    for(long i=0; i<in->nx*in->ny; i++){
	if(!in->p[i]->header){
	    in->p[i]->header=strdup(in->header);
	}
	map[i]=d2map(in->p[i]);
    }
    return map;
}

/**
 * Read map_t from file
 */
map_t *mapread(const char *format, ...){
    format2fn;
    file_t *fp=zfopen(fn,"rb");
    char *header=NULL;
    uint32_t magic=read_magic(fp, &header);
    map_t *map=NULL;
    if(magic==M_DBL){
	dmat *in=dreaddata(fp, magic);
	if(in){
	    in->header=header;
	}else{
	    free(header);
	}
	map=d2map(in);
	dfree(in);
    }else if(iscell(magic)){/*old format. */
	dcell *in=dcellreaddata(fp, magic);
	if(in){
	    in->header=header;
	}else{
	    free(header);
	}

	if(fabs(in->p[0]->p[0]-in->p[0]->p[1])>1.e-14){
	    error("Map should be square\n");
	}
	map=calloc(1, sizeof(map_t));
	map->dx=in->p[0]->p[0];
	map->ox=in->p[0]->p[2];
	map->oy=in->p[0]->p[3];
	map->h=in->p[0]->p[4];
    
	dmat *ampg=dref(in->p[1]);
	map->p=ampg->p;
	map->nx=ampg->nx;
	map->ny=ampg->ny;
	dcellfree(in);
	dfree_keepdata(ampg);
	if(map->ox/map->dx*2+map->nx > 2
	   ||map->oy/map->dx*2+map->ny > 2){
	    warning("Ampground %s is not centered.\n",fn);
	}
    }else{
	error("Invalid format. magic=%u\n", magic);
    }
    zfclose(fp);
    return map;
}

/**
 * Read map_t arr from file
 */

map_t **maparrread(int *nlayer, const char *format, ...){
    format2fn;
    dcell *in=dcellread("%s", fn);
    map_t **map=dcell2map(nlayer, in);
    dcellfree(in);
    return map;
}


/**
   convert a dmat to map_t.
*/
rectmap_t* d2rectmap(dmat *in){
    rectmap_t *map=realloc(dref(in), sizeof(rectmap_t));
    char *header=in->header;
    if(!in->header){
	error("this dmat has no header\n");
    }
    map->ox=search_header_num(header,"ox");
    map->oy=search_header_num(header,"oy");
    map->dx=search_header_num(header,"dx");
    map->dy=search_header_num(header,"dy");
    map->txdeg=search_header_num(header,"txdeg");
    map->tydeg=search_header_num(header,"tydeg");
    map->ftel=search_header_num(header,"ftel");
    map->fexit=search_header_num(header,"fexit");
    map->fsurf=search_header_num(header,"fsurf");
    if(isnan(map->ox) || isnan(map->oy) || isnan(map->dx) || isnan(map->fsurf)){
	error("header is needed to convert dmat to map_t\n");
    }
    return map;
}

/**
 * convert a mmap'ed dcell to map_t array
 */
rectmap_t **dcell2rectmap(int *nlayer, dcell *in){
    *nlayer=in->nx*in->ny;
    rectmap_t **map=calloc(in->nx*in->ny, sizeof(rectmap_t*));
    for(long i=0; i<in->nx*in->ny; i++){
	map[i]=d2rectmap(in->p[i]);
    }
    return map;
}
/**
 * Readrtectmap_t from file
 */
rectmap_t *rectmapread(const char *format, ...){
    format2fn;
    dmat *in=dread("%s", fn);
    rectmap_t *map=d2rectmap(in);
    dfree(in);
    return map;
}

/**
 * Read rectmap_t arr from file
 */
rectmap_t **rectmaparrread(int *nlayer, const char *format, ...){
    format2fn;
    dcell *in=dcellread("%s", fn);
    rectmap_t **map=dcell2rectmap(nlayer, in);
    dcellfree(in);
    return map;
}

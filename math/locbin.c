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

#include "mathdef.h"
#include "locbin.h"

/**
   Verify the magic, dimension and read in the loc_t by calling locreaddata2().
 */
loc_t *locreaddata(file_t *fp, header_t *header){
    const double tol=1e-7;
    header_t header2;
    if(!header){
	header=&header2;
	read_header(header, fp);
    }
    if((header->magic&M_DBL)!=M_DBL){
	error("magic=%u. Expect %x\n", header->magic, M_DBL);
    }
    double dx=fabs(search_header_num(header->str,"dx"));
    double dy=fabs(search_header_num(header->str,"dy"));
    free(header->str);header->str=0;
    long nx=header->nx;
    long ny=header->ny;
    loc_t *out;
    if(nx==0 || ny==0){
	out=NULL;
    }else{
	out=locnew(nx, 0, 0);
	zfread(out->locx, sizeof(double), nx, fp);
	zfread(out->locy, sizeof(double), nx, fp);
	double dxd=INFINITY, dyd=INFINITY;
	for(long i=0; i<out->nloc-1; i++){
	    double dxi=fabs(out->locx[i+1]-out->locx[i]);
	    double dyi=fabs(out->locy[i+1]-out->locy[i]);
	    if(dxi>tol && dxi+tol<dxd){
		dxd=dxi;
	    }
	    if(dyi>tol && dyi+tol<dyd){
		dyd=dyi;
	    }
	}
	if(fabs(dx)<tol || isnan(dx)){
	    dx=dxd;//use value derived from data
	}else if(fabs(dx-dxd)>tol && isfinite(dxd)){
	    error("Specified dx=%.15g doesn't agree with data: %.15g\n", dx, dxd);
	}

	if(fabs(dy)<tol || isnan(dy)){
	    dy=dyd;
	}else if(fabs(dy-dyd)>tol && isfinite(dyd)){
	    error("Specified dy=%.15g doesn't agree with data: %.15g\n", dy, dyd);
	}

	out->dx=fabs(dx);
	out->dy=fabs(dy);
	out->map=NULL;
    }
    return out;
}

void locwritedata(file_t *fp, const loc_t *loc){
    char str[120];
    snprintf(str,120,"dx=%.15g;\ndy=%.15g;iac=%.15g\n",loc->dx,loc->dy,loc->iac);
    header_t header={M_LOC64, 0, 0, str};
    if(loc){
	header.nx=loc->nloc;
	header.ny=2;
    }
    write_header(&header,fp);
    if(loc){
	zfwrite(loc->locx, sizeof(double),loc->nloc,fp);
	zfwrite(loc->locy, sizeof(double),loc->nloc,fp);
    }
}

void mapwritedata(file_t *fp, map_t *map){
    if(map){
	if(!map->header){
	    char header[1024];
	    snprintf(header,1024,"ox=%.15g\noy=%.15g\ndx=%.15g\ndy=%.15g\nh=%.15g\nvx=%.15g\nvy=%.15g\n",
		     map->ox, map->oy, map->dx, map->dy, map->h, map->vx, map->vy);
	    map->header=strdup(header);
	}
    }
    dmat *in=(dmat*) map;
    writebindata(fp, in);
}

/**
   convert a dmat to map_t.
*/
map_t* d2map(dmat *in){
    map_t *map=myrealloc(dref(in), 1, map_t);
    memset((char*)map+sizeof(dmat), 0, sizeof(map_t)-sizeof(dmat));
    char *header=in->header;
    map->iac=0;
    map->ox=search_header_num(header,"ox");
    map->oy=search_header_num(header,"oy");
    map->dx=search_header_num(header,"dx");
    map->dy=search_header_num(header,"dy");
    map->h =search_header_num(header,"h");
    map->vx=search_header_num(header,"vx");
    map->vy=search_header_num(header,"vy");
    if(isnan(map->dx)){
	warning("dx is not specified in header, set to 1/64.\n");
	map->dx=1./64.;
    }
    if(isnan(map->dy)){
	map->dy=map->dx;
    }
    if(isnan(map->ox) || isnan(map->oy)){
	map->ox=-map->nx/2*map->dx;
	map->oy=-map->ny/2*map->dy;
    }
    if(isnan(map->h)) map->h=0;
    if(isnan(map->vx)) map->vx=0;
    if(isnan(map->vy)) map->vy=0;
    return map;
}

/**
 * convert a mmap'ed dcell to map_t array
 */
mapcell *dcell2map(dcell *in){
    mapcell *map=(mapcell*)cellnew(in->nx, in->ny);
    for(long i=0; i<in->nx*in->ny; i++){
	if(!in->p[i]->header && in->header){
	    in->p[i]->header=strdup(in->header);
	}
	map->p[i]=d2map(in->p[i]);
    }
    return map;
}

/**
 * Read map_t from file
 */
map_t *mapreaddata(file_t *fp, header_t *header){
    header_t header2={0,0,0,0};
    if(!header){
	header=&header2;
	read_header(header, fp);
    }
    map_t *map=NULL;
    if(header->magic==M_DBL){
	dmat *in=dreaddata(fp, header);
	map=d2map(in);
	dfree(in);
    }else if(iscell(&header->magic)){/*old format. */
	dcell *in=dcellreaddata(fp, header);
	if(fabs(in->p[0]->p[0]-in->p[0]->p[1])>1.e-14){
	    error("Map should be square\n");
	}
	map=mapnew(in->p[1]->nx, in->p[1]->ny, 
		   in->p[0]->p[0], in->p[0]->p[0], in->p[1]->p);
	map->ox=in->p[0]->p[2];
	map->oy=in->p[0]->p[3];
	map->h=in->p[0]->p[4];
	dfree_keepdata(in->p[1]);
	dcellfree(in);
    }else{
	error("Invalid format. magic=%u\n", header->magic);
    }
    free(header->str);header->str=0;
    return map;
}

/**
   convert a dmat to rmap_t.
*/
rmap_t* d2rmap(dmat *in){
    rmap_t *map=myrealloc(dref(in), 1, rmap_t);
    memset((char*)map+sizeof(dmat), 0, sizeof(rmap_t)-sizeof(dmat));
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
    if(isnan(map->ox) || isnan(map->oy) || isnan(map->dx)|| isnan(map->dy) || isnan(map->fsurf)){
	error("header is needed to convert dmat to map_t\n");
    }
    return map;
}

/**
 * convert a mmap'ed dcell to map_t array
 */
rmap_t **dcell2rmap(int *nlayer, dcell *in){
    *nlayer=in->nx*in->ny;
    rmap_t **map=mycalloc(in->nx*in->ny,rmap_t*);
    for(long i=0; i<in->nx*in->ny; i++){
	map[i]=d2rmap(in->p[i]);
    }
    return map;
}
/**
 * Readrtectmap_t from file
 */
rmap_t *rmapread(const char *format, ...){
    format2fn;
    dmat *in=dread("%s", fn);
    rmap_t *map=d2rmap(in);
    dfree(in);
    return map;
}

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
#include "mathdef.h"

void *cellnew(long nx, long ny){
    cell *dc;
    if(nx<0) nx=0;
    if(ny<0) ny=0;
    dc=calloc(1, sizeof(cell));
    dc->id=MCC_ANY;
    dc->nx=nx;
    dc->ny=ny;
    dc->p=calloc(nx*ny, sizeof(void*));
    return dc;
}
void free_by_id(void *pix){
    if(!pix) return;
    long id=*((long*)(pix));
    switch(id & 0xFFFF){
    case MCC_ANY:
	cellfree_do(pix);break;
    case M_DBL:
	dfree(pix);break;
    case M_CMP:
	cfree(pix);break;
    case M_FLT:
	sfree(pix);break;
    case M_ZMP:
	zfree(pix);break;
    case M_LOC64:
	locfree(pix);break;
    default:
	error("Unknown id=%lx\n", id);
    }
}
/**
   Free a cell object.
*/
void cellfree_do(void *dc_){
    cell *dc=(cell*)dc_;
    if(dc && dc->id!=MCC_ANY) error("Invalid use\n");
    if(!dc_) return;
    if(dc->p){
	for(int ix=0; ix<dc->nx*dc->ny; ix++){
	    free_by_id(dc->p[ix]);
	}
	if(dc->mmap){
	    mmap_unref(dc->mmap);
	}else{
	    free(dc->header);
	}
	free(dc->p);dc->p=0;
    }
    if(dc->m) free_by_id(dc->m);
    if(dc->fft) dfft_free_plan(dc->fft);
    free(dc);
}
void write_by_id(file_t *fp, const void *pix, long id){
    if(!id && !pix){
	error("cannot deterine id\n");
    }if(!id){
	id=*((long*)(pix));
    }else if(pix){
	assert(id==*((long*)(pix)));
    }
    switch(id){
    case MCC_ANY:
	cellwritedata(fp, pix);break;
    case M_DBL:
	dwritedata(fp, pix);break;
    case M_CMP:
	cwritedata(fp, pix);break;
    case M_FLT:
	swritedata(fp, pix);break;
    case M_ZMP:
	zwritedata(fp, pix);break;
    case M_LOC64:
	locwritedata(fp, pix);break;
    case M_MAP64:
	mapwritedata(fp, (map_t*)pix);break;
    default:
	error("Unknown id=%lx\n", id);
    }
}

/**
   Free a cell object.
*/
void cellwritedata(file_t *fp, const void *dc_){
    cell* dc=(cell*)dc_;
    if(dc && dc->id!=MCC_ANY) error("Invalid use\n");
    uint64_t nx=0;
    uint64_t ny=0;
    if(dc){
	nx=dc->nx;
	ny=dc->ny;
    }
    header_t header={MCC_ANY, nx, ny, dc?dc->header:NULL};
    write_header(&header, fp);
    if(nx && ny){
	for(int ix=0; ix<dc->nx*dc->ny; ix++){
	    write_by_id(fp, dc->p[ix], 0);
	}
    }
}
void cellwrite(const void *dc, const char* format,...){
    format2fn;
    file_t *fp=zfopen(fn,"wb");
    cellwritedata(fp,dc);
    zfclose(fp);
}
void *readdata_by_id(file_t *fp, uint32_t id, header_t *header){
    header_t header2;
    if(!header){
	header=&header2;
	read_header(header, fp);
    }
    if(iscell(header->magic)){
	return cellreaddata(fp, id, header);
    }
    if(!id) id=header->magic;
    switch(id){
    case M_DBL: return dreaddata(fp, header);break;
    case M_FLT: return sreaddata(fp, header);break;
    case M_CMP: return creaddata(fp, header);break;
    case M_ZMP: return zreaddata(fp, header);break;
    case M_LOC64: return locreaddata(fp, header); break;
    case M_MAP64: return mapreaddata(fp, header); break;
    default:error("id=%ux\n", id);
    }
    return 0;
}
void* cellreaddata(file_t *fp, uint32_t id, header_t *header){
    header_t header2;
    if(!header){
	header=&header2;
	read_header(header, fp);
    }
    long nx,ny;
    header_t *headerc=check_cell(header, &nx, &ny);//return null if is cell.
    cell *out;
    if(!headerc || !zfisfits(fp)){/*genuine cell array or mat data.*/
	out=cellnew((long)nx,(long)ny);
	if(!headerc){
	    out->header=header->str; header->str=NULL;
	}
        for(long ix=0; ix<nx*ny; ix++){
	    out->p[ix]=readdata_by_id(fp, id, headerc);
	}
    }else{/*fits format. need to read extensions*/
	if(ny!=1) error("inconsistent\n");
	void **tmp=calloc(nx, sizeof(void*));
	tmp[nx-1]=readdata_by_id(fp, id, headerc);
	while(!read_header2(headerc, fp)){
	    nx++;
	    tmp=realloc(tmp, nx*sizeof(void*));
	    tmp[nx-1]=readdata_by_id(fp, id, headerc);
	}
	out=cellnew(nx, 1);
	memcpy(out->p, tmp, sizeof(void*)*nx);
	free(tmp);
    }
    return out;
}
void* cellread(uint32_t id, const char *format, ...){
    format2fn;
    file_t *fp=zfopen(fn,"rb");
    void *out=cellreaddata(fp, id, 0);
    zfeof(fp);
    zfclose(fp);
    return out;
}

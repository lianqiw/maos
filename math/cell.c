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
    if(nx*ny>0){
	dc->p=calloc(nx*ny, sizeof(void*));
    }
    return dc;
}
void cellresize(void *in, long nx, long ny){
    cell *A=(cell*)in;
    if(A->nx==nx || A->ny==1){
	A->p=realloc(A->p, sizeof(cell*)*nx*ny);
	if(nx*ny>A->nx*A->ny){
	    memset(A->p+A->nx*A->ny, 0, (nx*ny-A->nx*A->ny)*sizeof(cell*));
	}
    }else{
	cell **p=calloc(nx*ny, sizeof(cell*));
	long minx=A->nx<nx?A->nx:nx;
	long miny=A->ny<ny?A->ny:ny;
	for(long iy=0; iy<miny; iy++){
	    memcpy(p+iy*nx, A->p+iy*A->nx, sizeof(cell*)*minx);
	}
	free(A->p);
	A->p=p;
    }
    A->nx=nx;
    A->ny=ny;
}
/**
   free a mat or cell object.
*/
void cellfree_do(void *pix){
    if(!pix) return;
    long id=*((long*)(pix));
    switch(id & 0xFFFF){
    case MCC_ANY:{
	cell *dc=(cell*)pix;
	if(dc->p){
	    for(int ix=0; ix<dc->nx*dc->ny; ix++){
		cellfree_do(dc->p[ix]);
	    }
	    memset(dc->p, 0, sizeof(void*)*dc->nx*dc->ny);
	    free(dc->p);dc->p=0;
	}
	if(dc->mmap){
	    mmap_unref(dc->mmap);
	}else{
	    free(dc->header);
	}
	if(dc->m) cellfree_do(dc->m);
	if(dc->fft) dfft_free_plan(dc->fft);
	memset(dc, 0, sizeof(cell));
	free(dc);
    }break;
    case M_DBL:
	dfree_do(pix,0);break;
    case M_CMP:
	cfree_do(pix,0);break;
    case M_FLT:
	sfree_do(pix,0);break;
    case M_ZMP:
	zfree_do(pix,0);break;
    case M_LOC64:
	locfree_do(pix);break;
    case M_SP32:
    case M_SP64:
	spfree_do(pix);break;
    case M_SSP32:
    case M_SSP64:
	sspfree_do(pix);break;
    case M_CSP32:
    case M_CSP64:
	cspfree_do(pix);break;
    case M_ZSP64:
    case M_ZSP32:
	zspfree_do(pix);break;
    default:
	error("Unknown id=%lx\n", id);
    }
}

void writedata_by_id(file_t *fp, const void *pix, long id){
    if(pix){
	if(!id){
	    id=*((long*)(pix));
	}else if (id!=MCC_ANY){
	    long id2=*((long*)(pix));
	    if((id & id2)!=id && (id & id2) != id2){
		error("id=%ld, id2=%ld, mismatch\n", id, id2);
	    }
	}
    }else if(!id){
	id=MCC_ANY;//default for empty array
    }
    switch(id){
    case MCC_ANY:{
	const cell *dc=pix;
	uint64_t nx=0;
	uint64_t ny=0;
	if(dc){
	    nx=dc->nx;
	    ny=dc->ny;
	}
	id=0;/*determine id first for empty cell*/
	if(nx && ny){
	    for(int ix=0; ix<nx*ny; ix++){
		if(dc->p[ix]){
		    id=*((long*)(dc->p[ix]));
		    if(!id){
			error("id is not set\n");
		    }
		}
	    }
	    if(!id){
		nx=0; ny=0;
	    }
	}
	header_t header={MCC_ANY, nx, ny, dc?dc->header:NULL};
	write_header(&header, fp);
	if(id){
	    for(int ix=0; ix<dc->nx*dc->ny; ix++){
		writedata_by_id(fp, dc->p[ix], 0);
	    }
	}
    }
	break;
    case M_DBL:
	dwritedata(fp, (dmat*)pix);break;
    case M_CMP:
	cwritedata(fp, (cmat*)pix);break;
    case M_FLT:
	swritedata(fp, (smat*)pix);break;
    case M_ZMP:
	zwritedata(fp, (zmat*)pix);break;
    case M_LONG:
	lwritedata(fp, (lmat*)pix);break;
    case M_LOC64:
	locwritedata(fp, (loc_t*)pix);break;
    case M_MAP64:
	mapwritedata(fp, (map_t*)pix);break;
    case M_SP32:
    case M_SP64:
	spwritedata(fp, (dsp*)pix);break;
    case M_SSP32:
    case M_SSP64:
	sspwritedata(fp, (ssp*)pix);break;
    case M_CSP32:
    case M_CSP64:
	cspwritedata(fp, (csp*)pix);break;
    case M_ZSP64:
    case M_ZSP32:
	zspwritedata(fp, (zsp*)pix);break;	
    default:
	error("Unknown id=%lx\n", id);
    }
}

void write_by_id(const void *dc, long id, const char* format,...){
    format2fn;
    file_t *fp=zfopen(fn,"wb");
    writedata_by_id(fp, dc, id);
    zfclose(fp);
}
void *readdata_by_id(file_t *fp, long id, int level, header_t *header){
    header_t header2={0};
    if(!header){
	header=&header2;
	read_header(header, fp);
    }
    void *out=0;
    if(zfisfits(fp) || level==0){
	switch(level){
	case 0:/*read a mat*/
	    if(!id) id=header->magic;
	    switch(id){
	    case M_DBL: out=dreaddata(fp, header);break;
	    case M_FLT: out=sreaddata(fp, header);break;
	    case M_CMP: out=creaddata(fp, header);break;
	    case M_ZMP: out=zreaddata(fp, header);break;
	    case M_LONG: out=lreaddata(fp, header);break;
	    case M_LOC64: out=locreaddata(fp, header); break;
	    case M_MAP64: out=mapreaddata(fp, header); break;
	    case M_SP32:
	    case M_SP64:
		out=spreaddata(fp, header);break;
	    case M_SSP32:
	    case M_SSP64:
		out=sspreaddata(fp, header);break;
	    case M_CSP32:
	    case M_CSP64:
		out=cspreaddata(fp, header);break;
	    case M_ZSP64:
	    case M_ZSP32:
		out=zspreaddata(fp, header);break;	
	    default:error("id=%lx\n", id);
	    }
	    break;
	case 1:{/*read a cell from fits*/
	    int maxlen=10;
	    void **tmp=malloc(maxlen*sizeof(void*));
	    int nx=0;
	    do{
		if(nx>=maxlen){
		    maxlen*=2;
		    tmp=realloc(tmp, sizeof(void*)*maxlen);
		}
		tmp[nx++]=readdata_by_id(fp, id, 0, header);
		free(header->str);header->str=0;
	    }while(!read_header2(header, fp));
	    cell *dcout=cellnew(nx, 1);
	    memcpy(dcout->p, tmp, sizeof(void*)*nx);
	    free(tmp);
	    out=dcout;
	}
	    break;
	default:
	    error("Only support zero or one level of cell when reading fits file\n");
	}
    }else{
	if(!iscell(header->magic)){
	    error("Traying to read cell from non cell data\n");
	}else{
	    long nx=header->nx;
	    long ny=header->ny;
	    cell *dcout=cellnew(nx, ny);
	    for(long i=0; i<nx*ny; i++){
		dcout->p[i]=readdata_by_id(fp, id, level-1, 0);
	    }
	    out=dcout;
	}
    }
    free(header->str);header->str=0;
    return out;
}

void* read_by_id(long id, int level, const char *format, ...){
    format2fn;
    file_t *fp=zfopen(fn,"rb");
    void *out=readdata_by_id(fp, id, level, 0);
    zfeof(fp);
    zfclose(fp);
    return out;
}

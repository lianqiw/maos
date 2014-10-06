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
/**
   \file cell.c
*/

/**
   Create a generic cell array.
*/
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
/**
   Resize a generic cell array.
*/
void cellresize(void *in, long nx, long ny){
    cell *A=(cell*)in;
    if(A->nx==nx || A->ny==1){
	int nold=A->nx*A->ny;
	int nnew=nx*ny;
	if(nnew<nold && iscell(A->id)){
	    for(int i=nnew; i<nold; i++){
		cellfree_do(A->p[i]);
	    }
	}
	A->p=realloc(A->p, sizeof(cell*)*nnew);
	if(nnew>nold){
	    memset(A->p+nold, 0, (nnew-nold)*sizeof(cell*));
	}
    }else{
	cell **p=calloc(nx*ny, sizeof(cell*));
	long minx=A->nx<nx?A->nx:nx;
	long miny=A->ny<ny?A->ny:ny;
	for(long iy=0; iy<miny; iy++){
	    for(long ix=0; ix<minx; ix++){
		p[ix+iy*nx]=A->p[ix+iy*A->nx];
		A->p[ix+iy*A->nx]=0;
	    }
	}
	if(iscell(A->id)){
	    for(long i=0; i<A->nx*A->ny; i++){
		cellfree_do(A->p[i]);
	    }
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
void cellfree_do(void *A){
    if(!A) return;
    long id=*((long*)(A));
    switch(id & 0xFFFF){
    case MCC_ANY:{
	cell *dc=(cell*)A;
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
	dfree_do(A,0);break;
    case M_CMP:
	cfree_do(A,0);break;
    case M_FLT:
	sfree_do(A,0);break;
    case M_ZMP:
	zfree_do(A,0);break;
    case M_LOC64:
	locfree_do(A);break;
    case M_DSP:
	dspfree_do(A);break;
    case M_SSP:
	sspfree_do(A);break;
    case M_CSP:
	cspfree_do(A);break;
    case M_ZSP:
	zspfree_do(A);break;
    default:
	error("Unknown id=%lx\n", id);
    }
}

void writedata_by_id(file_t *fp, const void *A, long id){
    if(A){
	if(!id){
	    id=*((long*)(A));
	}else if (id!=MCC_ANY){
	    long id2=*((long*)(A));
	    if((id & id2)!=id && (id & id2) != id2){
		error("id=%ld, id2=%ld, mismatch\n", id, id2);
	    }
	}
    }else if(!id){
	id=MCC_ANY;//default for empty array
    }
    switch(id){
    case MCC_ANY:{
	const cell *dc=A;
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
	writebindata(fp, (dmat*)A);break;
    case M_CMP:
	writebindata(fp, (cmat*)A);break;
    case M_FLT:
	writebindata(fp, (smat*)A);break;
    case M_ZMP:
	writebindata(fp, (zmat*)A);break;
    case M_LONG:
	lwritedata(fp, (lmat*)A);break;
    case M_LOC64:
	locwritedata(fp, (loc_t*)A);break;
    case M_MAP64:
	mapwritedata(fp, (map_t*)A);break;
    case M_DSP:
	dspwritedata(fp, (dsp*)A);break;
    case M_SSP:
	sspwritedata(fp, (ssp*)A);break;
    case M_CSP:
	cspwritedata(fp, (csp*)A);break;
    case M_ZSP:
	zspwritedata(fp, (zsp*)A);break;	
    default:
	error("Unknown id=%lx\n", id);
    }
}

void write_by_id(const void *A, long id, const char* format,...){
    format2fn;
    file_t *fp=zfopen(fn,"wb");
    writedata_by_id(fp, A, id);
    zfclose(fp);
}
void *readdata_by_id(file_t *fp, long id, int level, header_t *header){
    header_t header2={0};
    if(!header){
	header=&header2;
	read_header(header, fp);
    }
    void *out=0;
    if(level<0 && !iscell(header->magic)){
	level=0;
    }
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
	    case M_DSP32: case M_DSP64: /**Possible to read mismatched integer*/
		out=dspreaddata(fp, header);break;
	    case M_SSP32: case M_SSP64:
		out=sspreaddata(fp, header);break;
	    case M_CSP32: case M_CSP64:
		out=cspreaddata(fp, header);break;
	    case M_ZSP64: case M_ZSP32:
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
	    error("Trying to read cell from non cell data\n");
	}else{
	    long nx=header->nx;
	    long ny=header->ny;
	    cell *dcout=cellnew(nx, ny);
	    dcout->header=header->str; header->str=0;
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

void* readbin(const char *format, ...){
    format2fn;
    return read_by_id(0, -1, "%s", fn);
}

void writebin(const void *A, const char *format, ...){
    format2fn;
    write_by_id(A, 0, "%s", fn);
}

void writebindata(file_t *fp, const void *A){
    writedata_by_id(fp, A, 0);
}

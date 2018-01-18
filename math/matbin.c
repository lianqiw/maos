/*
  Copyright 2009-2018 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "type.h"
#include "mathdef.h"
#include "defs.h"
/**
   basic routines
*/

/**
   Function to write dense matrix data into a file pointer. Generally used by
   library developer */
void X(writedata)(file_t *fp, const X(mat) *A){
    uint64_t nx=0, ny=0;
    if(A){
	nx=(uint64_t)A->nx;
	ny=(uint64_t)A->ny;
    }
    writearr(fp, 0, sizeof(T), M_T, A?A->header:NULL, A?A->p:NULL, nx, ny);
}

/**
   Function to read dense matrix into memory from file pointer. Generally used
   by library developer.  */
X(mat) *X(readdata)(file_t *fp, header_t *header){
    header_t header2={0,0,0,0};
    if(!header){
	header=&header2;
	read_header(header, fp);
    }
    uint64_t nx, ny;
    nx=header->nx;
    ny=header->ny;
    if(!nx || !ny) {
	free(header->str);
	return NULL;
    }
    X(mat) *out;
    out=X(new)((long)nx,(long)ny);
    out->header=header->str; header->str=NULL;
    if(nx>0 && ny>0){
	if(header->magic==M_T){
	    zfread(out->p,sizeof(T),nx*ny,fp);
	}else if(M_T==M_DBL && header->magic==M_FLT){
	    float *p=mymalloc(nx*ny,float);
	    zfread(p, sizeof(float), nx*ny, fp);
	    for(uint64_t i=0; i<nx*ny; i++){
		out->p[i]=(T)p[i];
	    }
	    free(p);
#ifdef USE_COMPLEX
	}else if(M_T==M_CMP && header->magic==M_ZMP){
	    fcomplex *p=mymalloc(nx*ny,fcomplex);
	    zfread(p, sizeof(fcomplex), nx*ny, fp);
	    for(uint64_t i=0; i<nx*ny; i++){
		out->p[i]=(T)p[i];
	    }
	    free(p);
#endif
	}else{
	    error("%s is not a X(mat) file. magic=%x. We want %x\n", 
		  zfname(fp), header->magic, M_T);
	}
    }
    return out;
}
#include "matmmap.c"

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
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <fcntl.h>
#include <unistd.h>
#include "type.h"
#include "mathdef.h"
#include "defs.h"
/**
   Contains routines to write/read dense/sparse matrix into/from file.
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
    do_write(fp, 0, sizeof(T), M_T, A?A->header:NULL, A?A->p:NULL, nx, ny);
}

/**
   Function to read dense matrix into memory from file pointer. Generally used
   by library developer.  */
X(mat) *X(readdata)(file_t *fp, header_t *header){
    header_t header2;
    if(!header){
	header=&header2;
	read_header(header, fp);
    }
    uint64_t nx, ny;
    nx=header->nx;
    ny=header->ny;
    if(!nx || !ny) return NULL;
    X(mat) *out;
    out=X(new)((long)nx,(long)ny);
    out->header=header->str; header->str=NULL;
    if(nx>0 && ny>0){
	if(header->magic==M_T){
	    zfread(out->p,sizeof(T),nx*ny,fp);
	}else if(M_T==M_DBL && header->magic==M_FLT){
	    float *p=malloc(nx*ny*sizeof(float));
	    zfread(p, sizeof(float), nx*ny, fp);
	    for(int i=0; i<nx*ny; i++){
		out->p[i]=(T)p[i];
	    }
	    free(p);
	}else if(M_T==M_CMP && header->magic==M_ZMP){
	    fcomplex *p=malloc(nx*ny*sizeof(fcomplex));
	    zfread(p, sizeof(fcomplex), nx*ny, fp);
	    for(int i=0; i<nx*ny; i++){
		out->p[i]=(T)p[i];
	    }
	    free(p);
	}else{
	    error("%s is not a X(mat) file. magic=%x. We want %x\n", 
		  zfname(fp), header->magic, M_T);
	}
    }
    return out;
}
/**
   Scale a dcell array and save to file.
*/
void X(cellswrite)(X(cell) *A, double scale, const char *format, ...){
    format2fn;
    X(cell) *tmp=NULL;
    if(scale<1.e-14){
	error("scale=%g\n",scale);
    }
    X(celladd)(&tmp, 0, A, scale);
    X(cellwrite)(tmp,"%s",fn);
    X(cellfree)(tmp);
}

/**
   Scale a dcell array and save to file.
*/
void X(swrite)(X(mat) *A, double scale, const char *format, ...){
    format2fn;
    X(mat) *tmp=NULL;
    if(scale<1.e-14){
	error("scale=%g\n",scale);
    }
    X(add)(&tmp, 0, A, scale);
    X(write)(tmp,"%s",fn);
    X(free)(tmp);
}

#include "matmmap.c"

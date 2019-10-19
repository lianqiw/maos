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
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>

#include <fcntl.h>

#include "type.h"
#include "spbin.h"
#include "defs.h"
/*
   Function to write sparse matrix data into file pointed using a file
   pointer. Generally used by library developer.  We do not convert data during
   saving, but rather do the conversion during reading.
*/
void X(spwritedata)(file_t *fp, const X(sp) *sp){
    header_t header={M_SPT, 0, 0, NULL};
    if(sp && sp->nzmax){
	header.nx=sp->nx;
	header.ny=sp->ny;
	header.str=sp->header;
    }
    write_header(&header, fp);
    if(sp && sp->nzmax){
	X(spsort)((X(sp)*)sp);/*sort the matrix to have the right order */
	uint64_t nzmax;
	nzmax=sp->p[sp->ny];/*don't use sp->nzmax, which maybe larger than actual */
	zfwrite(&nzmax, sizeof(uint64_t), 1, fp);
	zfwrite(sp->p, sizeof(spint), sp->ny+1, fp);
	zfwrite(sp->i, sizeof(spint), nzmax, fp);
	zfwrite(sp->x ,sizeof(T),nzmax,fp);  
    }
}


/**
   Function to read sparse matrix data from file pointer into memory. Used by
   library developer.
  */
X(sp) *X(spreaddata)(file_t *fp, header_t *header){
    header_t header2={0,0,0,0};
    if(!header){
	header=&header2;
	read_header(header, fp);
    }
    long m=header->nx;
    long n=header->ny;
    uint64_t nzmax;
    X(sp) *out=NULL;
    if(m!=0 && n!=0){
	uint32_t magic2=0;
	uint32_t magic1=0;
	switch(header->magic){
	case M_DSP64:
	    magic1=M_DBL;
	    magic2=M_INT64;
	    break;
	case M_SSP64:
	    magic1=M_FLT;
	    magic2=M_INT64;
	    break;
	case M_DSP32:
	    magic1=M_DBL;
	    magic2=M_INT32;
	    break;
	case M_SSP32:
	    magic1=M_FLT;
	    magic2=M_INT32;
	    break;
	default:
	    error("This is not a valid sparse matrix file. magic=%x\n", header->magic);
	}
	zfread(&nzmax,sizeof(uint64_t),1,fp);
	if(nzmax!=0){
	    out=X(spnew)(m,n,nzmax);
	    out->header=header->str; header->str=NULL;
	    //readspintdata(fp, magic2, out->p, n+1);
	    //readspintdata(fp, magic2, out->i, nzmax);
	    readvec(out->p, M_SPINT, magic2, sizeof(spint), n+1, fp);
	    readvec(out->i, M_SPINT, magic2, sizeof(spint), nzmax, fp);
	    readvec(out->x, M_T, magic1, sizeof(T), nzmax, fp);
	}
    }
    free(header->str);
    return out;
}

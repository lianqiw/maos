/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "defs.h"
/**
   Routines for input/output to bin/fits file.
*/

/**
   Function to write dense matrix data into a file pointer. Generally used by
   library developer */
void X(writedata)(file_t* fp, const X(mat)* A, long ncol){
	if(ncol==-1){//initialize async data
		if(A){
			if(!A->async){
				((X(mat)*)A)->async=async_init(fp, sizeof(T), M_T, A->keywords, P(A), A->nx, A->ny);
			}else{
				dbg("%s: async is already initialized or A=%p is empty\n", zfname(fp), A);
			}
		}
	}else if(A && ncol>0 && A->async){//write async data
		long nx=A->nx;
		long ny=A->ny;
		if(ny==1){//for row vectors, allow saving by row instead.
			ny=nx;
			nx=1;
		}
		if(nx && ny){
			if(ncol>ny){
				warning("ncol=%ld>ny=%ld\n", ncol, ny);
				ncol=ny;
			}
			async_write(A->async, nx*ncol*sizeof(T), 0);
		}
	}else if(fp){//normal writing
		writearr(fp, 0, sizeof(T), M_T, A?A->keywords:NULL, A?P(A):NULL, A?A->nx:0, A?A->ny:0);
	}else{
		dbg("writedata called with invalid fp(%p) or async (ncol=%ld) information, canceled.\n", fp, ncol);
		print_backtrace();
	}
}

/**
   Function to read dense matrix into memory from file pointer. Generally used
   by library developer.  */
X(mat)* X(readdata)(file_t* fp, header_t* header){
	header_t header2={0,0,0,0};
	if(!header){
		header=&header2;
	}
	if(header->magic==0){
		read_header(header, fp);
	}
	uint64_t nx, ny;
	nx=header->nx;
	ny=header->ny;

	X(mat)* out=X(new)((long)nx, (long)ny);
	out->keywords=header->str; header->str=NULL;
	if(nx&&ny){
		readvec(P(out), M_T, header->magic, sizeof(T), nx*ny, fp);
	}
	header->magic=0; header->nx=0; header->ny=0;//prevent reuse.
	return out;
}


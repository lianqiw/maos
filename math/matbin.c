/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>

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
void X(writedata)(file_t* fp, const X(mat)* A){
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
X(mat)* X(readdata)(file_t* fp, header_t* header){
	header_t header2={0,0,0,0};
	if(!header){
		header=&header2;
		read_header(header, fp);
	}
	uint64_t nx, ny;
	nx=header->nx;
	ny=header->ny;

	X(mat)* out=X(new)((long)nx, (long)ny);
	out->header=header->str; header->str=NULL;
	readvec(out->p, M_T, header->magic, sizeof(T), nx*ny, fp);
	return out;
}


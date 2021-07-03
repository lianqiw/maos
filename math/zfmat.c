/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "zfmat.h"
#include "mathdef.h"
/**
  Routines for streaming dmat data to file.
*/
struct zfmat{
	file_t* fp;     /**<pointer to file*/
	long cur;       /**<Current element*/
	long tot;       /**<Total number of elements*/
	//long curx;      /**<Current row (inner dimension)*/
	//long cury;	   /**<Current column*/
	//long nx;       /**<Total number of rows*/
	//long ny;       /**<Total number of columns*/
	uint32_t id;
};

/**
   Initializing an zfmat object that contains dmat
 */
zfmat* zfmat_init(long nx, long ny, const char* format, ...){
	format2fn;
	if(!fn) return NULL;
	if(nx<0) nx=0;
	if(ny<0) ny=0;
	zfmat* out=mycalloc(1, zfmat);
	out->fp=zfopen(fn, "wb");
	out->tot=nx*ny;
	out->id=M_REAL;
	header_t header={out->id, (uint64_t)nx, (uint64_t)ny, NULL};
	write_header(&header, out->fp);
	return out;
}

/**
   Append a A of type type into the zfmat ca, at location i.
*/
void zfmat_push(zfmat* ca, long count, real *p){
	if(!ca){
		warning_once("zfmat is NULL\n");
		return;
	}
	zfwrite(p, count*sizeof(p), 1, ca->fp);
	ca->cur+=count;
}


/**
   Close the zfmat.
*/
void zfmat_close(zfmat* ca){
	if(!ca) return;
	if(ca->tot && ca->cur>ca->tot){
		warning("zfmat %s is initialized with %ld elements, "
			"but %ld elements were written\n",
			zfname(ca->fp), ca->tot, ca->cur);
	} else if(ca->cur<ca->tot){
		zfseek(ca->fp, (ca->tot-ca->cur-1)*sizeof(real), SEEK_SET);
		real zero=0;
		zfmat_push(ca, 1, &zero);
	}
	zfclose(ca->fp);
	free(ca);
}
/**
   Close an array of zfmat
*/
void zfmat_close_n(zfmat** ca, int nc){
	if(!ca) return;
	for(int ic=0; ic<nc; ic++){
		zfmat_close(ca[ic]);
	}
	free(ca);
}

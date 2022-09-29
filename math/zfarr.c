/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#include "zfarr.h"
#include "mathdef.h"
/**
  Routines for streaming data to file.
*/
struct zfarr{
	file_t* fp;     /**<pointer to file*/
	long cur;       /**<Current element*/
	long tot;       /**<Total number of elements*/
	uint32_t id;
};

/**
   Initializing an zfarray object that contains arrays of dmat, cmat, dcell or ccell
 */
zfarr* zfarr_init(long nx, long ny, const char* format, ...){
	format2fn;
	if(!fn) return NULL;
	if(nx<0) nx=0;
	if(ny<0) ny=0;
	zfarr* ca=mycalloc(1, zfarr);
	ca->fp=zfopen(fn, "wb");
	ca->cur=0;
	if(!zfisfits(ca->fp)){
		ca->tot=nx*ny;
	}
	header_t header={MCC_ANY, (uint64_t)nx, (uint64_t)ny, NULL};
	write_header(&header, ca->fp);
	return ca;
}
/**
   Append a A of type type into the zfarr ca, at location i.
*/
void zfarr_push_cell(zfarr* ca, int i, const cell* A){
	if(!ca){
		warning_once("zfarr is NULL\n");
		return;
	}
	if(i>=0&&ca->cur>i){
		warning("Invalid. %s, cur=%ld, i=%x, skip.\n", zfname(ca->fp), ca->cur, i);
		print_backtrace();
		return;
	}
	/*
	if(ca->tot && i>=ca->tot){
		warning_once("zfarr %s overflow. Size is %ld, current position is %ld\n", zfname(ca->fp), ca->tot, ca->cur);
		print_backtrace();
		return;
	}*/
	M_ID id=0;
	if(A){
		id=A->id;
		if(!ca->id){
			ca->id=id;
		}else if(ca->id!=id){
			warning("Data mismatch: %s, existing data is %x, new data is %x", zfname(ca->fp), ca->id, id);
		}
	}
	while(ca->cur<i && ca->tot){//fill blank
		writedata_by_id(ca->fp, 0, ca->id, 0);
		ca->cur++;
	}
	if(A || ca->tot){
		writedata_by_id(ca->fp, A, ca->id, 0);
		ca->cur++;
	}
}

/**
   Close the zfarr.
*/
void zfarr_close(zfarr* ca){
	if(!ca) return;
	if(ca->tot){//total is specified
		if(ca->cur>ca->tot){//overflow
			warning("zfarr %s is initialized with %ld elements, "
				"but %ld elements were written\n",
				zfname(ca->fp), ca->tot, ca->cur);
		} else if(ca->cur<ca->tot){//under fill
			zfarr_push_cell(ca, ca->tot-1, NULL);
		}
	}else if(!zfisfits(ca->fp)){//total is not specified for bin file
		zfrewind(ca->fp);
		header_t header={MCC_ANY, (uint64_t)ca->cur, (uint64_t)1, NULL};
		write_header(&header, ca->fp);
	}
	zfclose(ca->fp);
	free(ca);
}
/**
   Close an array of zfarr
*/
void zfarr_close_n(zfarr** ca, int nc){
	if(!ca) return;
	for(int ic=0; ic<nc; ic++){
		zfarr_close(ca[ic]);
	}
	free(ca);
}

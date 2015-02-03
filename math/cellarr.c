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
#include <sys/file.h>
#include "cellarr.h"
#include "mathdef.h"

/**
   Initializing an cellarray object that contains arrays of dmat, cmat, dcell or ccell
 */
cellarr* cellarr_init(long nx, long ny,const char*format,...){
    format2fn;
    if(nx<0) nx=0;
    if(ny<0) ny=0;
    cellarr *out=calloc(1, sizeof(cellarr));
    out->fp=zfopen(fn,"wb");
    out->cur=0;
    out->tot=nx*ny;
    header_t header={MCC_ANY, nx, ny, NULL};
    write_header(&header, out->fp);
    return out;
}

/**
   Append a A of type type into the cellarr ca, at location i.
*/
void cellarr_push(cellarr *ca, int i, const void *p){
    if(!p) return;//nothing to be done
    if(!ca) error("cellarr is NULL\n");
    if(i>=0 && ca->cur>i) warning("Invalid. cur=%ld, i=%d\n", ca->cur, i);
    long id=*((long*)(p));
    while(ca->cur<i) {
	writedata_by_id(ca->fp, 0, id); 
	ca->cur++;
    }
    writedata_by_id(ca->fp, p, 0); 
    ca->cur++;
}

#define cellarr_cell(type)					\
    void cellarr_##type##cell(cellarr *ca, int i, const type##cell *A){ \
	cellarr_push(ca, i, A);}					\

#define cellarr_mat(type)					\
    void cellarr_##type##mat(cellarr *ca, int i, const type##mat *A){	\
	cellarr_push(ca, i, A);}

cellarr_cell(d);
cellarr_cell(s);
cellarr_cell(c);
cellarr_cell(z);
cellarr_mat(d);
cellarr_mat(s);
cellarr_mat(c);
cellarr_mat(z);

/**
   Close the cellarr.
*/
void cellarr_close(cellarr *ca){
    if(!ca) return;
    if(ca->cur !=ca->tot){
	warning2("cellarr %s is initialized with %ld elements, "
		 "but %ld elements are written\n",
		 zfname(ca->fp),ca->tot,ca->cur);
    }
    zfclose(ca->fp);
    free(ca);
}
/**
   Close an array of cellarr
*/
void cellarr_close_n(cellarr **ca, int nc){
    if(!ca) return;
    for(int ic=0; ic<nc; ic++){
	cellarr_close(ca[ic]);
    }
    free(ca);
}

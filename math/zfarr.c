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
#include <sys/file.h>
#include "zfarr.h"
#include "mathdef.h"

/**
   Initializing an zfarray object that contains arrays of dmat, cmat, dcell or ccell
 */
zfarr* zfarr_init(long nx, long ny,const char*format,...){
    format2fn;
    if(nx<0) nx=0;
    if(ny<0) ny=0;
    zfarr *out=mycalloc(1,zfarr);
    out->fp=zfopen(fn,"wb");
    out->cur=0;
    out->tot=nx*ny;
    header_t header={MCC_ANY, (uint64_t)nx, (uint64_t)ny, NULL};
    write_header(&header, out->fp);
    return out;
}

/**
   Append a A of type type into the zfarr ca, at location i.
*/
void zfarr_push(zfarr *ca, int i, const void *p){
    if(!p){
	warning_once("p is empty\n");
	return;//nothing to be done
    }
    if(!ca){
	warning_once("zfarr is NULL\n");
	return;
    }
    if(i>=0 && ca->cur>i) {
	warning("Invalid. cur=%ld, i=%d, skip.\n", ca->cur, i);
	return;
    }
    uint32_t id=((cell*)p)->id;
    while(ca->cur<i && !zfisfits(ca->fp)) {
	writedata_by_id(ca->fp, 0, id); 
	ca->cur++;
    }
    writedata_by_id(ca->fp, p, id); 
    ca->cur++;
}

#define zfarr_cell(type)					\
    void zfarr_##type##cell(zfarr *ca, int i, const type##cell *A){ \
	zfarr_push(ca, i, A);}					\

#define zfarr_mat(type)					\
    void zfarr_##type##mat(zfarr *ca, int i, const type##mat *A){	\
	zfarr_push(ca, i, A);}

zfarr_cell(d);
zfarr_cell(s);
zfarr_cell(c);
zfarr_cell(z);
zfarr_mat(d);
zfarr_mat(s);
zfarr_mat(c);
zfarr_mat(z);

/**
   Close the zfarr.
*/
void zfarr_close(zfarr *ca){
    if(!ca) return;
    /*if(ca->cur !=ca->tot){
	warning("zfarr %s is initialized with %ld elements, "
		 "but %ld elements are written\n",
		 zfname(ca->fp),ca->tot,ca->cur);
		 }*/
    zfclose(ca->fp);
    free(ca);
}
/**
   Close an array of zfarr
*/
void zfarr_close_n(zfarr **ca, int nc){
    if(!ca) return;
    for(int ic=0; ic<nc; ic++){
	zfarr_close(ca[ic]);
    }
    free(ca);
}

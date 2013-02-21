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
#include "dmat.h"
#include "smat.h"
#include "cmat.h"
#include "zmat.h"
/**
   \file cellarr.c
   
   cellarr is an object used to write arrays of dcell or ccell into file.
   Mainly used to output PSF into files.
 */

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
   Append a dcell A into the cellarr ca.
 */
void cellarr_dcell(cellarr *ca, const dcell *A){
    if(!ca) error("callarr is NULL\n");
    dcellwritedata(ca->fp,A);
    ca->cur++;
    zflush(ca->fp);
}

/**
   Append a dcell A into the cellarr ca.
 */
void cellarr_scell(cellarr *ca, const scell *A){
    if(!ca) error("callarr is NULL\n");
    scellwritedata(ca->fp,A);
    ca->cur++;
    zflush(ca->fp);
}

/**
   Append a ccell A into the cellarr ca.
 */
void cellarr_ccell(cellarr *ca, const ccell *A){
    if(!ca) error("callarr is NULL\n");
    ccellwritedata(ca->fp,A);
    ca->cur++;
    zflush(ca->fp);
}
/**
   Append a ccell A into the cellarr ca.
 */
void cellarr_zcell(cellarr *ca, const zcell *A){
    if(!ca) error("callarr is NULL\n");
    zcellwritedata(ca->fp,A);
    ca->cur++;
    zflush(ca->fp);
}

/**
   Append a dmat A into the cellarr ca.
 */
void cellarr_dmat(cellarr *ca, const dmat *A){
    if(!ca) error("callarr is NULL\n");
    dwritedata(ca->fp,A);
    ca->cur++;
    zflush(ca->fp);
}

/**
   Append a dmat A into the cellarr ca.
 */
void cellarr_smat(cellarr *ca, const smat *A){
    if(!ca) error("callarr is NULL\n");
    swritedata(ca->fp,A);
    ca->cur++;
    zflush(ca->fp);
}
/**
   Append a ccell A into the cellarr ca.
 */
void cellarr_cmat(cellarr *ca, const cmat *A){
    if(!ca) error("callarr is NULL\n");
    cwritedata(ca->fp,A);
    ca->cur++;
    zflush(ca->fp);
}

/**
   Append a ccell A into the cellarr ca.
 */
void cellarr_zmat(cellarr *ca, const zmat *A){
    if(!ca) error("callarr is NULL\n");
    zwritedata(ca->fp,A);
    ca->cur++;
    zflush(ca->fp);
}
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

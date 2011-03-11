/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "cmat.h"
#include "bin.h"
/**
   \file cellarr.c
   
   cellarr is an object used to write arrays of dcell or ccell into file.
   Mainly used to output PSF into files.
 */

/**
   Initializing an cellarray object that contains arrays of dmat, cmat, dcell or ccell
 */
cellarr* cellarr_init(uint64_t nx, uint64_t ny,const char*format,...){
    format2fn;
    cellarr *out=calloc(1, sizeof(cellarr));
    out->fp=zfopen(fn,"wb");
    out->cur=0;
    out->tot=nx*ny;
    write_magic(MCC_ANY, out->fp);
    zfwrite(&nx,sizeof(uint64_t),1,out->fp);
    zfwrite(&ny,sizeof(uint64_t),1,out->fp);
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
   Append a ccell A into the cellarr ca.
 */
void cellarr_ccell(cellarr *ca, const ccell *A){
    if(!ca) error("callarr is NULL\n");
    ccellwritedata(ca->fp,A);
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
   Append a ccell A into the cellarr ca.
 */
void cellarr_cmat(cellarr *ca, const cmat *A){
    if(!ca) error("callarr is NULL\n");
    cwritedata(ca->fp,A);
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

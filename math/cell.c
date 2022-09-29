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

#include <unistd.h>
#include "mathdef.h"

/**
   Create a generic cell array.
*/
cell* cellnew(long nx, long ny){
	cell* dc;
	if(nx<0) nx=0;
	if(ny<0) ny=0;
	dc=mycalloc(1, cell);
	dc->id=MCC_ANY;
	dc->nx=nx;
	dc->ny=ny;
	if(nx*ny>0){
		P(dc)=mycalloc(nx*ny, cell*);
	}
	return dc;
}

/**
   check the size of cell array if exist. Otherwise create it
*/
void cellinit(cell** A, long nx, long ny){
	if(*A){
		if(!(iscell(*A)&&(*A)->nx==nx&&(*A)->ny==ny)){
			error("type or size mismatch: id=%u, require (%ld, %ld), A is (%ld, %ld)\n",
				(*A)->id, nx, ny, (*A)->nx, (*A)->ny);
		}
	} else{
		*A=cellnew(nx, ny);
	}
}
/**
   check the size of cell array if exist. Otherwise create it. m is created for dcell
*/
void cellinit2(cell** A, const cell* B){
	if(!iscell(B)){
		error("Invalid usage. B must be cell\n");
	}if(*A){
		cellinit(A, B->nx, B->ny);
	} else{
		M_ID magic=0;
		for(int ib=0; ib<B->nx*B->ny; ib++){
			if(P(B, ib)){
				if(!magic){
					magic=P(B, ib)->id;
				} else{
					if(magic!=P(B, ib)->id){
						error("Mixed type is not supported: %u vs %u\n",
							magic, P(B, ib)->id);
					}
				}
			}
		}
		if(magic==M_REAL){
			*A=(cell*)dcellnew2((const dcell*)B);
		} else if(magic==M_COMP){
			*A=(cell*)ccellnew2((const ccell*)B);
		} else{
			*A=(cell*)cellnew(B->nx, B->ny);
		}
	}
}

/**
   Obtain the dimensions.
*/
void celldim(const cell* A, long* nx, long* ny, long** nxs, long** nys){
	*nx=0;
	*ny=0;
	if(!A || !A->nx || !A->ny){
		*nxs=*nys=NULL;
		return;
	} 
	*nxs=mycalloc(A->nx, long);
	*nys=mycalloc(A->ny, long);
	
	for(long ix=0; ix<A->nx; ix++){
		for(long iy=0; iy<A->ny; iy++){
			if(!isempty(P(A, ix, iy))){
				if(!(*nxs)[ix]){
					*nx+=P(A, ix, iy)->nx;
					(*nxs)[ix]=P(A, ix, iy)->nx;
				} else if(P(A, ix, iy)->nx != (*nxs)[ix]){
					error("Cells along a row has different number of rows\n");
				}
			}
		}
	}
	for(long iy=0; iy<A->ny; iy++){
		for(long ix=0; ix<A->nx; ix++){
			if(!isempty(P(A, ix, iy))){
				if(!(*nys)[iy]){
					*ny+=P(A, ix, iy)->ny;
					(*nys)[iy]=P(A, ix, iy)->ny;
				} else if((*nys)[iy]!=P(A, ix, iy)->ny){
					error("Cells along a column has different number of columns\n");
				}
			}
		}
	}
}
/**
   Resize a generic cell array.
*/
void cellresize_do(cell* A, long nx, long ny){
	if(!A || A->m){ 
		warning("Can only resize a valid cell array with m unset.\n");
		return;
	}
	if(A->nx==nx||A->ny==1){
		int nold=A->nx*A->ny;
		int nnew=nx*ny;
		if(nnew<nold&&iscell(A)){
			for(int i=nnew; i<nold; i++){
				cellfree_do(P(A,i));
			}
		}
		P(A)=myrealloc(P(A), nnew, cell*);
		if(nnew>nold){
			memset(&P(A, nold), 0, (nnew-nold)*sizeof(cell*));
		}
	} else{
		cell** p=mycalloc(nx*ny, cell*);
		long minx=A->nx<nx?A->nx:nx;
		long miny=A->ny<ny?A->ny:ny;
		for(long iy=0; iy<miny; iy++){
			for(long ix=0; ix<minx; ix++){
				p[ix+iy*nx]=P(A,ix,iy);
				P(A,ix,iy)=0;
			}
		}
		if(iscell(A)){
			for(long i=0; i<A->nx*A->ny; i++){
				cellfree_do(P(A,i));
			}
		}
		free(P(A));
		P(A)=p;
	}
	A->nx=nx;
	A->ny=ny;
}


/**
   free a mat or cell object.
*/
void cellfree_do(cell* A){
	if(!A) return;
	M_ID id=A->id;
	switch(id){
	case MCC_ANY:{
		cell* dc=A;
		//if(dc->fp) writebin_async(dc, 0);//don't do this. only cells need to finish sync
		//do not close for yet for individual cells may need to do synchronization.
		if(P(dc)){
			for(int ix=0; ix<dc->nx*dc->ny; ix++){
				cellfree_do(P(dc,ix));
			}
			free(P(dc));P(dc)=0;
		}
		if(dc->header){
			free(dc->header); dc->header=NULL;		
		}
		if(dc->m){
			cellfree_do(dc->m);dc->m=NULL;		
		}
		if(dc->fft){
			dfft_free_plan(dc->fft);
			dc->fft=NULL;
		}
		if(dc->fp) zfclose(dc->fp);
		free(dc);
	}break;
	case M_DBL:
		dfree_do((dmat*)A);break;
	case M_CMP:
		cfree_do((cmat*)A);break;
	case M_FLT:
		sfree_do((smat*)A);break;
	case M_ZMP:
		zfree_do((zmat*)A);break;
	case M_LONG:
		lfree_do((lmat*)A);break;
	case M_LOC:
		locfree_do((loc_t*)A);break;
	case M_DSP:
		dspfree_do((dsp*)A);break;
	case M_SSP:
		sspfree_do((ssp*)A);break;
	case M_CSP:
		cspfree_do((csp*)A);break;
	case M_ZSP:
		zspfree_do((zsp*)A);break;
	default:
		warning("Unknown id=%x, A=%p\n", id, A);
		print_backtrace();
	}
}
/**
 * ncol: 0: normal writing. -1: initialize async data. >0: async writing.
 * 	
 * */
void writedata_by_id(file_t* fp, const cell* A, M_ID id, long ncol){
	if(fp&&ncol>0){
		error("writedata_by_id should be called with either fp or ncol, but not both, aborted.\n");
		return;
	}
	if(A){
		if(!id){
			id=A->id;
		} else if(id!=MCC_ANY){
			M_ID id2=A->id;
			if((id&id2)!=id&&(id&id2)!=id2){
				warning("write as %x, data is %x, mismatch\n", id, id2);
				id=id2;
			}
		}
	} else if(!id){
		id=MCC_ANY;//default for empty array
	}
	switch(id){
	case MCC_ANY:{//write cell
		long nx=0;
		long ny=0;
		if(A){
			nx=A->nx;
			ny=A->ny;
		}
		M_ID id2=0;//determine id first for non-empty cell
		
		if(nx&&ny){
			for(long ix=0; ix<(nx*ny); ix++){
				if(P(A,ix)){
					if(!id2) id2=P(A,ix)->id;
					else if(id2!=P(A, ix)->id){
						warning("cell (%ld) has different id (%u) than others (%u)\n", ix, P(A,ix)->id, id2);
					}
				}
			}
			if(!id2){//empty cell
				nx=0; ny=0;
			}
		}
		if(ncol<=0){
			header_t header={MCC_ANY, nx, ny, A?A->header:NULL};
			write_header(&header, fp);
		}
		
		for(long ix=0; ix<(nx*ny); ix++){
			writedata_by_id(fp, P(A,ix), id2, ncol);
		}
	}
	break;
	case M_DBL:
		dwritedata(fp, (dmat*)A, ncol);break;
	case M_CMP:
		cwritedata(fp, (cmat*)A, ncol);break;
	case M_FLT:
		swritedata(fp, (smat*)A, ncol);break;
	case M_ZMP:
		zwritedata(fp, (zmat*)A, ncol);break;
	case M_LONG:
		lwritedata(fp, (lmat*)A, ncol);break;
	case M_LOC:
		locwritedata(fp, (loc_t*)A);break;
	case M_MAP:
		map_header((map_t*)A);
		dwritedata(fp, (dmat*)A, ncol);break;
	case M_RECTMAP:
		rmap_header((rmap_t*)A);
		dwritedata(fp, (dmat*)A, ncol);break;
	case M_DSP:
		dspwritedata(fp, (dsp*)A);break;
	case M_SSP:
		sspwritedata(fp, (ssp*)A);break;
	case M_CSP:
		cspwritedata(fp, (csp*)A);break;
	case M_ZSP:
		zspwritedata(fp, (zsp*)A);break;
	default:
		error("Unknown id=%u\n", id);
	}
}

void write_by_id(const cell* A, M_ID id, const char* format, ...){
	format2fn;
	if(!fn) return;
	file_t* fp=zfopen(fn, "wb");
	if(!fp) return;
	writedata_by_id(fp, A, id, 0);
	zfclose(fp);
}

/**
 * Calls writebin with build in fp
 * */
void writecell_async(const cell* A, long ncol){
	if(ncol==0){
		warning("writecell_async shall not be called with ncol=0, aborted.\n");
	}else{
		writedata_by_id(NULL, A, 0, ncol);
	}
}
/**
   A generic routine for write data to file with separate header
 */
void writebin_header(cell* Ac, const char* header, const char* format, ...){
	format2fn;
	if(Ac && header){
		free(Ac->header);
		Ac->header=strdup(header);
	}
	write_by_id(Ac, M_0, "%s", fn);
}
/**
 * Read basic array. Not cell. header is already read from file
 * */
static void *readdata_mat(file_t*fp, M_ID id, header_t *header){
	//dbg("calling readdata_mat with dimension %ldx%ld\n", header->nx, header->ny);
	void *out;
	if(!id) id=header->magic;
	switch(id){
	case M_DBL:
		out=dreaddata(fp, header);break;
	case M_FLT:
		out=sreaddata(fp, header);break;
	case M_CMP:
		out=creaddata(fp, header);break;
	case M_ZMP:
		out=zreaddata(fp, header);break;
	case M_LONG:
		out=lreaddata(fp, header);break;
	case M_LOC64: case M_LOC32:
		out=locreaddata(fp, header); break;
	case M_MAP64: case M_MAP32:
	{
		dmat *tmp=dreaddata(fp, header);
		out=d2map(tmp);
		dfree(tmp);
	}
	break;
	case M_RECTMAP64: case M_RECTMAP32:
	{
		dmat *tmp=dreaddata(fp, header);
		out=d2rmap(tmp);
		dfree(tmp);
	}
	break;
	case M_DSP32: case M_DSP64: /**Possible to read mismatched integer*/
		out=dspreaddata(fp, header);break;
	case M_SSP32: case M_SSP64:
		out=sspreaddata(fp, header);break;
	case M_CSP32: case M_CSP64:
		out=cspreaddata(fp, header);break;
	case M_ZSP64: case M_ZSP32:
		out=zspreaddata(fp, header);break;
	default:
		warning("data type id=%ux not supported\n", id);
		out=NULL;
	}
	return out;
}
/*
	read cell. do not scan.  header is already read from file
	id is intended only for the fundamental data
*/
static cell *readdata_cell(file_t *fp, M_ID id, header_t *header){
	long nx=header->nx;
	long ny=header->ny;
	//dbg("calling readdata_cell with dimension %ldx%ld\n", nx, ny);
	cell *dcout=cellnew(nx, ny);
	dcout->header=header->str; header->str=0;
	header_t headerc={0,0,0,0};
	int ans=0;
	
	for(long i=0; i<nx*ny; i++){
		if((ans=read_header(&headerc, fp))){
			dbg3("read_header failed:%s, err=%d, i=%ld\n", zfname(fp), ans, i);
			break;
		}
		//Avoid calling readdata_by_id to avoid scan. 
		if(iscell(&headerc.magic)){
			P(dcout, i)=readdata_cell(fp, id, &headerc);
		}else{
			P(dcout, i)=readdata_mat(fp, id, &headerc);
		}
	}
	return dcout;
}
//scan file for next set of data. first header is already read.
static cell *readdata_scan(file_t *fp, M_ID id, header_t *header){
	cell *dcout=0;
	int maxlen=10;
	void **tmp=mymalloc(maxlen, void *);
	int nx=0;
	do{
		void *tmp2=iscell(&header->magic)?readdata_cell(fp, id, header):readdata_mat(fp, id, header);
		free(header->str);header->str=0;
		if(tmp2){
			if(nx>=maxlen){
				maxlen*=2;
				tmp=myrealloc(tmp, maxlen, void *);
			}
			tmp[nx++]=tmp2;
		} else{
			break;
		}
	} while(!read_header(header, fp));
	
	if(nx>0){
		dcout=cellnew(nx, 1);
		memcpy(P(dcout), tmp, sizeof(void *)*nx);
	}
	free(tmp);
	return dcout;
}
/**
 * Read data from file. level indicate intension, which may not match the file, in which case conversion is done.
 * Usage: level=0: array of fundamental data. 
 * 		  level=1: cell of fundamental data. 
 * 		  level>1: cell of cell ...
 * 	      level==-1: automatically identify and read all data. top level non-cell data will be combined to a cell (e.g., fits)
 * 		  level<-1: only from recursive calling when the initial level is -1
 * 
 * id: magic number of request fundamental data. data is converted if it does not match magic number from the file. 
 * 
 * There are two cases the file may be scanned to real all data, only at the top level though
 * 1. Found array of fundamental data by level is not 0. Used to wrap mat to cell or read fits file extension
 * 2. Found cell array with zero dimension. 
 */
cell* readdata_by_id(file_t* fp, M_ID id, int level, header_t* header){
	//dbg("level=%d\n", level);
	header_t header2={0,0,0,0};
	if(!header){
		header=&header2;
	}
	//if header is not supplied, read from file
	if(header->magic==0 && read_header(header, fp)){
		return NULL;//read header failed
	}
	void* out=0;
	//scan file in automatic mode or when cell dimension is 0 or when request cell but data is not cell.
	
	if(!iscell(&header->magic)){//data is not cell
		if(level==0||level<-1){// and no scan is desired
			out=readdata_mat(fp, id, header);
		}else{//want to scan or want cell array
			out=readdata_scan(fp, id, header);
		}
	}else if(header->nx && header->ny){//cell with deterministic size
		out=readdata_cell(fp, id, header);
	}else if(!read_header(header, fp)){
		out=readdata_scan(fp, id, header);
		//error("unhandled: fn=%s level=%d, magic=%x, nx=%ld, ny=%ld\n", zfname(fp), level, header->magic, (long)header->nx, (long)header->ny);
	}
	
	if(level==0){//want mat, but found cell
	 	while(out && iscell(out)){
			cell *dcout=(cell*)out;
			if(PN(dcout)==0){
				out=NULL;
			}else if(PN(dcout)>=1){
				if(PN(dcout)>1){
					warning("%s: using the first cell when multiple is found\n", zfname(fp));
				}
				out=P(dcout,0);
				P(dcout,0)=0;
				cellfree(dcout);
			}
		}
	}else if (level==-1&&!iscell(&header->magic)&&out&&iscell(out)&&PN((cell*)out)==1){
		//automatic scan mode with only one element, do not wrap into array
		//dbg("level=%d, scan mode find only a single element.\n", level);
		cell* out2=P((cell*)out, 0);
		P((cell*)out, 0)=0;
		cellfree(out);
		out=out2;
	}
	free(header->str);header->str=0;
	header->magic=0; header->nx=0; header->ny=0;//prevent reuse.
	return (cell*)out;
}

cell* read_by_id(M_ID id, int level, const char* format, ...){
	format2fn;
	file_t* fp;
	if(!fn||!(fp=zfopen(fn, "rb"))){
		dbg("Unable to read from file %s (empty).\n", fn);
		return NULL;
	}
	
	cell* out=readdata_by_id(fp, id, level, 0);
	zfclose(fp);
	return out;
}

/**
   A generic routine for reading data from socket. User need to cast the result.
   We dup the fd to avoid close it after read.
 */
cell* readsock(int sock){
	file_t* fp=zfdopen(dup(sock));
	cell* out=fp?readdata_by_id(fp, 0, -1, 0):NULL;
	zfclose(fp);
	return out;
}
/**
   A generic routine for write data to socket.
   We dup the fd to avoid close it after read.
 */
void writesock(const cell* A, int sock){
	file_t* fp=zfdopen(dup(sock));
	if(fp) writedata_by_id(fp, A, 0, 0);
	zfclose(fp);
}

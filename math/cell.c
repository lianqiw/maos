/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
void cellinit(panyarray A_, long nx, long ny){
	cell** A=A_.c;
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
void cellinit2(panyarray A_, const_anyarray B_){
	const cell* B=B_.c;
	cell** A=A_.c;
	if(!iscell(B)){
		error("Invalid usage. B must be cell\n");
	}
	if(*A){
		cellinit(A, B->nx, B->ny);
	} else if(iscell(B)){
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
 * @brief Construct a reference base on type.
 * 
 */
cell* cellref(const_anyarray in){
	if(!in.c) return NULL;
	if(iscell(in.c)){
		const cell *A=in.c;
		cell *out=cellnew(A->nx, A->ny);
		for(long i=0; i<PN(A); i++){
			P(out, i)=cellref(P(A, i));
		}
		if(A->keywords) out->keywords=strdup(A->keywords);
		if(A->m) out->m=cellref(A->m);
		return out;
	}else switch(in.c->id){
		case M_DBL: return (cell*)dref(in.dm); break;
		case M_CMP: return (cell*)cref(in.cm); break;
		case M_FLT: return (cell*)sref(in.sm); break;
		case M_ZMP: return (cell*)zref(in.zm); break;
		case M_LONG: return (cell*)lref(in.lm); break;
		case M_DSP: return (cell*)dspref(in.ds); break;
		case M_CSP: return (cell*)cspref(in.cs); break;
		case M_SSP: return (cell*)sspref(in.ss); break;
		case M_ZSP: return (cell*)zspref(in.zs); break;
		default:error("Unknown id: %u\n", in.c->id);
		return NULL;
	}
}
uint32_t cellhash(const_anyarray in, uint32_t key){
	if(!in.c) return key;
	if(iscell(in.c)){
		const cell *A=in.c;
		for(long i=0; i<PN(A); i++){
			key=cellhash(P(A,i), key);
		}
	}else switch(in.c->id&0xFFFF){
		case M_DBL: key=dhash(in.dm, key); break;
		case M_CMP: key=chash(in.cm, key); break;
		case M_FLT: key=shash(in.sm, key); break;
		case M_ZMP: key=zhash(in.zm, key); break;
		default: 
			warning("Hash is skipped: to implement\n");
	}
	return key;
}
/**
 * @brief Convert elements of cell array
 * 
 */
cell *cellconvert(cell *A, cell* (*fun_convert)(cell*)){
	if(!A) return NULL;
	cell *out;
	if(iscell(A)){
		out=cellnew(A->nx, A->ny);
		for(long i=0; i<PN(A); i++){
			P(out, i)=cellconvert(P(A,i), fun_convert);
		}
	}else{
		out=fun_convert(A);
	}
	return out;
}

/**
   Obtain the dimensions.
*/
void celldim(const_anyarray A_, long* nx, long* ny, long** nxs, long** nys){
	const cell* A=A_.c;
	*nx=0;
	*ny=0;
	if(!A || !iscell(A) || !A->nx || !A->ny){
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
void cellresize(anyarray A_, long nx, long ny){
	cell* A=A_.c;
	if(!A || A->m){ 
		error("Can only resize a valid cell array with m unset.\n");
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
			memset(A->p+nold, 0, (nnew-nold)*sizeof(cell*));
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
		free(A->p);
		A->p=p;
	}
	A->nx=nx;
	A->ny=ny;
}

/**Check whether the cell has only diagonal entries.
 */
int cell_is_diag(const_anyarray A_){
	const cell *A=A_.c;
	if(!iscell(A)) return 0;
	for(int iy=0; iy<NY(A); iy++){
		for(int ix=0; ix<NX(A); ix++){
			if(ix!=iy){//off diagonal should be empty.
				if(P(A, ix, iy) && PN(A, ix, iy)){
					return 0;
				}
			}
		}
	}
	return 1;
}
/**
   free a mat or cell object.
*/
void cellfree_do(anyarray A_){
	cell* A=A_.c;
	if(!A) return;
	M_ID id=A->id;
	switch(id&0xFFFF){
	case MCC_ANY:{
		cell* dc=A;
		for(int ix=0; ix<PN(dc); ix++){
			cellfree_do(P(dc,ix));
		}
		free(P(dc));P(dc)=0;
		if(dc->keywords){
			free(dc->keywords); dc->keywords=NULL;		
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
		dfree_do(A_.dm);break;
	case M_CMP:
		cfree_do(A_.cm);break;
	case M_FLT:
		sfree_do(A_.sm);break;
	case M_ZMP:
		zfree_do(A_.zm);break;
	case M_LONG:
		lfree_do(A_.lm);break;
	case M_DSP:
		dspfree_do(A_.ds);break;
	case M_SSP:
		sspfree_do(A_.ss);break;
	case M_CSP:
		cspfree_do(A_.cs);break;
	case M_ZSP:
		zspfree_do(A_.zs);break;
	default:
		warning("Unknown id=%x, A=%p\n", id, A);
		print_backtrace();
	}
}
/**
 * ncol: 0: normal writing. -1: initialize async data. >0: async writing.
 * 	
 * */
void writedata(file_t* fp, const_anyarray A_, long ncol){
	if(fp&&ncol>0){
		error("writedata_by_id should be called with either fp or ncol, but not both, aborted.\n");
		return;
	}
	if(!A_.c || iscell(A_.c)){//empty data is treated as cell
		const cell* A=A_.c;
		if(ncol<=0){//not async write.
			header_t header={MCC_ANY, NX(A), NY(A), A?A->keywords:NULL};
			write_header(&header, fp);
		}
		for(long ix=0; ix<PN(A); ix++){
			writedata(fp, P(A,ix), ncol);
		}
	}else switch(A_.c->id&0xFFFF){
		case M_DBL:
			dwritedata(fp, A_.dm, ncol);break;
		case M_CMP:
			cwritedata(fp, A_.cm, ncol);break;
		case M_FLT:
			swritedata(fp, A_.sm, ncol);break;
		case M_ZMP:
			zwritedata(fp, A_.zm, ncol);break;
		case M_LONG:
			lwritedata(fp, A_.lm, ncol);break;
		case M_DSP:
			dspwritedata(fp, A_.ds);break;
		case M_SSP:
			sspwritedata(fp, A_.ss);break;
		case M_CSP:
			cspwritedata(fp, A_.cs);break;
		case M_ZSP:
			zspwritedata(fp, A_.zs);break;
		default:
			error("Unknown id=%u\n", A_.c->id);
	}
}

void writebin(const_anyarray A, const char* format, ...){
	format2fn;
	if(!fn) return;
	file_t* fp=zfopen(fn, "wb");
	if(!fp) return;
	writedata(fp, A, 0);
	zfclose(fp);
}

/**
 * Calls writebin with build in fp
 * */
void writecell_async(const_anyarray A, long ncol){
	if(ncol==0){
		dbg("writecell_async should not be called with ncol=0.\n");
	}else{
		writedata(NULL, A, ncol);
	}
}
/**
   A generic routine for write data to file with separate keywords
 */
void writebin_header(anyarray Ac_, const char* keywords, const char* format, ...){
	format2fn;
	cell* Ac=Ac_.c;
	if(Ac && keywords){
		free(Ac->keywords);
		Ac->keywords=strdup(keywords);
	}
	writebin(Ac, "%s", fn);
}

/*
	Read a block which may be cell array or basic arrys. do not scan.  id is
	intended only for the fundamental data. header is read first if header==0 or header->magic==0.
*/
cell *readdata(file_t *fp, M_ID id, header_t *header){
	header_t header2={0,0,0,0};
	if(!header){//header is not yet read
		header=&header2;
	}
	if(!header->magic){
		if(read_header(header, fp)) return NULL;
	}
	if(!header->nx || !header->ny){
		if(!header->str){
			return NULL;
		}
	}
	long nx=header->nx;
	long ny=header->ny;
	cell *res=NULL;
	//dbg("calling readdata_cell with dimension %ldx%ld\n", nx, ny);
	if(iscell(&header->magic)){
		cell *dcout=cellnew(nx, ny);
		dcout->keywords=header->str; header->str=0;
		header_t headerc={0,0,0,0};
		int ans=0;
		
		for(long i=0; i<nx*ny; i++){
			if((ans=read_header(&headerc, fp))){
				dbg3("read_header failed:%s, err=%d, i=%ld\n", zfname(fp), ans, i);
				break;
			}
			//Avoid calling readdata_auto to avoid scan. 
			P(dcout, i)=readdata(fp, id, &headerc);
		}
		res=dcout;
	}else{
		void *out=NULL;
		
		if(!id) id=header->magic;
		switch(id&0xFFFF){
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
	
		res=out;
	}
	return res;
}
//scan file for next set of data. first header is already read.
static cell *readdata_scan(file_t *fp, M_ID id, header_t *header){
	cell *dcout=0;
	int maxlen=10;
	cell **tmp=mymalloc(maxlen, cell*);
	int nx=0;
	do{
		cell *tmp2=readdata(fp, id, header);
		free(header->str);header->str=0;
		if(tmp2){
			if(nx>=maxlen){
				maxlen*=2;
				tmp=myrealloc(tmp, maxlen, cell*);
			}
			tmp[nx++]=tmp2;
		} else{
			break;
		}
	} while(!read_header(header, fp));
	
	if(nx>0){
		dcout=cellnew(nx, 1);
		memcpy(P(dcout), tmp, sizeof(cell*)*nx);
	}
	free(tmp);
	return dcout;
}
/**
 * Read data from file. level indicate intention, which may not match the file, in which case conversion is done.
 * Usage: level=0: array of fundamental data. 
 * 		  level=1: cell of fundamental data. 
 * 		  level>1: cell of cell ...
 * 	      level==-1: automatically identify and read all data. top level non-cell data will be combined to a cell (e.g., fits)
 * 		  level<-1: only from recursive calling when the initial level is -1
 * 
 * id: magic number of request fundamental data. data is converted if it does not match magic number from the file. 
 * 
 * There are two cases the file may be scanned to real all data, only at the top level though
 * 1. Found array of fundamental data but level is not 0. Used to wrap mat to cell or read fits file extension
 * 2. Found cell array with zero dimension. 
 * 
 * Function should be called at the top level only (right after opening file or socket)
 */
static cell* readdata_auto(file_t* fp, M_ID id, int level, header_t* header){
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
			out=readdata(fp, id, header);
		}else{//want to scan or want cell array
			out=readdata_scan(fp, id, header);
		}
	}else if(header->nx && header->ny){//cell with deterministic size
		out=readdata(fp, id, header);
	}else if(!read_header(header, fp)){//cell with unknown size (0). scan. 
		out=readdata_scan(fp, id, header);
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

cell* readbin_id(M_ID id, int level, const char* format, ...){
	format2fn;
	file_t* fp=NULL;
	if(!fn||!(fp=zfopen(fn, "rb"))){
		dbg("Unable to read from file %s (empty).\n", fn);
		return NULL;
	}
	cell* out=readdata_auto(fp, id, level, 0);
	zfclose(fp);
	return out;
}
cell* readbin(const char* format, ...){
	format2fn;
	return readbin_id(M_0, -1, "%s", fn);
}
/**
   A generic routine for reading data from socket. User need to cast the result.
   We dup the fd to avoid close it after read.
 */
cell* readsock(int sock){
	file_t* fp=zfdopen(dup(sock));
	cell* out=fp?readdata_auto(fp, 0, -1, 0):NULL;
	zfclose(fp);
	return out;
}
/**
   A generic routine for write data to socket.
   We dup the fd to avoid close it after read.
 */
void writesock(const_anyarray A, int sock){
	file_t* fp=zfdopen(dup(sock));
	if(fp) writedata(fp, A, 0);
	zfclose(fp);
}

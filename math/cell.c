/*
  Copyright 2009-2026 Lianqi Wang
  
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

   @param[in] nx Number of rows. Negative is treated as 0.
   @param[in] ny Number of columns. Negative is treated as 0.
   @return Newly allocated cell array with NULL elements; the caller frees it with cellfree().
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

   @param[in,out] A_ Pointer to the cell array. Created if NULL, otherwise its type and size are checked.
   @param[in] nx Number of rows.
   @param[in] ny Number of columns.
   @return 0 on success, -1 if the existing cell array does not match.
*/
int cellinit(panyarray A_, long nx, long ny){
	cell** A=A_.c;
	if(!A) return -1;
	if(*A){
		if(!(iscell(*A)&&(*A)->nx==nx&&(*A)->ny==ny)){
			error("type or size mismatch: id=%u, require (%ld, %ld), A is (%ld, %ld)\n",
				(*A)->id, nx, ny, (*A)->nx, (*A)->ny);
			return -1;
		}
	} else{
		*A=cellnew(nx, ny);
	}
	return 0;
}
/**
   check the size of cell array if exist. Otherwise create it with the same layout as B. m is created for
   dcell

   @param[in,out] A_ Pointer to the cell array. Created if NULL, otherwise its type and size are checked.
   @param[in] B_ The reference cell array. Must be a cell.
*/
void cellinit2(panyarray A_, const_anyarray B_){
	const cell* B=B_.c;
	cell** A=A_.c;
	if(!iscell(B)){
		error("Invalid usage. B must be cell\n");
		return;
	}
	if(*A){
		cellinit(A, B->nx, B->ny);
	} else if(iscell(B)){
		M_ID magic=(M_ID)0;
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
   Construct a reference base on type. Cell arrays are recreated with each element referenced.

   @param[in] in Array or cell array to reference.
   @return Newly allocated reference sharing the data with in; the caller frees it. NULL if in is NULL.
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
/**
   Compute a hash of the numerical content of array or cell array in (recursively)
   using key as the seed. Returns the updated key. Sparse arrays are skipped.

   @param[in] in Array or cell array.
   @param[in] key Seed of the hash.
   @return The updated hash key.
*/
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
   Convert elements of cell array in place.

   @param[in,out] A Array or cell array whose (non-cell) elements are replaced by the result of fun_convert.
   @param[in] fun_convert Conversion function. It takes ownership of its input and returns the converted
   array.
   @return The converted array (A itself if A is a cell array). NULL if A is NULL.
*/
cell *cellconvert(cell *A, cell* (*fun_convert)(cell*)){
	if(!A) return NULL;
	if(iscell(A)){
		for(long i=0; i<PN(A); i++){
			P(A, i)=cellconvert(P(A,i), fun_convert);
		}
	}else{
		A=fun_convert(A);
	}
	return A;
}

/**
   Obtain the dimensions of the dense matrix assembled from cell array A. nx, ny
   are the total rows and columns. *nxs and *nys are allocated (caller frees) to
   contain the number of rows of each cell row and columns of each cell column.
   Cells along a row (column) must have the same number of rows (columns).

   @param[in] A_ The cell array.
   @param[out] nx Total number of rows.
   @param[out] ny Total number of columns.
   @param[out] nxs Newly allocated number of rows of each cell row; the caller frees it. NULL if A is empty.
   @param[out] nys Newly allocated number of columns of each cell column; the caller frees it. NULL if A is
   empty.
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
					goto fail;
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
					goto fail;
				}
			}
		}
	}
	return;
fail:
	free(*nxs); *nxs=NULL;
	free(*nys); *nys=NULL;
	*nx=0; *ny=0;
}
/**
   Resize a generic cell array. Elements outside of the new size are freed and new elements are NULL.

   @param[in,out] A_ The cell array. Must be valid and with m unset.
   @param[in] nx New number of rows.
   @param[in] ny New number of columns.
*/
void cellresize(anyarray A_, long nx, long ny){
	cell* A=A_.c;
	if(!A || A->m){ 
		error("Can only resize a valid cell array with m unset.\n");
		return;
	}
	if(A->nx==nx||(A->ny==1&&(ny==1||nx>=A->nx))){//layout is preserved by linear realloc
		long nold=A->nx*A->ny;
		long nnew=nx*ny;
		if(nnew<nold&&iscell(A)){
			for(long i=nnew; i<nold; i++){
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

/**
   Check whether the cell has only diagonal entries.

   @param[in] A_ The cell array.
   @return 1 if all off diagonal elements are empty, 0 otherwise or if A_ is not a cell array.
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

   @param[in,out] A_ The object to free (recursively for cell arrays). May be NULL.
*/
void cellfree_do(anyarray A_){
	cell* A=A_.c;
	if(!A) return;
	M_ID id=A->id;
	switch(id&0xFFFF){
	case MCC_ANY:{
		cell* dc=A;
		if(dc->fp && zfisfits(dc->fp)){//async writing is not supported for fits. write all data now.
			writedata(dc->fp, dc, 0);
		}
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
   Write array or cell array to file.

   @param[in] fp File to write to. Must be NULL when ncol>0.
   @param[in] A_ Array or cell array to write.
   @param[in] ncol 0: normal writing. -1: initialize async data. >0: async writing of the first ncol columns.
*/
void writedata(file_t* fp, const_anyarray A_, long ncol){
	if(fp&&ncol>0){
		error("writedata_by_id should be called with either fp or ncol, but not both, aborted.\n");
		return;
	}
	if(!A_.c || iscell(A_.c)){//empty data is treated as cell
		const cell* A=A_.c;
		if(ncol==-1 && fp && zfisfits(fp)){
			//async writing is not supported for fits. data is written when A is freed (see cellfree_do).
			return;
		}
		if(ncol<=0){//not async write.
			header_t header={MCC_ANY, (uint64_t)NX(A), (uint64_t)NY(A), A?A->keywords:NULL};
			write_header(&header, fp);
		}
		uint32_t id=0;//find the type
		for(long ix=0; ix<PN(A); ix++){
			if(P(A,ix)){
				id=P(A,ix)->id;
				break;
			}
		}
		for(long ix=0; ix<PN(A); ix++){
			if(P(A,ix)){
				writedata(fp, P(A,ix), ncol);
			}else if(ncol<=0){//write empty (not for async writing, where fp is NULL)
				header_t header={id, (uint64_t)0, (uint64_t)0, NULL};
				write_header(&header, fp);
			}
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
		default:{
			const cell *A=A_.c;
			const int id=A->id&0xFFFF;
			if(id==M_INT64||id==M_INT32||id==M_INT16||id==M_UINT8){//other integer arrays, e.g., from python
				const size_t size=id==M_INT64?8:(id==M_INT32?4:(id==M_INT16?2:1));
				if(ncol<=0) writearr(fp, 0, size, A->id, A->keywords, P(A), NX(A), NY(A));
			}else{
				error("Unknown id=%u\n", A->id);
			}
		}
	}
}

/**
   Write array or cell array A to file with name given by the printf-style
   format. Does nothing if the file name is empty or cannot be opened.

   @param[in] A Array or cell array to write.
   @param[in] format printf-style format of the file name.
*/
void writebin(const_anyarray A, const char* format, ...){
	format2fn;
	if(!fn) return;
	file_t* fp=zfopen(fn, "wb");
	if(!fp) return;
	writedata(fp, A, 0);
	zfclose(fp);
}

/**
   Asynchronously write data to the file previously initialized with writedata() with ncol=-1.

   @param[in] A Array or cell array to write.
   @param[in] ncol Number of columns that are ready to write. Must not be 0.
*/
void writecell_async(const_anyarray A, long ncol){
	if(ncol==0){
		dbg("writecell_async should not be called with ncol=0.\n");
	}else{
		writedata(NULL, A, ncol);
	}
}
/**
   A generic routine for write data to file with separate keywords

   @param[in,out] Ac_ Array or cell array to write; its keywords are replaced by keywords if not NULL.
   @param[in] keywords Keywords to write.
   @param[in] format printf-style format of the file name.
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

/**
	Read a block which may be cell array or basic arrys. do not scan.  id is
	intended only for the fundamental data. header is read first if header==0 or header->id==0.

   @param[in] fp File to read from.
   @param[in] id Requested type of the fundamental data. 0 to use the type in the file
   (integer arrays of all sizes are read as lmat).
   @param[in,out] header Header of the block. May be NULL. Its str is taken by the returned cell array.
   @return Newly allocated data; the caller frees it. NULL if empty or if reading fails.
*/
cell *readdata(file_t *fp, M_ID id, header_t *header){
	header_t header2={0,0,0,0};
	if(!header){//header is not yet read
		header=&header2;
	}
	if(!header->id){
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
	if(iscell(&header->id)){
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
		
		if(!id){
			id=header->id;
			switch(id&0xFFFF){//integer arrays of all sizes are read as lmat
			case M_INT64: case M_INT32: case M_INT16: case M_UINT8:
				id=M_LONG;break;
			}
		}
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
	
		res=(cell*)out;
	}
	return res;
}
/**
   scan file for next set of data. first header is already read.

   @param[in] fp File to read from.
   @param[in] id Requested type of the fundamental data. 0 to use the type in the file.
   @param[in,out] header The already read header. Used to read the following headers.
   @return Newly allocated nx*1 cell array containing all the blocks; the caller frees it. NULL if none.
*/
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
   Read data from file. level indicate intention, which may not match the file, in which case conversion is done.

   There are two cases the file may be scanned to real all data, only at the top level though
   1. Found array of fundamental data but level is not 0. Used to wrap mat to cell or read fits file extension
   2. Found cell array with zero dimension.

   Function should be called at the top level only (right after opening file or socket)

   @param[in] fp File to read from.
   @param[in] id Magic number of request fundamental data. data is converted if it does not match magic
   number from the file.
   @param[in] level 0: array of fundamental data. 1: cell of fundamental data. >1: cell of cell ... -1:
   automatically identify and read all data; top level non-cell data will be combined to a cell (e.g.,
   fits). <-1: only from recursive calling when the initial level is -1.
   @param[in,out] header Header that is already read. May be NULL. It is reset on return.
   @return Newly allocated data; the caller frees it. NULL if reading fails.
*/
static cell* readdata_auto(file_t* fp, M_ID id, int level, header_t* header){
	//dbg("level=%d\n", level);
	header_t header2={0,0,0,0};
	if(!header){
		header=&header2;
	}
	//if header is not supplied, read from file
	if(header->id==0 && read_header(header, fp)){
		return NULL;//read header failed
	}
	void* out=0;
	//scan file in automatic mode or when cell dimension is 0 or when request cell but data is not cell.
	
	if(!iscell(&header->id)){//data is not cell
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
				cellfree(dcout);
			}else if(PN(dcout)>=1){
				if(PN(dcout)>1){
					warning("%s: using the first cell when multiple is found\n", zfname(fp));
				}
				out=P(dcout,0);
				P(dcout,0)=0;
				cellfree(dcout);
			}
		}
	}else if (level==-1&&!iscell(&header->id)&&out&&iscell(out)&&PN((cell*)out)==1){
		//automatic scan mode with only one element, do not wrap into array
		//dbg("level=%d, scan mode find only a single element.\n", level);
		cell* out2=P((cell*)out, 0);
		P((cell*)out, 0)=0;
		cellfree(out);
		out=out2;
	}
	free(header->str);header->str=0;
	header->id=0; header->nx=0; header->ny=0;//prevent reuse.
	return (cell*)out;
}

/**
   Read array or cell array from an opened file, socket, or pipe, like readbin_id().

   @param[in] fp File to read from, e.g., from zfdopen_as().
   @param[in] id Requested fundamental data type. 0 for native.
   @param[in] level Cell nesting level, see readdata_auto(). -1 to detect automatically.
   @return Newly allocated data; the caller frees it. NULL if empty or if reading fails.
*/
cell* readbin_fp(file_t* fp, M_ID id, int level){
	return fp?readdata_auto(fp, id, level, 0):NULL;
}
/**
   Read array or cell array from file with name given by the printf-style format.
   id is the requested fundamental data type (converted if necessary, 0 for
   native), level is the cell nesting level, see readdata_auto(). Returns NULL
   if the file cannot be opened. The caller owns the result.

   @param[in] id Requested fundamental data type. 0 for native.
   @param[in] level Cell nesting level, see readdata_auto().
   @param[in] format printf-style format of the file name.
   @return Newly allocated data; the caller frees it. NULL if the file cannot be opened.
*/
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
/**
   Read any data from file with automatic type and nesting detection. The
   caller needs to cast the result and free it.

   @param[in] format printf-style format of the file name.
   @return Newly allocated data; the caller frees it. NULL if the file cannot be opened.
*/
cell* readbin(const char* format, ...){
	format2fn;
	return readbin_id(M_0, -1, "%s", fn);
}
/**
   A generic routine for reading data from socket. User need to cast the result.
   We dup the fd to avoid close it after read.

   @param[in] sock The socket.
   @return Newly allocated data; the caller frees it. NULL on failure.
*/
cell* readsock(int sock){
	file_t* fp=zfdopen(dup(sock));
	cell* out=fp?readdata_auto(fp, (M_ID)0, -1, 0):NULL;
	zfclose(fp);
	return out;
}
/**
   A generic routine for write data to socket.
   We dup the fd to avoid close it after read.

   @param[in] A Array or cell array to write.
   @param[in] sock The socket.
*/
void writesock(const_anyarray A, int sock){
	file_t* fp=zfdopen(dup(sock));
	if(fp) writedata(fp, A, 0);
	zfclose(fp);
}

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
		dc->p=mycalloc(nx*ny, cell*);
	}
	return dc;
}

/**
   Allocate a new array of the same type
 */
cell* cellnew2(const void* A_){
	const cell* A=(const cell*)A_;
	if(!A){
		return 0;
	} else if(iscell(A)){
		return cellnew(A->nx, A->ny);
	} else if(A->id==M_REAL){
		return (cell*)dnew(A->nx, A->ny);
	} else{
		error("Invalid type: id=%u\n", A->id);
		return 0;
	}
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
		uint32_t magic=0;
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
void cellresize(void* A_, long nx, long ny){
	if(!A_) return;
	cell* A=cell_cast(A_);
	if(A->nx==nx||A->ny==1){
		int nold=A->nx*A->ny;
		int nnew=nx*ny;
		if(nnew<nold&&iscell(A)){
			for(int i=nnew; i<nold; i++){
				cellfree_do(P(A,i));
			}
		}
		A->p=myrealloc(A->p, nnew, cell*);
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
   free a mat or cell object.
*/
void cellfree_do(void* A){
	if(!A) return;
	uint32_t id=((cell*)A)->id;
	switch(id){
	case MCC_ANY:{
		cell* dc=(cell*)A;
		if(dc->p){
			for(int ix=0; ix<dc->nx*dc->ny; ix++){
				cellfree_do(P(dc,ix));
			}
			free(dc->p);dc->p=0;
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
		print_backtrace();
		warning("Unknown id=%x, A=%p\n", id, A);
	}
}

void writedata_by_id(file_t* fp, const void* A_, uint32_t id){
	if(!fp) return;
	const cell* A=(const cell*)A_;
	if(A){
		if(!id){
			id=A->id;
		} else if(id!=MCC_ANY){
			uint32_t id2=A->id;
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
		uint64_t nx=0;
		uint64_t ny=0;
		if(A){
			nx=A->nx;
			ny=A->ny;
		}
		id=0;/*determine id first for non-empty cell*/
		if(nx&&ny){
			for(uint64_t ix=0; ix<(nx*ny); ix++){
				if(P(A,ix)){
					id=P(A,ix)->id;
					if(!id){
						warning("id is not set\n");
					}else{
						break;
					}
				}
			}
			if(!id){//empty cell
				nx=0; ny=0;
			}
		}
		header_t header={MCC_ANY, nx, ny, A?A->header:NULL};
		write_header(&header, fp);
		if(id){
			for(long ix=0; ix<(A->nx*A->ny); ix++){
				writedata_by_id(fp, P(A,ix), 0);
			}
		}
	}
				break;
	case M_DBL:
		dwritedata(fp, (dmat*)A);break;
	case M_CMP:
		cwritedata(fp, (cmat*)A);break;
	case M_FLT:
		swritedata(fp, (smat*)A);break;
	case M_ZMP:
		zwritedata(fp, (zmat*)A);break;
	case M_LONG:
		lwritedata(fp, (lmat*)A);break;
	case M_LOC:
		locwritedata(fp, (loc_t*)A);break;
	case M_MAP:
		map_header((map_t*)A);
		dwritedata(fp, (dmat*)A);break;
	case M_RECTMAP:
		rmap_header((rmap_t*)A);
		dwritedata(fp, (dmat*)A);break;
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

void write_by_id(const void* A, uint32_t id, const char* format, ...){
	format2fn;
	file_t* fp=zfopen(fn, "wb");
	if(!fp) return;
	writedata_by_id(fp, A, id);
	zfclose(fp);
}
/**
 * Read data from file. 
 * Usage: level=0: array of fundamental data. 
 * 		  level=1: cell of fundamental data. 
 * 		  level>1: cell of cell ...
 * 	      level<0: automatically identify and read all data
 * 
 * id: magic number of request fundamental data. data is converted if it does not match magic number from the file. 
 * Cell dimension may be zero. In this case the dimension will be automatically detected.
 *  */
cell* readdata_by_id(file_t* fp, uint32_t id, int level, header_t* header){
	header_t header2={0,0,0,0};
	if(!header){
		header=&header2;
	}
	if(header->magic==0 && read_header(header, fp)){
		return NULL;
	}
	//info("level=%d, magic=%x\n", level, header->magic);
	void* out=0;
	//scan file in automatic mode or when cell dimension is 0 or when request cell but data is not cell.
	if((level==0||level<-1)&&!iscell(&header->magic)){
		//info("level==0 && !iscell\n");
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
			dmat* tmp=dreaddata(fp, header);
			out=d2map(tmp);
			dfree(tmp);
		}
		break;
		case M_RECTMAP64: case M_RECTMAP32:
		{
			dmat* tmp=dreaddata(fp, header);
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
			warning("data type id=%u not supported\n", id);
			out=NULL;
		}
	}else if(iscell(&header->magic) && header->nx>0 && header->ny>0){//cell array with known dimensions
		//info("iscell && nx, ny >0\n");
		long nx=header->nx;
		long ny=header->ny;
		cell* dcout=cellnew(nx, ny);
		dcout->header=header->str; header->str=0;
		header_t headerc={0,0,0,0};
		for(long i=0; i<nx*ny; i++){
			if(read_header(&headerc, fp)){
				break;
			}
			if(!headerc.str&&dcout->header){//copy str from cell to mat.
				headerc.str=strdup(dcout->header);
			}
			P(dcout,i)=readdata_by_id(fp, id, level-1, &headerc);
		}
		out=dcout;
	}else if(level>-2){//scan file only when auto.
		if(!iscell(&header->magic)||!read_header(header, fp)){
			int maxlen=10;
			void** tmp=mymalloc(maxlen, void*);
			int nx=0;
			do{
				if(nx>=maxlen){
					maxlen*=2;
					tmp=myrealloc(tmp, maxlen, void*);
				}
				void* tmp2=readdata_by_id(fp, id, level-1, header);
				free(header->str);header->str=0;
				if(tmp2){
					tmp[nx++]=tmp2;
				} else{
					break;
				}
			} while(!read_header(header, fp));
			cell* dcout=cellnew(nx, 1);
			memcpy(dcout->p, tmp, sizeof(void*)*nx);
			free(tmp);
			out=dcout;
		}
	}
	free(header->str);header->str=0;
	if(level==0){//take first none-cell element.
		while(out && iscell(out)){
			cell* out2=P((cell*)out,0);
			P((cell*)out,0)=0;
			cellfree(out);
			out=out2;
		}
	} else if(level<0&&!iscell(&header->magic)&&out&&iscell(out)&&((cell*)out)->nx==1&&((cell*)out)->ny==1){
		//automatic mode with only one element
		cell* out2=P((cell*)out, 0);
		P((cell*)out, 0)=0;
		cellfree(out);
		out=out2;
	}
	header->magic=0; header->nx=0; header->ny=0;//prevent reuse.
	return (cell*)out;
}

cell* read_by_id(uint32_t id, int level, const char* format, ...){
	format2fn;
	file_t* fp=zfopen(fn, "rb");
	if(!fp) return 0;
	cell* out=readdata_by_id(fp, id, level, 0);
	zfclose(fp);
	return out;
}
/**
   A generic routine for reading data from file. User need to cast the result.
 */
cell* readbin(const char* format, ...){
	format2fn;
	return read_by_id(0, -1, "%s", fn);
}
/**
   A generic routine for write data to file.
 */
void writebin(const void* A, const char* format, ...){
	format2fn;
	write_by_id(A, 0, "%s", fn);
}
/**
   A generic routine for write data to file with separate header
 */
void writebin2(void* A, const char* header, const char* format, ...){
	format2fn;
	cell* Ac=cell_cast(A);
	if(Ac){
		free(Ac->header);
		Ac->header=strdup(header);
	}
	write_by_id(Ac, 0, "%s", fn);
}
/**
   A generic routine for reading data from socket. User need to cast the result.
 */
cell* readsock(int sock){
	file_t* fp=zfdopen(sock, "rb");
	cell* out=readdata_by_id(fp, 0, -1, 0);
	zfclose(fp);
	return out;
}
/**
   A generic routine for write data to socket.
 */
void writesock(const void* A, int sock){
	file_t* fp=zfdopen(sock, "wb");
	writedata_by_id(fp, A, 0);
	zfclose(fp);
}

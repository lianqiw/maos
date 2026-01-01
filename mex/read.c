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
#ifdef __INTEL_COMPILER
#undef _GNU_SOURCE /*avoid compiling problem*/
#endif

static const int skip_unicell=0; //strip single cells: a{1}->a

#include "io.h"
mxArray* SKIPPED=(mxArray*)(1);
static mxArray* readdata(file_t* fp, mxArray** header, int start, int howmany){
	/*
	  if start!=0 || howmany!=0 and data is cell, will only read cell from start to start+howmany
	  if data is not cell and start==-1, will skip the data.
	  Only the first call to readdata will possibly have howmany!=0
	*/
	if(fp->eof) return NULL;
	header_t header2={0,0,0,0,0};
	if(read_header2(&header2, fp)){
		return NULL;
	}
	uint32_t magic=header2.magic&0xFFFF;
	if(magic==0){//end of file or empty file
		fp->eof=1;
		return NULL;
	}
	if(header){
		if(header2.str)
			*header=mxCreateString(header2.str);
		else
			*header=mxCreateString("");
	}
	long ntot=header2.ntot;
	int start_save=start;
	if(iscell(magic)){
		if(start!=-1&&howmany==0){
			if(start!=0){
				error("Invalid use\n");
			}
			howmany=ntot-start;
		}
	} else if(fp->isfits){
		if(search_keyword_int(header2.str, "Sparse=")){
			magic=M_DSP64;
		}
		if(howmany!=0){
			/*first read of fits file. determine if we need it*/
			if(start>0){
				start=-1;//skip this block.
			}
		}
	} else{
		if((start!=0&&start!=-1)||howmany!=0){
			error("invalid use");
		}
	}
	int iscell=0;
	if(fp->eof) return NULL;
	mxArray* out=NULL;
	mwSize byte=0;
	mxClassID id=(mxClassID)0;
	mxComplexity mxFLAG=mxREAL;
	int ibyte=8;
	int issp=0;
	switch(magic){
	case M_DSP64:  byte=8; id=mxDOUBLE_CLASS; issp=1; break;
	case M_DSP32:  byte=8; id=mxDOUBLE_CLASS; issp=1;ibyte=4; break;
	case M_DBL:   byte=8; id=mxDOUBLE_CLASS;break;

	case M_CSP64: byte=16;id=mxDOUBLE_CLASS;mxFLAG=mxCOMPLEX;issp=1;break;
	case M_CSP32: byte=16;id=mxDOUBLE_CLASS;mxFLAG=mxCOMPLEX;issp=1; ibyte=4; break;
	case M_CMP:   byte=16;id=mxDOUBLE_CLASS;mxFLAG=mxCOMPLEX;break;

	case M_SSP64:  byte=4; id=mxSINGLE_CLASS;issp=1;break;
	case M_SSP32:  byte=4; id=mxSINGLE_CLASS;issp=1;ibyte=4; break;
	case M_FLT:   byte=4; id=mxSINGLE_CLASS;break;

	case M_ZSP64: byte=8; id=mxSINGLE_CLASS;mxFLAG=mxCOMPLEX;issp=1;break;
	case M_ZSP32: byte=8; id=mxSINGLE_CLASS;mxFLAG=mxCOMPLEX;issp=1;ibyte=4; break;
	case M_ZMP:   byte=8; id=mxSINGLE_CLASS;mxFLAG=mxCOMPLEX;break;

	case M_INT64: byte=8; id=mxINT64_CLASS; break;
	case M_INT32: byte=4; id=mxINT32_CLASS; break;
	case M_INT16: byte=2; id=mxINT16_CLASS; break;
	case M_INT8:  byte=1; id=mxINT8_CLASS; break;
	default: id=(mxClassID)0;
	}
	if(issp&&start==-1){
		warning("Skipping sparse matrix not implemented yet.\n");
		start=0;
	}
	if(magic>=0x6410&&magic<=0x6424){//any array
		iscell=1;
		long ix;
		if(fp->eof) return NULL;
		/*if(strstr(header2.str,"type=struct"))
		  {
		  out=mxCreateStructMatrix(1,1,ntot);
		  }else*/
		{
			if(ntot>1||!skip_unicell){
				out=mxCreateCellArray(header2.ndim, header2.dims);
			}
		}
		mxArray* header0=mxCreateCellMatrix(ntot+1, 1);
		int nheader0=0;
		for(ix=0; ix<ntot; ix++){
			int start2=0;
			if(start==-1||ix<start||ix>=start+howmany){
				start2=-1;
			}
			mxArray* header3=NULL;
			mxArray* tmp=readdata(fp, &header3, start2, 0);
			if(fp->eof){
				break;
			}
			if(tmp&&tmp!=SKIPPED){
				if(mxIsStruct(out)){
					/*
					  char *key=strstr(header3.str, "struct_key=");
					  char *key2=strstr(key, ";\n");
					  char field[64];
					  if(!key){
					  warning("key name is not found.\n");
					  snprintf(field, 64, "key%d", ix);
					  }else{
					  key+=11;
					  strncpy(field, 64, key, key2-key);
					  }
					  mxSetField(out, ix, field, tmp);*/
				} else{
					if(ntot>1||!skip_unicell){
						mxSetCell(out, ix, tmp);
					} else{
						out=tmp;
					}
				}
				if(header3){
					mxSetCell(header0, ix, header3);
					if(mxGetNumberOfElements(header3)){
						nheader0++;
					}
				}
			}
		}
		if(nheader0){
			if(header){
				mxSetCell(header0, ntot, *header);
				*header=header0;
			} else{
				mxDestroyArray(header0);
			}
		} else{
			mxDestroyArray(header0);//just use the global header2.
		}
	} else if(start==-1){//skip read.
		if(zfseek(fp, byte*ntot, SEEK_CUR)){
			error("Seek failed\n");
		}
		out=SKIPPED;
	} else if(magic!=M_COMMENT){
		if(issp){//sparse matrix
			int64_t nzmax;
			int64_t nx, ny;
			if(header2.ndim!=2){
				error("Invalid dims\n");
			}

			nx=header2.dims[0];
			ny=header2.dims[1];
			if(nx!=0&&ny!=0){
				zfread(&nzmax, sizeof(uint64_t), 1, fp);
			} else{
				nzmax=0;
			}
			if(fp->eof) return NULL;
			out=mxCreateSparse(nx, ny, nzmax, mxFLAG);
			if(nx!=0&&ny!=0&&nzmax!=0){
				if(sizeof(mwIndex)==ibyte){/*Match*/
					zfread(mxGetJc(out), ibyte, ny+1, fp);
					zfread(mxGetIr(out), ibyte, nzmax, fp);
				} else{//Convert index.
					long i;
					mwIndex* Jc0=mxGetJc(out);
					mwIndex* Ir0=mxGetIr(out);
					void* Jc=malloc(ibyte*(ny+1));
					void* Ir=malloc(ibyte*nzmax);
					zfread(Jc, ibyte, ny+1, fp);
					zfread(Ir, ibyte, nzmax, fp);
					if(ibyte==4){
						uint32_t* Jc2=(uint32_t*)Jc;
						uint32_t* Ir2=(uint32_t*)Ir;
						for(i=0; i<ny+1; i++){
							Jc0[i]=Jc2[i];
						}
						for(i=0; i<(long)nzmax; i++){
							Ir0[i]=Ir2[i];
						}
						free(Jc);
						free(Ir);
					} else if(ibyte==8){
						uint64_t* Jc2=(uint64_t*)Jc;
						uint64_t* Ir2=(uint64_t*)Ir;
						for(i=0; i<ny+1; i++){
							Jc0[i]=Jc2[i];
						}
						for(i=0; i<nzmax; i++){
							Ir0[i]=Ir2[i];
						}
						free(Jc);
						free(Ir);
					} else{
						mexErrMsgTxt("Invalid sparse format\n");
					}
				}
			}
			ntot=nzmax;
		} else{
			out=mxCreateNumericArray(header2.ndim, header2.dims, id, mxFLAG);
		}
#if !MX_HAS_INTERLEAVED_COMPLEX
	//convert from interleaved complex to separate complex.
		if(mxFLAG==mxCOMPLEX){
			void* tmp0=malloc(ntot*byte);
			zfread(tmp0, byte, ntot, fp);
			void* Pr0=mxGetPr(out);
			void* Pi0=mxGetPi(out);
			if(id==mxDOUBLE_CLASS){
				dcomplex* tmp=(dcomplex*)tmp0;
				double* Pr=(double*)Pr0;
				double* Pi=(double*)Pi0;
				for(long i=0; i<ntot; i++){
					Pr[i]=tmp[i].x;
					Pi[i]=tmp[i].y;
				}
			} else if(id==mxSINGLE_CLASS){
				fcomplex* tmp=(fcomplex*)tmp0;
				float* Pr=(float*)Pr0;
				float* Pi=(float*)Pi0;
				for(long i=0; i<ntot; i++){
					Pr[i]=tmp[i].x;
					Pi[i]=tmp[i].y;
				}
			} else{
				error("Invalid\n");
			}
			free(tmp0);
		} else
#endif
			zfread(mxGetData(out), byte, ntot, fp);

	} else{
		fprintf(stderr, "magic=%x\n", magic);
		warning("Unrecognized file. Please recompile the mex routines in the newest code\n");
		out=NULL;
	}
	start=start_save;
	if(!iscell&&fp->isfits==1){/*fits file may contain extra extensions.*/
		fp->isfits++;
		int icell=0;
		mxArray** outarr=NULL;
		mxArray** headerarr=NULL;
		while(out){
			icell++;
			outarr=(mxArray**)realloc(outarr, icell*sizeof(mxArray*));
			headerarr=(mxArray**)realloc(headerarr, icell*sizeof(mxArray*));
			if(out==SKIPPED) out=NULL;
			outarr[icell-1]=out;
			if(header) headerarr[icell-1]=*header;
			int start2=0;
			if(howmany!=0){//selective reading.
				if(icell<start||icell+1>start+howmany){
					start2=-1;//don't read next data.
				}
			}
			out=readdata(fp, header, start2, 0);
		}
		if(icell>1){/*set output.*/
			out=mxCreateCellMatrix(icell, 1);
			if(header) *header=mxCreateCellMatrix(icell, 1);
			int i;
			for(i=0; i<icell; i++){
				mxSetCell(out, i, outarr[i]);
				if(header) mxSetCell(*header, i, headerarr[i]);
			}
			free(outarr);
			free(headerarr);
		} else{
			out=outarr[0];
			if(header) *header=headerarr[0];
		}
	}
	if(!out){
		out=mxCreateDoubleMatrix(0, 0, mxREAL);
	}
	free(header2.dims);
	free(header2.str); header2.str=NULL;
	return out;
}
static char* mx2str(const mxArray* A){
	int nlen=mxGetNumberOfElements(A)+1;
	char* fn=(char*)malloc(nlen);
	mxGetString(A, fn, nlen);
	return fn;
}
void usage(){
	mexErrMsgTxt("Usage: [var, [header]]=read('filename' [,howmany] [, start]). [] means optional\n");
}
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	file_t* fp;
	char* fn=NULL;
	int start=0, howmany=0;
	switch(nrhs){
	case 3:
		start=(long)mxGetScalar(prhs[2]);//starting block to read. matlab index.
		//fallthrough
	case 2:
		howmany=(long)mxGetScalar(prhs[1]);
		//fallthrough
	case 1:
		fn=mx2str(prhs[0]);
		break;
	default:
		usage();
	}
	if(howmany>0){
		if(start>0){
			start--;//convert to C index.
		}
		if(start<0){
			start=0;
		}
	}
	fp=zfopen(fn, "rb");
	if(!fp){
		error2("Unable to open file: %s (%s)\n", fn, strerror(errno));
		return;
	}
	free(fn);
	switch(nlhs){
	case 0:
	case 1:
		plhs[0]=readdata(fp, NULL, start, howmany); break;
	case 2:
		plhs[0]=readdata(fp, &plhs[1], start, howmany); break;
	default:
		usage();
	}
	if(start==0&&howmany==0){
		int res=zfeof(fp);
		if(res){
			warning("There is unread data: res=%d\n", res);
		}
	}
	zfclose(fp);
}

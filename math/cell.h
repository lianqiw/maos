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
/**
   \file cell.h

   Function for a generic cell. use CELL() to case to cell.

*/
#ifndef AOS_LIB_CELL_H
#define AOS_LIB_CELL_H
#include "type.h"
cell* cellnew(long nx, long ny)CHECK_UNUSED_RESULT;
static inline cell* cell_cast(const void* A){
	return iscell(A)?(cell*)A:0;
}
cell *cellref(const_anyarray in);
uint32_t cellhash(const_anyarray A, uint32_t key) CHECK_UNUSED_RESULT;
cell *cellconvert(cell *A, cell* (*fun_convert)(cell*));
void cellinit(panyarray A, long nx, long ny);
void cellinit2(panyarray A, const_anyarray B);
void celldim(const_anyarray A_, long* nx, long* ny, long** nxs, long** nys);
void cellresize(anyarray in, long nx, long ny);
int cell_is_diag(const_anyarray A);
/*!free a cell array and zero the pointer.*/
#define cellfree(A) if(A){cellfree_do(A); A=NULL;}

void cellfree_do(anyarray dc);
void writedata(file_t* fd, const_anyarray A, long ncol);
void writebin(const_anyarray dc, const char* format, ...) CHECK_ARG(2);
cell *readdata(file_t *fp, M_ID id, header_t *header);
cell* readbin_id(M_ID id, int level, const char* format, ...) CHECK_ARG(3);
cell* readbin(const char* format, ...) CHECK_ARG(1);

//(const void* dc, const char* format, ...) CHECK_ARG(2);
void writecell_async(const_anyarray A, long ncol);
#define writebin_async(A,ncol) writecell_async(A, ncol)
void writebin_header(anyarray dc, const char* keywords, const char* format, ...) CHECK_ARG(3);
cell* readsock(int sock);
void writesock(const_anyarray dc, int sock);
#endif

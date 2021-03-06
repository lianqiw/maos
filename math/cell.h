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
/**
   \file cell.h

   Function for a generic cell. The data those functions are passed
   as void * to accomodate any cell type.

*/
#ifndef AOS_LIB_CELL_H
#define AOS_LIB_CELL_H
#include "type.h"
cell* cellnew(long nx, long ny)CHECK_UNUSED_RESULT;
cell* cellnew2(const void*)CHECK_UNUSED_RESULT;
static inline cell* cell_cast(const void* A){
	return iscell(A)?(cell*)A:0;
}
void cellinit(cell** A, long nx, long ny);
void cellinit2(cell** A, const cell* B);
void celldim(const cell* A_, long* nx, long* ny, long** nxs, long** nys);
void cellresize(void* in, long nx, long ny);
static inline void cellreshape(void* in_, long nx, long ny){
	cell* in=cell_cast(in_);
	if(in->nx*in->ny!=nx*ny){
		error("Incorrect new dimension\n");
	} else{
		in->nx=nx;
		in->ny=ny;
	}
}
/*!free a cell array and zero the pointer.*/
#define cellfree(A) ({cellfree_do(A); A=NULL;})

void cellfree_do(void* dc);
void writedata_by_id(file_t* fd, const void* pix, uint32_t id);
void write_by_id(const void* dc, uint32_t id, const char* format, ...) CHECK_ARG(3);

cell* readdata_by_id(file_t* fp, uint32_t id, int level, header_t* header);
cell* read_by_id(uint32_t id, int level, const char* format, ...) CHECK_ARG(3);

cell* readbin(const char* format, ...) CHECK_ARG(1);
void writebin(const void* dc, const char* format, ...) CHECK_ARG(2);
void writebin_auto(const void* A);
void writebin2(void* dc, const char* header, const char* format, ...) CHECK_ARG(3);
cell* readsock(int sock);
void writesock(const void* dc, int sock);
#endif

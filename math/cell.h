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
/**
   \file cell.h 

   Function for a generic cell. The data those functions are passed
   as void * to accomodate any cell type.
   
*/
#ifndef AOS_LIB_CELL_H
#define AOS_LIB_CELL_H
#include "type.h"
void* cellnew(long nx, long ny)CHECK_UNUSED_RESULT;
cell *cell_cast(const void *A) CHECK_UNUSED_RESULT;
void cellinit(cell **A, long nx, long ny);

void cellresize(void *in, long nx, long ny);
void cellfree_do(void* dc);
void writedata_by_id(file_t *fd, const void* pix, long id);
void write_by_id(const void* dc, long id, const char* format,...);

void* readdata_by_id(file_t *fp, long id, int level, header_t *header);
void* read_by_id(long id, int level, const char *format, ...);

void* readbin(const char *format, ...);
void writebin(const void *dc, const char *format, ...);
void writebindata(file_t *fp, const void *dc);
#endif

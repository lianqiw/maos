/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#ifndef AOS_LIB_SHM_H
#define AOS_LIB_SHM_H
#if USE_POSIX_SHM == 1
#include <sys/mman.h>
#include <sys/file.h>
#include <unistd.h>
void shm_unmap(void *p, int fd);
void shm_free_older(long sec);
int  shm_free_unused(char *fnshm, int timeout);
long shm_get_avail(void);
extern int shm_keep_unused;
#endif
#endif

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
#ifndef AOS_SYS_ZSTD_H
#define AOS_SYS_ZSTD_H
#include <stdio.h>
#include <stddef.h>
/*
\file zstdio.h
	wrapper for zstd compression/decompression.
*/
#ifdef __cplusplus
extern "C" {
#endif

typedef struct zstd_t zstd_t;

zstd_t* zstd_open(int fd, int mode); // mode = 'r': read, 'w': write
size_t zstd_read(zstd_t *zf, void *buf, size_t size);
int zstd_write(zstd_t *zf, const void *buf, size_t size);
int zstd_flush(zstd_t *zf);
int zstd_close(zstd_t *zf, int close_fd);

#ifdef __cplusplus
}
#endif

#endif

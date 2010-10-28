/*
  Copyright 2009,2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#ifndef AOS_SOCKIO_H
#define AOS_SOCKIO_H
int swrite(int *fd, const void *p, size_t len);
int sread(int *fd, void *p, size_t len);
void swriteint(int *fd, int cmd);
void swriteintarr(int *fd, int* cmd, unsigned int len);
int sreadint(int *fd);
void sreadintarr(int *fd, int* cmd, unsigned int len);
void swritestr(int *fd, const char *str);
char *sreadstr(int *fd);
#endif

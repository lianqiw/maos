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
#ifndef AOS_IO_H
#define AOS_IO_H
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

void writestr(int fd, const char *str);
char * readstr(int fd);
void fwritestr(FILE *fp, const char *str);
char* freadstr(FILE *fp);

void writeint(int fd, int cmd);
int readint(int fd);
void fwriteint(FILE *fp, int cmd);
int   freadint(FILE *fp);
#endif

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

#ifndef AOS_LIB_MISC_H
#define AOS_LIB_MISC_H
//#include "common.h"
#include "sys/misc.h"
#undef strdup
#define strdup mystrdup //our strdup handles NULL correctly.
#undef strndup
#define strndup mystrndup //our strdup handles NULL correctly.

char *FF(const char *format,...) CHECK_ARG(1);
long factorial(long n);
char *mybasename(const char *fn);
int check_suffix(const char *fn, const char *suffix);
void copyfile(const char *dest, const char *src);
char *argv2str(int argc, char **argv);
void print_file(const char *fnin);
#endif


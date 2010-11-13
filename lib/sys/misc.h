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
#ifndef AOS_SYS_MISC_H
#define AOS_SYS_MISC_H
#include "common.h"
int myclocki(void);
double myclockd(void);
const char *myasctime(void);
char *strtime(void);
const char *myhostname(void);
char *mygetcwd(void);
char *myabspath(const char *path);
void mysymlink(const char *fn, const char *fnlink);
int exist(const char *fn);/*test file exist*/
void touch(const char *fn);
char *stradd(const char* a, ...) CHECK_NULL_TERMINATED;
void expand_filename(char **fnout, const char *fn);
char *mystrndup(const char *A, int len);
char *mystrdup(const char *A);
void maxmindbl(const double *restrict p, long N, 
	       double *restrict max, double *restrict min, double *restrict sum);
void maxminlong(const long *restrict p, long N,
		long *restrict max, long *restrict min, long *restrict sum);
void maxmincmp(const dcomplex *restrict p, long N,
	       double *restrict max, double *restrict min, double *restrict sum);
void remove_file_older(const char *fndir, long sec);
void mymkdir(const char *format,...) CHECK_ARG(1);
int mystrcmp(const char *a, const char *b);
#endif

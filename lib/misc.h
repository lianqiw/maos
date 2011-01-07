/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "common.h"
long factorial(long n);
char *mybasename(const char *fn);
int check_suffix(const char *fn, const char *suffix);
void copyfile(const char *dest, const char *src);
char *argv2str(int argc, char **argv);
void print_file(const char *fnin);
int myclocki(void);
double myclockd(void);
const char *myasctime(void);
char *strtime(void);
const char *myhostname(void);
char *mygetcwd(void);
char *myabspath(const char *path);
void mysymlink(const char *fn, const char *fnlink);
int exist(const char *fn);/*test file exist*/
int isdir(const char *fn);
int isfile(const char *fn);
void touch(const char *fn);
char *stradd(const char* a, ...) CHECK_NULL_TERMINATED;
void expand_filename(char **fnout, const char *fn);

void remove_file_older(const char *fndir, long sec);
void mymkdir(const char *format,...) CHECK_ARG(1);
int mystrcmp(const char *a, const char *b);
char *mystrndup(const char *A, int len);
void cloexec(int fd);
void mysleep(double sec);
#if USE_MEM == 1
char *mystrdup(const char *A);
#undef strdup
#define strdup mystrdup //our strdup handles NULL correctly.
#endif //USE_MEM=1
#endif


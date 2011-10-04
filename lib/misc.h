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
double factorial(long n);
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
int exist(const char *fn);
int isdir(const char *fn);
int isfile(const char *fn);
int islink(const char*fn);
size_t flen(const char *fn);
time_t fmtime(const char *fn);
char *stradd(const char* a, ...) CHECK_NULL_TERMINATED;
char *strnadd(int argc, char **argv, const char *delim);
void expand_filename(char **fnout, const char *fn);

void remove_file_older(const char *fndir, long sec);
void mymkdir(const char *format,...) CHECK_ARG(1);
int mystrcmp(const char *a, const char *b);
char *mystrndup(const char *A, int len);
void cloexec(int fd);
void mysleep(double sec);
void mypause(void);
long available(const char *path);
/*strdup for static pointers or in constructor should be using this strdup0 to
  avoid reporting of unfreed mem.*/
extern char* (*strdup0)(const char *);
#if USE_MEM == 1
char *mystrdup(const char *A);
#undef strdup
#define strdup mystrdup /*our strdup handles NULL correctly, and talk to mem.c */
#endif /*USE_MEM=1 */

typedef struct ARGOPT_T{
    char *name; /**<The long name*/
    char key;   /**<The short name*/
    int type;   /**<The type of result expected*/
    int opt;    /**<0: the key does not need value, like -d. 
		   1:The key needs value, like -s 1.
		   2:The key does not need value, and the supplied val is a function pointer.
		   3:The key needs value, and the supplied val is a function pointer.
		*/
    void *val;  /**<The the address to put the return result.*/
    int *nval;  /**<If val is array, this is the counter.*/
}ARGOPT_T;
char *parse_argopt(int argc, char **argv, ARGOPT_T *options);
#endif

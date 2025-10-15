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
#ifndef AOS_SYS_MISC_H
#define AOS_SYS_MISC_H
#include <time.h>
#include "common.h"
/**
   \file sys/misc.h
   misc routines.
*/
//const char *mybasename(const char *fn);
char *mydirname(const char *fn);
int check_suffix(const char *fn, const char *suffix);
int copyfile(const char* src, const char* dest);
char *argv2str(int argc, const char *argv[], const char *delim);
void print_file(const char *fnin);
time_t myclocki(void);
double myclockd(void);
const char *myasctime(time_t at);
char *strtime_pid(void);
const char *myhostname(void);
char *mygetcwd(void);
char *myabspath(const char *path);
int mysymlink(const char *source, const char *dest);
int mylink(const char* source, const char* dest);
int exist(const char *fn);
void touch(const char* format, ...) CHECK_ARG(1);
int isdir(const char *fn);
int isfile(const char *fn);
int islink(const char*fn);
int issock(int fd);
size_t flen(const char *fn);
time_t fmtime(const char *fn);
char *stradd(const char* a, ...) CHECK_NULL_TERMINATED;
char *strnadd(int argc, const char *argv[], const char *delim);
char *expand_filename(const char *fn);

void remove_file_older(const char *fndir, int level, long sec);
void mymkdir(const char *format,...) CHECK_ARG(1);
int mystrcmp(const char *a, const char *b);
char *mystrndup(const char *A, size_t len);
void cloexec(int fd);
void mysleep(double sec);
int mypause(int fd1, int fd2);
long available_space(const char* path);
/*strdup for static pointers or in constructor should be using this strdup0 to
  avoid reporting of unfreed mem.*/
extern char* (*strdup0)(const char *);
char *mystrdup(const char *A);
int mysnprintf(char *restrict str, size_t size, const char *restrict format, ...) CHECK_ARG(3);
#undef strdup
#define strdup mystrdup /*our strdup handles NULL correctly, and talk to mem.c */
#undef strndup
#define strndup mystrndup
#undef snprintf
#define snprintf mysnprintf /*our snprintf avoids the truncation error by checking the returned size. */
void mystrrep(char*str, const char *prefix, const char *substitute);
typedef struct argopt_t{
    const char *name;    /**<The long name*/
    char key;      /**<The short name*/
    int type;      /**<The type of result expected*/
    int valtype;   /**<The type of input expected:
		      0: The key does not need value;
		      1: The key needs value;
		      2: The key accepts an array of values.*/
    int isfun;  /**< Whether val indicate a function pointer
		      0: The address to put the return results.
		      1: A function pointer to call. */
    void *val;  /**<The address to put the return result.*/
    int *nval;  /**<If val is array, this is the counter.*/
}argopt_t;
void parse_argopt(char *cmds, argopt_t *options);
//sem_lock is removed in favor of lock_file which can recover when a process exits.
//int sig_block(int block);
//int sem_lock(const char *key, int lock);
void set_realtime(int icpu, int niceness);
void free_strarr(char **str, int n);
extern const int default_color_table[];
#define default_color(i) default_color_table[i%11]
void print_version(void);
int sec2str(char*tmp, long stmp, double sec);
void rename_log(int sig, const char *exe);
#endif

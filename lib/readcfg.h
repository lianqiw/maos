/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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


#ifndef __READ_CONFIG_H
#define __READ_CONFIG_H
/**
   \file readcfg.h

   Routines to read .conf type config files. Each entry is composed of a key and
   a value. The entries are maintained in a hash table. Each entry can be
   retrieved from the key.
*/
#include "../math/mathdef.h"
#define format2key				\
    char key[512];				\
    va_list ap;					\
    va_start(ap,format);			\
    vsnprintf(key,sizeof(key), format, ap);	\
    va_end(ap)

void   open_config (const char*fn, const char *prefix, int priority);
void   close_config(const char*format,...) CHECK_ARG(1);

int    readcfg_peek(const char*format,...) CHECK_ARG(1);
int    readcfg_peek_n(const char *format, ...) CHECK_ARG(1);
int    readcfg_peek_override(const char *format,...) CHECK_ARG(1);
void   readcfg_ignore(const char *format, ...) CHECK_ARG(1);

char*  readcfg_str (const char*format,...) CHECK_ARG(1);
int    readcfg_int (const char*format,...) CHECK_ARG(1);
real   readcfg_dbl (const char*format,...) CHECK_ARG(1);

int    readcfg_strarr(char ***res, int len, int relax, const char *format, ...) CHECK_ARG(4);
int    readcfg_intarr(int **ret, int len, int relax, const char *format,...) CHECK_ARG(4);
int    readcfg_dblarr(real **ret,int len, int relax, const char *format,...) CHECK_ARG(4);

dmat*  readstr_dmat(int n, int relax, const char* key, const char *str);
dmat*  readcfg_dmat(int n, int relax, const char *format, ...) CHECK_ARG(3);
lmat*  readcfg_lmat(int n, int relax, const char *format,...) CHECK_ARG(3);

dcell* readcfg_dcell(const char *format,...) CHECK_ARG(1);
#endif

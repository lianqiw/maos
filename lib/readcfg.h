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

#ifndef __READ_CONFIG_H
#define __READ_CONFIG_H

#include "misc.h"
#define format2key				\
    char key[512];				\
    va_list ap;					\
    va_start(ap,format);			\
    vsnprintf(key,sizeof(key), format, ap);	\
    va_end(ap)
/**
   \file readcfg.h

   Routines to read .conf type config files. Each entry is composed of a key and
   a value. The entries are maintained a hash table. Each entry can be
   retrieved from the key.
*/
void   open_config (const char*fn, long protect);
void   close_config(const char*format,...) CHECK_ARG(1);

int    readcfg_peek(const char*format,...) CHECK_ARG(1);
int    readcfg_override(const char *format,...) CHECK_ARG(1);
char*  readcfg_str (const char*format,...) CHECK_ARG(1);
/**
   Read an int number from config hash table*/
int    readcfg_int (const char*format,...) CHECK_ARG(1);
/**
   Read a double number from config hash table
*/
double readcfg_dbl (const char*format,...) CHECK_ARG(1);
/**
   Read an int array from config hash table
*/
int    readcfg_intarr(int **ret,const char*format,...) CHECK_ARG(2);
/**
   Read a double array from config hash table
*/
int    readcfg_dblarr(double **ret,const char*format,...) CHECK_ARG(2);
int    readcfg_strarr(char ***res, const char *format,...) CHECK_ARG(2);
#endif

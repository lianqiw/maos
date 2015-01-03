/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "common.h"
/**
   \file readstr.h

   Group all routines that are used to parse values from string that contain
   key=value pairs.
*/
/*Check whether the char is space. We only treat real space, \t and \v as space. */
INLINE int is_space(char c){
    if(c==' ' || c=='\t' || c=='\v'){
	return 1;
    }else{
	return 0;
    }
}
/*Check whether the current express ends. String end or newline is treated as end.*/
INLINE int is_end(char c){
    if(c=='\0' || c=='\n' || c=='\r'){
	return 1;
    }else{
	return 0;
    }
}
int readstr_strarr(char ***res, int len, const char *sdata);
double readstr_num(const char *data, char **endptr0);
int readstr_numarr(void **ret, int len, int *nrow0, int *ncol0, int type, const char *data);
int readstr_intarr(int**ret, int len, const char *data);
void readstr_intarr_nmax(int **ret, int len, const char *data);
void readstr_intarr_relax(int **ret, int len, const char *data);

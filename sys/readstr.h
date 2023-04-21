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

#include "common.h"
#include "bin.h"
/**
   \file readstr.h

   Group all routines that are used to parse values from string that contain
   key=value pairs.
*/
/*Check whether the char is space. We only treat real space, \t and \v as space. */
/*
static inline int is_space(char c){
    if(c==' ' || c=='\t' || c=='\v'){
	return 1;
    }else{
	return 0;
    }
}
*/
/*Check whether the current express ends. String end or newline is treated as end.*/
static inline int is_end(char c){
    if(c=='\0' || c=='\n' || c=='\r'){
	return 1;
    }else{
	return 0;
    }
}
int readstr_strarr(char ***res, int len, int relax, const char *key, const char *sdata);
double readstr_num(const char *key, const char *data, char **endptr0);
int readstr_numarr(void **ret, int *nrow0, int *ncol0, int len, int relax, int type, const char *key, const char *data);
void trim_string(const char** pheader, const char** pend);
const char* search_keyword(const char* keywords, const char* key);
double search_keyword_num(const char* keywords, const char *key);
double search_keyword_num_valid(const char* keywords, const char *key);
double search_keyword_num_default(const char* keywords, const char *key, double value0);

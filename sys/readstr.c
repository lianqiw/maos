/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#include "readstr.h"
#include "misc.h"

/**
   Group all routines that are used to parse values from string that contain
   key=value pairs.
*/

/**
   Obtain a string array value from the key.
 */
int readstr_strarr(char ***res, /**<[out] Result*/
		   int len,     /**<[out] Length of array*/
		   const char *sdata /**<[in] Input string*/
		   ){
    int count=0;
    int maxcount=5;
    if(len && *res){
	memset(*res, 0, sizeof(char*)*len);
    }else{
	if(len) maxcount=len;
	*res=calloc(maxcount,sizeof(char*));
    }
    const char *sdataend=sdata+strlen(sdata)-1;
    const char *sdata2;
    if(sdata[0]!='[' || sdataend[0]!=']'){
	/*spaces are already trimmed out.*/
	warning2("{%s}: Should start with [ and end with ]\n", sdata);
	sdata2=sdata;
    }else{
	sdata2=sdata+1;
    }
    /*skip spaces*/
    while(sdata2<sdataend && sdata2[0]==' '){
	sdata2++;
    }
    while(sdata2<sdataend){
	char mark=' ';
	if(sdata2[0]=='"' || sdata2[0]=='\''){
	    mark=sdata2[0];
	    sdata2++;
	}
	const char *sdata4=strchr(sdata2, mark);
	if(!sdata4){
	    if(mark!=' '){
		error("{%s}: Unmatched string\n", sdata);
	    }else{
		sdata4=sdataend;
	    }
	}
	if(sdata4>sdata2 && sdata4<=sdataend){/*found non-empty str*/
	    if(!len && count>=maxcount){
		maxcount*=2;
		*res=realloc(*res,sizeof(char*)*maxcount);
	    }
	    (*res)[count]=mystrndup(sdata2, sdata4-sdata2);
	}else{/*found empty str*/
	    (*res)[count]=NULL;
	}
	count++;
	sdata2=sdata4+1;
	while(sdata2<sdataend && sdata2[0]==' '){
	    sdata2++;
	}
    }
    if(!len){
	if(count>0){
	    *res=realloc(*res,sizeof(char*)*count);
	}else{
	    free(*res); *res=NULL;
	}
    }
    return count;
}

/**
   Read in a number from the value string. Will interpret * and / operators. *endptr0 will
   be updated to point to the next valid entry, or at separator like coma
   (spaced are skipped).  */
double readstr_num(const char *data, /**<[in] Input string*/
		   char **endptr0    /**<[out] Location in Input string after readed number.*/
		   ){
    if(!data || strlen(data)==0){
	error("{%s}: Unable to parse for a number\n", data);
	return NAN;
    }
    char *endptr;
    double res=strtod(data, &endptr);
    if(data==endptr){
	error("{%s}: Unable to parse for a number\n", data);
	return NAN;
    }
    while(is_space(endptr[0])) endptr++;
    while(endptr[0]=='/' || endptr[0]=='*'){
	int power=1;
	if(endptr[0]=='/'){
	    power=-1;
	}
	endptr++;
	while(is_space(endptr[0])) endptr++;
	data=endptr;
	double tmp=strtod(data, &endptr);
	if(data==endptr){
	    error("{%s}: Failed to parse for a number\n", data);
	}
	if(power==1){
	    res*=tmp;
	}else{
	    res/=tmp;
	}
	while(is_space(endptr[0])) endptr++;
    }
    if(endptr0){
	*endptr0=endptr;
    }
    return res;
}


/**
   Read numerical array from a string. if len is nonzero, *ret should be already allocated. 
   Can read in the following formats:
   [1 2 3]   as 3 col, 1 row [1 2 3] in memory
   [1 2 3]+2 as [3 4 5] in memory
   2*[1 2 3] or [1 2 3]*2 as [2 4 6] in memory
   [1 2 3]/2 as [0.5 1 1.5] in memory
   2/[1 2 3] as [2 1 2/3] in memory
   2/[1 2 3]*2+1 as [5 3 7/3] in memory

   2-d arrays are groups as columns (in contrast to matlab that uses rows); semi-colon is used to separate columns.
   [1 2 3; 4 5 6] is 2 d array, 3 columns, 2 rows. In memory is stored as [1 2 3 4 5 6]
   transpose is supported:
   [1 2 3; 4 5 6]' is 2 d array, 2 columns, 3 rows. In memory is stored as [1 4 2 5 3 6]

   \return Number of values actually read.
*/
int readstr_numarr(void **ret, /**<[out] Result*/
		   int len,    /**<[in]  Max number of values to read.*/
		   int *nrow0, /**<[out] Number of rows (nx)*/
		   int *ncol0, /**<[out] Number of columns (ny)*/
		   int type,   /**<[in]  Data type*/
		   const char *data /**<[in] Input string*/
		   ){
    if(!data || strlen(data)==0){
	return 0;
    }
    size_t nmax=10;
    size_t size=0;/*size of each number */
    switch(type){
    case T_INT:
	size=sizeof(int);
	break;
    case T_DBL:
	size=sizeof(double);
	break;
    default:
	error("Invalid type");
    }
    if(len==0){
	if(!(*ret=malloc(size*nmax))){
	    error("Failed to allocate memory for ret\n");
	}
    }else{
	nmax=len;
	if(*ret){
	    memset(*ret, 0, size*len);
	}else{
	    *ret=calloc(size, len);
	}
    }
    const char *startptr=data;
    char *endptr, *startptr2;
    double fact=1;
    int power=1;
    int addition=0; double addval=0; 
    int trans=0;
    /*process possible numbers before the array. */
    if(strchr(startptr,'[')){/*there is indeed '[' */
	while(startptr[0]!='['){/*parse number before [*/
	    double fact1=strtod(startptr, &endptr);/*get the number */
	    if(startptr==endptr){
		error("{%s}: Invalid entry to parse for numerical array\n", data);
	    }else{/*valid number */
		if(power==1){
		    fact*=fact1;
		}else{
		    fact/=fact1;
		}
		while(is_space(endptr[0])) endptr++;
		if(endptr[0]=='/'){
		    power=-1;
		}else if(endptr[0]=='*'){
		    power=1;
		}else{
		    error("{%s}: Invalid entry to parse for numerical array.\n", data);
		}
		startptr=endptr+1;
	    }
	}
	if(startptr[0]!='['){
	    error("{%s}: Invalid entry to parse for numerical array\n", data);
	}
	startptr++;/*points to the beginning of the array in [] */
	/*process possible numbers after the array. do not use startptr here.*/
    }else{
	/*warning2("Expecting array: {%s} should start with [\n", data); */
    }
    if(strchr(startptr,']')){/*there is indeed ']' */
	endptr=strchr(startptr,']')+1;
	while(is_space(endptr[0]) || endptr[0]=='\''){
	    if(endptr[0]=='\'') trans=1-trans;
	    endptr++;
	}
	while(endptr[0]=='/' || endptr[0]=='*' || endptr[0]=='+' || endptr[0]=='-'){
	    int power2=1;
	    if(endptr[0]=='/'){
		power2=-1;
	    }else if(endptr[0]=='*'){
		power2=1;
	    }else if(endptr[0]=='+' || endptr[0]=='-'){
		power2=0;
		if(endptr[0]=='+'){
		    addition=1;
		}else{
		    addition=-1;
		}
	    }
	    if(addition && power2){
		error("We do not yet support * or / after + or - \n");
	    }
	    endptr++;
	    while(is_space(endptr[0])) endptr++;
	    startptr2=endptr;
	    double fact2=strtod(startptr2, &endptr);
	    if(startptr2==endptr){
		error("{%s}: Invalid entry to parse for numerical array\n", data);
	    }
	    while(is_space(endptr[0])) endptr++;
	    if(addition){/*addition*/
		if(addition==1){
		    addval+=fact2;
		}else{
		    addval-=fact2;
		}
	    }else{
		if(power2==1){
		    fact*=fact2;
		}else{
		    fact/=fact2;
		}
	    }
	}
	if(!is_end(endptr[0])){
	    error("{%s}: There is garbage in the end of the string.\n", data);
	}
    }
    int count=0;
    int nrow=0;/*number of rows*/ 
    int ncol=0;/*number of columns*/
    int rowbegin=0;/*beginning of this row*/
    /*Read in the array */
    while(startptr[0]!=']' && !is_end(startptr[0])){
	if(count>=nmax){
	    if(len){
		error("{%s}: Needs %d numbers, but more are supplied.\n", data, len);
	    }else{
		nmax*=2;
		*ret=realloc(*ret, size*nmax);
	    }
	}
	/*parse the string for a floating point number.  */
	double res=readstr_num(startptr, &endptr);
	startptr=endptr;
	/*apply the factors appear before or after [] */
	if(power==1){
	    res=fact*res;
	}else{
	    res=fact/res;
	}
	if(addition){
	    res+=addval;
	}
	/*assign the value to appropriate array. convert to int if necessary. */
	switch(type){
	case T_INT:
	    ((int*)(*ret))[count]=(int)res;
	    break;
	case T_DBL:
	    ((double*)(*ret))[count]=res;
	    break;
	default:
	    error("Invalid type");
	}
	count++;
	/*Skip the number separators. */
	while(startptr[0]==' '||startptr[0]==';'){
	    if(startptr[0]==';'){
		ncol++;
		if(nrow==0){
		    nrow=count-rowbegin;
		}else if(nrow!=count-rowbegin){
		    error("{%s}: last row has %d numbers while new row has %d numbers\n", data, nrow, count-rowbegin);
		}
		rowbegin=count;
	    }
	    startptr++;
	}
    }
    if(rowbegin<count){/*if last row is not ended with ;*/
	ncol++;
	if(nrow==0){
	    nrow=count-rowbegin;
	}else if(nrow!=count-rowbegin){
	    error("{%s}: last row has %d numbers while new row has %d numbers\n", data, nrow, count-rowbegin);
	}
    }
    if(nrow*ncol!=count){
	error("nrow=%d, ncol=%d, count=%d\n", nrow, ncol, count);
    }

    if(trans && count>0){
	info("Transposing %dx%d array\n", ncol, nrow);
	void *newer=malloc(size*count);
	switch(type){
	case T_INT:{
	    int *from=(*ret);
	    int *to=newer;
	    for(int icol=0; icol<ncol; icol++){
		for(int irow=0; irow<ncol; irow++){
		    to[icol+ncol*irow]=from[irow+nrow*icol];
		}
	    }
	}
	    break;
	case T_DBL:{
	    double *from=(*ret);
	    double *to=newer;
	    for(int icol=0; icol<ncol; icol++){
		for(int irow=0; irow<ncol; irow++){
		    to[icol+ncol*irow]=from[irow+nrow*icol];
		}
	    }
	}
	    break;
	default:
	    error("Invalid type");
	}
	int tmp=ncol; ncol=nrow; nrow=tmp;
	free(*ret);
	*ret=newer;
    }else if(!len){
	if(count>0){
	    *ret=realloc(*ret, size*count);
	}else{
	    free(*ret);
	    *ret=NULL;
	}
    }
    if(nrow0) *nrow0=nrow;
    if(ncol0) *ncol0=ncol;
    return count;
}
/**
   Read an integer array. Duplicate if only one number is present.
 */
void readstr_intarr_nmax(int **ret, /**<[out] Result*/
			 int len,   /**<[in]  Max number of values to read.*/
			 const char *data /**<[in] Input string*/
			 ){
    int len2=readstr_numarr((void**)ret, len,NULL,NULL, T_INT, data);
    if(len2==1){
	for(int i=1; i<len; i++){
	    (*ret)[i]=(*ret)[0];
	}
    }else if(len2!=0 && len2!=len){
	error("{%s}: Require %d numbers, but got %d\n", data, len, len2);
    }
}

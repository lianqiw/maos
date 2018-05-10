/*
  Copyright 2009-2018 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include <tgmath.h>
#include <ctype.h> /*isspace */
#include "readstr.h"
#include "misc.h"
#include "bin.h"

/**
   Group all routines that are used to parse values from string that contain
   key=value pairs.
*/

/**
   Obtain a string array value from the key. String entries must be separated by
   space, coma [,] or semi-column [;]. Enclosing with quote is optional, but
   necessary to protect spaces.
   Examples:
   []: one empty string
   [,]: two empty strings
   [a b]: two strings: "a" and "b";
   [a,b]: same as [a b]
   ['a b']: one string "a b".
 */
int readstr_strarr(char ***res, /**<[out] Result*/
		   int len,     /**<[in] Length of array*/
		   const char *sdata /**<[in] Input string*/
		   ){
    int count=0;
    int maxcount=len;
    if(len && *res){//pre allocatrd data.
	memset(*res, 0, sizeof(char*)*len);
    }else{
	if(!len) maxcount=5;
	*res=mycalloc(maxcount,char*);
    }
    if(!sdata) return count;
    const char *sdataend=sdata+strlen(sdata)-1;
    const char *sdata2=sdata;
    if(sdata[0]=='['){
	sdata2++;
	if(sdataend[0]!=']'){
	    error("{%s}: Does not end in ].\n", sdata);
	}
    }
    /*skip spaces*/
    while(sdata2<sdataend && sdata2[0]==' '){
	sdata2++;
    }
    int end_coma=0;//end with coma. append an additional element.
    while(sdata2<sdataend || end_coma){
	const char *sdata4=sdataend;
	const char *sdata3;
	end_coma=0;
	//If entry is quoted
	if(sdata2[0]=='"' || sdata2[0]=='\''){
	    char sep=sdata2[0];
	    sdata2++;
	    sdata4=strchr(sdata2, sep);//find matching quote.
	    if(!sdata4){
		error("{%s}: Unmatched quote in input\n", sdata);
	    }
	    sdata3=sdata4+1;
	    //Skip spaces
	    while(sdata3<sdataend && sdata3[0]==' '){
		sdata3++;
	    }
	    //Ignore separator following the quote.
	    if(sdata3[0]==',' || sdata3[0]==';'){
		sdata3++;
		end_coma=1;
	    }
	}else{
	    //Find the next space, coma or semi-colon.
	    char sep[]=" ,;";
	    for(size_t is=0; is<sizeof(sep); is++){
		const char *tmp=strchr(sdata2, sep[is]);
		if(tmp && tmp<sdata4 && tmp<sdataend){
		    sdata4=tmp;
		}
	    }
	    if(sdata4[0]==',' || sdata4[0]==';'){
		end_coma=1;
	    }
	    //skip the separator
	    sdata3=sdata4+1;
	}
	if(count>=maxcount){//check memory
	    if(!len){
		maxcount*=2;
		*res=myrealloc(*res,maxcount,char*);
	    }else{
		error("{%s}: need %d, got more than %d elements\n", sdata, len, count);
	    }
	}
	if(sdata4>sdata2){/*found non-empty str*/
	    (*res)[count++]=mystrndup(sdata2, sdata4-sdata2);
	}else{/*found empty str*/
	    (*res)[count++]=NULL;
	}
	sdata2=sdata3;
	while(sdata2<sdataend && (sdata2[0]==' ')){//skip space
	    sdata2++;
	}
    }
    if(!len){
	if(count>0){
	    *res=myrealloc(*res,count,char*);
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
    while(isspace(endptr[0])) endptr++;
    while(endptr[0]=='/' || endptr[0]=='*'){
	int power=1;
	if(endptr[0]=='/'){
	    power=-1;
	}
	endptr++;
	while(isspace(endptr[0])) endptr++;
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
	while(isspace(endptr[0])) endptr++;
    }
    if(endptr0){
	*endptr0=endptr;
    }
    return res;
}


/**
   Read numerical array from a string. if len is nonzero, *ret should be already allocated. 
   NOTICe that continuous numbers are readin as column vector, not row vector(as in matlab)
   Can read in the following formats:
   [1 2 3]   as 3 rows, 1 col [1 2 3] in memory
   [1 2 3]+2 as [3 4 5] in memory
   2*[1 2 3] or [1 2 3]*2 as [2 4 6] in memory
   [1 2 3]/2 as [0.5 1 1.5] in memory
   2/[1 2 3] as [2/1 2/2 2/3] in memory
   2/[1 2 3]*2+1 as [5 3 7/3] in memory

   2-d arrays are groups as column vecotrs (in contrast to matlab that uses rows); semi-colon is used to separate columns.
   [1 2 3; 4 5 6] is 2 d array, 3 rows, 2 columns. In memory is stored as [1 2 3 4 5 6]
   transpose is supported:
   [1 2 3; 4 5 6]' is 2 d array, 2 rows, 3 columns. In memory is stored as [1 4 2 5 3 6]

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
    case M_INT:
	size=sizeof(int);
	break;
    case M_LONG:
	size=sizeof(long);
	break;
    case M_DBL:
	size=sizeof(double);
	break;
    default:
	error("Invalid type");
    }
    if(len==0){
	if(!(*ret=calloc(nmax, size))){
	    error("Failed to allocate memory for ret\n");
	}
    }else{
	nmax=len;
	if(*ret){
	    memset(*ret, 0, size*len);
	}else{
	    *ret=calloc(len, size);
	}
    }
    const char *startptr=data;
    const char *endptr, *startptr2;
    double fact=1;
    int power=1;
    int addition=0; double addval=0; 
    int trans=0;
    /*process possible numbers before the array. */
    if(strchr(startptr,'[')){/*there is indeed '[' */
	while(startptr[0]!='['){/*parse number before [*/
	    double fact1=strtod(startptr, (char**)&endptr);/*get the number */
	    if(startptr==endptr){
		error("{%s}: Invalid entry to parse for numerical array\n", data);
	    }else{/*valid number */
		if(power==1){
		    fact*=fact1;
		}else{
		    fact/=fact1;
		}
		while(isspace(endptr[0])) endptr++;
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
	/*warning("Expecting array: {%s} should start with [\n", data); */
    }
    if(strchr(startptr,']')){/*there is indeed ']'. Handle operations after ] */
	endptr=strchr(startptr,']')+1;
	while(isspace(endptr[0]) || endptr[0]=='\''){
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
		error("{%s}: We do not yet support * or / after + or - \n", data);
	    }
	    endptr++;
	    while(isspace(endptr[0])) endptr++;
	    startptr2=endptr;
	    double fact2=strtod(startptr2, (char**)&endptr);
	    if(startptr2==endptr){
		error("{%s}: Invalid entry to parse for numerical array\n", data);
	    }
	    while(isspace(endptr[0])) endptr++;
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
	if(!is_end(endptr[0]) && endptr[0]!=';'){
	    error("{%s}: There is garbage in the end of the string.\n", data);
	    _Exit(1);
	}
    }
    size_t count=0;
    size_t nrow=0;/*number of rows*/ 
    size_t ncol=0;/*number of columns*/
    size_t rowbegin=0;/*beginning of this row*/
    /*Read in the array */
    while(startptr[0]!=']' && !is_end(startptr[0])){
	if(count>=nmax){
	    if(len){
		error("{%s}: Needs %d numbers, but more are supplied.\n", data, len);
	    }else{
		*ret=myrealloc(*ret, size*nmax*2,char);
		memset((char*)*ret+size*nmax, 0, size*nmax);
		nmax*=2;
	    }
	}
	/*parse the string for a floating point number.  */
	double res=readstr_num(startptr, (char**)&endptr);
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
	case M_INT:
	    if(fabs(res-(int)res)>EPS){
		warning("Floating point number supplied while integer is needed: {%s}\n", data);
	    }
	    ((int*)(*ret))[count]=(int)res;
	    break;
	case M_LONG:
	    if(fabs(res-(long)res)>EPS){
		warning("Floating point number supplied while long integer is needed: {%s}\n", data);
	    }
	    ((long*)(*ret))[count]=(long)res;
	    break;
	case M_DBL:
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
		    error("{%s}: last row has %zu numbers while new row has %zu numbers\n", data, nrow, count-rowbegin);
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
	    error("{%s}: last row has %zu numbers while new row has %zu numbers\n", data, nrow, count-rowbegin);
	}
    }
    if(nrow*ncol!=count){
	error("{%s}: nrow=%zu, ncol=%zu, count=%zu\n", data, nrow, ncol, count);
    }

    if(trans && count>0){
	dbg("Transposing %zux%zu array\n", ncol, nrow);
	void *newer=calloc(count, size);
	switch(type){
	case M_INT:{
	    int *from=(int*)(*ret);
	    int *to=(int*)newer;
	    for(size_t icol=0; icol<ncol; icol++){
		for(size_t irow=0; irow<ncol; irow++){
		    to[icol+ncol*irow]=from[irow+nrow*icol];
		}
	    }
	}
	    break;
	case M_LONG:{
	    long *from=(long*)(*ret);
	    long *to=(long*)newer;
	    for(size_t icol=0; icol<ncol; icol++){
		for(size_t irow=0; irow<ncol; irow++){
		    to[icol+ncol*irow]=from[irow+nrow*icol];
		}
	    }
	}
	    break;
	case M_DBL:{
	    double *from=(double*)(*ret);
	    double *to=(double*)newer;
	    for(size_t icol=0; icol<ncol; icol++){
		for(size_t irow=0; irow<ncol; irow++){
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
	    *ret=myrealloc(*ret, size*count,char);
	}else{
	    free(*ret);
	    *ret=NULL;
	}
    }
    if(nrow0) *nrow0=nrow;
    if(ncol0) *ncol0=ncol;
    return count;
}
int readstr_intarr(int**ret, int len, const char *data){
    return readstr_numarr((void**)ret, len,NULL,NULL, M_INT, data);
}
/**
   Read an integer array. Duplicate if only one number is present.
 */
void readstr_intarr_nmax(int **ret, /**<[out] Result*/
			 int len,   /**<[in]  Max number of values to read.*/
			 const char *data /**<[in] Input string*/
			 ){
    int len2=readstr_intarr(ret, len, data);
    if(len2==1){
	for(int i=1; i<len; i++){
	    (*ret)[i]=(*ret)[0];
	}
    }else if(len2!=0 && len2!=len){
	error("{%s}: Require %d numbers, but got %d\n", data, len, len2);
    }
}
void readstr_intarr_relax(int **ret, /**<[out] Result*/
			  int len,   /**<[in]  Max number of values to read.*/
			  const char *data /**<[in] Input string*/
    ){
    int *ret0=0;
    int len2=readstr_intarr(&ret0, 0, data);
    int *ret2=ret0;
    if(len2==1){
	for(int i=0; i<len; i++){
	    (*ret)[i]=ret2[0];
	}
    }else if(len2>=len){
	for(int i=0; i<len; i++){
	    (*ret)[i]=ret2[i];
	}
    }else{
	error("{%s}: Require %d numbers, but got %d\n", data, len, len2);
    }
    free(ret2);
}
/**
   Search and return the value correspond to key. NULL if not found. Do not free the
   returned pointer. The key must be preceeded by space, semicolon, coma or new line (isspace),
   and succeeded by = sign. */
const char *search_header(const char *header, const char *key){
    if(!header) return NULL;
    const char *ans=NULL;
    //const char *ans_bak=NULL;
    const char *val=header;
    while(val[0]!='\0' && (val=strstr(val, key))){
	if(val>header){
	    char prev=*(val-1);
	    if(!isspace((int)prev) && prev!=';' && prev !=','){
		//ans_bak=val;
		val=val+strlen(key);
		continue;/*Invalid */
	    }
	}
	val=val+strlen(key);
	while(val[0]==' ') val++;
	if(val[0] == '='){
	    val++;
	}else{
	    continue;//invalid key
	}
	while(val[0]==' ') val++;
	ans=val;
	break;
    }
    //if(!ans) ans=ans_bak;
    return ans;
}
/**
   Read a number from the header with key
*/
double search_header_num(const char *header, const char *key){
    if(!header) return NAN;
    const char *val=search_header(header, key);
    if(val){
	return readstr_num(val, NULL);
    }else{
	return NAN;/*not found. */
    }
}
/**
   Read a number from the header and verify.
*/
double search_header_num_valid(const char *header, const char *key){
    double val=search_header_num(header, key);
    if(!(val==val)){
	error("{%s}: Unable to parse a number for %s from %s\n", header, key, search_header(header, key));
    }
    return val;
}

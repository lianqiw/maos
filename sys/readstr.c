/*
  Copyright 2009-2026 Lianqi Wang
  
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
#include <ctype.h> /*isspace */
#include <strings.h> //strncasecmp
#include "readstr.h"
#include "misc.h"
#include "bin.h"
#include "scheduler_client.h"
/**
   Group all routines that are used to parse values from string that contain
   key=value pairs.
*/

/**
   Obtain a string array value from the key. String entries must be separated by
   coma [,]. Enclosing with quote is optional, but necessary to protect spaces.
   Examples:
   []: one empty string
   [,]: two empty strings
   [a,b]: two strings: "a" and "b";
   [a,b,]: three strings, "a", "b", and ""
   ['a b']: one string "a b".
 */
int readstr_strarr(char*** res, /**<[out] Result*/
	int len,     /**<[in] max number f values to read*/
	int relax,	 /**<[in] Whether fewer entries are permitted. If true, will copy from last.*/
	const char *key, /**<[in] the key that needs the value.*/
	const char *sdata /**<[in] Input string*/
){
	int count=0;
	int maxcount=len;
	if(len&&*res){//pre allocatrd data.
		memset(*res, 0, sizeof(char*)*len);
	} else{
		if(!len) maxcount=5;
		*res=mycalloc(maxcount, char*);
	}
	if(!sdata) return count;
	const char* sdataend=0;//sdata+strlen(sdata);
	const char* sdata2=sdata;
	trim_string(&sdata2, &sdataend);
	if(sdata[0]=='['){
		sdata2++;
		sdataend--;
		if(sdataend[0]!=']'){
			error("%s=%s: unmatched [] in input.\n", key, sdata);
		}
	}
	/*skip spaces*/
	while(sdata2<sdataend&&isspace(sdata2[0])){
		sdata2++;
	}
	int end_coma=0;//end with coma. append an additional element.
	while(sdata2<sdataend||end_coma){
		//sdata2 is current value star
		//sdata4 is current value end
		//sdata3 is next value start		
		const char* sdata4=sdataend;
		const char* sdata3;
		end_coma=0;
		//If entry is quoted
		if(sdata2[0]=='"'||sdata2[0]=='\''){
			char sep=sdata2[0];
			sdata2++;
			sdata4=strchr(sdata2, sep);//find matching quote.
			if(!sdata4){
				error("%s=%s: unmatched quote in input\n", key, sdata);
			}
			sdata3=sdata4+1;
			//Skip spaces
			while(sdata3<sdataend&&isspace(sdata3[0])){
				sdata3++;
			}
			//separator following the quote.
			if(sdata3[0]==','/*||sdata3[0]==';'*/){
				sdata3++;
				end_coma=1;
			}
		} else{
			//Find the next separator
			char sep[]=",";//entries are separated by ,
			for(size_t is=0; is<sizeof(sep); is++){
				const char* tmp=strchr(sdata2, sep[is]);
				if(tmp&&tmp<sdata4&&tmp<sdataend){
					sdata4=tmp;
					break;
				}
			}
			if(sdata4[0]==','/*||sdata4[0]==';'*/){
				end_coma=1;
			}
			//skip the separator
			sdata3=sdata4+1;
		}
		if(count>=maxcount){//check memory
			if(!len){
				maxcount*=2;
				*res=myrealloc(*res, maxcount, char*);
			} else{
				error("%s=%s: need %d numbers, but got more than %d.\n", key, sdata, len, count);
			}
		}
		if(sdata4>sdata2){/*found non-empty str*/
			trim_string(&sdata2, &sdata4);
			(*res)[count++]=mystrndup(sdata2, sdata4-sdata2);
		} else{/*found empty str*/
			(*res)[count++]=NULL;
		}
		sdata2=sdata3;
		while(sdata2<sdataend&&(sdata2[0]==' ')){//skip space
			sdata2++;
		}
	}
	if(count && count<len){//not enough entry
		if(!relax||(relax==1&&count>1)){
			error("%s=%s: require %d numbers, but got %d.\n", key, sdata, len, count);
		}else{
			dbg3("%s=%s: Fill %d to %d by value in %d\n", key, sdata, count, len, count-1);
			for(int i=count; i<len; i++){
				(*res)[i]=(*res)[count-1]?strdup((*res)[count-1]):NULL;
			}
		}
	}else if(!len){
		if(count>0){
			*res=myrealloc(*res, count, char*);
		} else{
			free(*res); *res=NULL;
		}
	}
	return count;
}
/**
 * Parse the extra expression to the form _*mul+add. It supports +,-,*,/ operators with no space before the operator.
 * */
static void parse_expr(double* vmul, double* vadd, char** pendptr){
	*vmul=1;
	*vadd=0;
	if(!pendptr || !*pendptr) {
		dbg("parse_expr: string is empty.\n");
		return;
	}
	char *endptr=*pendptr;
	double valpriv=1;
	int optpriv=0;

	while(endptr[0]=='/'||endptr[0]=='*'||endptr[0]=='+'||endptr[0]=='-'){
		char op=endptr[0];
		endptr++;
		while(isspace(endptr[0])) endptr++;
		char *data=endptr;
		double tmp=strtod(data, (char**)&endptr);
		if(data==endptr){
			error("%s: Unable to read a number after operator %c\n", *pendptr, op);
			break;
		}
		//convert - to + and / to * for easy handling
		if(op=='-'){
			tmp=-tmp; 
			op='+';
		}else if(op=='/'){
			tmp=1./tmp; 
			op='*';
		}
		//For + operator, commit previous stage if any, stage it and check the next operator
		//For * operator, apply to the stage or the result.
		if(op=='+'){
			if(optpriv){
				*vadd+=valpriv;//add previous val to result
			}
			valpriv=tmp;//record addition for later use
			optpriv=1;
		}else if(op=='*'){
			if(optpriv){//multiply val to previous value and keep it
				valpriv*=tmp;
			}else{//multiply to result
				*vmul*=tmp;
			}
		}
	}
	if(optpriv){
		*vadd+=valpriv;
	}
	*pendptr=endptr;
}

/**
   Read in a number from the value string. Will interpret +,0,*,/ operators if
   there is nospace in between. *endptr0 will be updated to point to the next
   valid entry, or at separator like coma (spaced are skipped).  */
double readstr_num(const char *key, /**<[in] the key that needs the value.*/
	const char *data, /**<[in] Input string*/
	char** endptr0    /**<[out] Location in Input string after readed number.*/
){
	char* endptr;
	while(isspace(data[0]) && data[0]!='\0') data++;
	if(!data || data[0]=='\0'){
		warning("%s: cannot parse a number from an empty string. Assume 0.\n", key);
		return 0;
	}
	double res=strtod(data, &endptr);
	if(data==endptr){
		error("%s=%s: unable to parse for a number, assume nan.\n", key, data);
		return NAN;
	}
	double vmul, vadd;
	parse_expr(&vmul, &vadd, &endptr);
	if(endptr0) *endptr0=endptr;
	return res*vmul+vadd;
}
/**
 * @brief Append the array with new numbers or duplicate the last entry ndup times.
 * 
 * @param ret 	The vector for containing the numbers
 * @param len 	Total number of elements requested (ret is pre-allocated). 0 means no limit (ret is no pre-allocated)
 * @param type 	The date type of the vector in memory
 * @param size	The data element size
 * @param count The existing number of elements in the vector.
 * @param res 	The number to append
 * @param ndup 	Number of times to duplicate the last entry.
 * @return -1 if invalid entry, -2 if more numbers are supplied, 0 if success.
 */
int fill_num(void **ret, int len, int *nmax, int type, int size, int *count, double res, int nfill){
	if(len>0 && *count+nfill>len){
		//warning("Needs %d numbers, but more are supplied.\n", len);
		*count+=nfill;
		return -2;
	}
	for(int i=0; i<nfill; i++){
		if(*count>=*nmax){
			*ret=myrealloc(*ret, size*(*nmax)*2, char);
			memset((char *)*ret+size*(*nmax), 0, size*(*nmax));
			*nmax*=2;
		}
		/*assign the value to appropriate array. convert to int if necessary. */
		switch(type){
		case M_INT:
			if(fabs(res-(int)res)>EPS){
				warning("Floating point number supplied while integer is required.\n") ;
			}
			((int *)(*ret))[*count]=(int)res;
			break;
		case M_LONG:
			if(fabs(res-(long)res)>EPS){
				warning("Floating point number supplied while integer is required.\n");
			}
			((long *)(*ret))[*count]=(long)res;
			break;
		case M_DBL:
			((double *)(*ret))[*count]=res;
			break;
		case M_FLT:
			((float *)(*ret))[*count]=(float)res;
			break;
		default:
			error("Invalid type");return -1;
		}
		(*count)++;
	}
	return 0;
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

   when there is no [], only read until there is ';', '\n', or end of str.
   \return Number of values actually read not including duplicated when relax is set.
*/
int readstr_numarr(void **ret, /**<[out] Result*/
	int *nrow0, /**<[out] Number of rows (nx)*/
	int *ncol0, /**<[out] Number of columns (ny)*/
	int len,    /**<[in]  Max number of values to read.*/
	int relax,  /**<[in] Whether fewer entries are permitted. If true, will copy from last.*/
	int type,   /**<[in]  Data type*/
	const char *key, /**<[in] the key that needs the value.*/
	const char *data /**<[in] Input string*/
){
	if(!ret){
		warning("%s=%s: ret is not set\n", key, data);
		return 0;
	}
	int size=0;/*size of each number */

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
	case M_FLT:
		size=sizeof(float);
		break;
	default:
		error("%s=%s: invalid type", key, data); return -1;
	}
	int nmax=10;
	if(len==0){
		if(!(*ret=calloc(nmax, size))){
			error("%s=%s: failed to allocate memory for ret\n", key, data);return -1;
		}
	} else{
		nmax=len;
		if(*ret){
			memset(*ret, 0, size*len);
		} else{
			*ret=calloc(len, size);
		}
	}
	if(!data||strlen(data)==0){
		return 0;
	}
	const char *startptr=data;
	const char *endptr=0;//temporary pointer
	const char *endarr=0;//indicate end of array parsing
	const char *bopen=0;//bracket open
	const char *bclose=0;//bracket close
	double fact=1;
	double addval=0;
	int power=1;
	int trans=0;
	
	/*process possible numbers before the array. */
	bopen=strchr(startptr, '[');
	const char *sep=strchr(startptr, ';');
	if(sep && bopen && sep<bopen){
		bopen=NULL;//; appears before [, no array.
	}
	if(bopen){/*there is indeed '[' */
		while(startptr[0]!='['){/*parse expression before [. only * and / are allowed*/
			double fact1=strtod(startptr, (char**)&endptr);/*get the number */
			if(startptr==endptr){
				error("%s=%s: invalid entry to parse for numerical array\n", key, data);return -1;
			} else{/*valid number */
				if(power==1){
					fact*=fact1;
				} else{
					fact/=fact1;
				}
				while(isspace(endptr[0])) endptr++;
				if(endptr[0]=='/'){
					power=-1;
				} else if(endptr[0]=='*'){
					power=1;
				} else{
					error("%s=%s: Invalid entry to parse for numerical array.\n", key, data);return -1;
				}
				startptr=endptr+1;
			}
		}
		if(startptr[0]!='['){
			error("%s=%s: Invalid entry to parse for numerical array\n", key, data);return -1;
		}
		startptr++;/*points to the beginning of the array in [] */
		/*process possible numbers after the array. do not use startptr here.*/
	} else{
		/*warning("Expecting array: {%s} should start with [\n", data); */
		endarr=strchr(startptr, ';');
		if(!endarr) endarr=strchr(startptr, '\n');
		if(!endarr) endarr=startptr+strlen(startptr);
	}
	bclose=strchr(startptr, ']');
	if(!bopen && bclose && bclose>endarr){
		bclose=NULL;//] appears after endarr, ignore it.
	}
	if(bclose){/*there is indeed ']'. Handle operations after ] */
		endarr=bclose;
		endptr=bclose+1;
		while(isspace(endptr[0])||endptr[0]=='\''){
			if(endptr[0]=='\'') trans=1-trans;
			endptr++;
		}
		double vmul, vadd;
		parse_expr(&vmul, &vadd, (char**)&endptr);
		addval+=vadd;
		fact*=vmul;
		while(isspace(endptr[0])) endptr++;
		if(!is_end(endptr[0])&&endptr[0]!=';'&&endptr[0]!=','){
			error("%s=%s: There is garbage in the end of the string.\n", key, data);return -1;
		}
	}
	if((bopen!=NULL) != (bclose!=NULL)){
		error("%s=%s: unmatched [] in input.\n", key, data);return -1;
	}
	int count=0;
	int nrow=0;/*number of rows*/
	int ncol=0;/*number of columns*/
	int rowbegin=0;/*beginning of this row*/
	double res=0;
	/*Read in the array between [ and ]*/
	while(startptr<endarr){
		/*parse the string for a floating point number.  */
		res=readstr_num(key, startptr, (char**)&endptr);
		if(startptr==endptr){
			warning("%s=%s: unable to parse a number at %s\n", key, data, startptr);
			break;
		}else startptr=endptr;
		/*apply the factors appear before or after [] */
		if(power==1){
			res=fact*res;
		} else{
			res=fact/res;
		}
		res+=addval;
		
		if(fill_num(ret, len, &nmax, type, size, &count, res, 1)==-1){
			warning("%s=%s: read failed\n", key, data);
			return -1;
		}
		
		if(startptr<endarr){//more data to read
			if(!(isspace(startptr[0])||startptr[0]=='@'|| startptr[0]==';'||startptr[0]==',')){
				error("%s=%s: garbage is found: {%s}\n", key, data, startptr);return -1;
			}
		}
		/*process the number separators. */
		while(startptr<endarr && isspace(startptr[0])) startptr++;//continuous spaces are ignored
		/*if the number is followed by @n, it is repeated n times. */
		if(startptr[0]=='@'){
			startptr++;
			int nrep=strtol(startptr, (char **)&endptr, 10);
			if(startptr==endptr){
				nrep=0;
				error("Failed to read an integer after @: %s\n", startptr); return -1;
			} else{
				startptr=endptr;
				if(fill_num(ret, len, &nmax, type, size, &count, res, nrep-1)==-1){
					warning("%s=%s: read failed\n", key, data);
					return -1;
				}
			}
		}
		while(startptr<endarr&&isspace(startptr[0])) startptr++;//continuous spaces are ignored
		if(startptr<endarr && startptr[0]==','){
			startptr++; //a single coma is permitted
		} else if(startptr<endarr && startptr[0]==';'){//; is used to separate into a new row
			ncol++;
			if(nrow==0){
				nrow=count-rowbegin;
			} else if(nrow!=count-rowbegin){
				error("%s=%s: previous row has %d numbers while the current row has %d numbers\n", key, data, nrow, count-rowbegin); return -1;
			}
			rowbegin=count;
			startptr++;
		}
		while(startptr<endarr && isspace(startptr[0])) startptr++;//continuous spaces are ignored
	}
	if(startptr!=endarr){
		error("%s=%s: garbage is found: {%s}\n", key, data, startptr);return -1;
	}
	//postprocessing to satisfy request.
	if(rowbegin<count){/*if last row is not ended with ;*/
		ncol++;
		if(nrow==0){
			nrow=count-rowbegin;
		} else if(nrow!=count-rowbegin){
			error("%s=%s: previous row has %d numbers while the current row has %d numbers\n", key, data, nrow, count-rowbegin); return -1;
		}
	}
	if(nrow*ncol!=count){
		error("%s=%s: nrow=%d, ncol=%d, count=%d\n", key, data, nrow, ncol, count); return -1;
	}
	//count is gauranteed to not exceed len
	if(count>0&&count<len){//not enough values are read
		if(!relax||(relax==1&&count>1)||trans){
			error("%s=%s: require %d numbers, but got %d\n", key, data, len, count);return -1;
		} else{//Fill the array with the last number
			dbg3("%s=%s: fill %d to %d by value in %d\n", key, data, count, len, count-1);
			if(fill_num(ret, len, &nmax, type, size, &count, res, len-count)==-1){
				warning("%s=%s: read failed\n", key, data);
				return -1;
			}
			if(nrow*ncol!=len){//vector
				ncol=1;
				nrow=len;
			}
		}
	}else if(trans&&count>0){
		dbg("%s=%s: transposing %dx%d array\n", key, data, ncol, nrow);
		void* newer=calloc(count, size);
#define DO_TRANS(T)						\
		{							\
			T *from=(T*)(*ret);					\
			T *to=(T*)newer;					\
			for(int icol=0; icol<ncol; icol++){		\
			for(int irow=0; irow<ncol; irow++){		\
				to[icol+ncol*irow]=from[irow+nrow*icol];	\
			}						\
			}							\
		}

		switch(type){
		case M_INT:
			DO_TRANS(int);
			break;
		case M_LONG:
			DO_TRANS(long);
			break;
		case M_DBL:
			DO_TRANS(double);
			break;
		case M_FLT:
			DO_TRANS(float);
			break;
		default:
			error("%s=%s: invalid type", key, data);return -1;
		}
#undef DO_TRANS
		int tmp=ncol; ncol=nrow; nrow=tmp;
		free(*ret);
		*ret=newer;
	} else {
		if(!len){//automatic read
			if(count>0){
				*ret=myrealloc(*ret, size*count, char);
			} else{
				free(*ret);
				*ret=NULL;
			}
		}else{//keep data at zero.
			nrow=len;
			ncol=1;
			if(count>len){
				warning("%s=%s: need %d, but got %d values, ignore the rest.\n", key, data, len, count);
				count=len;
			}
		}
	}
	if(nrow0) *nrow0=nrow;
	if(ncol0) *ncol0=ncol;
	return count;
}
/**
	update header and end to point to valid region. Does not modify the string
*/
void trim_string(const char **pstart, const char **pend){
	if(!pstart) return;
	const char* start=*pstart;
	const char* end=(pend && *pend)?*pend:(start+strlen(start));
repeat:
	while(isspace((int)start[0]) && start<end) start++;
	while(isspace((int)end[-1])&& start<end) end--;
	if(start>=end){ 
		start=NULL;
		end=NULL;
	} else{
		if(start+1<end && (start[0]=='\''||start[0]=='"')){
			if(end[-1]!=start[0]){
				error("Quote is not matched: {%s}\n", *pstart);
			} else{
				end--;
			}
			start++;
			goto repeat;
		}
	}
	*pstart=start;
	if(pend) *pend=end;
}
/**
   Search and return the value correspond to key. Case is ignored; 
   NULL if not found. Do not free the returned pointer. 
   The key must be preceeded by space, semicolon, coma or new line (isspace),
   and succeeded by = sign. */
const char* search_keyword(const char* keywords, const char* key){
	if(!keywords) return NULL;
	const char* ans=NULL;
	//const char* val=header;
	const char* end=NULL;
	trim_string(&keywords, &end);
	if(keywords&&end){
		const int nkey=strlen(key);
		int was_space=1;
		for(const char* p=keywords; p<end; p++){
			const char c=*p;
			if(!isspace((int)c)&&c!=';'&&c!=','){
				if(was_space){//start of key
					if(!strncasecmp(p, key, nkey)){//match regardless of case.
						p+=nkey;
						while(isspace((int)*p)&&p<end) p++;
						if(*p=='='){
							p++;
							while(isspace((int)*p)&&p<end) p++;
							ans=p;
							break;
						} else{
							p--;
						}
					}
				}
				was_space=0;
			} else{
				was_space=1;
			}
		}
	}
	return ans;
}
/**
   Read a number from the header with key
*/
double search_keyword_num(const char* keywords, const char* key){
	if(!keywords) return NAN;
	const char* val=search_keyword(keywords, key);
	if(val){
		return readstr_num(key, val, NULL);
	} else{
		return NAN;/*not found. */
	}
}
/**
   Read a number from the header and verify.
*/
double search_keyword_num_valid(const char* keywords, const char* key){
	double val=search_keyword_num(keywords, key);
	if(!(val==val)){
		error("{%s}: Unable to parse a number for %s from %s\n", keywords, key, search_keyword(keywords, key));
	}
	return val;
}
/**
   Read a number from the header and use value0 if not found.
*/
double search_keyword_num_default(const char* keywords, const char* key, double value0){
	double val=search_keyword_num(keywords, key);
	if(!(val==val)){
		val=value0;
	}
	return val;
}

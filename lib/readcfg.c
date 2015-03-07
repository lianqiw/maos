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

/*
  read the config file, parse the tokens and store into a binary tree. 
  then user can query tokens using read_int, read_double, or read_array
*/



#include <search.h>
#include <ctype.h>


#include <limits.h>
#include "readcfg.h"
/**
   Routines to read .conf type config files. Each entry is composed of a key and
   a value. The entries are maintained in a tree. Each entry can be
   retrieved from the key.
*/
/**
   Compatibility mode: old keys are automatically renamed to new keys.
*/
#define COMPATIBILITY 1

#if COMPATIBILITY == 1
#define RENAME(old,new)							\
    if(!strcmp(var,#old)){						\
	warning2("Deprecated: please change %s to %s.\n",#old,#new);	\
	var=#new;/*strcpy may overflow. just reference the char*/	\
    }
#define IGNORE(old)				\
    if(!strcmp(var,#old)){			\
	ssline[0]='\0';				\
	continue;				\
    }
#else
#define RENAME(old,new) /*do nothing. */
#define IGNORE(old)
#endif
static void *MROOT=NULL;

typedef struct STORE_T{
    char *key;
    char *data;
    long protect;
    long count;
}STORE_T;
static int key_cmp(const void *a, const void *b){
    return strcmp(((STORE_T*)a)->key, ((STORE_T*)b)->key);
}

/*Keep record of all the strings so that we can check whether they have all been used. */
static long nstore=0;/*number of total records */
static long nused=0;/*number of read records */

#define STRICT 1

/**
   trim the spaces or coma before and after string.*/
static void strtrim(char **str){
    if(!*str) return;
    int iend;
    /*turn non-printable characters, coma, and semicolumn, to space */
    for(char *tmp=*str; !is_end(*tmp); tmp++){
	if(!isgraph((int)*tmp) || *tmp==',' || is_space(*tmp)){
	    *tmp=' ';
	}
    }
    /*remove leading spaces. */
    while(!is_end(**str) && is_space((*str)[0])) (*str)++;
    iend=strlen(*str)-1;
    /*remove tailing spaces. */
    while((is_space((*str)[iend]) || (*str)[iend]==';') && iend>=0){
	(*str)[iend]='\0';
	iend--;
    }
    /*remove duplicated spaces */
    char *tmp2=*str;
    int issp=0;
    for(char *tmp=*str; !is_end(*tmp); tmp++){
	if(*tmp==' '){
	    if(issp){/*there is already a space, skip. */
		continue;
	    }else{
		issp=1;
	    }
	}else{
	    issp=0;
	}
	/*copy the string. */
	*tmp2=*tmp;
	tmp2++;
    }
    *tmp2='\0';/*properly end the string */
    if(is_end(**str)) *str=NULL;
}
/**
   remove "" or '' around a string from config file. Accept whole string if
quotes are not found.  */
static char* strextract(const char *data){
    if(!data || !strlen(data)) return NULL;
    if(data[0]=='"' || data[0] =='\''){
	if(data[strlen(data)-1]!=data[0]){
	    error("Record is {%s}, quotes must come in pair\n", data);
	    return NULL;
	}else{
	    char *res=strdup(data+1);
	    res[strlen(res)-1]='\0';
	    if(is_end(res[0])){
		free(res);
		res=NULL;
	    }
	    return res;
	}
    }else{
	return strdup(data);
    }
}
/**
   Remove comment (after #), leading spaces, and trailing line break and spaces from a
   string. 
 */
static char *squeeze(char *line){
    char *sline, *comment;
    if(!line) {
	sline=NULL;
    }else{
	/*Remove comment*/
	comment=strchr(line,'#');
	if(comment){
	    comment[0]='\0';
	}
	/*Remove trailing linebreak, semi-colon, and spaces*/
	int nread=strlen(line)-1;
	while(nread>=0 && (is_space(line[nread]) || is_end(line[nread]))){
	    line[nread]='\0';
	    nread--;
	}
	/*Remove leading spaces*/
	sline=line;
	while(is_space(sline[0])) sline++;
	if(is_end(sline[0]))  sline=NULL;
    }
    return sline;
}
static FILE *fpout=0;
static void print_key(const void *key, VISIT which, int level){
    const STORE_T *store=*((const STORE_T**)key);
    (void)level;
    if(which==leaf || which==postorder){
	if(fpout){
	    if(!store->protect){
		fprintf(fpout, "#");
	    }
	    fprintf(fpout, "%s=", store->key);
	    if(store->data && strcmp(store->data, "ignore")){
		fprintf(fpout, "%s\n", store->data);
	    }else{
		fprintf(fpout, "\n");
	    }
	}
	if(store->count==0){
	    error("key \"%s\" is not recognized, value is %s\n", store->key, store->data);
	}else if(store->count!=1){
	    error("Key %s is used %ld times\n", store->key, store->count);
	}
    }
}
static void delete_leaf(const void *key, VISIT which, int level){
    (void)level;
    if(which==leaf){
	STORE_T *store=*((STORE_T**)key);
	tdelete(store, &MROOT, key_cmp);
	free(store->key);
	free(store->data);
	free(store);
    }
}
/**
   Save all configs to file and check for unused config options.
 */
void close_config(const char *format, ...){
    format2fn;
    const char *fnout=format?fn:NULL;
    if(MROOT){
	info2("Used %ld of %ld supplied keys\n",nused, nstore);
	if(fnout && strlen(fnout)>0 && !disable_save) fpout=fopen(fnout, "w");
	twalk(MROOT, print_key);
	if(fpout) fclose(fpout); fpout=0;
    }
    while(MROOT){
	twalk(MROOT, delete_leaf);
    }
    nused=0;
}
/**
   Start the read config process by opening .conf files and fill the entries in
   a hash table. A key can be protected. If a key is not protected, newer
   entries with the same key will override a previous entry.
 */
void open_config(char* config_file, /**<[in]The .conf file to read*/
		 const char* prefix,/**<[in]if not NULL, prefix the key with this.*/
		 long protect       /**<[in]whether we protect the value*/
		 ){
    if(!config_file) return;
    FILE *fd=NULL;
    char *fn=NULL;
    int print_override=1;
    if(check_suffix(config_file, ".conf")){
	if(exist(config_file)){
	    fn=strdup(config_file);
	}else{
	    fn=find_file(config_file);
	    print_override=0;
	}
	if(!fn || !(fd=fopen(fn,"r"))){
	    perror("fopen");
	    error("Cannot open file %s for reading.\n",fn);
	}
    }else{
	parse_argopt(config_file, NULL);
	char *end; 
	//Remove trailing space
	for(end=config_file+strlen(config_file); end>=config_file; end--){
	    if(isspace((int)*end)) *end='\0'; 
	    else break;
	}
	if(end<config_file) return;
	fn=config_file;
    }
    
    char *sline=NULL;
    char *var=NULL, *value=NULL;
    int countnew=0;
    int countold=0;
    
#define MAXLN 40960
    char ssline[MAXLN];
    ssline[0]='\0';/*stores the available line. */
    char line[MAXLN];
    while(1){
	if(fd){/*read from file*/
	    if(!fgets(line, MAXLN, fd)) break;
	}else{/*read from string*/
	    if(!config_file) break;
	    char *p0=strchr(config_file, '\n');
	    int len;
	    if(p0){
		len=p0-config_file;
	    }else{
		len=strlen(config_file);
	    }
	    if(len+1>MAXLN){
		error("Input line is too long. Please make MAXLN larger to accomodate.\n");
	    }
	    strncpy(line, config_file, len);
	    line[len]='\0';
	    if(p0){
		config_file=p0+1;
	    }else{
		config_file=0;
	    }
	}
	sline=squeeze(line);
	if(!sline || is_end(sline[0]))/*skip empty lines. */
	    continue;
	if(sline){
	    if((strlen(sline)+strlen(ssline))>=MAXLN){
		error("Input line is too long. Please make MAXLN larger to accomodate.\n");
	    }
	    strcat(ssline, sline);
	}
	/*lines ended with \ will continue on next line */
	if(ssline[strlen(ssline)-1]=='\\'){
	    ssline[strlen(ssline)-1]='\0';
	    continue;
	}
	char *eql=index(ssline,'=');
	if(!eql){//no equal sign
	    if(check_suffix(ssline, ".conf")){
		char *embeded=strextract(ssline);
		if(embeded){
		    open_config(embeded, prefix, protect);
		    free(embeded);
		}
	    }else if(!strcmp(ssline, "__protect_start")){
		protect+=10000;
	    }else if(!strcmp(ssline, "__protect_end")){
		protect-=10000; 
		if(protect<0){
		    error("__protect_end must appear after __protect_start, in the same file");
		}
	    }else{
		error("Input (%s) is not valid\n", ssline);
	    }
	    ssline[0]='\0';
	    continue;
	}
	int append=0;
	if(eql[-1]=='+'){
	    append=1;//append to key
	    eql[-1]='\0';
	}
	eql[0]='\0';
	var=ssline;
	value=eql+1;
	strtrim(&var);
	strtrim(&value);
	if(!var || strlen(var)==0){
	    error("Line '%s' is invalid\n",line);
	}else if(!strcmp(var,"path") || !strcmp(var, "PATH")){
	    char *val2=strextract(value);
	    addpath(val2);
	    free(val2);
	}else if(!strcmp(var,"include")){
	    /*info("Opening embeded config file %s\n",value); */
	    char *embeded=strextract(value);
	    if(embeded){
		open_config(embeded,prefix,protect);
		free(embeded);
	    }
	}else{
#if COMPATIBILITY == 1	    
	    /*
	      Compatibility mode: rename old key names to new key names. Will
	      remove in the future.
	     */
	    RENAME(atm.zadeg, sim.zadeg);
	    /*Added on 2011-04-28 */
	    RENAME(dbg.noatm, sim.noatm);
	    RENAME(dbg.fitonly, sim.fitonly);
	    RENAME(dbg.evlol, sim.evlol);
	    RENAME(sim.fitonly, sim.idealfit);
	    IGNORE(sim.servotype_hi);
	    IGNORE(sim.servotype_lo);
	    IGNORE(sim.epfocus);
	    IGNORE(dbg.dxonedge);
	    RENAME(sim.gtypeII_lo, sim.eplo);
	    RENAME(sim.epngs, sim.eplo);
	    RENAME(sim.apngs, sim.aplo);
	    RENAME(dbg.splitlrt, tomo.splitlrt);
	    RENAME(tomo.split_wt, tomo.ahst_wt);
	    RENAME(tomo.split, recon.split);
	    RENAME(evl.opdcov, evl.cov);
	    RENAME(evl.psfpttr, evl.pttr);
#endif
	    if(value && !strcmp(value, "ignore")){
		ssline[0]='\0';
		continue;
	    }
	    STORE_T *store=calloc(1, sizeof(STORE_T));
	    if(prefix){
		store->key=stradd(prefix,var,NULL);
	    }else{
		store->key=strdup(var);
	    }

	    if(value && strlen(value)>0){
		store->data=strdup(value);
	    }else{
		store->data=NULL;
	    }
	    store->protect=protect;
	    store->count=0;

	    void **entryfind=tfind(store, &MROOT, key_cmp);
	    if(entryfind){ 
		/*same key found */
		STORE_T *oldstore=*entryfind;
		if(append){
		    /*concatenate new value with old value for arrays. both have to start/end with [/]*/
		    const char *olddata=oldstore->data;
		    const char *newdata=store->data;
		    const int nolddata=strlen(olddata);
		    const int nnewdata=strlen(newdata);
		    if(olddata[0]!='[' || olddata[nolddata-1]!=']'){
			error("olddata='%s' should be encapsulated by bare [].\n", olddata);
		    }else if(newdata[0]!='[' || newdata[nnewdata-1]!=']'){
			error("newdata='%s' should be encapsulated by bare [].\n", newdata);
		    }else{
			oldstore->data=realloc(oldstore->data, (nolddata+nnewdata));
			oldstore->data[nolddata-1]=' ';
			strncat(oldstore->data, newdata+1, nnewdata-1);
		    }
		}else{
		    if(oldstore->protect<=protect){
			if(print_override && 
			   (((oldstore->data==NULL || store->data==NULL)
			     &&(oldstore->data != store->data))||
			    ((oldstore->data!=NULL && store->data!=NULL)
			     &&strcmp(oldstore->data, store->data)))){
			    info2("Overriding %-20s\t{%s}-->{%s}\n", 
				  store->key, oldstore->data, store->data);
			}
			/*free old value */
			free(oldstore->data);
			/*move pointer of new value. */
			oldstore->data=store->data; store->data=0;
			oldstore->protect=store->protect;
			oldstore->count=store->count;
		    }
		}
		countold++;
		free(store->data);
		free(store->key);
		free(store);
	    }else{
		/*new key */
		if(!tsearch(store, &MROOT, key_cmp)){
		    error("Error inserting to tree\n");
		}
		countnew++;
		nstore++;
	    }
	}
	ssline[0]='\0';
    }
    info2("loaded %3d (%3d new) records from '%s'\n",countnew+countold,countnew, fn);
    if(fd){
	fclose(fd);
	free(fn);
    }
#undef MAXLN
}
/**
   Get the record number of a key.
 */
static const STORE_T* getrecord(char *key, int mark){
    STORE_T store;
    void **found=0;
    strtrim(&key);
    store.key=key;
    if((found=tfind(&store, &MROOT, key_cmp))){
	if(mark){
	    if((*(STORE_T**)found)->count){
		error("This record %s is already read\n",key);
	    }
	    (*(STORE_T**)found)->count++;
	    nused++;
	}
    }else if(mark){
	print_file("change.log");
	error("Record %s not found\n",key);
    }
    return found?(*found):0;
}
/**
   Check whether a have a record of a key.
 */
int readcfg_peek(const char *format,...){
    /*Check whether key exists */
    format2key;
    return getrecord(key, 0)?1:0;
}
/**
   Check the size of an array input
*/
int readcfg_peek_n(const char *format, ...){
    format2key;
    const STORE_T *store=getrecord(key, 0);
    if(!store) return 0;
    const char *sdata=store->data;
    const char *startptr=strchr(sdata,'[');
    const char *endptr=strchr(sdata,']');
    if(!startptr) startptr=sdata;
    if(!endptr) endptr=startptr+strlen(sdata)-1;
    const char *quote=strchr(startptr,'"');
    int count=0;
    if(quote && quote<endptr){/*this is string array */
	char **ret=NULL;
	count=readstr_strarr(&ret, 0, sdata);
	for(int i=0; i<count; i++){
	    free(ret[i]);
	}
	free(ret);
    }else{/*this is numerical array */
	void *ret;
	count=readstr_numarr(&ret, 0, NULL,NULL,T_DBL, sdata);
	free(ret);
    }
    return count;
}
/**
   Check whether the record is overriden by user supplied conf files.
*/
int readcfg_peek_override(const char *format,...){
    /*Check whether key exists */
    format2key;
    const STORE_T *store=getrecord(key, 0);
    if(!store){
	return 0;
    }else{
	return store->protect;
    }
}
/**
   Obtain a string value from the key.
 */
char *readcfg_str(const char *format,...){
    format2key;
    char *data;
    const STORE_T *store=getrecord(key, 1);
    if(store){
	const char *sdata=store->data;
	if(sdata && strlen(sdata)>0){
	    data=strextract(sdata);
	}else{
	    data=NULL;
	}
    }else{
	error("key '%s' not found\n", key);
	data=NULL;
    }
    return data;
}
int readcfg_strarr(char ***res, const char *format,...){
   /*Read str array. */
    format2key;
    *res=NULL;/*initialize */
    const STORE_T *store=getrecord(key, 1);
    if(!store){/*record not found. */
	error("key '%s' not found\n", key);
	return 0;
    }else{
	const char *sdata=store->data;
	if(!sdata){/*record is empty. */
	    return 0;
	}
	return readstr_strarr(res, 0, sdata);
    }
}

/**
   Read integer array
*/
int readcfg_intarr(int **ret, const char *format,...){
    format2key;
    return readstr_numarr((void**)ret, 0,NULL,NULL, T_INT, getrecord(key, 1)->data);
}
/**
   Read as an lmat.
 */
lmat *readcfg_lmat_do(int n, char *key){
    long *val=0;
    long **ret=&val;
    int nx, ny;
    readstr_numarr((void**)ret, n, &nx, &ny, T_LONG, getrecord(key, 1)->data);
    lmat *res=0;
    if(!nx || !ny){
	free(val); val=0;
    }
    res=lnew_data(nx,ny, val);
    return res;
}
/**
   Read as an lmat. 
 */
lmat *readcfg_lmat(const char *format,...){
    format2key;
    return readcfg_lmat_do(0, key);
}
/**
   Read as an lmat. Exactly n numbers if n>0
 */
lmat *readcfg_lmat_n(int n, const char *format,...){
    format2key;
    lmat *out=readcfg_lmat_do(n, key);
    int nread=out?(out->nx*out->ny):0;
    if(n!=0 && nread!=n){
	error("Need %d elements, got %d\n", n, nread);
    }
    return out;
}
/**
   Read as an lmat. A max of n numbers
 */
lmat *readcfg_lmat_nmax(int n, const char *format,...){
    format2key;
    lmat *out=readcfg_lmat_do(n, key);
    long nread=out?(out->nx*out->ny):0;
    if(nread<=1){
	lresize(out, n, 1);
	if(nread==1){
	    for(int i=1; i<n; i++){
		out->p[i]=out->p[0];
	    }
	}
    }else if(nread!=n){
	error("Need %d elements, got %ld\n", n, nread);
    }
    return out;
}
/**
   Read double array
*/
int readcfg_dblarr(double **ret, const char *format,...){
    format2key;
    return readstr_numarr((void**)ret, 0,NULL,NULL,T_DBL, getrecord(key, 1)->data);
}
/**
   Read as a dmat. It can be a file name or an array.
 */
dmat *readcfg_dmat_do(int n, char *key){
    double *val=NULL;
    char *str=getrecord(key, 1)->data;
    if(!str){
	return 0;
    }else if((str[0]<='Z' && str[0]>='A')
       || (str[0]<='z' && str[0]>='a')
       ||(str[0]=='"' || str[0]=='\'')){
	char *fn=strextract(str);
	return dread("%s", fn);
	free(fn);
    }else{
        double **pval=&val;
	int nx, ny;
	readstr_numarr((void**)pval, n, &nx, &ny,T_DBL, str);
	dmat *res=0;
	if(!nx || !ny) {
	    free(val); val=0;
	}
	res=dnew_data(nx, ny, val);
	return res;
    }
}
/**
   Read as a dmat. It can be a file name or an array.
 */
dmat *readcfg_dmat(const char *format,...){
    format2key;
    return readcfg_dmat_do(0, key);
}
/**
   Read as a dmat. It can be a file name or an array.
 */
dmat *readcfg_dmat_n(int n, const char *format,...){
    format2key;
    dmat *out=readcfg_dmat_do(n, key);
    long nread=out?(out->nx*out->ny):0;
    if(n!=0 && nread!=n){
	error("Need %d elements, got %ld\n", n, nread);
    }
    return out;
}
/**
   Read as a dmat. It can be a file name or an array.
 */
dmat *readcfg_dmat_nmax(int n, const char *format,...){
    format2key;
    dmat *out=readcfg_dmat_do(n, key);
    long nread=out?(out->nx*out->ny):0;
    if(nread<=1){
	dresize(out, n, 1);
	if(nread==1){
	    dset(out, out->p[0]);
	}
    }else if(nread!=0 && nread!=n){
	error("Need %d elements, got %ld\n", n, nread);
    }
    return out;
}
/**
   Read string array of len elements
*/
void readcfg_strarr_n(char ***ret, int len, const char *format,...){
    format2key;
    int len2;
    if(len!=(len2=readstr_strarr((char***)ret, len, getrecord(key, 1)->data))){
	error("%s: Require %d elements, but got %d\n", key, len, len2);
    }
}
/**
   Read str array of upto len elements
*/
void readcfg_strarr_nmax(char ***ret, int len, const char *format,...){
    format2key;
    int len2=readstr_strarr((char***)ret, len, getrecord(key, 1)->data);
    if(len2==1){
	for(int i=1; i<len; i++){
	    (*ret)[i]=(*ret)[0]?strdup((*ret)[0]):NULL;
	}
    }else if(len2!=0 && len2!=len){
	error("%s: Require %d elements, but got %d\n", key, len, len2);
    }
}
/**
   Read integer array of len elements
*/
void readcfg_intarr_n(int **ret, int len, const char *format,...){
    format2key;
    int len2;
    if(len!=(len2=readstr_numarr((void**)ret, len, NULL,NULL,T_INT, getrecord(key, 1)->data))){
	error("%s: Need %d, got %d integers\n", key, len, len2);
    }
}
/**
   Read integer array of maximum of len elements
*/
void readcfg_intarr_nmax(int **ret, int len, const char *format,...){
    format2key;
    int len2=readstr_numarr((void**)ret, len,NULL,NULL, T_INT, getrecord(key, 1)->data);
    if(len2==1){
	for(int i=1; i<len; i++){
	    (*ret)[i]=(*ret)[0];
	}
    }else if(len2!=0 && len2!=len){
	error("%s: Require %d numbers, but got %d\n", key, len, len2);
    }
}
/**
   Read double array of len elements
*/
void readcfg_dblarr_n(double **ret, int len, const char *format,...){
    format2key;
    int len2;
    if(len!=(len2=readstr_numarr((void**)ret, len,NULL,NULL, T_DBL, getrecord(key, 1)->data))){
	error("%s: Need %d, got %d double\n", key, len, len2);
    }
}
/**
   Read double array of len elements
*/
void readcfg_dblarr_nmax(double **ret, int len, const char *format,...){
    format2key;
    int len2=readstr_numarr((void**)ret, len, NULL,NULL,T_DBL, getrecord(key, 1)->data);
    if(len2==1){
	for(int i=1; i<len; i++){
	    (*ret)[i]=(*ret)[0];
	}
    }else if(len2!=0 && len2!=len){
	error("%s: Require %d numbers, but got %d\n", key, len, len2);
    }
}
/**
   Read integer
*/
int readcfg_int( const char *format,...){
    format2key;
    char *val=getrecord(key, 1)->data;
    char *endstr;
    double ans=readstr_num(val, &endstr);
    if(fabs(ans-(int)ans)>EPS){
	warning("Floating point number supplied while integer is needed: %s=%s\n", key, val);
    }
    if(endstr[0]!='\0'){
	error("Garbage found in %s=%s.\n", key, val);
    }
    return (int)ans;
}
/**
   Read double
*/
double readcfg_dbl(const char *format,...){
    format2key;
    char *val=getrecord(key, 1)->data;
    char *endstr;
    double ans=readstr_num(val, &endstr);
    if(endstr[0]!='\0'){
	error("Garbage found in %s=%s.\n", key, val);
    }
    return ans;
}

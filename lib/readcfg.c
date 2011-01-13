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

/*
  read the config file, parse the tokens and store into a hash table. 
  then user can query tokens using read_int, read_double, or read_array
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <search.h>
#include <ctype.h>
#include <stdarg.h>
#include <signal.h>
#include <limits.h>
#include "readcfg.h"
#include "path.h"
#include "common.h"
/**
   file readcfg.c

   Routines to read .conf type config files. Each entry is composed of a key and
   a value. The entries are maintain ed a hash table. Each entry can be
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
#else
#define RENAME(old,new) //do nothing.
#endif
static int MAX_ENTRY=1000;
typedef struct STORE_T{
    char *key;
    char *data;
    long protect;
    long count;
}STORE_T;
//Keep record of all the strings so that we can check whether they have all been used.
static STORE_T *store=NULL;//store the pointers.
static long nstore=0;//number of total records
static long nused=0;//number of read records
static int print_override=1;

#define STRICT 1
/*
   We only use one hash table in the beginning of code. so the non-reentrant hsearch is fine.
*/
/**
trim the spaces, ", ', before and after string.*/
static void strtrim(char **str){
    if(!*str) return;
    int iend;
    //turn non-printable characters, coma, and semicolumn, to space
    for(char *tmp=*str; *tmp!='\0'; tmp++){
      if(!isgraph((int)*tmp) || *tmp==';' || *tmp==',' || isspace((int)*tmp)){
	    *tmp=' ';
	}
    }
    //remove leading spaces.
    while((*str)[0]!='\0' && isspace((int)(*str)[0])) (*str)++;
    iend=strlen(*str)-1;
    //remove tailing spaces.
    while(isspace((int)(*str)[iend]) && iend>=0){
	(*str)[iend]='\0';
	iend--;
    }
    //remove duplicated spaces
    char *tmp2=*str;
    int issp=0;
    for(char *tmp=*str; *tmp!='\0'; tmp++){
	if(*tmp==' '){
	    if(issp){//there is already a space, skip.
		continue;
	    }else{
		issp=1;
	    }
	}else{
	    issp=0;
	}
	//copy the string.
	*tmp2=*tmp;
	tmp2++;
    }
    *tmp2='\0';//properly end the string
    if((*str)[0]=='\0') *str=NULL;
}
/**
   remove "" or '' around a string from config file. Error if quotes are not
found.  */
static char* strextract(const char *data){
    if(!data || strlen(data)<=2) return NULL;
    if((data[0]!='"' && data[0]!='\'') || 
       (data[strlen(data)-1]!='"' && data[strlen(data)-1]!='\'')){
	error("Record is {%s}. \nThis is not a string constant. "
	      "Strings must be embrased by \"\" or \'\'\n",data);
    }
    char *res=strdup(data+1);
    res[strlen(res)-1]='\0';
    return res;
}
/**
   Remove comment (after #), leading spaces, and trailing line break and spaces from a
   string. Convert ' to " for easy processing later.
 */
static char *squeeze(char *line){
    char *sline, *comment;
    if(!line) {
	sline=NULL;
    }else{
	/*Remove comment*/
	comment=strchr(line,'#');
	if(comment)
	    comment[0]='\0';
	/*Remove trailing linebreak and spaces*/
	int nread=strlen(line)-1;
	while(nread>=0 && line[nread]==' '){
	    line[nread]='\0';
	    nread--;
	}
	/*Remove leading spaces*/
	sline=line;
	while(isspace((int)sline[0])) sline++;
	if(sline[0]=='\0')  sline=NULL;
	//Convert single quote ' to double quotes" if any
	if(sline){
	    const int nx=strlen(sline);
	    for(int ix=0; ix<nx; ix++){
		if(sline[ix]=='\''){
		    sline[ix]='"';
		}
	    }
	}
    }
    return sline;
}
/**
   End the read config process and output the effective config hash table into a
   file for later reference.
 */
void close_config(const char *format, ...){
    format2fn;
    const char *fnout=fn;
    if(store){
	info2("Used %ld of %ld supplied keys\n",nused,nstore);
	if(fnout && strlen(fnout)>0){
	    FILE *fp=fopen(fnout,"w");
	    for(int i=0; i<nstore; i++){
		if(store[i].data)
		    fprintf(fp,"%s=%s\n",store[i].key,store[i].data);
		else
		    fprintf(fp,"%s=\n",store[i].key);
	    }
	    fclose(fp);
	}
	for(int i=0; i<nstore; i++){
	    if(store[i].count==0){
		print_file("change.log");
		error("key \"%s\" is not recognized, value is %s\n", store[i].key,store[i].data);
	    }else if(store[i].count!=1){
		//this should not happen.
		error("Key %s is used %ld times\n", store[i].key,store[i].count);
	    }
#if defined(__linux__)
	    //Linux doesn't free the keys. but Mac does.
	    free(store[i].key);//free key;
#endif
	    if(store[i].data){
		free(store[i].data);//free value
	    }
	}
	hdestroy();
	free(store);
	store=NULL;
    }else{
	warning("Store is empty\n");
    }
}
/**
   Start the read config process by opening .conf files and fill the entries in
   a hash table. A key can be protected. If a key is not protected, newer
   entries with the same key will override a previous entry.
 */
void open_config(const char* config_file, /**<The .conf file to read*/
		 const char* prefix,      /**<if not NULL, prefix the key with this.*/
		 long protect             /**<whether we protect the value*/
		 ){
    if(!config_file) return;
    if(!check_suffix(config_file, ".conf")){
	error("config file '%s' doesn't end with .conf.\n", config_file);
    }
    FILE *fd=NULL;
    char *sline=NULL;
    char *var=NULL, *value=NULL;
    ENTRY entry;
    ENTRY *entryfind=NULL;
    int recordcount=0;
    if(!store){
	hcreate(MAX_ENTRY);
	store=calloc(MAX_ENTRY,sizeof(STORE_T));
	nstore=0;
	nused=0;
    }
    char *fn=find_file(config_file);
    if(!(fd=fopen(fn,"r"))){
	error("File %s doesn't exist\n",fn);
    }
#define MAXLN 4096
    char ssline[MAXLN];
    ssline[0]='\0';//stores the available line.
    char line[MAXLN];
    while (fgets(line, MAXLN, fd)){
	sline=squeeze(line);
	//lines ended with \ will continue on next line
	if(!sline || sline[0]=='\0')
	    continue;
	if(sline){
	    if((strlen(sline)+strlen(ssline))>=MAXLN){
		error("Please make ssline longer\n");
	    }
	    strcat(ssline, sline);
	}
	if(ssline[strlen(ssline)-1]=='\\'){
	    ssline[strlen(ssline)-1]=' ';
	    continue;
	}
	var=strtok(ssline, "=");
	value=strtok(NULL, "=");
	strtrim(&var);
	strtrim(&value);
	if(strtok(NULL,"=") || !var || strlen(var)==0){
	    error("Line '%s' is invalid\n",line);
	}else if(var && !strcmp(var,"include")){
	    //info("Opening embeded config file %s\n",value);
	    char *embeded=strextract(value);
	    if(embeded){
		print_override=0;
		open_config(embeded,prefix,protect);
		print_override=1;
		free(embeded);
	    }
	}else{
#if COMPATIBILITY == 1	    
	    /*
	      Compatibility mode: rename old key names to new key names. Will
	      remove in the future.
	     */
	    RENAME(atm.zadeg, sim.zadeg);
	    RENAME(powfs.msa, powfs.order);
	    RENAME(fit.tik_cstr, fit.tikcr);
	    RENAME(tomo.tik_cstr, tomo.tikcr);
	    RENAME(tomo.split_wt, tomo.ahst_wt);
	    RENAME(tomo.split_idealngs, tomo.ahst_idealngs);
	    RENAME(tomo.split_rtt, tomo.ahst_rtt);
	    RENAME(evl.wvl, evl.psfwvl);
	    RENAME(cn2.nhtrecon, cn2.nhtomo);
#endif
	    if(prefix){
		entry.key=stradd(prefix,var,NULL);
	    }else{
		entry.key=strdup(var);
	    }
	    store[nstore].key=entry.key;
	    if(value && strlen(value)>0){
		store[nstore].data=strdup(value);
	    }else{
		store[nstore].data=NULL;
	    }
	    store[nstore].protect=protect;
	    entry.data=(void*)nstore;//record entry.
	    if((entryfind=hsearch(entry,FIND))){
		//same key found
		long istore=(long)(entryfind->data);
		if(store[istore].protect && !protect){
		    info2("%s is protected, has {%s}, will not override by {%s}\n",
			  (char*)entry.key, (char*)store[istore].data,
			  (char*)store[nstore].data);
		    free(store[nstore].key);
		    free(store[nstore].data);
		}else{
		    if(print_override && 
		       (((store[istore].data==NULL || store[nstore].data==NULL)
			 &&(store[istore].data != store[nstore].data))||
			((store[istore].data!=NULL && store[nstore].data!=NULL)
			 &&strcmp((char*)store[istore].data,(char*)store[nstore].data)))){
			info2("Overriding %s:\t{%s}-->{%s}\n", 
			      entry.key, (char*)store[istore].data,
			      (char*)store[nstore].data);
		    }
		    //free old value
		    if(store[istore].data)
			free(store[istore].data);
		    //copy pointer of new value.
		    store[istore].data=store[nstore].data;
		    store[istore].protect=store[nstore].protect;
		    //free key.
		    free(store[nstore].key);
		}
		store[nstore].key=NULL;
		store[nstore].data=NULL;
	    }else{
		//new key
		entryfind=hsearch(entry,ENTER);
		recordcount++;
		nstore++;
		if(nstore>MAX_ENTRY-2){
		    MAX_ENTRY*=2;
		    store=realloc(store,MAX_ENTRY*sizeof(STORE_T));
		}
	    }
	}
	ssline[0]='\0';
    }
    fclose(fd);
    info2("loaded %3d new records from %s\n",recordcount,fn);
    free(fn);
#undef MAXLN
}
/**
   Get the record number of a key.
 */
static long getrecord(char *key, int mark){
    long irecord;
    ENTRY entry, *entryfind;
    strtrim(&key);
    entry.key=key;
    if((entryfind=hsearch(entry,FIND))){
	irecord=(long)entryfind->data;
	if(mark){
	    if(store[irecord].count){
		error("This record %s is already read\n",key);
	    }
	    store[irecord].count++;//record read
	    nused++;
	}
    }else{
	irecord=-1;
	print_file("change.log");
	error("Record %s not found\n",key);
    }
    return irecord;
}
/**
   Check whether a have a record of a key.
 */
int readcfg_peek(const char *format,...){
    //Check whether key exists
    format2key;
    long irecord=getrecord(key, 0);
    if(irecord==-1){
	return 0;
    }else{
	return 1;
    }
}
/**
   Check the size of an array input
*/
int readcfg_peak_n(const char *format, ...){
    format2key;
    long irecord=getrecord(key, 0);
    const char *sdata=store[irecord].data;
    const char *startptr=strchr(sdata,'[');
    const char *endptr=strchr(sdata,']');
    const char *quote=strchr(startptr,'"');
    int count=0;
    if(quote && quote<endptr){//this is string array
	char **ret=NULL;
	count=readstr_strarr(&ret, 0, sdata);
	for(int i=0; i<count; i++){
	    free(ret[i]);
	}
	free(ret);
    }else{//this is numerical array
	double *ret;
	count=readstr_numarr((void**)&ret, 0, T_DBL, sdata);
	free(ret);
    }
    return count;
}
/**
   Check whether the record is overriden by user supplied conf files.
*/
int readcfg_override(const char *format,...){
    //Check whether key exists
    format2key;
    long irecord=getrecord(key, 0);
    if(irecord==-1){
	return 0;
    }else{
	return store[irecord].protect;
    }
}
/**
   Obtain a string value from the key.
 */
char *readcfg_str(const char *format,...){
    format2key;
    char *data;
    long irecord=getrecord(key, 1);
    if(irecord!=-1){
	const char *sdata=store[irecord].data;
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
   //Read str array.
    format2key;
    *res=NULL;//initialize
    long irecord=getrecord(key, 1);
    if(irecord==-1){//record not found.
	error("key '%s' not found\n", key);
	return 0;
    }else{
	const char *sdata=store[irecord].data;
	if(!sdata){//record is empty.
	    return 0;
	}
	return readstr_strarr(res, 0, sdata);
    }
}
/**
   Obtain a string array value from the key. revised on 2010-10-16 to be more
   strict on the entry to avoid mistakes.
 */
int readstr_strarr(char ***res, int len, const char *sdata){
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
	info("sdata[0]=%c, sdataend[0]=%c\n",sdata[0], sdataend[0]);
	warning2("Entry {%s} should start with [ and end with ]\n", sdata);
	sdata2=sdata;
    }else{
	sdata2=sdata+1;
    }
    //sdataend--;
    //find each string.
    while(sdata2<sdataend && sdata2[0]==' '){
	sdata2++;
    }
    while(sdata2<sdataend){
	if(sdata2[0]!='"'){
	    error("Unable to parse {%s} for str array\n", sdata);
	}
	const char *sdata3=sdata2+1;
	const char *sdata4=strchr(sdata3,'"');
	if(!sdata4) error("Unmatched ""\n");
	if(sdata4>sdata3){
	    if(!len && count>=maxcount){
		maxcount*=2;
		*res=realloc(*res,sizeof(char*)*maxcount);
	    }
	    (*res)[count]=mystrndup(sdata3, sdata4-sdata3);
	}else{
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
   Read integer array
*/
int readcfg_intarr(int **ret, const char *format,...){
    format2key;
    return readstr_numarr((void**)ret, 0, T_INT, store[getrecord(key, 1)].data);
}

/**
   Read double array
*/
int readcfg_dblarr(double **ret, const char *format,...){
    format2key;
    return readstr_numarr((void**)ret, 0,T_DBL, store[getrecord(key, 1)].data);
}
/**
   Read integer array of len elements
*/
void readcfg_strarr_n(char ***ret, int len, const char *format,...){
    format2key;
    int len2;
    if(len!=(len2=readstr_strarr((char***)ret, len, store[getrecord(key, 1)].data))){
	error("%s: Need %d, got %d integers\n", key, len, len2);
    }
}
/**
   Read integer array of len elements
*/
void readcfg_intarr_n(int **ret, int len, const char *format,...){
    format2key;
    int len2;
    if(len!=(len2=readstr_numarr((void**)ret, len, T_INT, store[getrecord(key, 1)].data))){
	error("%s: Need %d, got %d integers\n", key, len, len2);
    }
}
/**
   Read integer array of maximum of len elements
*/
void readcfg_intarr_nmax(int **ret, int len, const char *format,...){
    format2key;
    int len2=readstr_numarr((void**)ret, len, T_INT, store[getrecord(key, 1)].data);
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
    if(len!=(len2=readstr_numarr((void**)ret, len, T_DBL, store[getrecord(key, 1)].data))){
	error("%s: Need %d, got %d double\n", key, len, len2);
    }
}
/**
   Read double array of len elements
*/
void readcfg_dblarr_nmax(double **ret, int len, const char *format,...){
    format2key;
    int len2=readstr_numarr((void**)ret, len, T_DBL, store[getrecord(key, 1)].data);
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
    return (int)readstr_num(store[getrecord(key, 1)].data, NULL);
}
/**
   Read double
*/
double readcfg_dbl(const char *format,...){
    format2key;
    return readstr_num(store[getrecord(key, 1)].data, NULL);
}
/**
   Read in a number from the value string. Will interpret * and / operators. *endptr0 will
   be updated to point to the next valid entry, or at separator like coma
   (spaced are skipped).  */
double readstr_num(const char *data, char **endptr0){
    if(!data || strlen(data)==0){
	error("Unable to parse {%s} for a number\n", data);
	return 0;
    }
    char *endptr;
    double res=strtod(data, &endptr);
    if(data==endptr){
	error("Unable to parse {%s} for a number\n", data);
	return 0;
    }
    while(isspace((int)endptr[0])) endptr++;
    while(endptr[0]=='/' || endptr[0]=='*'){
	int power=1;
	if(endptr[0]=='/'){
	    power=-1;
	}
	endptr++;
	while(isspace((int)endptr[0])) endptr++;
	data=endptr;
	double tmp=strtod(data, &endptr);
	if(data==endptr){
	    error("Failed to parse {%s} for a number\n", data);
	}
	if(power==1){
	    res*=tmp;
	}else{
	    res/=tmp;
	}
	while(isspace((int)endptr[0])) endptr++;
    }
    if(endptr0){
	*endptr0=endptr;
    }
    return res;
}
/**
   Read numerical array from a string. if len is nonzero, *ret should be already allocated.
*/
int readstr_numarr(void **ret, int len, int type, const char *data){
    if(!data || strlen(data)==0){
	return 0;
    }
    size_t nmax=10;
    size_t size=0;//size of each number
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
    //process possible numbers before the array.
    if(strchr(startptr,'[')){//there is indeed '['
	while(startptr[0]!='['){
	    double fact1=strtod(startptr, &endptr);//get the number
	    if(startptr==endptr){
		error("Invalid entry to parse for numerical array: {%s}\n", data);
	    }else{//valid number
		if(power==1){
		    fact*=fact1;
		}else{
		    fact/=fact1;
		}
		while(isspace((int)endptr[0])) endptr++;
		if(endptr[0]=='/'){
		    power=-1;
		}else if(endptr[0]=='*'){
		    power=1;
		}else{
		    error("Invalid entry to parse for numerical array: {%s}\n", data);
		}
		startptr=endptr+1;
	    }
	}
	if(startptr[0]!='['){
	    error("Invalid entry to parse for numerical array: {%s}\n", data);
	}
	startptr++;//points to the beginning of the array in []
	/*process possible numbers after the array. do not use startptr here.*/
    }else{
	warning2("Expecting array: {%s} should start with [\n", data);
    }
    if(strchr(startptr,']')){//there is indeed ']'
	endptr=strchr(startptr,']')+1;
	while(isspace((int)endptr[0])) endptr++;
	while(endptr[0]=='/' || endptr[0]=='*'){
	    int power2=1;
	    if(endptr[0]=='/'){
		power2=-1;
	    }else{
		power2=1;
	    }
	    endptr++;
	    while(isspace((int)endptr[0])) endptr++;
	    startptr2=endptr;
	    double fact2=strtod(startptr2, &endptr);
	    if(startptr2==endptr){
		error("Invalid entry to parse for numerical array: {%s}\n", data);
	    }
	    while(isspace((int)endptr[0])) endptr++;
	    if(power2==1){
		fact*=fact2;
	    }else{
		fact/=fact2;
	    }
	}
	if(endptr[0]!='\0'){
	    error("There is garbage in the end of the string: {%s}\n", data);
	}
    }
    int count=0;
    //Read in the array
    while(startptr[0]!=']' && startptr[0]!='\0'){
	if(count>=nmax){
	    if(len){
		error("We want %d numbers, but more are supplied\n", len);
	    }else{
		nmax*=2;
		*ret=realloc(*ret, size*nmax);
	    }
	}
	//parse the string for a floating point number. 
	double res=readstr_num(startptr, &endptr);

	startptr=endptr;
	//apply the factors appear before or after []
	if(power==1){
	    res=fact*res;
	}else{
	    res=fact/res;
	}
	//assign the value to appropriate array. convert to int if necessary.
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
	//Skip the number separators.
	while(startptr[0]==' '){
	    startptr++;
	}
    }
    if(!len){
	if(count>0){
	    *ret=realloc(*ret, size*count);
	}else{
	    free(*ret);
	    *ret=NULL;
	}
    }
    return count;
}


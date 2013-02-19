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

#include "../sys/sys.h"
#include "mathmisc.h"
#include "readstr.h"
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
#define IGNORE(old)				\
    if(!strcmp(var,#old)){			\
	ssline[0]='\0';				\
	continue;				\
    }
#else
#define RENAME(old,new) /*do nothing. */
#define IGNORE(old)
#endif

static int MAX_ENTRY=1000;
typedef struct STORE_T{
    char *key;
    char *data;
    long protect;
    long count;
}STORE_T;
/*Keep record of all the strings so that we can check whether they have all been used. */
static STORE_T *store=NULL;/*store the pointers. */
static long nstore=0;/*number of total records */
static long nused=0;/*number of read records */
static int print_override=1;

#define STRICT 1
/*
   We only use one hash table in the beginning of code. so the non-reentrant hsearch is fine.
*/

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
    while(is_space((*str)[iend]) && iend>=0){
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
	    const char *pre;
	    for(int i=0; i<nstore; i++){
		if(store[i].protect){
		    pre="";
		}else{
		    pre="#";//comment out default parameters
		}
		if(store[i].data){
		    if(strcmp(store[i].data, "ignore")){
			fprintf(fp,"%s%s=%s\n",pre,store[i].key,store[i].data);
		    }
		}else{
		    fprintf(fp,"%s%s=\n",pre,store[i].key);
		}
	    }
	    fclose(fp);
	}
	for(int i=0; i<nstore; i++){
	    if(store[i].count==0){
		//print_file("change.log");
		error("key \"%s\" is not recognized, value is %s\n", store[i].key,store[i].data);
	    }else if(store[i].count!=1){
		/*this should not happen. */
		error("Key %s is used %ld times\n", store[i].key,store[i].count);
	    }
#if defined(__linux__)
	    /*Linux doesn't free the keys. but Mac does. */
	    free(store[i].key);/*free key; */
#endif
	    if(store[i].data){
		free(store[i].data);/*free value */
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
    int countnew=0;
    int countold=0;
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
#define MAXLN 40960
    char ssline[MAXLN];
    ssline[0]='\0';/*stores the available line. */
    char line[MAXLN];
    while (fgets(line, MAXLN, fd)){
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
	if(!index(ssline, '=')){
	    if(check_suffix(ssline, ".conf")){
		char *embeded=strextract(ssline);
		if(embeded){
		    open_config(embeded, prefix, protect);
		    free(embeded);
		}
	    }else{
		error("Input (%s) is not valid\n", ssline);
	    }
	    ssline[0]='\0';
	    continue;
	}
	var=strtok(ssline, "=");
	value=strtok(NULL, "=");
	strtrim(&var);
	strtrim(&value);
	if(!strcmp(var,"path") || !strcmp(var, "PATH")){
	    addpath(value);
	}else if(strtok(NULL,"=") || !var || strlen(var)==0){
	    error("Line '%s' is invalid\n",line);
	}else if(var && !strcmp(var,"include")){
	    /*info("Opening embeded config file %s\n",value); */
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
	    RENAME(atm.wsdeg, atm.wddeg);
	    RENAME(atm.zadeg, sim.zadeg);
	    RENAME(powfs.msa, powfs.order);
	    RENAME(fit.tik_cstr, fit.tikcr);
	    RENAME(tomo.tik_cstr, tomo.tikcr);
	    RENAME(tomo.split_wt, tomo.ahst_wt);
	    RENAME(tomo.split_idealngs, tomo.ahst_idealngs);
	    RENAME(tomo.split_rtt, tomo.ahst_ttr);
	    RENAME(tomo.ahst_rtt, tomo.ahst_ttr);
	    RENAME(evl.psfwvl, evl.wvl);
	    RENAME(cn2.nhtrecon, cn2.nhtomo);
	    /*Added on 2011-04-28 */
	    RENAME(dbg.noatm, sim.noatm);
	    RENAME(dbg.fitonly, sim.fitonly);
	    RENAME(dbg.evlol, sim.evlol);
	    RENAME(sim.fitonly, sim.idealfit);
	    RENAME(evl.ht, evl.hs);
	    RENAME(sim.recon, recon.alg);
	    RENAME(sim.glao, recon.glao);
	    RENAME(dbg.parallel, sim.parallel);
	    RENAME(tomo.split, recon.split);
	    RENAME(llt.fn, llt.fnprof);
	    RENAME(aper.dx, evl.dx);
	    IGNORE(sim.servotype_hi);
	    IGNORE(sim.servotype_lo);
	    IGNORE(sim.epfocus);
	    IGNORE(dbg.dxonedge);
	    RENAME(sim.gtypeII_lo, sim.eplo);
	    RENAME(sim.epngs, sim.eplo);
	    RENAME(sim.apngs, sim.aplo);
	    RENAME(dbg.splitlrt, tomo.splitlrt);
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
	    entry.data=(void*)nstore;/*record entry. */
	    if(store[nstore].data && !strcmp(store[nstore].data, "ignore")){
		store[nstore].count=1;/*Mark it as already consumed. */
	    }else{
		store[nstore].count=0;
	    }
	    if((entryfind=hsearch(entry,FIND))){
		/*same key found */
		long istore=(long)(entryfind->data);
		if(store[istore].protect && !protect){
		    info2("%s={%s} is protected. Will not be overriden by {%s}\n",
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
		    /*free old value */
		    if(store[istore].data)
			free(store[istore].data);
		    /*copy pointer of new value. */
		    store[istore].data=store[nstore].data;
		    store[istore].protect=store[nstore].protect;
		    store[istore].count=store[nstore].count;
		    /*free key. */
		    free(store[nstore].key);
		}
		store[nstore].key=NULL;
		store[nstore].data=NULL;
		countold++;
	    }else{
		/*new key */
		entryfind=hsearch(entry,ENTER);
		countnew++;
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
    info2("loaded %3d (%3d new) records from %s\n",countnew+countold,countnew,fn);
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
	    store[irecord].count++;/*record read */
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
    /*Check whether key exists */
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
int readcfg_peek_n(const char *format, ...){
    format2key;
    long irecord=getrecord(key, 0);
    const char *sdata=store[irecord].data;
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
   /*Read str array. */
    format2key;
    *res=NULL;/*initialize */
    long irecord=getrecord(key, 1);
    if(irecord==-1){/*record not found. */
	error("key '%s' not found\n", key);
	return 0;
    }else{
	const char *sdata=store[irecord].data;
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
    return readstr_numarr((void**)ret, 0,NULL,NULL, T_INT, store[getrecord(key, 1)].data);
}

/**
   Read double array
*/
int readcfg_dblarr(double **ret, const char *format,...){
    format2key;
    return readstr_numarr((void**)ret, 0,NULL,NULL,T_DBL, store[getrecord(key, 1)].data);
}
/**
   Read as a dmat. It can be a file name or an array.
 */
dmat *readcfg_dmat(const char *format,...){
    format2key;
    double *val=NULL;
    int nx, ny;
    char *str=store[getrecord(key, 1)].data;
    if((str[0]<='Z' && str[0]>='A')
       || (str[0]<='z' && str[0]>='a')
       ||(str[0]=='"' || str[0]=='\'')){
	char *fn=strextract(str);
	return dread("%s", fn);
	free(fn);
    }else{
        double **pval=&val;
	readstr_numarr((void**)pval, 0, &nx, &ny,T_DBL, str);
	return dnew_data(nx, ny, val);
    }
}
/**
   Read string array of len elements
*/
void readcfg_strarr_n(char ***ret, int len, const char *format,...){
    format2key;
    int len2;
    if(len!=(len2=readstr_strarr((char***)ret, len, store[getrecord(key, 1)].data))){
	error("%s: Require %d elements, but got %d\n", key, len, len2);
    }
}
/**
   Read str array of upto len elements
*/
void readcfg_strarr_nmax(char ***ret, int len, const char *format,...){
    format2key;
    int len2=readstr_strarr((char***)ret, len, store[getrecord(key, 1)].data);
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
    if(len!=(len2=readstr_numarr((void**)ret, len, NULL,NULL,T_INT, store[getrecord(key, 1)].data))){
	error("%s: Need %d, got %d integers\n", key, len, len2);
    }
}
/**
   Read integer array of maximum of len elements
*/
void readcfg_intarr_nmax(int **ret, int len, const char *format,...){
    format2key;
    int len2=readstr_numarr((void**)ret, len,NULL,NULL, T_INT, store[getrecord(key, 1)].data);
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
    if(len!=(len2=readstr_numarr((void**)ret, len,NULL,NULL, T_DBL, store[getrecord(key, 1)].data))){
	error("%s: Need %d, got %d double\n", key, len, len2);
    }
}
/**
   Read double array of len elements
*/
void readcfg_dblarr_nmax(double **ret, int len, const char *format,...){
    format2key;
    int len2=readstr_numarr((void**)ret, len, NULL,NULL,T_DBL, store[getrecord(key, 1)].data);
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

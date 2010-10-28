/*
  Copyright 2009,2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
    while((*str)[0]!='\0' && isspace((*str)[0])) (*str)++;
    iend=strlen(*str)-1;
    while(isspace((*str)[iend]) && iend>=0){
	(*str)[iend]='\0';
	iend--;
    }
    if((*str)[0]=='\0') *str=NULL;
}
/**
   remove "" or '' around a string from config file. Error if quotes are not
found.  */
static char* strextract(const char *data){
    if(!data || strlen(data)<=2) return NULL;
    if((data[0]!='"' && data[0]!='\'') || 
       (data[strlen(data)-1]!='"' && data[strlen(data)-1]!='\'')){
	error("Record is (%s). \nThis is not a string constant. "
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
	while(nread>=0 && (isspace(line[nread])||line[nread]==';')){
	    line[nread]='\0';
	    nread--;
	}
	/*Remove leading spaces*/
	sline=line;
	while(isspace(sline[0])) sline++;
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
void open_config(const char* config_file, long protect){
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
    const size_t sslineln=4096;
    char ssline[sslineln];
    ssline[0]='\0';//stores the available line.
    const int linemax=4096;
    char line[linemax];
    while (fgets(line, linemax, fd)){
	sline=squeeze(line);
	//lines ended with \ will continue on next line
	if(!sline || sline[0]=='\0')
	    continue;
	if(sline){
	    if((strlen(sline)+strlen(ssline))>=sslineln){
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
		open_config(embeded,protect);
		print_override=1;
		free(embeded);
	    }
	}else{
	    store[nstore].key=entry.key=strdup(var);
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
}
/**
   Get the record number of a key.
 */
static long getrecord(char *key){
    long irecord;
    ENTRY entry, *entryfind;
    strtrim(&key);
    entry.key=key;
    if((entryfind=hsearch(entry,FIND))){
	irecord=(long)entryfind->data;
	if(store[irecord].count){
	    error("This record %s is already read\n",key);
	}
	store[irecord].count++;//record read
    }else{
	irecord=-1;
	print_file("change.log");
	error("Record %s not found\n",key);
    }
    nused++;
    return irecord;
}
/**
   Check whether a have a record of a key.
 */
int readcfg_peek(const char *format,...){
    //Check whether key exists
    format2key;
    ENTRY entry;
    char *key2=key;
    strtrim(&key2);
    entry.key=key2;
    if(hsearch(entry,FIND)){
	return 1;
    }else{
	return 0;
    }
}
/**
   Obtain a string value from the key.
 */
char *readcfg_str(const char *format,...){
    format2key;
    char *data;
    long irecord=getrecord(key);
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
/**
   Obtain a string array value from the key. revised on 2010-10-16 to be more
   strict on the entry to avoid mistakes.
 */
int readcfg_strarr(char ***res, const char *format,...){
    //Read str array.
    format2key;
    long irecord=getrecord(key);
    int count=0, maxcount=5;
    *res=calloc(maxcount,sizeof(char*));
    if(irecord!=-1){
	const char *sdata=store[irecord].data;
	const char *sdataend=sdata+strlen(sdata)-1;
	const char *sdata2, *sdata3, *sdata4;
	if(sdata[0]!='[' || sdataend[0]!=']'){
	    info("sdata[0]=%c, sdataend[0]=%c\n",sdata[0], sdataend[0]);
	    error("key %s: Entry (%s) should start with [ and end with ]\n",key, sdata);
	}
	sdata2=sdata+1;
	sdataend--;
	//find each string.
	while(sdata2<sdataend && (sdata2[0]==','||sdata2[0]==';'||!isgraph(sdata2[0]))){
	    sdata2++;
	}
	while(sdata2<sdataend){
	    if(sdata2[0]!='"'){
		error("Unable to parse (%s) for str array\n", sdata);
	    }
	    sdata3=sdata2+1;
	    sdata4=strchr(sdata3,'"');
	    if(!sdata4) error("Unmatched ""\n");
	    if(sdata4>sdata3){
		if(count>=maxcount){
		    maxcount*=2;
		    *res=realloc(*res,sizeof(char*)*maxcount);
		}
		(*res)[count]=strndup(sdata3,sdata4-sdata3);
	    }else{
		(*res)[count]=NULL;
	    }
	    count++;
	    sdata2=sdata4+1;
	    while(sdata2<sdataend && (sdata2[0]==','||sdata2[0]==';'||!isgraph(sdata2[0]))){
		sdata2++;
	    }
	}
	/*if(*sdata2!=']'){
	    error("Entry (%s) does not have the right format\n", sdata);
	    }*/
	*res=realloc(*res,sizeof(char*)*count);
    }else{
	error("key '%s' not found\n", key);
	*res=NULL;
    }
    return count;
}
#define TYPE int
#define TYPEFUN1 readcfg_intarr
#define TYPEFUN2 readcfg_int
#define TYPENAME "int"
#define TYPECFUN(A,B) strtol(A,B,10)
#include "readcfg_numarr.c"
#undef TYPE
#undef TYPEFUN1
#undef TYPEFUN2
#undef TYPENAME
#undef TYPECFUN

#define TYPE double
#define TYPEFUN1 readcfg_dblarr
#define TYPEFUN2 readcfg_dbl
#define TYPENAME "double"
#define TYPECFUN(A,B) strtod(A,B)
#include "readcfg_numarr.c"
#undef TYPE
#undef TYPEFUN1
#undef TYPEFUN2
#undef TYPENAME
#undef TYPECFUN

#ifdef TEST
int main(int argc, char*argv[]){
    char *config_file;
    char *file;
    int nwfs,iwfs;
    int *wfs;
    int nwvl,iwvl;
    double *wvl;
    if(argc>1){
	config_file=argv[1];
    }else{
	config_file=NULL;
    }
    fprintf(stderr, "%p=%p\n", &wfs,wfs);
    
    open_config(config_file,0
    file=readcfg_str("file");
    nwfs=readcfg_intarr("wfs", &wfs);
    printf("file is %s\n",file);
    printf("nwfs is %d, wfs is ",nwfs);
    for(iwfs=0; iwfs<nwfs; iwfs++){
	printf("%d ",wfs[iwfs]);
    }
    printf("\n");
    
    nwvl=readcfg_dblarr("wvl",&wvl);
    printf("nwvl is %d\n", nwvl);
    for(iwvl=0; iwvl<nwvl; iwvl++){
	printf("%f ",wvl[iwvl]);
    }
    printf("atm.r0=%f\n",readcfg_dbl("atm.r0"));
    printf("deg=%ld\n",readcfg_int("deg"));
    close_config();
    return 0;
}

#endif

/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
  then user can query tokens using read_int, read_real, or read_array
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
DEF_ENV_FLAG(PERMISSIVE, 0)

static void* MROOT=NULL;
static char* default_config=NULL;
typedef struct STORE_T{
	char* key;    //Name of the entry
	char* data;   //Value of the entry
	int index;	  //Input order.
	int priority;//Priority of the entry. -1: a replacement entry. 0: default, larger: higher priority
	int flag;    //-1: a replacement entry. 0: not used. 1: used once. Others: error happens.
	int prefix;  //whether this is prefix'ed entry.
}STORE_T;
static int key_cmp(const void* a, const void* b){
	return strcmp(((STORE_T*)a)->key, ((STORE_T*)b)->key);
}

/*Keep record of all the strings so that we can check whether they have all been used. */
static long nstore=0;/*number of total records */
static long nused=0;/*number of read records */
static int changes_loaded=0;//whether changes.conf has been loaded.
#define STRICT 1

/**
   trim the spaces before and after string. str is modified to skip leading spaces.*/
static void strtrim(char** str){
	if(!*str) return;
	int iend;
	/*turn non-printable characters, coma, and semicolumn, to space */
	for(char* tmp=*str; !is_end(*tmp); tmp++){
		if(!isgraph((int)*tmp)||isspace(*tmp)){
			*tmp=' ';
		}
	}
	/*remove leading spaces. */
	while(!is_end(**str)&&isspace((*str)[0])) (*str)++;
	iend=strlen(*str)-1;
	/*remove tailing spaces. */
	while((isspace((*str)[iend])||(*str)[iend]==';')&&iend>=0){
		(*str)[iend]='\0';
		iend--;
	}
	/*remove continuous spaces */
	char* tmp2=*str;
	int issp=0;
	for(char* tmp=*str; !is_end(*tmp); tmp++){
		if(*tmp=='['){
			issp=2; //skip space after [
		} else if(*tmp==' '){
			if(issp){//skip multiple space
				continue;
			} else{
				issp=1;
			}
		} else if(*tmp==']'){
			if(issp==1){//remove space before ]
				tmp2--;
			} else{
				issp=2;//skip space after ]
			}
		} else{
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
static char* strextract(const char* data){
	const char *start=data, *end=0;
	trim_string(&start, &end);
	return mystrndup(start, end-start);
}
/**
   Remove comment (after #), leading spaces, and trailing line break and spaces from a
   string.
 */
static char* squeeze(char* line){
	char* sline, * comment;
	if(!line){
		sline=NULL;
	} else{
	/*Remove comment*/
		comment=strchr(line, '#');
		if(comment){
			comment[0]='\0';
		}
		/*Remove trailing linebreak, semi-colon, and spaces*/
		int nread=strlen(line)-1;
		while(nread>=0&&(isspace(line[nread])||is_end(line[nread]))){
			line[nread]='\0';
			nread--;
		}
		/*Remove leading spaces*/
		sline=line;
		while(isspace(sline[0])) sline++;
		if(is_end(sline[0]))  sline=NULL;
	}
	return sline;
}
static FILE* fpout=NULL;
static void print_key(const void* key, VISIT which, int level){
	const STORE_T* store=*((const STORE_T**)key);
	(void)level;
	if((which==leaf||which==postorder) && store->flag!=-1){
		if(fpout){
			fprintf(fpout, "%s%s=%s\n", store->priority>0?"":"#", store->key, (store->data&&strcmp(store->data, "ignore"))?store->data:"");
			
		}
		if(store->flag!=1&&(!store->data||strcmp(store->data, "ignore"))){
			if(store->flag==0){
				if(PERMISSIVE||store->priority>0||store->prefix){
					warning("key \"%s\" is not recognized, value is %s\n", store->key, store->data);
				} else{
					error("key \"%s\" is not recognized, value is %s. Set env MAOS_PERMISSIVE=1 to ignore the error.\n", store->key, store->data);
				}
			} else if(store->flag>1){
				if(PERMISSIVE){
					warning("Key %s is used %d times\n", store->key, store->flag);
				} else{
					error("Key %s is used %d times. Set env MAOS_PERMISSIVE=1 to ignore the error.\n", store->key, store->flag);
				}
			}
		}
	}
}
static void delete_leaf(const void* key, VISIT which, int level){
	(void)level;
	if(which==leaf){
		STORE_T* store=*((STORE_T**)key);
		tdelete(store, &MROOT, key_cmp);
		free(store->key);
		free(store->data);
		free(store);
	}
}
void erase_config(){
	while(MROOT){
		twalk(MROOT, delete_leaf);
	}
	nused=0;
	nstore=0;
	changes_loaded=0;
	free(default_config); default_config=NULL;
}
/**
   Save all configs to file and check for unused config options.
 */
void close_config(const char* format, ...){
	format2fn;
	const char* fnout=format?fn:NULL;
	if(MROOT){
		dbg("Used %ld of %ld supplied keys\n", nused, nstore);
		if(fnout&&strlen(fnout)>0&&!disable_save) fpout=fopen(fnout, "w");
		if(fpout && default_config) fprintf(fpout, "include=%s\n", default_config);
		twalk(MROOT, print_key);
		if(fpout){ fclose(fpout); fpout=NULL; }
	}
	erase_config();
}

/**
   Start the read config process by opening .conf files and fill the entries in
   a hash table. A key has a priority or 0 or higher. A new key with same or
   higher priority can override previous entry. priority is 0 for default configurations. -1 indicates optional configuration file that fails silently.
 */
void open_config_full(const char* config_in, /**<[in]The .conf file to read*/
	const char* prefix,    /**<[in]if not NULL, prefix the key with this.*/
	int priority,    /**<[in]Priorities of keys.*/
	const int index		   /**<[in]Only when prefix is set. Initial entry index. */
){
	if(!config_in) return;
	if(!changes_loaded){//on first call, load changes.conf
		changes_loaded=1;
		open_config_full("changes.conf", NULL, -1, 0);
	}
	FILE* fd=NULL;
	char* config_file=NULL;
	char* config_dir=NULL;//directory for config_file
	static int current_index=0;
	if(check_suffix(config_in, ".conf")){
		config_file=search_file(config_in, 0);
		if(!config_file||!(fd=fopen(config_file, "r"))){
			if(priority<0){
				dbg("Cannot open file %s for reading.\n", config_in);
			}else if(!prefix){
				error("Cannot open file %s for reading.\n", config_in);
			}else if(default_config){
				warning("Cannot open file %s for reading. Ignored (prefix=%s)\n", config_in, prefix);
			}
		}
		config_dir=mydirname(config_file);
		if(priority!=-1 && !default_config){//used in close_config to reproduce the simulation
			default_config=strdup(config_file);//use full path to avoid recursive inclusion when maos_recent.conf is used as the initial file.
			addpath2(-1, "%s", config_dir);
			if(priority==1) priority=0;//default config has priority=0
		}
	} else{
		config_file=strdup(config_in);
		parse_argopt(config_file, NULL);
		char* end;
		//Remove trailing space
		for(end=config_file+strlen(config_file); end>=config_file; end--){
			if(isspace((int)*end)) *end='\0';
			else break;
		}
		if(end<config_file){
			warning("config %s is empty\n", config_in);
			free(config_file);
			return;
		}
	}

	char* sline=NULL;
	int countnew=0;
	int countold=0;
	//int countskip=0;
#define MAXLN 10000
	char ssline[MAXLN];
	ssline[0]='\0';/*stores the available line. */
	char line[MAXLN];//line concatenates to ssline if ends with /
	char* config_line=config_file;//contains config lines
	while(1){
		if(fd){/*read from file*/
			if(!fgets(line, MAXLN, fd)) break;
		} else{/*read from string*/
			if(!config_line) break;
			char* p0=strchr(config_line, '\n');
			const int len=p0?(p0-config_line):(MAXLN-1);
			if(len+1>MAXLN){
				error("Input line is too long. Please make MAXLN larger to accomodate.\n");break;
			}else{
				strncpy(line, config_line, len); line[len]='\0';
			}
			if(p0){
				config_line=p0+1;//start of next line
			} else{
				config_line=NULL;//no more
			}
		}
		sline=squeeze(line);
		if(!sline||is_end(sline[0]))/*skip empty lines. */
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
		char* eql=strchr(ssline, '=');
		if(!eql){//no equal sign
			if(strlen(ssline)>2 && ssline[0]=='_'&&ssline[1]=='_'){
				ssline[0]='\0';
				continue;
			}else if(!check_suffix(ssline, ".conf")){
				error("Input (%s) is not valid.\n", ssline);
			}
		}else if(eql==ssline){
			error("Input (%s) should not start with =", ssline);
		}
		int append=0;
		int replace=0;//this indicates a replacement entry
		if(eql){
			if(eql[-1]=='+'){
				append=1;//append to key
				eql[-1]='\0';
			}else if(eql[-1]=='-'){
				append=-1;//remove from key
				eql[-1]='\0';
			}
			if(eql[1]=='='){//equivalent keys:old==new keys
				replace=1;
				eql[1]=' ';
				if(eql[2]=='>'){
					eql[2]=' ';
					replace=2;
				}
			}else if(eql[1]=='>'){//old=>new keys
				replace=2;
				eql[1]=' ';
			}
			eql[0]='\0';
		}
		char* var0=ssline;
		strtrim(&var0);
		const char* var=eql?var0:"include";//key name
		char* value=eql?(eql+1):var0; //key value
		strtrim(&value);

		if(!var||strlen(var)==0){
			error("Line '%s' is invalid\n", line);
		} else if(!strcmp(var, "path")||!strcmp(var, "PATH")){
			char* val2=strextract(value);
			if(config_dir&&val2[0]!='/'&&val2[0]!='~'){//not absolute path, make it relative to config_dir
				char *tmp=val2;
				val2=stradd(config_dir, "/", tmp, NULL);
				free(tmp);
			}
			dbg("addpath %s.\n", val2);
			addpath2(99, "%s", val2);
			free(val2);
		} else if(!strncmp(var, "MAOS_", 5)){
			//This is environment variable.
			if(!value){
				unsetenv(var);//not effective for read_sys_env()
			}else{
				setenv(var, value, 1);
				read_sys_env();
			}
		} else if(!strcmp(var, "include")){
			/*dbg("Opening embeded config file %s\n",value); */
			char* embeded=strextract(value);
			if(embeded){
				if(strcmp(embeded, config_in)){
					open_config_full(embeded, prefix, priority, index);
				}else{
					warning("Recursive inclusion: %s includes %s. Ignored.\n", config_in, embeded);
				}
				free(embeded);
			}
		} else{
			STORE_T* store=mycalloc(1, STORE_T);
			if(prefix){
				store->key=stradd(prefix, var, NULL);
				store->prefix=1;
			} else{
				store->key=strdup(var);
			}

			if(value&&strlen(value)>0){
				store->data=strdup(value);
			} else{
				store->data=NULL;
			}
			if(replace){
				store->flag=-1;
			}else{
				store->flag=0;
			}
			store->priority=priority;
			if(prefix){
				store->index=index;
			}else if(priority!=-1){
				store->index=++current_index;
			}
			void* entryfind=tfind(store, &MROOT, key_cmp);
			if(entryfind && priority!=-1){/*key found, check whether it is renamed.*/
				STORE_T* oldstore=*(STORE_T**)entryfind;
				if(oldstore->flag==-1){//a replacement entry
					if(oldstore->data&&oldstore->data[0]!=0&&oldstore->data[0]!='('){
						if(replace==2) warning("%s has been renamed to %s.\n", store->key, oldstore->data);
						free(store->key);
						store->key=strdup(oldstore->data);
						entryfind=tfind(store, &MROOT, key_cmp);//search again
					}else{
						warning("%s is no longer needed. %s\n", store->key, oldstore->data?oldstore->data:"");
						free(store->key); store->key=NULL;
						entryfind=NULL;
					}
				}
			}
			if(entryfind){
				STORE_T *oldstore=*(STORE_T **)entryfind;
			 	if(append){
					/*append (append=1) or remove (append=-1) new value with/from old value for arrays. 
					both have to start/end with [/] */
					const int nolddata=strlen(oldstore->data);
					const int nnewdata=strlen(store->data);
					if(nolddata>0 && (oldstore->data[0]!='['||oldstore->data[nolddata-1]!=']')){
						error("olddata='%s' should be entirely encapsulated by [].\n", oldstore->data);
					} else if(store->data[0]!='['||store->data[nnewdata-1]!=']'){
						error("newdata='%s' should be entirely encapsulated by [].\n", store->data);
					} else{
						if(append==1){//append
							oldstore->data=myrealloc(oldstore->data, (nolddata+nnewdata), char);
							oldstore->data[nolddata-1]=oldstore->data[nolddata-2]=='['?0:' ';
							strncat(oldstore->data, store->data+1, nnewdata-1);
						}else if(append==-1){//remove
							if(!strncmp(oldstore->data+(nolddata-nnewdata+1), store->data+1, nnewdata-1)){
								char *ptmp=oldstore->data+(nolddata-nnewdata+1);
								while(ptmp>oldstore->data && ptmp[-1]==' ') ptmp--;
								ptmp[0]=']';
								ptmp[1]=0;
							}else{
								error("Cannot remove %s from the end of %s\n", store->data, oldstore->data);
							}
						}else{
							error("Invalid value: append=%d\n", append);
						}
					}
				} else {//check if the values are identical
					if(((oldstore->data==NULL||store->data==NULL)&&(oldstore->data!=store->data))||
							((oldstore->data!=NULL&&store->data!=NULL)&&strcmp(oldstore->data, store->data))){
						if(oldstore->priority>priority || oldstore->index>store->index){//Entry with higher priority or entered layer prevails.
							dbg("Not overriding %-20s\t%10s by %s\n", store->key, oldstore->data, store->data);
						} else{//Replace value if different.
							if(priority>0){//Print if not default
								info("Overriding %-20s\t%10s --> %s\n", store->key, oldstore->data, store->data);
							}
							free(oldstore->data);/*free old value */
							oldstore->data=store->data; store->data=NULL;/*move pointer of new value. */
							oldstore->priority=store->priority;
							oldstore->index=store->index;
						}
					}
				}
				countold++;
			}else if(store->key){/*key may be NULL if it is declared obsolete.*/
				if(!tsearch(store, &MROOT, key_cmp)){
					error("Error inserting to tree\n");
				}
				countnew++;
				if(store->flag!=-1){
					nstore++;
				}
				store=NULL;//consumed by the tree.
			}
			if(store){
				free(store->data);
				free(store->key);
				free(store);
			}
		}
		ssline[0]='\0';
	}
	if(strcmp(config_in, "changes.conf")){
		info("loaded %3d (%3d new) records from '%s'\n", countnew+countold, countnew, fd?config_file:"command line");
	}else{
		dbg("loaded %3d (%3d new) records from '%s'\n", countnew+countold, countnew, fd?config_file:"command line");
	}
	if(fd){
		fclose(fd);
	}
	free(config_file);
	free(config_dir);
#undef MAXLN
}
void open_config(const char *config_in){
	open_config_full(config_in, NULL, 1, 0);
}
void open_config_prefix(const char *config_in, const char *prefix, const char *prekey){
	open_config_full(config_in, prefix, readcfg_peek_priority("%s",prekey), readcfg_peek_index("%s",prekey));
}
/**
   Get the record number of a key.

   Value of mark:
   1: mark the record as read.
   0: do not mark or check record existance
   -1:do not mark but error if record is not found.
 */
static const STORE_T* getrecord(char* key, int mark){
	STORE_T store;
	void* found=NULL;
	strtrim(&key);
	store.key=key;
	if((found=tfind(&store, &MROOT, key_cmp))){
		if(mark>0){
			if((*(STORE_T**)found)->flag>0){
				error("This record %s is already read.\n", key);
			} else if((*(STORE_T **)found)->flag<0){
				error("This record %s should not be read.\n", key);//it is a change record.
			}
			(*(STORE_T**)found)->flag=1;
			nused++;
		}
		const char* data=(*(STORE_T**)found)->data;
		if(!mystrcmp(data, "ref:")){
			//This key references another key's value
			char* key2=mystrdup(data+4);
			const STORE_T* refed=getrecord(key2, -1);
			free(key2);
			return refed;
		}
	} else if(mark){
		warning("Record %s not found\n", key);
	}
	return found?(*(STORE_T**)found):0;
}
static void *getrecord_data(char* key, int mark){
	const STORE_T *record=getrecord(key, mark);
	return record?record->data:0;
}
/**
   Check whether a have a record of a key.
 */
int readcfg_peek(const char* format, ...){
	/*Check whether key exists */
	format2key;
	return getrecord(key, 0)?1:0;
}
/**
   Check the size of an array input
*/
int readcfg_peek_n(const char* format, ...){
	format2key;
	const STORE_T* store=getrecord(key, 0);
	if(!store) return 0;
	const char* sdata=store->data;
	const char* startptr=strchr(sdata, '[');
	const char* endptr=strchr(sdata, ']');
	if(!startptr) startptr=sdata;
	if(!endptr) endptr=startptr+strlen(sdata)-1;
	const char* quote=strchr(startptr, '"');
	int count=0;
	if(quote&&quote<endptr){/*this is string array */
		char** ret=NULL;
		count=readstr_strarr(&ret, 0, 0, key, sdata);
		for(int i=0; i<count; i++){
			free(ret[i]);
		}
		free(ret);
	} else{/*this is numerical array */
		void* ret;
		count=readstr_numarr(&ret, NULL, NULL, 0, 0, M_REAL, key, sdata);
		free(ret);
	}
	return count;
}
/**
 * Ignore an entry and do not warn if it does not exist.
 * */
void readcfg_ignore(const char *format, ...){
	format2key;
	if(getrecord(key,0)){
		getrecord(key, 1);//mark as used.
	}//else: silently ignore non-existant input.
}
/**
   Return priority value of the record
*/
int readcfg_peek_priority(const char *format, ...){
	format2key;
	const STORE_T *store=getrecord(key, 0);
	return store?store->priority:0;
}/**
   Return index value of the record
*/
int readcfg_peek_index(const char *format, ...){
	format2key;
	const STORE_T *store=getrecord(key, 0);
	return store?store->index:0;
}
/**
   Obtain a string value from the key.
 */
char* readcfg_str(const char* format, ...){
	format2key;
	char* data=NULL;
	const STORE_T* store=getrecord(key, 1);
	const char* sdata=store?store->data:0;
	if(sdata&&strlen(sdata)>0){
		data=strextract(sdata);
	}
	return data;
}

/**
   Read as an lmat.
 */
lmat *readcfg_lmat(int n, int relax, const char *format, ...){
	format2key;
	long* val=NULL;
	int nx, ny;
	char *data;
	if(!(data=getrecord_data(key, 1))){
		warning("Record for %s not found, assume 0.\n", key);
		return NULL;
	}
	readstr_numarr((void**)&val, &nx, &ny, n, relax, M_LONG, key, data);
	lmat* res=NULL;
	if(nx && ny){
		res=lnew_do(nx, ny, val, mem_new(val));
	}
	return res;
}
/**
   Read as a dmat. It can be a file name or an array.
 */
dmat* readstr_dmat(int n, ///<[in]Number of elements requested
				int relax, ///<[in]1: allow fewer values and fill the rest
				const char *key, ///<[in] the key that needs the value.
				const char *str///<[in]input
				){
	if(!str){
		return NULL;
	}
	dmat* res=NULL;
	char* fn=strextract(str);
	if(check_suffix(fn, ".bin.gz")||check_suffix(fn, ".fits.gz")||
		check_suffix(fn, ".bin")||check_suffix(fn, ".fits")){
		if(!(res=dread("%s", fn))){
			error("%s: read file %s failed\n", key, fn);
		}
		if(n>0 && PN(res)!=n){
			if((relax==2&& n>0) || (relax==1 && n==1)){//resize and fill remaining values
				long n2=PN(res);
				dresize(res, n, 1);
				for(long i=n2; i<n; i++){
					P(res,i)=P(res,n2-1);
				}
			}else{
				error("%s=%s: need %d values, got %ld instead.\n", key, fn, n, PN(res));
			}
		}
	} else{
		int nx, ny;
		real* val=NULL;
		readstr_numarr((void**)&val, &nx, &ny, n, relax, M_REAL, key, fn);
		if(nx&&ny){
			res=dnew_do(nx, ny, val, mem_new(val));
		}
	}
	free(fn);
	return res;
}

/**
   Read as a dmat. It can be a file name or an array.
 */
dmat* readcfg_dmat(int n, int relax, const char* format, ...){
	format2key;
	return readstr_dmat(n, relax, key, getrecord_data(key, 1));
}

/**
   Read string array of len elements
*/
int readcfg_strarr(char*** ret, int len, int relax, const char* format, ...){
	format2key;
	return readstr_strarr((char ***)ret, len, relax, key, getrecord_data(key, 1));
}

/**
   Read integer array of len elements
*/
int readcfg_intarr(int** ret, int len, int relax, const char* format, ...){
	format2key;
	return readstr_numarr((void**)ret,  NULL, NULL, len, relax,  M_INT, key, getrecord_data(key, 1));

}
/**
   Read integer array of len elements
*/
int readcfg_dblarr(real **ret, int len, int relax, const char *format, ...){
	format2key;
	return readstr_numarr((void **)ret, NULL, NULL, len, relax, M_REAL, key, getrecord_data(key, 1));
}
/**
   Read integer
*/
int readcfg_int(const char* format, ...){
	format2key;
	char* val;
	char* endstr;
	real ans=0;
	if(!(val=getrecord_data(key, 1))){
		warning("Record for %s not found, assume 0.\n", key);
		ans=0;
	}else if(isnan(ans=readstr_num(key, val, &endstr))||endstr[0]!='\0'||(ans-(int)ans)!=0){
		error("Invalid data: %s=%s\n", key, val);
	}
	return (int)ans;
}
/**
   Read real
*/
real readcfg_dbl(const char* format, ...){
	format2key;
	char *data; 
	char *endstr;
	real ans=0;
	if(!(data=getrecord_data(key, 1))){
		warning("Record for %s not found, assume 0.\n", key);
		ans=0;
	}else if(isnan(ans=readstr_num(key, data, &endstr))||endstr[0]!='\0'){
		error("Invalid data: %s=%s\n", key, data);
	}
	return ans;
}

/**
   Read dcell
*/
dcell* readcfg_dcell(const char* format, ...){
	format2key;
	const char* str=getrecord_data(key, 1);
	if(str){
		return dcellread("%s", str);
	} else{
		return NULL;
	}
}

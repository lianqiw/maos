/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#include <sys/stat.h>
#include <sys/types.h>
#include <limits.h>

#include "path.h"
#include "common.h"
/**
   The linked list to store path where we look for files.
*/
typedef struct PATH_T{
    char *path;
    struct PATH_T *next;
}PATH_T;

/**
   \file path.c

   This file contains function managing the searching path for files, behaves as
   PATH in POSIX systemsx */
static PATH_T *PATH=NULL;//privately maintained path to locate config files.

/**
   Add a directory to path.
*/
void addpath(const char*path){
    char *abspath=myabspath(path);
    PATH_T *node=calloc(1,sizeof(PATH_T));
    node->path=abspath;
    node->next=PATH;
    PATH=node;
}
/**
   Remove a directory from path.
 */
void rmpath(const char *path){
    char *abspath=myabspath(path);
    PATH_T *ia,*ib=NULL;
    for(ia=PATH;ia;ia=ia->next){
	if(!strcmp(ia->path,abspath)){
	    if(ib){
		ib->next=ia->next;
	    }else{
		PATH=ia->next;
	    }
	    free(ia->path);
	    break;
	}
	ib=ia;
    }
    if(!ia){
	warning("%s not foun idn PATH\n",path);
    }
    free(abspath);
}
/**
   Print current path
*/
void printpath(void){
    info2("PATH is :\n");
    for(PATH_T *ia=PATH;ia;ia=ia->next){
	info2("%s\n",ia->path);
    }
}
/**
   Empty the path
*/
void freepath(void){
    PATH_T *ib;
    for(PATH_T *ia=PATH;ia;ia=ib){
	ib=ia->next;
	free(ia);
    }
    PATH=NULL;
}
/**
   Free the path at exit.
 */
static __attribute__((destructor)) void deinit(){
    freepath();
}
/**
   Try to find a file in path and return its absolute filename of exist, NULL
otherwise.  */
char *search_file(const char *fn){
    if(!fn) return NULL;
    char *fnout=NULL;
    if(exist(fn)){
	fnout=strdup(fn);
    }else{
	PATH_T *ia;
	char fntmp[PATH_MAX];
	for(ia=PATH;ia;ia=ia->next){
	    strcpy(fntmp,ia->path);
	    if(strlen(ia->path)>0 && fntmp[strlen(ia->path)-1]!='/') 
		strcat(fntmp,"/");
	    strcat(fntmp,fn);
	    if(exist(fntmp)){
		fnout=strdup(fntmp);
		break;
	    }
	}
    }
    return fnout;
}
/**
   Locate a file in path and return its absolute filename. Will emit error if
   not found.
 */
char *find_file(const char *fn){
    char *fnout=search_file(fn);
    if(!fnout || !exist(fnout)){
	info("Looking for %s, found %s\n",fn, fnout);
	printpath();
	error("Unable to find file %s.\n",fn);
	return NULL;
    }
    return fnout;
}
/**
   Find the directory that hold the configuration files. name is "maos" for
maos, "skyc" for skyc.  */
char *find_config(const char *name){
    const char *aos_config_path=getenv("AOS_CONFIG_PATH");
    char *config_path=NULL;
    if(aos_config_path){
	config_path=stradd(aos_config_path,"/",name,NULL);
    }
    if(!exist(config_path) && exist(SRCDIR)){
	/*If not specified, assume it is in the source tree*/
	if(config_path) free(config_path);
	config_path=stradd(SRCDIR,"/config/",name,NULL);
    }
 
    if(!exist(config_path)){
	info("Configuration files not found in %s\n", config_path);
	free(config_path);
	char *cwd=mygetcwd();
	config_path=stradd(cwd,"/config/",name,NULL);
	free(cwd);
    }

    if(!exist(config_path)){
	free(config_path);
	info("Configuration files not found in %s\n", config_path);
	config_path=stradd(HOME,"/.aos/config-",PACKAGE_VERSION,"/",name,NULL);
    }
    if(!exist(config_path)){
	info("Configuration files not found in %s\n", config_path);
#ifdef __CYGWIN__
	error("Please download the configuration files from github.com/lianqiw/maos/downloads and extract to current folder.\n");
#else
	warning("Unable to determine the path to the configuration files.\n");
	warning("Will download a copy from the website and put in %ss/.aos/config\n",HOME);
	char cmd[400];
	char FN[80];
	snprintf(FN,40,"config-%s.tar.bz2", PACKAGE_VERSION);
	snprintf(cmd,400,"wget %s/%s -O %s/%s && tar axvf %s/%s -C %s/.aos/ && rm -rf %s/%s",
		 BASEURL,FN,TEMP,FN,TEMP,FN,HOME,TEMP,FN);
	if(system(cmd)){
	    error("Unable to download the configuration files from the internet "
		  "and extract to %s/.aos/config", HOME);
	}
#endif
    }
    if(!exist(config_path)){
	error("Unable to determine the path to the configuration files.\n");
    }
    return config_path;
}

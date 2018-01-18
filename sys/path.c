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

#include <sys/stat.h>
#include <sys/types.h>


#include <limits.h>
#include "common.h"
#include "process.h"
#include "path.h"
#include "misc.h"
#include "thread.h"
/**
   This file contains function managing the searching path for files, behaves
   like PATH in POSIX systems.*/

/**
   The linked list to store path where we look for files.
*/
typedef struct PATH_T{
    char *path;
    struct PATH_T *next;
}PATH_T;
static PATH_T *PATH=NULL;/*privately maintained path to locate config files. */
PNEW(mutex_path);
/**
   Add a directory to path.
*/
void addpath(const char*path){
    char *abspath=myabspath(path);
    if(!path || !abspath){
	warning2("Path not found: path=%s; abspath=%s; pwd=%s. Ignored.\n", path, abspath,mygetcwd());
	return;
    }
    PATH_T *node=mycalloc(1,PATH_T);
    node->path=abspath;
    LOCK(mutex_path);
    node->next=PATH;
    PATH=node;
    UNLOCK(mutex_path);
}

/**
   Remove a directory from path.
 */
void rmpath(const char *path){
    char *abspath=myabspath(path);
    PATH_T *ia,*ib=NULL;
    LOCK(mutex_path);
    for(ia=PATH;ia;ia=ia->next){
	if(!strcmp(ia->path,abspath)){/*found */
	    if(ib){/*there is parent node */
		ib->next=ia->next;
	    }else{
		PATH=ia->next;
	    }
	    free(ia->path);
	    free(ia);
	    ia=ib;
	}else{
	    ib=ia;
	}
    }
    UNLOCK(mutex_path);
    free(abspath);
}
/**
   Print current path.
*/
void printpath(void){
    info2("PATH is :\n");
    for(PATH_T *ia=PATH;ia;ia=ia->next){
	info2("%s\n",ia->path);
    }
}
/**
   Empty the path.
*/
void freepath(void){
    LOCK(mutex_path);
    for(PATH_T *ia=PATH;ia;ia=PATH){
	PATH=ia->next;
	free(ia->path);
	free(ia);
    }
    UNLOCK(mutex_path);
}
/**
   Try to find a file in path and return its absolute filename if exist, NULL
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
    const char *maos_config_path=getenv("MAOS_CONFIG_PATH");
    char *config_path=NULL;
    if(!exist(config_path) && EXEP[0]){
	/*If not found, try the folder that contains the exe*/
	free(config_path);
	config_path=stradd(EXEP,"/config/",name,NULL);
    }
    if(!exist(config_path) && maos_config_path){
	config_path=stradd(maos_config_path,"/config/",name,NULL);
	if(!exist(config_path)){
	    free(config_path);
	    config_path=stradd(maos_config_path,"/",name,NULL);
	}
    }
    if(!exist(config_path) && exist(SRCDIR)){
	/*If not specified, assume it is in the source tree*/
	free(config_path);
	config_path=stradd(SRCDIR,"/config/",name,NULL);
    }

    if(!exist(config_path) && HOME){
	/*If not found, try .aos folder*/
	free(config_path);
	config_path=stradd(HOME,"/.aos/config-",PACKAGE_VERSION,"/",name,NULL);
    }
    if(!exist(config_path)){
	free(config_path);
	config_path=NULL;
	warning2("Unable to determine the path to the configuration files.\n");
	warning2("Tried %s/config/%s\n", SRCDIR, name);
	warning2("Tried %s/config/%s\n", EXEP, name);
	warning2("Tried %s/.aos/config-%s/%s\n", HOME, PACKAGE_VERSION, name);
	warning2("Please set env MAOS_CONFIG_PATH=/path/to/config");
    }
    return config_path;
}

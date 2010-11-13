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

/**
   Collection all functions with respect to shared memory.
*/

#include <sys/mman.h>
#include <sys/stat.h>        /* For mode constants */
#include <sys/file.h>
#include <sys/statvfs.h>
#include <sys/types.h>
#include <sys/time.h>
#include <fcntl.h>     
#include <unistd.h>
#include <dirent.h>
#include <stdio.h>
#include <signal.h>
#include <errno.h>
#include <limits.h>
#include "misc.h"
#include "common.h"
#include "config.h"
#include "shm.h"
#if USE_POSIX_SHM == 1
int shm_keep_unused=0;
/**
   unmap a shared segment opened by fd.
*/

void shm_unmap(void *p, int fd){
    if(fd<=0) return;
    futimes(fd,NULL);//set access, modification time to current.
    struct stat buf;
    if(fstat(fd, &buf)){
	perror("fstat");
	error("fstat failed");
    }else{
	munmap(p, buf.st_size);
	flock(fd, LOCK_UN);
	close(fd);
	info2("Unmapping shared memory of size %ld.\n", (long)buf.st_size);
    }
}
void shm_free_older(long sec){
    /**
       Unlink shared segments that are mmaped more than sec ago.
    */
    remove_file_older("/dev/shm",sec);  
}

/**
   Free unused shared segments. All used segments are locked via shared
   lock. All shared segments that are not locked, and not what we are going to
   use, are considered not used and unlinked if the modification time is more
   than timeout seconds old.*/

int shm_free_unused(char *fnshm, int timeout){
    DIR *dir=opendir("/dev/shm");
    if(!dir){
	error("Unable to open directory /dev/shm\n");
    }
    struct dirent *dp;
    char fnshm2[NAME_MAX];
    int already_exist=0;
    struct stat buf;
    long sec2=myclocki()-timeout;
    while((dp=readdir(dir))){
	snprintf(fnshm2,NAME_MAX,"/%s", dp->d_name);
	if(!mystrcmp(fnshm2, "/maos_atm_")){//this is a maos atm
	    if(!fnshm || strcmp(fnshm, fnshm2)){
		//not equal to fnshm. 
		if(!shm_keep_unused){
		    int fd=shm_open(fnshm2, O_RDONLY, 00777);
		    if(fd<0){
			perror("open");
			warning("unable to open shared segment %s\n", fnshm2);
		    }
		    if(!flock(fd, LOCK_EX|LOCK_NB) ){//nobody is using this shm.
			fstat(fd, &buf);
			if(buf.st_mtime<sec2){
			    info2("\nUnlink %s ...", fnshm2);
			    shm_unlink(fnshm2);
			}
		    }else{
			//info2("\nKeep %s ...", fnshm2);
		    }
		    close(fd);
		}
	    }else{
		already_exist=1;
	    }
	}
    }
    closedir(dir);//don't forget to close it.
    return already_exist;
}

long shm_get_avail(void){
    /**
       returns available shm in bytes.
     */

    struct statvfs buf;
    statvfs("/dev/shm", &buf);
    return (long)buf.f_bsize * (long)buf.f_bavail;
}

#endif

/*
  Copyright 2009-2016 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
  Linux specific routines.
 */
#if defined(__linux__)




#include <limits.h>
#include <sys/types.h>
#include <dirent.h>


#include "common.h"
#include "process.h"
#include "misc.h"

/**
   Get the executable name of current process with full path.
*/
int get_job_progname(char *res, int nres, int pid){
    char path[PATH_MAX];
    /*readlink doesn't append \0 */
    char procpath[PATH_MAX];
    snprintf(procpath, PATH_MAX, "/proc/%d/exe", pid>0?pid:(int)getpid());
    int nprog=readlink(procpath,path,PATH_MAX);
    if(nprog>0){
	path[nprog]='\0';
	char *delt=strstr(path, "(deleted)");
	if(delt){
	    delt[0]='\0';
	}
	strncpy(res, path, nres); res[nres-1]=0;
    }else{
	warning("Failed to readlink %s\n", procpath);
    }
    return nprog>0?0:1;
}
/**
   Get the memory usage of current process.
 */
int get_job_mem(void){
    int mem;
    FILE* pfile;
    static int pagesize=0;//in kB
    if(!pagesize){
	pagesize=getpagesize()/1024;
    }
    if ((pfile=fopen("/proc/self/statm","r"))){
	if(fscanf(pfile, "%*d %d", &mem)!=1){
	    warning("failed to read statm\n");
	}
	mem*=pagesize;
	fclose(pfile);
    } else {
	mem=0;
    }
    return mem;
}
/**
   Get the launch time of current process.
*/
double get_job_launchtime(int pid){
    double starttime;
    char fnjob[64];
    snprintf(fnjob,64,"/proc/%d/stat",pid);
    FILE *fp=fopen(fnjob,"r");
    if(!fp){
	perror("fopen");
	warning("Unable to open file %s\n",fnjob);
	return 0;
    }
    unsigned long long starttime0;
    int nread;
    nread=fscanf(fp,"%*d %*s %*c %*d %*d %*d %*d %*d "/*tpgid */
		 "%*u %*u %*u %*u %*u %*u %*u "/*stime */
		 "%*d %*d %*d %*d %*d %*d %llu "/*starttime */
		 , &starttime0);
    fclose(fp);
    if(nread<0){
	warning("Unable to read %s, nread=%d\n", fnjob,nread);
	return 0;
    }
    starttime=(double)starttime0/(double)sysconf(_SC_CLK_TCK);
    return starttime;
}
static int proc_read_status(char *fnjob,  char *exename, char *cstat,long* nthread){
    FILE *fp=fopen(fnjob,"r");
    if(!fp){
	/*warning("Unable to open file %s\n",fnjob);//job exited */
	return -1;
    }
    int nread=fscanf(fp,"%*d %s %c %*d %*d %*d %*d %*d "/*tpgid */
		     "%*u %*u %*u %*u %*u %*u %*u "/*stime */
		     "%*d %*d %*d %*d %ld" , exename, cstat, nthread);
    fclose(fp);
    if(nread<0){
	/*warning("Unable to read %s, nread=%d\n", fnjob,nread); */
	return -1;
    }else{
	return 0;
    }
}
/**
   Return the number of running jobs.
*/
int get_usage_running(void){
  
    struct dirent *dp, *dpsub;
    char fnjob[64],fnsub[64];
    char cstat;
    char exename[256];
    int nrunning=0;
    long nthread=0;
    
    DIR *dir=opendir("/proc");
    if(!dir){
	warning("Unable to open /proc folder\n");
	return NCPU;
    }
    char mypid[32];
    snprintf(mypid,32,"%d",(int)getpid());
    while((dp=readdir(dir))){
	if(dp->d_name[0]>'9' || dp->d_name[0]<'0' || !strcmp(dp->d_name, mypid)){
	    continue;
	}
	snprintf(fnjob,64,"/proc/%s/stat",dp->d_name);
	if(proc_read_status(fnjob,  exename, &cstat,&nthread)){
	    continue;
	}
	if(strcmp(exename,"(maos)") && strcmp(exename,"(skyc)") && strcmp(exename,"(MATLAB)")){
	    continue;/*we only count maos or MATLAB. */
	}
	if(cstat!='R' && nthread==1){
	    continue;/*Job is not running; */
	}
	if(nthread>1 && strcmp(exename,"(MATLAB)")){
	    /*Only count 1 thread of MATLAB*/
	    nthread=1;
	}else if(nthread>1){/*There are many threads. Check each one. */
	    snprintf(fnsub,64,"/proc/%s/task",dp->d_name);
	    DIR* dirsub=opendir(fnsub);
	    nthread=0;
		
	    if(!dirsub){
		warning("Unable to read subdir %s\n", fnsub);
	    }else{
		char stat2;
		long nthread2;
		while((dpsub=readdir(dirsub))){
		    if(dpsub->d_name[0]<='9' && dpsub->d_name[0]>='0'){
			snprintf(fnsub,64,"/proc/%s/task/%s/stat",dp->d_name,dpsub->d_name);
			if(proc_read_status(fnsub, exename, &stat2, &nthread2)){
			    continue;
			}
			if(stat2=='R'){
			    nthread++;
			}
		    }
		}
		if(nthread>0){
		    cstat='R';
		}else{/*check the process again to avoid racing condition */
		    if(proc_read_status(fnjob,  exename, &cstat,&nthread2)){
			nthread=0;/*failed. */
		    }else{
			if(cstat=='R'){
			    nthread=1;
			    warning("Racing condition found\n");
			}else{
			    nthread=0;
			}
		    }
		}
		closedir(dirsub);
	    }
	}
	nrunning+=nthread;
    }
    closedir(dir);
    return nrunning;
}
/**
   Get the system load.
 */
double get_usage_load(void){
    double load=0;
    FILE *fp;
    fp=fopen("/proc/loadavg","r");
    if(!fp || fscanf(fp,"%lf",&load)!=1){
	warning("failed to read loadavg\n");
	load=0;
    }
    fclose(fp);
    return load;
}
/**
   Get the system memory usage.
*/
double get_usage_mem(void){
    double mem=0;
    FILE *fp;
    long memtot;
    long memfree;
    long membuf;
    long memcache;
    fp=fopen("/proc/meminfo","r");
    char line[1024];
    if(!fp) return 0;
    int count=0;
    while(fgets(line, sizeof(line), fp) && count<4){
	if(!mystrcmp(line, "MemTotal:")){
	    if(sscanf(line, "%*s %ld %*s", &memtot)!=1) break;
	    count++;
	}else if(!mystrcmp(line, "MemFree:")){
	    if(sscanf(line, "%*s %ld %*s", &memfree)!=1) break;
	    count++;
	}else if(!mystrcmp(line, "Buffers:")){
	    if(sscanf(line, "%*s %ld %*s", &membuf)!=1) break;
	    count++;
	}else if(!mystrcmp(line, "Cached:")){
	    if(sscanf(line, "%*s %ld %*s", &memcache)!=1) break;
	    count++;
	}
    }
    fclose(fp);
    if(count!=4){
	warning("failed to read meminfo\n");
	return 0;
    }else{
	mem=(double)(memtot-(memfree+membuf+memcache))/(double)memtot;
	return mem;
    }
}

/**
   Read the cpu counter.*/
int read_cpu_counter(long *user, long *tot){
    FILE *fp;
    long t_usr,t_nice, t_sys, t_idle;
    fp=fopen("/proc/stat","r");
    if(!fp || fscanf(fp, "%*s %ld %ld %ld %ld %*d %*d %*d %*d",
	      &t_usr,&t_nice,&t_sys,&t_idle)<0){
	fclose(fp);
	return -1;
    }
    fclose(fp);
    *user=t_usr+t_nice+t_sys;
    *tot=*user+t_idle;
    return 0;
}
/**
   return CPU usage of current process */
double read_self_cpu(void){
 
    static double t_last=0,s_last=0;
    long stime,utime;
    long t_usr,t_nice, t_sys, t_idle;
    FILE *fp=fopen("/proc/self/stat","r");
    if(!fp || fscanf(fp,"%*d %*s %*c %*d %*d %*d %*d %*d %*u %*u %*u %*u %*u %ld %ld",
	      &stime, &utime)!=2){
	fclose(fp);
	return 0;
    }
    fclose(fp);
    fp=fopen("/proc/stat","r");
    if(!fp || fscanf(fp, "%*s %ld %ld %ld %ld %*d %*d %*d %*d",
	      &t_usr,&t_nice,&t_sys,&t_idle)<0){
	fclose(fp); return 0;
    }
    fclose(fp);
    double t_tot=(double)(t_usr+t_nice+t_sys+t_idle);
    double s_tot=(double)(stime+utime);
    /*info("t_tot/TCK=%g\n",(t_tot-t_last)/TCK); == NCPU */
    double frac=(s_tot-s_last)/(t_tot-t_last)*NCPU;
    t_last=t_tot;
    s_last=s_tot;
    return frac;
}
/**
   parsing /proc/cpuinfo to find real number of physical cores, while
   hyper-threading may be cheating us.*/
int get_ncpu(void){
    int ncpu0 = sysconf( _SC_NPROCESSORS_ONLN );
    FILE *fp=fopen("/proc/cpuinfo","r");
#define nmax 4096
    char line[1024];
    int phyid[nmax];/*record the list of physical cpu ids */
    int coreid[nmax];/*record the list of core ids */
    int iphy=0, icore=0;
    const char *s_phy="physical id"; /*records number of CPUs */
    const char *s_core="core id";    /*should record number of cores per cpu */
    const char *s_cores="cpu cores"; /*should record number of cores per cpu. */
    int ncore=0;
    while(fp && fgets(line,1024,fp)){
	if(!mystrcmp(line,s_phy)){/*contains physical id */
	    int kphy;
	    char *ss=index(line,':');
	    int jphy=strtol(ss+1, NULL, 10);
	    for(kphy=0; kphy<iphy; kphy++){
		if(phyid[kphy]==jphy){
		    break;
		}
	    }
	    if(kphy == iphy){
		phyid[iphy]=jphy;
		iphy++;
		if(iphy>nmax){
		    error("Over flow. Please make nmax bigger\n");
		}
	    }
	}
	if(!mystrcmp(line,s_core)){/*contains core id */
	    int kcore;
	    char *ss=index(line,':');
	    int jcore=strtol(ss+1, NULL, 10);
	    for(kcore=0; kcore<icore; kcore++){
		if(coreid[kcore]==jcore){
		    break;
		}
	    }
	    if(kcore == icore){
		coreid[icore]=jcore;
		icore++;
		if(icore>nmax){
		    error("Over flow\n");
		}
	    }
	}
	if(!mystrcmp(line,s_cores)){/*contains cpu cores */
	    int mcore=strtol(index(line,':')+1,NULL,10);
	    if(ncore==0 || ncore == mcore){
		ncore=mcore;
	    }else{
		warning("'cpu_core' Does not have a uniq number\n");
	    }
	}
    }
    int ncpu1=iphy*(ncore>icore?icore:ncore);
    fclose(fp);
    int ncpu=ncpu0<ncpu1?ncpu0:ncpu1;
    if(ncpu==0){
	warning("Unable to obtain CPU count\n");
	ncpu=4;
    }
    return ncpu;
#undef nmax
}

#endif

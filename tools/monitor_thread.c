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
  All routines in this thread runs in a separate threads from the main
  thread. In order to call any gdk/gtk routines, it calls to go through
  gdk_threads_add_idle
*/

#include <unistd.h>
#include "../sys/sys.h"
#include "monitor.h"
PROC_T **pproc;
int *nproc;
extern int* hsock;
static double *htime;//last time having signal from host.
int nhostup=0;
PNEW(mhost);
int sock_main[2]={0,0}; /*Use to talk to the thread that blocks in select()*/
static fd_set active_fd_set;

extern double *usage_cpu, *usage_cpu2;
extern double *usage_mem, *usage_mem2;

static PROC_T *proc_get(int id,int pid){
    PROC_T *iproc;
    if(id<0 || id>=nhost){
	error("id=%d is invalid\n", id);
    }
    LOCK(mhost);
    for(iproc=pproc[id]; iproc; iproc=iproc->next){
	if(iproc->pid==pid){
	    break;
	}
    }
    UNLOCK(mhost);
    return iproc;
}

static PROC_T *proc_add(int id,int pid){
    PROC_T *iproc;
    if((iproc=proc_get(id,pid))) return iproc;
    iproc=calloc(1, sizeof(PROC_T));
    iproc->iseed_old=-1;
    iproc->pid=pid;
    iproc->hid=id;
    LOCK(mhost);
    iproc->next=pproc[id];
    pproc[id]=iproc;
    nproc[id]++;
    gdk_threads_add_idle((GSourceFunc)update_title, GINT_TO_POINTER(id));
    UNLOCK(mhost);
    return iproc;
}

static void proc_remove_all(int id){
    PROC_T *iproc,*jproc=NULL;
    LOCK(mhost);
    for(iproc=pproc[id]; iproc; iproc=jproc){
	jproc=iproc->next;
	gdk_threads_add_idle((GSourceFunc)remove_entry, iproc);//frees iproc
    }
    nproc[id]=0;
    UNLOCK(mhost);
    pproc[id]=NULL;
    gdk_threads_add_idle((GSourceFunc)update_title, GINT_TO_POINTER(id));
}

static void proc_remove(int id,int pid){
    PROC_T *iproc,*jproc=NULL;
    for(iproc=pproc[id]; iproc; jproc=iproc,iproc=iproc->next){
	if(iproc->pid==pid){
	    LOCK(mhost);
	    if(jproc){
		jproc->next=iproc->next;
	    }else{
		pproc[id]=iproc->next;
	    }
	    nproc[id]--;
	    UNLOCK(mhost);
	    gdk_threads_add_idle((GSourceFunc)remove_entry, iproc);
	    gdk_threads_add_idle((GSourceFunc)update_title, GINT_TO_POINTER(id));
	    break;
	}
    }
}


static int host_from_sock(int sock){
    if(sock<0) return -1;
    for(int ihost=0; ihost<nhost; ihost++){
	if(hsock[ihost]==sock){
	    return ihost;
	}
    }
    return -1;
}


void add_host_wrap(int ihost){
    int cmd[3]={CMD_ADDHOST, 0, 0};
    cmd[1]=ihost;
    stwriteintarr(sock_main[1], cmd, 3);
}

/* Record the host upon connection */
static void host_added(int ihost, int sock){
    htime[ihost]=myclockd();
    proc_remove_all(ihost);/*remove all entries. */
    LOCK(mhost);
    nhostup++;
    hsock[ihost]=sock;
    FD_SET(sock, &active_fd_set);
    UNLOCK(mhost);
    add_host_wrap(-1);//wakes up listen_host().
    info2("connected to %s\n", hosts[ihost]);
    gdk_threads_add_idle(host_up, GINT_TO_POINTER(ihost));
}

/*remove the host upon disconnection*/
static void host_removed(int sock){
    int ihost=host_from_sock(sock);
    if(ihost==-1) return;
    close(sock);
    LOCK(mhost);
    nhostup--;
    hsock[ihost]=-1;
    FD_CLR(sock, &active_fd_set);
    UNLOCK(mhost);
    add_host_wrap(-1);
    gdk_threads_add_idle(host_down, GINT_TO_POINTER(ihost));
    info2("disconnected from %s\n", hosts[ihost]);
}

static void add_host(gpointer data){
    int ihost=GPOINTER_TO_INT(data);
    int todo=0;
    LOCK(mhost);
    if(hsock[ihost]==-1){
	hsock[ihost]--;//make it -2 so no concurrent access.
	todo=1;
    }
    UNLOCK(mhost);
    if(todo){
	int sock=connect_port(hosts[ihost], PORT, 0, 0);
	if(sock>-1){
	    int cmd[2];
	    cmd[0]=CMD_MONITOR;
	    cmd[1]=scheduler_version;
	    if(stwriteintarr(sock, cmd, 2)){//write failed.
		warning("Rare event: Failed to write to scheduler at %s\n", hosts[ihost]);
		close(sock);
		LOCK(mhost);
		hsock[ihost]=-1;
		UNLOCK(mhost);
	    }else{
		host_added(ihost, sock);
	    }
	}else{
	    LOCK(mhost);
	    hsock[ihost]=-1;
	    UNLOCK(mhost);
	}
    }
}

static int respond(int sock){
    int ihost=host_from_sock(sock);
    htime[ihost]=myclockd();
    int cmd[3];
    if(streadintarr(sock, cmd, 3)){
	return -1;//failed
    }
    int pid=cmd[2];
    switch(cmd[0]){
    case CMD_VERSION:
	break;
    case CMD_STATUS:
	{
	    PROC_T *p=proc_get(ihost,pid);
	    if(!p){
		p=proc_add(ihost,pid);
	    }
	    if(stread(sock, &p->status, sizeof(STATUS_T))){
		return -1;
	    }
	    if(p->status.info==S_REMOVE){
		proc_remove(ihost, pid);
	    }else{
		if(p->status.info==S_WAIT && p->status.iseed<0){
		    p->pid=-p->status.iseed;
		}
		gdk_threads_add_idle((GSourceFunc)refresh, p);
	    }
	}
	break;
    case CMD_PATH:
	{
	    PROC_T *p=proc_get(ihost,pid);
	    if(!p){
		p=proc_add(ihost,pid);
	    }
	    if(streadstr(sock, &p->path)){
		return -1;
	    }
	    char *tmp=NULL;
	    while((tmp=strchr(p->path, '\n'))){
		tmp[0]=' ';
	    }
	}
	break;
    case CMD_LOAD:
	{
	    usage_cpu[ihost]=(double)((pid>>16) & 0xFFFF)/100.;
	    usage_mem[ihost]=(double)(pid & 0xFFFF)/100.;
	    usage_cpu[ihost]=MAX(MIN(1,usage_cpu[ihost]),0);
	    usage_mem[ihost]=MAX(MIN(1,usage_mem[ihost]),0);
	    gdk_threads_add_idle((GSourceFunc)update_progress, GINT_TO_POINTER(ihost));
	}
	break;
    case CMD_ADDHOST:
	if(cmd[1]>-1 && cmd[1]<nhost){
	    pthread_t tmp;
	    pthread_create(&tmp, NULL, (void*(*)(void*))add_host, GINT_TO_POINTER(cmd[1]));
	}else if(cmd[2]==-2){
	    return -2;
	}
	break;
    default:
	warning3("Invalid cmd %d\n",cmd[0]);
	return -1;
    }
    return 0;
}
/**
   listen_host() live in a separate thread, it has the following resposibilities:
   1) listening commands from the main thread to initiate connection to servers
   2) listening to connected servers for maos status event and update the display
   3) monitor connected servers for activity. Disable pages when server is disconnected.

   write to sock_main[1] will be caught by select in listen_host(). This wakes it up.*/
void listen_host(){
    htime=calloc(nhost, sizeof(double));
    FD_ZERO(&active_fd_set);
    FD_SET(sock_main[0], &active_fd_set);
    int keep_listen=1;
    while(keep_listen){
	fd_set read_fd_set = active_fd_set;
	if(select(FD_SETSIZE, &read_fd_set, NULL, NULL, NULL)<0){
	    perror("select");
	    continue;
	}
	for(int i=0; i<FD_SETSIZE; i++){
	    if(FD_ISSET(i, &read_fd_set)){
		int res;
		res=respond(i);
		if(res==-2){//quit
		    keep_listen=0;
		    break;
		}else if(res==-1){//remove host
		    host_removed(i);
		}
	    }
	}
	double ntime=myclockd();
	for(int ihost=0; ihost<nhost; ihost++){
	    if(hsock[ihost]>-1){
		if(htime[ihost]+10<ntime){
		    //10 seconds grace period
		    info2("10 seconds no respond. disconnect\n");
		    host_removed(hsock[ihost]);
		}
	    }
	}
    }
    for(int i=0; i<FD_SETSIZE; i++){
	if(FD_ISSET(i, &active_fd_set)){
	    shutdown(i, SHUT_RDWR);
	    usleep(100);
	    close(i);
	    FD_CLR(i, &active_fd_set);
	}
    }
}

/**
   called by monitor to talk to scheduler.
*/
int scheduler_cmd(int host,int pid, int command){
    int sock=hsock[host];
    if(sock==-1){
	add_host_wrap(host);
	sleep(1);
    }
    if(sock==-1) return 1;
    int cmd[2];
    cmd[0]=command;
    cmd[1]=pid;/*pid */
    int ans=stwriteintarr(sock,cmd,2);
    if(ans){/*communicated failed.*/
	close(sock);
    }
    return ans;
}

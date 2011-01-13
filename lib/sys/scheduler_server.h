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
#ifndef AOS_SCHEDULER_SERVER_H
#define AOS_SCHEDULER_SERVER_H
#define MAXMSG  512
/*
#include <sys/types.h>       
#include <sys/socket.h>
*/
#include <netinet/in.h>
/*
  Struct to hold status information of a running process.
 */
#if defined(__INTEL_COMPILER)
/*with htons defined in glibc 2.4, intel compiler complains
  about conversion from in to uint16. THis is an ugly workaround*/
#undef htons
#define htons myhtons
static inline uint16_t myhtons(uint16_t port){
    uint16_t ans;
#if __BYTE_ORDER == __BIG_ENDIAN
    ans=(port);
#else
    ans=(unsigned short int)
	((((port) >> 8) & 0xff) | (unsigned short int)(((port) & 0xff) << 8));
#endif
    return ans;
}
#endif
/**
   Status about each process. Timing, wavefront error, etc.
 */
typedef struct STATUS_T{
    //Individual timing
    double wfs;
    double recon;
    double cache;
    double eval;
    double tot;
    double scale;
    double mean;
    long rest;
    long laps;
    //status report
    int iseed;
    int nseed;
    int nthread;
    int isim;
    int simstart;
    int simend;
    int done;
    double clerrlo;
    double clerrhi;
    int info;
    time_t timstart;
    time_t timend;
}STATUS_T;
/**
   Struct to hold information of queued jobs waited to start.
*/
typedef struct QUEUE_T{
    int pid;
    int sock;
    int nthread;
}QUEUE_T;
/**
   Struct to hold information of running jobs.
*/
typedef struct RUN_T{
    struct RUN_T *next;
    STATUS_T status;
    double started;//started execution.
    double launchtime;
    int pid;
    int sock;
    int nthread;
    int time;
    char *path;
}RUN_T;

/**
  Struct to hold available monitors waiting for information.
*/
typedef struct MONITOR_T{
    int sock;
    int load;//handle machine load information.
    struct MONITOR_T *next;
}MONITOR_T;

int myhostid(const char *host);
int make_socket (uint16_t port, int retry);
MONITOR_T *monitor_get(int hostid);
void monitor_remove(int hostid);
MONITOR_T *monitor_add(int hostid);
void monitor_send(RUN_T *run,char*path);
void monitor_send_initial(MONITOR_T *ic);
void monitor_send_load(void);
void print_backtrace(int sig);
void print_backtrace_symbol(void *const *buffer, int size);

enum{
    CMD_START=1,
    CMD_FINISH,
    CMD_STATUS,
    CMD_CRASH,//4
    CMD_MONITOR,
    CMD_PATH,//6
    CMD_KILL,
    CMD_TRACE,//8
    CMD_DRAW,
    CMD_BROADCAST,
    CMD_REMOVE,
    CMD_HOSTUP,
    CMD_VERSION,//13
    CMD_LOAD,//14
    CMD_SHUTRD,//15
    CMD_SHUTWR,//16
};
enum{
    S_RUNNING=1,
    S_WAIT,
    S_START,
    S_FINISH=11,
    S_CRASH,
    S_TOKILL,
    S_REMOVE,
    S_KILLED,
};
extern uint16_t PORTMON;
extern uint16_t PORT;
//#define MAXNHOST 10
extern int nhost;
extern char** hosts;
extern int hid;
#define scheduler_version 0x19
#endif

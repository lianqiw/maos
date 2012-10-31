/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include <netinet/in.h>

/**
   Status about each process. Timing, wavefront error, etc.
 */
typedef struct STATUS_T{
    /*Individual timing */
    double wfs;
    double recon;
    double cache;
    double eval;
    double tot;
    double scale;
    double mean;
    long rest;
    long laps;
    /*status report */
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

enum{
    CMD_START=1,
    CMD_FINISH,
    CMD_STATUS,
    CMD_CRASH,/*4 */
    CMD_MONITOR,
    CMD_PATH,/*6 */
    CMD_KILL,
    CMD_TRACE,/*8 */
    CMD_DRAW,
    CMD_BROADCAST,
    CMD_REMOVE,
    CMD_HOSTUP,
    CMD_VERSION,/*13 */
    CMD_LOAD,/*14 */
    CMD_SHUTRD,/*15 */
    CMD_SHUTWR,/*16 */
    CMD_ADDHOST,/*17*/
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

#define scheduler_version 0x24
#define BACKTRACE_CMD_LEN 200


void print_backtrace(int sig);
void print_backtrace_symbol(void *const *buffer, int size);
void socket_tcp_keepalive(int sock);
void scheduler(void);


#endif

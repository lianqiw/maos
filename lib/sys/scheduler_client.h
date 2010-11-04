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
#ifndef AOS_SCHEDULER_CLIENT_H
#define AOS_SCHEDULER_CLIENT_H
#define MAXMSG  512
#include "scheduler_server.h"

#if defined (__CYGWIN__)
#define start_scheduler(A...)
#define scheduler_start(A...)
#define scheduler_wait(A...) 0
#define scheduler_finish(A...)
#define scheduler_report(A...)
#define scheduler_get_drawdaemon(A...) NULL
#define scheduler()
#define print_backtrace(A...)
#else
/*
#include <sys/types.h>   
#include <sys/socket.h>
#include <netinet/in.h>*/
void start_scheduler(int argc, char **argv);
int init_sockaddr (struct sockaddr_in *name, const char *hostname, uint16_t port);
int scheduler_connect(int ihost, int block,int mode);
int scheduler_connect_self(int block,int mode);
//called by maos
//call scheduler_start, and then scheduler_wait if waiting is 1.
int scheduler_start(char *path, int nthread, int waiting);
int scheduler_wait(void);
void scheduler_finish(int status);
void scheduler_report(STATUS_T *status);
//called by monitor
int scheduler_kill_job(int ihost,int pid);
int scheduler_remove_job(int ihost, int pid);

char* scheduler_get_drawdaemon(int pid);
void  scheduler_shutdown(int *sock, int mode);
#endif//cygwin
#endif

/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "scheduler.h"
#include "thread.h"
/**
   \file scheduler_client.h
   Contains routines that will be used to talk to the scheduler.
*/
#include <stdint.h>
extern int nhost;
extern char **hosts;
extern char **hostshort;
/*hosts*/
void parse_host(const char *line);
void free_hosts();
void init_hosts();
const char *lookup_hostaddr(const char *hostname);
const char *lookup_hostname(const char *hostaddr);
/*called by maos */
int scheduler_connect(const char* hostname);
pthread_t scheduler_listen(thread_fun fun);
void scheduler_report_path(const char* path);
void scheduler_start(int nthread, int ngpu, int waiting);
void scheduler_finish(int status);
void scheduler_report(status_t *status);
/*called by monitor */
int scheduler_launch_exe(char *hostarr, int argc, const char *argv[]);
/*save a socket for draw()*/
int scheduler_socket(int dir, int *sfd, int id);
/*Handling backtrace*/
int call_addr2line(char *ans, int nans, const char *cmd);
int print_backtrace_symbol(void *const *buffer, int size);

void print_backtrace();

#endif

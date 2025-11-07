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
#ifndef AOS_SCHEDULER_H
#define AOS_SCHEDULER_H
#ifdef HAVE_CONFIG_H
#include "config.h" 
#endif
#include <poll.h>
typedef struct ws_proxy_t{
	int fd_remote;//socket passed to maos
	int fd;//socket to read from maos
	int ismonitor;//this is a direct monitor server
	int (*send)(char* buf, int nlen, int mode, void *userdata);
	int (*forward)(int fd, int nlen, int mode, void *userdata);
	void *userdata;
}ws_proxy_t;
void runned_remove(int pid);
void running_kill(int pid);
int  maos_command(int pid, int sock, int cmd);
int send_draw_sock(int sock, int pid);
void monitor_add(int sock, int flag, int (*func)(char*buf, int nlen, int mode, void *userdata), void* userdata);
void monitor_remove(int sock);
void ws_proxy_command(char* in, size_t len, ws_proxy_t ws);
void ws_proxy_remove(void *userdata, int toclose);
#if HAS_LWS
int  start_lws(short port);
#endif
extern const int http_padding;
int http_handshake(struct pollfd *pfd, int flag);
#endif

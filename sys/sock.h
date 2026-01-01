/*
  Copyright 2009-2026 Lianqi Wang
  
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
#ifndef AOS_SOCKET_H
#define AOS_SOCKET_H
/**
   \file sock.h
   Routines to establish socket connection
*/
#include <netinet/in.h>
#include <poll.h>
int socket_recv_timeout(int sock, long sec);
int socket_send_timeout(int sock, long sec);
int connect_port(const char *hostname, int port, int block, int nodelay);
typedef struct listen_opt_t{
	int port;
	int nodelay;
	const char *localpath;//for local unix socket if set
	const char *ipaddr;//ip address to listen for port if set
	int (*responder)(struct pollfd*, int);
	int (*http_responder)(struct pollfd*, int);
	void (*timeout_fun)();
	double timeout_sec;
}listen_opt_t;
void listen_port(listen_opt_t opt);//uint16_t port, char *localpath, int (*responder)(struct pollfd*, int), double timeout_sec, void (*timeout_fun)(), int nodelay);
void listen_port_add(int sock, short events, int(*handler)(struct pollfd*, int), const char *caller);
void listen_port_del(int sock, int toclose, const char *caller);
int socket_nopipe(int sock);
int socket_block(int sock, int block);
int bind_socket(int protocol, const char* ipaddr, uint16_t port);
int socket_port(int sock);
in_addr_t socket_peer(int sock);
const char *addr2name(in_addr_t s_addr);
#endif

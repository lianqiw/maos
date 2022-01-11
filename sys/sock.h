/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
int socket_recv_timeout(int sock, double sec);
int socket_send_timeout(int sock, double sec);
int connect_port(const char *hostname, int port, int block, int nodelay);
void listen_port(uint16_t port, char *localpath, int (*respond)(int), double timeout_sec, void (*timeout_fun)(), int nodelay);
int socket_nopipe(int sock);
int socket_block(int sock, int block);
int bind_socket(int protocol, char* ip, uint16_t port);
int socket_port(int sock);
in_addr_t socket_peer(int sock);
const char *addr2name(in_addr_t s_addr);
#endif

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
#ifndef AOS_SCHEDULER_WS_H
#define AOS_SCHEDULER_WS_H
typedef struct l_message{
	char* payload;
	size_t len;
	struct l_message* next;
}l_message;
int  ws_start(short port);
void ws_end();
void ws_service();
void ws_append(const char* in, size_t len);
int  ws_proxy_get_fd(void* userdata);
void ws_proxy_add(int sock, void* userdata);
void ws_proxy_remove(void* userdata, int toclose);
int  ws_proxy_read(struct pollfd *pfd, int flag);
void scheduler_push_ws(l_message** head, long prepad);
void scheduler_receive_ws(char* in, size_t len, void* userdata);
#endif

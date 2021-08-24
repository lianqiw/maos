/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

int ws_service(int waiting);
int ws_start(short port);
struct l_message{
	char* payload;
	size_t len;
	struct l_message* next;
};
typedef struct l_message l_message;
void ws_push(const char* in, size_t len);
void html_push_all(l_message** head, l_message** tail, long prepad, long postpad);
void scheduler_handle_ws(char* in, size_t len);
#endif

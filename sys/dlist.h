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
#ifndef AOS_DLIST_H
#define AOS_DLIST_H
typedef struct dlist{
	struct dlist *next;
	struct dlist *prev;
}dlist;
void add_head(dlist **phead, dlist **ptail, dlist *pdata);
dlist *remove_head(dlist **phead, dlist **ptail);
void add_tail(dlist **phead, dlist **ptail, dlist *pdata);
dlist *remove_tail(dlist **phead, dlist **ptail);
#endif

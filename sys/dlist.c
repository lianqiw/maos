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
/**
	\file dhlist.c Implements a double linked list. Not thread safe.
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "dlist.h"

void add_head(dlist **phead, dlist **ptail, dlist *pdata){
	if(!phead||!ptail||!pdata) return;
	pdata->next=*phead;
	pdata->prev=NULL;
	if(*phead){
		(*phead)->prev=pdata;
	}
	*phead=pdata;
	if(!*ptail){
		*ptail=*phead;
	}
}
dlist *remove_head(dlist **phead, dlist **ptail){
	if(!phead||!ptail||!*phead) return NULL;
	dlist *pdata=*phead;
	*phead=pdata->next;
	if(*phead){
		(*phead)->prev=NULL;
	} else{
		*ptail=NULL;
	}
	pdata->next=NULL;
	return pdata;
}
void add_tail(dlist **phead, dlist **ptail, dlist *pdata){
	if(!phead||!ptail||!pdata) return;
	if(!*ptail){
		add_head(phead, ptail, pdata);
	} else{
		(*ptail)->next=pdata;
		pdata->prev=*ptail;
		pdata->next=NULL;
		*ptail=pdata;
	}
}
dlist *remove_tail(dlist **phead, dlist **ptail){
	if(!phead||!ptail||!*ptail) return NULL;
	dlist *pdata=*ptail;
	*ptail=pdata->prev;
	pdata->prev=NULL;
	if(*ptail){
		(*ptail)->next=NULL;
	} else{
		*phead=NULL;
	}
	return pdata;
}

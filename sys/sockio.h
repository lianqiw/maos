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
#ifndef AOS_SOCKIO_H
#define AOS_SOCKIO_H
int stwrite(int sfd, const void *p, size_t len);
int stread(int sfd, void *p, size_t len);
int stwriteint(int sfd, int cmd);
int stwriteintarr(int sfd, int* cmd, unsigned int len);
int streadint(int sfd, int *res);
int streadintarr(int sfd, int* cmd, unsigned int len);
int stwritestr(int sfd, const char *str);
int streadstr(int sfd, char **str);
#endif

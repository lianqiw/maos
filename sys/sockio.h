/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#ifdef HAVE_CONFIG_H
#include "config.h" 
#endif
int stwrite(int sfd, const void *p, size_t len);
int stread(int sfd, void *p, size_t len);
INLINE int stwriteint(int sfd, int cmd){
    return stwrite(sfd, &cmd, sizeof(int));
}
INLINE int stwriteintarr(int sfd, int* cmd, unsigned int len){
    return stwrite(sfd,cmd,len*sizeof(int));
}
INLINE int streadint(int sfd, int *cmd){
    return stread(sfd, cmd, sizeof(int));
}
INLINE int streadintarr(int sfd, int* cmd, unsigned int len){
    return stread(sfd,cmd,len*sizeof(int));
}
int stwritestr(int sfd, const char *str);
int streadstr(int sfd, char **str);
int stwritestrarr(int sfd, const char *const *str, int nstr);
int streadstrarr(int sfd, char ***str, int *nstr);
int stwritefd(int sfd, int fd);
int streadfd(int sfd, int *fd);
#endif

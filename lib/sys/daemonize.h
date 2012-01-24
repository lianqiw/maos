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
/*
  System specific routines to interact with the system to daemonize a program.
*/
#ifndef AOS_PROCESS_H
#define AOS_PROCESS_H
void single_instance_daemonize(const char *lockfolder_in, 
			       const char *progname,long version,
			       void(*daemon_func)(void*), 
			       void* daemon_arg);
int lock_file(const char *fn, long block, long version);
void daemonize(void);
extern int detached;
#endif


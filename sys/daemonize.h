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

/*
  System specific routines to interact with the system to daemonize a program.
*/
#ifndef AOS_PROCESS_H
#define AOS_PROCESS_H
/**
   \file daemonize.h
   Handling process creation.
*/
int single_instance_daemonize(const char *lockfolder_in, 
			       const char *progname, int version,
			       void(*daemon_func)(void*), 
			       void* daemon_arg);
int lock_file(const char *fn, int block);
int lock_file_version(const char *fn, int block, int version);
void daemonize(void);
void redirect(void);
pid_t launch_exe(const char *exepath, const char *cmd);
char* find_exe(const char *name);
int spawn_process(const char *exename, const char *arg, const char *path);
extern int detached;
#endif

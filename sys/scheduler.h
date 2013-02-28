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
#ifndef AOS_SCHEDULER_SERVER_H
#define AOS_SCHEDULER_SERVER_H
/**
   Status about each process. Timing, wavefront error, etc.
 */
typedef struct STATUS_T{
    /*Individual timing */
    double wfs;
    double recon;
    double cache;
    double eval;
    double tot;
    double scale;
    double mean;
    long rest;
    long laps;
    /*status report */
    int iseed;
    int nseed;
    int nthread;
    int isim;
    int simstart;
    int simend;
    int done;
    double clerrlo;
    double clerrhi;
    int info;
    time_t timstart;
    time_t timend;
}STATUS_T;

enum{
    CMD_START=1,
    CMD_FINISH,
    CMD_STATUS,
    CMD_CRASH,/*4 */
    CMD_MONITOR,
    CMD_PATH,/*6 */
    CMD_KILL,
    CMD_TRACE,/*8 */
    CMD_DRAW,
    CMD_SOCK,  /**<10 We are pass a socket*/
    CMD_REMOVE,/**<11 Remove a job*/
    CMD_DISPLAY,/**<12 Remote display for telemetry*/
    CMD_VERSION,/*13 */
    CMD_LOAD,/*14 */
    CMD_UNUSED2,/*15*/
    CMD_UNUSED3,/*16*/
    CMD_ADDHOST,/*17*/
    CMD_LAUNCH,/*18*/
    CMD_PAUSE,
    CMD_RESUME,
};
enum{
    S_RUNNING=1,/*1*/
    S_WAIT,  /*2*/
    S_START,/*3*/
    S_QUEUED,/*4*/
    S_FINISH=11,/*11*/
    S_CRASH,/*12*/
    S_TOKILL,/*13*/
    S_REMOVE,/*14*/
    S_KILLED,/*15*/
};

#define scheduler_version 0x26
#endif

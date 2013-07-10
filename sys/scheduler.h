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
	char done;
	char warning;
	char unused1;
	char unused2;
    double clerrlo;
    double clerrhi;
    int info;
    time_t timstart;
    time_t timend;
}STATUS_T;
/*For scheduler()*/
enum{
    CMD_START=1,
    CMD_FINISH,
    CMD_STATUS,
    CMD_CRASH,/*4 */
    CMD_MONITOR,
    CMD_PATH,/*6 */
    CMD_KILL,
    CMD_TRACE,/*8 */
    CMD_UNUSED0,
    CMD_SOCK,  /**<10 We are pass a socket*/
    CMD_REMOVE,/**<11 Remove a job*/
    CMD_DISPLAY,/**<12 Remote display for telemetry*/
    CMD_UNUSED1,/*13 */
    CMD_UNUSED2,/*14 */
    CMD_RESTART,/*15*/
    CMD_UNUSED3,/*16*/
    CMD_UNUSED4,/*17*/
    CMD_LAUNCH,/*18*/
    CMD_PAUSE,
    CMD_RESUME,
};
/*command from scheduler etc to maos*/
enum{
    MAOS_SERVER=1,/*tell maos to act as server*/
    MAOS_DRAW=10,/*tell maos to start drawing*/
    MAOS_ASSIGN_WFS=100,
    MAOS_ASSIGN_EVL,/*101*/
    MAOS_ASSIGN_RECON,/*102*/
    MAOS_ASSIGN_DONE=199,/*199*/
};
/*command from scheduler etc and monitor*/
enum{
    MON_STATUS=3,
    MON_PATH=6,
    MON_VERSION=13,
    MON_LOAD=14,
    MON_ADDHOST=17,
};
/*For job status*/
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
    S_NONEXIST,/*16*/
};
#define scheduler_version 0x27
#endif

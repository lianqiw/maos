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
#ifndef AOS_SCHEDULER_SERVER_H
#define AOS_SCHEDULER_SERVER_H
/**
   Status about each process. Timing, wavefront error, etc.
 */
typedef struct status_t{
    /*Individual timing */
    double wfs;
    double recon;
    double other;
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
    time_t timstart;//Job launch time
    time_t timlast; //last status change
}status_t;
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
    CMD_PROBE,/*9 just for probe connection*/
    CMD_SOCK,  /**<10 We are pass a socket*/
    CMD_REMOVE,/**<11 Remove a job*/
    CMD_DISPLAY,/**<12 Remote display for telemetry*/
    CMD_MAOS, /*13 */
    CMD_MAOSDAEMON, /*14 */
    CMD_RESTART,/*15*/
    CMD_KILLED,/*16*/
    CMD_UNUSED4,/*17*/
    CMD_LAUNCH,/*18*/
    CMD_PAUSE,
    CMD_RESUME,
};
/*command from scheduler etc to maos*/
enum{
    MAOS_SERVER=1,/*not used*/
    MAOS_DRAW=10, /*tell maos to start drawing with the received fd*/
    MAOS_VAR=20,  /*tell maos to enable variable access*/
    MAOS_PAUSE,   /*tell maos to pause or resume*/
    MAOS_ASSIGN_WFS=100,
    MAOS_ASSIGN_EVL,/*101*/
    MAOS_ASSIGN_RECON,/*102*/
    MAOS_ASSIGN_DONE=199,/*199*/
};
/*command from scheduler etc and monitor*/
enum{
    MON_CMD=1, /*called by monitor main thread to relay cmd to scheduler*/
    MON_DRAWDAEMON=2,/*start drawdaemon*/
    MON_STATUS=3,
    MON_PATH=6,
    
    MON_CLEARJOB=7,
    MON_KILLJOB=8,
    MON_SAVEJOB=9,

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
    S_UNKNOWN,//5
    S_FINISH=11,/*11*/
    S_CRASH,/*12*/
    S_TOKILL,/*13*/
    S_REMOVE,/*14*/
    S_KILLED,/*15*/
    S_NONEXIST,/*16*/
};
//For tools/drawdaemon and lib/draw.c
enum{
    DRAW_START=0, /*Mark the starting of data stream. */
    DRAW_DATA,//1
    DRAW_SHM,//2
    DRAW_POINTS,//3
    DRAW_STYLE,//4
    DRAW_CIRCLE,//5
    DRAW_LIMIT,//6
    DRAW_FIG,//7
    DRAW_NAME,//8
    DRAW_TITLE,//9
    DRAW_XLABEL,//10
    DRAW_YLABEL,//11
    DRAW_ZLIM,//12
    DRAW_LEGEND,//13
    DRAW_XYLOG,//14
    DRAW_FIGFN,//15
    DRAW_PAUSE,//16
    DRAW_RESUME,//17
    DRAW_FINAL,//18 this client is done.
    DRAW_FLOAT,//19 number of bytes for float
    DRAW_FRAME,//20 information about UDP frame, following this is frame index, total bytes, sub-frame index, sub-frame bytes/
    DRAW_SINGLE,//21 toggle single/continuous drawing
    DRAW_UDPPORT,//22 send udp port of the client
    DRAW_END=100,
    DRAW_ENTRY=9999 /*A new entry*/
};
#define scheduler_version 51
#endif

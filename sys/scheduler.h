/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
/*Commands handled by scheduler server. Do not modify existing values*/
enum{
    CMD_START=1,/*1: called by maos to check whether it can start*/
    CMD_FINISH, /*2: called by maos to report finishes*/
    CMD_STATUS, /*3: called by maos to report status*/
    CMD_CRASH,	/*4: called by maos te report crashes */
    CMD_MONITOR,/*5: called by monitor upon connection*/
    CMD_PATH,	/*6: called by maos to report job path */
    CMD_KILL, 	/*7: called by monitor to kill a job*/
    CMD_TRACE,	/*8: called by maos to do a backtrace */
    CMD_PROBE,	/*9: called by monitor to probe connection*/
    CMD_DRAWCLI,/*10: called by maos to save or request a drawdaemon socket*/
    CMD_REMOVE, /*11: called by monitor to remove a job */
    CMD_DRAWSER,/*12: called by monitor or drawdaemon to provide drawdaemon connection*/
	CMD_MAOSCLI,/*13: for a maos client to create a client link to maos. */
	CMD_MAOSSER,/*14: called by maos to save a port to run maos_command */
	CMD_RESTART,/*15: called by monitor to restart a job*/
	CMD_KILLED, /*16: called by maos to indicate that job is cancelled or killed*/
    CMD_DUMMY,	/*17. no action. for probing connection with no reply anticipated*/
	CMD_LAUNCH,	/*18: called from maos from another machine to start a job in this machine*/
};
/*Command handled by maos built-in server. Do not modify existing values.*/
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
/*Commands handled by monitor. May modify by restarting scheduler*/
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
/*For job status. Do not modify existing values*/
enum{
    //below 10 is active
    S_RUNNING=1,/*1*/
    S_WAIT,  /*2*/
    S_START,/*3*/
    S_QUEUED,/*4*/
    S_UNKNOWN,//5
    //above 10 is no longer active
    S_FINISH=11,/*11*/
    S_CRASH,/*12*/
    S_TOKILL,/*13*/
    S_REMOVE,/*14*/
    S_KILLED,/*15*/
    S_NONEXIST,/*16*/
};
//Commands exchanged between drawdaemon and draw() from maos or drawres. Do not modify existing values
enum{
    DRAW_START=0, /*Mark the starting of data stream. */
    DRAW_DATA,//1
    DRAW_HEARTBEAT,//2
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
    DRAW_INIT,//23 set the time to delete old plots that are not updated when draw_final is called
    DRAW_PID,//24 followed by pid
    DRAW_ZLOG,//25 log scale for image.
	DRAW_PATH,//26 directory used for saving files
	DRAW_EXENAME,//27 executable that calls draw().
    DRAW_END=100,//data send is over. start drawing
    
    DRAW_ENTRY=9999 /*A new entry*/
};
enum {
	DRAW_ID_MAOS=1,
	DRAW_ID_RES,
	DRAW_ID_BIN,
	DRAW_ID_TOT
};
#define scheduler_version 53
#endif

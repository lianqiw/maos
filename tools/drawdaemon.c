/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#include "../lib/draw.h"
#include "process.h"
#include "misc.h"
#include "mathmisc.h"
#include "drawdaemon.h"

char *fifo=NULL;
int fifopid=0;
int main(int argc, char *argv[])
{
    if(!g_thread_supported()){
	g_thread_init(NULL);
	gdk_threads_init();
    }
    gtk_init(&argc, &argv);
    if(argc>1){
	fifo=argv[1];
	char *tmp=strstr(fifo,"drawdaemon_");
	if(tmp){
	    fifopid=strtol(tmp+11,NULL,10);
	}
    }else{
	fifo=calloc(80,sizeof(char));
	fifopid=(int)getppid();
	snprintf(fifo,80,"%s/drawdaemon_%d.fifo",TEMP,fifopid);
    }
    int ppid;
    const char *fifo2=strstr(fifo,"drawdaemon");
    if(!fifo2){
	warning("drawdaemon not found in string\n");
	ppid=getpid();
    }else{
	sscanf(fifo2, "drawdaemon_%d.fifo", &ppid);
    }
    {//record the drawdaemon pid. and redirect output
	char fnpid[PATH_MAX];
	snprintf(fnpid, PATH_MAX,"%s/drawdaemon_%d.pid", TEMP, ppid);
	FILE *fp=fopen(fnpid, "w");
	fprintf(fp, "%d", (int)getpid());
	fclose(fp);
    
	char fnlog[PATH_MAX];
	snprintf(fnlog, PATH_MAX,"%s/drawdaemon_%d.log", TEMP, ppid);
	if(!freopen(fnlog, "w", stdout)) {
	    perror("freopen");
	    warning("Error redirect stdout\n");
	}
	if(!freopen(fnlog, "w", stderr)) {
	    perror("freopen");
	    warning("Error redirect stderr\n");
	}
	setbuf(stdout,NULL);//disable buffering.
	setbuf(stderr,NULL);
    }
    create_window();
    g_thread_create((GThreadFunc)open_fifo, NULL, 0, NULL);
    gdk_threads_enter();
    gtk_main();
    gdk_threads_leave();
    {
	char fn[PATH_MAX];
	snprintf(fn, PATH_MAX,"%s/drawdaemon_%d.pid", TEMP, ppid);
	remove(fn);
	snprintf(fn, PATH_MAX,"%s/drawdaemon_%d.log", TEMP, ppid);
	remove(fn);
    }
    remove(fifo);
}//main

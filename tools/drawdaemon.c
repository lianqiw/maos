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

#include "../lib/aos.h"
#include "drawdaemon.h"
int sock;
int main(int argc, char *argv[]){
#if GLIB_MAJOR_VERSION<3 && GLIB_MINOR_VERSION<32
    if(!g_thread_supported()){
	g_thread_init(NULL);
	gdk_threads_init();
    }
#endif
    gtk_init(&argc, &argv);
    if(argc<2){
	error("Must call drawdaemon with at least one argument\n");
    }
    sock=strtol(argv[1], NULL, 10);
    if(sock<0){
	error("sock=%d\n", sock);
    }
    info("sock=%d\n", sock);
    {
	char fnlog[PATH_MAX];
	snprintf(fnlog, PATH_MAX,"%s/drawdaemon.log", TEMP);
	if(!freopen(fnlog, "w", stdout)) {
	    perror("freopen");
	    warning("Error redirect stdout\n");
	}
	if(!freopen(fnlog, "w", stderr)) {
	    perror("freopen");
	    warning("Error redirect stderr\n");
	}
	setbuf(stdout,NULL);
	setbuf(stderr,NULL);
    }
    g_thread_new("listen_draw", (GThreadFunc)listen_draw, NULL);
    create_window();
    gtk_main();
}/*main */

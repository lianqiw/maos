/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>

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

#include "drawdaemon.h"
#include "icon-draw.h"
int sock;
int sock_idle=0;
int cumu=0;
#if MAC_INTEGRATION
#include <gtkosxapplication.h>
#endif
GdkPixbuf* icon_main=NULL;
int main(int argc, char* argv[]){
	{
		char fnlog[PATH_MAX];
		snprintf(fnlog, PATH_MAX, "%s/drawdaemon.log", TEMP);
		if(!freopen(fnlog, "w", stdout)){
			perror("freopen");
			warning("Error redirect stdout\n");
		}
		if(!dup2(fileno(stdout), fileno(stderr))){
			perror("freopen");
			warning("Error redirect stderr\n");
		}
		setbuf(stdout, NULL);
	}
	info("drawdaemon is launched with %s %s\n", argv[0], argv[1]);
#if GLIB_MAJOR_VERSION<3 && GLIB_MINOR_VERSION<32
	if(!g_thread_supported()){
		g_thread_init(NULL);
		gdk_threads_init();
	}
#endif
#if GTK_MAJOR_VERSION<4
	gtk_init(&argc, &argv);
#else
	gtk_init();
#endif
#if MAC_INTEGRATION
	GtkosxApplication* theApp=g_object_new(GTKOSX_TYPE_APPLICATION, NULL);
#endif

	icon_main=gdk_pixbuf_new_from_inline(-1, icon_inline_draw, FALSE, NULL);
	if(argc<2){
		error("Must call drawdaemon with the socket fd.\n");
	}
	sock=strtol(argv[1], NULL, 10);
	if(sock<0){
		error("sock=%d\n", sock);
	}
	//dbg("sock=%d\n", sock);
	socket_block(sock, 0);

#if MAC_INTEGRATION
	gtkosx_application_set_dock_icon_pixbuf(theApp, icon_main);
	gtkosx_application_ready(theApp);
#endif
	thread_new(listen_draw, NULL);
	//g_thread_new("listen_draw", (GThreadFunc)listen_draw, NULL);
	create_window();
	gtk_main();
}/*main */

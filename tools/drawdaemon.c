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


#include "drawdaemon.h"
//int sock;
const char *host=0;
#if MAC_INTEGRATION
#include <gtkosxapplication.h>
#endif
GdkPixbuf* icon_main=NULL;
#if GTK_MAJOR_VERSION>=4
static void
activate (GtkApplication *app,
          gpointer        user_data)
{
	(void)user_data;
	GtkWidget *window = gtk_application_window_new (app);
	create_window(window);
}
#endif

#if MAC_INTEGRATION
void mac_terminate(GtkosxApplication *app, gpointer psock){
	(void)app;
	if(psock){
		int sock0=*(int*)psock;
		info("close %d socket in mac_terminate\n", sock0);
		if(sock0!=-1) close(sock0);
	}
	sleep(1);
	gtk_main_quit();
}
#endif
int main(int argc, char* argv[]){
	info("drawdaemon is launched with %s %s\n", argv[0], argv[1]);
	if(getenv("MAOS_DIRECT_LAUNCH")){
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
		register_signal_handler(NULL);
	}
	
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

	icon_main=gdk_pixbuf_new_from_resource("/maos/icon-draw.png", NULL);
	if(argc<2){
		error("Usage: %s socket or hostname.\n", argv[0]);
	}

#if MAC_INTEGRATION
	GtkosxApplication* theApp=g_object_new(GTKOSX_TYPE_APPLICATION, NULL);
	gtkosx_application_set_dock_icon_pixbuf(theApp, icon_main);
	gtkosx_application_ready(theApp);
	g_signal_connect(theApp, "NSApplicationWillTerminate", G_CALLBACK(mac_terminate), &sock);
#endif
	thread_new(listen_draw, argv[1]);
	//g_thread_new("listen_draw", (GThreadFunc)listen_draw, NULL);
#if GTK_MAJOR_VERSION<4
	create_window(NULL);
	gtk_main();
#else
	GtkApplication *app;
  	int status;
  	app = gtk_application_new ("maos.drawdaemon", G_APPLICATION_FLAGS_NONE);
  	g_signal_connect (app, "activate", G_CALLBACK (activate), NULL);
  	status = g_application_run (G_APPLICATION (app), argc, argv);
  	g_object_unref (app);
	return status;
#endif
}/*main */

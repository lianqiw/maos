/*
  Copyright 2009-2026 Lianqi Wang
  
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
#include <unistd.h>
//int sock;
const char *host=0;
#if MAC_INTEGRATION
#include <gtkosxapplication.h>
#endif
#if GTK_MAJOR_VERSION<4 || MAC_INTEGRATION
GdkPixbuf* icon_main=NULL;
#endif
GdkPixbuf *icon_log=NULL;
GdkPixbuf *icon_avg=NULL;
/**
   Parse drawdaemon command line arguments and start the listen_draw() thread.
   argv[1] is the socket or hostname to connect to (localhost if absent), argv[2]
   is the optional draw_id. When MAOS_DIRECT_LAUNCH is set, stdout and stderr are
   redirected to TEMP/drawdaemon.log.

   @param[in] argc Number of command line arguments.
   @param[in] argv Command line arguments.
*/
void handle_args(int argc, char *argv[]){
	if(argc==1){
		info("Usage: %s socket or hostname. Assume localhost if not supplied.\n", argv[0]);
	}
	const char *target=argc>1?argv[1]:"localhost";
	info("drawdaemon is launched with %s %s\n", argv[0], target);
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
		setbuf(stderr, NULL);
		register_signal_handler(NULL);
	}
	if(argc>2){//2nd argument is draw_id
		char *end;
		extern int draw_id;
		draw_id=strtol(argv[2], &end, 10);
		if(argv[2]==end){
			draw_id=0;
		}
		info_time("draw_id=%d\n",draw_id);
	}
	thread_new(listen_draw, (void*)target);
}
#if USE_APPLICATION
/**
   Handler of the GtkApplication activate signal. Creates the main window.

   @param[in] app The application.
   @param[in] user_data Unused.
*/
static void
activate (GtkApplication *app,
          gpointer        user_data)
{
	(void)user_data;
	GtkWidget *window = gtk_application_window_new (app);
	create_window(window);
}
/**
   Handler of the GApplication command-line signal. Parses the arguments and
   activates the application.

   @param[in] application The application.
   @param[in] cmdline The command line.
   @return 0.
*/
static int
command_line(GApplication *application,
			  GApplicationCommandLine *cmdline){
	gchar **argv;
	gint argc;
	argv=g_application_command_line_get_arguments(cmdline, &argc);
	handle_args(argc, argv);
	g_application_activate(application);//this activates the GUI which calls create_window
	return 0;
}
#endif

#if MAC_INTEGRATION
/**
   Handler of the macOS NSApplicationWillTerminate signal. Closes the socket
   pointed to by psock (if valid) and quits the GTK main loop.

   @param[in] app Unused.
   @param[in] psock Pointer to the socket (int) to close. May be NULL.
*/
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
/**
   Drawdaemon: GUI that receives and plots data sent by maos, drawres, drawbin
   or other clients via draw(). Usage: drawdaemon [socket|hostname] [draw_id].

   @param[in] argc Number of command line arguments.
   @param[in] argv Command line arguments.
   @return Exit status.
*/
int main(int argc, char* argv[]){
	
#if GLIB_MAJOR_VERSION<3 && GLIB_MINOR_VERSION<32
	if(!g_thread_supported()){
		g_thread_init(NULL);
		gdk_threads_init();
	}
#endif
#if USE_APPLICATION==0
	gtk_init_check(&argc, &argv);
#endif
	handle_args(argc, argv);
#if USE_APPLICATION==0 || MAC_INTEGRATION
	icon_main=gdk_pixbuf_new_from_resource("/maos/icon-draw.png", NULL);
#endif
	icon_log=gdk_pixbuf_new_from_resource_at_scale("/maos/icon-log.png", -1, 16, 1, NULL);
	icon_avg=gdk_pixbuf_new_from_resource_at_scale("/maos/icon-avg.png", -1, 16, 1, NULL);

#if MAC_INTEGRATION
	GtkosxApplication* theApp=(GtkosxApplication*)g_object_new(GTKOSX_TYPE_APPLICATION, NULL);
	gtkosx_application_set_dock_icon_pixbuf(theApp, icon_main);
	gtkosx_application_ready(theApp);
	g_signal_connect(theApp, "NSApplicationWillTerminate", G_CALLBACK(mac_terminate), &sock);
#endif
	//g_thread_new("listen_draw", (GThreadFunc)listen_draw, NULL);
#if USE_APPLICATION==0
	create_window(NULL);
	gtk_main();
#else
	GtkApplication *app=gtk_application_new("maos.drawdaemon", (GApplicationFlags)0);
	g_signal_connect(app, "activate", G_CALLBACK(activate), NULL);
	g_signal_connect(app, "command-line", G_CALLBACK(command_line), NULL);
  	int status = g_application_run (G_APPLICATION (app), 0, NULL);
	info("status=%d\n", status);
  	g_object_unref (app);
	return status;
#endif
}/*main */

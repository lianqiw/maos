/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
  A monitor that monitors the progress of all running processes in all machines.

  The widgets duplicates and manages the string.
  Set NOTIFY_IS_SUPPORTED to 0 is libnotify is not installed.

  Todo:
  1) Log the activities into a human readable file
  2) Store the activity list into file when quit and reload later
*/

#include <errno.h>
//#include <glib.h>
//#include <gtk/gtk.h>
//#include <glib/gprintf.h>
#include <pthread.h>
#include <fcntl.h>

#include <ctype.h>
#include <sys/types.h>
#include <sys/socket.h>
#include "monitor.h"


#if WITH_NOTIFY
#include <libnotify/notify.h>
static int notify_daemon=1;
#endif
#ifndef GTK_WIDGET_VISIBLE
#define GTK_WIDGET_VISIBLE gtk_widget_get_visible
#endif
int sock_main[2]={0,0}; /*Write by GUI thread. Read by listen_host*/
GdkPixbuf* icon_main=NULL;
GdkPixbuf *icon_draw=NULL;
GdkPixbuf* icon_finished=NULL;
GdkPixbuf* icon_failed=NULL;
GdkPixbuf* icon_running=NULL;
GdkPixbuf* icon_waiting=NULL;
GdkPixbuf* icon_cancel=NULL;
GdkPixbuf* icon_save=NULL;
GdkPixbuf* icon_skip=NULL;
GdkPixbuf* icon_clear=NULL;
GdkPixbuf* icon_connect=NULL;
GtkWidget* notebook=NULL;
GtkWidget** pages=NULL;
#if GTK_MAJOR_VERSION<4
GtkWidget** buttons=NULL;
#endif
GtkWidget* window=NULL;
//static GtkWidget** tabs;
static GtkWidget** titles;
//static GtkWidget** cmdconnect;
//GtkWidget *textscroll=NULL;
//GtkTextBuffer* textbuffer=NULL;
double* usage_cpu, * usage_cpu2;
//double *usage_mem, *usage_mem2;
static GtkWidget** prog_cpu;
//static GtkWidget **prog_mem;
PangoAttrList* pango_active, * pango_down;
#if GTK_MAJOR_VERSION<3
static GtkStatusIcon* status_icon=0;
GdkColor blue;
GdkColor green;
GdkColor red;
GdkColor yellow;
GdkColor white;
GdkColor color_even;
GdkColor color_odd;
#endif
GtkWidget* toolbar;
#if GTK_MAJOR_VERSION>=3
#include "gtk3-css.h"
#endif
#define MAX_HOST 20
const char *mailto=NULL;
int headless=0;//set to 1 when GTK cannot create window
/*
#define DIALOG_MSG(A...) {				\
	GtkWidget *dialog0=gtk_message_dialog_new	\
	    (GTK_WINDOW(window),		\
	     GTK_DIALOG_DESTROY_WITH_PARENT,		\
	     GTK_MESSAGE_INFO,				\
	     GTK_BUTTONS_CLOSE,				\
	     A);					\
	gtk_dialog_run(GTK_DIALOG(dialog0));		\
	gtk_widget_destroy(dialog0);			\
    }
	*/

void scheduler_cmd_wrap(int ihost, int pid, int command){
	int cmd[4]={MON_CMD, ihost, pid, command};
	stwriteintarr(sock_main[1], cmd, 4);
	//the following may block the GUI. 
	//int ans;
	//streadint(sock_main[1], &ans);
	//return ans;
}
void add_host_wrap(int ihost){
	int cmd[3]={MON_ADDHOST, ihost, 0};
	stwriteintarr(sock_main[1], cmd, 3);
}
void clear_job_wrap(int ihost, int flag){
	int cmd[3]={MON_CLEARJOB, ihost, flag};
	stwriteintarr(sock_main[1], cmd, 3);
}
void kill_job_wrap(int ihost, int pid){
	int cmd[3]={MON_KILLJOB, ihost, pid};
	stwriteintarr(sock_main[1], cmd, 3);
}
void save_job_wrap(){
	int cmd[3]={MON_SAVEJOB, 0, 0};
	stwriteintarr(sock_main[1], cmd, 3);
}
/**
   The number line pattern determines how dash is drawn for gtktreeview. the
   first number is the length of the line, and second number of the length
   of blank after the line.

   properties belonging to certain widget that are descendent of GtkWidget, need
to specify the classname before the key like GtkTreeView::allow-rules */
#if GTK_MAJOR_VERSION<3
static const gchar* rc_string_widget=
{
"style \"widget\" {\n"
	"font_name = \"Sans 12\""
"}\n"
"class \"GtkWidget\" style \"widget\" \n"
};
static const gchar* rc_string_treeview=
{
"style \"solidTreeLines\"{                       \n"
/*" GtkTreeView::grid-line-pattern=\"\111\111\"\n" */
" GtkTreeView::grid-line-width   = 1             \n"
" GtkTreeView::allow-rules       = 1             \n"
/*" GtkTreeView::odd-row-color     = \"#EFEFEF\"   \n"
  " tkTreeView::even-row-color    = \"#FFFFFF\"   \n"*/
" GtkTreeView::horizontal-separator = 0          \n"
" GtkTreeView::vertical-separator = 0            \n"
"}\n                                             \n"
"class \"GtkTreeView\" style \"solidTreeLines\"  \n"
};

static const gchar* rc_string_entry=
{
"style \"entry\" {               \n"
/*"base[NORMAL] = \"white\"        \n"
  "bg[SELECTED] = \"#0099FF\"         \n"*/
"xthickness   = 0                \n"
"ythickness   = 0                \n"
"GtkEntry::inner-border    = {1,1,1,1}     \n"
"GtkEntry::progress-border = {0,0,0,0}     \n"
"GtkEntry::has-frame       = 0 \n"
"}\n"
"class \"GtkEntry\" style \"entry\" \n"
};
#endif
/**
   Enables the notebook page for a host*/
gboolean host_up(gpointer data){
	int ihost=GPOINTER_TO_INT(data);
	//gtk_widget_set_sensitive(cmdconnect[ihost], 0);
	//gtk_widget_hide(cmdconnect[ihost]);
	gtk_widget_set_sensitive(pages[ihost], 1);
	gtk_label_set_attributes(GTK_LABEL(titles[ihost]), pango_active);
	return 0;
}
/**
   Disables the notebook page for a host*/
gboolean host_down(gpointer data){
	int ihost=GPOINTER_TO_INT(data);
	char status[100];
	snprintf(status, sizeof(status), "Disconnected at %s", myasctime(0));
	/*gtk_button_set_label(GTK_BUTTON(cmdconnect[ihost]), status);
	gtk_widget_show(cmdconnect[ihost]);
	gtk_widget_set_sensitive(cmdconnect[ihost], 1);*/
	gtk_widget_set_sensitive(pages[ihost], 0);
	gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(prog_cpu[ihost]), 0);
	//gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(prog_mem[ihost]), 0);
	gtk_label_set_attributes(GTK_LABEL(titles[ihost]), pango_down);
	return 0;
}
/**
* determine host index from name
*/
int host2i(const char *hostn){
	for(int ihost=0; ihost<nhost; ihost++){
		if(!strcmp(hostshort[ihost], hostn)){
			return ihost;
		}
	}
	warning("host %s is not found\n", hostn);
	return -1;
}
/**
   modifies the color of progress bar*/
   /*
static void modify_bg(GtkWidget* widget, int type){
#if GTK_MAJOR_VERSION>=3
	//GtkCssProvider* provider;
	const char *class[]={NULL,NULL};
	switch(type){
	case 1:
		class[0]="blue";
		break;
	case 2:
		class[0]="red";
		break;
	default:
		class[0]=NULL;
	}
	if(class[0]){
#if GTK_MAJOR_VERSION>=4
		gtk_widget_set_css_classes(widget, class);
#else
		//gtk_widget_set_name(widget, class[0]);
		GtkStyleContext *context=gtk_widget_get_style_context(widget);
		gtk_style_context_add_class(context, class[0]);
#endif
	}
#else
	GdkColor* color;
	switch(type){
	case 1:
		color=&blue;
		break;
	case 2:
		color=&red;
		break;
	case 3:
		color=&green;
		break;
	default:
		color=NULL;
	}
	gtk_widget_modify_bg(widget, GTK_STATE_SELECTED, color);
	gtk_widget_modify_bg(widget, GTK_STATE_PRELIGHT, color);
#endif
}*/
/**
   updates the progress bar for a job*/
gboolean update_progress(gpointer input){
	int ihost=GPOINTER_TO_INT(input);
	//double last_cpu=usage_cpu2[ihost];
	//double last_mem=usage_mem2[ihost];
	usage_cpu2[ihost]=usage_cpu[ihost];
	//usage_mem2[ihost]=usage_mem[ihost];
	if(GTK_WIDGET_VISIBLE(window)){
		gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(prog_cpu[ihost]), usage_cpu[ihost]);
		//gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(prog_mem[ihost]), usage_mem[ihost]);
		/*if(usage_cpu[ihost]>=0.8&&last_cpu<0.8){
			modify_bg(prog_cpu[ihost], 2);
		} else if(usage_cpu[ihost]<0.8&&last_cpu>=0.8){
			modify_bg(prog_cpu[ihost], 1);
		}*/
			/*
		if(usage_mem[ihost]>=0.8 && last_mem<0.8){
			modify_bg(prog_mem[ihost],2);
		}else if(usage_mem[ihost]<0.8 && last_mem>=0.8){
			modify_bg(prog_mem[ihost],1);
				}*/
	}
	return 0;
}

void notify_user(proc_t* p){
	if(p->status.info==p->oldinfo || p->status.done) return;
	p->oldinfo=p->status.info;
#if WITH_NOTIFY
	if(!notify_daemon) return;

	static NotifyNotification* notify_urgent=NULL, * notify_normal=NULL, * notify_low=NULL;
	if(!notify_urgent){
#if !defined(NOTIFY_CHECK_VERSION) || !NOTIFY_CHECK_VERSION(0,7,0)
	/*newer versions doesnot have _new_with_status_icon */
		notify_low=notify_notification_new_with_status_icon("Low", NULL, NULL, status_icon);
		notify_normal=notify_notification_new_with_status_icon("Normal", NULL, NULL, status_icon);
		notify_urgent=notify_notification_new_with_status_icon("Urgent", NULL, "error", status_icon);
#endif
		if(!notify_low){
			notify_low=notify_notification_new("Low", NULL, NULL);
			notify_normal=notify_notification_new("Normal", NULL, NULL);
			notify_urgent=notify_notification_new("Urgent", NULL, NULL);
		}
		notify_notification_set_icon_from_pixbuf(notify_low, icon_main);
		notify_notification_set_timeout(notify_low, NOTIFY_EXPIRES_DEFAULT);
		notify_notification_set_urgency(notify_low, NOTIFY_URGENCY_LOW);

		notify_notification_set_icon_from_pixbuf(notify_normal, icon_main);
		notify_notification_set_timeout(notify_normal, NOTIFY_EXPIRES_DEFAULT);
		notify_notification_set_urgency(notify_normal, NOTIFY_URGENCY_NORMAL);

		notify_notification_set_timeout(notify_urgent, NOTIFY_EXPIRES_NEVER);
		notify_notification_set_urgency(notify_urgent, NOTIFY_URGENCY_CRITICAL);

	}
	static char summary[80];
	NotifyNotification* notify;
	switch(p->status.info){
	case S_START:
		notify=notify_low;
		snprintf(summary, 80, "Job Started on %s", hosts[p->hid]);
		break;
	case S_FINISH:
		notify=notify_normal;
		snprintf(summary, 80, "Job Finished on %s", hosts[p->hid]);
		break;
	case S_CRASH:
		notify=notify_urgent;
		snprintf(summary, 80, "Job Crashed on %s!", hosts[p->hid]);
		break;
	case S_KILLED:
		notify=notify_urgent;
		snprintf(summary, 80, "Job is Killed on %s!", hosts[p->hid]);
		break;
	default:
		return;
	}
	notify_notification_update(notify, summary, p->path, NULL);
	notify_notification_show(notify, NULL);

#endif
}
/**
   quite the program
*/
static void quitmonitor(GtkWidget* widget, gpointer data){
	(void)widget;
	(void)data;
#if WITH_NOTIFY
	if(notify_daemon)
		notify_uninit();
#endif
	add_host_wrap(-2);
	/*if(gtk_main_level()>0){
		gtk_main_quit();
	}*/
	exit(0);
}
/**
   update the job count
*/
gboolean update_title(gpointer data){
	size_t tmp=GPOINTER_TO_SIZE(data);
	int hid=tmp&0xFF;//lowest 1 byte
	int nproc=tmp>>8&0xFFF;//2-3 byte
	int npending=tmp>>20;
	char tit[40];
	if(hid<nhost){
		if(nproc>0){
			snprintf(tit, 40, "%s (%d/%d)", hostshort[hid], npending, nproc);
		}else{
			snprintf(tit, 40, "%s", hostshort[hid]);
		}
	} else{
		snprintf(tit, 40, "All");
		gtk_label_set_attributes(GTK_LABEL(titles[hid]), pango_active);
	}
	gtk_label_set_text(GTK_LABEL(titles[hid]), tit);
	return 0;//always return 0 to avoid repeat itself when use g_idle_add
}
/**
   respond to the kill job event
*/
/*
void kill_job(int hid, int pid){
	GtkWidget* dia=gtk_message_dialog_new
	(GTK_WINDOW(window), GTK_DIALOG_DESTROY_WITH_PARENT,
		GTK_MESSAGE_QUESTION,
		GTK_BUTTONS_NONE,
		"Kill job %d on server %s?",
		pid, hosts[hid]);
	if(p->status.info==S_WAIT){
		gtk_dialog_add_buttons(GTK_DIALOG(dia), "Kill Job", 0, "Cancel", 2, NULL);
	} else{
		gtk_dialog_add_buttons(GTK_DIALOG(dia), "Kill Job", 0, "Remote Display", 1, "Cancel", 2, NULL);
	}
	int result=gtk_dialog_run(GTK_DIALOG(dia));
	gtk_widget_destroy(dia);
	switch(result){
	case 0:
		if(scheduler_cmd(p->hid, p->pid, CMD_KILL)){
			warning("Failed to kill the job\n");
		}

		refresh(p);
		break;
	case 1:
	{
		scheduler_display(hid, pid);
	}
	break;
	}
}*/
void kill_all_jobs_do(const int this_host){
	for(int ihost=0; ihost<nhost; ihost++){
		if(this_host!=nhost&&ihost!=this_host){
			continue;
		}
		kill_job_wrap(ihost, 0);
	}
}
#if GTK_MAJOR_VERSION>=4
void kill_all_job_callback(GtkDialog *dialog, int result, gpointer data){
	gtk_window_destroy(GTK_WINDOW(dialog));
	int this_host=GPOINTER_TO_INT(data);
	if(result){
		kill_all_jobs_do(this_host);
	}
}
#endif
void kill_all_jobs(GtkWidget* btn, gpointer data){
	(void)btn;
	(void)data;
	int this_host=gtk_notebook_get_current_page(GTK_NOTEBOOK(notebook));

	GtkWidget* dia=gtk_message_dialog_new
	(GTK_WINDOW(window), GTK_DIALOG_DESTROY_WITH_PARENT,
		GTK_MESSAGE_QUESTION,
		GTK_BUTTONS_NONE,
		"Kill all jobs on %s?", this_host==nhost?"all servers":hostshort[this_host]);
	gtk_dialog_add_buttons(GTK_DIALOG(dia), "Kill all", 1, "Cancel", 0, NULL);
#if GTK_MAJOR_VERSION>=4
	g_signal_connect(GTK_DIALOG(dia), "response", G_CALLBACK(kill_all_job_callback), GINT_TO_POINTER(this_host));
	gtk_window_present(GTK_WINDOW(dia));
#else
	int result=gtk_dialog_run(GTK_DIALOG(dia));
	gtk_widget_destroy(dia);
	if(result){
		kill_all_jobs_do(this_host);
	}
#endif
}
#if GTK_MAJOR_VERSION < 3
static void status_icon_on_click(void* widget,
	gpointer data){
	(void)widget;
	(void)data;
	static int cx=0, cy=0;
	static int x=0, y=0;
	int force_show=GPOINTER_TO_INT(data);
	if(GTK_WIDGET_VISIBLE(window)&&!force_show){
		gtk_window_get_size(GTK_WINDOW(window), &x, &y);
		gtk_window_get_position(GTK_WINDOW(window), &cx, &cy);
		gtk_widget_hide(window);
	} else{
		if(x&&y){
			gtk_window_set_default_size(GTK_WINDOW(window), x, y);
			gtk_window_move(GTK_WINDOW(window), cx, cy);
		}
		gtk_widget_show(window);
		//gtk_window_deiconify(GTK_WINDOW(window));
	}
}
static void status_icon_on_popup(GtkStatusIcon* status_icon0, guint button,
	guint32 activate_time, gpointer user_data){
	(void)status_icon0;
	(void)button;
	(void)activate_time;
	/*gtk_menu_popup(GTK_MENU(user_data),NULL,NULL,
		   gtk_status_icon_position_menu,
		   status_icon0,button,activate_time);*/
	status_icon_on_click(NULL, user_data);
}
#endif
static void create_status_icon(){
#if GTK_MAJOR_VERSION < 3
	status_icon=gtk_status_icon_new_from_pixbuf(icon_main);
	//g_signal_connect(GTK_STATUS_ICON(status_icon),"activate", G_CALLBACK(status_icon_on_click),0);//useless.
	g_signal_connect(GTK_STATUS_ICON(status_icon), "popup-menu", G_CALLBACK(status_icon_on_popup), 0);
#endif
	/*GtkWidget *menu, *menuItemShow,*menuItemExit;
	menu=gtk_menu_new();
	menuItemShow=gtk_menu_item_new_with_label("Show/Hide");
	menuItemExit=gtk_menu_item_new_with_mnemonic("_Exit");
	g_signal_connect(G_OBJECT(menuItemShow),"activate",
			 G_CALLBACK(status_icon_on_click),NULL);
	g_signal_connect(G_OBJECT(menuItemExit),"activate",
			 G_CALLBACK(quitmonitor),NULL);
	gtk_menu_shell_append(GTK_MENU_SHELL(menu),menuItemShow);
	gtk_menu_shell_append(GTK_MENU_SHELL(menu),menuItemExit);
	gtk_widget_show_all(menu);
	g_signal_connect(GTK_STATUS_ICON(status_icon),
			 "popup-menu",G_CALLBACK(trayIconPopup),menu);
#if GTK_MAJOR_VERSION >=3 || GTK_MINOR_VERSION>=16
	gtk_status_icon_set_tooltip_text(status_icon, "MAOS Job Monitoring");
#endif
	gtk_status_icon_set_visible(status_icon, TRUE);*/
}

/*
static gboolean delete_window(void){
	if(status_icon && gtk_status_icon_is_embedded(status_icon)){
	gtk_widget_hide(window);
	return TRUE;//do not quit
	}else{
		return FALSE;//quit.
	}
}*/
/*
static gboolean
window_state_event(GtkWidget *widget,GdkEventWindowState *event,gpointer data){
	(void)data;
	if(event->changed_mask == GDK_WINDOW_STATE_ICONIFIED
	   && (event->new_window_state == GDK_WINDOW_STATE_ICONIFIED
	   || event->new_window_state ==
	   (GDK_WINDOW_STATE_ICONIFIED | GDK_WINDOW_STATE_MAXIMIZED)))
	{
		gtk_widget_hide (GTK_WIDGET(widget));
	}
	return TRUE;
	}*/

void clear_jobs(GtkWidget* btn, gpointer flag){
	(void)btn;
	int this_host=gtk_notebook_get_current_page(GTK_NOTEBOOK(notebook));
	if(this_host==nhost){
		for(int ihost=0; ihost<nhost; ihost++){
			clear_job_wrap(ihost, GPOINTER_TO_INT(flag));
		}
	}else{
		clear_job_wrap(this_host, GPOINTER_TO_INT(flag));
	}
}

void save_all_jobs(GtkWidget* btn, gpointer data){
	(void) btn;
	(void) data;
	save_job_wrap();
}
GtkWidget* monitor_new_entry_progress(void){
	GtkWidget* prog=gtk_entry_new();
	gtk_editable_set_editable(GTK_EDITABLE(prog), FALSE);
#if GTK_MAJOR_VERSION>=4
	gtk_editable_set_width_chars(GTK_EDITABLE(prog), 12);
#else
	gtk_entry_set_width_chars(GTK_ENTRY(prog), 12);
#endif
#if GTK_VERSION_AFTER(2,18)
	gtk_widget_set_can_focus(prog, FALSE);//2.18
#else
	g_object_set(prog, "can-focus", FALSE, NULL);
#endif

	/*gtk_widget_modify_base(prog,GTK_STATE_NORMAL, &white);
	gtk_widget_modify_bg(prog,GTK_STATE_SELECTED, &blue);
	*/
#if GTK_VERSION_AFTER(2,10) && !GTK_VERSION_AFTER(3,4)
	gtk_entry_set_inner_border(GTK_ENTRY(prog), 0);//after 2.10, beore 3.4
#endif
	gtk_entry_set_has_frame(GTK_ENTRY(prog), 0);
	gtk_entry_set_alignment(GTK_ENTRY(prog), 0.5);//since 2.4
	return prog;
}
GtkWidget* monitor_new_progress(int vertical, int length){
	GtkWidget* prog=gtk_progress_bar_new();
	gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(prog), 0);
	if(vertical){
#if GTK_MAJOR_VERSION>=3
		gtk_orientable_set_orientation(GTK_ORIENTABLE(prog),
			GTK_ORIENTATION_VERTICAL);
		gtk_progress_bar_set_inverted(GTK_PROGRESS_BAR(prog), TRUE);
		gtk_progress_bar_set_show_text(GTK_PROGRESS_BAR(prog), FALSE);
#else
		gtk_progress_bar_set_orientation(GTK_PROGRESS_BAR(prog),
			GTK_PROGRESS_BOTTOM_TO_TOP);
		//g_object_set(G_OBJECT(prog), "show-text", FALSE, NULL);
#endif
		gtk_widget_set_size_request(prog, 4,4);
	} else{
		gtk_widget_set_size_request(prog, length, 12);
	}
	return prog;
}
/*
  Try to connect to hosts in response to button click.
 */
static void add_host_event(GtkWidget* button, gpointer data){
	(void)button;
	int ihost=GPOINTER_TO_INT(data);
	if(ihost>-1&&ihost<nhost){
		add_host_wrap(ihost);//only this host
	} else{
		for(ihost=0; ihost<nhost; ihost++){
			add_host_wrap(ihost);//all host
		}
	}
}
#if USE_TITLEBAR
static void switch_page(GtkWidget* button, gpointer data){
	int ihost=GPOINTER_TO_INT(data);
	int ihost_old=gtk_notebook_get_current_page(GTK_NOTEBOOK(notebook));
	if(ihost==ihost_old){//keep active
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
	}else if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button))){
		gtk_notebook_set_current_page(GTK_NOTEBOOK(notebook), ihost);
	#if GTK_MAJOR_VERSION<4
		for(int jhost=0; jhost<=nhost; jhost++){
			if(jhost!=ihost){
				if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(buttons[jhost]))){
					gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(buttons[jhost]), FALSE);
				}
			}
		}
	#endif	
	}
}
#endif
gboolean host_added_wrap(gpointer data){
	(void) data;
	add_host_wrap(-1);//just wakeup respond() with noop.
	return FALSE;//just run once.
}
void parse_icons(){
#if GDK_MAJOR_VERSION < 3
	gdk_color_parse("#EE0000", &red);
	gdk_color_parse("#00CC00", &green);
	gdk_color_parse("#FFFF33", &yellow);
	gdk_color_parse("#FFFFFF", &white);
	gdk_color_parse("#0099FF", &blue);
	gdk_color_parse("#FFFFAB", &color_even);
	gdk_color_parse("#FFFFFF", &color_odd);
#endif
	icon_main=gdk_pixbuf_new_from_resource("/maos/icon-monitor.png", NULL);
	icon_draw=gdk_pixbuf_new_from_resource("/maos/icon-draw-small.png", NULL);
	icon_finished=gdk_pixbuf_new_from_resource("/maos/icon-finished.png", NULL);
	icon_failed=gdk_pixbuf_new_from_resource("/maos/icon-error.png", NULL);
	icon_running=gdk_pixbuf_new_from_resource("/maos/icon-play.png", NULL);
	icon_skip=gdk_pixbuf_new_from_resource("/maos/icon-skip.png", NULL);
	icon_waiting=gdk_pixbuf_new_from_resource("/maos/icon-waiting.png", NULL);
	icon_cancel=gdk_pixbuf_new_from_resource("/maos/icon-cancel.png", NULL);
	icon_save=gdk_pixbuf_new_from_resource("/maos/icon-save.png", NULL);
	icon_clear=gdk_pixbuf_new_from_resource("/maos/icon-clear.png", NULL);
	icon_connect=gdk_pixbuf_new_from_resource("/maos/icon-connect.png", NULL);
}
void parse_provider(){

#if GTK_MAJOR_VERSION<3
	gtk_rc_parse_string(rc_string_widget);
	gtk_rc_parse_string(rc_string_treeview);
	gtk_rc_parse_string(rc_string_entry);
	//    GtkStyle *style=gtk_widget_get_style(window);
#else
	//	copied from drawdaemon_gui.c
	GtkCssProvider*provider_default=gtk_css_provider_new();
#if GTK_MAJOR_VERSION>=4
	gtk_css_provider_load_from_string(provider_default, all_style);
	gtk_style_context_add_provider_for_display(gdk_display_get_default(), GTK_STYLE_PROVIDER(provider_default),
		GTK_STYLE_PROVIDER_PRIORITY_USER);
	GtkSettings *settings=gtk_settings_get_for_display(gdk_display_get_default());
	g_object_set(G_OBJECT(settings), "gtk-hint-font-metrics", TRUE, NULL);//fix gtk4 subpixel positioning blurry.
#else
	gtk_css_provider_load_from_data(provider_default, all_style, strlen(all_style), NULL);

	GtkStyleContext *all_context=gtk_widget_get_style_context(window);
	GdkScreen *screen=gtk_style_context_get_screen(all_context);
	gtk_style_context_add_provider_for_screen(screen, GTK_STYLE_PROVIDER(provider_default),
		GTK_STYLE_PROVIDER_PRIORITY_USER);
#endif
#endif
}

void create_window(
#if USE_APPLICATION
	GtkApplication* app,
	gpointer        user_data
#endif
		  ){
#if USE_APPLICATION
	window=gtk_application_window_new(app);
	(void)user_data;
#else
	window=gtk_window_new(GTK_WINDOW_TOPLEVEL);
	gtk_window_set_position(GTK_WINDOW(window), GTK_WIN_POS_CENTER);
#endif	
#if GTK_VERSION_AFTER(3, 14)
#if GTK_VERSION_AFTER(4, 0)
	GtkIconTheme *icon_theme=gtk_icon_theme_get_for_display(gtk_widget_get_display(window));
#else	
	GtkIconTheme *icon_theme=gtk_icon_theme_get_default();
#endif
	gtk_icon_theme_add_resource_path(icon_theme, "/maos/");
	/*if(gtk_icon_theme_has_icon(icon_theme, "icon-monitor")!=1){
		warning("gtk_icon_theme does not have icon-monitor\n");
	}
	if(gtk_icon_theme_has_icon(icon_theme, "video-display")!=1){
		warning("gtk_icon_theme does not have video-display\n");
	}*/
	//gtk_window_set_icon_name(GTK_WINDOW(window), "icon-monitor");
	gtk_window_set_default_icon_name("icon-monitor");//this does not work in macos as quartz backend does not support it. 
	//info_time("gtk_icon_theme_has_icon returned %d for icon-monitor\n", gtk_icon_theme_has_icon(icon_theme, "icon-monitor"));//test ok.
	//gtk_window_set_icon_name(GTK_WINDOW(window), "video-display");
	//icon_name is ignored on all non-x11 systems (wayland, windows, macos).
#else
	gtk_window_set_icon(GTK_WINDOW(window), icon_main);
#endif
	gtk_window_set_resizable(GTK_WINDOW(window), TRUE);
	parse_provider();//requires window to be set
	gtk_window_set_title(GTK_WINDOW(window), "MAOS Monitor");
#if USE_TITLEBAR
	toolbar=gtk_header_bar_new();
#if GTK_MAJOR_VERSION<4
	gtk_header_bar_set_show_close_button(GTK_HEADER_BAR(toolbar), TRUE);
	gtk_header_bar_set_title(GTK_HEADER_BAR(toolbar), NULL);
#else
	gtk_header_bar_set_show_title_buttons(GTK_HEADER_BAR(toolbar), TRUE);
#endif
#else
	//toolbar=gtk_toolbar_new();
	//gtk_toolbar_set_icon_size(GTK_TOOLBAR(toolbar),GTK_ICON_SIZE_MENU);
	toolbar=gtk_hbox_new(FALSE, 0);
#endif
	new_toolbar_item(toolbar, NULL, 0,"computer", icon_connect, "Connect", add_host_event, GINT_TO_POINTER(-1));
	new_toolbar_item(toolbar, NULL, 0, NULL, NULL, NULL, NULL, NULL);
	new_toolbar_item(toolbar, NULL, 0, "object-select", icon_finished, "Clear finished jobs", clear_jobs, GINT_TO_POINTER(-1));
	new_toolbar_item(toolbar, NULL, 0, "media-skip-forward", icon_skip, "Clear skipped jobs", clear_jobs, GINT_TO_POINTER(-2));
	new_toolbar_item(toolbar, NULL, 0, "dialog-error", icon_failed, "Clear crashed jobs", clear_jobs, GINT_TO_POINTER(-3));
	new_toolbar_item(toolbar, NULL, 0, "edit-clear-all", icon_clear, "Clear all jobs", clear_jobs, GINT_TO_POINTER(-4));
	new_toolbar_item(toolbar, NULL, 0, NULL, NULL, NULL, NULL, NULL);
	new_toolbar_item(toolbar, NULL, 0, "process-stop", icon_cancel, "Kill all jobs", kill_all_jobs, GINT_TO_POINTER(-1));
	new_toolbar_item(toolbar, NULL, 0, "media-floppy", icon_save, "Save jobs to file", save_all_jobs, GINT_TO_POINTER(-1));
	new_toolbar_item(toolbar, NULL, 0, NULL, NULL, NULL, NULL, NULL);
#if GTK_MAJOR_VERSION<4
	gtk_widget_show_all(toolbar);
#endif
	notebook=gtk_notebook_new();
#if USE_TITLEBAR
	gtk_notebook_set_show_tabs(GTK_NOTEBOOK(notebook), FALSE);
#else
	gtk_notebook_set_tab_pos(GTK_NOTEBOOK(notebook), GTK_POS_TOP);
#endif
	//gtk_notebook_set_scrollable(GTK_NOTEBOOK(notebook), TRUE);
#if USE_TITLEBAR
	gtk_window_set_titlebar(GTK_WINDOW(window), toolbar);
#else
	gtk_notebook_set_action_widget(GTK_NOTEBOOK(notebook), toolbar, GTK_PACK_START);
#endif
	//gtk_widget_show(notebook);
	//box_append(GTK_BOX(vbox), notebook, TRUE, TRUE, 0);

	//GtkWidget *vbox=gtk_vbox_new(FALSE, 0);
#if GTK_MAJOR_VERSION<4
	gtk_container_add(GTK_CONTAINER(window), notebook);
#else
	gtk_window_set_child(GTK_WINDOW(window), notebook);
#endif

	/*{//text are to show details job information
		//GtkWidget* seperator=gtk_hseparator_new();
		//box_append(GTK_BOX(pages[ihost]), seperator, FALSE, FALSE, 0);
		GtkWidget* textview=gtk_text_view_new();
		textbuffer=gtk_text_view_get_buffer(GTK_TEXT_VIEW(textview));
		textscroll=gtk_scrolled_window_new(
#if GTK_MAJOR_VERSION<4
				NULL, NULL
#endif
		);
		gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(textscroll),
			GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
		gtk_widget_set_size_request(textscroll, -1, 50);
#if GTK_MAJOR_VERSION<4
		gtk_container_add(GTK_CONTAINER(textscroll), textview);//page is scrollable.
#else
		gtk_widget_set_hexpand(textscroll, TRUE);//important
		gtk_widget_set_vexpand(textscroll, TRUE);
		gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(textscroll), textview);
#endif
		box_append(GTK_BOX(vbox), textscroll, TRUE, TRUE, 0);
		gtk_text_buffer_set_text(textbuffer, "", -1);
		gtk_text_view_set_editable(GTK_TEXT_VIEW(textview), 0);
		gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(textview), GTK_WRAP_WORD);
	}*/

	//g_signal_connect(window, "delete_event", G_CALLBACK (delete_window), NULL);
	g_signal_connect(window, "destroy", G_CALLBACK(quitmonitor), NULL);
	//g_signal_connect(G_OBJECT (window), "window-state-event", G_CALLBACK (window_state_event), NULL);

	gtk_window_set_default_size(GTK_WINDOW(window), 1200, 800);

	//tabs=mycalloc(nhost+1, GtkWidget*);
	pages=mycalloc(nhost+1, GtkWidget*);
	titles=mycalloc(nhost+1, GtkWidget*);
#if GTK_MAJOR_VERSION<4
	buttons=mycalloc(nhost+1, GtkWidget*);
#endif	
	//cmdconnect=mycalloc(nhost+1, GtkWidget*);
	//buffers=mycalloc(nhost+1, GtkTextBuffer*);

	
	prog_cpu=mycalloc(nhost, GtkWidget*);
	//prog_mem=mycalloc(nhost,GtkWidget *);

	pango_active=pango_attr_list_new();
	pango_down=pango_attr_list_new();
	pango_attr_list_insert(pango_down, pango_attr_foreground_new(0x88FF, 0x88FF, 0x88FF));
	pango_attr_list_insert(pango_active, pango_attr_foreground_new(0x0000, 0x0000, 0x0000));
#if GTK_MAJOR_VERSION>3 && USE_TITLEBAR
	GtkWidget *lastbtn=NULL;
#endif	
	for(int ihost=0; ihost<=nhost; ihost++){//ihost==nhost is the include all tab
		//char tit[40];
		//snprintf(tit,40,"%s(0)",hosts[ihost]);
		GtkWidget *eventbox=0;
		//create notebook header
		titles[ihost]=gtk_label_new(NULL);
		gtk_label_set_attributes(GTK_LABEL(titles[ihost]), pango_down);
		GtkWidget* hbox0=gtk_hbox_new(FALSE, 0);
		if(ihost<nhost){
#if GTK_MAJOR_VERSION>=3
			prog_cpu[ihost]=monitor_new_progress(1, 4);
			//prog_mem[ihost]=monitor_new_progress(1,4);
#else
			prog_cpu[ihost]=monitor_new_progress(1, 16);
			//prog_mem[ihost]=monitor_new_progress(1,16);
#endif
			box_append(GTK_BOX(hbox0), prog_cpu[ihost], FALSE, FALSE, 1);
			//modify_bg(prog_cpu[ihost], 2);
			//modify_bg(prog_mem[ihost], 1);
		}
		gtk_box_set_homogeneous(GTK_BOX(hbox0), 0);
		box_append(GTK_BOX(hbox0), titles[ihost], FALSE, TRUE, 0);
#if USE_TITLEBAR
		eventbox=gtk_toggle_button_new();
#if GTK_MAJOR_VERSION>3
		if(lastbtn){
			gtk_toggle_button_set_group(GTK_TOGGLE_BUTTON(eventbox), GTK_TOGGLE_BUTTON(lastbtn));
		}
		lastbtn=eventbox;
#endif		
	#if GTK_MAJOR_VERSION<4
		gtk_button_set_relief(GTK_BUTTON(eventbox), GTK_RELIEF_NONE);
		gtk_container_add(GTK_CONTAINER(eventbox), hbox0);
	#else
		gtk_button_set_has_frame(GTK_BUTTON(eventbox), FALSE);
		gtk_button_set_child(GTK_BUTTON(eventbox), hbox0);
	#endif
#else
#if GTK_MAJOR_VERSION>=4
		eventbox=hbox0;
#else		
		eventbox=gtk_event_box_new();//eventbox allows click but does not show click animation
		gtk_event_box_set_above_child(GTK_EVENT_BOX(eventbox), TRUE);
		gtk_event_box_set_visible_window(GTK_EVENT_BOX(eventbox), FALSE);
		gtk_container_add(GTK_CONTAINER(eventbox), hbox0);
#endif		
#endif			
		/*Put the box below an eventbox so that clicking on the progressbar
		switches tab also.*/
#if GTK_MAJOR_VERSION<4
		gtk_widget_show_all(eventbox);
		buttons[ihost]=eventbox;
#endif
		//pages[ihost]=gtk_paned_new(GTK_ORIENTATION_VERTICAL);//gtk_vbox_new(FALSE, 0);
		//gtk_paned_set_position(GTK_PANED(pages[ihost]), 500);
		//area for showing list of jobs
		GtkWidget* page=new_page(ihost);
		GtkWidget* scroll=gtk_scrolled_window_new(
#if GTK_MAJOR_VERSION<4
			NULL, NULL
#endif
			);
		gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroll),
			GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
		gtk_widget_set_size_request(scroll, -1, 300);
#if GTK_MAJOR_VERSION<4
		gtk_container_add(GTK_CONTAINER(scroll), page);//page is scrollable.
#else
		gtk_widget_set_hexpand(scroll, TRUE);//important
		gtk_widget_set_vexpand(scroll, TRUE);
		gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scroll), page);
#endif
		//box_append(GTK_BOX(pages[ihost]), scroll, TRUE, TRUE, 0);
		//gtk_paned_add1(GTK_PANED(pages[ihost]), scroll);
		pages[ihost]=scroll;

		/*if(ihost<nhost){
			cmdconnect[ihost]=gtk_button_new_with_label("Click to connect");
			g_signal_connect(cmdconnect[ihost], "clicked", G_CALLBACK(add_host_event), GINT_TO_POINTER(ihost));
			// button for reconnection
			GtkWidget* hbox=gtk_hbox_new(FALSE, 0);
			box_append(GTK_BOX(hbox), cmdconnect[ihost], TRUE, FALSE, 0);
			box_append(GTK_BOX(pages[ihost]), hbox, FALSE, FALSE, 0);
		}*/
#if USE_TITLEBAR
		new_toolbar_item(toolbar, eventbox, 1, NULL, NULL, "notebook", switch_page, GINT_TO_POINTER(ihost));
		eventbox=0;
#endif	
		//if(ihost<nhost){
		gtk_notebook_append_page(GTK_NOTEBOOK(notebook), pages[ihost], eventbox);
		//} else{
		//	gtk_notebook_insert_page(GTK_NOTEBOOK(notebook), pages[ihost], eventbox, 0);
		//}
		update_title(GINT_TO_POINTER(ihost));
		//gtk_widget_show_all(pages[ihost]);
	}
#if GTK_MAJOR_VERSION<4
	gtk_widget_show_all(window);
#else
	gtk_window_present(GTK_WINDOW(window));
#endif
	gtk_notebook_set_current_page(GTK_NOTEBOOK(notebook), nhost);//need to be after show_all
	//do it here to prevent memory error.
	usage_cpu=mycalloc(nhost, double);
	//usage_mem=mycalloc(nhost,double);
	usage_cpu2=mycalloc(nhost, double);
	//usage_mem2=mycalloc(nhost,double);
	for(int ihost=0; ihost<nhost; ihost++){
		add_host_wrap(ihost);
	}
	//gtk_widget_hide(textscroll);
}

void print_help(const char *cmd){
	fprintf(stderr, "%s [options] [port] [host1] [host2]\n"
	"\nOptions:\n\t\t--mailto someone@somehost.com (requires sendmail)\n"
	"\t\t--disable-plot or --disable-draw Disable plotting capability (drawdaemon)\n"
	"\t\t--help (-h) print this help\n"
	"\ndefault hosts can be specified in ~/.aos/hosts, one per line\n"
	"default port can be specified in ~/.aos/port\n", cmd
	);
}
int plot_enabled=1;
void handle_args(int argc, char *argv[]){
	if(argc==1 && nhost==1){
		print_help(argv[0]);
	}
	for(int i=1; i<argc; i++){
		if(argv[i][0]=='-'){//start with -
			const char *s=argv[i]+1;
			if(s[0]=='-') s++;
			const char *key="mailto";
			const char *key2a="disable-plot";
			const char *key2b="disable-draw";
			if(!mystrcmp(s, key)){
				s+=strlen(key);
				while(s[0]==' ') s++;
				if(s[0]=='=') s++;
				while(s[0]==' ') s++;
				if(s[0]){
					mailto=mystrdup(s);
				} else if(i+1<argc){
					mailto=mystrdup(argv[i+1]);
					i++;
				} else{
					warning("Invalid option: %s\n", argv[i]);
				}
				if(mailto){
					info("Will send mail to %s for host disconnection or job crashing.\n", mailto);
				}
			} else if(!mystrcmp(s, key2a)||!mystrcmp(s, key2b)){
				plot_enabled=0;
				warning("plot disabled\n");
			} else if((s[0]=='h'&&(!s[1]||isspace(s[1])))
					||(!mystrcmp(s, "help")&&(!s[4]||isspace(s[4])))
					){
				print_help(argv[0]);
				return;
			} else{
				warning("Unknown options: %s\n", argv[i]);
			}
		} else {
			char *endl;
			const char *endargi=argv[i]+strlen(argv[i]);
			int temp=strtol(argv[i], &endl, 10);
			if(endl==endargi){//argument is port
				extern int PORT;
				PORT=temp;
			}else{
				parse_host(argv[i]);
			}
		}
	}
	register_signal_handler(NULL);
	if(socketpair(AF_UNIX, SOCK_STREAM, 0, sock_main)){
		error("failed to create socketpair\n");
	}
	socket_nopipe(sock_main[0]);
	socket_nopipe(sock_main[1]);
	thread_new(listen_host, GINT_TO_POINTER(sock_main[0]));//must be after parse_host
}
#if USE_APPLICATION
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
int main(int argc, char* argv[]){
	if(0){
		char* fnlog=stradd(TEMP, "/monitor.log", NULL);
		//info("Check %s for log.\n", fnlog);
		if(!freopen(fnlog, "w", stdout)){
			warning("Unable to redirect output to %s\n", fnlog);
		} else{
			setbuf(stdout, NULL);
		}
		free(fnlog);
	}
#if	!GLIB_VERSION_AFTER(2, 32)
	if(!g_thread_supported()){
		g_thread_init(NULL);
		gdk_threads_init();
	}
#endif
#if USE_APPLICATION==0
	if(!gtk_init_check(&argc, &argv)){//check whether we can enable window
		headless=1;
		plot_enabled=0;
		warning("Running in headless mode\n");
	}
#endif
#if WITH_NOTIFY
	if(!notify_init("AOS Notification")){
		notify_daemon=0;
	}
#endif
	parse_icons();

	if(!headless){
		create_status_icon();
#if MAC_INTEGRATION
		GtkosxApplication* theApp=g_object_new(GTKOSX_TYPE_APPLICATION, NULL);
		gtkosx_application_set_dock_icon_pixbuf(theApp, icon_main);
	//g_signal_connect(theApp,"NSApplicationDidBecomeActive", G_CALLBACK(status_icon_on_click), GINT_TO_POINTER(1));//useless
		gtkosx_application_ready(theApp);
#endif
	}
#if USE_APPLICATION==0
	handle_args(argc, argv);
	if(!headless){
		create_window();
		gtk_main();
	}else{
		pause();
	}
#else
	GtkApplication* app;
	int status;
	//application id must contain a .
	app=gtk_application_new("maos.monitor", G_APPLICATION_HANDLES_OPEN|G_APPLICATION_HANDLES_COMMAND_LINE|G_APPLICATION_NON_UNIQUE);
	g_signal_connect(app, "command-line", G_CALLBACK(command_line), NULL);
	g_signal_connect(app, "activate", G_CALLBACK(create_window), NULL);
	status=g_application_run(G_APPLICATION(app), argc, argv);
	g_object_unref(app);
	return status;
#endif
}

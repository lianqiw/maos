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
/*
  A monitor that monitors the progress of all running processes in all machines.

  The widgets duplicates and manages the string.
  Set NOTIFY_IS_SUPPORTED to 0 is libnotify is not installed.

  Todo:
  1) Log the activities into a human readable file
  2) Store the activity list into file when quit and reload later
*/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <glib.h>
#include <gtk/gtk.h>
#include <glib/gprintf.h>
#include <pthread.h>
#include <fcntl.h>
#include <unistd.h>
#ifndef GTK_WIDGET_VISIBLE
#define GTK_WIDGET_VISIBLE gtk_widget_get_visible
#endif
#if WITH_NOTIFY
#include <libnotify/notify.h>
static int notify_daemon=1;
#endif
#include "../sys/sys.h"
#include "monitor.h"
#include "icon-monitor.h"
#include "icon-finished.h"
#include "icon-running.h"
#include "icon-failed.h"
GdkPixbuf *icon_main=NULL;
GdkPixbuf *icon_finished=NULL;
GdkPixbuf *icon_failed=NULL;
GdkPixbuf *icon_running=NULL;

static GtkStatusIcon *status_icon;
GtkWidget *notebook=NULL;
GtkWidget **pages;
static GtkWidget *window=NULL;
static GtkWidget **tabs;
static GtkWidget **titles;
static GtkWidget **cmdconnect;
double *usage_cpu, *usage_cpu2;
double *usage_mem, *usage_mem2;
static GtkWidget **prog_cpu;
static GtkWidget **prog_mem;
PangoAttrList *pango_active, *pango_down;
GdkColor blue;
GdkColor green;
GdkColor red;
GdkColor yellow;
GdkColor white;
GdkColor color_even;
GdkColor color_odd;

GtkActionGroup *topgroup;
GtkWidget *toptoolbar;
#if GTK_MAJOR_VERSION>=3 
//GtkCssProvider *provider_prog;
GtkCssProvider *provider_red;
GtkCssProvider *provider_blue;
#endif
int *hsock;

/**
   The number line pattern determines how dash is drawn for gtktreeview. the
   first number is the length of the line, and second number of the length
   of blank after the line.

   properties belonging to certain widget that are descendent of GtkWidget, need
to specify the classname before the key like GtkTreeView::allow-rules */
#if GTK_MAJOR_VERSION<3
static const gchar *rc_string_widget =
    {
	
	"style \"widget\" {               \n"
	"font_name = \"Sans 8\""
     	"}\n"
	"class \"GtkWidget\" style \"widget\" \n"
	};
static const gchar *rc_string_treeview = 
    {
	"style \"solidTreeLines\"{                       \n" 
	/*" GtkTreeView::grid-line-pattern=\"\111\111\"\n" */
	" GtkTreeView::grid-line-width   = 1             \n"
	" GtkTreeView::allow-rules       = 1             \n"
	" GtkTreeView::odd-row-color     = \"#EFEFEF\"   \n"
	" GtkTreeView::even-row-color    = \"#FFFFFF\"   \n"
	" GtkTreeView::horizontal-separator = 0          \n"
	" GtkTreeView::vertical-separator = 0            \n"
	"}\n                                             \n" 
	"class \"GtkTreeView\" style \"solidTreeLines\"  \n" 
    };

static const gchar *rc_string_entry =
    {
	"style \"entry\" {               \n"
	"base[NORMAL] = \"white\"        \n"
	"bg[SELECTED] = \"#0099FF\"         \n"
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
    gtk_widget_set_sensitive(cmdconnect[ihost],0);
    gtk_widget_hide(cmdconnect[ihost]);
    gtk_label_set_attributes(GTK_LABEL(titles[ihost]), pango_active);
    return 0;
}
/**
   Disables the notebook page for a host*/
gboolean host_down(gpointer data){
    int ihost=GPOINTER_TO_INT(data);
    gtk_button_set_label(GTK_BUTTON(cmdconnect[ihost]),"Click to connect");
    gtk_widget_show_all(cmdconnect[ihost]);
    gtk_widget_set_sensitive(cmdconnect[ihost],1);
    gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(prog_cpu[ihost]), 0);
    gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(prog_mem[ihost]), 0);
    gtk_label_set_attributes(GTK_LABEL(titles[ihost]), pango_down);
    return 0;
}

/**
   modifies the color of progress bar*/
static void modify_bg(GtkWidget *widget, int type){
#if GTK_MAJOR_VERSION>=3 
    GtkCssProvider *provider;
    switch(type){
    case 1:
	provider=provider_blue;
	break;
    case 2:
	provider=provider_red;
	break;
    default:
	provider=NULL;
    }
    GtkStyleContext *context=gtk_widget_get_style_context(widget);
    gtk_style_context_add_provider(context, GTK_STYLE_PROVIDER(provider), 
				   GTK_STYLE_PROVIDER_PRIORITY_USER);
#else
    GdkColor *color;
    switch(type){
    case 1:
	color=&blue;
	break;
    case 2:
	color=&red;
	break;
    default:
	color=NULL;
    }
    gtk_widget_modify_bg(widget, GTK_STATE_SELECTED, color);
    gtk_widget_modify_bg(widget, GTK_STATE_PRELIGHT, color);
#endif
}
/**
   updates the progress bar for a job*/
gboolean update_progress(gpointer input){
    int ihost=GPOINTER_TO_INT(input);
    double last_cpu=usage_cpu2[ihost];
    double last_mem=usage_mem2[ihost];
    usage_cpu2[ihost]=usage_cpu[ihost];
    usage_mem2[ihost]=usage_mem[ihost];
    if(GTK_WIDGET_VISIBLE(window)){
	gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(prog_cpu[ihost]), usage_cpu[ihost]);
	gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(prog_mem[ihost]), usage_mem[ihost]);
	if(usage_cpu[ihost]>=0.8 && last_cpu<0.8){
	    modify_bg(prog_cpu[ihost],2);
	}else if(usage_cpu[ihost]<0.8 && last_cpu>=0.8){
	    modify_bg(prog_cpu[ihost],1);
	}
	if(usage_mem[ihost]>=0.8 && last_mem<0.8){
	    modify_bg(prog_mem[ihost],2);
	}else if(usage_mem[ihost]<0.8 && last_mem>=0.8){
	    modify_bg(prog_mem[ihost],1);
	}
    }
    return 0;
}

void notify_user(PROC_T *p){
    if(p->status.done || p->status.info==S_START) return;
#if WITH_NOTIFY
    if(!notify_daemon) return;
    if(p->status.info==p->oldinfo) return;
    static NotifyNotification *notify_urgent=NULL, *notify_normal=NULL, *notify_low=NULL;
    if (!notify_urgent){
#if defined(NOTIFY_CHECK_VERSION) 
#if NOTIFY_CHECK_VERSION(0,7,0) /*newer versions doesnot have _new_with_status_icon */
	notify_low=notify_notification_new("Low", NULL, NULL);
	notify_normal=notify_notification_new("Low", NULL, NULL);
	notify_urgent=notify_notification_new("Low", NULL, NULL);
#else
	notify_low=notify_notification_new_with_status_icon ("Low",NULL,NULL,status_icon);
	notify_normal=notify_notification_new_with_status_icon ("Normal",NULL,NULL,status_icon);
	notify_urgent=notify_notification_new_with_status_icon ("Urgent",NULL,"error",status_icon);
#endif
#else
	notify_low=notify_notification_new_with_status_icon ("Low",NULL,NULL,status_icon);
	notify_normal=notify_notification_new_with_status_icon ("Normal",NULL,NULL,status_icon);
	notify_urgent=notify_notification_new_with_status_icon ("Urgent",NULL,"error",status_icon);
#endif
	notify_notification_set_icon_from_pixbuf(notify_low,icon_main);
	notify_notification_set_timeout(notify_low,NOTIFY_EXPIRES_DEFAULT);
	notify_notification_set_urgency(notify_low,NOTIFY_URGENCY_LOW);

	notify_notification_set_icon_from_pixbuf(notify_normal,icon_main);
	notify_notification_set_timeout(notify_normal,NOTIFY_EXPIRES_DEFAULT);
	notify_notification_set_urgency(notify_normal,NOTIFY_URGENCY_NORMAL);

	notify_notification_set_timeout(notify_urgent,NOTIFY_EXPIRES_NEVER);
	notify_notification_set_urgency(notify_urgent,NOTIFY_URGENCY_CRITICAL);

    }
    static char summary[80];
    NotifyNotification *notify;
    switch(p->status.info){
    case S_START:
	notify=notify_low;
	snprintf(summary,80,"Job Started on %s",hosts[p->hid]);
	break;
    case S_FINISH:
	notify=notify_normal;
	snprintf(summary,80,"Job Finished on %s",hosts[p->hid]);
	break;
    case S_CRASH:
	notify=notify_urgent;
	snprintf(summary,80,"Job Crashed on %s!",hosts[p->hid]);
	break;
    case S_KILLED:
	notify=notify_urgent;
	snprintf(summary,80,"Job is Killed on %s!",hosts[p->hid]);
	break;
    default:
	warning("Invalid status\n");
	return;
    }
    notify_notification_update(notify,summary,p->path,NULL);
    notify_notification_show(notify,NULL);
    p->oldinfo=p->status.info;
#endif
}
/**
   quite the program
*/
static void quitmonitor(GtkWidget *widget, gpointer data){
    (void)widget;
    (void)data;
#if WITH_NOTIFY
    if(notify_daemon)
	notify_uninit();
#endif
    add_host_wrap(-2);
    if(gtk_main_level ()>0){
	gtk_main_quit();
    }
    exit(0);
}
/**
   update the job count
*/
gboolean update_title(gpointer data){
    int id=GPOINTER_TO_INT(data);
    char tit[40];
    snprintf(tit,40,"%s (%d)",hosts[id], nproc[id]);
    gtk_label_set_text(GTK_LABEL(titles[id]),tit);
    return 0;
}
/**
   respond to the kill job event
*/
static void kill_job(PROC_T *p){
    GtkWidget *dia=gtk_message_dialog_new
	(NULL, GTK_DIALOG_DESTROY_WITH_PARENT,
	 GTK_MESSAGE_QUESTION,
	 GTK_BUTTONS_NONE,
	 "Kill job %d on server %s?",
	 p->pid,hosts[p->hid]);
    if(p->status.info==S_WAIT){
	gtk_dialog_add_buttons(GTK_DIALOG(dia), "Kill Job", 0, "Cancel", 2, NULL);
    }else{
	gtk_dialog_add_buttons(GTK_DIALOG(dia), "Kill Job", 0, "Remote Display", 1, "Cancel", 2, NULL);
    }
    int result=gtk_dialog_run(GTK_DIALOG(dia));
    gtk_widget_destroy (dia);
    switch(result){
    case 0:
	if(scheduler_cmd(p->hid,p->pid,CMD_KILL)){
	    warning("Failed to kill the job\n");
	}
	p->status.info=S_TOKILL;
	refresh(p);
    	break;
    case 1:
	{
	    /*connect to scheduler with a new port. Pass the port to
	      drawdaemon.*/
	    int sock=connect_port(hosts[p->hid], PORT, 0, 0);
	    int cmd[2]={CMD_DISPLAY, p->pid};
	    if(stwriteintarr(sock, cmd, 2)){
		warning("Failed to communicate to scheduler\n");
		close(sock);
	    }
	    if(streadintarr(sock, cmd, 1)){
		warning("Read response failed\n");
	    }else if(cmd[0]){
		warning("The scheduler failed to talk to maos\n");
	    }else{
		if(spawn_drawdaemon(sock)){
		    warning("spwn drawdaemon failed\n");
		}
		close(sock);
	    }
	}
	break;
    }
}
void kill_job_event(GtkWidget *btn, GdkEventButton *event, PROC_T *p){
    (void)btn;
    if(event->button==1){
	kill_job(p);
    }
}
static void kill_all_jobs(GtkAction *btn){
    (void)btn;
    int ihost=gtk_notebook_get_current_page (GTK_NOTEBOOK(notebook));
    GtkWidget *dia=gtk_message_dialog_new
	(NULL, GTK_DIALOG_DESTROY_WITH_PARENT,
	 GTK_MESSAGE_QUESTION,
	 GTK_BUTTONS_YES_NO,
	 "Kill all jobs on server %s?", hosts[ihost]);
    int result=gtk_dialog_run(GTK_DIALOG(dia));
    gtk_widget_destroy (dia);
    if(result==GTK_RESPONSE_YES){
	for(PROC_T *iproc=pproc[ihost]; iproc; iproc=iproc->next){
	    if(iproc->hid == ihost 
	       && (iproc->status.info==S_WAIT 
		   || iproc->status.info==S_RUNNING
		   || iproc->status.info==S_START)){
		if(scheduler_cmd(iproc->hid,iproc->pid,CMD_KILL)){
		    warning("Failed to kill the job\n");
		}
		iproc->status.info=S_TOKILL;
		refresh(iproc);
	    }
	}
    }
}
static void status_icon_on_click(GtkStatusIcon *status_icon0, 
			       gpointer data){
    (void)status_icon0;
    (void)data;
    static int cx=0, cy=0;
    static int x, y;
    if(GTK_WIDGET_VISIBLE(window)){
	gtk_window_get_size(GTK_WINDOW(window), &x, &y);
	gtk_window_get_position(GTK_WINDOW(window), &cx, &cy);
	gtk_widget_hide(window);
    }else{
	gtk_window_set_default_size(GTK_WINDOW(window), x, y);
	gtk_window_move(GTK_WINDOW(window), cx, cy);
	gtk_widget_show(window);
	gtk_window_deiconify(GTK_WINDOW(window));
	gtk_window_stick(GTK_WINDOW(window));
    }
}
static void trayIconPopup(GtkStatusIcon *status_icon0, guint button, 
			  guint32 activate_time, gpointer popUpMenu){
    gtk_menu_popup(GTK_MENU(popUpMenu),NULL,NULL,
		   gtk_status_icon_position_menu,
		   status_icon0,button,activate_time);
}
static void create_status_icon(){
    status_icon = gtk_status_icon_new_from_pixbuf(icon_main);
    g_signal_connect(G_OBJECT(status_icon),"activate",
		     G_CALLBACK(status_icon_on_click), NULL);
    GtkWidget *menu, *menuItemShow,*menuItemExit;
    menu=gtk_menu_new();
    menuItemShow=gtk_menu_item_new_with_label("Show/Hide");
    menuItemExit=gtk_image_menu_item_new_from_stock(GTK_STOCK_QUIT,NULL);
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
    gtk_status_icon_set_visible(status_icon, TRUE);
}


static gboolean delete_window(void){
    if(gtk_status_icon_is_embedded(status_icon)){
	gtk_widget_hide(window);
	return TRUE;/*do not quit */
    }else{
	return FALSE;/*quit. */
    }
}

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
}
static int test_skipped(int status, double frac){
    return status==S_FINISH && frac<1;
}
static int test_finished(int status, double frac){
    (void)frac;
    return status==S_FINISH;
}
static int test_crashed(int status, double frac){
    (void)frac;
    return status==S_CRASH || status==S_KILLED || status==S_TOKILL;
}
static void clear_jobs(GtkAction *btn, int (*testfun)(int, double)){
    (void)btn;
    int ihost=gtk_notebook_get_current_page (GTK_NOTEBOOK(notebook));
    if(!pproc[ihost]) return;
    int sock=hsock[ihost];
    if(sock==-1) return;
    int cmd[2];
    cmd[0]=CMD_REMOVE;
    for(PROC_T *iproc=pproc[ihost]; iproc; iproc=iproc->next){
	if(testfun(iproc->status.info, iproc->frac)){
	    cmd[1]=iproc->pid;
	    if(stwriteintarr(sock,cmd,2)){
		warning("write to socket %d failed\n",sock);
		break;
	    }
	}
    }
}
static void clear_jobs_finished(GtkAction *btn){
    clear_jobs(btn, test_finished);
}
static void clear_jobs_crashed(GtkAction *btn){
    clear_jobs(btn, test_crashed);
}
static void clear_jobs_skipped(GtkAction *btn){
    clear_jobs(btn, test_skipped);
}
GtkWidget *monitor_new_entry_progress(void){
    GtkWidget *prog=gtk_entry_new();
    gtk_editable_set_editable(GTK_EDITABLE(prog),FALSE);
    gtk_entry_set_width_chars(GTK_ENTRY(prog),12);
#if GTK_MAJOR_VERSION>=3 || GTK_MINOR_VERSION >= 18
    gtk_widget_set_can_focus(prog,FALSE);
#else
    g_object_set(prog,"can-focus",FALSE,NULL);
#endif
    
    /*gtk_widget_modify_base(prog,GTK_STATE_NORMAL, &white);
    gtk_widget_modify_bg(prog,GTK_STATE_SELECTED, &blue);
    */
#if GTK_MAJOR_VERSION==2 && GTK_MINOR_VERSION >= 10 || GTK_MAJOR_VERSION==3 && GTK_MINOR_VERSION < 4
    gtk_entry_set_inner_border(GTK_ENTRY(prog), 0);
#endif
    gtk_entry_set_has_frame (GTK_ENTRY(prog), 0);
    gtk_entry_set_alignment(GTK_ENTRY(prog),0.5);
    return prog;
}
GtkWidget *monitor_new_progress(int vertical, int length){
    GtkWidget *prog=gtk_progress_bar_new();
    if(vertical){
#if GTK_MAJOR_VERSION>=3 
	gtk_orientable_set_orientation(GTK_ORIENTABLE(prog),
				       GTK_ORIENTATION_VERTICAL);
	gtk_progress_bar_set_inverted(GTK_PROGRESS_BAR(prog), TRUE);

#else
	gtk_progress_bar_set_orientation(GTK_PROGRESS_BAR(prog),
					 GTK_PROGRESS_BOTTOM_TO_TOP);
#endif
	gtk_widget_set_size_request(prog, 8,length);
	g_object_set(G_OBJECT(prog), "show-text", FALSE, NULL);
    }else{
	gtk_widget_set_size_request(prog, length,12);
    }
    return prog;
}

static void add_host_event(GtkButton *button, gpointer data){
    (void)button;
    int ihost=GPOINTER_TO_INT(data);
    if(ihost>-1 && ihost<nhost){
	add_host_wrap(ihost);
    }else{
	for(ihost=0; ihost<nhost; ihost++){
	    add_host_wrap(ihost);
	}
    }
}

int main(int argc, char *argv[])
{
#if GLIB_MAJOR_VERSION<3 && GLIB_MINOR_VERSION<32
    if(!g_thread_supported()){
	g_thread_init(NULL);
	gdk_threads_init();
    }
#endif
    gtk_init(&argc, &argv);
#if WITH_NOTIFY
    if(!notify_init("AOS Notification")){
	notify_daemon=0;
    }
#endif
    if(argc>1){
	hosts=realloc(hosts, sizeof(char*)*(nhost+(argc-1)));
	for(int i=1; i<argc; i++){
	    hosts[nhost]=strdup(argv[i]);
	    nhost++;
	}
    }else if(nhost==1){
	info2("Using '%s host1 host2' to monitor other machines, or put their hostnames in ~/.aos/hosts\n", argv[0]);
    }
    gdk_color_parse("#EE0000",&red);
    gdk_color_parse("#00CC00",&green);
    gdk_color_parse("#FFFF33",&yellow);
    gdk_color_parse("#FFFFFF",&white);
    gdk_color_parse("#0099FF",&blue);
    gdk_color_parse("#FFFFAB",&color_even);
    gdk_color_parse("#FFFFFF",&color_odd);


    icon_main=gdk_pixbuf_new_from_inline(-1,icon_monitor, FALSE, NULL);
    icon_finished=gdk_pixbuf_new_from_inline(-1,icon_inline_finished,FALSE,NULL);
    icon_failed=gdk_pixbuf_new_from_inline(-1,icon_inline_failed,FALSE,NULL);
    icon_running=gdk_pixbuf_new_from_inline(-1,icon_inline_running,FALSE,NULL);
    /*GtkWidget *image_finished=gtk_image_new_from_stock(GTK_STOCK_APPLY, GTK_ICON_SIZE_MENU);
    gtk_widget_show(image_finished);
    icon_finished=gtk_image_get_pixbuf(GTK_IMAGE(image_finished));
    */
    create_status_icon();
    window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(window),"MAOS Monitor");

    gtk_window_set_icon(GTK_WINDOW(window),icon_main);


#if GTK_MAJOR_VERSION<3
    gtk_rc_parse_string(rc_string_widget); 
    gtk_rc_parse_string(rc_string_treeview);
    gtk_rc_parse_string(rc_string_entry);
    //    GtkStyle *style=gtk_widget_get_style(window);
#else
    /*trough is the main background. progressbar is the sizable bar.*/

    const gchar *prog_blue=".progressbar{"
	"background-image:-gtk-gradient(linear,left bottom, right bottom, from (#0000FF), to (#0000FF));\n}" ;
    const gchar *prog_red=".progressbar{"
	"background-image:-gtk-gradient(linear,left bottom, right bottom, from (#FF0000), to (#FF0000));\n}" ;

    //    provider_prog=gtk_css_provider_new();
    //gtk_css_provider_load_from_data(provider_prog, prog_style, strlen(prog_style), NULL);
    provider_blue=gtk_css_provider_new();
    gtk_css_provider_load_from_data(provider_blue, prog_blue, strlen(prog_blue), NULL);
    provider_red=gtk_css_provider_new();
    gtk_css_provider_load_from_data(provider_red, prog_red, strlen(prog_red), NULL);
    /*Properties not belonging to GtkWidget need to begin with -WidgetClassName*/
    const gchar *all_style="*{"
	"padding:0;\n"
	"border-radius:0;\n"
	"border-width:1;\n"
	"font:Sans 8;\n"
	"-GtkToolbar-button-releaf:0\n"
	"}"
	".notebook{"
	"-GtkNotebook-tab-overlap: 2;"
	"-GtkNotebook-tab-curvature: 2;"
	"-GtkWidget-focus-line-width: 0;"
	"background-color:@theme_bg_color;"
        "base:@theme_bg_color;"
	"}"
	".notebook tab{"
	"-adwaita-focus-border-radius: 0;"
	"border-width: 0;"
	"padding: 0 0 0;"
	"background-color:#CCCCCC;"
	"}"
	".notebook tab:active{"
	"-adwaita-focus-border-radius: 1;"
	"border-width: 2;"
	"padding: 0 0 0;"
        "background-color:@theme_bg_color;"
	"}"
	".entry.progressbar{"
	"-GtkEntry-has-frame:1;\n"
	"-GtkEntry-progress-border:0,0,0,0;\n"
	"-GtkEntry-inner-border:0,0,0,0;\n"
	"background-image:-gtk-gradient(linear,left bottom, right bottom, from (#0000FF), to (#0000FF));\n"
	"border-width:1; border-style:none; \n"
	"border-color:#000000;\n"
	"}"
	"GtkProgressBar {\n"
	"-GtkProgressBar-min-horizontal-bar-height:6;\n"
	"-GtkProgressBar-min-horizontal-bar-width:6;\n"
	"-GtkProgressBar-min-vertical-bar-height:6;\n"
	"-GtkProgressBar-min-vertical-bar-width:6;\n"
	"-GtkProgressBar-xspacing:0;\n"
	"-GtkProgressBar-yspacing:0;\n"
	"padding:1;\n" //this is the gap betwen trough and bar.
	"}"
	"GtkProgressBar.progressbar{"
	"background-image:-gtk-gradient(linear,left bottom, right bottom, from (#0000FF), to (#0000FF));\n" //this is the bar
	"}"
	"GtkProgressBar.progressbar,.progressbar{"
	"border-width:0;"
	"padding:0;"
	"border-radius:0;" //make the bar square
	"-adwaita-progressbar-pattern: none;"
	"}"
	".progressbar.vertical{"
	"-adwaita-progressbar-pattern: none;"
	"border-radius:0;"
	"padding:0;"
	"}"
	".trough, .trough row, .trough row:hover{"
	"border-width: 1;\n"
	"-GtkProgressBar-xspacing: 0;"
	"-GtkProgressBar-yspacing: 0;"
	"border-radius:0;\n" //make the trough square
	//"background-image:-gtk-gradient(linear, left bottom, left top, from (#0000FF), to (#FF0000);\n"
	"}";
    GtkCssProvider *provider_default=gtk_css_provider_new();
    gtk_css_provider_load_from_data(provider_default, all_style, strlen(all_style), NULL);
    GtkStyleContext *all_context=gtk_widget_get_style_context(window);
    GdkScreen *screen=gtk_style_context_get_screen(all_context);
    gtk_style_context_add_provider_for_screen(screen, GTK_STYLE_PROVIDER(provider_default),
					      GTK_STYLE_PROVIDER_PRIORITY_USER);
    
#endif


    GtkWidget *vbox=gtk_vbox_new(FALSE, 0);
    gtk_container_add(GTK_CONTAINER(window), vbox);
    {
	/*first create actions*/
	topgroup=gtk_action_group_new("topgroup");
	GtkAction *action_clear_finish=gtk_action_new("act-clear-finish", "Clear finished", "Clear finished jobs", GTK_STOCK_APPLY);
	g_signal_connect(action_clear_finish, "activate", G_CALLBACK(clear_jobs_finished),NULL);
	gtk_action_group_add_action(topgroup, action_clear_finish);

	GtkAction *action_clear_skip=gtk_action_new("act-clear-skip", "Clear skipped", "Clear skipped jobs", GTK_STOCK_MEDIA_NEXT);
	g_signal_connect(action_clear_skip, "activate", G_CALLBACK(clear_jobs_skipped),NULL);
	gtk_action_group_add_action(topgroup, action_clear_skip);

	GtkAction *action_clear_crash=gtk_action_new("act-clear-crash", "Clear crashed", "Clear crashed jobs", GTK_STOCK_DELETE);
	g_signal_connect(action_clear_crash, "activate", G_CALLBACK(clear_jobs_crashed),NULL);
	gtk_action_group_add_action(topgroup, action_clear_crash);

	GtkAction *action_connect_all=gtk_action_new("act-connect-all", "Connect all", "Connect to all hosts", GTK_STOCK_CONNECT);
	g_signal_connect(action_connect_all, "activate", G_CALLBACK(add_host_event),GINT_TO_POINTER(-1));
	gtk_action_group_add_action(topgroup, action_connect_all);

	GtkAction *action_kill_all=gtk_action_new("act-kill-all", "Kill all jobs", "Kill al jobs", GTK_STOCK_CANCEL);
	g_signal_connect(action_kill_all, "activate", G_CALLBACK(kill_all_jobs),NULL);
	gtk_action_group_add_action(topgroup, action_kill_all);
	
	/*set toolbar*/
	toptoolbar=gtk_toolbar_new();
	gtk_toolbar_insert(GTK_TOOLBAR(toptoolbar), GTK_TOOL_ITEM(gtk_action_create_tool_item(action_clear_finish)), -1);
	gtk_toolbar_insert(GTK_TOOLBAR(toptoolbar), GTK_TOOL_ITEM(gtk_action_create_tool_item(action_clear_skip)), -1);
	gtk_toolbar_insert(GTK_TOOLBAR(toptoolbar), GTK_TOOL_ITEM(gtk_action_create_tool_item(action_clear_crash)), -1);
	gtk_toolbar_insert(GTK_TOOLBAR(toptoolbar), GTK_TOOL_ITEM(gtk_action_create_tool_item(action_connect_all)), -1);
	gtk_toolbar_insert(GTK_TOOLBAR(toptoolbar), GTK_TOOL_ITEM(gtk_action_create_tool_item(action_kill_all)), -1);
	//gtk_toolbar_insert(GTK_TOOLBAR(toptoolbar), gtk_separator_tool_item_new(), -1);
	gtk_widget_show_all(toptoolbar);
	gtk_box_pack_start(GTK_BOX(vbox), toptoolbar, FALSE, FALSE, 0);
	gtk_toolbar_set_icon_size(GTK_TOOLBAR(toptoolbar), GTK_ICON_SIZE_MENU);
    }

    notebook=gtk_notebook_new();
    gtk_widget_show(notebook);
    gtk_box_pack_start(GTK_BOX(vbox), notebook, TRUE, TRUE, 0);       

    gtk_notebook_set_scrollable(GTK_NOTEBOOK(notebook), TRUE);
    gtk_notebook_set_tab_pos(GTK_NOTEBOOK(notebook),GTK_POS_LEFT);


    g_signal_connect(window, "delete_event", G_CALLBACK (delete_window), NULL);
    g_signal_connect(window, "destroy", G_CALLBACK (quitmonitor), NULL);
    g_signal_connect(G_OBJECT (window), "window-state-event", G_CALLBACK (window_state_event), NULL);
    gtk_window_set_position(GTK_WINDOW(window),GTK_WIN_POS_CENTER);
    gtk_window_set_default_size(GTK_WINDOW(window), 840, 400);
    gtk_widget_show_all(window);

    tabs=calloc(nhost,sizeof(GtkWidget*));
    pages=calloc(nhost,sizeof(GtkWidget*));
    titles=calloc(nhost,sizeof(GtkWidget*));
    pproc=calloc(nhost,sizeof(PROC_T*));
    nproc=calloc(nhost,sizeof(int));
    cmdconnect=calloc(nhost,sizeof(GtkWidget*));
    hsock=calloc(nhost, sizeof(int));
    usage_cpu=calloc(nhost,sizeof(double));
    usage_mem=calloc(nhost,sizeof(double));
    usage_cpu2=calloc(nhost,sizeof(double));
    usage_mem2=calloc(nhost,sizeof(double));
    prog_cpu=calloc(nhost, sizeof(GtkWidget *));
    prog_mem=calloc(nhost, sizeof(GtkWidget *));

    pango_active=pango_attr_list_new();
    pango_down=pango_attr_list_new();
    pango_attr_list_insert(pango_down, pango_attr_foreground_new(0x88FF, 0x88FF, 0x88FF));
    pango_attr_list_insert(pango_active, pango_attr_foreground_new(0x0000, 0x0000, 0x0000));

    for(int ihost=0; ihost<nhost; ihost++){
	char tit[40];
	snprintf(tit,40,"%s(0)",hosts[ihost]);
    	titles[ihost]=gtk_label_new(tit);
	gtk_label_set_attributes(GTK_LABEL(titles[ihost]), pango_down);
	GtkWidget *hbox0=gtk_hbox_new(FALSE,0);
	prog_cpu[ihost]=monitor_new_progress(1,16);
	prog_mem[ihost]=monitor_new_progress(1,16);
	gtk_box_pack_start(GTK_BOX(hbox0),prog_cpu[ihost], FALSE, FALSE, 1);
	gtk_box_pack_start(GTK_BOX(hbox0),prog_mem[ihost], FALSE, FALSE, 1);
	modify_bg(prog_cpu[ihost], 1);
	modify_bg(prog_mem[ihost], 1);
	gtk_box_pack_start(GTK_BOX(hbox0),titles[ihost], FALSE, TRUE, 0);
	gtk_widget_show_all(hbox0);
	GtkWidget *eventbox=gtk_event_box_new();
	gtk_container_add(GTK_CONTAINER(eventbox),hbox0);
	/*Put the box below an eventbox so that clicking on the progressbar
	  switches tab also.*/
	gtk_widget_show_all(eventbox);
	gtk_event_box_set_above_child(GTK_EVENT_BOX(eventbox),TRUE);
	gtk_event_box_set_visible_window(GTK_EVENT_BOX(eventbox),FALSE);

	cmdconnect[ihost]=gtk_button_new_with_label("Click to connect");
	g_signal_connect(cmdconnect[ihost],"clicked", G_CALLBACK(add_host_event), GINT_TO_POINTER(ihost));
	GtkWidget *hbox=gtk_hbox_new(FALSE,0);
	gtk_box_pack_start(GTK_BOX(hbox),cmdconnect[ihost],TRUE,FALSE,0);

	pages[ihost]=gtk_vbox_new(FALSE,0);
	gtk_box_pack_end(GTK_BOX(pages[ihost]),hbox,FALSE,FALSE,0);
	GtkWidget *page=new_page(ihost);
	if(page){
	    gtk_box_pack_start(GTK_BOX(pages[ihost]),page,FALSE,FALSE,0);
	}
	tabs[ihost]=gtk_scrolled_window_new(NULL,NULL);
	gtk_scrolled_window_set_policy
	    (GTK_SCROLLED_WINDOW(tabs[ihost]),
	     GTK_POLICY_AUTOMATIC,GTK_POLICY_AUTOMATIC);
	gtk_scrolled_window_add_with_viewport
	    (GTK_SCROLLED_WINDOW(tabs[ihost]),pages[ihost]);
	gtk_notebook_append_page
	    (GTK_NOTEBOOK(notebook),tabs[ihost],eventbox);
	update_title(GINT_TO_POINTER(ihost));
	gtk_widget_show_all(tabs[ihost]);
	gtk_notebook_set_current_page(GTK_NOTEBOOK(notebook), ihost);
    }
    extern int pipe_main[2];
    if(pipe(pipe_main)<0){
	error("failed to create pipe\n");
    }
    pthread_t tmp;
    pthread_create(&tmp, NULL, (void*(*)(void*))listen_host, NULL);
    for(int ihost=0; ihost<nhost; ihost++){
	add_host_wrap(ihost);
    }
    gtk_main();
}

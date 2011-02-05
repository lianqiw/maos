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
#include <netdb.h>
#include <unistd.h>
#include <sys/types.h>
#include <fcntl.h> 
#include <errno.h>
#include <arpa/inet.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <string.h>
#include <glib.h>
#include <gtk/gtk.h>
#include <glib/gprintf.h>
#include <pthread.h>
#if WITH_NOTIFY
#include <libnotify/notify.h>
static int notify_daemon=1;
#endif
#include "common.h"
#include "misc.h"
#include "daemonize.h"
#include "scheduler_client.h"
#include "io.h"
#include "monitor.h"
#include "icon-monitor.h"
#include "icon-finished.h"
#include "icon-running.h"
#include "icon-failed.h"
GdkPixbuf *icon_main=NULL;
GdkPixbuf *icon_finished=NULL;
GdkPixbuf *icon_failed=NULL;
GdkPixbuf *icon_running=NULL;

static const char *ProgName="Job Monitoring";
static GtkStatusIcon *status_icon;
GtkWidget *notebook=NULL;
GtkWidget **pages;
static GtkWidget *window=NULL;
static GtkWidget **tabs;
static GtkWidget **titles;
static GtkWidget **cmdconnect;
static int *hsock;
static double *usage_cpu;
static double *usage_mem;
static GtkWidget **prog_cpu;
static GtkWidget **prog_mem;
GdkColor blue;
GdkColor green;
GdkColor red;
GdkColor yellow;
GdkColor white;
GdkColor *bg;
GdkColor color_even;
GdkColor color_odd;
static int scheduler_crashed;

/**
   The number line pattern determines how dash is drawn for gtktreeview. the
   first number is the length of the line, and second number of the length
   of blank after the line.
*/
static const gchar *rc_string_treeview = 
    {
	"style \"solidTreeLines\"{                       \n" 
	" GtkTreeView::allow-rules       = 1             \n"
	" GtkTreeView::grid-line-pattern = \"\255\1\"    \n" 
	" GtkTreeView::grid-line-width   = 1             \n"
	" GtkTreeView::odd-row-color     = \"#EFEFFF\"   \n"
	" GtkTreeView::even-row-color    = \"#FFFFFF\"   \n"
	" GtkTreeView::horizontal-separator = 0          \n"
	" GtkTreeView::vertical-separator = 0            \n"
	"}\n                                             \n" 
	"class \"GtkTreeView\" style \"solidTreeLines\"  \n" 
    };
static const gchar *rc_string_progress = 
    {
	"style \"progressbar\" {                       \n"
	"fg[SELECTED] = \"red\"                            \n"
	"bg[SELECTED] = \"yellow\"                         \n"
	"base[SELECTED] = \"green\"                        \n"
	"text[SELECTED] = \"purple\"                       \n"
	"}\n                                           \n"
	"class \"GtkCellRenderer\" style \"progressbar\" \n"
    };

PROC_T **pproc;
int *nproc;
static int nhostup=0;
static int quitall=0;
pthread_cond_t pcond;
pthread_mutex_t pmutex;
static PROC_T *proc_get(int id,int pid);
static PROC_T *proc_add(int id,int pid);
static void add_host_wakeup(void);
static int host_from_sock(int sock){
    for(int ihost=0; ihost<nhost; ihost++){
	if(hsock[ihost]==sock){
	    return ihost;
	}
    }
    return -1;
}
static void host_up(int host){
    pthread_mutex_lock(&pmutex);
    nhostup++;
    pthread_mutex_unlock(&pmutex);
    gtk_widget_set_sensitive(cmdconnect[host],0);
    gtk_widget_hide(cmdconnect[host]);
}
static void host_down(int host, int info){
    pthread_mutex_lock(&pmutex);
    nhostup--;
    pthread_mutex_unlock(&pmutex);
    add_host_wakeup();
    proc_remove_all(host);//remove all entries.
    static const char *infotext[]={"Connection is lost. Click to reconnect.",
				   "Scheduler version is too old",
				   "Scheduler verison is too new, plase update monitor"};

    gtk_button_set_label(GTK_BUTTON(cmdconnect[host]),infotext[info]);
    gtk_widget_show_all(cmdconnect[host]);
    gtk_widget_set_sensitive(cmdconnect[host],1);
    gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(prog_cpu[host]), 0);
}

static void channel_removed(gpointer data){
    /*
      The socket seems to be already shutdown and closed, so don't do it again.
    */
    int sock=GPOINTER_TO_INT(data);
    int ihost=host_from_sock(sock);
    if(ihost!=-1){
	host_down(ihost,0);
	hsock[ihost]=0;
    }
}
static gboolean respond(GIOChannel *source, GIOCondition cond, gpointer data){
    /*
      Return false causes GIOChannel to be removed from watch list.
    */
    int sock=GPOINTER_TO_INT(data);
    gsize nread;
    int cmd[3];
    if(cond&G_IO_HUP || cond&G_IO_ERR || cond&G_IO_NVAL){
	warning2("Lost connection to %s\n", hosts[host_from_sock(sock)]);
	return FALSE;
    }
    GIOStatus status;
    
    /*
      g_io_channel_read_to_end only returns at EOF. which
      couldn't happen unless socket closes g_io_channel_read_line
      returns at newline \n.
    */

    status=g_io_channel_read_chars(source,(gchar*)cmd,
				   3*sizeof(int),&nread,NULL);
	
    if(status==G_IO_STATUS_EOF){
	return FALSE;
    }
    if(status==G_IO_STATUS_EOF|| status==G_IO_STATUS_ERROR
       || status==G_IO_STATUS_AGAIN||nread!=3*sizeof(int)){
	warning("Error encountered. Disconnect\n");
	return FALSE;//disconnect.
    }
    int host=-1;
    for(int ihost=0; ihost<nhost; ihost++){
	if(hsock[ihost]==sock){
	    host=ihost;
	    break;
	}
    }
    if(host<0){
	warning("sock not found\n");
	return FALSE;
    }

    int pid=cmd[2];
    switch(cmd[0]){
    case CMD_VERSION:
	break;
    case CMD_STATUS:
	{
	    PROC_T *p=proc_get(host,pid);
	    if(!p){
		p=proc_add(host,pid);
	    }
	    if(g_io_channel_read_chars(source,(gchar*)&p->status,sizeof(STATUS_T),&nread,NULL)
	       !=G_IO_STATUS_NORMAL){
		warning("Error reading status\n");
		return FALSE;
	    }
	    refresh(p);
	}
	break;
    case CMD_PATH:
	{
	    PROC_T *p=proc_get(host,pid);
	    if(!p){
		p=proc_add(host,pid);
	    }
	    p->path=readstr(sock);
	}
	break;
    case CMD_LOAD:
	{
	    double last_cpu=usage_cpu[host];
	    double last_mem=usage_mem[host];
	    usage_cpu[host]=(double)((pid>>16) & 0xFFFF)/100.;
	    usage_mem[host]=(double)(pid & 0xFFFF)/100.;

	    if(GTK_WIDGET_VISIBLE(window)){
		gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(prog_cpu[host]), usage_cpu[host]);
		gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(prog_mem[host]), usage_mem[host]);
		if(usage_cpu[host]>=0.8 && last_cpu<0.8){
		    gtk_widget_modify_bg(prog_cpu[host],GTK_STATE_SELECTED,&red);
		    //somehow the XFCE engines only respect PRELIGHT
		    gtk_widget_modify_bg(prog_cpu[host],GTK_STATE_PRELIGHT,&red);
		}else if(usage_cpu[host]<0.8 && last_cpu>=0.8){
		    gtk_widget_modify_bg(prog_cpu[host],GTK_STATE_SELECTED,&blue);
		    gtk_widget_modify_bg(prog_cpu[host],GTK_STATE_PRELIGHT,&blue);
		}
		if(usage_mem[host]>=0.8 && last_mem<0.8){
		    gtk_widget_modify_bg(prog_mem[host],GTK_STATE_SELECTED,&red);
		    gtk_widget_modify_bg(prog_mem[host],GTK_STATE_PRELIGHT,&red);
		}else if(usage_mem[host]<0.8 && last_mem>=0.8){
		    gtk_widget_modify_bg(prog_mem[host],GTK_STATE_SELECTED,&blue);
		    gtk_widget_modify_bg(prog_mem[host],GTK_STATE_PRELIGHT,&blue);
		}
	    }
	}
	break;
    default:
	warning3("Invalid cmd %d\n",cmd[0]);
	return FALSE;
    }
    return TRUE;
}

int init_sockaddr (struct sockaddr_in *name,
                    const char *hostname, uint16_t port){
    struct hostent *hostinfo;
    
    name->sin_family = AF_INET;
    name->sin_port = htons(port);
    hostinfo = gethostbyname (hostname);
    if (hostinfo == NULL){
	perror("gethostbyname");
	fprintf (stderr, "Unknown host %s.\n", hostname);
	return -1;
    }else{
	name->sin_addr = *(struct in_addr *) hostinfo->h_addr;
	return 0;
    }
}
/**
   To open a port and connect to scheduler
*/
int scheduler_connect(int ihost, int block, int mode){
    /*
      mode=0: read/write
      mode=1: read only by the client. the server won't read
      mode=2: write only by the client. the server won't write
     */
    const char *host;
    if(ihost==-1){
	ihost=hid;
    }
    host=hosts[ihost];
    int sock;
 
    if(scheduler_crashed) {
	return -1;
    }

    struct sockaddr_in servername;
 
    /* Create the socket. */
    sock = socket (PF_INET, SOCK_STREAM, 0);
    socket_tcp_keepalive(sock);
    if (sock < 0) {
	perror ("socket (scheduler)");
	scheduler_crashed=1; 
	return sock;
    }

    int oldflag=fcntl(sock,F_GETFD,0);
    oldflag |= FD_CLOEXEC;
    fcntl(sock, F_SETFD, oldflag);//close on exec.
    
    if(init_sockaddr (&servername, host, PORT)){
	warning3("Unable to init_sockaddr.");
	scheduler_crashed=1; 
	close(sock);
	return -1;
    }
    int count=0;
    while(connect(sock, (struct sockaddr *)&servername, sizeof (servername))<0){
	if(!block){
	    close(sock);
	    return -1;
	}
	sleep(4);
	count++;
	if(count>1000){
	    close(sock);
	    sock=-1;
	    error("Failed to connect to scheduer\n");
	}
    }
    scheduler_shutdown(&sock,mode);
    return sock;
}

/**
   called by monitor to kill a process
*/
int scheduler_kill_job(int host,int pid){
    int sock=scheduler_connect(host,0,2);
    if(sock==-1) return 1;
    int cmd[2];
    cmd[0]=CMD_KILL;
    cmd[1]=pid;//pid
    if(write(sock,cmd,2*sizeof(int))!=sizeof(int)*2){
	warning("write to socket %d failed\n",sock);
    }
    close(sock);
    return 0;//success
}

/**
   Ask scheduler to remove run by monitor
*/
int scheduler_remove_job(int host, int pid){
    int sock=scheduler_connect(host,0,2);
    if(sock==-1) return 1;
    int cmd[2];
    cmd[0]=CMD_REMOVE;
    cmd[1]=pid;//pid
    if(write(sock,cmd,2*sizeof(int))!=sizeof(int)*2){
	warning("write to socket %d failed\n",sock);
    }
    close(sock);
    return 0;//success
}

/**
   This thread lives in a separate thread, so need to use gdk_thread_enter
   before calls gtk functions.
*/
static void add_host_thread(void){
    /*From GTK Manual: g_io_channel_get_flags (): Gets the current flags for a
      GIOChannel, including read-only flags such as G_IO_FLAG_IS_READABLE.  The
      values of the flags G_IO_FLAG_IS_READABLE and G_IO_FLAG_IS_WRITEABLE are
      cached for internal use by the channel when it is created. If they should
      change at some later point (e.g. partial shutdown of a socket with the
      UNIX shutdown() function), the user should immediately call
      g_io_channel_get_flags() to update the internal values of these flags.*/
    while(!quitall){
	for(int ihost=0; ihost<nhost; ihost++){
	    if(!hsock[ihost]){
		int sock=scheduler_connect(ihost,0,0);
		if(sock==-1){
		    hsock[ihost]=0;
		    warning2("Unable to connect to %s\n",hosts[ihost]);
		}else{
		    int cmd[2];
		    cmd[0]=CMD_MONITOR;
		    cmd[1]=scheduler_version;
		    if(write(sock,cmd,sizeof(int)*2)!=sizeof(int)*2){
			hsock[ihost]=0;
		    }else{
			//2010-07-03:we don't write. to detect remote close.
			shutdown(sock,SHUT_WR);
			gdk_threads_enter();
			GIOChannel *channel=g_io_channel_unix_new(sock);
			g_io_channel_set_encoding(channel,NULL,NULL);
			//must not be buffered
			g_io_channel_set_buffered(channel,FALSE);
			g_io_channel_set_close_on_unref (channel,1);
			g_io_add_watch_full
			    (channel,0, (GIOCondition)
			     (G_IO_IN|G_IO_HUP|G_IO_ERR|G_IO_PRI|G_IO_NVAL), 
			     respond,GINT_TO_POINTER(sock),channel_removed);
		    
			g_io_channel_unref(channel);
			gdk_threads_leave();
			hsock[ihost]=sock;
		    }
		}
		gdk_threads_enter();
		if(hsock[ihost]){
		    host_up(ihost);
		}
		gdk_threads_leave();
	    }
	}
	pthread_mutex_lock(&pmutex);
	if(nhostup==nhost){//all host is up
	    pthread_cond_wait(&pcond, &pmutex);
	}else{//sleep 5 seconds before retry.
	    struct timespec abstime;
	    abstime.tv_sec=myclockd()+5;
	    abstime.tv_nsec=0;
	    pthread_cond_timedwait(&pcond, &pmutex, &abstime);
	}
	pthread_mutex_unlock(&pmutex);
    }
}
static void add_host_wakeup(void){
    //Wakeup add_host_thread
    pthread_cond_signal(&pcond);
}
void notify_user(PROC_T *p){
    if(p->status.done) return;
#if WITH_NOTIFY
    if(!notify_daemon) return;
    if(p->status.info==p->oldinfo) return;
    static NotifyNotification *notify_urgent=NULL, *notify_normal=NULL, *notify_low=NULL;
    if (!notify_urgent){
	notify_low=notify_notification_new_with_status_icon ("Low",NULL,NULL,status_icon);
	notify_notification_set_icon_from_pixbuf(notify_low,icon_main);
	notify_notification_set_timeout(notify_low,NOTIFY_EXPIRES_DEFAULT);
	notify_notification_set_urgency(notify_low,NOTIFY_URGENCY_LOW);

	notify_normal=notify_notification_new_with_status_icon ("Normal",NULL,NULL,status_icon);
	notify_notification_set_icon_from_pixbuf(notify_normal,icon_main);
	notify_notification_set_timeout(notify_normal,NOTIFY_EXPIRES_DEFAULT);
	notify_notification_set_urgency(notify_normal,NOTIFY_URGENCY_NORMAL);

	notify_urgent=notify_notification_new_with_status_icon ("Urgent",NULL,"error",status_icon);
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
static void quitmonitor(GtkWidget *widget, gpointer data){
    (void)widget;
    (void)data;
#if WITH_NOTIFY
    if(notify_daemon)
	notify_uninit();
#endif
    quitall=1;
    add_host_wakeup();
    if(gtk_main_level ()>0){
	gtk_main_quit();
    }
    exit(0);
}
PROC_T *proc_get(int id,int pid){
    PROC_T *iproc;
    for(iproc=pproc[id]; iproc; iproc=iproc->next){
	if(iproc->pid==pid){
	    break;
	}
    }
    return iproc;
}
static void update_title(int id){
    char tit[40];
    snprintf(tit,40,"%s (%d)",hosts[id], nproc[id]);
    gtk_label_set_text(GTK_LABEL(titles[id]),tit);
}

PROC_T *proc_add(int id,int pid){
    PROC_T *iproc;
    if((iproc=proc_get(id,pid))) return iproc;
    iproc=calloc(1, sizeof(PROC_T));
    iproc->iseed_old=-1;
    iproc->pid=pid;
    iproc->hid=id;
    iproc->next=pproc[id];
    pproc[id]=iproc;
    nproc[id]++;
    update_title(id);
    return iproc;
}
void proc_remove_all(int id){
    PROC_T *iproc,*jproc=NULL;
    for(iproc=pproc[id]; iproc; iproc=jproc){
	jproc=iproc->next;
	remove_entry(iproc);
	if(iproc->path) free(iproc->path);
	free(iproc);
    }
    nproc[id]=0;
    pproc[id]=NULL;
    update_title(id);
}
void proc_remove(int id,int pid){
    PROC_T *iproc,*jproc=NULL;
    for(iproc=pproc[id]; iproc; jproc=iproc,iproc=iproc->next){
	if(iproc->pid==pid){
	    if(jproc){
		jproc->next=iproc->next;
	    }else{
		pproc[id]=iproc->next;
	    }
	    remove_entry(iproc);
	    if(iproc->path) free(iproc->path);
	    free(iproc);
	    nproc[id]--;
	    update_title(id);
	    break;
	}
    }
}

void kill_job(PROC_T *p){
	GtkWidget *dia=gtk_message_dialog_new
	    (NULL, GTK_DIALOG_DESTROY_WITH_PARENT,
	     GTK_MESSAGE_QUESTION,
	     GTK_BUTTONS_YES_NO,
	     "Kill job %d on server %s?",
	     p->pid,hosts[p->hid]);
	int result=gtk_dialog_run(GTK_DIALOG(dia));
	gtk_widget_destroy (dia);
	if(result==GTK_RESPONSE_YES){
	    if(scheduler_kill_job(p->hid,p->pid)){
		warning("Failed to kill the job\n");
	    }
	    p->status.info=S_TOKILL;
	    refresh(p);
	}
}
void kill_job_event(GtkWidget *btn, GdkEventButton *event, PROC_T *p){
    (void)btn;
    if(event->button==1){
	kill_job(p);
    }
}

static void status_icon_on_click(GtkStatusIcon *status_icon0, 
			       gpointer data){
    (void)status_icon0;
    (void)data;
    if(GTK_WIDGET_VISIBLE(window)){
	gtk_widget_hide(window);
    }else{
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
		     "popup-menu",GTK_SIGNAL_FUNC(trayIconPopup),menu);
    gtk_status_icon_set_tooltip(status_icon, ProgName);
    gtk_status_icon_set_visible(status_icon, TRUE);
}


static gboolean delete_window(void){
    if(gtk_status_icon_is_embedded(status_icon)){
	gtk_widget_hide(window);
	return TRUE;//do not quit
    }else{
	return FALSE;//quit.
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
static void clear_jobs_finished(GtkWidget *btn){
    (void)btn;
    int ihost=gtk_notebook_get_current_page (GTK_NOTEBOOK(notebook));
    PROC_T *iproc,*jproc;
    if(!pproc[ihost]) return;

    int sock=scheduler_connect(pproc[ihost]->hid,0,0);
    if(sock==-1) return;
    int cmd[2];
    cmd[0]=CMD_REMOVE;
    
    for(iproc=pproc[ihost]; iproc; iproc=jproc){
	jproc=iproc->next;
	if(iproc->status.info==S_FINISH){
	    cmd[1]=iproc->pid;
	    if(write(sock,cmd,2*sizeof(int))!=sizeof(int)*2){
		warning("write to socket %d failed\n",sock);
	    }
	    while (gtk_events_pending ())
		gtk_main_iteration ();
	    
	}
    }
    close(sock);
}
static void clear_jobs_crashed(GtkWidget *btn){
    (void)btn;
    int ihost=gtk_notebook_get_current_page (GTK_NOTEBOOK(notebook));
    PROC_T *iproc,*jproc;

    if(!pproc[ihost]) return;
    int sock=scheduler_connect(pproc[ihost]->hid,0,0);
    if(sock==-1) return;
    int cmd[2];
    cmd[0]=CMD_REMOVE;
    
    for(iproc=pproc[ihost]; iproc; iproc=jproc){
	jproc=iproc->next;
	if(iproc->status.info==S_CRASH || iproc->status.info==S_KILLED 
	   ||iproc->status.info==S_TOKILL){
	    cmd[1]=iproc->pid;
	    if(write(sock,cmd,2*sizeof(int))!=sizeof(int)*2){
		warning("write to socket %d failed\n",sock);
	    }
	    while (gtk_events_pending ())
		gtk_main_iteration ();
	    
	}
    }
    close(sock);
}
GtkWidget *monitor_new_entry_progress(void){
    GtkWidget *prog=gtk_entry_new();
    gtk_entry_set_editable(GTK_ENTRY(prog),FALSE);
    gtk_entry_set_has_frame(GTK_ENTRY(prog),FALSE);
    gtk_entry_set_width_chars(GTK_ENTRY(prog),12);
    static GtkBorder border={1,1,1,1};
    gtk_entry_set_inner_border(GTK_ENTRY(prog),&border);
#if GTK_MAJOR_VERSION>=3 || GTK_MINOR_VERSION >= 18
    gtk_widget_set_can_focus(prog,FALSE);
#else
    g_object_set(prog,"can-focus",FALSE,NULL);
#endif
    gtk_widget_modify_base(prog,GTK_STATE_NORMAL, &white);
    gtk_widget_modify_bg(prog,GTK_STATE_SELECTED, &blue);

    GtkRcStyle *style=gtk_widget_get_modifier_style(prog);
    if(!style){
	style=gtk_rc_style_new();
    }
    style->xthickness=0;//use the names specified in gtkrc.
    style->ythickness=0;
    gtk_widget_modify_style(prog, style);
    gtk_entry_set_alignment(GTK_ENTRY(prog),0.5);
    return prog;
}
GtkWidget *monitor_new_progress(int vertical, int length){
    GtkWidget *prog=gtk_progress_bar_new();
    if(vertical){
	gtk_progress_bar_set_orientation(GTK_PROGRESS_BAR(prog),
					 GTK_PROGRESS_BOTTOM_TO_TOP);
	gtk_widget_set_size_request(prog, 8,length);
    }else{
	gtk_widget_set_size_request(prog, length,12);
    }
    return prog;
}

int main(int argc, char *argv[])
{
    if(!g_thread_supported()){
	g_thread_init(NULL);
	gdk_threads_init();
    }
    gtk_init(&argc, &argv);
#if WITH_NOTIFY
    if(!notify_init("AOS Notification")){
	notify_daemon=0;
    }
#endif
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
    
 
    gtk_rc_parse_string(rc_string_treeview);
    gtk_rc_parse_string(rc_string_progress);
    create_status_icon();
    window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(window),"MAOS Monitor");
    GtkStyle *style=gtk_widget_get_style(window);
    bg=gdk_color_copy(&style->bg[GTK_STATE_NORMAL]);
    gtk_window_set_icon(GTK_WINDOW(window),icon_main);
    GtkWidget *toolbar=gtk_hbox_new(FALSE,0);
    GtkWidget *tool_clear_crashed=gtk_button_new_with_label("Clear crashed");
    GtkWidget *tool_clear_finished=gtk_button_new_with_label("Clear finished");
    GtkWidget *im=gtk_image_new_from_stock(GTK_STOCK_NO,GTK_ICON_SIZE_MENU);
    gtk_button_set_image(GTK_BUTTON(tool_clear_crashed), im);
    im=gtk_image_new_from_stock(GTK_STOCK_YES,GTK_ICON_SIZE_MENU);
    gtk_button_set_image(GTK_BUTTON(tool_clear_finished), im);

    g_signal_connect(tool_clear_crashed, "clicked", G_CALLBACK(clear_jobs_crashed),NULL);
    g_signal_connect(tool_clear_finished,"clicked", G_CALLBACK(clear_jobs_finished),NULL);
   
    gtk_box_pack_start(GTK_BOX(toolbar),tool_clear_crashed, FALSE,FALSE,0);
    gtk_box_pack_start(GTK_BOX(toolbar),tool_clear_finished,FALSE,FALSE,0);

    notebook=gtk_notebook_new();
    gtk_widget_show(notebook);
    gtk_notebook_set_scrollable((GtkNotebook*)notebook, TRUE);
    GtkWidget *vbox=gtk_vbox_new(FALSE,0);
    gtk_widget_show_all(toolbar);
    gtk_container_add(GTK_CONTAINER(window), vbox);
#if GTK_MAJOR_VERSION<=2 && GTK_MINOR_VERSION <20
    gtk_box_pack_start(GTK_BOX(vbox), toolbar,FALSE,FALSE,0);
#else
    //Newer GTK_NOTEBOOK has action widgets.
    gtk_notebook_set_action_widget(GTK_NOTEBOOK(notebook), toolbar, GTK_PACK_END);
#endif
    gtk_box_pack_start(GTK_BOX(vbox), notebook, TRUE,TRUE,0);
    g_signal_connect(window, "delete_event", G_CALLBACK (delete_window), NULL);
    g_signal_connect(window, "destroy", 
		     G_CALLBACK (quitmonitor), NULL);
    g_signal_connect(G_OBJECT (window), "window-state-event", 
		     G_CALLBACK (window_state_event), NULL);
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
    prog_cpu=calloc(nhost, sizeof(GtkWidget *));
    prog_mem=calloc(nhost, sizeof(GtkWidget *));
   
    for(int ihost=0; ihost<nhost; ihost++){
	char tit[40];
	snprintf(tit,40,"%s(0)",hosts[ihost]);
    	titles[ihost]=gtk_label_new(tit);

	
	GtkWidget *hbox0=gtk_hbox_new(FALSE,0);
	prog_cpu[ihost]=monitor_new_progress(1,16);
	prog_mem[ihost]=monitor_new_progress(1,16);
	gtk_box_pack_start(GTK_BOX(hbox0),prog_cpu[ihost], TRUE, TRUE, 1);
	gtk_box_pack_start(GTK_BOX(hbox0),prog_mem[ihost], TRUE, TRUE, 1);
	gtk_box_pack_start(GTK_BOX(hbox0),titles[ihost], TRUE, TRUE, 0);
	gtk_widget_show_all(hbox0);
	GtkWidget *eventbox=gtk_event_box_new();
	gtk_container_add(GTK_CONTAINER(eventbox),hbox0);
	/*Put the box below an eventbox so that clicking on the progressbar
	  switches tab also.*/
	gtk_widget_show_all(eventbox);
	gtk_event_box_set_above_child(GTK_EVENT_BOX(eventbox),TRUE);
	gtk_event_box_set_visible_window(GTK_EVENT_BOX(eventbox),FALSE);

	cmdconnect[ihost]
	    =gtk_button_new_with_label("Click to make connection.");
	g_signal_connect(cmdconnect[ihost],"clicked", G_CALLBACK(add_host_wakeup), NULL);
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
	update_title(ihost);
	gtk_widget_show_all(tabs[ihost]);
	gtk_notebook_set_current_page(GTK_NOTEBOOK(notebook), ihost);
    }
    pthread_mutex_init(&pmutex, NULL);
    pthread_cond_init(&pcond, NULL);
    if(NULL==g_thread_create((GThreadFunc)add_host_thread, NULL, FALSE, NULL)){
	warning("Thread creation failed.\n");
    }
    gdk_threads_enter();
    gtk_main();
    gdk_threads_leave();
    gdk_color_free(bg);
}

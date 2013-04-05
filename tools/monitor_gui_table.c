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
  Implementing Monitor using GtkTable. Deprecated.
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
#include "monitor.h"
GtkWidget **tables;
gint *nrows;
#define ncol 7
static void delete_hbox_event(GtkWidget *btn, GdkEventButton *event,PROC_T *p){
    (void)btn;
    if(event->button==1){
	scheduler_cmd(p->hid,p->pid,CMD_REMOVE);
    }
}
static GtkWidget *new_button(void){
    GtkWidget *ev=gtk_event_box_new();
    gtk_widget_set_events(ev, GDK_BUTTON_PRESS_MASK);
    gtk_widget_show_all(ev);
    return ev;
    }
static void change_button(PROC_T *p, const gchar *STOCKID, void (*func)){
    if(p->btn){
	g_signal_handler_disconnect(p->btn, p->btnhandler);
    }else{
	p->btn=new_button();
    }
    if(p->btnimage){
	gtk_image_clear(GTK_IMAGE(p->btnimage));
	gtk_image_set_from_stock(GTK_IMAGE(p->btnimage),STOCKID,GTK_ICON_SIZE_MENU);
    }else{
	p->btnimage=gtk_image_new_from_stock(STOCKID,GTK_ICON_SIZE_MENU);
	gtk_container_add(GTK_CONTAINER(p->btn), p->btnimage);
    }
    p->btnhandler=g_signal_connect(p->btn,"button_press_event",G_CALLBACK(func),(gpointer)p);
    gtk_widget_show_all(p->btn);
}
static GtkWidget* new_entry(const char *text, int width, gfloat align){
    GtkWidget *prog=monitor_new_entry_progress();
    if(!text) text="Not Available";
    gtk_entry_set_text(GTK_ENTRY(prog),text);
    gtk_entry_set_width_chars(GTK_ENTRY(prog),width);
    gtk_entry_set_alignment(GTK_ENTRY(prog),align);
    return prog;
}
static GtkWidget *new_label(const char *text, int width,float align){
    GtkWidget *prog=gtk_label_new(text);
    gtk_label_set_width_chars(GTK_LABEL(prog),width);
    gtk_misc_set_alignment(GTK_MISC(prog),align,0.5);
    return prog;
}

#define WIDTH_PID 24
#define WIDTH_PATH 20
#define WIDTH_ISEED 5
#define WIDTH_TIMING 25
#define WIDTH_ERRLO 7
#define WIDTH_ERRHI 7

static void create_entry(PROC_T *p){
    p->hbox=gtk_hbox_new(FALSE,0);
    char lb[12];
    char stime[80];
    snprintf(lb,12," %5d",p->pid);
    struct tm *tim=localtime(&(p->status.timstart));
    strftime(stime,80,"[%F %k:%M:%S]",tim);
    strcat(stime,lb);
    p->entry_pid=new_label(stime,WIDTH_PID,0);
    p->entry_path=new_label(p->path,WIDTH_PATH,1);
    gtk_label_set_selectable(GTK_LABEL(p->entry_path), TRUE);
    gtk_label_set_ellipsize(GTK_LABEL(p->entry_path),PANGO_ELLIPSIZE_START);
#if GTK_MAJOR_VERSION>=3 || GTK_MINOR_VERSION >= 12
    gtk_widget_set_tooltip_text(p->entry_path, p->path);
#endif
    p->entry_errlo=new_label("Lo (nm)",WIDTH_ERRLO,1);
    p->entry_errhi=new_label("Hi (nm)",WIDTH_ERRHI,1);
    p->entry_iseed=new_entry("iSEED",WIDTH_ISEED,0.5);
    p->entry_timing=new_entry("Timing",WIDTH_TIMING,1);
    /*kill_button_new(p); */
    change_button(p, GTK_STOCK_STOP, kill_job_event);
    int irow=nrows[p->hid];
    nrows[p->hid]++;
    gtk_table_resize(GTK_TABLE(tables[p->hid]), nrows[p->hid],ncol);
    grid_attach(tables[p->hid], p->entry_pid, 0,1,irow,irow+1,0);
    grid_attach(tables[p->hid], p->entry_path, 1,2,irow,irow+1,7);
    grid_attach(tables[p->hid], p->entry_errlo, 2,3,irow,irow+1,0);
    grid_attach(tables[p->hid], p->entry_errhi, 3,4,irow,irow+1,0);
    grid_attach(tables[p->hid], p->entry_iseed, 4,5,irow,irow+1,0);
    grid_attach(tables[p->hid], p->entry_timing, 5,6,irow,irow+1,0);
    grid_attach(tables[p->hid], p->btn, 6,7,irow,irow+1,0);
    gtk_widget_show_all(tables[p->hid]);
    gtk_notebook_set_current_page(GTK_NOTEBOOK(notebook),p->hid);
}
static void update_prog(PROC_T *p){
    /*Update the progress bar*/
    if(p->status.nseed>0){
	if(p->status.laps>0){
	    p->frac=(double)p->status.laps/(double)(p->status.rest+p->status.laps);
	}
    
	const double tkmean=p->status.scale;
	const long rest=p->status.rest;
	const long laps=p->status.laps;
	const long resth=rest/3600;
	const long restm=(rest-resth*3600)/60;
	const long lapsh=laps/3600;
	const long lapsm=(laps-lapsh*3600)/60;
	const double tot=p->status.tot*tkmean;
	if(p->iseed_old!=p->status.iseed){
	    char tmp[64];
	    snprintf(tmp,64,"%d of %d",p->status.iseed+1,p->status.nseed);
	    gtk_entry_set_text(GTK_ENTRY(p->entry_iseed), tmp);
	    p->iseed_old=p->status.iseed;
#if GTK_MAJOR_VERSION>=3 || GTK_MINOR_VERSION >= 16
	    gtk_entry_set_progress_fraction
		(GTK_ENTRY(p->entry_iseed),
		 (double)(p->status.iseed+1)/(double)p->status.nseed);
#endif
	}
	char tmp[64];
	snprintf(tmp,64, "%d of %d %ld:%02ld %ld:%02ld %.2fs",
		 p->status.isim+1,p->status.simend,
		 lapsh,lapsm,resth,restm,tot);
	gtk_entry_set_text(GTK_ENTRY(p->entry_timing),tmp);
	snprintf(tmp,64,"%.2f",p->status.clerrlo);
	gtk_label_set_text(GTK_LABEL(p->entry_errlo),tmp);
	snprintf(tmp,64,"%.2f",p->status.clerrhi);
	gtk_label_set_text(GTK_LABEL(p->entry_errhi),tmp);
#if GTK_MAJOR_VERSION>=3 || GTK_MINOR_VERSION >= 16
	gtk_entry_set_progress_fraction(GTK_ENTRY(p->entry_timing),p->frac);
#endif	
    }
}
gboolean remove_entry(PROC_T *iproc){
    /*Delete widget; */
    gtk_widget_destroy(iproc->entry_pid);
    gtk_widget_destroy(iproc->entry_path);
    gtk_widget_destroy(iproc->entry_errlo);
    gtk_widget_destroy(iproc->entry_errhi);
    gtk_widget_destroy(iproc->entry_iseed);
    gtk_widget_destroy(iproc->entry_timing);
    gtk_widget_destroy(iproc->btn);
    nrows[iproc->hid]--;
    free(iproc->path);
    free(iproc);
    return 0;
}
gboolean refresh(PROC_T *p){
    if(!p->entry_iseed) create_entry(p);

    switch(p->status.info){
    case 0:
	break;
    case S_RUNNING:
	break;
    case S_WAIT:
	/*waiting to start */
	gtk_entry_set_text(GTK_ENTRY(p->entry_timing),"Waiting");
	break;
    case S_START:
	/*just started. */
	gtk_entry_set_text(GTK_ENTRY(p->entry_timing),"Started");
	{
	    char lb[12];
	    char stime[80];
	    snprintf(lb,12," %5d",p->pid);
	    struct tm *tim=localtime(&(p->status.timstart));
	    strftime(stime,80,"[%F %k:%M:%S]",tim);
	    strcat(stime,lb);
	    gtk_label_set_text(GTK_LABEL(p->entry_pid), stime);
	}
	notify_user(p);
	break;
    case S_QUEUED:
	/*queued in scheduler */
	gtk_entry_set_text(GTK_ENTRY(p->entry_timing),"Queued");
	break;
    case S_FINISH:/*Finished */
	p->frac=1;
	change_button(p,GTK_STOCK_APPLY,delete_hbox_event);
	/*progress bar color. */
	gtk_widget_modify_bg(p->entry_timing,GTK_STATE_SELECTED,&green);
	notify_user(p);
	break;
    case S_CRASH:/*Error */
	gtk_entry_set_text(GTK_ENTRY(p->entry_timing),"Error");
	change_button(p,GTK_STOCK_CLOSE,delete_hbox_event);
	gtk_widget_modify_base(p->entry_timing,GTK_STATE_NORMAL,&red);
	notify_user(p);
	break;
    case S_TOKILL:/*kill command sent */
	gtk_entry_set_text(GTK_ENTRY(p->entry_timing),"Kill command sent");
	change_button(p,GTK_STOCK_CLOSE,delete_hbox_event);
	gtk_widget_modify_base(p->entry_timing,GTK_STATE_NORMAL,&yellow);
	break;
    case S_KILLED:
	gtk_entry_set_text(GTK_ENTRY(p->entry_timing),"Killed");
	change_button(p,GTK_STOCK_CLOSE,delete_hbox_event);
	gtk_widget_modify_base(p->entry_timing,GTK_STATE_NORMAL,&red);
	notify_user(p);
	break;
    default:
	warning("Unknown info\n");
    }
    update_prog(p);
    return 0;
}
GtkWidget *new_page(int ihost){
    if(!tables){
	tables=calloc(nhost, sizeof(GtkWidget*));
	nrows=calloc(nhost, sizeof(gint));
    }
    nrows[ihost]=0;
    tables[ihost]=gtk_table_new(nrows[ihost],ncol,0);
#if GTK_MAJOR_VERSION <3 || GTK_MINOR_VERSION < 4
    gtk_table_set_row_spacings(GTK_TABLE(tables[ihost]), 2);
    gtk_table_set_col_spacings(GTK_TABLE(tables[ihost]), 2);
#endif
    return tables[ihost];
}
void kill_selected_jobs(GtkAction *btn){
    (void)btn;
}

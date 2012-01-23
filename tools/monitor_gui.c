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
#include "scheduler_client.h"
#include "common.h"
#include "monitor.h"

#define WIDTH_START 13
#define WIDTH_PID 6
#define WIDTH_PATH 20
#define WIDTH_ISEED 4
#define WIDTH_TIMING 20
#define WIDTH_ERRLO 7
#define WIDTH_ERRHI 7

static void delete_hbox_event(GtkWidget *btn, GdkEventButton *event,PROC_T *p){
    (void)btn;
    if(event->button==1){
	scheduler_remove_job(p->hid,p->pid);
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

static void create_entry(PROC_T *p){
    p->hbox=gtk_hbox_new(FALSE,0);
    char lb[12];
    snprintf(lb,12," %5d",p->pid);
    struct tm *tim=localtime(&(p->status.timstart));
    char stime[80];
    strftime(stime,80,"[%a %k:%M:%S]",tim);
    p->entry_start=new_label(stime,WIDTH_START,0.5);
    p->entry_pid=new_label(lb,WIDTH_PID,1);
    p->entry_path=new_label(p->path,WIDTH_PATH,0);
    gtk_label_set_selectable(GTK_LABEL(p->entry_path), TRUE);
    gtk_label_set_ellipsize(GTK_LABEL(p->entry_path),PANGO_ELLIPSIZE_START);
#if GTK_MAJOR_VERSION>=3 || GTK_MINOR_VERSION >= 12
    gtk_widget_set_tooltip_text(p->entry_path, p->path);
#endif
    p->entry_errlo=new_label("Lo (nm)",WIDTH_ERRLO,1);
    p->entry_errhi=new_label("Hi (nm)",WIDTH_ERRHI,1);
    gtk_label_set_max_width_chars(GTK_LABEL(p->entry_errlo),WIDTH_ERRLO);
    gtk_label_set_max_width_chars(GTK_LABEL(p->entry_errhi),WIDTH_ERRHI);
    p->entry_iseed=new_entry("seed",WIDTH_ISEED,0.5);
    p->entry_timing=new_entry("Timing",WIDTH_TIMING,1);
    gtk_box_pack_start(GTK_BOX(p->hbox),gtk_vseparator_new(),FALSE,FALSE,0);
    gtk_box_pack_start(GTK_BOX(p->hbox),p->entry_start,FALSE,FALSE,0);
    gtk_box_pack_start(GTK_BOX(p->hbox),gtk_vseparator_new(),FALSE,FALSE,0);
    gtk_box_pack_start(GTK_BOX(p->hbox),p->entry_pid,FALSE,FALSE,0);
    gtk_box_pack_start(GTK_BOX(p->hbox),gtk_vseparator_new(),FALSE,FALSE,0);
    gtk_box_pack_start(GTK_BOX(p->hbox),p->entry_path,TRUE,TRUE,0);
    gtk_box_pack_start(GTK_BOX(p->hbox),gtk_vseparator_new(),FALSE,FALSE,0);
    gtk_box_pack_start(GTK_BOX(p->hbox),p->entry_errlo,FALSE,FALSE,0);
    gtk_box_pack_start(GTK_BOX(p->hbox),gtk_vseparator_new(),FALSE,FALSE,0);
    gtk_box_pack_start(GTK_BOX(p->hbox),p->entry_errhi,FALSE,FALSE,0);
    gtk_box_pack_start(GTK_BOX(p->hbox),gtk_vseparator_new(),FALSE,FALSE,0);
    gtk_box_pack_start(GTK_BOX(p->hbox),p->entry_iseed,FALSE,FALSE,0);
    gtk_box_pack_start(GTK_BOX(p->hbox),gtk_vseparator_new(),FALSE,FALSE,0);
    gtk_box_pack_start(GTK_BOX(p->hbox),p->entry_timing,FALSE,FALSE,0);
    gtk_box_pack_start(GTK_BOX(p->hbox),gtk_vseparator_new(),FALSE,FALSE,0);
    change_button(p, GTK_STOCK_STOP, kill_job_event);
    gtk_box_pack_start(GTK_BOX(p->hbox),p->btn,FALSE,FALSE,0);
    p->vbox=gtk_vbox_new(FALSE,0);
    if(nproc[p->hid]==1){
	gtk_box_pack_start(GTK_BOX(p->vbox),gtk_hseparator_new(),FALSE,FALSE,0);
    }
    gtk_box_pack_start(GTK_BOX(p->vbox),p->hbox,FALSE,FALSE,0);
    gtk_box_pack_start(GTK_BOX(p->vbox),gtk_hseparator_new(),FALSE,FALSE,0);
    gtk_box_pack_start(GTK_BOX(pages[p->hid]),p->vbox,FALSE,FALSE,0);
    
    /*gtk_notebook_set_current_page(GTK_NOTEBOOK(notebook),p->hid); */
    gtk_widget_show_all(p->vbox);
}
static void update_prog(PROC_T *p){
    /*Update the progress bar*/
    if(p->status.nseed>0){
	if(p->status.laps>0){
	    p->frac=(double)p->status.laps/(double)(p->status.rest+p->status.laps);
	}
    
	const double tkmean=p->status.scale;
	const long tot=p->status.rest+p->status.laps;/*show total time. */
	const long toth=tot/3600;
	const long totm=(tot-toth*3600)/60;
	const long tots=tot-toth*3600-totm*60;
	const long rest=p->status.rest;
	const long resth=rest/3600;
	const long restm=(rest-resth*3600)/60;
	const long rests=rest-resth*3600-restm*60;
	const double step=p->status.tot*tkmean;
	if(p->iseed_old!=p->status.iseed){
	    char tmp[64];
	    snprintf(tmp,64,"%d/%d",p->status.iseed+1,p->status.nseed);
	    gtk_entry_set_text(GTK_ENTRY(p->entry_iseed), tmp);
	    p->iseed_old=p->status.iseed;
#if GTK_MAJOR_VERSION>=3 || GTK_MINOR_VERSION >= 16
	    gtk_entry_set_progress_fraction
		(GTK_ENTRY(p->entry_iseed),
		 (double)(p->status.iseed+1)/(double)p->status.nseed);
#endif
	}
	char tmp[64];
	/*snprintf(tmp,64, "%5.2fs %d/%d %2ld:%02ld/%ld:%02ld",
	  step,p->status.isim+1,p->status.simend,
	  restm,rests,totm,tots);*/
	if(toth>99){
	    snprintf(tmp,64, "%d %5.2fs %ldh/%ldh",p->status.simend,
		     step, resth,toth);
	}else if(toth>0){
	    snprintf(tmp,64, "%d %5.2fs %ldh%02ld/%ldh%02ld",p->status.simend,
		     step, resth,restm,toth,totm);
	}else{
	    snprintf(tmp,64, "%d %5.2fs %2ld:%02ld/%ld:%02ld",p->status.simend,
		     step, restm,rests,totm,tots);	
	}
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
void remove_entry(PROC_T *iproc){
    /*Delete widget; */
    if(iproc->vbox){
	/*warning3("destroy hbox\n"); */
	gtk_widget_destroy(iproc->vbox);
	iproc->vbox=NULL;
    }else{
	/*warning("hbox is empty\n"); */
    }
}
void refresh(PROC_T *p){
    if(p->status.info==S_REMOVE){
	proc_remove(p->hid,p->pid);
	return;
    }
    if(!p->entry_iseed) create_entry(p);
    if(p->done) return;
    switch(p->status.info){
    case S_RUNNING:
	break;
    case S_WAIT: /*waiting to start */
	gtk_entry_set_text(GTK_ENTRY(p->entry_timing),"Waiting to start");
	break;
    case S_START: /*just started. */
	gtk_entry_set_text(GTK_ENTRY(p->entry_timing),"Started");
	notify_user(p);
	break;
    case S_FINISH:/*Finished */
	p->done=1;
	change_button(p,GTK_STOCK_APPLY,delete_hbox_event);
	gtk_widget_modify_bg(p->entry_timing,GTK_STATE_SELECTED,&green);/*progress bar color.(filled bkgrnd) */
	gtk_widget_modify_base(p->entry_timing,GTK_STATE_NORMAL,&green);/*progress bar color.(empty bkgrnd) */
	gtk_entry_set_text(GTK_ENTRY(p->entry_timing),"Skipped");
	notify_user(p);
	break;
    case S_CRASH:/*Error */
	p->done=1;
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
	p->done=1;
	gtk_entry_set_text(GTK_ENTRY(p->entry_timing),"Killed");
	change_button(p,GTK_STOCK_CLOSE,delete_hbox_event);
	gtk_widget_modify_base(p->entry_timing,GTK_STATE_NORMAL,&red);
	notify_user(p);
	break;
    default:
	warning("Unknown info\n");
    }
    update_prog(p);
}
GtkWidget *new_page(int ihost){
    (void)ihost;
    return NULL;
}

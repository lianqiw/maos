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
#define NOTIFY_IS_SUPPORTED 1

#include <stdio.h>
#include <stdlib.h>
#include <netdb.h>
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
#include <gdk/gdk.h>
#include <glib/gprintf.h>
#include "common.h"
#include "monitor.h"
#include "misc.h"
//#include "gtkcellrendererprogressnew.h" //modify appearance of progress
static GtkListStore **lists=NULL;
static void list_get_iter(PROC_T *p, GtkTreeIter *iter){
    GtkListStore *list=lists[p->hid];
    GtkTreePath* tpath=gtk_tree_row_reference_get_path (p->row);
    gtk_tree_model_get_iter(GTK_TREE_MODEL(list),iter,tpath);
    gtk_tree_path_free(tpath);
}
static void list_modify_icon(PROC_T *p, GdkPixbuf *newicon){
    GtkTreeIter iter;
    GtkListStore *list=lists[p->hid];
    list_get_iter(p, &iter);
    gtk_list_store_set(list, &iter, COL_ACTION, newicon,-1);
}
static void list_modify_color(PROC_T *p, const char *color){
    GtkTreeIter iter;
    GtkListStore *list=lists[p->hid];
    list_get_iter(p, &iter);
    gtk_list_store_set(list, &iter, COL_COLOR, color,-1);
}
static void list_modify_status(PROC_T *p, const char *status){
    GtkTreeIter iter;
    GtkListStore *list=lists[p->hid];
    list_get_iter(p, &iter);
    gtk_list_store_set(list, &iter, COL_TIMING, status,-1);
}
/*static void list_modify_progress(PROC_T *p, double prog){
    GtkTreeIter iter;
    GtkListStore *list=lists[p->hid];
    list_get_iter(p, &iter);
    gtk_list_store_set(list, &iter, COL_STEPP, prog,-1);
    }*/
static void list_update_progress(PROC_T *p){
    if(p->status.nseed==0) return;
    GtkListStore *list=lists[p->hid];
    double tot=(double)(p->status.rest+p->status.laps);
    if(fabs(tot)>1.e-10){
	p->frac=(double)p->status.laps/tot;
    }else{
	p->frac=0;
    }
    GtkTreeIter iter;
    list_get_iter(p, &iter);
	    
    const double tkmean=p->status.scale;
    const long rest=p->status.rest;
    const long laps=p->status.laps;
    const long resth=rest/3600;
    const long restm=(rest-resth*3600)/60;
    const long lapsh=laps/3600;
    const long lapsm=(laps-lapsh*3600)/60;
    tot=p->status.tot*tkmean;
    if(p->status.iseed!=p->iseed_old-1){
	char tmp[64];
	snprintf(tmp,64,"%d/%d",p->status.iseed+1,p->status.nseed);
	gtk_list_store_set(list, &iter, 
			   COL_SEED, tmp, 
			   COL_SEEDP,
			   100*(double)(p->status.iseed+1.)/(double)p->status.nseed,
			   -1);
	p->iseed_old=p->status.iseed+1;
    }
    char tmp[64];
    snprintf(tmp,64,"%d/%d",p->status.isim+1,p->status.simend);
    gtk_list_store_set(list, &iter, 
		       COL_STEP,tmp, 
		       COL_STEPP,p->frac*100, -1);
    snprintf(tmp,64,"%.2f",tot);
    gtk_list_store_set(list, &iter, COL_TIMING,tmp, -1);
    
    snprintf(tmp,64,"%ld:%02ld", lapsh,lapsm);
    gtk_list_store_set(list, &iter, COL_LAPS, tmp, -1);
    
    snprintf(tmp,64,"%ld:%02ld", resth,restm);
    gtk_list_store_set(list, &iter, COL_REST, tmp, -1);
    
    snprintf(tmp,64,"%.2f",p->status.clerrlo);
    gtk_list_store_set(list, &iter, COL_ERRLO,tmp, -1);
    snprintf(tmp,64,"%.2f",p->status.clerrhi);
    gtk_list_store_set(list, &iter, COL_ERRHI,tmp, -1);
}
void remove_entry(PROC_T *p){
    GtkTreePath *path=gtk_tree_row_reference_get_path(p->row);
    GtkTreeIter iter;
    if(!gtk_tree_model_get_iter(GTK_TREE_MODEL(lists[p->hid]),&iter,path)){
	warning("Unable to find entry");
    }else{
	gtk_list_store_remove (lists[p->hid],&iter);
    }
}
void refresh(PROC_T *p){
    if(p->status.info==S_REMOVE){
	proc_remove(p->hid,p->pid);
	return;
    }
    if(p->done) return;//we are done with it (finished or crashed)
    if(!p->row){
	char sdate[80];
	char spid[12];
	snprintf(spid,12," %d ",p->pid);
	time_t t;
	t=myclocki();
	struct tm *tim=localtime(&t);
	strftime(sdate,80,"%m-%d %k:%M:%S",tim);
	char *spath=p->path;
	if(!spath) spath="Unknown";
	GtkListStore *list=lists[p->hid];
	GtkTreeIter iter;
	gtk_list_store_append(list, &iter);
	gtk_list_store_set(list,&iter,
			   COL_DATE,sdate,
			   COL_PID, spid,
			   COL_PATH,spath,
			   COL_SEED,"",
			   COL_STEP,"",
			   -1);
	GtkTreePath *tpath=gtk_tree_model_get_path(GTK_TREE_MODEL(list), &iter);
	p->row=gtk_tree_row_reference_new(GTK_TREE_MODEL(list),tpath);
	gtk_tree_path_free(tpath);
	list_update_progress(p);
    }
    switch(p->status.info){
    case S_RUNNING:
	list_update_progress(p);
	list_modify_icon(p,icon_running);
	break;
    case S_WAIT:
	//waiting to start
	//list_modify_status(p, "Queued");
	break;
    case S_START:
	//just started.
	list_modify_status(p, "Started");
	notify_user(p);
	break;
    case S_FINISH://Finished
	//p->frac=1;
	//p->done=1;
	//list_modify_progress(p,100);
	list_update_progress(p);
	//list_modify_status(p, "Finished");
	list_modify_icon(p, icon_finished);
	list_modify_color(p,"#00DD00");
	//gtk_widget_modify_base(p->prog3,GTK_STATE_NORMAL,&green);
	notify_user(p);
	break;
    case S_CRASH://Error
	p->done=1;
	list_modify_status(p, "Error");
	list_modify_icon(p,icon_failed);
	list_modify_color(p,"#CC0000");
	notify_user(p);
	break;
    case S_TOKILL://kill command sent
	list_modify_status(p, "Kill command sent");
	list_modify_icon(p,icon_failed);
	list_modify_color(p,"#CCCC00");
	break;
    case S_KILLED:
	p->done=1;
	list_modify_status(p, "Killed");
	list_modify_icon(p,icon_failed);
	list_modify_color(p,"#CC0000");
	notify_user(p);
	break;
    default:
	warning("Unknown info: %d\n",p->status.info);
    }
}

static void action_clicked(GtkTreeViewColumn *col){
    (void)col;
    info("clicked\n");
}
 
GtkWidget *new_page(int ihost){
    if(!lists){
	lists=calloc(nhost,sizeof(GtkListStore*));
    }

    lists[ihost]=gtk_list_store_new(COL_TOT,
				    G_TYPE_STRING,//DATE
				    G_TYPE_STRING,//PID
				    G_TYPE_STRING,//PATH
				    G_TYPE_STRING,//SEED
				    G_TYPE_FLOAT, //SEEDP
				    G_TYPE_STRING,//STEP
				    G_TYPE_FLOAT, //STEPP
				    G_TYPE_STRING,//TIMING
				    G_TYPE_STRING,//LAPS
				    G_TYPE_STRING,//REST
				    G_TYPE_STRING,//ERRLO
				    G_TYPE_STRING,//ERRHI
				    GDK_TYPE_PIXBUF,//ACTION
				    G_TYPE_STRING //COLOR
				    );
    GtkWidget *view;
    view=gtk_tree_view_new_with_model(GTK_TREE_MODEL(lists[ihost]));
    g_object_unref(lists[ihost]);
    GtkCellRenderer *render;
    GtkTreeViewColumn *col;
    GtkTreeSelection *viewsel=gtk_tree_view_get_selection(GTK_TREE_VIEW(view));
    //gtk_tree_selection_set_select_function(viewsel, treeselfun,NULL,NULL);
    gtk_tree_selection_set_mode(viewsel,GTK_SELECTION_MULTIPLE);
    /*
      The implementation of GtkTreeView hardcoded GDK_LINE_ON_OFF_DASH in
      gtk_tree_view_set_grid_lines, which makes it impossible to make solid
      lines. In rc_string (monitor.c) "\255\256" in grid-line-pattern means 255
      pixels of on and 1 pixel of off. I can not have \0 for the second part
      because it terminates the string.
    */
    int spacing=0;
    float align=0.5;
    //gtk_tree_view_set_grid_lines(GTK_TREE_VIEW(view), GTK_TREE_VIEW_GRID_LINES_VERTICAL);
    gtk_tree_view_set_enable_search(GTK_TREE_VIEW(view), TRUE);
    //g_object_set(G_OBJECT(view),"rules-hint", TRUE, NULL);
    gtk_tree_view_set_rules_hint(GTK_TREE_VIEW(view), TRUE);

    render=gtk_cell_renderer_text_new();
    g_object_set(G_OBJECT(render), "ypad", 0, NULL);
    col=gtk_tree_view_column_new();
    gtk_tree_view_column_set_spacing(col, spacing);
    gtk_tree_view_column_set_alignment(col,align);
    gtk_tree_view_column_pack_start(col,render,FALSE);
    gtk_tree_view_column_add_attribute(col,render,"text",COL_DATE);
    gtk_tree_view_column_set_title(col,"Date");
    gtk_tree_view_append_column(GTK_TREE_VIEW(view),col);

    render=gtk_cell_renderer_text_new();
    g_object_set(G_OBJECT(render), "ypad", 0, NULL);
    col=gtk_tree_view_column_new();
    gtk_tree_view_column_set_spacing(col, spacing);
    gtk_tree_view_column_set_alignment(col,align);
    gtk_tree_view_column_pack_start(col,render,FALSE);
    gtk_tree_view_column_add_attribute(col,render,"text",COL_PID);
    gtk_tree_view_column_set_title(col,"PID");
    gtk_tree_view_append_column(GTK_TREE_VIEW(view),col);

    render=gtk_cell_renderer_text_new();
    g_object_set(G_OBJECT(render), "ypad", 0, NULL);
    col=gtk_tree_view_column_new();
    gtk_tree_view_column_set_spacing(col, spacing);
    gtk_tree_view_column_set_alignment(col,align);
    gtk_tree_view_column_pack_start(col,render,TRUE);
    gtk_tree_view_column_add_attribute(col,render,"text",COL_PATH);
    gtk_tree_view_column_set_title(col,"Path");
    gtk_tree_view_column_set_min_width(col,20);
    gtk_tree_view_column_set_expand(col,TRUE);
    g_object_set(G_OBJECT(render),"ellipsize",PANGO_ELLIPSIZE_START,NULL);
    gtk_tree_view_append_column(GTK_TREE_VIEW(view),col);


    render=gtk_cell_renderer_text_new();
    g_object_set(G_OBJECT(render), "ypad", 0, NULL);
    col=gtk_tree_view_column_new();
    gtk_tree_view_column_set_spacing(col, spacing);
    gtk_tree_view_column_set_alignment(col,align);
    gtk_tree_view_column_pack_start(col,render,TRUE);
    gtk_tree_view_column_add_attribute(col,render,"text",COL_ERRLO);
    gtk_tree_view_column_set_title(col,"ErrLo");
    gtk_tree_view_append_column(GTK_TREE_VIEW(view),col);

    render=gtk_cell_renderer_text_new();
    g_object_set(G_OBJECT(render), "ypad", 0, NULL);
    col=gtk_tree_view_column_new();
    gtk_tree_view_column_set_spacing(col, spacing);
    gtk_tree_view_column_set_alignment(col,align);
    gtk_tree_view_column_pack_start(col,render,TRUE);
    gtk_tree_view_column_add_attribute(col,render,"text",COL_ERRHI);
    gtk_tree_view_column_set_title(col,"ErrHi");
    gtk_tree_view_append_column(GTK_TREE_VIEW(view),col);
	
    render=gtk_cell_renderer_progress_new();
    g_object_set(G_OBJECT(render), "ypad", 0, NULL);
    col=gtk_tree_view_column_new();
    gtk_tree_view_column_set_spacing(col, spacing);
    gtk_tree_view_column_set_alignment(col,align);
    gtk_tree_view_column_pack_start(col,render,TRUE);
    gtk_tree_view_column_add_attribute(col,render,"text",COL_SEED);
    gtk_tree_view_column_add_attribute(col,render,"value",COL_SEEDP);
    gtk_tree_view_column_add_attribute(col,render,"cell-background",COL_COLOR);
    gtk_tree_view_column_set_title(col,"Seed");
    gtk_tree_view_append_column(GTK_TREE_VIEW(view),col);
   
    render=gtk_cell_renderer_progress_new();
    g_object_set(G_OBJECT(render), "ypad", 0, NULL);
    gtk_cell_renderer_set_fixed_size(render,30,5);
    col=gtk_tree_view_column_new();
    gtk_tree_view_column_set_spacing(col, spacing);
    gtk_tree_view_column_set_alignment(col,align);
    gtk_tree_view_column_pack_start(col,render,TRUE);
    gtk_tree_view_column_add_attribute(col,render,"text",COL_STEP);
    gtk_tree_view_column_add_attribute(col,render,"value",COL_STEPP);
    gtk_tree_view_column_add_attribute(col,render,"cell-background",COL_COLOR);
    gtk_tree_view_column_set_title(col,"Progress");
    gtk_tree_view_column_set_min_width(col,80);
    gtk_tree_view_column_set_max_width(col,80);
    gtk_tree_view_append_column(GTK_TREE_VIEW(view),col);
    //gtk_tree_view_column_set_min_width(GTK_TREE_VIEW_COLUMN(col),200);

    render=gtk_cell_renderer_text_new();
    g_object_set(G_OBJECT(render), "ypad", 0, NULL);
    col=gtk_tree_view_column_new();
    gtk_tree_view_column_set_spacing(col, spacing);
    gtk_tree_view_column_set_alignment(col,align);
    gtk_tree_view_column_pack_start(col,render,TRUE);
    gtk_tree_view_column_add_attribute(col,render,"text",COL_TIMING);
    //gtk_tree_view_column_add_attribute(col,render,"foreground",COL_COLOR);
    //gtk_tree_view_column_add_attribute(col,render,"cell-background",COL_COLOR);
    gtk_tree_view_column_set_title(col,"Time");
    gtk_tree_view_append_column(GTK_TREE_VIEW(view),col);

    render=gtk_cell_renderer_text_new();
    g_object_set(G_OBJECT(render), "ypad", 0, NULL);
    col=gtk_tree_view_column_new();
    gtk_tree_view_column_set_spacing(col, spacing);
    gtk_tree_view_column_set_alignment(col,align);
    gtk_tree_view_column_pack_start(col,render,TRUE);
    gtk_tree_view_column_add_attribute(col,render,"text",COL_REST);
    //gtk_tree_view_column_add_attribute(col,render,"foreground",COL_COLOR);
    //gtk_tree_view_column_add_attribute(col,render,"cell-background",COL_COLOR);
    gtk_tree_view_column_set_title(col,"Left");
    gtk_tree_view_append_column(GTK_TREE_VIEW(view),col);
  
    render=gtk_cell_renderer_text_new();
    g_object_set(G_OBJECT(render), "ypad", 0, NULL);
    col=gtk_tree_view_column_new();
    gtk_tree_view_column_set_spacing(col, spacing);
    gtk_tree_view_column_set_alignment(col,align);
    gtk_tree_view_column_pack_start(col,render,TRUE);
    gtk_tree_view_column_add_attribute(col,render,"text",COL_LAPS);
    //gtk_tree_view_column_add_attribute(col,render,"foreground",COL_COLOR);
    //gtk_tree_view_column_add_attribute(col,render,"cell-background",COL_COLOR);
    gtk_tree_view_column_set_title(col,"Used");
    gtk_tree_view_append_column(GTK_TREE_VIEW(view),col);

 
    render=gtk_cell_renderer_pixbuf_new();
    g_object_set(G_OBJECT(render), "ypad", 0, NULL);
    col=gtk_tree_view_column_new();
    gtk_tree_view_column_set_spacing(col, spacing);
    gtk_tree_view_column_set_alignment(col,align);
    gtk_tree_view_column_pack_start(col,render,TRUE);
    gtk_tree_view_column_add_attribute(col,render,"pixbuf",COL_ACTION);
    gtk_tree_view_column_set_title(col,"   ");
    //gtk_tree_view_column_set_clickable(col,TRUE);
    gtk_tree_view_append_column(GTK_TREE_VIEW(view),col);
    g_signal_connect(col, "clicked",G_CALLBACK(action_clicked),NULL);


    return view;
}

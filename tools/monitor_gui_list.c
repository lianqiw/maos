/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "monitor.h"
enum{
    COL_DATE,
    COL_PID,
    COL_PATH,
    COL_SEED,
    COL_SEEDP,/*float */
    COL_STEP,
    COL_STEPP,/*float */
    COL_TIMING,
    COL_ALL,
    COL_REST,
    COL_ERRLO,
    COL_ERRHI,
    COL_ACTION,
    COL_COLOR,
    COL_TOT,
};
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
    double total=(double)(p->status.rest+p->status.laps);
    if(fabs(total)>1.e-10){
	p->frac=(double)p->status.laps/total;
    }else{
	p->frac=0;
    }
    GtkTreeIter iter;
    list_get_iter(p, &iter);
	    
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
    if(p->status.iseed!=p->iseed_old){
	char tmp[64];
	snprintf(tmp,64,"%d/%d",p->status.iseed+1,p->status.nseed);
	gtk_list_store_set(list, &iter, 
			   COL_SEED, tmp, 
			   COL_SEEDP,
			   100*(double)(p->status.iseed+1.)/(double)p->status.nseed,
			   -1);
	p->iseed_old=p->status.iseed;
    }
    char tmp[64];
    snprintf(tmp,64,"%5.2f", step); gtk_list_store_set(list, &iter, COL_TIMING,tmp, -1);
    if(toth>99){
	snprintf(tmp,64,"%ldh", resth);  gtk_list_store_set(list, &iter, COL_REST,tmp, -1);
	snprintf(tmp,64,"%ldh", toth);   gtk_list_store_set(list, &iter, COL_ALL,tmp, -1);
    }else if(toth>0){
	snprintf(tmp,64,"%ldh%02ld", resth, restm);  gtk_list_store_set(list, &iter, COL_REST,tmp, -1);
	snprintf(tmp,64,"%ldh%02ld", toth, totm);    gtk_list_store_set(list, &iter, COL_ALL,tmp, -1);
    }else{
 	snprintf(tmp,64,"%02ld:%02ld", restm, rests);  gtk_list_store_set(list, &iter, COL_REST,tmp, -1);
	snprintf(tmp,64,"%02ld:%02ld", totm, tots);    gtk_list_store_set(list, &iter, COL_ALL,tmp, -1);
     }
    
    if(toth>99){
	snprintf(tmp,64, "%d/%d %5.2fs %ldh/%ldh",p->status.isim+1,p->status.simend, step, resth,toth);
    }else if(toth>0){
	snprintf(tmp,64, "%d/%d %5.2fs %ldh%02ld/%ldh%02ld",p->status.isim+1,p->status.simend, step, resth,restm,toth,totm);
    }else{
	snprintf(tmp,64, "%d/%d %5.2fs %2ld:%02ld/%ld:%02ld",p->status.isim+1,p->status.simend, step, restm,rests,totm,tots);	
    }
    //snprintf(tmp,64,"%d/%d",p->status.isim+1,p->status.simend);
    gtk_list_store_set(list, &iter, COL_STEP,tmp, COL_STEPP,p->frac*100, -1);
    
    snprintf(tmp,64,"%.2f",p->status.clerrlo);
    gtk_list_store_set(list, &iter, COL_ERRLO,tmp, -1);
    snprintf(tmp,64,"%.2f",p->status.clerrhi);
    gtk_list_store_set(list, &iter, COL_ERRHI,tmp, -1);
    
}
gboolean remove_entry(PROC_T *p){
    GtkTreePath *path=gtk_tree_row_reference_get_path(p->row);
    GtkTreeIter iter;
    if(!gtk_tree_model_get_iter(GTK_TREE_MODEL(lists[p->hid]),&iter,path)){
	warning("Unable to find entry");
    }else{
	gtk_list_store_remove (lists[p->hid],&iter);
    }
    free(p->path);
    free(p);
    return 0;
}
gboolean refresh(PROC_T *p){
    if(p->done) return 0;/*we are done with it (finished or crashed) */
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
	list_update_progress(p);
	gtk_tree_path_free(tpath);
    }
    switch(p->status.info){
    case S_RUNNING:
	list_update_progress(p);
	list_modify_icon(p,icon_running);
	break;
    case S_WAIT: /*waiting to start */
	break;
    case S_START: /*just started. */
	list_modify_status(p, "Started");
	notify_user(p);
	break;
    case S_FINISH:/*Finished */
	list_update_progress(p);
	list_modify_icon(p, icon_finished);
	list_modify_color(p,"#00DD00");
	notify_user(p);
	break;
    case S_CRASH:/*Error */
	p->done=1;
	list_modify_status(p, "Error");
	list_modify_icon(p,icon_failed);
	list_modify_color(p,"#CC0000");
	notify_user(p);
	break;
    case S_TOKILL:/*kill command sent */
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
    return 0;
}

static  GtkTreeViewColumn *new_column(int type, int width, const char *title, ...){
    GtkTreeViewColumn *col;
    GtkCellRenderer *render=NULL;
    col=gtk_tree_view_column_new();
    switch(type){
    case 0:
	render=gtk_cell_renderer_text_new();
#if GTK_MAJOR_VERSION>=3 || GTK_MINOR_VERSION >= 18
	gtk_cell_renderer_set_padding(render, 1, 1);
	gtk_cell_renderer_set_alignment(render, 1, 0.5);
#endif
	break;
    case 1:
	render=gtk_cell_renderer_progress_new();
	break;
    case 2:
	render=gtk_cell_renderer_pixbuf_new();
	break;
    default:
	error("Invalid\n");
    }
    gtk_tree_view_column_set_title(col, title);
    gtk_tree_view_column_set_spacing(col, 2);
    gtk_tree_view_column_set_alignment(col, 1);
    if(width>0){/*minimum width*/
	gtk_tree_view_column_set_min_width(col, width);
	gtk_tree_view_column_set_expand(col,TRUE);
	if(type==0){
	    g_object_set(G_OBJECT(render),"ellipsize",PANGO_ELLIPSIZE_START,NULL);
	}
    }else if(width<0){/*exact width*/
	gtk_tree_view_column_set_min_width(col,-width);
	gtk_tree_view_column_set_max_width(col,-width);
    }
    gtk_tree_view_column_pack_start(col,render,TRUE);
    va_list ap;
    va_start(ap, title);
    const char *att=NULL;
    while((att=va_arg(ap, const char *))){
	gtk_tree_view_column_add_attribute(col, render, att, va_arg(ap,int));
    }
    va_end(ap);
    return col;
}
GtkWidget *new_page(int ihost){
    if(!lists){
	lists=calloc(nhost,sizeof(GtkListStore*));
    }

    lists[ihost]=gtk_list_store_new(COL_TOT,
				    G_TYPE_STRING,/*DATE */
				    G_TYPE_STRING,/*PID */
				    G_TYPE_STRING,/*PATH */
				    G_TYPE_STRING,/*SEED */
				    G_TYPE_FLOAT, /*SEEDP */
				    G_TYPE_STRING,/*STEP */
				    G_TYPE_FLOAT, /*STEPP */
				    G_TYPE_STRING,/*TIMING */
				    G_TYPE_STRING,/*TOT */
				    G_TYPE_STRING,/*REST */
				    G_TYPE_STRING,/*ERRLO */
				    G_TYPE_STRING,/*ERRHI */
				    GDK_TYPE_PIXBUF,/*ACTION */
				    G_TYPE_STRING /*COLOR */
				    );
    GtkWidget *view;
    view=gtk_tree_view_new_with_model(GTK_TREE_MODEL(lists[ihost]));
    g_object_unref(lists[ihost]);
    GtkTreeSelection *viewsel=gtk_tree_view_get_selection(GTK_TREE_VIEW(view));
    /*gtk_tree_selection_set_select_function(viewsel, treeselfun,NULL,NULL); */
    gtk_tree_selection_set_mode(viewsel,GTK_SELECTION_MULTIPLE);
    /*
      The implementation of GtkTreeView hardcoded GDK_LINE_ON_OFF_DASH in
      gtk_tree_view_set_grid_lines, which makes it impossible to make solid
      lines. In rc_string (monitor.c) "\255\256" in grid-line-pattern means 255
      pixels of on and 1 pixel of off. I can not have \0 for the second part
      because it terminates the string.
    */
    /*gtk_tree_view_set_grid_lines(GTK_TREE_VIEW(view), GTK_TREE_VIEW_GRID_LINES_VERTICAL); */
    gtk_tree_view_set_enable_search(GTK_TREE_VIEW(view), TRUE);
    /*g_object_set(G_OBJECT(view),"rules-hint", TRUE, NULL); */
    gtk_tree_view_set_rules_hint(GTK_TREE_VIEW(view), TRUE);
    
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 0, "Date", "text", COL_DATE, NULL));
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 0, "PID" , "text", COL_PID, NULL));
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 20,"Path", "text", COL_PATH, NULL));
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, -46, "Low" , "text", COL_ERRLO, NULL));
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, -46, "High", "text", COL_ERRHI, NULL));
    /*gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 0, "Step", "text", COL_TIMING, NULL));
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 0, "Left", "text", COL_REST, NULL));
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 0, "Tot", "text", COL_ALL, NULL));*/
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(1, 0, "Seed", "text", COL_SEED, "value",COL_SEEDP,NULL));
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(1, 0, "Progress", "text", COL_STEP, "value",COL_STEPP,"cell-background",COL_COLOR,NULL));
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(2, 0, " ", "pixbuf", COL_ACTION, NULL));
    
    return view;
}

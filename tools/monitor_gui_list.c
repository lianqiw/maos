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
static GtkWidget **views=NULL;
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
/*static void list_modify_color(PROC_T *p, const char *color){
    GtkTreeIter iter;
    GtkListStore *list=lists[p->hid];
    list_get_iter(p, &iter);
    gtk_list_store_set(list, &iter, COL_COLOR, color,-1);
    }*/
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
    if(p->row){
	GtkTreePath *path=gtk_tree_row_reference_get_path(p->row);
	GtkTreeIter iter;
	if(!gtk_tree_model_get_iter(GTK_TREE_MODEL(lists[p->hid]),&iter,path)){
	    warning("Unable to find entry");
	}else{
	    gtk_list_store_remove (lists[p->hid],&iter);
	}
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
	//time_t t;
	//t=myclocki();
	struct tm *tim=localtime(&p->status.timstart);
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
    case 0:
	break;
    case S_RUNNING:
	list_update_progress(p);
	list_modify_icon(p,icon_running);
	break;
    case S_WAIT: /*waiting to start */
	break;
    case S_START: /*just started. */
	list_modify_status(p, "Started");
	{
	    GtkTreeIter iter;
	    GtkListStore *list=lists[p->hid];
	    list_get_iter(p, &iter);
	    char spid[12];
	    snprintf(spid,12," %d ",p->pid);
	
	    char sdate[80];
	    struct tm *tim=localtime(&p->status.timstart);
	    strftime(sdate,80,"%m-%d %k:%M:%S",tim);
	    gtk_list_store_set(list, &iter, 
			       COL_PID,spid, 
			       COL_DATE, sdate,
			       -1);
	}
	notify_user(p);
	break;
    case S_QUEUED:
	break;
    case S_FINISH:/*Finished */
	list_update_progress(p);
	list_modify_icon(p, icon_finished);
	//list_modify_color(p,"#00DD00");
	notify_user(p);
	break;
    case S_CRASH:/*Error */
	p->done=1;
	list_modify_status(p, "Error");
	list_modify_icon(p,icon_failed);
	//list_modify_color(p,"#CC0000");
	notify_user(p);
	break;
    case S_TOKILL:/*kill command sent */
	list_modify_status(p, "Kill command sent");
	list_modify_icon(p,icon_failed);
	//list_modify_color(p,"#CCCC00");
	break;
    case S_KILLED:
	p->done=1;
	list_modify_status(p, "Killed");
	list_modify_icon(p,icon_failed);
	//list_modify_color(p,"#CC0000");
	notify_user(p);
	break;
    case S_REMOVE:
	break;
    default:
	warning("Unknown info: %d\n",p->status.info);
    }
    return 0;
}

static GtkTreeViewColumn *new_column(int type, int width, const char *title, ...){
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
/**
   Clear the clipboard so we can append to it multiple times.
 */
static void clipboard_clear(){
    GdkAtom atoms[2]={GDK_SELECTION_CLIPBOARD, GDK_SELECTION_PRIMARY};
    for(int iatom=0; iatom<2; iatom++){
	GtkClipboard *clip=gtk_clipboard_get(atoms[iatom]);
	gtk_clipboard_set_text(clip,"",-1);
    }
}
/**
   Append text to clipboard (primary and default).
*/
static void clipboard_append(const char *jobinfo){
    GdkAtom atoms[2]={GDK_SELECTION_CLIPBOARD, GDK_SELECTION_PRIMARY};
    for(int iatom=0; iatom<2; iatom++){
	GtkClipboard *clip=gtk_clipboard_get(atoms[iatom]);
	gchar *old=gtk_clipboard_wait_for_text(clip);
	gchar *newer=stradd(old, jobinfo, "\n", NULL);
	gtk_clipboard_set_text(clip, newer, -1);
	g_free(newer);
	g_free(old);
    }
}
static void handle_selected(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer user_data, int cmd, char *action){
    gint ihost=GPOINTER_TO_INT(user_data);
    GValue value=G_VALUE_INIT;
    gtk_tree_model_get_value(model, iter, COL_PID, &value);
    int pid=strtol(g_value_get_string(&value), NULL, 10);
    g_value_unset(&value);
    if(cmd<0){
	switch(cmd){
	case -1:{
	   
	    gtk_tree_model_get_value(model, iter, COL_PATH, &value);
	    gchar *jobinfo=g_strdup(g_value_get_string(&value));
	    g_value_unset(&value);
	    if(!strcmp(action,"CopyPath")){
		char *pos=NULL;
		pos=strstr(jobinfo, "/maos ");
		if(!pos){
		    pos=strstr(jobinfo, "/skyc ");
		}
		if(pos){
		    pos[1]='\0';
		    gchar *tmp=jobinfo;
		    jobinfo=stradd("cd ", tmp, NULL);
		    g_free(tmp);
		}
	    }
	    clipboard_append(jobinfo);
	    g_free(jobinfo);
	}
	    break;
	}
    }else{
	if(scheduler_cmd(ihost,pid,cmd)){
	    warning("Failed to %s the job\n", action);
	}
    }
}
static void kill_selected(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer user_data){
    handle_selected(model, path, iter, user_data, CMD_KILL, "Kill");
}
static void restart_selected(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer user_data){
    handle_selected(model, path, iter, user_data, CMD_RESTART, "Restart");
}
static void clear_selected(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer user_data){
    handle_selected(model, path, iter, user_data, CMD_REMOVE, "Remove");
}
static void copy_selected(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer user_data){
    handle_selected(model, path, iter, user_data, -1, "Copy");
}
static void copy_selectedpath(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer user_data){
    handle_selected(model, path, iter, user_data, -1, "CopyPath");
}

static void handle_selected_event(GtkMenuItem *menuitem, gpointer user_data, GtkTreeSelectionForeachFunc func, char *action){
    gint ihost=GPOINTER_TO_INT(user_data);
    GtkWidget *view=views[ihost];
    GtkTreeSelection *selection=gtk_tree_view_get_selection(GTK_TREE_VIEW(view));
    int nsel=gtk_tree_selection_count_selected_rows(selection);
    if(nsel<1) {
	warning("nsel=%d\n", nsel);
	return;
    }
    int result=GTK_RESPONSE_YES;
    if(!strcmp(action, "Kill") || !strcmp(action, "Restart")){
	GtkWidget *dia=gtk_message_dialog_new
	    (NULL, GTK_DIALOG_DESTROY_WITH_PARENT,
	     GTK_MESSAGE_QUESTION,
	     GTK_BUTTONS_YES_NO,
	     "%s %d jobs?", action, nsel);
	result=gtk_dialog_run(GTK_DIALOG(dia));
	gtk_widget_destroy (dia);
    }
    if(result==GTK_RESPONSE_YES){
	gtk_tree_selection_selected_foreach(selection, func, GINT_TO_POINTER(ihost));
    }
}

static void kill_selected_event(GtkMenuItem *menuitem, gpointer user_data){
    handle_selected_event(menuitem, user_data, kill_selected, "Kill");
}
static void restart_selected_event(GtkMenuItem *menuitem, gpointer user_data){
    handle_selected_event(menuitem, user_data, restart_selected, "Restart");
}
static void clear_selected_event(GtkMenuItem *menuitem, gpointer user_data){
    handle_selected_event(menuitem, user_data, clear_selected, "Clear");
}
static void copy_selected_event(GtkMenuItem *menuitem, gpointer user_data){
    clipboard_clear();
    handle_selected_event(menuitem, user_data, copy_selected, "Copy");
}
static void copy_selectedpath_event(GtkMenuItem *menuitem, gpointer user_data){
    clipboard_clear();
    handle_selected_event(menuitem, user_data, copy_selectedpath, "Copy Path of ");
}

static gboolean view_popup_menu(GtkWidget *view, gpointer user_data){
    GtkTreeSelection *selection=gtk_tree_view_get_selection(GTK_TREE_VIEW(view));
    int nsel=gtk_tree_selection_count_selected_rows(selection);

    GtkWidget *menu=gtk_menu_new();
    GtkWidget *menuitem;
    char text[40];
    snprintf(text, 40, "%d selected", nsel);
    menuitem=gtk_menu_item_new_with_label(text);
    gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);
    menuitem=gtk_separator_menu_item_new();
    gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);
    if(nsel>0){
	menuitem=gtk_image_menu_item_new_from_stock(GTK_STOCK_CANCEL, NULL);
	gtk_menu_item_set_label(GTK_MENU_ITEM(menuitem), "Kill selected jobs");
	g_signal_connect(menuitem, "activate", G_CALLBACK(kill_selected_event), user_data);
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);

	menuitem=gtk_image_menu_item_new_from_stock(GTK_STOCK_MEDIA_PLAY, NULL);
	gtk_menu_item_set_label(GTK_MENU_ITEM(menuitem), "Restart selected jobs");
	g_signal_connect(menuitem, "activate", G_CALLBACK(restart_selected_event), user_data);
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);

	menuitem=gtk_image_menu_item_new_from_stock(GTK_STOCK_CLEAR, NULL);
	gtk_menu_item_set_label(GTK_MENU_ITEM(menuitem), "Clear selected jobs");
	g_signal_connect(menuitem, "activate", G_CALLBACK(clear_selected_event), user_data);
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);

	menuitem=gtk_image_menu_item_new_from_stock(GTK_STOCK_COPY, NULL);
	gtk_menu_item_set_label(GTK_MENU_ITEM(menuitem), "Copy selected jobs");
	g_signal_connect(menuitem, "activate", G_CALLBACK(copy_selected_event), user_data);
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);

	menuitem=gtk_image_menu_item_new_from_stock(GTK_STOCK_COPY, NULL);
	gtk_menu_item_set_label(GTK_MENU_ITEM(menuitem), "Copy path of selected jobs");
	g_signal_connect(menuitem, "activate", G_CALLBACK(copy_selectedpath_event), user_data);
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);
    }
    gtk_widget_show_all(menu);
    gtk_menu_popup(GTK_MENU(menu), NULL, NULL, NULL, NULL, 3, gtk_get_current_event_time());
    return TRUE;
}
static gboolean view_click_event(GtkWidget *view, GdkEventButton *event, gpointer user_data){
    if(event->button==3){//right menu
	GtkTreeSelection *selection=gtk_tree_view_get_selection(GTK_TREE_VIEW(view));
	int nsel=gtk_tree_selection_count_selected_rows(selection);
	if (nsel < 1){
	    GtkTreePath *path;
	    if(gtk_tree_view_get_path_at_pos(GTK_TREE_VIEW(view),
					     (gint)event->x,
					     (gint)event->y,
					     &path, NULL, NULL, NULL)){
		gtk_tree_selection_select_path(selection, path);
		gtk_tree_path_free(path);
	    }
	}
	view_popup_menu(view, user_data);
	return TRUE;
    }else{
	return FALSE;
    }
}
GtkWidget *new_page(int ihost){
    if(!lists){
	lists=calloc(nhost,sizeof(GtkListStore*));
	views=calloc(nhost, sizeof(GtkWidget*));
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
    views[ihost]=view=gtk_tree_view_new_with_model(GTK_TREE_MODEL(lists[ihost]));
    g_object_unref(lists[ihost]);
    g_signal_connect(view, "button-press-event", 
		     G_CALLBACK(view_click_event), GINT_TO_POINTER(ihost));
    g_signal_connect(view, "popup-menu", 
		     G_CALLBACK(view_popup_menu), GINT_TO_POINTER(ihost));
    GtkTreeSelection *viewsel=gtk_tree_view_get_selection(GTK_TREE_VIEW(view));
    /*gtk_tree_selection_set_select_function(viewsel, treeselfun,NULL,NULL); */
    gtk_tree_selection_set_mode(viewsel,GTK_SELECTION_MULTIPLE);
    gtk_tree_view_set_tooltip_column(GTK_TREE_VIEW(view), COL_PATH);
    /*
      The implementation of GtkTreeView hardcoded GDK_LINE_ON_OFF_DASH in
      gtk_tree_view_set_grid_lines, which makes it impossible to make solid
      lines. In rc_string (monitor.c) "\255\256" in grid-line-pattern means 255
      pixels of on and 1 pixel of off. I can not have \0 for the second part
      because it terminates the string.
    */
    gtk_tree_view_set_grid_lines(GTK_TREE_VIEW(view), GTK_TREE_VIEW_GRID_LINES_VERTICAL);
    gtk_tree_view_set_enable_search(GTK_TREE_VIEW(view), TRUE);
    /*g_object_set(G_OBJECT(view),"rules-hint", TRUE, NULL); */
    gtk_tree_view_set_rules_hint(GTK_TREE_VIEW(view), TRUE);
    
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 0, "Date", "text", COL_DATE, NULL));
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 0, "PID" , "text", COL_PID, NULL));
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 20,"Path", "text", COL_PATH, NULL));
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, -50, "Low" , "text", COL_ERRLO, NULL));
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, -50, "High", "text", COL_ERRHI, NULL));
    /*gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 0, "Step", "text", COL_TIMING, NULL));
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 0, "Left", "text", COL_REST, NULL));
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 0, "Tot", "text", COL_ALL, NULL));*/
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(1, 0, "Seed", "text", COL_SEED, "value",COL_SEEDP,NULL));
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(1, 0, "Progress", "text", COL_STEP, "value",COL_STEPP,"cell-background",COL_COLOR,NULL));
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(2, 0, " ", "pixbuf", COL_ACTION, NULL));
    
    return view;
}


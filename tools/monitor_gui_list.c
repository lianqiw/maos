/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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



#include <tgmath.h>
#include <netdb.h>
#include <netdb.h>
#include <fcntl.h> 
#include <errno.h>
#include <arpa/inet.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>

#include <glib.h>
#include <gtk/gtk.h>
#include <gdk/gdk.h>
#include <glib/gprintf.h>
#include "monitor.h"
/*DO not modify enum without modify list store*/
enum{
    COL_DATE=0,
    COL_PID,
    COL_FULL,/*full path+arguments+output*/
    COL_START,/*starting path*/
    COL_ARGS,/*arguments*/
    COL_OUT, /*output folder*/
    COL_SEED,
    COL_SEEDP,/*gint */
    COL_STEP,
    COL_STEPP,/*gint */
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
static void list_modify_color(PROC_T *p, const char *color){
    GtkTreeIter iter;
    GtkListStore *list=lists[p->hid];
    list_get_iter(p, &iter);
    gtk_list_store_set(list, &iter, COL_COLOR, color,-1);
}

static void list_update_progress(PROC_T *p){
    if(p->status.nseed==0) return;
    GtkListStore *list=lists[p->hid];
    double total=(double)(p->status.rest+p->status.laps);
    if(fabs(total)>1.e-10){
	p->frac=(double)p->status.laps/total;
    }else{
	p->frac=0;
    }
    if(p->frac>1){
	p->frac=1;
    }else if(p->frac<0){
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
			   COL_SEEDP, (100*(p->status.iseed+1)/p->status.nseed),
			   -1);
	p->iseed_old=p->status.iseed;
    }
    char tmp[64];
  
    if(toth>99){
	snprintf(tmp,64, "%d/%d %.3fs %ldh/%ldh",p->status.isim+1,p->status.simend, step, resth,toth);
    }else if(toth>0){
	snprintf(tmp,64, "%d/%d %.3fs %ldh%02ld/%ldh%02ld",p->status.isim+1,p->status.simend, step, resth,restm,toth,totm);
    }else{
	snprintf(tmp,64, "%d/%d %.3fs %2ld:%02ld/%ld:%02ld",p->status.isim+1,p->status.simend, step, restm,rests,totm,tots);	
    }
    gtk_list_store_set(list, &iter, 
		       COL_STEP,tmp, 
		       COL_STEPP,(gint)(p->frac*100), 
		       -1);
    
    snprintf(tmp,64,"%.2f",p->status.clerrlo);
    gtk_list_store_set(list, &iter, COL_ERRLO,tmp, -1);
    snprintf(tmp,64,"%.2f",p->status.clerrhi);
    gtk_list_store_set(list, &iter, COL_ERRHI,tmp, -1);
    
}
static void list_modify_reset(PROC_T *p){
    GtkTreeIter iter;
    GtkListStore *list=lists[p->hid];
    list_get_iter(p, &iter);
    char spid[12];
    snprintf(spid,12," %d ",p->pid);
	
    char sdate[80];
    struct tm *tim=localtime(&p->status.timstart);
    strftime(sdate,80,"%m-%d %k:%M:%S",tim);
    gtk_list_store_set(list, &iter, 
		       COL_DATE, sdate,
		       COL_PID, spid, 
		       COL_SEED," ",
		       COL_SEEDP, 0, //Don't use 0.
		       COL_STEP," ",
		       COL_STEPP, 0,
		       COL_ERRHI," ",
		       COL_ERRLO," ",
		       -1);
    p->iseed_old=-1;
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
    if(!p->row){
	char sdate[80];
	char spid[12];
	snprintf(spid,12," %d ",p->pid);
	struct tm *tim=localtime(&p->status.timstart);
	strftime(sdate,80,"%m-%d %k:%M:%S",tim);
	char *spath=p->path;
	char *sstart=NULL, *sout=NULL, *sargs=NULL;
	if(spath){
	    const char *pos=NULL;
	    pos=strstr(spath, "/maos ");
	    if(!pos){
		pos=strstr(spath, "/skyc ");
	    }
	    if(pos){
		sstart=(char*)malloc(pos-spath+1);
		memcpy(sstart, spath, pos-spath);
		sstart[pos-spath]='\0';
		sargs=strdup(pos+1);
		char *pos2=NULL;
		for(char *tmp=sargs; (tmp=strstr(tmp, " -o ")); tmp+=4){
		    pos2=tmp;
		}
		if(pos2){
		    pos2+=4;
		    char *pos3=strchr(pos2, ' ');
		    if(!pos3) pos3=strchr(pos2, '\0');
		    if(pos3){
			sout=(char*)malloc(pos3-pos2+1);
			memcpy(sout, pos2, pos3-pos2);
			sout[pos3-pos2]='\0';
			memmove(pos2-4, pos3, strlen(pos3)+1);
		    }
		}
	    }
	}
	GtkListStore *list=lists[p->hid];
	GtkTreeIter iter;
	gtk_list_store_append(list, &iter);
	gtk_list_store_set(list,&iter,
			   COL_DATE,  sdate,
			   COL_PID,   spid,
			   COL_FULL,  spath?spath:"Unknown",
			   COL_START, sstart?sstart:" ",
			   COL_ARGS,  sargs?sargs:" ",
			   COL_OUT,   sout?sout:" ",
			   COL_ERRHI, " ",
			   COL_ERRLO, " ",
			   COL_SEED,  " ",
			   COL_SEEDP, 0,
			   COL_STEP,  " ",
			   COL_STEPP, 0,
			   -1);
	free(sstart); free(sout); free(sargs);
	GtkTreePath *tpath=gtk_tree_model_get_path(GTK_TREE_MODEL(list), &iter);
	p->row=gtk_tree_row_reference_new(GTK_TREE_MODEL(list),tpath);
	list_update_progress(p);
	gtk_tree_path_free(tpath);
	gtk_tree_view_columns_autosize(GTK_TREE_VIEW(views[p->hid]));
    }
    switch(p->status.info){
    case 0:
	break;
    case S_RUNNING:
	list_update_progress(p);
	list_modify_icon(p,icon_running);
	break;
    case S_WAIT: /*waiting to start */
	//list_modify_status(p, "Waiting");
	list_modify_reset(p);
	list_modify_icon(p,icon_waiting);
	break;
    case S_START: /*just started. */
	//list_modify_status(p, "Started");
	list_modify_reset(p);
	list_modify_icon(p,icon_running);
	notify_user(p);
	break;
    case S_QUEUED:
	//list_modify_status(p, "Queued");
	list_modify_reset(p);
	list_modify_icon(p,icon_waiting);
	break;
    case S_FINISH:/*Finished */
	list_update_progress(p);
	list_modify_icon(p, icon_finished);
	//list_modify_color(p,"#00DD00");
	notify_user(p);
	break;
    case S_CRASH:/*Error */
	//list_modify_status(p, "Error");
	list_modify_icon(p,icon_failed);
	//list_modify_color(p,"#CC0000");
	notify_user(p);
	break;
    case S_TOKILL:/*kill command sent */
	//list_modify_status(p, "Kill command sent");
	list_modify_icon(p,icon_failed);
	//list_modify_color(p,"#CCCC00");
	break;
    case S_KILLED:
	//list_modify_status(p, "Killed");
	list_modify_icon(p,icon_failed);
	//list_modify_color(p,"#CC0000");
	notify_user(p);
	break;
    case S_REMOVE:
	break;
    default:
	warning("Unknown info: %d\n",p->status.info);
    }
    if(p->status.warning){
	list_modify_color(p, "#FF0000");
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
    //resizeable makes the column very small if not expand.
    gtk_tree_view_column_set_resizable(col, (width)?TRUE:FALSE);
    gtk_tree_view_column_set_expand(col,(width==-1)?TRUE:FALSE);

    //column only hides text when 1) maxwidth is set, 2) ellipsize is set
    if(width>0){
	gtk_tree_view_column_set_min_width(col,width);
	//gtk_tree_view_column_set_max_width(col,width*5);
    }
    if(type==0 && width && width!=-2){
	//set ellipsize makes it prefer to shrink
	g_object_set(G_OBJECT(render),"ellipsize",PANGO_ELLIPSIZE_START,NULL);
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
	if(!clip) continue;
	gchar *old=gtk_clipboard_wait_for_text(clip);
	gchar *newer=NULL;
	if(old){
	    newer=stradd(old, jobinfo, "\n", NULL);
	    g_free(old);
	}else{
	    newer=stradd(jobinfo, "\n", NULL);
	}
	gtk_clipboard_set_text(clip, newer, -1);
	free(newer);
    }
}
/*A general routine handle actions to each item.*/
static void handle_selected(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer user_data, int cmd, const char *action){
    (void)path;
    gint ihost=GPOINTER_TO_INT(user_data);
    GValue value=G_VALUE_INIT;
    gtk_tree_model_get_value(model, iter, COL_PID, &value);
    int pid=strtol(g_value_get_string(&value), NULL, 10);
    g_value_unset(&value);
    if(cmd<0){
	switch(cmd){
	case -1:{
	    if(!strcmp(action, "CopyPath")){
		gtk_tree_model_get_value(model, iter, COL_START, &value);
	    }else if(!strcmp(action, "CopyOutPath")){
		gtk_tree_model_get_value(model, iter, COL_OUT, &value);
	    }else if(!strcmp(action, "Copy")){
		gtk_tree_model_get_value(model, iter, COL_FULL, &value);
	    }
	    gchar *jobinfo=g_strdup(g_value_get_string(&value));
	    g_value_unset(&value);
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
/*wrapper routines for each item.*/
static void kill_selected(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer user_data){
    handle_selected(model, path, iter, user_data, CMD_KILL, "Kill");
}
static void restart_selected(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer user_data){
    handle_selected(model, path, iter, user_data, CMD_RESTART, "Restart");
}
static void plot_selected(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer user_data){
    handle_selected(model, path, iter, user_data, CMD_DISPLAY, "Plot");
}
static void clear_selected(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer user_data){
    handle_selected(model, path, iter, user_data, CMD_REMOVE, "Remove");
}
static void copy_selected(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer user_data){
    handle_selected(model, path, iter, user_data, -1, "Copy");
}
static void copy_selected_path(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer user_data){
    handle_selected(model, path, iter, user_data, -1, "CopyPath");
}
static void copy_selected_outpath(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer user_data){
    handle_selected(model, path, iter, user_data, -1, "CopyOutPath");
}
/*Handles menu item clicks in a general way.*/
static void handle_selected_event(GtkMenuItem *menuitem, gpointer user_data, GtkTreeSelectionForeachFunc func, const char *action){
    (void)menuitem;
    gint ihost=GPOINTER_TO_INT(user_data);
    GtkWidget *view=views[ihost];
    GtkTreeSelection *selection=gtk_tree_view_get_selection(GTK_TREE_VIEW(view));
    int nsel=gtk_tree_selection_count_selected_rows(selection);
    if(nsel<1) {
	warning("nsel=%d\n", nsel);
	return;
    }
    int result=GTK_RESPONSE_YES;
    /*Alert user in Kill or Restart event*/
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
/* The following routines respond to menu item clicks.*/
static void kill_selected_event(GtkMenuItem *menuitem, gpointer user_data){
    handle_selected_event(menuitem, user_data, kill_selected, "Kill");
}
static void restart_selected_event(GtkMenuItem *menuitem, gpointer user_data){
    handle_selected_event(menuitem, user_data, restart_selected, "Restart");
}
static void plot_selected_event(GtkMenuItem *menuitem, gpointer user_data){
    handle_selected_event(menuitem, user_data, plot_selected, "Plot");
}
static void clear_selected_event(GtkMenuItem *menuitem, gpointer user_data){
    handle_selected_event(menuitem, user_data, clear_selected, "Clear");
}
static void copy_selected_event(GtkMenuItem *menuitem, gpointer user_data){
    clipboard_clear();
    handle_selected_event(menuitem, user_data, copy_selected, "Copy");
}
static void copy_selected_path_event(GtkMenuItem *menuitem, gpointer user_data){
    clipboard_clear();
    handle_selected_event(menuitem, user_data, copy_selected_path, "Copy Path of ");
}
static void copy_selected_outpath_event(GtkMenuItem *menuitem, gpointer user_data){
    clipboard_clear();
    handle_selected_event(menuitem, user_data, copy_selected_outpath, "Copy Output Path of ");
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
	menuitem=gtk_menu_item_new_with_label("Kill selected jobs");
	g_signal_connect(menuitem, "activate", G_CALLBACK(kill_selected_event), user_data);
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);

	menuitem=gtk_menu_item_new_with_label("Restart selected jobs");
	g_signal_connect(menuitem, "activate", G_CALLBACK(restart_selected_event), user_data);
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);

	menuitem=gtk_menu_item_new_with_label("Plot selected jobs");
	g_signal_connect(menuitem, "activate", G_CALLBACK(plot_selected_event), user_data);
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);

	menuitem=gtk_separator_menu_item_new();
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);

	menuitem=gtk_menu_item_new_with_label("Clear selected jobs");
	g_signal_connect(menuitem, "activate", G_CALLBACK(clear_selected_event), user_data);
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);

	menuitem=gtk_menu_item_new_with_label("Copy cmdline selected jobs");
	g_signal_connect(menuitem, "activate", G_CALLBACK(copy_selected_event), user_data);
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);

	menuitem=gtk_menu_item_new_with_label("Copy path of selected jobs");
	g_signal_connect(menuitem, "activate", G_CALLBACK(copy_selected_path_event), user_data);
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);

	menuitem=gtk_menu_item_new_with_label("Copy output path of selected jobs");
	g_signal_connect(menuitem, "activate", G_CALLBACK(copy_selected_outpath_event), user_data);
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);
	
    }
    gtk_widget_show_all(menu);
    gtk_menu_popup(GTK_MENU(menu), NULL, NULL, NULL, NULL, 3, gtk_get_current_event_time());
    return TRUE;
}
/* Handle click on treeview */
static gboolean view_click_event(GtkWidget *view, GdkEventButton *event, gpointer user_data){
    if(event->button!=3) return FALSE;//let system handle left click
    GtkTreePath *path=NULL;
    GtkTreeViewColumn *column=NULL;
    GtkTreeSelection *selection=gtk_tree_view_get_selection(GTK_TREE_VIEW(view));
    if(gtk_tree_view_get_path_at_pos(GTK_TREE_VIEW(view),
				     (gint)event->x,
				     (gint)event->y,
				     &path, &column, NULL, NULL)){
	/*clicks within treeview.*/
	if(!gtk_tree_selection_path_is_selected(selection, path)){
	    /*the path clicked is not selected. move selection to this path*/
	    gtk_tree_selection_unselect_all(selection);
	    gtk_tree_selection_select_path(selection, path);
	}
	gtk_tree_path_free(path);
    }else{//clicks outsize of treeview. unselect all.
	gtk_tree_selection_unselect_all(selection);
	path=NULL; column=NULL;
    }
    //right click popup menu
    view_popup_menu(view, user_data);
    return TRUE;
}
static gboolean view_release_event(GtkWidget *view, GdkEventButton *event, gpointer user_data){
    if(event->button!=1) return FALSE;
    GtkTreePath *path=NULL;
    GtkTreeViewColumn *column=NULL;
    if(gtk_tree_view_get_path_at_pos(GTK_TREE_VIEW(view),
				     (gint)event->x,
				     (gint)event->y,
				     &path, &column, NULL, NULL)){
	const gchar *title=gtk_tree_view_column_get_title(column);
	if(!strcmp(title," ")){
	    GtkTreeIter iter;
	    gint ihost=GPOINTER_TO_INT(user_data);
	    GtkTreeModel *model=GTK_TREE_MODEL(lists[ihost]);
	    gtk_tree_model_get_iter(model, &iter, path);
	    GValue value=G_VALUE_INIT;
	    gtk_tree_model_get_value(model, &iter, COL_PID, &value);
	    int pid=strtol(g_value_get_string(&value), NULL, 10);
	    g_value_unset(&value);
	    PROC_T *p=proc_get(ihost, pid);
	    if(p){
		if(p->status.info<10){
		    kill_job(p);
		}else{
		    scheduler_cmd(ihost, pid, CMD_REMOVE);
		}
	    }
	}
    }
    return TRUE;
}
GtkWidget *new_page(int ihost){
    if(!lists){
	lists=mycalloc(nhost,GtkListStore*);
	views=mycalloc(nhost,GtkWidget*);
    }

    lists[ihost]=gtk_list_store_new(COL_TOT,
				    G_TYPE_STRING,/*DATE */
				    G_TYPE_STRING,/*PID */
				    G_TYPE_STRING,/*FULL */
				    G_TYPE_STRING,/*PATH */
				    G_TYPE_STRING,/*ARGS */
				    G_TYPE_STRING,/*OUT */
				    G_TYPE_STRING,/*SEED */
				    G_TYPE_INT, /*SEEDP */
				    G_TYPE_STRING,/*STEP */
				    G_TYPE_INT, /*STEPP */
				    G_TYPE_STRING,/*TIMING */
				    G_TYPE_STRING,/*ALL */
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
    g_signal_connect(view, "button-release-event", 
		     G_CALLBACK(view_release_event), GINT_TO_POINTER(ihost));
    g_signal_connect(view, "popup-menu", 
		     G_CALLBACK(view_popup_menu), GINT_TO_POINTER(ihost));
    GtkTreeSelection *viewsel=gtk_tree_view_get_selection(GTK_TREE_VIEW(view));
    /*gtk_tree_selection_set_select_function(viewsel, treeselfun,NULL,NULL); */
    gtk_tree_selection_set_mode(viewsel,GTK_SELECTION_MULTIPLE);
    gtk_tree_view_set_tooltip_column(GTK_TREE_VIEW(view), COL_FULL);
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
    //gtk_tree_view_set_headers_clickable(GTK_TREE_VIEW(view), TRUE);
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 0, "Date", "text", COL_DATE, NULL));
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 0, "PID" , "text", COL_PID, NULL));
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 100,"Path", "text", COL_START, NULL));
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, -1,"Args", "text", COL_ARGS, NULL));
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, -2,"Out", "text",  COL_OUT, NULL));
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 0, "Low" , "text", COL_ERRLO,"foreground",COL_COLOR, NULL));
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 0, "High", "text", COL_ERRHI,"foreground",COL_COLOR, NULL));
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(1, 0, "Seed", "text", COL_SEED, "value",COL_SEEDP,NULL));
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(1, 0, "Progress", "text", COL_STEP, "value",COL_STEPP,NULL));
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(2, 0, " ", "pixbuf", COL_ACTION, NULL));
    
    return view;
}


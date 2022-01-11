/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
	This file is a supplement to monitor.c It defines and operations an treeview to organize jobs.
*/
#define NOTIFY_IS_SUPPORTED 1
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
	COL_DATE=0,//date in string
	COL_TIME,//time in time_t. //only set once, to keep task order stable.
	COL_PID,
	COL_FULL,/*full path+arguments+output*/
	COL_PATH,/*starting path*/
	COL_ARGS,/*arguments*/
	COL_OUT, /*output folder*/
	//COL_SEED,
	//COL_SEEDP,/*gint */
	COL_STEPT,
	COL_STEP,
	COL_STEPP,/*gint */
	COL_TIMING,
	COL_ALL,
	COL_REST,
	COL_ERRLO,
	COL_ERRHI,
	COL_ACTION,
	COL_COLOR,
	COL_HOST,
	COL_TOT,
};
static GtkListStore* listall=NULL;
static GtkTreeModel** lists=NULL;
static GtkWidget** views=NULL;
static int list_get_iter(GtkTreeRowReference* row, GtkTreeIter* iter){
	GtkTreePath* tpath=NULL;
	if(row && (tpath=gtk_tree_row_reference_get_path(row))){
		gtk_tree_model_get_iter(GTK_TREE_MODEL(listall), iter, tpath);
		gtk_tree_path_free(tpath);
		return 0;
	}
	return 1;
}
static void list_modify_icon(GtkTreeRowReference* row, GdkPixbuf* newicon){
	GtkTreeIter iter;
	if(!list_get_iter(row, &iter)){
		gtk_list_store_set(listall, &iter, COL_ACTION, newicon, -1);
	}
}
static void list_modify_color(GtkTreeRowReference* row, const char* color){
	GtkTreeIter iter;
	if(!list_get_iter(row, &iter)){
		gtk_list_store_set(listall, &iter, COL_COLOR, color, -1);
	}
}

static void list_proc_update(proc_t* p){
	if(!p || p->status.nseed==0) return;
	GtkTreeIter iter;
	if(list_get_iter(p->row, &iter)) return;//failure

	double total=(double)(p->status.rest+p->status.laps);
	if(fabs(total)>1.e-10){
		p->frac=(double)p->status.laps/total;
	} else{
		p->frac=0;
	}
	if(p->frac>1){
		p->frac=1;
	} else if(p->frac<0){
		p->frac=0;
	}

	const long tot=p->status.rest+p->status.laps;/*show total time. */
	const long toth=tot/3600;
	const long totm=(tot-toth*3600)/60;
	const long tots=tot-toth*3600-totm*60;
	const long rest=p->status.rest;
	const long resth=rest/3600;
	const long restm=(rest-resth*3600)/60;
	const long rests=rest-resth*3600-restm*60;
	const double step=p->status.tot*p->status.scale;//scale is normally 1
	/*if(p->status.iseed!=p->iseed_old){
		char tmp[64];
		snprintf(tmp, 64, "%d/%d", p->status.iseed+1, p->status.nseed);
		gtk_list_store_set(listall, &iter,
			COL_SEED, tmp,
			COL_SEEDP, (100*(p->status.iseed+1)/p->status.nseed),
			-1);
		p->iseed_old=p->status.iseed;
	}*/
	char stmp[8];
	if(snprintf(stmp, 8, "%.2f", step)>=8){
		warning_once("format is truncated\n");
	}
	char tmp[64];
	size_t np=0;
	if(toth>99){
		np=snprintf(tmp, 64, "%d/%d %d/%d %ldh/%ldh", p->status.iseed+1, p->status.nseed, p->status.isim+1, p->status.simend, resth, toth);
	} else if(toth>0){
		np=snprintf(tmp, 64, "%d/%d %d/%d %ldh%02ld/%ldh%02ld", p->status.iseed+1, p->status.nseed, p->status.isim+1, p->status.simend, resth, restm, toth, totm);
	} else{
		np=snprintf(tmp, 64, "%d/%d %d/%d %02ld:%02ld/%02ld:%02ld", p->status.iseed+1, p->status.nseed, p->status.isim+1, p->status.simend, restm, rests, totm, tots);
	}
	if(np>=sizeof(tmp)){
		warning_once("format is truncated\n");
	}
	gtk_list_store_set(listall, &iter,
		COL_STEPT, stmp,
		COL_STEP, tmp,
		COL_STEPP, (gint)(p->frac*100),
		-1);

	snprintf(tmp, 64, "%.1f", p->status.clerrlo);
	gtk_list_store_set(listall, &iter, COL_ERRLO, tmp, -1);
	snprintf(tmp, 64, "%.1f", p->status.clerrhi);
	gtk_list_store_set(listall, &iter, COL_ERRHI, tmp, -1);

}
static void list_proc_reset(proc_t* p){
	GtkTreeIter iter;
	if(!p || !p->row || list_get_iter(p->row, &iter)) return;
	char spid[12];
	snprintf(spid, 12, " %d ", p->pid);

	char sdate[80];
	struct tm* tim=localtime(&p->status.timlast);
	strftime(sdate, 80, "%m-%d %H:%M:%S", tim);
	gtk_list_store_set(listall, &iter,
		COL_DATE, sdate,
		COL_PID, spid,
		//COL_SEED, " ",
		//COL_SEEDP, 0, //Don't use 0.
		COL_STEPT," ",
		COL_STEP, " ",
		COL_STEPP, 0,
		COL_ERRHI, " ",
		COL_ERRLO, " ",
		-1);
	//p->iseed_old=-1;
}
/**
 * Append proc to the list.
*/
static void list_proc_append(proc_t *p){
	if(p->row) return;
	char sdate[80];
	char spid[12];
	snprintf(spid, 12, " %d ", p->pid);
	struct tm* tim=localtime(&p->status.timlast);
	strftime(sdate, 80, "%m-%d %H:%M:%S", tim);
	char* spath=p->path;
	char* sstart=NULL, * sout=NULL, * sargs=NULL;
	if(spath){
		const char* pos=NULL;
		pos=strstr(spath, "/maos ");
		if(!pos){
			pos=strstr(spath, "/skyc ");
		}
		if(pos){
			sstart=(char*)malloc(pos-spath+1);
			memcpy(sstart, spath, pos-spath);
			sstart[pos-spath]='\0';
			sargs=strdup(pos+1);
			char* pos2=NULL;
			for(char* tmp=sargs; (tmp=strstr(tmp, " -o ")); tmp+=4){
				pos2=tmp;
			}
			if(pos2){
				pos2+=4;
				char* pos3=strchr(pos2, ' ');
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
	GtkTreeIter iter;
	gtk_list_store_append(listall, &iter);
	gtk_list_store_set(listall, &iter,
		COL_DATE, sdate,
		COL_TIME, p->status.timstart,
		COL_PID, spid,
		COL_FULL, spath?spath:" ",
		COL_PATH, sstart?sstart:" ",
		COL_ARGS, sargs?sargs:" ",
		COL_OUT, sout?sout:" ",
		COL_ERRHI, " ",
		COL_ERRLO, " ",
		COL_STEPT, " ",
		COL_STEP, " ",
		COL_STEPP, 0,
		COL_HOST, hosts[p->hid],
		-1);
	free(sstart); free(sout); free(sargs);
	GtkTreePath* tpath=gtk_tree_model_get_path(GTK_TREE_MODEL(listall), &iter);
	p->row=gtk_tree_row_reference_new(GTK_TREE_MODEL(listall), tpath);
	gtk_tree_path_free(tpath);
}
//calls by monitor_thread
gboolean remove_entry(GtkTreeRowReference* row){
	GtkTreeIter iter;
	if(!list_get_iter(row, &iter)){
		gtk_list_store_remove(listall, &iter);
	}
	return 0;//must return false
}
gboolean refresh(proc_t* p){
	if(!p) return 0;
	if(!p->row){
		list_proc_append(p);
		list_proc_update(p);
	}
	switch(p->status.info){
	case 0:
		break;
	case S_RUNNING:
		list_proc_update(p);
		list_modify_icon(p->row, icon_running);
		break;
	case S_WAIT: /*waiting to start */
	//list_modify_status(p, "Waiting");
		list_proc_reset(p);
		list_modify_icon(p->row, icon_waiting);
		break;
	case S_START: /*just started. */
	//list_modify_status(p, "Started");
		list_proc_reset(p);
		list_modify_icon(p->row, icon_running);
		notify_user(p);
		break;
	case S_QUEUED:
	//list_modify_status(p, "Queued");
		list_proc_reset(p);
		list_modify_icon(p->row, icon_waiting);
		break;
	case S_UNKNOWN:
		list_modify_icon(p->row, icon_cancel);
		break;
	case S_FINISH:/*Finished */
		list_proc_update(p);
		list_modify_icon(p->row, p->frac==0?icon_skip:icon_finished);
		//list_modify_icon(p->row, icon_finished);
		//list_modify_color(p->row,"#00DD00");
		notify_user(p);
		break;
	case S_CRASH:/*Error */
	//list_modify_status(p, "Error");
		list_modify_icon(p->row, icon_failed);
		//list_modify_color(p->row,"#CC0000");
		notify_user(p);
		break;
	case S_TOKILL:/*kill command sent */
	//list_modify_status(p, "Kill command sent");
		list_modify_icon(p->row, icon_failed);
		//list_modify_color(p->row,"#CCCC00");
		break;
	case S_REMOVE:
		break;
	case S_KILLED:
	//list_modify_status(p, "Killed");
		list_modify_icon(p->row, icon_failed);
		//list_modify_color(p->row,"#CC0000");
		notify_user(p);
		break;
	default:
		warning("Unknown info: %d\n", p->status.info);
	}
	if(p->status.warning){
		list_modify_color(p->row, "#FF0000");
	}
	return 0;//must return false
}

static GtkTreeViewColumn* new_column(int type, int width, const char* title, ...){
	GtkTreeViewColumn* col;
	GtkCellRenderer* render=NULL;
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
	gtk_tree_view_column_set_expand(col, (width<0)?TRUE:FALSE);

	//column only hides text when 1) maxwidth is set, 2) ellipsize is set
	if(width){
		gtk_tree_view_column_set_min_width(col, abs(width));
	}
	if(type==0&&width){//set ellipsize makes it shrinkable
		g_object_set(G_OBJECT(render), "ellipsize", PANGO_ELLIPSIZE_START, NULL);
	}

	gtk_tree_view_column_pack_start(col, render, TRUE);
	va_list ap;
	va_start(ap, title);
	const char* att=NULL;
	while((att=va_arg(ap, const char*))){
		gtk_tree_view_column_add_attribute(col, render, att, va_arg(ap, int));
	}
	va_end(ap);
	return col;
}
#if GTK_MAJOR_VERSION>3	
//to be used with callback
void scheduler_cmd_wrap2(GtkDialog *dialog, int response_id, int* cmds){
	gtk_window_destroy(GTK_WINDOW(dialog));
	int ans;
	extern int *sock_main;
	if(response_id==GTK_RESPONSE_YES){
		stwriteintarr(sock_main[1], cmds, 4);
		streadint(sock_main[1], &ans);
	}
	free(cmds);
}
CHECK_ARG(2)
void dialog_confirm(int *cmds, const char* format, ...){
#else
CHECK_ARG(1)
gboolean dialog_confirm(const char* format, ...){
#endif
	format2fn;
	GtkWidget* dia=gtk_message_dialog_new
	(GTK_WINDOW(window), GTK_DIALOG_DESTROY_WITH_PARENT|GTK_DIALOG_MODAL,
		GTK_MESSAGE_QUESTION,
		GTK_BUTTONS_YES_NO,
		"%s", fn);
	
#if GTK_MAJOR_VERSION>3
	g_signal_connect(GTK_DIALOG(dia), "response", G_CALLBACK(scheduler_cmd_wrap2), cmds);
	gtk_widget_show(dia);
#else
	int result=gtk_dialog_run(GTK_DIALOG(dia));
	gtk_widget_destroy(dia);
	return result==GTK_RESPONSE_YES;
#endif	
}
#if GTK_MAJOR_VERSION<4
/**
   Append text to clipboard (primary and default).
*/
static void clipboard_append(const char* jobinfo){
#if GTK_MAJOR_VERSION > 3 
//GTK4 replaces gtkclipboard by gdkclipboard
	GdkClipboard *clipboard = gtk_widget_get_clipboard(window);
	GdkContentProvider *provider=gdk_clipboard_get_content(clipboard);
	GValue value=G_VALUE_INIT;
	g_value_init(&value, G_TYPE_STRING);
	if(!gdk_content_provider_get_value(provider, &value, NULL)) return;
	const char *old=g_value_get_string(&value);
	char *newer=NULL;
	if(old){
		newer=stradd(old, jobinfo, "\n", NULL);
		g_value_reset(&value);
	}else{
		newer=stradd(jobinfo, "\n", NULL);
	}
	g_value_take_string(&value, newer);
	gdk_clipboard_set_value(clipboard, &value);
	g_value_unset(&value);
#else
	GdkAtom atoms[2]={GDK_SELECTION_CLIPBOARD, GDK_SELECTION_PRIMARY};
	for(int iatom=0; iatom<2; iatom++){
		GtkClipboard* clip=gtk_clipboard_get(atoms[iatom]);
		if(!clip) continue;
		gchar* old=gtk_clipboard_wait_for_text(clip);
		gchar* newer=NULL;
		if(old){
			newer=stradd(old, jobinfo, "\n", NULL);
			g_free(old);
		} else{
			newer=stradd(jobinfo, "\n", NULL);
		}
		gtk_clipboard_set_text(clip, newer, -1);
		free(newer);
	}
#endif
}

typedef struct{
	const char* menu;
	const char* action;
	GdkPixbuf** icon;
	int command;
}menudata_t;
gint cur_page=-1;
menudata_t menudata[]={
	{"Plot selected jobs", "Plot", NULL, CMD_DISPLAY},
	{"Clear selected jobs", "Remove", NULL, CMD_REMOVE},
	{NULL, NULL, NULL, -1},
	{"Kill selected jobs", "Kill", &icon_failed, CMD_KILL},
	{"Restart selected jobs", "Restart", NULL, CMD_RESTART},
	{NULL, NULL, NULL, -1},
	{"Copy cmdline selected jobs", "Copy", NULL, -1},
	{"Copy path of selected jobs", "CopyPath", NULL, -1},
	{"Copy output path of selected jobs", "CopyOutPath", NULL, -1}
};
/*A general routine handle actions to each item.*/
static void handle_selection(GtkTreeModel* model, GtkTreePath* path, GtkTreeIter* iter, gpointer user_data){
	(void)path;
	menudata_t* data=user_data;
	int cmd=data->command;
	const char* action=data->action;
	GValue value=G_VALUE_INIT;
	if(cmd<0){
		if(cmd==-1){
			if(!strcmp(action, "CopyPath")){
				gtk_tree_model_get_value(model, iter, COL_PATH, &value);
			} else if(!strcmp(action, "CopyOutPath")){
				gtk_tree_model_get_value(model, iter, COL_OUT, &value);
			} else if(!strcmp(action, "Copy")){
				gtk_tree_model_get_value(model, iter, COL_FULL, &value);
			}
			const gchar* jobinfo=g_value_get_string(&value);
			clipboard_append(jobinfo);
			g_value_unset(&value);
		}
	} else{
		gtk_tree_model_get_value(model, iter, COL_PID, &value);
		int pid=strtol(g_value_get_string(&value), NULL, 10);
		g_value_unset(&value);
		gtk_tree_model_get_value(model, iter, COL_HOST, &value);
		const char* hostn=g_value_get_string(&value);
		int ihost=host2i(hostn);
		g_value_unset(&value);
		if(scheduler_cmd_wrap(ihost, pid, cmd)){
			warning("Failed to %s the job\n", action);
		}
	}
}

/*Handles menu item clicks in a general way.*/
static void handle_menu_event(GtkMenuItem* menuitem, gpointer user_data){
	(void)menuitem;
	menudata_t* data=user_data;
	GtkWidget* view=views[cur_page];
	GtkTreeSelection* selection=gtk_tree_view_get_selection(GTK_TREE_VIEW(view));
	int nsel=gtk_tree_selection_count_selected_rows(selection);
	if(nsel<1){
		warning("nsel=%d\n", nsel);
		return;
	}
	gboolean ans=1;
	const char* action=data->action;
	/*Alert user in Kill or Restart event*/
	if(!strcmp(action, "Kill")||!strcmp(action, "Restart")){
		ans=dialog_confirm("%s %d jobs?", action, nsel);
		/*GtkWidget* dia=gtk_message_dialog_new
		(GTK_WINDOW(window), GTK_DIALOG_DESTROY_WITH_PARENT,
			GTK_MESSAGE_QUESTION,
			GTK_BUTTONS_YES_NO,
			"%s %d jobs?", action, nsel);
		result=gtk_dialog_run(GTK_DIALOG(dia));
		gtk_widget_destroy(dia);*/
	}
	if(ans){
		gtk_tree_selection_selected_foreach(selection, handle_selection, user_data);
	}
}

static gboolean view_popup_menu(GtkWidget* view, gpointer user_data){
	(void)user_data;
	GtkTreeSelection* selection=gtk_tree_view_get_selection(GTK_TREE_VIEW(view));
	int nsel=gtk_tree_selection_count_selected_rows(selection);

	GtkWidget* menu=gtk_menu_new();
	GtkWidget* menuitem;
	char text[40];
	snprintf(text, 40, "%d selected", nsel);
	menuitem=gtk_menu_item_new_with_label(text);
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);
	menuitem=gtk_separator_menu_item_new();
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);

	if(nsel>0){
		cur_page=GPOINTER_TO_INT(user_data);
		for(size_t i=0; i<sizeof(menudata)/sizeof(menudata_t); i++){
			if(menudata[i].menu){
				if(menudata[i].icon && *menudata[i].icon){
					menuitem=gtk_menu_item_new();
#if GTK_MAJOR_VERSION >=3
					GtkWidget* box=gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
#else
					GtkWidget* box=gtk_hbox_new(0, 6);
#endif
					GtkWidget* icon=gtk_image_new_from_pixbuf(*menudata[i].icon);
					GtkWidget* label=gtk_label_new(menudata[i].menu);
					gtk_container_add(GTK_CONTAINER(box), label);
					gtk_box_pack_end(GTK_BOX(box), icon, TRUE, FALSE, 0);
					//gtk_container_add(GTK_CONTAINER(box), icon);
					gtk_container_add(GTK_CONTAINER(menuitem), box);
				} else{
					menuitem=gtk_menu_item_new_with_label(menudata[i].menu);
					//gtk_widget_set_margin_start(menuitem, 16);
				}
				g_signal_connect(menuitem, "activate", G_CALLBACK(handle_menu_event), menudata+i);
			} else{
				menuitem=gtk_separator_menu_item_new();
			}
			gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);
		}
	}
	gtk_widget_show_all(menu);
#if GTK_MAJOR_VERSION>=3 		
	gtk_menu_popup_at_pointer(GTK_MENU(menu), NULL);
#else
	gtk_menu_popup(GTK_MENU(menu), NULL, NULL, NULL, NULL, 3, gtk_get_current_event_time());
#endif
	return TRUE;
}
#endif
/* Handle click on treeview */
#if GTK_MAJOR_VERSION<4
static gboolean view_click_event(GtkWidget* view, GdkEventButton* event, gpointer user_data){
	gdouble x=event->x;
	gdouble y=event->y;
	int button=event->button;
#else
static gboolean view_click_event(GtkGestureClick* self, gint n_press, gdouble x, gdouble y, gpointer user_data){
	int button=gtk_gesture_single_get_current_button(GTK_GESTURE_SINGLE(self));
	GtkWidget *view=user_data;
	(void)n_press;
#endif

	GtkTreePath* path=NULL;
	GtkTreeViewColumn* column=NULL;
	GtkTreeSelection* selection=gtk_tree_view_get_selection(GTK_TREE_VIEW(view));
	
	if(gtk_tree_view_get_path_at_pos(GTK_TREE_VIEW(view),
		(gint)x,(gint)y,&path, &column, NULL, NULL)){
		if(button==3){//right click, select if nothing is selected. 
			if(gtk_tree_selection_count_selected_rows(selection)==0){
				//the path clicked is not selected. move selection to this path
				gtk_tree_selection_unselect_all(selection);
				gtk_tree_selection_select_path(selection, path);
			}
		}
	} else{//clicks outsize of treeview. unselect all.
		if(button==1){
			gtk_tree_selection_unselect_all(selection);
			column=NULL;
		}
	}
	if(path) {
		gtk_tree_path_free(path); path=NULL;
	}
	//right click popup menu
	if(button==3){
#if GTK_MAJOR_VERSION<4	
		view_popup_menu(view, user_data);
#endif
		return TRUE;//return TRUE to avoid modifying selection
	}
	return FALSE;//default handler handles multi selection
}
/* Handle click on icon column to handle individual jobs*/
#if GTK_MAJOR_VERSION<4
static gboolean view_release_event(GtkWidget* view, GdkEventButton* event, gpointer user_data){
	int button=event->button;
	gdouble x=event->x;
	gdouble y=event->y;
	(void) user_data;
#else
static gboolean view_release_event(GtkGestureClick*self, gint n_press, gdouble x, gdouble y, gpointer user_data){
	int button=gtk_gesture_single_get_current_button(GTK_GESTURE_SINGLE(self));
	GtkWidget* view=user_data;
	(void)n_press;
#endif
	;
	if(button!=1) return FALSE;//only handle left click
	GtkTreePath* path=NULL;
	GtkTreeViewColumn* column=NULL;
	if(gtk_tree_view_get_path_at_pos(GTK_TREE_VIEW(view),
		(gint)x, (gint)y, &path, &column, NULL, NULL)){
		const gchar* title=gtk_tree_view_column_get_title(column);
		if(!strcmp(title, " ")){
			GtkTreeIter iter;
			//gint ihost=GPOINTER_TO_INT(user_data);
			GtkTreeModel* model=gtk_tree_view_get_model(GTK_TREE_VIEW(view));
			//GtkTreeModel* model=GTK_TREE_MODEL(lists[ihost]);
			gtk_tree_model_get_iter(model, &iter, path);
			GValue value=G_VALUE_INIT;
			gtk_tree_model_get_value(model, &iter, COL_PID, &value);
			int pid=strtol(g_value_get_string(&value), NULL, 10);
			g_value_unset(&value);
			gtk_tree_model_get_value(model, &iter, COL_HOST, &value);
			const char* hostn=g_value_get_string(&value);
			int ihost=host2i(hostn);
			g_value_unset(&value);
			GdkPixbuf*status2=0;
			gtk_tree_model_get(model, &iter, COL_ACTION, &status2, -1);
			
			if(status2==icon_running||status2==icon_waiting){
#if GTK_MAJOR_VERSION>3				
				int* cmds=malloc(sizeof(int)*4);
				cmds[0]=MON_CMD; cmds[1]=ihost; cmds[2]=pid; cmds[3]=CMD_KILL;
				dialog_confirm(cmds, "Kill %d?", pid);
#else				
				if(dialog_confirm("Kill %d?", pid)){
					
					scheduler_cmd_wrap(ihost, pid, CMD_KILL);
				}
#endif
			}else{
				scheduler_cmd_wrap(ihost, pid, CMD_REMOVE);
			}
		}
	}
	return TRUE;
}
static void concat_selected_path(GtkTreeModel* model, GtkTreePath* path, GtkTreeIter* iter, gpointer user_data){
	(void)path;
	gchar** buf=(gchar**)user_data;
	if(!*buf){
		GValue value=G_VALUE_INIT;
		gtk_tree_model_get_value(model, iter, COL_FULL, &value);
		const gchar* val=g_value_get_string(&value);
		*buf=g_strdup(val);
		g_value_unset(&value);
	} else{
		//*buf=g_strjoin("\n", *buf, val, NULL);
	}
}
/*
  handle section change event to update the GtkTextBuffer buffers[ihost]
 */
void view_selection_event(GtkTreeSelection* selection, gpointer user_data){
	int ipage=GPOINTER_TO_INT(user_data);
	int nsel=gtk_tree_selection_count_selected_rows(selection);
	if(nsel){
		gchar* buf=0;
		gtk_tree_selection_selected_foreach(selection, concat_selected_path, &buf);
		gtk_text_buffer_set_text(buffers[ipage], buf, -1);
		g_free(buf);
	}
}
static gboolean filter_host(GtkTreeModel* model, GtkTreeIter* iter, gpointer host){
	gchar* host2;
	gtk_tree_model_get(model, iter, COL_HOST, &host2, -1);
	int ans=host2&&(!strcmp(host2, (gchar*)host));
	g_free(host2);
	return ans;
}
static gboolean filter_status(GtkTreeModel* model, GtkTreeIter* iter, gpointer status){
	GdkPixbuf* status2=0;
	(void)status;
	gtk_tree_model_get(model, iter, COL_ACTION, &status2, -1);
	return (status2==icon_running || status2==icon_finished);
}
GtkWidget* new_page(int ihost){
	if(!listall){
		listall=gtk_list_store_new(COL_TOT,
			G_TYPE_STRING,/*DATE */
			G_TYPE_LONG, /*TIME*/
			G_TYPE_STRING,/*PID */
			G_TYPE_STRING,/*FULL */
			G_TYPE_STRING,/*PATH */
			G_TYPE_STRING,/*ARGS */
			G_TYPE_STRING,/*OUT */
			//G_TYPE_STRING,/*SEED */
			//G_TYPE_INT, /*SEEDP */
			G_TYPE_STRING,/*STEPT*/
			G_TYPE_STRING,/*STEP */
			G_TYPE_INT, /*STEPP */
			G_TYPE_STRING,/*TIMING */
			G_TYPE_STRING,/*ALL */
			G_TYPE_STRING,/*REST */
			G_TYPE_STRING,/*ERRLO */
			G_TYPE_STRING,/*ERRHI */
			GDK_TYPE_PIXBUF,/*ACTION */
			G_TYPE_STRING, /*COLOR */
			G_TYPE_STRING /*HOST*/
		);
		lists=mycalloc(nhost+1, GtkTreeModel*);
		views=mycalloc(nhost+1, GtkWidget*);
		gtk_tree_sortable_set_sort_column_id(GTK_TREE_SORTABLE(listall), COL_TIME, GTK_SORT_ASCENDING);
	}
	if(views[ihost]){
		return views[ihost];
	}
	lists[ihost]=gtk_tree_model_filter_new(GTK_TREE_MODEL(listall), NULL);
	if(ihost<nhost){
		gtk_tree_model_filter_set_visible_func(GTK_TREE_MODEL_FILTER(lists[ihost]), filter_host, hosts[ihost], NULL);
	} else{
		gtk_tree_model_filter_set_visible_func(GTK_TREE_MODEL_FILTER(lists[ihost]), filter_status, NULL, NULL);
	}
	
	GtkWidget* view;
	views[ihost]=view=gtk_tree_view_new_with_model(GTK_TREE_MODEL(lists[ihost]));
	//g_object_unref(lists[ihost]);
	{
#if GTK_MAJOR_VERSION < 4			
		g_signal_connect(view, "button-press-event",
			G_CALLBACK(view_click_event), GINT_TO_POINTER(ihost));
		g_signal_connect(view, "button-release-event",
			G_CALLBACK(view_release_event), GINT_TO_POINTER(ihost));

		g_signal_connect(view, "popup-menu",
			G_CALLBACK(view_popup_menu), GINT_TO_POINTER(ihost));
#else
		GtkGesture* gesture=gtk_gesture_click_new();
		g_signal_connect(gesture, "pressed", 
			G_CALLBACK(view_click_event), view);
		g_signal_connect(gesture, "released",
		G_CALLBACK(view_release_event), view);
		gtk_widget_add_controller(view, GTK_EVENT_CONTROLLER(gesture));
#endif
	}
	GtkTreeSelection* viewsel=gtk_tree_view_get_selection(GTK_TREE_VIEW(view));
	/*gtk_tree_selection_set_select_function(viewsel, treeselfun,NULL,NULL); */
	gtk_tree_selection_set_mode(viewsel, GTK_SELECTION_MULTIPLE);
	g_signal_connect(viewsel, "changed", G_CALLBACK(view_selection_event), GINT_TO_POINTER(ihost));
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
	//gtk_tree_view_set_rules_hint(GTK_TREE_VIEW(view), TRUE);
	//gtk_tree_view_set_headers_clickable(GTK_TREE_VIEW(view), TRUE);
	gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(2, 0, " ", "pixbuf", COL_ACTION, NULL));
	gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 0, "Date", "text", COL_DATE, NULL));
	if(ihost==nhost){
		gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 0, "Host", "text", COL_HOST, NULL));
	} else{
		gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 0, "PID", "text", COL_PID, NULL));
	}
	//gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 50, "Path", "text", COL_PATH, NULL));
	//gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 50, "Args", "text", COL_ARGS, NULL));
	gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, -100, "Out", "text", COL_OUT, NULL));
	gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 0, "Low", "text", COL_ERRLO, "foreground", COL_COLOR, NULL));
	gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 0, "High", "text", COL_ERRHI, "foreground", COL_COLOR, NULL));
	//gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(1, 0, "Seed", "text", COL_SEED, "value", COL_SEEDP, NULL));
	gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 0, "Step", "text", COL_STEPT, NULL));
	gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(1, 0, "Progress", "text", COL_STEP, "value", COL_STEPP, NULL));
	
	gtk_tree_view_columns_autosize(GTK_TREE_VIEW(view));
	return view;
}


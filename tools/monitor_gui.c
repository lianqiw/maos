/*
  Copyright 2009-2026 Lianqi Wang
  
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
#include <ctype.h>
#include <glib.h>
#include <gtk/gtk.h>
#include <gdk/gdk.h>
#include <glib/gprintf.h>
#include "monitor.h"
/*DO not modify enum without modify list store*/
enum{
	COL_DATE=0,//date in string
	COL_TIME,//time in time_t. //only set once, to keep task order stable.
	COL_HOST,
	COL_PID,
	COL_FULL,/*full path+arguments+output*/
	COL_TOOLTIP,/*reformatted from COL_FULL for tooltip*/
	COL_PATH,/*starting path*/
	COL_ARGS,/*arguments*/
	COL_OUT, /*output folder*/
	//COL_SEED,
	//COL_SEEDP,/*gint */
	COL_STEPT,
	COL_STEP,
	COL_STEPP,/*gint */
	COL_TIMING,
	//COL_ALL,
	//COL_REST,
	COL_ERRLO,
	COL_ERRHI,
	COL_ACTION,
	COL_COLOR,
	
	COL_TOT,
};
#define SHOW_SOUT 1 //Show output path
#define SHOW_STEPT 1 //Show step time separately
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
static void list_modify_icon(proc_t *proc, GdkPixbuf* newicon){
	if(!proc || !newicon || proc->oldinfo==proc->status.info) return;
	GtkTreeIter iter;
	if(!list_get_iter(proc->row, &iter)){
		gtk_list_store_set(listall, &iter, COL_ACTION, newicon, -1);
	}
}
static void list_modify_color(proc_t *proc, const char *color){
	if(!proc||!color||proc->oldinfo==proc->status.info) return;
	GtkTreeIter iter;
	if(!list_get_iter(proc->row, &iter)){
		gtk_list_store_set(listall, &iter, COL_COLOR, color, -1);
	}
}

static void list_proc_update(proc_t* p){
	if(!p || p->status.nseed==0 || p->status.simend==0) return;//skipped simulation
	GtkTreeIter iter;
	if(list_get_iter(p->row, &iter)) return;//failure
	status_t *ps=&p->status;
	double total=(double)(ps->rest+ps->laps);
	if(fabs(total)>1.e-10){
		p->frac=MIN(1, (double)ps->laps/total);
	} else{
		p->frac=1;//finish too soon.
	}
	const double step=ps->tot*ps->scale;//scale is normally 1
	/*if(ps->iseed!=p->iseed_old){
		char tmp[64];
		snprintf(tmp, 64, "%d/%d", ps->iseed+1, ps->nseed);
		gtk_list_store_set(listall, &iter,
			COL_SEED, tmp,
			COL_SEEDP, (100*(ps->iseed+1)/ps->nseed),
			-1);
		p->iseed_old=ps->iseed;
	}*/
	char stmp[8];
	sec2str(stmp, sizeof(stmp), step);

	char tmp[64];
	snprintf(tmp, 64, "%d/%d %d/%d ", ps->iseed+1, ps->nseed, ps->isim+1, ps->simend);
	unsigned int offset=strlen(tmp);
	
	offset+=sec2str(tmp+offset, sizeof(tmp)-offset, ps->rest);
	if(offset+1<sizeof(tmp)){
		tmp[offset]='/';offset++;
	}
	offset+=sec2str(tmp+offset, sizeof(tmp)-offset, ps->rest+ps->laps);
#if SHOW_STEPT==0	
	if(offset<sizeof(tmp)){
		tmp[offset]=' ';offset++;
	}
	offset+=sec2str(tmp+offset, sizeof(tmp)-offset, step);
#endif
	gtk_list_store_set(listall, &iter,
		COL_STEPT, stmp,
		COL_STEP, tmp,
		COL_STEPP, (gint)(p->frac*100),
		-1);
#define ERR_MAX 9999.9
	snprintf(tmp, 64, "%.1f", MIN(ps->clerrlo, ERR_MAX));
	gtk_list_store_set(listall, &iter, COL_ERRLO, tmp, -1);
	snprintf(tmp, 64, "%.1f", MIN(ps->clerrhi, ERR_MAX));
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
	char* stooltip=NULL;//for tooltip
	char* sstart=NULL;//starting directory
	char* sout=NULL;  //output directory
	char* sargs=NULL; //job arguments
	if(spath){
		const char* pos=NULL;
		pos=strstr(spath, "/maos ");
		if(!pos){
			pos=strstr(spath, "/skyc ");
		}
		if(pos){
			sstart=(char*)malloc(pos-spath+1);//job start directory
			memcpy(sstart, spath, pos-spath);
			sstart[pos-spath]='\0';
			sargs=strdup(pos+1);//job arguments. 
			//Move output directory to the end of sargs
			char* tmp=sargs;
			while((tmp=strstr(tmp, " -o"))){
				char* pos2=tmp+3;//start of output directory
				while(isspace(pos2[0])) pos2++;//skip space
				char* pos3=strchr(pos2, ' ');//end of output directory
				if(!pos3) pos3=strchr(pos2, 0);//end of string
				if(pos3){
					if(sout) free(sout);
					sout=(char*)malloc(pos3-pos2+1);
					memcpy(sout, pos2, pos3-pos2);
					sout[pos3-pos2]='\0';
					memmove(tmp, pos3, strlen(pos3)+1);
				}
			}
#if SHOW_SOUT == 0			
			if(sout){
				if(strstr(spath, " -o ")){//there is room for space after -o
					strcat(sargs, " -o ");
				}else{
					strcat(sargs, " -o");
				}
				strcat(sargs, sout);//append sout back to sargs. Cannot overflow
			}
#endif			
		}
		{
			stooltip=strdup(sargs);
			parse_argopt(stooltip, NULL);
		}
	}
	GtkTreeIter iter;
	gtk_list_store_append(listall, &iter);
	gtk_list_store_set(listall, &iter,
		COL_DATE, sdate,
		COL_TIME, p->status.timstart,
		COL_HOST, hostshort[p->hid],
		COL_PID, spid,
		COL_FULL, spath?spath:" ",
		COL_TOOLTIP, stooltip?stooltip:"",
		COL_PATH, sstart?sstart:" ",
		COL_ARGS, sargs?sargs:" ",//arguments
		COL_OUT, sout?sout:" ",//output directory
		COL_ERRHI, " ",
		COL_ERRLO, " ",
		COL_STEPT, " ",
		COL_STEP, " ",
		COL_STEPP, 0,
		-1);
	free(sstart); free(sout); free(sargs); free(stooltip);
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
		list_modify_icon(p, icon_running);
		break;
	case S_WAIT: /*waiting to start */
	//list_modify_status(p, "Waiting");
		list_proc_reset(p);
		list_modify_icon(p, icon_waiting);
		break;
	case S_START: /*just started. */
	//list_modify_status(p, "Started");
		list_proc_reset(p);
		list_modify_icon(p, icon_running);
		notify_user(p);
		break;
	case S_QUEUED:
	//list_modify_status(p, "Queued");
		list_proc_reset(p);
		list_modify_icon(p, icon_waiting);
		break;
	case S_UNKNOWN:
		list_modify_icon(p, icon_cancel);
		break;
	case S_FINISH:/*Finished */
		list_proc_update(p);
		list_modify_icon(p, p->frac==0?icon_skip:icon_finished);
		//list_modify_icon(p, icon_finished);
		//list_modify_color(p,"#00DD00");
		notify_user(p);
		break;
	case S_CRASH:/*Error */
	//list_modify_status(p, "Error");
		list_modify_icon(p, icon_failed);
		//list_modify_color(p,"#CC0000");
		notify_user(p);
		break;
	case S_TOKILL:/*kill command sent */
	//list_modify_status(p, "Kill command sent");
		list_modify_icon(p, icon_failed);
		//list_modify_color(p,"#CCCC00");
		break;
	case S_REMOVE:
		break;
	case S_KILLED:
	//list_modify_status(p, "Killed");
		list_modify_icon(p, icon_failed);
		//list_modify_color(p,"#CC0000");
		notify_user(p);
		break;
	default:
		warning("Unknown info: %d\n", p->status.info);
	}
	if(p->status.warning){
		list_modify_color(p, "#FF0000");
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
		break;
	case 1:
		render=gtk_cell_renderer_progress_new();
		g_object_set(G_OBJECT(render), "text-xalign", 0.0, NULL);
		break;
	case 2:
		render=gtk_cell_renderer_pixbuf_new();
		break;
	default:
		error("Invalid\n");
	}
#if GTK_VERSION_AFTER(2, 18)
	gtk_cell_renderer_set_padding(render, 1, 1);
	gtk_cell_renderer_set_alignment(render, 1, 0.5);
#endif
	gtk_tree_view_column_set_title(col, title);
	gtk_tree_view_column_set_spacing(col, 2);
	gtk_tree_view_column_set_alignment(col, 1);
	//column only hides text when 1) maxwidth is set, 2) ellipsize is set
	if(width){
		gtk_tree_view_column_set_min_width(col, abs(width));
		if(type==0 && width<0){//set ellipsize makes it shrinkable than minimal needed size.
			g_object_set(G_OBJECT(render), "ellipsize", PANGO_ELLIPSIZE_START, NULL);
		}
	}
	const int expandable=(width<0)?TRUE:FALSE;
	const int resizable=width?TRUE:FALSE;//resizeable makes the column very small if not expand.

	gtk_tree_view_column_set_resizable(col, resizable);//allow user to grab and resize
	gtk_tree_view_column_set_expand(col, expandable);//cell is expandable within column
	gtk_tree_view_column_pack_start(col, render, expandable);//column is expandable
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
static void dialog_confirm(int *cmds, const char* format, ...){
#else
CHECK_ARG(1)
static gboolean dialog_confirm(const char* format, ...){
#endif
	format2fn;
	GtkWidget* dia=gtk_message_dialog_new
	(GTK_WINDOW(window), (GtkDialogFlags)(GTK_DIALOG_DESTROY_WITH_PARENT|GTK_DIALOG_MODAL),
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
int reset_clipboard=0;
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
		gchar* old=reset_clipboard?NULL:gtk_clipboard_wait_for_text(clip);
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
	if(reset_clipboard){
		reset_clipboard=0;
	}
#endif
}

typedef struct{
	const char* menu;
	const char* action;
	GdkPixbuf** icon;//use address for compile time constant
	int command;
}menudata_t;
gint cur_page=-1;
menudata_t menudata[]={
	{"Plot", "Plot", &icon_draw, CMD_DRAWSER},
	{"Clear", "Clear", &icon_clear, CMD_REMOVE},
	{NULL, NULL, NULL, -1},
	{"Kill", "Kill", &icon_failed, CMD_KILL},
	{"Restart", "Restart", &icon_running, CMD_RESTART},
	{NULL, NULL, NULL, -1},
	{"Copy command line", "Copy", NULL, -1},
	{"Copy starting path", "CopyPath", NULL, -1},
	{"Copy output path", "CopyOutPath", NULL, -1}
};
/*A general routine handle actions to each item.*/
static void handle_selection(GtkTreeModel* model, GtkTreePath* path, GtkTreeIter* iter, gpointer user_data){
	(void)path;
	menudata_t* data=(menudata_t*) user_data;
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
		scheduler_cmd_wrap(ihost, pid, cmd);
	}
}

/*Handles menu item clicks in a general way.*/
static void handle_menu_event(GtkMenuItem* menuitem, gpointer user_data){
	(void)menuitem;
	menudata_t* data=(menudata_t*)user_data;
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
		reset_clipboard=1;
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
	snprintf(text, 40, "%d jobs selected:", nsel);
	menuitem=gtk_menu_item_new_with_label(text);
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);
	menuitem=gtk_separator_menu_item_new();
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);

	if(nsel>0){
		cur_page=GPOINTER_TO_INT(user_data);
		for(size_t i=0; i<sizeof(menudata)/sizeof(menudata_t); i++){
			if(menudata[i].menu){
				if(1){
					menuitem=gtk_menu_item_new();
#if GTK_MAJOR_VERSION >=3
					GtkWidget* box=gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
#else
					GtkWidget* box=gtk_hbox_new(0, 6);
#endif
					GtkWidget* label=gtk_label_new(menudata[i].menu);
					gtk_container_add(GTK_CONTAINER(box), label);
					if(menudata[i].icon&&*menudata[i].icon){
						GtkWidget *icon=gtk_image_new_from_pixbuf(*menudata[i].icon);
						gtk_box_pack_end(GTK_BOX(box), icon, FALSE, FALSE, 0);
					}
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
	GtkWidget *view=GTK_WIDGET(user_data);
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
	}/* else{//clicks outsize of treeview. unselect all.
		if(button==1){
			gtk_tree_selection_unselect_all(selection);
			column=NULL;
		}
	}*/
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
	GtkWidget* view=GTK_WIDGET(user_data);
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
				int* cmds=mymalloc(4, int);
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
/*static void concat_selected_path(GtkTreeModel* model, GtkTreePath* path, GtkTreeIter* iter, gpointer user_data){
	(void)path;
	gchar** buf=(gchar**)user_data;
	if(!*buf){
		GValue value=G_VALUE_INIT;
		gtk_tree_model_get_value(model, iter, COL_FULL, &value);
		const gchar* val=g_value_get_string(&value);
		*buf=g_strdup(val);
		parse_argopt(*buf, NULL);
		g_value_unset(&value);
	} else{
		// *buf=g_strjoin("\n", *buf, val, NULL);
	}
}*/
/*
  handle section change event to update the GtkTextBuffer buffers[ihost]
 */
/*void view_selection_event(GtkTreeSelection* selection, gpointer user_data){
	(void)user_data;
	//int ipage=GPOINTER_TO_INT(user_data);
	int nsel=gtk_tree_selection_count_selected_rows(selection);
	if(nsel){
		gchar* buf=0;
		gtk_tree_selection_selected_foreach(selection, concat_selected_path, &buf);
		gtk_text_buffer_set_text(textbuffer, buf, -1);
		gtk_widget_show(textscroll);
		g_free(buf);
	}else{
		gtk_widget_hide(textscroll);
	}
}
static gboolean view_unselect_all(GtkTreeView *view, gpointer user_data){
	(void)view;
	(void)user_data;
	gtk_widget_hide(textscroll);
	return 0;
}*/
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
	return (status2==icon_running||status2==icon_finished||status2==icon_waiting);
}
GtkWidget* new_page(int ihost){
	if(!listall){
		listall=gtk_list_store_new(COL_TOT,
			G_TYPE_STRING,/*COL_DATE */
			G_TYPE_LONG, /*COL_TIME*/
			G_TYPE_STRING, /*COL_HOST*/
			G_TYPE_STRING,/*COL_PID */
			G_TYPE_STRING,/*COL_FULL */
			G_TYPE_STRING,/*COL_TOOLTIP */
			G_TYPE_STRING,/*COL_PATH */
			G_TYPE_STRING,/*COL_ARGS */
			G_TYPE_STRING,/*COL_OUT */
			//G_TYPE_STRING,/*COL_SEED */
			//G_TYPE_INT, /*COL_SEEDP */
			G_TYPE_STRING,/*COL_STEPT*/
			G_TYPE_STRING,/*COL_STEP */
			G_TYPE_INT, /*COL_STEPP */
			G_TYPE_STRING,/*COL_TIMING */
			//G_TYPE_STRING,/*COL_ALL */
			//G_TYPE_STRING,/*COL_REST */
			G_TYPE_STRING,/*COL_ERRLO */
			G_TYPE_STRING,/*COL_ERRHI */
			GDK_TYPE_PIXBUF,/*COL_ACTION */
			G_TYPE_STRING /*COL_COLOR */

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
		gtk_tree_model_filter_set_visible_func(GTK_TREE_MODEL_FILTER(lists[ihost]), filter_host, hostshort[ihost], NULL);
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
		//g_signal_connect(view, "unselect-all", G_CALLBACK(view_unselect_all), GINT_TO_POINTER(ihost));
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
	//g_signal_connect(viewsel, "changed", G_CALLBACK(view_selection_event), GINT_TO_POINTER(ihost));
	gtk_tree_view_set_tooltip_column(GTK_TREE_VIEW(view), COL_TOOLTIP);
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

	gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 0, "Date", "text", COL_DATE, NULL));
	if(ihost==nhost){
		gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 0, "Host", "text", COL_HOST, NULL));
	} else{
		gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 0, "PID", "text", COL_PID, NULL));
	}
	gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 100, "Start Dir", "text", COL_PATH, NULL));
	gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, -100, "Arguments", "text", COL_ARGS, NULL));
#if SHOW_SOUT 
	gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 50, "Out Dir", "text", COL_OUT, NULL));
#endif
	gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 0, "Low", "text", COL_ERRLO, "foreground", COL_COLOR, NULL));
	gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 0, "High", "text", COL_ERRHI, "foreground", COL_COLOR, NULL));
	//gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(1, 0, "Seed", "text", COL_SEED, "value", COL_SEEDP, NULL));
#if SHOW_STEPT	
	gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(0, 0, "Step", "text", COL_STEPT, NULL));
#endif	
	gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(2, 0, " ", "pixbuf", COL_ACTION, NULL));
	gtk_tree_view_append_column(GTK_TREE_VIEW(view), new_column(1, 0, "Progress", "text", COL_STEP, "value", COL_STEPP, NULL));
	
	gtk_tree_view_columns_autosize(GTK_TREE_VIEW(view));
	return view;
}


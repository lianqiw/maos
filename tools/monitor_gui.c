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
	COL_STEPT,/*str; step progressing and timing*/
	COL_STEP, /*str; step timing*/
	COL_STEPP,/*gint; percentage done */
	COL_TIMING,//remaining time if running or total time
	//COL_ALL,
	//COL_REST,
	COL_ERRLO,
	COL_ERRHI,
	COL_ACTION,
	COL_COLOR,
	
	COL_TOT,
};
#define SHOW_START 0 //Show start path separately
#define SHOW_SOUT 0 //Show output path separately
#define SHOW_STEPT 1 //Show step time separately
#if USE_COLUMN_VIEW
#define MON_TYPE_JOB (mon_job_get_type())
G_DECLARE_FINAL_TYPE(MonJob, mon_job, MON, JOB, GObject)
/**
   A row of the job list. It holds the same columns as the GtkListStore used
   with GtkTreeView. Emits "changed" when any column is modified.
*/
struct _MonJob{
	GObject parent;
	char* str[COL_TOT];//string columns. NULL for the others.
	long time;//COL_TIME
	int stepp;//COL_STEPP
	GdkPixbuf* icon;//COL_ACTION. Points to a global icon. Not referenced.
};
G_DEFINE_FINAL_TYPE(MonJob, mon_job, G_TYPE_OBJECT)
static guint job_changed_signal=0;
static GListStore* listall=NULL;
static GtkFilter** filters=NULL;
static GtkWidget** views=NULL;
/**
   Free the string columns of a job row.

   @param[in] obj The job row.
*/
static void mon_job_finalize(GObject* obj){
	MonJob* job=MON_JOB(obj);
	for(int i=0; i<COL_TOT; i++){
		g_free(job->str[i]);
	}
	G_OBJECT_CLASS(mon_job_parent_class)->finalize(obj);
}
/**
   Initialize the MonJob class and register its "changed" signal.

   @param[in] klass The class.
*/
static void mon_job_class_init(MonJobClass* klass){
	G_OBJECT_CLASS(klass)->finalize=mon_job_finalize;
	job_changed_signal=g_signal_new("changed", G_TYPE_FROM_CLASS(klass), G_SIGNAL_RUN_LAST, 0, NULL, NULL, NULL, G_TYPE_NONE, 0);
}
/**
   Initialize a job row. All fields are already zero.

   @param[in] job Unused.
*/
static void mon_job_init(MonJob* job){
	(void)job;
}
/**
   Set columns of a job row and notify the views. Refilters the "All" page if
   the status icon changes.

   @param[in,out] job The job row. May be NULL.
   @param[in] ... Pairs of column index and value, terminated by -1, as with
   gtk_list_store_set(). COL_TIME takes long, COL_STEPP int, COL_ACTION
   GdkPixbuf*, and the others strings.
   @return 0 on success, 1 if the row is invalid.
*/
static int list_set(job_row_t* job, ...){
	if(!job) return 1;
	int icon_changed=0;
	int col;
	va_list ap;
	va_start(ap, job);
	while((col=va_arg(ap, int))!=-1){
		switch(col){
		case COL_TIME:
			job->time=va_arg(ap, long);
			break;
		case COL_STEPP:
			job->stepp=va_arg(ap, int);
			break;
		case COL_ACTION:{
			GdkPixbuf* icon=va_arg(ap, GdkPixbuf*);
			if(job->icon!=icon){
				job->icon=icon;
				icon_changed=1;
			}
		}
			break;
		default:
			g_free(job->str[col]);
			job->str[col]=g_strdup(va_arg(ap, const char*));
		}
	}
	va_end(ap);
	g_signal_emit(job, job_changed_signal, 0);
	if(icon_changed&&filters&&filters[nhost]){
		gtk_filter_changed(filters[nhost], GTK_FILTER_CHANGE_DIFFERENT);
	}
	return 0;
}
#else
static GtkListStore* listall=NULL;
static GtkTreeModel** lists=NULL;
static GtkWidget** views=NULL;
/**
   Get the iterator in the list store of a row reference.

   @param[in] row The row reference. May be NULL.
   @param[out] iter The iterator.
   @return 0 on success, 1 if the row is invalid.
*/
static int list_get_iter(GtkTreeRowReference* row, GtkTreeIter* iter){
	GtkTreePath* tpath=NULL;
	if(row && (tpath=gtk_tree_row_reference_get_path(row))){
		gtk_tree_model_get_iter(GTK_TREE_MODEL(listall), iter, tpath);
		gtk_tree_path_free(tpath);
		return 0;
	}
	return 1;
}
/**
   Set columns of a job row.

   @param[in] row The row reference. May be NULL.
   @param[in] ... Pairs of column index and value, terminated by -1, as with
   gtk_list_store_set().
   @return 0 on success, 1 if the row is invalid.
*/
static int list_set(job_row_t* row, ...){
	GtkTreeIter iter;
	if(list_get_iter(row, &iter)) return 1;
	va_list ap;
	va_start(ap, row);
	gtk_list_store_set_valist(listall, &iter, ap);
	va_end(ap);
	return 0;
}
#endif
/**
   Set the status icon of the job if its status has changed.

   @param[in] proc The job.
   @param[in] newicon The new icon.
*/
static void list_modify_icon(proc_t *proc, GdkPixbuf* newicon){
	if(!proc || !newicon || proc->oldinfo==proc->status.info) return;
	list_set(proc->row, COL_ACTION, newicon, -1);
}
/**
   Set the text color of the job if its status has changed.

   @param[in] proc The job.
   @param[in] color The color string, e.g. "#FF0000".
*/
static void list_modify_color(proc_t *proc, const char *color){
	if(!proc||!color||proc->oldinfo==proc->status.info) return;
	list_set(proc->row, COL_COLOR, color, -1);
}

/**
   Update the progress, step time, remaining or total time, and errors of the job
   in the list. Also updates p->frac.

   @param[in,out] p The job.
*/
static void list_proc_update(proc_t* p){
	if(!p || p->status.nseed==0 || p->status.simend==0) return;//skipped simulation
	if(!p->row) return;//failure
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
		list_set(p->row,
			COL_SEED, tmp,
			COL_SEEDP, (100*(ps->iseed+1)/ps->nseed),
			-1);
		p->iseed_old=ps->iseed;
	}*/
	char stmp[8];
	sec2str(stmp, sizeof(stmp), step);
	char tmp[64];
	snprintf(tmp, 64, "%d/%d %d/%d ", ps->iseed+1, ps->nseed, ps->isim+1, ps->simend);

	char timing[64]; timing[0]=0;
	if(ps->info==S_RUNNING){
		sec2str(timing, sizeof(timing), ps->rest);
	}else if(ps->info==S_FINISH){
		sec2str(timing, sizeof(timing), ps->laps);
	}

	list_set(p->row,
		COL_STEPT, stmp,
		COL_STEP, tmp,
		COL_TIMING, timing,
		COL_STEPP, (gint)(p->frac*100),
		-1);
#define ERR_MAX 9999.9
	snprintf(tmp, 64, "%.1f", MIN(ps->clerrlo, ERR_MAX));
	list_set(p->row, COL_ERRLO, tmp, -1);
	snprintf(tmp, 64, "%.1f", MIN(ps->clerrhi, ERR_MAX));
	list_set(p->row, COL_ERRHI, tmp, -1);

}
/**
   Reset the progress and error columns of the job and update its date and PID.

   @param[in] p The job.
*/
static void list_proc_reset(proc_t* p){
	if(!p || !p->row) return;
	char spid[12];
	snprintf(spid, 12, " %d ", p->pid);

	char sdate[80];
	struct tm* tim=localtime(&p->status.timlast);
	strftime(sdate, 80, "%m-%d %H:%M:%S", tim);
	list_set(p->row,
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
   Append proc to the list. Does nothing if it already has a row.

   @param[in,out] p The job. p->row is set and the output directory is moved to the end of p->path.
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
		char* pos=strstr(spath, "/maos ");
		if(!pos){
			pos=strstr(spath, "/skyc ");
		}
		if(pos){
			sstart=spath;//start of path
#if SHOW_START
			pos[0]=0;\\terminate string
			sargs=pos+1;//start of arguments
#else
			sargs=sstart;
#endif
			int has_space=0;
			//Extract output directory 
			char* tmp=pos+5;
			while((tmp=strstr(tmp, " -o"))){
				char* pos2=tmp+3;//start of output directory
				while(isspace((unsigned char)pos2[0])){
					has_space=1;
					pos2++;//skip space
				}
				char* pos3=strchr(pos2, ' ');//end of output directory
				if(!pos3) pos3=strchr(pos2, 0);//end of string
				if(pos3){
					if(sout) free(sout);
					sout=(char*)malloc((pos3-pos2)+1);
					memcpy(sout, pos2, pos3-pos2);
					sout[pos3-pos2]='\0';
					memmove(tmp, pos3, strlen(pos3)+1);
				}
			}
			{
				stooltip=strdup(pos+6);
				char *dc=strstr(stooltip, "-c ");
				while(dc){//remove space between -c and the following conf to avoid breaking them up
					dc[0]=' '; dc[1]='-'; dc[2]='c';
					dc=strstr(stooltip, "-c ");
				}
				parse_argopt(stooltip, NULL);
			}
#if SHOW_SOUT == 0			
			if(sout){
				int len=strlen(sout);
				memmove(pos+5+len+3+has_space, pos+5, strlen(pos+5)+1);
				memcpy(pos+5, " -o ", 4); 
				memcpy(pos+8+has_space, sout, len);
			}
#endif			
		}
	}
#if USE_COLUMN_VIEW
	p->row=(job_row_t*)g_object_new(MON_TYPE_JOB, NULL);
#else
	GtkTreeIter iter;
	gtk_list_store_append(listall, &iter);
	GtkTreePath* tpath=gtk_tree_model_get_path(GTK_TREE_MODEL(listall), &iter);
	p->row=gtk_tree_row_reference_new(GTK_TREE_MODEL(listall), tpath);
	gtk_tree_path_free(tpath);
#endif
	list_set(p->row,
		COL_DATE, sdate,
		COL_TIME, (long)p->status.timstart,
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
	free(sout); free(stooltip);
#if USE_COLUMN_VIEW
	g_list_store_append(listall, p->row);//p->row keeps its own reference
#endif
}
/**
   Remove the job at row from the job list. Called through g_idle_add() by the
   monitor_thread so that the GUI is modified in the GUI thread.

   @param[in] row Row of the job. Freed.
   @return 0 so that it runs only once.
*/
gboolean remove_entry(job_row_t* row){
	if(!row) return 0;
#if USE_COLUMN_VIEW
	guint pos;
	if(listall&&g_list_store_find(listall, row, &pos)){
		g_list_store_remove(listall, pos);
	}
	g_object_unref(row);
#else
	GtkTreeIter iter;
	if(!list_get_iter(row, &iter)){
		gtk_list_store_remove(listall, &iter);
	}
	gtk_tree_row_reference_free(row);
#endif
	return 0;//must return false
}
/**
   Update the job list entry of p (adding it if needed) and its status icon
   according to p->status, and notify the user of status changes. Called through
   g_idle_add() by the monitor_thread.

   @param[in,out] p The job.
   @return 0 so that it runs only once.
*/
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

#if USE_COLUMN_VIEW
/**
   Definition of a job list column.
*/
typedef struct{
	int type;//0: text, 1: progress, 2: icon
	int width;//minimum width. Negative to make the column expandable (and text ellipsized). 0 is not resizable.
	int col;//column of the text, progress text, or icon
	int col2;//column of the text color (type 0) or progress value (type 1). -1 for none.
}coldef_t;
/**
   Get the paintable of an icon. It is created on first use and kept with the
   icon.

   @param[in] icon The icon. May be NULL.
   @return Reference to the paintable (do not free), or NULL.
*/
static GdkPaintable* icon_paintable(GdkPixbuf* icon){
	if(!icon) return NULL;
	GdkPaintable* paintable=(GdkPaintable*)g_object_get_data(G_OBJECT(icon), "paintable");
	if(!paintable){
		paintable=GDK_PAINTABLE(gdk_texture_new_for_pixbuf(icon));
		g_object_set_data_full(G_OBJECT(icon), "paintable", paintable, g_object_unref);
	}
	return paintable;
}
/**
   Update the widget of a cell from its job row. Handler of the "changed" signal
   of the job row.

   @param[in] job The job row.
   @param[in] item The list item of the cell.
*/
static void cell_update(MonJob* job, GtkListItem* item){
	GtkWidget* child=gtk_list_item_get_child(item);
	const coldef_t* def=(const coldef_t*)g_object_get_data(G_OBJECT(child), "coldef");
	const char* text=job->str[def->col]?job->str[def->col]:"";
	switch(def->type){
	case 0:
		if(def->col2>=0&&job->str[def->col2]){
			char* markup=g_markup_printf_escaped("<span foreground=\"%s\">%s</span>", job->str[def->col2], text);
			gtk_label_set_markup(GTK_LABEL(child), markup);
			g_free(markup);
		} else{
			gtk_label_set_text(GTK_LABEL(child), text);
		}
		break;
	case 1:
		gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(child), job->stepp*0.01);
		gtk_progress_bar_set_text(GTK_PROGRESS_BAR(child), text);
		break;
	case 2:
		gtk_image_set_from_paintable(GTK_IMAGE(child), icon_paintable(job->icon));
		break;
	}
	gtk_widget_set_tooltip_text(child, job->str[COL_TOOLTIP]);
}
static void job_icon_clicked(const char* host, const char* spid, GdkPixbuf* icon);
/**
   Handle left button release on the icon cell to act on an individual job.

   @param[in] gesture Unused.
   @param[in] n_press Unused.
   @param[in] x Unused.
   @param[in] y Unused.
   @param[in] item The list item of the cell.
*/
static void cell_icon_released(GtkGestureClick* gesture, int n_press, double x, double y, GtkListItem* item){
	(void)gesture; (void)n_press; (void)x; (void)y;
	MonJob* job=(MonJob*)gtk_list_item_get_item(item);
	if(job){
		job_icon_clicked(job->str[COL_HOST], job->str[COL_PID], job->icon);
	}
}
/**
   Handle right button press on a cell: select its row if nothing is selected.

   @param[in] gesture The click gesture.
   @param[in] n_press Unused.
   @param[in] x Unused.
   @param[in] y Unused.
   @param[in] item The list item of the cell.
*/
static void cell_right_pressed(GtkGestureClick* gesture, int n_press, double x, double y, GtkListItem* item){
	(void)n_press; (void)x; (void)y;
	GtkWidget* view=gtk_widget_get_ancestor(gtk_event_controller_get_widget(GTK_EVENT_CONTROLLER(gesture)), GTK_TYPE_COLUMN_VIEW);
	if(!view) return;
	GtkSelectionModel* model=gtk_column_view_get_model(GTK_COLUMN_VIEW(view));
	GtkBitset* selected=gtk_selection_model_get_selection(model);
	if(gtk_bitset_is_empty(selected)){
		gtk_selection_model_select_item(model, gtk_list_item_get_position(item), TRUE);
	}
	gtk_bitset_unref(selected);
	gtk_gesture_set_state(GTK_GESTURE(gesture), GTK_EVENT_SEQUENCE_CLAIMED);//keep the selection
}
/**
   Create the widget of a cell. Handler of the "setup" signal of the column factory.

   @param[in] factory Unused.
   @param[in] item The list item of the cell.
   @param[in] def The column definition.
*/
static void cell_setup(GtkSignalListItemFactory* factory, GtkListItem* item, const coldef_t* def){
	(void)factory;
	GtkWidget* child=NULL;
	GtkGesture* gesture;
	switch(def->type){
	case 0:
		child=gtk_label_new(NULL);
		gtk_label_set_xalign(GTK_LABEL(child), 1);
		if(def->width<0){//set ellipsize makes it shrinkable than minimal needed size.
			gtk_label_set_ellipsize(GTK_LABEL(child), PANGO_ELLIPSIZE_END);
		}
		break;
	case 1:
		child=gtk_progress_bar_new();
		gtk_progress_bar_set_show_text(GTK_PROGRESS_BAR(child), TRUE);
		break;
	case 2:
		child=gtk_image_new();
		gesture=gtk_gesture_click_new();
		gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(gesture), 1);
		g_signal_connect(gesture, "released", G_CALLBACK(cell_icon_released), item);
		gtk_widget_add_controller(child, GTK_EVENT_CONTROLLER(gesture));
		break;
	default:
		error("Invalid\n");
	}
	if(def->width){
		gtk_widget_set_size_request(child, abs(def->width), -1);
	}
	gesture=gtk_gesture_click_new();
	gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(gesture), 3);
	g_signal_connect(gesture, "pressed", G_CALLBACK(cell_right_pressed), item);
	gtk_widget_add_controller(child, GTK_EVENT_CONTROLLER(gesture));
	g_object_set_data(G_OBJECT(child), "coldef", (gpointer)def);
	gtk_list_item_set_child(item, child);
}
/**
   Show the job row in a cell and follow its changes. Handler of the "bind"
   signal of the column factory.

   @param[in] factory Unused.
   @param[in] item The list item of the cell.
   @param[in] def Unused.
*/
static void cell_bind(GtkSignalListItemFactory* factory, GtkListItem* item, const coldef_t* def){
	(void)factory; (void)def;
	MonJob* job=(MonJob*)gtk_list_item_get_item(item);
	g_signal_connect(job, "changed", G_CALLBACK(cell_update), item);
	cell_update(job, item);
}
/**
   Stop following the job row of a cell. Handler of the "unbind" signal of the
   column factory.

   @param[in] factory Unused.
   @param[in] item The list item of the cell.
   @param[in] def Unused.
*/
static void cell_unbind(GtkSignalListItemFactory* factory, GtkListItem* item, const coldef_t* def){
	(void)factory; (void)def;
	MonJob* job=(MonJob*)gtk_list_item_get_item(item);
	if(job){
		g_signal_handlers_disconnect_by_data(job, item);
	}
}
/**
   Create a column of the job list and append it to the view.

   @param[in] view The column view.
   @param[in] type Cell type: 0 text, 1 progress, 2 icon.
   @param[in] width Minimum width. Negative to make the column expandable (and text ellipsized). 0 is not resizable.
   @param[in] title Column title.
   @param[in] col Column of the text (type 0), progress text (type 1), or icon (type 2).
   @param[in] col2 Column of the text color (type 0) or progress value (type 1). -1 for none.
*/
static void append_column(GtkWidget* view, int type, int width, const char* title, int col, int col2){
	coldef_t* def=g_new0(coldef_t, 1);
	def->type=type;
	def->width=width;
	def->col=col;
	def->col2=col2;
	GtkListItemFactory* factory=gtk_signal_list_item_factory_new();
	g_object_set_data_full(G_OBJECT(factory), "coldef", def, g_free);
	g_signal_connect(factory, "setup", G_CALLBACK(cell_setup), def);
	g_signal_connect(factory, "bind", G_CALLBACK(cell_bind), def);
	g_signal_connect(factory, "unbind", G_CALLBACK(cell_unbind), def);
	GtkColumnViewColumn* column=gtk_column_view_column_new(title, factory);
	gtk_column_view_column_set_resizable(column, width?TRUE:FALSE);//allow user to grab and resize
	gtk_column_view_column_set_expand(column, width<0?TRUE:FALSE);
	gtk_column_view_append_column(GTK_COLUMN_VIEW(view), column);
	g_object_unref(column);
}
#else
/**
   Create a treeview column and append it to the view.

   @param[in] view The treeview.
   @param[in] type Renderer type: 0 text, 1 progress, 2 pixbuf.
   @param[in] width Minimum width. Negative to make the column expandable (and text ellipsized). 0 is not resizable.
   @param[in] title Column title.
   @param[in] icol Column of the text (type 0), progress text (type 1), or pixbuf (type 2).
   @param[in] icol2 Column of the text color (type 0) or progress value (type 1). -1 for none.
*/
static void append_column(GtkWidget* view, int type, int width, const char* title, int icol, int icol2){
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
			g_object_set(G_OBJECT(render), "ellipsize", PANGO_ELLIPSIZE_END, NULL);
		}
	}
	const int expandable=(width<0)?TRUE:FALSE;
	const int resizable=width?TRUE:FALSE;//resizeable makes the column very small if not expand.

	gtk_tree_view_column_set_resizable(col, resizable);//allow user to grab and resize
	gtk_tree_view_column_set_expand(col, expandable);//cell is expandable within column
	gtk_tree_view_column_pack_start(col, render, expandable);//column is expandable
	gtk_tree_view_column_add_attribute(col, render, type==2?"pixbuf":"text", icol);
	if(icol2>=0){
		gtk_tree_view_column_add_attribute(col, render, type==1?"value":"foreground", icol2);
	}
	gtk_tree_view_append_column(GTK_TREE_VIEW(view), col);
}
#endif
#if GTK_MAJOR_VERSION>3
/**
   Callback of dialog_confirm() with the user answer.

   @param[in] yes Whether the user answered yes.
   @param[in] data The user data passed to dialog_confirm().
*/
typedef void (*confirm_cb_t)(int yes, gpointer data);
/**
   A pending confirmation of dialog_confirm().
*/
typedef struct{
	confirm_cb_t cb;
	gpointer data;
}confirm_t;
#if USE_ALERT_DIALOG
/**
   Completion handler of the confirmation dialog of dialog_confirm(). Calls the
   callback with the user answer.

   @param[in] dialog The alert dialog.
   @param[in] result The result of the dialog.
   @param[in] user_data The confirm_t. Freed.
*/
static void dialog_confirm_response(GObject* dialog, GAsyncResult* result, gpointer user_data){
	confirm_t* conf=(confirm_t*)user_data;
	int button=gtk_alert_dialog_choose_finish(GTK_ALERT_DIALOG(dialog), result, NULL);
	conf->cb(button==1, conf->data);
	free(conf);
}
#else
/**
   Response handler of the confirmation dialog of dialog_confirm(). Destroys the
   dialog and calls the callback with the user answer.

   @param[in] dialog The dialog.
   @param[in] response_id The response.
   @param[in] user_data The confirm_t. Freed.
*/
static void dialog_confirm_response(GtkDialog* dialog, int response_id, gpointer user_data){
	confirm_t* conf=(confirm_t*)user_data;
	gtk_window_destroy(GTK_WINDOW(dialog));
	conf->cb(response_id==GTK_RESPONSE_YES, conf->data);
	free(conf);
}
#endif
/**
   Callback of the job kill confirmation. Relays cmds to the monitor thread if
   the answer is yes.

   @param[in] yes Whether the user answered yes.
   @param[in] data The 4 command integers to send through sock_main. Freed.
*/
static void scheduler_cmd_wrap2(int yes, gpointer data){
	int* cmds=(int*)data;
	if(yes){
		scheduler_cmd_wrap(cmds[1], cmds[2], cmds[3]);
	}
	free(cmds);
}
#else
typedef void* confirm_cb_t;
#endif
/**
   Show a yes/no confirmation dialog with the formatted message. With GTK4 it
   returns 0 immediately and cb(yes, data) is called with the user response.

   @param[in] cb Response callback (GTK4 only), e.g., scheduler_cmd_wrap2().
   @param[in] data User data for cb (GTK4 only). Ownership is transferred to cb.
   @param[in] format printf style format of the message.
   @return Whether the user answered yes. Always 0 with GTK4.
*/
CHECK_ARG(3)
static gboolean dialog_confirm(confirm_cb_t cb, gpointer data, const char* format, ...){
	format2fn;
#if GTK_MAJOR_VERSION>3
	confirm_t* conf=mycalloc(1, confirm_t);
	conf->cb=cb;
	conf->data=data;
#if USE_ALERT_DIALOG
	GtkAlertDialog* dia=gtk_alert_dialog_new("%s", fn);
	const char* buttons[]={"No", "Yes", NULL};
	gtk_alert_dialog_set_buttons(dia, buttons);
	gtk_alert_dialog_set_cancel_button(dia, 0);
	gtk_alert_dialog_set_default_button(dia, 1);
	gtk_alert_dialog_choose(dia, GTK_WINDOW(window), NULL, dialog_confirm_response, conf);
	g_object_unref(dia);
#else
	GtkWidget* dia=gtk_message_dialog_new
	(GTK_WINDOW(window), (GtkDialogFlags)(GTK_DIALOG_DESTROY_WITH_PARENT|GTK_DIALOG_MODAL),
		GTK_MESSAGE_QUESTION,
		GTK_BUTTONS_YES_NO,
		"%s", fn);
	g_signal_connect(GTK_DIALOG(dia), "response", G_CALLBACK(dialog_confirm_response), conf);
	gtk_window_present(GTK_WINDOW(dia));
#endif
	return 0;
#else
	GtkWidget* dia=gtk_message_dialog_new
	(GTK_WINDOW(window), (GtkDialogFlags)(GTK_DIALOG_DESTROY_WITH_PARENT|GTK_DIALOG_MODAL),
		GTK_MESSAGE_QUESTION,
		GTK_BUTTONS_YES_NO,
		"%s", fn);
	(void)cb; (void)data;
	int result=gtk_dialog_run(GTK_DIALOG(dia));
	gtk_widget_destroy(dia);
	return result==GTK_RESPONSE_YES;
#endif
}
int reset_clipboard=0;

/**
   Append text to clipboard (primary and default).

   @param[in] jobinfo The text to append.
*/
static void clipboard_append(const char* jobinfo){
#if GTK_MAJOR_VERSION > 3
    GdkClipboard *clipboard = gtk_widget_get_clipboard(window);

    char *old = NULL; /*gdk_clipboard_read_text_finish(
        clipboard,
        gdk_clipboard_read_text_async(clipboard, NULL, NULL, NULL),
        NULL
    );
	*/
    char *newer = NULL;
    if (old) {
        newer = stradd(old, jobinfo, "\n", NULL);
        g_free(old);
    } else {
        newer = stradd(jobinfo, "\n", NULL);
    }

    gdk_clipboard_set_text(clipboard, newer);
    free(newer);
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
/**
   Apply an action to a job.

   @param[in] data The menudata_t of the action.
   @param[in] host Short host name of the job.
   @param[in] spid PID of the job as a string.
   @param[in] icon Status icon of the job.
   @param[in] path Starting path of the job (COL_PATH).
   @param[in] out Output folder of the job (COL_OUT).
   @param[in] full Full path and arguments of the job (COL_FULL).
*/
static void handle_job(const menudata_t* data, const char* host, const char* spid, const GdkPixbuf* icon,
	const char* path, const char* out, const char* full){
	int cmd=data->command;
	const char* action=data->action;
	if(cmd<0){
		if(cmd==-1){
			const char* jobinfo=NULL;
			if(!strcmp(action, "CopyPath")){
				jobinfo=path;
			} else if(!strcmp(action, "CopyOutPath")){
				jobinfo=out;
			} else if(!strcmp(action, "Copy")){
				jobinfo=full;
			}
			if(jobinfo){
				clipboard_append(jobinfo);
			}
		}
	} else{
		if(!(cmd==CMD_DRAWSER&&icon!=icon_running)){
			int pid=strtol(spid, NULL, 10);
			int ihost=host2i(host);
			scheduler_cmd_wrap(ihost, pid, cmd);
		}
	}
}
#if !USE_COLUMN_VIEW
/**
   A general routine handle actions to each item.

   @param[in] model The tree model.
   @param[in] path Unused.
   @param[in] iter Iterator of the selected row.
   @param[in] user_data The menudata_t of the action.
*/
static void handle_selection(GtkTreeModel* model, GtkTreePath* path, GtkTreeIter* iter, gpointer user_data){
	(void)path;
	gchar *host=NULL, *spid=NULL, *spath=NULL, *sout=NULL, *sfull=NULL;
	GdkPixbuf* icon=NULL;
	gtk_tree_model_get(model, iter, COL_HOST, &host, COL_PID, &spid, COL_PATH, &spath, COL_OUT, &sout,
		COL_FULL, &sfull, COL_ACTION, &icon, -1);
	handle_job((menudata_t*)user_data, host, spid, icon, spath, sout, sfull);
	g_free(host); g_free(spid); g_free(spath); g_free(sout); g_free(sfull);
	if(icon) g_object_unref(icon);
}
#endif
/**
   Count the selected jobs of the current page.

   @return The number of selected jobs.
*/
static int count_selected(void){
#if USE_COLUMN_VIEW
	GtkSelectionModel* model=gtk_column_view_get_model(GTK_COLUMN_VIEW(views[cur_ihost]));
	GtkBitset* selected=gtk_selection_model_get_selection(model);
	int nsel=(int)gtk_bitset_get_size(selected);
	gtk_bitset_unref(selected);
	return nsel;
#else
	GtkTreeSelection* selection=gtk_tree_view_get_selection(GTK_TREE_VIEW(views[cur_ihost]));
	return gtk_tree_selection_count_selected_rows(selection);
#endif
}
/**
   Apply an action to each selected job of the current page.

   @param[in] data The menudata_t of the action.
*/
static void handle_selected(menudata_t* data){
	reset_clipboard=1;
#if USE_COLUMN_VIEW
	GtkSelectionModel* model=gtk_column_view_get_model(GTK_COLUMN_VIEW(views[cur_ihost]));
	GtkBitset* selected=gtk_selection_model_get_selection(model);
	GtkBitsetIter it;
	guint pos;
	for(gboolean ok=gtk_bitset_iter_init_first(&it, selected, &pos); ok; ok=gtk_bitset_iter_next(&it, &pos)){
		MonJob* job=(MonJob*)g_list_model_get_item(G_LIST_MODEL(model), pos);
		if(job){
			handle_job(data, job->str[COL_HOST], job->str[COL_PID], job->icon, job->str[COL_PATH], job->str[COL_OUT], job->str[COL_FULL]);
			g_object_unref(job);
		}
	}
	gtk_bitset_unref(selected);
#else
	GtkTreeSelection* selection=gtk_tree_view_get_selection(GTK_TREE_VIEW(views[cur_ihost]));
	gtk_tree_selection_selected_foreach(selection, handle_selection, data);
#endif
}
#if GTK_MAJOR_VERSION>3
/**
   Callback of the confirmation in handle_menu_event(). Applies the action to
   the selected jobs if the answer is yes.

   @param[in] yes Whether the user answered yes.
   @param[in] data The menudata_t of the action.
*/
static void handle_menu_event_response(int yes, gpointer data){
	if(yes){
		handle_selected((menudata_t*)data);
	}
}
#endif
/**
   Handles menu item clicks in a general way.

   @param[in] menuitem Unused.
   @param[in] user_data The menudata_t of the action.
*/
void handle_menu_event(GtkWidget* menuitem, gpointer user_data){
	(void)menuitem;
	menudata_t* data=(menudata_t*)user_data;
	int nsel=count_selected();
	if(nsel<1){
		warning("nsel=%d\n", nsel);
		return;
	}
	gboolean ans=1;
	const char* action=data->action;
	/*Alert user in Kill or Restart event*/
	if(!strcmp(action, "Kill")||!strcmp(action, "Restart")){
#if GTK_MAJOR_VERSION>3
		ans=dialog_confirm(handle_menu_event_response, data, "%s %d jobs?", action, nsel);
#else
		ans=dialog_confirm(NULL, NULL, "%s %d jobs?", action, nsel);
#endif
	}
	if(ans){
		handle_selected(data);
	}
}
/**
   Apply action imenu (index into menudata: plot, clear, kill, restart, or copy
   command line / paths) to the selected jobs.

   @param[in] imenu Index into menudata. Ignored if out of range.
*/
void handle_job_action(int imenu){
	if(imenu>-1 && imenu<9){
		handle_menu_event(NULL, menudata+imenu);
	}
}
/**
   Act on an individual job when its status icon is clicked: kill a running or
   waiting job after confirmation, or clear it otherwise.

   @param[in] host Short host name of the job.
   @param[in] spid PID of the job as a string.
   @param[in] icon Status icon of the job.
*/
static void job_icon_clicked(const char* host, const char* spid, GdkPixbuf* icon){
	int pid=strtol(spid, NULL, 10);
	int ihost=host2i(host);
	if(icon==icon_running||icon==icon_waiting){
#if GTK_MAJOR_VERSION>3
		int* cmds=mymalloc(4, int);
		cmds[0]=MON_CMD; cmds[1]=ihost; cmds[2]=pid; cmds[3]=CMD_KILL;
		dialog_confirm(scheduler_cmd_wrap2, cmds, "Kill %d?", pid);
#else
		if(dialog_confirm(NULL, NULL, "Kill %d?", pid)){
			scheduler_cmd_wrap(ihost, pid, CMD_KILL);
		}
#endif
	}else{
		scheduler_cmd_wrap(ihost, pid, CMD_REMOVE);
	}
}
#if GTK_MAJOR_VERSION<4
/**
   Pop up the job menu for the selected jobs of the treeview.

   @param[in] view The treeview.
   @param[in] user_data Unused.
   @return TRUE.
*/
static gboolean view_popup_menu(GtkWidget* view, gpointer user_data){
	(void) user_data;
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
		//cur_ihost=GPOINTER_TO_INT(wrap_ihost);
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
#if !USE_COLUMN_VIEW
#if GTK_MAJOR_VERSION<4
/**
   Handle mouse button press on the treeview. A right click selects the clicked
   row if nothing is selected and pops up the job menu.

   @param[in] view The treeview.
   @param[in] event The button event.
   @param[in] wrap_ihost Host index of the page, as GINT_TO_POINTER().
   @return TRUE for right click to keep the selection, FALSE otherwise.
*/
static gboolean view_click_event(GtkWidget* view, GdkEventButton* event, gpointer wrap_ihost){
	gdouble x=event->x;
	gdouble y=event->y;
	int button=event->button;
#else
/**
   Handle mouse button press on the treeview. A right click selects the clicked
   row if nothing is selected.

   @param[in] self The click gesture.
   @param[in] n_press Unused.
   @param[in] x Horizontal position of the click in the treeview.
   @param[in] y Vertical position of the click in the treeview.
   @param[in] wrap_ihost Unused.
   @return TRUE for right click to keep the selection, FALSE otherwise.
*/
static gboolean view_click_event(GtkGestureClick* self, gint n_press, gdouble x, gdouble y, gpointer wrap_ihost){
	int button=gtk_gesture_single_get_current_button(GTK_GESTURE_SINGLE(self));
	GtkWidget *view=gtk_event_controller_get_widget(GTK_EVENT_CONTROLLER(self));
	(void)n_press;
	(void)wrap_ihost;
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
		view_popup_menu(view, wrap_ihost);
#endif
		return TRUE;//return TRUE to avoid modifying selection
	}
	return FALSE;//default handler handles multi selection
}
#if GTK_MAJOR_VERSION<4
/**
   Handle left button release on the icon column to act on an individual job
   (see job_icon_clicked()).

   @param[in] view The treeview.
   @param[in] event The button event.
   @param[in] user_data Unused.
   @return FALSE if not a left click, TRUE otherwise.
*/
static gboolean view_release_event(GtkWidget* view, GdkEventButton* event, gpointer user_data){
	int button=event->button;
	gdouble x=event->x;
	gdouble y=event->y;
#else
/**
   Handle left button release on the icon column to act on an individual job
   (see job_icon_clicked()).

   @param[in] self The click gesture.
   @param[in] sequence The event sequence used to get the position.
   @param[in] user_data Unused.
   @return FALSE if not a left click, TRUE otherwise.
*/
static gboolean view_release_event(GtkGesture*self, GdkEventSequence *sequence, gpointer user_data){
	int button=gtk_gesture_single_get_current_button(GTK_GESTURE_SINGLE(self));
	GtkWidget *view=gtk_event_controller_get_widget(GTK_EVENT_CONTROLLER(self));
	gdouble x=0, y=0;
	gtk_gesture_get_point(self, sequence, &x, &y);
#endif
	(void) user_data;
	if(button!=1) return FALSE;//only handle left click
	GtkTreePath* path=NULL;
	GtkTreeViewColumn* column=NULL;
	if(gtk_tree_view_get_path_at_pos(GTK_TREE_VIEW(view),
		(gint)x, (gint)y, &path, &column, NULL, NULL)){
		const gchar* title=gtk_tree_view_column_get_title(column);
		if(!strcmp(title, " ")){
			GtkTreeIter iter;
			GtkTreeModel* model=gtk_tree_view_get_model(GTK_TREE_VIEW(view));
			gtk_tree_model_get_iter(model, &iter, path);
			gchar *host=NULL, *spid=NULL;
			GdkPixbuf* icon=NULL;
			gtk_tree_model_get(model, &iter, COL_HOST, &host, COL_PID, &spid, COL_ACTION, &icon, -1);
			job_icon_clicked(host, spid, icon);
			g_free(host); g_free(spid);
			if(icon) g_object_unref(icon);
		}
		gtk_tree_path_free(path);
	}
	return TRUE;
}
#endif
#if USE_COLUMN_VIEW
/**
   Visibility filter of a host page: show jobs of the host.

   @param[in] item The job row.
   @param[in] host Short host name to match.
   @return TRUE if the row belongs to the host.
*/
static gboolean filter_host(gpointer item, gpointer host){
	const char* host2=MON_JOB(item)->str[COL_HOST];
	return host2&&(!strcmp(host2, (gchar*)host));
}
/**
   Visibility filter of the "All" page: show running, waiting, and finished jobs.

   @param[in] item The job row.
   @param[in] status Unused.
   @return TRUE if the row is shown.
*/
static gboolean filter_status(gpointer item, gpointer status){
	(void)status;
	GdkPixbuf* status2=MON_JOB(item)->icon;
	return (status2==icon_running||status2==icon_finished||status2==icon_waiting);
}
/**
   Order job rows by their start time.

   @param[in] a The first job row.
   @param[in] b The second job row.
   @param[in] user_data Unused.
   @return Negative, zero, or positive if a starts before, with, or after b.
*/
static int sort_time(gconstpointer a, gconstpointer b, gpointer user_data){
	(void)user_data;
	long ta=((const MonJob*)a)->time;
	long tb=((const MonJob*)b)->time;
	return (ta>tb)-(ta<tb);
}
#else
/**
   Visibility filter of a host page: show jobs of the host.

   @param[in] model The tree model.
   @param[in] iter The row.
   @param[in] host Short host name to match.
   @return TRUE if the row belongs to the host.
*/
static gboolean filter_host(GtkTreeModel* model, GtkTreeIter* iter, gpointer host){
	gchar* host2;
	gtk_tree_model_get(model, iter, COL_HOST, &host2, -1);
	int ans=host2&&(!strcmp(host2, (gchar*)host));
	g_free(host2);
	return ans;
}
/**
   Visibility filter of the "All" page: show running, waiting, and finished jobs.

   @param[in] model The tree model.
   @param[in] iter The row.
   @param[in] status Unused.
   @return TRUE if the row is shown.
*/
static gboolean filter_status(GtkTreeModel* model, GtkTreeIter* iter, gpointer status){
	GdkPixbuf* status2=0;
	(void)status;
	gtk_tree_model_get(model, iter, COL_ACTION, &status2, -1);
	if(status2) g_object_unref(status2);//only the address is compared
	return (status2==icon_running||status2==icon_finished||status2==icon_waiting);
}
#endif
/**
   Get the job list view for host ihost, creating it (and the shared list
   store) on first call. ihost equal to nhost is the "All" page that shows
   running, waiting, and finished jobs of all hosts.

   @param[in] ihost Host index, or nhost for the "All" page.
   @return The treeview (GtkColumnView with GTK 4.10 or later). Owned by the module; do not free.
*/
GtkWidget* new_page(int ihost){
#if USE_COLUMN_VIEW
	if(!listall){
		listall=g_list_store_new(MON_TYPE_JOB);
		filters=mycalloc(nhost+1, GtkFilter*);
		views=mycalloc(nhost+1, GtkWidget*);
	}
	if(views[ihost]){
		return views[ihost];
	}
	//listall -> sort by start time -> filter by host or status -> multiple selection
	GtkSorter* sorter=GTK_SORTER(gtk_custom_sorter_new(sort_time, NULL, NULL));
	GtkSortListModel* sorted=gtk_sort_list_model_new(G_LIST_MODEL(g_object_ref(listall)), sorter);
	GtkCustomFilter* filter;
	if(ihost<nhost){
		filter=gtk_custom_filter_new(filter_host, hostshort[ihost], NULL);
	} else{
		filter=gtk_custom_filter_new(filter_status, NULL, NULL);
	}
	filters[ihost]=GTK_FILTER(filter);//owned by the filter model
	GtkFilterListModel* filtered=gtk_filter_list_model_new(G_LIST_MODEL(sorted), GTK_FILTER(filter));
	GtkMultiSelection* selection=gtk_multi_selection_new(G_LIST_MODEL(filtered));
	GtkWidget* view;
	views[ihost]=view=gtk_column_view_new(GTK_SELECTION_MODEL(selection));
	gtk_column_view_set_show_column_separators(GTK_COLUMN_VIEW(view), TRUE);
#else
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
		gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(gesture), 0);
		g_signal_connect(gesture, "pressed", G_CALLBACK(view_click_event), GINT_TO_POINTER(ihost));
		g_signal_connect(gesture, "released", G_CALLBACK(view_release_event), GINT_TO_POINTER(ihost));
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
#endif
	append_column(view, 0, 0, "Date", COL_DATE, -1);
	if(ihost==nhost){
		append_column(view, 0, 0, "Host", COL_HOST, -1);
	} else{
		append_column(view, 0, 0, "PID", COL_PID, -1);
	}
#if SHOW_START
	append_column(view, 0, 100, "Start Dir", COL_PATH, -1);
#endif
	append_column(view, 0, -100, "Arguments", COL_ARGS, -1);
#if SHOW_SOUT
	append_column(view, 0, 50, "Out Dir", COL_OUT, -1);
#endif
	append_column(view, 0, 0, "Low", COL_ERRLO, COL_COLOR);
	append_column(view, 0, 0, "High", COL_ERRHI, COL_COLOR);
	//append_column(view, 1, 0, "Seed", COL_SEED, COL_SEEDP);

	append_column(view, 0, 0, "Step", COL_STEPT, -1);

	append_column(view, 0, 0, "Tot", COL_TIMING, -1);

	append_column(view, 2, 0, " ", COL_ACTION, -1);
	append_column(view, 1, 0, "Progress", COL_STEP, COL_STEPP);
#if !USE_COLUMN_VIEW
	gtk_tree_view_columns_autosize(GTK_TREE_VIEW(view));
#endif
	return view;
}


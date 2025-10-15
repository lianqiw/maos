/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include <gdk/gdkkeysyms.h>
#include "drawdaemon.h"
#ifndef GTK_WIDGET_HAS_FOCUS
#define GTK_WIDGET_HAS_FOCUS gtk_widget_has_focus
#endif
#if GTK_MAJOR_VERSION>=3
#define GDK_Left GDK_KEY_Left
#define GDK_Right GDK_KEY_Right
#define GDK_Up GDK_KEY_Up
#define GDK_Down GDK_KEY_Down
#define GDK_plus GDK_KEY_plus
#define GDK_minus GDK_KEY_minus
#define GDK_equal GDK_KEY_equal
#define GDK_0 GDK_KEY_0
#define GDK_1 GDK_KEY_1
#endif
#define DRAWAREA_MIN_WIDTH 400
#define DRAWAREA_MIN_HEIGHT 320
#define MAX_ZOOM 1000 //larger zoom triggers problems in cairo (points missing)
#define MIN_ZOOM 1
static GtkWidget *fpslabel=NULL;
static GSList* windows=NULL;
static GtkWidget* curwindow=NULL;
static GtkWidget* curtopnb=NULL;
#if GTK_MAJOR_VERSION<4
//static GtkWidget* contexmenu=NULL;
#endif
static GtkWidget* cur_menu_cumu=NULL;
static GtkWidget *cur_menu_icumu=NULL;
static GtkWidget *cur_menu_zlog=NULL;
static drawdata_t* drawdata_dialog=NULL;
//static GtkToolItem* toggle_cumu=NULL;
#if GTK_MAJOR_VERSION<3
static GtkRcStyle* btn_rcstyle=NULL;
#endif
PangoFontDescription* desc=NULL;
int desc_font_size=0;
static void subnb_page_removed(GtkNotebook *subnb, GtkWidget *child, guint n, GtkWidget *topnb);
int font_name_version=0;
char* font_name=NULL;
const char *font_name_default="Arial Regular 12";
char **pfilename_gif=NULL;
float font_size;
//cairo_font_slant_t font_style=CAIRO_FONT_SLANT_NORMAL;
//cairo_font_weight_t font_weight=CAIRO_FONT_WEIGHT_NORMAL;
float lpf=1;//low pass filter the update. Set using dialog.
int window_proportional=1;//resize window proportionally
/*
  Routines in this file are about the GUI.
*/

GdkCursor* cursors[2]={0,0};
//GdkPixbuf *pix_hand=NULL;
//GdkPixbuf *pix_arrow=NULL;
#if GTK_MAJOR_VERSION < 3
static const char* rc_string_notebook={
	"style \"noborder\"{                      \n"
	"GtkNoteBook::draw-border={0 0 0 0}       \n"
	"}                                        \n"
	"class \"GtkNoteBook\" style \"noborder\" \n"
};
#endif
#if GTK_MAJOR_VERSION>=3
#include "gtk3-css.h"
#endif
#if GTK_VERSION_AFTER(4,10)
#define error_msg(A...) {\
	GtkAlertDialog *dialog0=gtk_alert_dialog_new(A);\
	gtk_alert_dialog_set_modal(dialog0, true);\
	const char*buttons[]={"OK", NULL};\
	gtk_alert_dialog_set_buttons(dialog0, buttons);\
	gtk_alert_dialog_show(dialog0, GTK_WINDOW(curwindow));\
}
#else
#define error_msg_create(A...) 		\
	GtkWidget *dialog0=gtk_message_dialog_new	\
	    (GTK_WINDOW(curwindow),		\
	     GTK_DIALOG_DESTROY_WITH_PARENT,		\
	     GTK_MESSAGE_ERROR,				\
	     GTK_BUTTONS_CLOSE,				\
	     A)
#if GTK_VERSION_AFTER(4,0)
#define error_msg(A...)	{error_msg_create(A);\
	gtk_window_present(GTK_WINDOW(dialog0));}
#else
#define error_msg(A...) {error_msg_create(A);\
	gtk_dialog_run(GTK_DIALOG(dialog0));\
	gtk_widget_destroy(dialog0);}
#endif

#endif

static GtkWidget *get_topnb(GtkWidget *window){
#if USE_TITLEBAR==0
	GtkWidget *vbox=gtk_bin_get_child(GTK_BIN(window));
	GList *list=gtk_container_get_children(GTK_CONTAINER(vbox));
	GtkWidget* topnb=(GtkWidget *)list->next->data;
	g_list_free(list);
#elif GTK_MAJOR_VERSION<4
	/*GtkWidget* vbox=gtk_bin_get_child(GTK_BIN(window));
	GList* list=gtk_container_get_children(GTK_CONTAINER(vbox));
	GtkWidget* topnb=(GtkWidget*)list->next->data;
	g_list_free(list);*/
	GtkWidget *topnb=gtk_bin_get_child(GTK_BIN(window));
#else
	GtkWidget *topnb=gtk_window_get_child(GTK_WINDOW(window));
#endif
	return topnb;
}
static GtkWidget *get_toolbar(GtkWidget *window){
#if USE_TITLEBAR==0
	GtkWidget *vbox=gtk_bin_get_child(GTK_BIN(window));
	GList *list=gtk_container_get_children(GTK_CONTAINER(vbox));
	GtkWidget *toolbar=(GtkWidget *)list->data;
	g_list_free(list);
#else
	/*GtkWidget* vbox=gtk_bin_get_child(GTK_BIN(window));
	GList* list=gtk_container_get_children(GTK_CONTAINER(vbox));
	GtkWidget* toolbar=(GtkWidget*)list->data;
	g_list_free(list);*/
	GtkWidget *toolbar=gtk_window_get_titlebar(GTK_WINDOW(window));
#endif
	return toolbar;
}

static const char *topnb_label_get(GtkWidget *topnb, GtkWidget *subnb){
	GtkWidget *label=gtk_notebook_get_tab_label(GTK_NOTEBOOK(topnb), subnb);
	return gtk_label_get_text(GTK_LABEL(label));
}

static void window_changed(GtkWidget* window){
	//info("window_changed: %p\n", window);
	curwindow=window;
	cur_menu_cumu=(GtkWidget*)g_object_get_data(G_OBJECT(window), "menu_cumu");
	cur_menu_icumu=(GtkWidget*)g_object_get_data(G_OBJECT(window), "menu_icumu");
	cur_menu_zlog=(GtkWidget*)g_object_get_data(G_OBJECT(window), "menu_zlog");
	curtopnb=get_topnb(window);
}
/*Get the current page for notebook*/
static GtkWidget* get_current_page(GtkWidget* notebook){
	if(!notebook){
		dbg_time("notebook is null\n");
		return NULL;
	}
	int n=gtk_notebook_get_current_page(GTK_NOTEBOOK(notebook));
	return gtk_notebook_get_nth_page(GTK_NOTEBOOK(notebook), n);
}
static drawdata_t* get_current_drawdata(void){
	if(!curtopnb) {
		dbg_time("curtopnb is null\n");
		return NULL;
	}
	GtkWidget* topnb=curtopnb;
	GtkWidget* subnb=get_current_page(topnb);
	GtkWidget* drawarea=get_current_page(subnb);

	return drawarea?(drawdata_t*)g_object_get_data(G_OBJECT(drawarea), "drawdata"):NULL;
}

/**
 * Move a tab from one notebook topnb to another one topnb2
*/
static void move_tab_page(GtkWidget* topnb, int ipage, GtkWidget* topnb2){
	GtkWidget *page=gtk_notebook_get_nth_page(GTK_NOTEBOOK(topnb), ipage);
	GtkWidget *label=gtk_notebook_get_tab_label(GTK_NOTEBOOK(topnb), page);
	g_signal_handlers_disconnect_by_func(page, (gpointer)subnb_page_removed, topnb);
	g_object_ref(page);
	g_object_ref(label);
	const gchar *text=topnb_label_get(topnb, page);
	int itab;
	for(itab=0; itab<gtk_notebook_get_n_pages(GTK_NOTEBOOK(topnb2)); itab++){
		GtkWidget *page2=gtk_notebook_get_nth_page(GTK_NOTEBOOK(topnb2), itab);
		const gchar *text2=topnb_label_get(topnb2, page2);
		if(strcmp(text2, text)>0){
			break;
		}
	}
	gtk_notebook_remove_page(GTK_NOTEBOOK(topnb), ipage);
	gtk_notebook_insert_page(GTK_NOTEBOOK(topnb2), page, label, itab);
	g_signal_connect(page, "page-removed", G_CALLBACK(subnb_page_removed), topnb2);
#if GTK_MAJOR_VERSION<4
	gtk_notebook_set_tab_detachable(GTK_NOTEBOOK(topnb2), page, TRUE);
#endif
	gtk_notebook_set_tab_reorderable(GTK_NOTEBOOK(topnb2), page, TRUE);
	g_object_unref(page);
	g_object_unref(label);
}

static void topnb_detach_btn(GtkWidget *btn, gpointer data){
	GtkWidget *topnb=GTK_WIDGET(data);
	(void)btn;
	if(!curtopnb){
		dbg_time("curtopnb is null\n");
		return ;
	}
	//GtkWidget *topnb=curtopnb;
	if(gtk_notebook_get_n_pages(GTK_NOTEBOOK(topnb))>1){
		int n=gtk_notebook_get_current_page(GTK_NOTEBOOK(topnb));
		GtkWidget *window=create_window(NULL);/*create a new window. */
		GtkWidget *topnb2=get_topnb(window);/*another topnb */
		move_tab_page(topnb, n, topnb2);
	}
}
static void update_toolbar(){//always use current page
	drawdata_t *drawdata=get_current_drawdata();
	if(!drawdata) return;
	int cumu_supported=!drawdata->p&&drawdata->npts>0&&!drawdata->square;
	if(cur_menu_cumu){
		toggle_button_set_active(cur_menu_cumu, drawdata->cumu);
		gtk_widget_set_sensitive(cur_menu_cumu, cumu_supported);
	}
	if(cur_menu_icumu){
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(cur_menu_icumu), drawdata->cumu?drawdata->icumu:0);
		gtk_widget_set_sensitive(cur_menu_icumu, cumu_supported);
	}
	if(cur_menu_zlog){
		toggle_button_set_active(cur_menu_zlog, drawdata->zlog);
		gtk_widget_set_sensitive(cur_menu_zlog, drawdata->p?TRUE:FALSE);
	}
}
gboolean finalize_gif(void*data){
	(void)data;
	if(pfilename_gif && *pfilename_gif){
		info_time("finalizing gif for %s\n", *pfilename_gif);
		const char *convert[]={"/usr/bin/convert", "/usr/local/bin/convert", "/opt/homebrew/bin/convert"};
		for(unsigned int i=0; i<sizeof(convert)/sizeof(convert[0]); i++){
			if(exist(convert[i])){
				char args[PATH_MAX];
				snprintf(args, PATH_MAX, "%s -delay 5 %s/???.png %s.gif &", convert[i], *pfilename_gif, *pfilename_gif);
				if(system(args)==-1){
					warning("run %s failed\n", args);
				}else{
					break;
				}
			}
		}
		free(*pfilename_gif);
		*pfilename_gif=NULL;
	}
	return FALSE;
}
/*
	drawdaemon_draw draws to pixmap2 and then to pixmap to avoid flickering
*/
static void update_pixmap(drawdata_t* drawdata){
	if(!drawdata) return;
	/*no more pending updates, do the updating. */
	if(__atomic_load_n(&drawdata->recycle, __ATOMIC_SEQ_CST)){
		warning_time("recycle is set, do not draw\n");
		return;
	}
	drawdata->iframe++;
	gboolean hasdata=(drawdata->nx&&drawdata->ny)||drawdata->npts||drawdata->ncir;
	if(hasdata){
		int locked=0;
		gint width=drawdata->width;
		gint height=drawdata->height;
		if(drawdata->pixmap){
			if(width!=drawdata->pwidth||height!=drawdata->pheight){//need a new size
				//need to lock mutex to synchronize with drawdaemon_draw which uses pixmap
				//destroy of pixmap need to happen in Gtk main thread because it is used to paint drawarea
				pthread_mutex_lock(&drawdata->mutex); locked=1;
				cairo_surface_destroy(drawdata->pixmap);
				drawdata->pixmap=NULL;
			}
		}
		if(!drawdata->pixmap){
			/*Create a new server size pixmap and then draw on it. */
	#if GTK_MAJOR_VERSION>=3
			int gtk_scale=gtk_widget_get_scale_factor(curwindow);//since 3.10
	#else
			int gtk_scale=1;
	#endif
			if(!locked){
				drawdata->pwidth=width;
				drawdata->pheight=height;
				pthread_mutex_lock(&drawdata->mutex); locked=1;
			}
			drawdata->pixmap=cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width*gtk_scale, height*gtk_scale);
	#if GTK_MAJOR_VERSION>=3
			cairo_surface_set_device_scale(drawdata->pixmap, gtk_scale, gtk_scale);
	#endif
			pthread_mutex_unlock(&drawdata->mutex); 
		}
		cairo_draw(drawdata);//cairo_draw calls drawarea_refresh()
		//g_idle_add((GSourceFunc)drawarea_refresh, drawdata->drawarea);
		//info("%d %d %d\n", drawdata->frame_io, drawdata->frame_draw, drawdata->frame_gif);
	}else{
		//dbg_time("hasdata=0\n");
	}
}
/// @brief called by drawdaemon_draw to refresh the widget
/// @param drawdata 
/// @return FALSE to run only once
gboolean drawarea_refresh(GtkWidget *drawarea){
	if(drawarea){
		gtk_widget_queue_draw(drawarea);
	}
	return FALSE;//run only once
}


#if  GTK_MAJOR_VERSION>=4
static gboolean
on_resize_event(GtkDrawingArea* widget, int width, int height, drawdata_t *drawdata){
	(void)widget;
	//info("on_resize_event with %dx%d\n", width, height);
	if(width>1&&height>1){
		drawdata->width=width;
		drawdata->height=height;
	}
	if(drawdata==get_current_drawdata()){
		drawdata->drawn=0;
	}else{
		update_pixmap(drawdata);
	}
	return FALSE;
}
#else
/*
   Expose event happens when the widget is first configured and whenever it is
   resized.  Create a new backup pixmap of the appropriate size and draw there.*/
static gboolean
on_configure_event(GtkWidget* widget, GdkEventConfigure* event, drawdata_t *drawdata){
	(void)event;
	(void)widget;
	drawdata->width=event->width;
	drawdata->height=event->height;
	if(drawdata!=get_current_drawdata()){
		drawdata->drawn=0;//mark for future drawing
	}else if(event->width>1&&event->height>1){
		update_pixmap(drawdata);
	}
	return FALSE;
}
#endif
#if GTK_MAJOR_VERSION>=4
void drawarea_draw_func(GtkDrawingArea *widget, cairo_t *cr, int width, int height, gpointer pdata){
	(void)widget; (void)width; (void)height;
	drawdata_t *drawdata=(drawdata_t *)pdata;
	//info("draw_func called with %dx%d. pixmap=%p\n", width, height, drawdata->pixmap);
	drawdata->width=width;
	drawdata->height=height;
#elif GTK_MAJOR_VERSION>=3
static gboolean on_draw_event(GtkWidget* widget, cairo_t* cr, gpointer pdata){
	drawdata_t *drawdata=(drawdata_t *)pdata;
	(void)widget;
#else
static gboolean on_expose_event(GtkWidget*widget, GdkEventExpose*event, gpointer pdata){
	(void) event;
	drawdata_t *drawdata=(drawdata_t *)pdata;
	if(!drawdata) return;
	if(drawdata->cumu!=drawdata->cumu_last || drawdata->zlog!=drawdata->zlog_last){
		update_toolbar();
	}
	cairo_t *cr=gdk_cairo_create(widget->window);
#endif
	if(drawdata->pixmap){
		//pthread_mutex_lock(&drawdata->mutex);
		//info_time("painting drawarea: ok\n");
		cairo_set_source_surface(cr, drawdata->pixmap, 0, 0);
		cairo_paint(cr);
		//pthread_mutex_unlock(&drawdata->mutex);

		if(drawdata->draw_rect){
			cairo_set_source_rgba(cr, 0, 0, 1, 0.1);
			cairo_set_line_width(cr, 1);
			cairo_rectangle(cr, drawdata->mxdown, drawdata->mydown, drawdata->dxdown, drawdata->dydown);
			cairo_fill_preserve(cr);
			cairo_set_source_rgba(cr, 0, 0, 0, 1);
			cairo_stroke(cr);
			drawdata->draw_rect=0;
		}
	}
	if(drawdata->font_name_version!=font_name_version||!drawdata->drawn||drawdata->cumu!=drawdata->cumulast){
		update_pixmap(drawdata);//it queues another draw in the background.
	}
#if GTK_MAJOR_VERSION<3
	cairo_destroy(cr);
#endif

#if GTK_MAJOR_VERSION<=3
	return FALSE;
#endif
}

/**
   Delete a figure page.
*/
int delete_page(drawdata_t *drawdata){
	//dbg_time("deleting page: for %p\n",drawdata);
	__atomic_store_n(&drawdata->recycle, 1, __ATOMIC_SEQ_CST);
	GtkWidget *subnb=drawdata->subnb;
	if(subnb){
		int ipage=gtk_notebook_page_num(GTK_NOTEBOOK(subnb), drawdata->drawarea);
		if(ipage!=-1){
			//dbg_time("delete_page:call page removal for %p (drawarea=%p) at index %d\n",drawdata, drawdata->drawarea, ipage);
			gtk_notebook_remove_page(GTK_NOTEBOOK(subnb), ipage);
			//dbg_time("delete_page:removed\n");
			//drawdata->drawarea=NULL;
			//drawdata->subnb=NULL;
			
		}else{
			warning_time("delete_page: page not found\n");
		}
	}else{
		warning_time("delete_page: %p subnb is NULL\n", drawdata);
	}
	return 0;
}
void drawarea_deleted(gpointer data){
	drawdata_t *drawdata=(drawdata_t*)data;
	if(drawdata){
		//dbg_time("%p: drawarea (%p) is deleted. delete drawdata\n", drawdata, drawdata->drawarea);
		pthread_mutex_lock(&drawdata->mutex);

		if(drawdata->image){
			cairo_surface_destroy(drawdata->image);	drawdata->image=NULL;
		}
		if(drawdata->pixmap){
			cairo_surface_destroy(drawdata->pixmap);drawdata->pixmap=NULL;
		}
		if(drawdata->pixmap2){
			cairo_surface_destroy(drawdata->pixmap2);drawdata->pixmap2=NULL;
		}
		if(drawdata->cacheplot){
			cairo_surface_destroy(drawdata->cacheplot);	drawdata->cacheplot=NULL;
		}
		//dbg_time("%p: pixmap deleted\n",drawdata);
		__atomic_store_n(&drawdata->recycle, 2,__ATOMIC_SEQ_CST);
		pthread_mutex_unlock(&drawdata->mutex);
	}else{
		info_time("drawdata is null\n");
	}
}
gboolean update_title(gpointer window){
	if(!window) window=curwindow;
	if(!window || !GTK_IS_WINDOW(window)) return FALSE;//deleted
	//static int client_pid_last=-2;
	//if(client_pid_last!=client_pid){
	//	client_pid_last=client_pid;
	char title[200];
	extern int client_pidold;
	const char *name=client_exename?client_exename:"Drawdaemon";
	snprintf(title, sizeof(title), "%s (%s:%d %s)", name, client_hostname, client_pid>0?client_pid:client_pidold, client_path?client_path:"");
	gtk_window_set_title(GTK_WINDOW(window), title);
		//update_fpslabel(fpslabel);
	//}else if(client_pid>0){
	//	update_fpslabel(fpslabel);//update frame rate.
	//}
	//return client_pid<0?FALSE:TRUE;//returns true to run repeatedly.
	return FALSE;//when called by g_idle_add, must return FALSE;
}
static GtkWidget* subnb_label_new(const gchar *str, drawdata_t *drawdata){
	GtkWidget* out;
	out=gtk_hbox_new(FALSE, 0);
	/* create label for tab */
	GtkWidget* label=gtk_label_new(str);

#if GTK_MAJOR_VERSION >= 3
	gtk_widget_set_halign(label, GTK_ALIGN_START);
	gtk_widget_set_margin_start(label, 5);//since 3.12
	//gtk_widget_set_margin_end(label, 0);
	//gtk_widget_set_margin_top(label, 0);//no use
	//gtk_widget_set_margin_bottom(label, 0);//no use

	gtk_widget_set_hexpand(label, 1);
	gtk_label_set_xalign(GTK_LABEL(label), 0);//works
#else
	gtk_misc_set_padding(GTK_MISC(label), 0, 0);
	gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
#endif
	gtk_label_set_width_chars(GTK_LABEL(label), 12);
	gtk_label_set_justify(GTK_LABEL(label), GTK_JUSTIFY_LEFT);//this is the default

	box_append(GTK_BOX(out), label, TRUE, TRUE, 0);

	GtkWidget* close_btn=gtk_button_new();
#if GTK_MAJOR_VERSION>=4
	gtk_button_set_icon_name(GTK_BUTTON(close_btn), "window-close");
#else
	GtkWidget *image=gtk_image_new_from_icon_name("window-close", GTK_ICON_SIZE_BUTTON);
	gtk_button_set_image(GTK_BUTTON(close_btn), image);
	gtk_image_set_pixel_size(GTK_IMAGE(image),10);//works
	//gtk_button_set_always_show_image(GTK_BUTTON(close_btn), TRUE);//>=gtk3
    gtk_button_set_relief(GTK_BUTTON(close_btn), GTK_RELIEF_NONE);
#endif

#if GTK_MAJOR_VERSION<3
	gtk_widget_modify_style(close_btn, btn_rcstyle);//helps a tiny bit
	//gtk_widget_set_size_request(close_btn, 24, 24);//just crops faceplate. does not reduce icon.
	//gtk_container_set_border_width(GTK_CONTAINER(close_btn), 0);//useless
#endif
	g_signal_connect_swapped(close_btn, "clicked", G_CALLBACK(delete_page), drawdata);
	box_append(GTK_BOX(out), close_btn, FALSE, FALSE, 0);
#if GTK_MAJOR_VERSION < 4
    gtk_widget_show_all(out);
#endif
	return out;
}
/*Get the label of itab'th child of notebook subnb*/
static const char* subnb_label_get(GtkWidget* subnb, GtkWidget* page){
	GtkWidget* label=gtk_notebook_get_tab_label(GTK_NOTEBOOK(subnb), page);
#if GTK_MAJOR_VERSION>=4
	return gtk_label_get_text(GTK_LABEL(gtk_widget_get_first_child(label)));
#else
	GList* list=gtk_container_get_children(GTK_CONTAINER(label));
	return gtk_label_get_text(GTK_LABEL(g_list_first(list)->data));
#endif
}
static void do_move(drawdata_t* drawdata, float xdiff, float ydiff){
	//dbg_time("do_move: %g %g\n", xdiff, ydiff);
	drawdata->offx+=xdiff;
	drawdata->offy+=ydiff;
	update_pixmap(drawdata);/*no need delay since motion notify already did it. */
}
#if GTK_MAJOR_VERSION>=4
static gboolean drawarea_drag_update(GtkGestureDrag *drag, gdouble dx, gdouble dy, drawdata_t *drawdata){
	gint button=gtk_gesture_single_get_current_button(GTK_GESTURE_SINGLE(drag));
	//dbg_time("drag_update at %g %g\n", dx, dy);
#else
static gboolean drawarea_motion_notify(GtkWidget* widget, GdkEventMotion* event, drawdata_t* drawdata){
	gdouble x=event->x;
	gdouble y=event->y;
	gdouble dx=x-drawdata->mxdown;
	gdouble dy=y-drawdata->mydown;
	gint button=(event->state&GDK_BUTTON1_MASK)?1:((event->state&GDK_BUTTON3_MASK)?3:0);
	(void)widget;
#endif
	if(button==1 && drawdata->region==2){//inside legend box
		drawdata->legendoffx+=drawdata->legbox_rx*(dx-drawdata->dxdown); CLIP(drawdata->legendoffx, 0, 1);
		drawdata->legendoffy+=drawdata->legbox_ry*(dy-drawdata->dydown); CLIP(drawdata->legendoffy, 0, 1);
		update_pixmap(drawdata);
		drawdata->dxdown=dx;
		drawdata->dydown=dy;
	}else if(button&&drawdata->region==1){
		//gdouble dt=myclockd()-drawdata->mtdown;
		/*move with left cursor */
		if((fabs(dx)>3 || fabs(dy)>3)){
			if(button==1){
				do_move(drawdata, (dx-drawdata->dxdown), -(dy-drawdata->dydown));/*notice the reverse sign. */
				drawdata->dxdown=dx;
				drawdata->dydown=dy;
			} else if(button==3){/*select and zoom. */
				if(drawdata->square){/*enforce aspect ratio*/
					float ratio=1;
					if(drawdata->p){
						ratio=(float)drawdata->nx/(float)drawdata->ny;
					}
					if(fabs(dx)<fabs(dy)*ratio){
						dy*=fabs(dx/(dy*ratio));
					} else{
						dx*=fabs(dy*ratio/dx);
					}
				}
				drawdata->dxdown=dx;
				drawdata->dydown=dy;
				drawdata->draw_rect=1;
				if(drawdata->pixmap){/*force a refresh to remove previous rectangule */
					gtk_widget_queue_draw(drawdata->drawarea);
				}
			}
		}
	}
#if GTK_MAJOR_VERSION>=4 //a separate controller for cursor motion
	return FALSE;
}
static void drawarea_motion_event(GtkEventControllerMotion *self, gdouble x, gdouble y, drawdata_t *drawdata){
	(void)self;
#endif
	/*we set the cursor */
	if(x>drawdata->xoff&&x < drawdata->xoff+drawdata->widthim
		&&y > drawdata->yoff&&y<drawdata->yoff+drawdata->heightim){
		if(!drawdata->cursorinside){
			drawdata->cursorinside=1;
			//set_cursor(drawdata->drawarea, cursors[0]);
			gtk_widget_set_has_tooltip(drawdata->drawarea, 1);
		}
	
		gdouble x2=(x-drawdata->xoff)/(gdouble)(drawdata->widthim);
		gdouble y2=1-(y-drawdata->yoff)/(gdouble)(drawdata->heightim);
		x2=(1.-x2)*drawdata->limit0[0]+x2*drawdata->limit0[1];
		y2=(1.-y2)*drawdata->limit0[2]+y2*drawdata->limit0[3];
		if(drawdata->xylog[0]!='n') x2=pow(10, x2);
		if(drawdata->xylog[1]!='n') y2=pow(10, y2);
		
		if(drawdata->p0){//2d image plot
			int xi=(int)((x2-drawdata->limit[0])/(drawdata->limit[1]-drawdata->limit[0])*drawdata->nx);
			int yi=(int)((y2-drawdata->limit[2])/(drawdata->limit[3]-drawdata->limit[2])*drawdata->ny);
			//convert to pixel coordinate
			x2=(xi+0.5)*(drawdata->limit[1]-drawdata->limit[0])/drawdata->nx+drawdata->limit[0];
			y2=(yi+0.5)*(drawdata->limit[3]-drawdata->limit[2])/drawdata->ny+drawdata->limit[2];
			
			if(xi>=0 && xi<drawdata->nx && yi>=0 && yi<drawdata->ny){
				float val=drawdata->p0[drawdata->nx*yi+xi];
				snprintf(drawdata->tooltip, sizeof(drawdata->tooltip), "(%g, %g): %g", x2, y2, val);
			}else{
				drawdata->tooltip[0]=0;
				//dbg_time("point (%d, %d) is invalid for array (%d,%d)\n", xi, yi, drawdata->nx, drawdata->ny);
			}
		}else{//line plot
			snprintf(drawdata->tooltip, sizeof(drawdata->tooltip), "(%g, %g)", x2, y2);
		}
		gtk_widget_set_tooltip_text(drawdata->drawarea, drawdata->tooltip);
	
	} else{
		if(drawdata->cursorinside){
			drawdata->cursorinside=0;
			//set_cursor(drawdata->drawarea, NULL);
			gtk_widget_set_has_tooltip(drawdata->drawarea, 0);
		}
	}
#if GTK_MAJOR_VERSION<4
	return FALSE;
#endif
}

static void do_zoom(drawdata_t* drawdata, float xdiff, float ydiff, int mode){
	if(!drawdata) return;
	float old_zoomx=drawdata->zoomx;
	float old_zoomy=drawdata->zoomy;
	if(mode==1){/*zoom in */
		drawdata->zoomx*=1.2;
		drawdata->zoomy*=1.2;
	} else if(mode==-1){
		drawdata->zoomx/=1.2;
		drawdata->zoomy/=1.2;
	} else if(mode==0){/*reset everything */
		drawdata->zoomx=1;
		drawdata->zoomy=1;
		drawdata->offx=0;
		drawdata->offy=0;
		drawdata->drawn=0;
	}
	if(drawdata->zoomx<MIN_ZOOM){
		drawdata->zoomx=MIN_ZOOM;
		drawdata->offx=0;
	}else if(drawdata->zoomx>MAX_ZOOM){
		drawdata->zoomx=MAX_ZOOM;
	}
	if(drawdata->zoomy<MIN_ZOOM){
		drawdata->zoomy=MIN_ZOOM;
		drawdata->offy=0;
	}else if(drawdata->zoomy>MAX_ZOOM){
		drawdata->zoomy=MAX_ZOOM;
	}
	//preserve visual center
	drawdata->offx=(drawdata->offx-xdiff)*drawdata->zoomx/old_zoomx+xdiff;
	drawdata->offy=(drawdata->offy-ydiff)*drawdata->zoomy/old_zoomy+ydiff;
	
	update_pixmap(drawdata);
}
#if GTK_MAJOR_VERSION>=4
static gboolean drawarea_scroll_event(GtkEventControllerScroll *scroll, gdouble dx, gdouble dy, drawdata_t *drawdata){
	GdkEvent *event=gtk_event_controller_get_current_event(GTK_EVENT_CONTROLLER(scroll));
	double x=drawdata->mxdown;
	double y=drawdata->mydown;

	unsigned int time=gdk_event_get_time(event);
	int zoom_dir=dy!=0?dy:dx;
	//info("dx=%g, dy=%g, x=%g, y=%g\n", dx, dy, x, y);
#else
static gboolean drawarea_scroll_event(GtkWidget* widget, GdkEventScroll* event, drawdata_t* drawdata){
	(void)widget;
	double x=event->x;
	double y=event->y;
	unsigned int time=event->time;
	int zoom_dir=(event->direction==GDK_SCROLL_UP)?1:-1;
	//dbg_time("scroll_event with dir %d\n", zoom_dir);
#endif
#define DO_ZOOM 1
#if DO_ZOOM
	static unsigned int last_time=0;

	if(time>last_time+100){//prevent fast scroll
		float xdiff=x?(x-drawdata->centerx):0;
		float ydiff=y?-(y-drawdata->centery):0;
		do_zoom(drawdata, xdiff, ydiff, zoom_dir);
		last_time=time;
	}
#else
	float xscale=0;
	float yscale=0;
	switch(event->direction){
	case GDK_SCROLL_UP:
		yscale=0.01;break;
	case GDK_SCROLL_DOWN:
		yscale=-0.01;break;
	case GDK_SCROLL_LEFT:
		xscale=-0.01;break;
	case GDK_SCROLL_RIGHT:
		xscale=0.01;break;
	}
	do_move(drawdata, drawdata->widthim*xscale, drawdata->heightim*yscale);
#endif
	return FALSE;
}

#if GTK_MAJOR_VERSION<4
static gboolean focus_in_handler(GtkWidget* widget, GdkEvent* event, drawdata_t* drawdata){
	(void)event;
	(void)widget;
	drawdata->region=0;
	//dbg_time("focus_in_handler.\n");
	return FALSE;
}

static gboolean drawarea_button_press(GtkWidget* widget, GdkEventButton* event, drawdata_t* drawdata){
	gdouble x=event->x;
	gdouble y=event->y;
	/*Grab focus so the keys work */
	if(!GTK_WIDGET_HAS_FOCUS(widget)) gtk_widget_grab_focus(widget);
#else
static gboolean drawarea_drag_begin(GtkGestureDrag*self, gdouble x, gdouble y, drawdata_t* drawdata){
	(void)self;
	//always triggered when mouse is clicked.
	//dbg_time("drag_begin at %g %g\n", x, y);
#endif
	if(x>drawdata->legbox_ox&&x<drawdata->legbox_ox+drawdata->legbox_width
		&&y>drawdata->legbox_oy&&y<drawdata->legbox_oy+drawdata->legbox_height){
		drawdata->region=2;
	}else if(x>drawdata->xoff&&x < drawdata->xoff+drawdata->widthim
		&&y > drawdata->yoff&&y<drawdata->yoff+drawdata->heightim){
		drawdata->region=1;
	} else{
		drawdata->region=0;
	}
	if(drawdata->region){
		drawdata->mxdown=x;
		drawdata->mydown=y;
		drawdata->dxdown=0;
		drawdata->dydown=0;
		drawdata->mtdown=myclockd();
	}
	//dbg_time("drawarea_button_press %g %g.\n", x, y);
	return FALSE;
}

#if GTK_MAJOR_VERSION<4
static gboolean drawarea_button_release(GtkWidget* widget, GdkEventButton* event, drawdata_t* drawdata){
	gdouble x=event->x;
	gdouble y=event->y;
	gdouble dx=x-drawdata->mxdown;
	gdouble dy=y-drawdata->mydown;
	gint button=event->button;
	(void)widget;
#else
static gboolean drawarea_drag_end(GtkGestureDrag *drag, gdouble dx, gdouble dy, drawdata_t* drawdata){
	gint button=gtk_gesture_single_get_current_button(GTK_GESTURE_SINGLE(drag));
	gdouble x=drawdata->mxdown+dx;
	gdouble y=drawdata->mydown+dy;
	//dbg_time("drag_end at %g %g\n", x, y);
#endif

	if(!drawdata->region) return FALSE;

	gdouble dt=myclockd()-drawdata->mtdown;
	dbg2_time("drawarea_button_release %g %g. dx=%g %g. button is %d.dt is %g\n",x, y, dx, dy, button, dt);
	if((fabs(dx)<3&&fabs(dy)<3)||dt<0.16){
		drawdata->draw_rect=0;
		gtk_widget_queue_draw(drawdata->drawarea);
		//dbg_time("Ignore accidental click\n");
	} else if(button==1){/*move only on left button */
		do_move(drawdata, (dx-drawdata->dxdown), -(dy-drawdata->dydown));
	} else if(button==3){/*right button select and zoom. */
		float xx=drawdata->mxdown;
		float yy=drawdata->mydown;
		if(drawdata->square){
			float ratio=1;
			if(drawdata->p){
				ratio=(float)drawdata->nx/(float)drawdata->ny;
			}
			if(fabs(dx)<fabs(dy)*ratio){
				dy*=fabs(dx/(dy*ratio));
			} else{
				dx*=fabs(dy*ratio/dx);
			}
		}
		if(dx<0) xx+=dx;
		if(dy>0) yy+=dy;
		float diffx=(drawdata->limit0[1]-drawdata->limit0[0])/drawdata->widthim;
		float diffy=(drawdata->limit0[3]-drawdata->limit0[2])/drawdata->heightim;
		drawdata->limit0[0]+=diffx*(xx-drawdata->xoff);
		drawdata->limit0[1]=drawdata->limit0[0]+diffx*fabs(dx);
		drawdata->limit0[2]+=diffy*(drawdata->yoff+drawdata->heightim-yy);
		drawdata->limit0[3]=drawdata->limit0[2]+diffy*fabs(dy);
		update_zoom(drawdata);
		update_pixmap(drawdata);
	}
	drawdata->region=0;
	return FALSE;

}


static void switch_tab(int lr, int ud){
	if(!curtopnb){
		dbg_time("curtopnb is null\n");
		return;
	}
	GtkWidget* topnb=curtopnb;
	if(lr){//switch top page
		gtk_notebook_set_current_page(GTK_NOTEBOOK(topnb),
			gtk_notebook_get_current_page(GTK_NOTEBOOK(topnb))+lr);
	} else if(ud){//switch sub page
		GtkWidget* subnb=gtk_notebook_get_nth_page(GTK_NOTEBOOK(topnb),
			gtk_notebook_get_current_page(GTK_NOTEBOOK(topnb)));
		gtk_notebook_set_current_page(GTK_NOTEBOOK(subnb),
			gtk_notebook_get_current_page(GTK_NOTEBOOK(subnb))+ud);
	}
}
#if GTK_MAJOR_VERSION<4
static gboolean drawarea_key_press(GtkWidget* widget, GdkEventKey* event, drawdata_t* drawdata){
	(void)widget;
	guint keyval=event->keyval;
	guint state=event->state;
#else
static gboolean drawarea_key_press(GtkEventControllerKey *ec, guint keyval, guint keycode, GdkModifierType state, drawdata_t*drawdata){
	(void) ec; (void) keycode;
	//info("drawarea_key_press with keyval %ud\n", keyval);
#endif
	if(state&GDK_CONTROL_MASK){
		switch(keyval){
		case GDK_Left:
			switch_tab(-1, 0);break;
		case GDK_Right:
			switch_tab(1, 0);break;
		case GDK_Up:
			switch_tab(0, -1);break;
		case GDK_Down:
			switch_tab(0, 1);break;
		default:
			return FALSE;
		}
	} else{
		switch(keyval){
		case GDK_plus:
		case GDK_equal:
			do_zoom(drawdata, 0, 0, 1); break;
		case GDK_minus:
			do_zoom(drawdata, 0, 0, -1);break;
		case GDK_0:
		case GDK_1:
			do_zoom(drawdata, 0, 0, 0);break;
		case GDK_Left:
			do_move(drawdata, -10, 0);break;
		case GDK_Right:
			do_move(drawdata, 10, 0);break;
		case GDK_Up:
			do_move(drawdata, 0, 10);break;
		case GDK_Down:
			do_move(drawdata, 0, -10);break;
		default:
			return FALSE;
		}
	}
	return TRUE;
}
//handles action when either top of sub label note book changed
//-1 means current page
static void page_changed(int topn, int subn){
	if(!curtopnb){
		dbg_time("page_changed: curtopnb is NULL\n");
		return;
	}
	GtkWidget* topnb=curtopnb;
	GtkWidget* subnb, * subpage;
	if(topn>-1){
		//info("page_changed 1, topnb=%p, topn=%d, subn=%d\n", topnb, topn, subn);
		subnb=gtk_notebook_get_nth_page(GTK_NOTEBOOK(topnb), topn);
	} else{
		subnb=get_current_page(topnb);
	}
	if(!subnb) return;

	if(subn>-1){
		//info("page_changed 2, subnb=%p\n", subnb);
		subpage=gtk_notebook_get_nth_page(GTK_NOTEBOOK(subnb), subn);
	} else{
		subpage=get_current_page(subnb);
	}
	if(!subpage) return;
	if(topn!=-1 || subn!=-1){
		drawdata_t* drawdata=(drawdata_t *)g_object_get_data(G_OBJECT(subpage), "drawdata");
		if(drawdata && !drawdata->drawn){
			update_pixmap(drawdata);
		}
		update_toolbar();
	}else{
		static int client_pid_last=0;
		if(client_pid==client_pid_last){
			return;//avoid sending repeatedly to the same client
		}
		client_pid_last=client_pid;
	}
	if(sock!=-1&&client_pid>0){
		const char *fig=topnb_label_get(topnb, subnb);
		const char *fn=subnb_label_get(subnb, subpage);
		dbg("send fig=%s, fn=%s\n", fig, fn);
		if(stwriteint(sock, DRAW_FIGFN)||
			stwritestr(sock, fig)||
			stwritestr(sock, fn)){
			dbg_time("Talk to client failed\n");
			close(sock);
			sock=-1;
		}
	}else{
		//dbg("cannot send fig=%s, fn=%s, sock=%d, client_pid=%d\n", fig, fn, sock, client_pid);
	}
	//info("done\n");
	io_time2=0;
}
/*These signal handlers are called before the notebook page switch is done.*/
static void topnb_page_switch(GtkNotebook* topnb, GtkWidget* page, guint n, GtkWidget* toolbar){
	(void)topnb; (void)page; (void) toolbar;
	//gtk_widget_set_sensitive(toolbar, TRUE);
	//dbg_time("topnb_page_switch: %u\n", n);
	page_changed(n, -1);
}
static void subnb_page_switch(GtkNotebook* subnb, GtkWidget* page, guint n, gpointer dummy){
	//dbg_time("subnb_page_switch: %u\n", n);
	(void)subnb; (void)page;
	(void)dummy;
	page_changed(-1, n);
}
static void subnb_page_removed(GtkNotebook *subnb, GtkWidget *drawarea, guint n, GtkWidget *topnb){
	//dbg_time("subnb_page_removed, subnb=%p, child=%p, n=%u\n",  subnb, drawarea, n);
	(void)n; (void)drawarea;
	if(gtk_notebook_get_n_pages(subnb)==0){
		//GtkWidget *topnb=gtk_widget_get_parent(GTK_WIDGET(subnb));//not work in gtk4
		if(topnb){
			int ipage=gtk_notebook_page_num(GTK_NOTEBOOK(topnb), GTK_WIDGET(subnb));
			if(ipage==-1){
				//dbg_time("page not found\n");//happens if topnb is destroyed already
			}else{
				gtk_notebook_remove_page(GTK_NOTEBOOK(topnb), ipage);
			}
		}
	}
}
//close additional window if there are no more pages left
static void topnb_page_removed(GtkNotebook* topnb, GtkWidget* child, guint n, GtkWidget* toolbar){
	//info("topnb_page_removed with n=%d for %p\n", n, topnb);
	(void)child;
	(void)n;
	(void)toolbar;
	int npage=gtk_notebook_get_n_pages(topnb);
	if(npage==0){/*no more pages left. */
		//gtk_widget_set_sensitive(toolbar, FALSE);
		//info("window list length is %d\n", g_slist_length(windows));
		if(g_slist_length(windows)>1){
			GtkWidget* window=gtk_widget_get_parent(gtk_widget_get_parent(GTK_WIDGET(topnb)));
			if(window){
				window_destroy(window);
			}
		}
	}
	/*gtk_notebook_set_show_tabs(topnb, npage!=1); */
}
static void topnb_page_added(GtkNotebook *topnb, GtkWidget *child, guint n, GtkWidget *toolbar){
	//info("topnb_page_added with n=%d for %p\n", n, topnb);
	(void)topnb; (void)child;
	gtk_widget_set_sensitive(toolbar, TRUE);
	page_changed(n, -1);
	/*gtk_notebook_set_show_tabs(topnb, npage!=1); */
}

/**
   2012-10-27: GTK 3.6 deprecated gdk_threads_enter(). So it is hard to call
   addpage from child threads. Modify the routine so that addpage is called by
   the main thread when gtk is idle.
*/
gboolean addpage(drawdata_t *drawdata){
	if(!drawdata->fig || !drawdata->name){
		dbg_time("Must set fig and name before calling addpage");
		return 0;
	}
	if(__atomic_load_n(&drawdata->recycle, __ATOMIC_SEQ_CST)){
		//this is caused by race condition
		//When the GUI is removing the page, the drawdaemon_io was receiving the data and then queued a addpage
		dbg_time("unexpected: drawdata(%p). recycle is set\n", drawdata);
		return 0;
	}
	//info("Adding page %s, %s\n", drawdata->fig, drawdata->name);
	GSList* subnbs=NULL;
	int nsubnb=0;
	GtkWidget* window=0;
	GtkWidget* topnb=0;
	int jtab=-1;//default: insert to the end
	if(!windows){//not windows remaining
		create_window(NULL);
	}
	for(GSList *p=windows; p; p=p->next){
		//scan through all window to find all topnb page that has the same "fig"
		window=GTK_WIDGET(p->data);
		topnb=get_topnb(window);
		for(int itab=0; itab<gtk_notebook_get_n_pages(GTK_NOTEBOOK(topnb)); itab++){
			GtkWidget* subnb=gtk_notebook_get_nth_page(GTK_NOTEBOOK(topnb), itab);
			int res=strcmp(drawdata->fig, topnb_label_get(topnb, subnb));
			if(!res){//found
				subnbs=g_slist_append(subnbs, subnb);
				nsubnb++;/*number of subnbs find. */
			} else if(res<0 && jtab==-1){//not found.
				jtab=itab;//mark insert location
			}
		}
	}

	if(!nsubnb){/*subnb not found. create one. */
		GtkWidget* subnb=gtk_notebook_new();
		subnbs=g_slist_append(subnbs, subnb);
		nsubnb++;
		//gtk_container_set_border_width(GTK_CONTAINER(subnb),0);
#if GTK_VERSION_AFTER(2, 24)
		gtk_notebook_set_group_name(GTK_NOTEBOOK(subnb), "secondlevel");//since 2.24
#elif GTK_VERSION_AFTER(2, 12)
		gtk_notebook_set_group(GTK_NOTEBOOK(subnb), "secondlevel");//since 2.12
#endif
		gtk_notebook_set_tab_pos(GTK_NOTEBOOK(subnb), GTK_POS_RIGHT);
		gtk_notebook_set_scrollable(GTK_NOTEBOOK(subnb), TRUE);
		//we use connect_after so that get_current_page gives right answer
		g_signal_connect_after(subnb, "switch-page", G_CALLBACK(subnb_page_switch), NULL);
		g_signal_connect(subnb, "page-removed", G_CALLBACK(subnb_page_removed), topnb);
		GtkWidget* label=gtk_label_new(drawdata->fig);
		gtk_notebook_insert_page(GTK_NOTEBOOK(topnb), subnb, label, jtab);//todo: restore tab_button_cb
		//gtk_notebook_set_tab_detachable(GTK_NOTEBOOK(topnb), subnb, TRUE);
		gtk_notebook_set_tab_reorderable(GTK_NOTEBOOK(topnb), subnb, TRUE);
	}
	GtkWidget* drawarea=NULL;
	GtkWidget* subnb=(GtkWidget*)subnbs->data;
	jtab=-1;
	for(GSList* p=subnbs; p; p=p->next){
		/*scan through all the subnb pages with same label */
		//info("subnb=%p\n", subnb);
		for(int itab=0; itab<gtk_notebook_get_n_pages(GTK_NOTEBOOK(p->data)); itab++){
			GtkWidget* tmp=gtk_notebook_get_nth_page(GTK_NOTEBOOK(p->data), itab);
			drawdata_t* drawdata_old=(drawdata_t*)g_object_get_data(G_OBJECT(tmp), "drawdata");
			if(drawdata_old==drawdata){//We match using drawdata pointer, not ->name to avoid race condition when tab is deleted.
				drawarea=drawdata->drawarea;
				break;
			}else{//locate a insert position
				const gchar* labeltext=subnb_label_get(GTK_WIDGET(p->data), tmp);
				int res=strcmp(drawdata->name, labeltext);
				if(res<0 && jtab==-1){
					jtab=itab;
					subnb=(GtkWidget*)p->data;
				}
			}
		}
	}
	if(drawarea){
		update_pixmap(drawdata);/*update pixmap */
		if(get_current_drawdata()!=drawdata){
			page_changed(-1, -1);/*notify client that it is not drawing to active page */
		}
	} else{
		//dbg_time("%p: addpage\n", drawdata);
		/*new tab inside the fig to contain the plot. */
		drawarea=drawdata->drawarea=gtk_drawing_area_new();
		drawdata->subnb=subnb;
		g_object_set_data_full(G_OBJECT(drawarea), "drawdata", drawdata, drawarea_deleted);

#if GTK_VERSION_AFTER(3, 22)
		gtk_widget_set_can_focus(drawarea, TRUE);//2.18
		gtk_widget_set_sensitive(drawarea, TRUE);
		//gtk_widget_set_focusable(drawarea, TRUE);//needed to receive keyboard events
		gtk_widget_set_focus_on_click(drawarea, TRUE);//3.20
#else
		GTK_WIDGET_SET_FLAGS(drawarea, GTK_CAN_FOCUS);
		GTK_WIDGET_SET_FLAGS(drawarea, GTK_SENSITIVE);
#endif
#if GTK_MAJOR_VERSION>=4
		GtkGesture *drag=gtk_gesture_drag_new();
		gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(drag), 0);
		gtk_widget_add_controller(drawarea, GTK_EVENT_CONTROLLER(drag));
		g_signal_connect(drag, "drag-begin", G_CALLBACK(drawarea_drag_begin), drawdata);
		g_signal_connect(drag, "drag-update", G_CALLBACK(drawarea_drag_update), drawdata);
		g_signal_connect(drag, "drag-end", G_CALLBACK(drawarea_drag_end), drawdata);

		GtkEventController *ec=gtk_event_controller_key_new();
		gtk_widget_add_controller(drawarea, GTK_EVENT_CONTROLLER(ec));
		g_signal_connect(ec, "key-pressed", G_CALLBACK(drawarea_key_press), drawdata);//not useful

		GtkEventController *scroll=gtk_event_controller_scroll_new(GTK_EVENT_CONTROLLER_SCROLL_VERTICAL);
		gtk_widget_add_controller(drawarea, GTK_EVENT_CONTROLLER(scroll));
		g_signal_connect(scroll, "scroll", G_CALLBACK(drawarea_scroll_event), drawdata);

		GtkEventController* motion=gtk_event_controller_motion_new();
		gtk_widget_add_controller(drawarea, GTK_EVENT_CONTROLLER(motion));
		g_signal_connect(motion, "motion", G_CALLBACK(drawarea_motion_event), drawdata);
#else
		gtk_widget_add_events(drawarea, GDK_BUTTON_PRESS_MASK|
			GDK_BUTTON_RELEASE_MASK|
			GDK_POINTER_MOTION_MASK|GDK_POINTER_MOTION_HINT_MASK|
			GDK_BUTTON1_MOTION_MASK|
			GDK_KEY_PRESS_MASK|
			GDK_KEY_RELEASE_MASK);
		g_signal_connect(drawarea, "motion-notify-event",
			G_CALLBACK(drawarea_motion_notify), drawdata);

		//notice that there is no drawarea_button_press event when click through to activate a window.
		//so we also connect focus-in-event
		g_signal_connect(drawarea, "button-press-event",
			G_CALLBACK(drawarea_button_press), drawdata);
		g_signal_connect(drawarea, "button-release-event",
			G_CALLBACK(drawarea_button_release), drawdata);
		g_signal_connect(drawarea, "focus-in-event",
			G_CALLBACK(focus_in_handler), drawdata);

		g_signal_connect(drawarea, "scroll-event",
			G_CALLBACK(drawarea_scroll_event), drawdata);

   		/*g_signal_connect(drawarea,"button-release-event",
	 	G_CALLBACK(drawarea_button_release),drawdata);*/
		g_signal_connect(drawarea, "key-press-event",
			G_CALLBACK(drawarea_key_press), drawdata);

#endif
		GtkWidget* button=subnb_label_new(drawdata->name, drawdata);
		gtk_notebook_insert_page(GTK_NOTEBOOK(subnb), drawarea, button, jtab);
		gtk_notebook_set_tab_reorderable(GTK_NOTEBOOK(subnb), drawarea, TRUE);
#if GTK_MAJOR_VERSION>=4
		gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(drawarea), drawarea_draw_func, drawdata, NULL);
#elif GTK_MAJOR_VERSION==3
		g_signal_connect(drawarea, "draw",G_CALLBACK(on_draw_event), drawdata);
#else
		g_signal_connect(drawarea, "expose-event",G_CALLBACK(on_expose_event), drawdata);
#endif
	/*handles zooming. */
#if GTK_MAJOR_VERSION>=4
		g_signal_connect(drawarea, "resize",G_CALLBACK(on_resize_event), drawdata);
#else
		g_signal_connect(drawarea, "configure-event", G_CALLBACK(on_configure_event), drawdata);

#endif
		gtk_widget_set_size_request(drawarea, DRAWAREA_MIN_WIDTH, DRAWAREA_MIN_HEIGHT);
#if GTK_MAJOR_VERSION<4
		gtk_widget_show_all(subnb);
#endif
		if(get_current_drawdata()!=drawdata){
			//The new page is not activated. notify client that it is not drawing to active page
			page_changed(-1, -1);
		}
	}
	g_slist_free(subnbs);
	return 0;//return 0 cause the function to be removed from g_idle_add()
}
static void save_file(drawdata_t *drawdata);
#if GTK_VERSION_AFTER(4,10)
static void file_saved(GObject *dialog, GAsyncResult *result, void *data){
	drawdata_t *drawdata=(drawdata_t*)data;
	(void) result;
	if(drawdata->file) g_object_unref(drawdata->file);
	drawdata->file=gtk_file_dialog_save_finish(GTK_FILE_DIALOG(dialog), result, NULL);//do not free the result
	if(drawdata->file){
		if(drawdata->filename) g_free(drawdata->filename);
		drawdata->filename=g_file_get_path(drawdata->file);
		save_file(drawdata);
	}
}
#endif
static void tool_save(GtkWidget* button, gpointer data){
	drawdata_t *drawdata=get_current_drawdata();
	if(!drawdata){
		dbg_time("drawdata is empty\n"); return;
	}
	(void)button;(void)data;
#if GTK_VERSION_AFTER(4,10)
	GtkFileDialog *dialog=gtk_file_dialog_new();
	gtk_file_dialog_set_title(GTK_FILE_DIALOG(dialog), "Select a file to save");
	gtk_file_dialog_set_modal(GTK_FILE_DIALOG(dialog), true);
	if(!drawdata->file){
		drawdata->file=g_file_new_for_path(client_path_full);
	}
	gtk_file_dialog_set_initial_file(GTK_FILE_DIALOG(dialog), drawdata->file);
	gtk_file_dialog_save(GTK_FILE_DIALOG(dialog), GTK_WINDOW(curwindow), NULL, file_saved, drawdata);

#else
	GtkWidget *dialog=gtk_file_chooser_dialog_new
	("Select a file to save",
		GTK_WINDOW(curwindow),
		GTK_FILE_CHOOSER_ACTION_SAVE,
		"_Cancel", GTK_RESPONSE_CANCEL,
		"_Save", GTK_RESPONSE_ACCEPT,
		NULL);
	gtk_file_chooser_set_do_overwrite_confirmation(GTK_FILE_CHOOSER(dialog), TRUE);
	if(drawdata->filename){
		gtk_file_chooser_set_filename(GTK_FILE_CHOOSER(dialog), drawdata->filename);
	} else if(client_path_full){
		gtk_file_chooser_set_current_folder
		(GTK_FILE_CHOOSER(dialog), client_path_full);
	} else{
		char curpath[PATH_MAX];
		if(!getcwd(curpath, PATH_MAX)){
			strcpy(curpath, HOME);
		}
		gtk_file_chooser_set_current_folder
		(GTK_FILE_CHOOSER(dialog), curpath);
	}

	if(gtk_dialog_run(GTK_DIALOG(dialog))==GTK_RESPONSE_ACCEPT){
		if(client_path_full) free(client_path_full);
		client_path_full=gtk_file_chooser_get_current_folder(GTK_FILE_CHOOSER(dialog));
		if(drawdata->filename) g_free(drawdata->filename);
		drawdata->filename=gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
		window_destroy(dialog);
	} else{
		window_destroy(dialog);
		return;
	}
	save_file(drawdata);
#endif
}
static void save_file(drawdata_t *drawdata){
	if(!drawdata){
		dbg_time("drawdata is empty\n");
		return;
	}
		const char *filename=drawdata->filename;
	if(!filename){
		dbg_time("filename is empty\n");
		return;
	}else{
		info("filename is %s\n", filename);
	}

	const char* suffix=strrchr(filename, '.');
	if(!suffix){
		char* filename2=stradd(filename, ".png", NULL);
		g_free(drawdata->filename);
		filename=drawdata->filename=filename2;
		suffix=strrchr(filename, '.');
	}

	if(strcmp(suffix, ".eps")==0){
#if CAIRO_HAS_PS_SURFACE == 1
		drawdata->surface=cairo_ps_surface_create(filename, drawdata->width, drawdata->height);
		cairo_ps_surface_set_eps(drawdata->surface, TRUE);
#endif
	} else if(strcmp(suffix, ".pdf")==0){
#if CAIRO_HAS_PDF_SURFACE == 1
		drawdata->surface=cairo_pdf_surface_create(filename, drawdata->width, drawdata->height);
#endif
	} else if(strcmp(suffix, ".svg")==0){
#if CAIRO_HAS_SVG_SURFACE == 1
		drawdata->surface=cairo_svg_surface_create(filename, drawdata->width, drawdata->height);
#endif
	} else if(strcmp(suffix, ".png")==0){
#if CAIRO_HAS_PNG_FUNCTIONS
#if GTK_MAJOR_VERSION>=3
		cairo_surface_write_to_png(drawdata->pixmap, filename);
		return;
#else
		float scale=2;
		drawdata->surface=cairo_image_surface_create((cairo_format_t)CAIRO_FORMAT_ARGB32, drawdata->width*scale, drawdata->height*scale);
		cairo_surface_set_device_scale(drawdata->surface, scale, scale);
#endif		
#endif		
	}else if(strcmp(suffix, ".gif")==0){
#if CAIRO_HAS_PNG_FUNCTIONS
#if GTK_MAJOR_VERSION>=3
		drawdata->filename_gif=strndup(filename, strlen(filename)-4);
		mymkdir("%s", drawdata->filename_gif);
		pfilename_gif=&drawdata->filename_gif;
		//drawdata->zlim_manual=1;
		drawdata->frame_io=0;
		return;
#endif
#endif
	}
	if(!drawdata->surface){
		error_msg("%s: file type is not supported.\n", filename);
		return;
	}
	cairo_draw(drawdata);
}

static void tool_zoom(GtkWidget* button, gpointer data){
	(void)button;
	drawdata_t* drawdata=get_current_drawdata();
	int mode=GPOINTER_TO_INT(data)-1;
	do_zoom(drawdata, 0, 0, mode);
}

static void limit_change(GtkSpinButton* spin, gfloat* val){
	*val=gtk_spin_button_get_value(spin);
	update_zoom(drawdata_dialog);
	update_pixmap(drawdata_dialog);
}

static void zlim_changed(GtkSpinButton* spin, gfloat* val){
	*val=gtk_spin_button_get_value(spin);
	drawdata_dialog->zlim_changed=1;
	drawdata_dialog->zlim_manual=2;//1: set by draw(). 2: set by gui.
	drawdata_dialog->drawn=0;
	update_pixmap(drawdata_dialog);
}
static void checkbtn_toggle(GtkWidget* btn, gint* key){
	*key=check_button_get_active(btn);
	drawdata_dialog->drawn=0;
	update_pixmap(drawdata_dialog);
}
static void checkbtn_toggle_inv(GtkWidget *btn, gint *key){
	*key=!check_button_get_active(btn);
	drawdata_dialog->drawn=0;
	update_pixmap(drawdata_dialog);
}
static void checkbtn_toggle_char(GtkWidget *btn, char *key){
	*key=check_button_get_active(btn)?'y':'n';
	drawdata_dialog->drawn=0;
	update_pixmap(drawdata_dialog);
}
/*static void range_changed(GtkRange* range, gfloat* val){
	*val=gtk_range_get_value(range);
	update_pixmap(drawdata_dialog);
}*/
static void entry_changed(GtkEditable* entry, char** key){
	free(*key);
	*key=gtk_editable_get_chars(entry, 0, -1);
	drawdata_dialog->drawn=0;
	update_pixmap(drawdata_dialog);
}
static void spin_changed(GtkSpinButton* spin, gfloat* val){
	*val=gtk_spin_button_get_value(spin);
	drawdata_dialog->drawn=0;
	update_pixmap(drawdata_dialog);
}
static void spin_icumu(GtkSpinButton* spin){
	drawdata_t *drawdata=get_current_drawdata();
	if(drawdata && !drawdata->p&&!drawdata->square && drawdata->cumu){
		drawdata->icumu=gtk_spin_button_get_value(spin);
		update_pixmap(drawdata);
		//dbg("set %p to %d\n", &drawdata->cumu, drawdata->cumu);
	}
}
static void togglebutton_toggle(GtkWidget* btn, void *data){
	int key=GPOINTER_TO_INT(data);
	int val=gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(btn));
	drawdata_t *drawdata=get_current_drawdata();
	if(!drawdata) return;
	switch(key){
      case 1:
		drawdata->zlog=val;
		break;
	  case 2:
		drawdata->cumu=val;
		break;
	}
	update_toolbar();
	update_pixmap(drawdata);
}

static void togglebutton_stop(GtkWidget* btn, gpointer data){
	(void)btn; (void) data;
	if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(btn))){
		keep_listen=0;//disable reconnection
		close(sock); sock=-1;//close connection
		mysleep(0.001);//wait for listen_draw to quit.
	}else{
		if(!keep_listen){
			keep_listen=1;//only modify keep_listen in this thread to avoid race condition.
			thread_new(listen_draw, NULL);
		}
	}
}
static void togglebutton_pause(GtkWidget* btn, gpointer data){
	(void) data;
	if(sock==-1) return;
	int cmd=gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(btn))?DRAW_PAUSE:DRAW_RESUME;
	if(stwriteint(sock, cmd)){
		close(sock);
		sock=-1;
	}
}
static void toolbutton_refresh(GtkWidget* btn, gpointer data){
	(void)btn; (void) data;
	if(sock==-1){
		
	}else if(stwriteint(sock, DRAW_SINGLE)){//toggles DRAW_SINGLE in draw()
		close(sock);
		sock=-1;
	}
}
static void set_window_size(GtkWidget *spin, gpointer data){
	int type=GPOINTER_TO_INT(data);
	int width, height;
#if GTK_MAJOR_VERSION<4	
	gtk_window_get_size(GTK_WINDOW(curwindow), &width, &height);
#else
	width=gtk_widget_get_width(curwindow);
	height=gtk_widget_get_height(curwindow);
#endif	
	double ratio=(double)width/height;
	switch(type){
		case 1:
			width=gtk_spin_button_get_value(GTK_SPIN_BUTTON(spin));	
			if(window_proportional){
				height=(int)round(width/ratio);
			}
			break;
		case 2:
			height=gtk_spin_button_get_value(GTK_SPIN_BUTTON(spin)); 
			if(window_proportional){
				width=(int)round(height*ratio);
			}
			break;
		default:
		break;
	}
	
#if GTK_MAJOR_VERSION<4	
	gtk_window_resize(GTK_WINDOW(curwindow), width, height);
#else
	gtk_widget_set_size_request(curwindow, width, height);
#endif	
}
/**
   Response to the quest to set the zaxis limit (the range of the color bar)
*/
static void tool_property(GtkWidget* button, gpointer data){
	(void)button;
	(void)data;
	drawdata_t* drawdata=get_current_drawdata();
	if(!drawdata){
		return;
	}
	drawdata_dialog=drawdata;
#if GTK_VERSION_AFTER(4,10)
	GtkWidget* dialog=window_new();
	gtk_window_set_modal(GTK_WINDOW(dialog), TRUE);
	gtk_window_set_transient_for(GTK_WINDOW(dialog), GTK_WINDOW(curwindow));//behaves like a dialog
	//GtkWidget *popover=gtk_popover_new();
	//gtk_popover_set_autohide(GTK_POPOVER(popover), TRUE);
#else
	GtkWidget* dialog=gtk_dialog_new_with_buttons("Figure Properties",
		GTK_WINDOW(curwindow),
		(GtkDialogFlags)(GTK_DIALOG_DESTROY_WITH_PARENT),
		"_OK", GTK_RESPONSE_ACCEPT,
		NULL);
	GtkWidget *content_area=gtk_dialog_get_content_area(GTK_DIALOG(dialog));
#endif
	GtkWidget* vbox=gtk_vbox_new(FALSE, 0);
	GtkWidget* hbox;
	GtkWidget* checkbtn, * label, * entry, * spin;
	int n;
	GtkWidget* spins[6];
	float diff[3];
	diff[0]=(drawdata->limit0[1]-drawdata->limit0[0]);
	diff[1]=(drawdata->limit0[3]-drawdata->limit0[2]);
	float *zlim=drawdata->zlim+(drawdata->zlog?2:0);
	if(zlim[0]||zlim[1]){
		diff[2]=(zlim[1]-zlim[0]);
		n=6;
	} else{
		n=4;
	}
	for(int i=0; i<n; i++){
	/*divide the separation to 100 steps. */
		if(diff[i/2]<EPS){
			diff[i/2]=1;
		}
		float step=pow(10, round(log10(fabs(diff[i/2])))-2);
		spins[i]=gtk_spin_button_new_with_range(drawdata->limit0[i]-1000*step, drawdata->limit0[i]+1000*step, step);
		//info("diff is %g, step is %g, limit is %g\n", diff[i/2], step, drawdata->limit0[i]);
		float* val;
		if(i<4){
			val=&drawdata->limit0[i];
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(spins[i]), *val);
			g_signal_connect(spins[i], "value-changed", G_CALLBACK(limit_change), val);
		} else{
			val=zlim+(i-4);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(spins[i]), *val);
			g_signal_connect(spins[i], "value-changed", G_CALLBACK(zlim_changed), val);
		}
	}
	drawdata->spins=spins;
	hbox=gtk_hbox_new(FALSE, 0);
	checkbtn=gtk_check_button_new_with_label("Make image square");
	check_button_set_active(checkbtn, drawdata->square);
	g_signal_connect(checkbtn, "toggled", G_CALLBACK(checkbtn_toggle), &drawdata->square);
	box_append(GTK_BOX(hbox), checkbtn, FALSE, FALSE, 0);

	checkbtn=gtk_check_button_new_with_label("Show grids");
	check_button_set_active(checkbtn, drawdata->grid);
	g_signal_connect(checkbtn, "toggled", G_CALLBACK(checkbtn_toggle), &drawdata->grid);
	box_append(GTK_BOX(hbox), checkbtn, FALSE, FALSE, 0);

	checkbtn=gtk_check_button_new_with_label("Put tic inside");
	check_button_set_active(checkbtn, drawdata->ticinside);
	g_signal_connect(checkbtn, "toggled", G_CALLBACK(checkbtn_toggle), &drawdata->ticinside);
	box_append(GTK_BOX(hbox), checkbtn, FALSE, FALSE, 0);

	box_append(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
	



	hbox=gtk_hbox_new(FALSE, 0);

	checkbtn=gtk_check_button_new_with_label("Legend");
	check_button_set_active(checkbtn, drawdata->legendbox);
	gtk_widget_set_sensitive(checkbtn, drawdata->legend!=NULL);
	g_signal_connect(checkbtn, "toggled", G_CALLBACK(checkbtn_toggle), &drawdata->legendbox);
	box_append(GTK_BOX(hbox), checkbtn, FALSE, FALSE, 0);

	checkbtn=gtk_check_button_new_with_label("On curve");
	check_button_set_active(checkbtn, drawdata->legendcurve);
	g_signal_connect(checkbtn, "toggled", G_CALLBACK(checkbtn_toggle), &drawdata->legendcurve);
	//box_append(GTK_BOX(vbox), checkbtn, FALSE, FALSE, 0);
	box_append(GTK_BOX(hbox), checkbtn, FALSE, FALSE, 0);
	
	checkbtn=gtk_check_button_new_with_label("Ellipsis");
	check_button_set_active(checkbtn, !noellipsis);
	gtk_widget_set_sensitive(checkbtn, drawdata->legend!=NULL);
	g_signal_connect(checkbtn, "toggled", G_CALLBACK(checkbtn_toggle_inv), &noellipsis);
	box_append(GTK_BOX(hbox), checkbtn, FALSE, FALSE, 0);

	box_append(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
	/*checkbtn=gtk_hscale_new_with_range(0, 1, 0.01);
	//gtk_widget_set_size_request(checkbtn, 80, 20);
	gtk_scale_set_draw_value(GTK_SCALE(checkbtn), 0);
	gtk_range_set_value(GTK_RANGE(checkbtn), drawdata->legendoffx);
	g_signal_connect(checkbtn, "value-changed", G_CALLBACK(range_changed), &drawdata->legendoffx);
	box_append(GTK_BOX(hbox), gtk_label_new("Position: H  "), FALSE, FALSE, 0);
	box_append(GTK_BOX(hbox), checkbtn, TRUE, TRUE, 0);

	checkbtn=gtk_hscale_new_with_range(0, 1, 0.01);
	//gtk_widget_set_size_request(checkbtn, 80, 20);
	gtk_scale_set_draw_value(GTK_SCALE(checkbtn), 0);
	gtk_range_set_value(GTK_RANGE(checkbtn), drawdata->legendoffy);
	g_signal_connect(checkbtn, "value-changed", G_CALLBACK(range_changed), &drawdata->legendoffy);
	box_append(GTK_BOX(hbox), gtk_label_new("  V  "), FALSE, FALSE, 0);
	box_append(GTK_BOX(hbox), checkbtn, TRUE, TRUE, 0);
	box_append(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
*/
	//set log scale
	hbox=gtk_hbox_new(FALSE, 0);
	//label=gtk_label_new("Set log scale:");
	//box_append(GTK_BOX(hbox), label, FALSE, FALSE, 0);
	checkbtn=gtk_check_button_new_with_label("log(X)");
	check_button_set_active(checkbtn, drawdata->xylog[0]=='y'?1:0);
	gtk_widget_set_sensitive(checkbtn, !drawdata->square);
	g_signal_connect(checkbtn, "toggled", G_CALLBACK(checkbtn_toggle_char), &drawdata->xylog[0]);
	box_append(GTK_BOX(hbox), checkbtn, FALSE, FALSE, 0);

	checkbtn=gtk_check_button_new_with_label("log(Y)");
	check_button_set_active(checkbtn, drawdata->xylog[1]=='y'?1:0);
	gtk_widget_set_sensitive(checkbtn, !drawdata->square);
	g_signal_connect(checkbtn, "toggled", G_CALLBACK(checkbtn_toggle_char), &drawdata->xylog[1]);
	box_append(GTK_BOX(hbox), checkbtn, FALSE, FALSE, 0);

	checkbtn=gtk_check_button_new_with_label("log(Value)");
	check_button_set_active(checkbtn, drawdata->zlog);
	gtk_widget_set_sensitive(checkbtn, drawdata->p?TRUE:FALSE);
	g_signal_connect(checkbtn, "toggled", G_CALLBACK(checkbtn_toggle), &drawdata->zlog);
	box_append(GTK_BOX(hbox), checkbtn, FALSE, FALSE, 0);

	checkbtn=gtk_check_button_new_with_label("Show colorbar");
	check_button_set_active(checkbtn, !hide_colorbar&&drawdata->p);
	gtk_widget_set_sensitive(checkbtn, drawdata->p?1:0);
	g_signal_connect(checkbtn, "toggled", G_CALLBACK(checkbtn_toggle_inv), &hide_colorbar);
	box_append(GTK_BOX(hbox), checkbtn, FALSE, FALSE, 0);
	
	box_append(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

	hbox=gtk_hbox_new(FALSE, 0);
	checkbtn=gtk_check_button_new_with_label("Cumulative average (");
	gtk_widget_set_sensitive(checkbtn, (drawdata->npts>0));
	check_button_set_active(checkbtn, drawdata->cumu);
	g_signal_connect(checkbtn, "toggled", G_CALLBACK(checkbtn_toggle), &drawdata->cumu);
	box_append(GTK_BOX(hbox), checkbtn, FALSE, FALSE, 0);

	checkbtn=gtk_check_button_new_with_label("overlay ");
	gtk_widget_set_sensitive(checkbtn, (drawdata->npts>0));
	check_button_set_active(checkbtn, drawdata->cumuover);
	g_signal_connect(checkbtn, "toggled", G_CALLBACK(checkbtn_toggle), &drawdata->cumuover);
	box_append(GTK_BOX(hbox), checkbtn, FALSE, FALSE, 0);

	checkbtn=gtk_check_button_new_with_label("quadrature) from ");
	gtk_widget_set_sensitive(checkbtn, (drawdata->npts>0));
	check_button_set_active(checkbtn, drawdata->cumuquad);
	g_signal_connect(checkbtn, "toggled", G_CALLBACK(checkbtn_toggle), &drawdata->cumuquad);
	box_append(GTK_BOX(hbox), checkbtn, FALSE, FALSE, 0);

	spin=gtk_spin_button_new_with_range(0, drawdata->limit[1], 1);
	gtk_widget_set_sensitive(spin, (drawdata->npts>0));
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin), drawdata->icumu);
	g_signal_connect(spin, "value-changed", G_CALLBACK(spin_changed), &drawdata->icumu);
	box_append(GTK_BOX(hbox), spin, TRUE, TRUE, 0);
	box_append(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

	//Set update low pass filter
	int lwidth=10;
	hbox=gtk_hbox_new(FALSE, 0);
	label=gtk_label_new("Update LPF");gtk_label_set_width_chars(GTK_LABEL(label), lwidth);
	//gtk_label_set_xalign(GTK_LABEL(label), 0);//gtk>=3.16
	spin=gtk_spin_button_new_with_range(0, 1, 0.01);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin), lpf);
	g_signal_connect(spin, "value-changed", G_CALLBACK(spin_changed), &lpf);
	box_append(GTK_BOX(hbox), label, FALSE, FALSE, 0);
	box_append(GTK_BOX(hbox), spin, TRUE, TRUE, 0);
	box_append(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

	hbox=gtk_hbox_new(FALSE, 0);
	entry=gtk_entry_new();
	checkbtn=gtk_check_button_new_with_label("Title");
	check_button_set_active(checkbtn, !hide_title);
	g_signal_connect(checkbtn, "toggled", G_CALLBACK(checkbtn_toggle_inv), &hide_title);
	box_append(GTK_BOX(hbox), checkbtn, FALSE, FALSE, 0);
	box_append(GTK_BOX(hbox), entry, TRUE, TRUE, 0);
	entry_set_text(GTK_ENTRY(entry), drawdata->title);
	g_signal_connect(GTK_EDITABLE(entry), "changed", G_CALLBACK(entry_changed), &drawdata->title);
	box_append(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

	hbox=gtk_hbox_new(FALSE, 0);
	entry=gtk_entry_new();
	checkbtn=gtk_check_button_new_with_label("X label");
	check_button_set_active(checkbtn, !hide_xlabel);
	g_signal_connect(checkbtn, "toggled", G_CALLBACK(checkbtn_toggle_inv), &hide_xlabel);
	box_append(GTK_BOX(hbox), checkbtn, FALSE, FALSE, 0);
	box_append(GTK_BOX(hbox), entry, TRUE, TRUE, 0);
	entry_set_text(GTK_ENTRY(entry), drawdata->xlabel);
	g_signal_connect(GTK_EDITABLE(entry), "changed", G_CALLBACK(entry_changed), &drawdata->xlabel);
	box_append(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

	hbox=gtk_hbox_new(FALSE, 0);
	entry=gtk_entry_new();
	checkbtn=gtk_check_button_new_with_label("Y label");
	check_button_set_active(checkbtn, !hide_ylabel);
	g_signal_connect(checkbtn, "toggled", G_CALLBACK(checkbtn_toggle_inv), &hide_ylabel);
	box_append(GTK_BOX(hbox), checkbtn, FALSE, FALSE, 0);
	box_append(GTK_BOX(hbox), entry, TRUE, TRUE, 0);
	entry_set_text(GTK_ENTRY(entry), drawdata->ylabel);
	g_signal_connect(GTK_EDITABLE(entry), "changed", G_CALLBACK(entry_changed), &drawdata->ylabel);
	box_append(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

	hbox=gtk_hbox_new(FALSE, 0);
	label=gtk_label_new("X min");gtk_label_set_width_chars(GTK_LABEL(label), 5);
	box_append(GTK_BOX(hbox), label, FALSE, FALSE, 0);
	box_append(GTK_BOX(hbox), spins[0], TRUE, TRUE, 0);

	label=gtk_label_new("Y min");gtk_label_set_width_chars(GTK_LABEL(label), 5);
	box_append(GTK_BOX(hbox), label, FALSE, FALSE, 0);
	box_append(GTK_BOX(hbox), spins[2], TRUE, TRUE, 0);

	box_append(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

	hbox=gtk_hbox_new(FALSE, 0);
	label=gtk_label_new("X max");gtk_label_set_width_chars(GTK_LABEL(label), 5);
	box_append(GTK_BOX(hbox), label, FALSE, FALSE, 0);
	box_append(GTK_BOX(hbox), spins[1], TRUE, TRUE, 0);
	
		
	label=gtk_label_new("Y max");gtk_label_set_width_chars(GTK_LABEL(label), 5);
	box_append(GTK_BOX(hbox), label, FALSE, FALSE, 0);
	box_append(GTK_BOX(hbox), spins[3], TRUE, TRUE, 0);
	
	box_append(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

	if(n>4){
		hbox=gtk_hbox_new(FALSE, 0);
		label=gtk_label_new("Z min");gtk_label_set_width_chars(GTK_LABEL(label), 5);
		box_append(GTK_BOX(hbox), label, FALSE, FALSE, 0);
		box_append(GTK_BOX(hbox), spins[4], TRUE, TRUE, 0);
		label=gtk_label_new("Z max");gtk_label_set_width_chars(GTK_LABEL(label), 5);
		box_append(GTK_BOX(hbox), label, FALSE, FALSE, 0);
		box_append(GTK_BOX(hbox), spins[5], TRUE, TRUE, 0);
		box_append(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
	}
	if(curwindow){//window size
		int width, height;
#if GTK_MAJOR_VERSION<4	
		gtk_window_get_size(GTK_WINDOW(curwindow), &width, &height);
#else
		width=gtk_widget_get_width(curwindow);
		height=gtk_widget_get_height(curwindow);
#endif	
		//dbg_time("width=%d, height=%d\n", width, height);
		hbox=gtk_hbox_new(FALSE, 0);
		label=new_spin_button("Width", width, 100, 2800, 10, set_window_size, GINT_TO_POINTER(1));
		box_append(GTK_BOX(hbox), label, TRUE,  TRUE, 0);
		label=new_spin_button("Height", height, 100, 8000, 10, set_window_size, GINT_TO_POINTER(2));
		box_append(GTK_BOX(hbox), label, TRUE,  TRUE, 0);

		checkbtn=gtk_check_button_new_with_label("Proportional");
		check_button_set_active(checkbtn, window_proportional);
		g_signal_connect(checkbtn, "toggled", G_CALLBACK(checkbtn_toggle), &window_proportional);
		box_append(GTK_BOX(hbox), checkbtn, FALSE, FALSE, 0);

		box_append(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
	}
#if GTK_VERSION_AFTER(4,10)
	gtk_window_set_child(GTK_WINDOW(dialog), vbox);
	gtk_window_present(GTK_WINDOW(dialog));
	//gtk_popover_set_child(GTK_POPOVER(popover),vbox);
	//gtk_popover_present(GTK_POPOVER(popover));
#else
	gtk_widget_show_all(vbox);
	gtk_container_add(GTK_CONTAINER(content_area), vbox);
	gtk_dialog_run(GTK_DIALOG(dialog));
	window_destroy(dialog);
#endif
	drawdata->spins=NULL;
}


void desc_updated(){
	font_name_version++;/*tell every expose event to update the figure. */
	if(font_name) free(font_name);
	font_name=strdup(pango_font_description_get_family(desc));

	/*PangoStyle style=pango_font_description_get_style(desc);
	switch(style){
	case PANGO_STYLE_NORMAL:
		font_style=CAIRO_FONT_SLANT_NORMAL;
		break;
	case PANGO_STYLE_OBLIQUE:
		font_style=CAIRO_FONT_SLANT_OBLIQUE;
		break;
	case PANGO_STYLE_ITALIC:
		font_style=CAIRO_FONT_SLANT_ITALIC;
		break;
	}*/
	/*
	PangoWeight weight=pango_font_description_get_weight(desc);
	if(weight<=500){
		font_weight=CAIRO_FONT_WEIGHT_NORMAL;
	} else{
		font_weight=CAIRO_FONT_WEIGHT_BOLD;
	}*/
	desc_font_size=pango_font_description_get_size(desc);
	/*get font_size in device unit (dots, pixels); */
	if(pango_font_description_get_size_is_absolute(desc)){
		font_size=(float)desc_font_size/(float)PANGO_SCALE;
		fprintf(stderr, "absolute: font_size=%g\n", font_size);
	} else{
#if GTK_MAJOR_VERSION<4
		gfloat dpi=gdk_screen_get_resolution(gdk_screen_get_default());
		if(dpi<0) dpi=72;
#else
		gfloat dpi=72;
#endif
		font_size=(float)desc_font_size/(float)PANGO_SCALE*(float)dpi/72.;
	}
	update_pixmap(get_current_drawdata());
}

static void tool_font_set(GtkWidget *btn){
	if(!desc) desc=pango_font_description_new();
	PangoFontDescription* desc2;
#if GTK_MAJOR_VERSION > 3
	//do not free desc. it belongs to the widget
	desc2=gtk_font_dialog_button_get_font_desc(GTK_FONT_DIALOG_BUTTON(btn));
#elif GTK_VERSION_AFTER(3,2)
	desc2=gtk_font_chooser_get_font_desc(GTK_FONT_CHOOSER(btn));
#else
	desc2=pango_font_description_from_string(gtk_font_button_get_font_name(GTK_FONT_BUTTON(btn)));
#endif	
	pango_font_description_set_family(desc, pango_font_description_get_family(desc2));
	pango_font_description_set_size(desc, pango_font_description_get_size(desc2));
#if !GTK_VERSION_AFTER(3,2)
	pango_font_description_free(desc2);
#endif	
	desc_updated();
}

#if GTK_MAJOR_VERSION>=3
static gboolean close_window(GtkWidget* window, GdkEvent* event)
#else
static gboolean close_window(GtkObject* window, GdkEvent* event)
#endif
{
	(void)event;
	windows=g_slist_remove(windows, window);
	info("close_window called\n");
	if(windows){
		window_changed(GTK_WIDGET(windows->data));
	}else{
		info("Close sock %d after last window is about to close.\n", sock);
		keep_listen=0;
		if(sock!=-1){
			close(sock);
			sock=-1;
		}
	}
	GtkWidget *topnb=get_topnb(GTK_WIDGET(window));
	GtkWidget *toolbar=get_toolbar(GTK_WIDGET(window));
	g_signal_handlers_disconnect_by_func(topnb, (gpointer)topnb_page_added, toolbar);
	g_signal_handlers_disconnect_by_func(topnb, (gpointer)topnb_page_removed, toolbar);
	g_signal_handlers_disconnect_by_func(topnb, (gpointer)topnb_page_switch, toolbar);
	if(windows){//there are existing windows. move tabs over
		GtkWidget* topnb2=get_topnb(GTK_WIDGET(windows->data));
		int npages=gtk_notebook_get_n_pages(GTK_NOTEBOOK(topnb));
		for(int ipage=npages-1; ipage>=0; ipage--){
			move_tab_page(topnb, ipage, topnb2);
		}
	}


	if(!windows){
#if GTK_MAJOR_VERSION < 4
		gtk_main_quit();
		info("call main_quit()\n");
#endif
		free(font_name); font_name=NULL;
		free(client_hostname);client_hostname=NULL;
		free(client_path_full); client_path_full=NULL;
		free(client_exename); client_exename=NULL;
	}
	/*don't call exit here. */
	return FALSE;
}
#if GTK_MAJOR_VERSION < 4
static gboolean window_state(GtkWidget* window, GdkEvent* event){
	if(event->type==GDK_FOCUS_CHANGE){
		GdkEventFocus* focus=(GdkEventFocus*)event;
		if(focus->in){
			window_changed(window);
		}
	}
	return FALSE;
}
#else
static gboolean window_activate_focus(GtkWidget *window){
	//info("window_activate_focus\n");
	window_changed(window);
	return FALSE;
}
#endif
gboolean update_fpslabel(gpointer data){
	gint repeat=GPOINTER_TO_INT(data);//TRUE: keep running, false: delete.
	char newtext[64];
	if(client_pid==0){
		snprintf(newtext, sizeof(newtext), "idle since %s", myasctime(0)+11);
		gtk_label_set_text(GTK_LABEL(fpslabel), newtext);
		repeat=FALSE;
	}else if(client_pid<0){
		snprintf(newtext, sizeof(newtext), "disconnected at %s", myasctime(0)+11);
		gtk_label_set_text(GTK_LABEL(fpslabel), newtext);
		repeat=FALSE;
	}else{
		float thistime=myclockd();
		extern float io_time1;//receiving time for latest frame
		extern float io_time2;//receiving time for previous frame or 0 if a different plot is received
		float fps=0;
		if(io_time1+10>thistime && io_time2+100>io_time1 && io_time1!=io_time2){//continuous update
			fps=1./(io_time1-io_time2);
			snprintf(newtext, sizeof(newtext), "%.1f Hz", fps);
		}else{
			newtext[0]=0;
			fps=0;
		}
		gtk_label_set_text(GTK_LABEL(fpslabel), newtext);
	}
	/*if(io_time1+60<thistime){//60 seconds no update, check connectivity.
		if(stwriteint(sock, DRAW_RESUME)){
			close(sock);
			sock=-1;
		}
	}*/
	return repeat;
}
//Create a new toolbar item.
/**
 * Create a new window with a topnb
*/
GtkWidget* create_window(GtkWidget* window){
	//GtkWidget* old_window=curwindow;
	if(!cursors[1]){
		cursors[0]=cursor_new_from_name("default");
		cursors[1]=cursor_new_from_name("crosshair");
#if GTK_MAJOR_VERSION < 3
		gtk_rc_parse_string(rc_string_notebook);
#endif
	}
#if GTK_MAJOR_VERSION < 4
	/*if(!contexmenu){
		contexmenu=gtk_menu_new();
		GtkWidget* item=gtk_menu_item_new_with_label("Detach");
		g_signal_connect(item, "activate", G_CALLBACK(topnb_detach), NULL);
		gtk_menu_shell_append(GTK_MENU_SHELL(contexmenu), item);
		gtk_widget_show_all(contexmenu);
		g_object_ref(contexmenu);
	}*/
#endif
	if(!window){
		window=window_new();
		/*if(old_window){
			gtk_window_set_transient_for(GTK_WINDOW(window), GTK_WINDOW(old_window));
		}*/
	}
	gtk_window_set_default_size(GTK_WINDOW(window), 720, 600);//1050, 900);
	gtk_window_set_resizable(GTK_WINDOW(window), TRUE);
	//gtk_window_stick(GTK_WINDOW(window));
	curwindow=window;
	windows=g_slist_append(windows, window);
	//g_signal_connect(window, "close-request", G_CALLBACK(window_close_request));
#if GTK_MAJOR_VERSION<4
	extern GdkPixbuf* icon_main;
	gtk_window_set_icon(GTK_WINDOW(window), icon_main);
#else
	gtk_window_set_icon_name(GTK_WINDOW(window), "computer");
#endif
	
#if GTK_MAJOR_VERSION>=3
	GtkCssProvider* provider_default=gtk_css_provider_new();
#if GTK_MAJOR_VERSION>=4
	gtk_css_provider_load_from_string(provider_default, all_style);
	gtk_style_context_add_provider_for_display(gdk_display_get_default(), GTK_STYLE_PROVIDER(provider_default),
		GTK_STYLE_PROVIDER_PRIORITY_USER);
	GtkSettings *settings=gtk_settings_get_for_display(gdk_display_get_default());
	g_object_set(G_OBJECT(settings), "gtk-hint-font-metrics", TRUE, NULL);//fix gtk4 subpixel positioning blurry.
#else
	gtk_css_provider_load_from_data(provider_default, all_style, strlen(all_style), NULL);

	GtkStyleContext* all_context=gtk_widget_get_style_context(window);
	GdkScreen* screen=gtk_style_context_get_screen(all_context);
	gtk_style_context_add_provider_for_screen(screen, GTK_STYLE_PROVIDER(provider_default),
		GTK_STYLE_PROVIDER_PRIORITY_USER);
#endif
#endif
#if GTK_MAJOR_VERSION<3
	if(!btn_rcstyle){
		btn_rcstyle=gtk_rc_style_new();
	//This option is used in many places to determine padding between text and border of widgets
		btn_rcstyle->xthickness=btn_rcstyle->ythickness=0;
	}
#endif
	GtkWidget *topnb=gtk_notebook_new();
#if USE_TITLEBAR==0
/*	GtkWidget* toolbar=gtk_toolbar_new();
	gtk_toolbar_set_icon_size(GTK_TOOLBAR(toolbar), GTK_ICON_SIZE_MENU);
	gtk_toolbar_set_style(GTK_TOOLBAR(toolbar), GTK_TOOLBAR_ICONS);*/
	GtkWidget *toolbar=gtk_hbox_new(0, 0);
#else
	GtkWidget* toolbar=gtk_header_bar_new();
#if GTK_MAJOR_VERSION<4
	gtk_header_bar_set_show_close_button(GTK_HEADER_BAR(toolbar), TRUE);
	gtk_header_bar_set_title(GTK_HEADER_BAR(toolbar), NULL);
#else
	gtk_header_bar_set_show_title_buttons(GTK_HEADER_BAR(toolbar), TRUE);
#endif
#endif
	//gtk_widget_set_sensitive(toolbar, FALSE);
	new_toolbar_item(toolbar, NULL, 0, "document-save-as", NULL, NULL, tool_save, NULL);
	new_toolbar_item(toolbar, NULL, 0, NULL, NULL, NULL, NULL, NULL); //separator
	new_toolbar_item(toolbar, NULL, 0, "zoom-in", NULL, NULL, tool_zoom, GINT_TO_POINTER(2));
	new_toolbar_item(toolbar, NULL, 0, "zoom-fit-best", NULL, NULL, tool_zoom, GINT_TO_POINTER(1));
	new_toolbar_item(toolbar, NULL, 0, "zoom-out", NULL, NULL, tool_zoom, GINT_TO_POINTER(0));
	new_toolbar_item(toolbar, NULL, 0, NULL, NULL, NULL, NULL, NULL); //separator
	new_toolbar_item(toolbar, NULL, 0, "document-properties", NULL, NULL, tool_property, NULL);

	//Log scale
	GtkWidget *menu_zlog=new_toolbar_item(toolbar, NULL, 1, NULL, icon_log, NULL, togglebutton_toggle, GINT_TO_POINTER(1));
	g_object_set_data(G_OBJECT(window), "menu_zlog", menu_zlog);
	//Cumulative
	GtkWidget *menu_cumu=new_toolbar_item(toolbar, NULL, 1, NULL, icon_avg, NULL, togglebutton_toggle, GINT_TO_POINTER(2));
	g_object_set_data(G_OBJECT(window), "menu_cumu", menu_cumu);

	//Cumulative step
	GtkWidget *menu_icumu=gtk_spin_button_new_with_range(0, 100000, 1);
	g_signal_connect(menu_icumu, "value-changed", G_CALLBACK(spin_icumu), NULL);
	new_toolbar_item(toolbar, menu_icumu, 0, "menu_icumu", NULL, NULL, NULL, NULL);
	g_object_set_data(G_OBJECT(window), "menu_icumu", menu_icumu);

	new_toolbar_item(toolbar, NULL, 0, NULL, NULL, NULL, NULL, NULL); //separator
#if GTK_VERSION_AFTER(4,10)
	GtkFontDialog *fontdiag=gtk_font_dialog_new();
	GtkWidget *fontsel=gtk_font_dialog_button_new(fontdiag);
	gtk_font_dialog_button_set_font_features(GTK_FONT_DIALOG_BUTTON(fontsel), font_name_default);
	//gtk_font_dialog_button_set_use_font(GTK_FONT_DIALOG_BUTTON(fontsel), TRUE);
	g_signal_connect(GTK_FONT_BUTTON(fontsel), "notify::font-desc", G_CALLBACK(tool_font_set), fontsel);
#else
	GtkWidget* fontsel=gtk_font_button_new_with_font(font_name_default);
	gtk_font_button_set_use_font(GTK_FONT_BUTTON(fontsel), TRUE);
	g_signal_connect(GTK_FONT_BUTTON(fontsel), "font-set", G_CALLBACK(tool_font_set), NULL);
#endif
	tool_font_set(fontsel);
	new_toolbar_item(toolbar, fontsel, 0, "font-set", NULL, NULL, NULL, NULL);
	new_toolbar_item(toolbar, NULL, 0, "view-refresh-symbolic", NULL, NULL, toolbutton_refresh, NULL);
	new_toolbar_item(toolbar, NULL, 1, "media-playback-pause", NULL, NULL, togglebutton_pause, NULL);
	new_toolbar_item(toolbar, NULL, 1, "media-playback-stop", NULL, NULL, togglebutton_stop, NULL);
	new_toolbar_item(toolbar, NULL, 0, "edit-copy", NULL, NULL, topnb_detach_btn, topnb);
	new_toolbar_item(toolbar, NULL, 1, NULL, NULL, NULL, NULL, NULL); //separator

	if(!fpslabel){//create only for the main window.
		fpslabel=gtk_label_new("");
		/*
	#if GTK_MAJOR_VERSION>=3	
		gtk_widget_set_hexpand(fpslabel, TRUE);
		gtk_widget_set_halign(fpslabel, GTK_ALIGN_END);
	#else
		gtk_misc_set_padding(GTK_MISC(fpslabel), 0, 0);
		gtk_misc_set_alignment(GTK_MISC(fpslabel), 1, 0.5);
	#endif*/
		new_toolbar_item(toolbar, fpslabel, 0, "fps", NULL, NULL, NULL, NULL);
	}
	//gtk_container_set_border_width(GTK_CONTAINER(topnb),2);
#if GTK_VERSION_AFTER(2,24)
	gtk_notebook_set_group_name(GTK_NOTEBOOK(topnb), "toplevel");//2.24
#elif GTK_VERSION_AFTER(2,12)
	gtk_notebook_set_group(GTK_NOTEBOOK(topnb), "toplevel");//2.12
#endif
	gtk_notebook_set_scrollable(GTK_NOTEBOOK(topnb), TRUE);
	//gtk_notebook_set_show_border(GTK_NOTEBOOK(topnb), FALSE);
	g_signal_connect(GTK_NOTEBOOK(topnb), "page-added", G_CALLBACK(topnb_page_added), toolbar);
	g_signal_connect(GTK_NOTEBOOK(topnb), "page-removed", G_CALLBACK(topnb_page_removed), toolbar);
	g_signal_connect_after(GTK_NOTEBOOK(topnb), "switch-page", G_CALLBACK(topnb_page_switch), toolbar);
#if USE_TITLEBAR==0
	GtkWidget *vbox=gtk_vbox_new(FALSE, 0);
	gtk_box_set_homogeneous(GTK_BOX(vbox), FALSE);
	box_append(GTK_BOX(vbox), toolbar, FALSE, FALSE, 0);
	box_append(GTK_BOX(vbox), topnb, TRUE, TRUE, 0);
	gtk_container_add(GTK_CONTAINER(window), vbox);
#else
	gtk_widget_set_vexpand(topnb, TRUE);
	gtk_window_set_titlebar(GTK_WINDOW(window), toolbar);
#if GTK_MAJOR_VERSION<4
	gtk_container_add(GTK_CONTAINER(window), topnb);
#else
	gtk_window_set_child(GTK_WINDOW(window), topnb);
#endif
#endif	
	
#if GTK_MAJOR_VERSION<4
	gtk_window_set_position(GTK_WINDOW(window), GTK_WIN_POS_CENTER);
	gtk_widget_add_events(GTK_WIDGET(window), GDK_FOCUS_CHANGE_MASK);/*in case this is not on. */
	g_signal_connect(GTK_WIDGET(window), "event", G_CALLBACK(window_state), NULL);
	gtk_widget_show_all(window);
	//tool_font_set(GTK_FONT_BUTTON(fontsel));/*initialize. */
#else
	g_signal_connect(GTK_WINDOW(window), "activate-focus", G_CALLBACK(window_activate_focus), NULL);

	//dnd
	//GtkDragSource *drag_src=gtk_drag_source_new();
	//g_signal_connect(drag_src, "prepare", G_CALLBACK(topnb_drag_prepare), topnb);
	//g_signal_connect(drag_src, "drag-begin", G_CALLBACK(topnb_drag_begin), topnb);
	//gtk_widget_add_controller(GTK_WIDGET(topnb), GTK_EVENT_CONTROLLER(drag_src));

#endif
	//font_set(font_name_default);

#if GTK_MAJOR_VERSION<4
	g_signal_connect(window, "delete-event", G_CALLBACK(close_window), NULL);
#else
	g_signal_connect(window, "close-request", G_CALLBACK(close_window), NULL);
#endif

	window_changed(window);
	gtk_window_present(GTK_WINDOW(window));

	return window;
}

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
#include <tgmath.h>
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
#define DRAWAREA_MIN_WIDTH 440
#define DRAWAREA_MIN_HEIGHT 320
#define MAX_ZOOM 10000
#define MIN_ZOOM 0.01
static GSList* windows=NULL;
static int iwindow=0;
static GtkWidget* curwindow=NULL;
static GtkWidget* curtopnb=NULL;
static GtkWidget* contexmenu=NULL;
static drawdata_t* drawdata_dialog=NULL;
//static GtkToolItem* toggle_cumu=NULL;
#if GTK_MAJOR_VERSION<3
static GtkRcStyle* btn_rcstyle=NULL;
#endif
PangoFontDescription* desc=NULL;
int font_name_version=0;
char* font_name=NULL;
float font_size;
cairo_font_slant_t font_style=CAIRO_FONT_SLANT_NORMAL;
cairo_font_weight_t font_weight=CAIRO_FONT_WEIGHT_NORMAL;
float lpf=1;//low pass filter the update. Set using dialog.
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

#define error_msg(A...) {				\
	GtkWidget *dialog0=gtk_message_dialog_new	\
	    (GTK_WINDOW(curwindow),		\
	     GTK_DIALOG_DESTROY_WITH_PARENT,		\
	     GTK_MESSAGE_ERROR,				\
	     GTK_BUTTONS_CLOSE,				\
	     A);					\
	gtk_dialog_run(GTK_DIALOG(dialog0));		\
	gtk_widget_destroy(dialog0);			\
    }

static void set_cur_window(GtkWidget* window){
	curwindow=window;
#if GTK_MAJOR_VERSION>=4 		
	GtkWidget* vbox=gtk_window_get_child(GTK_WINDOW(window));

	curtopnb=gtk_widget_get_next_sibling(gtk_widget_get_first_child(vbox));
#else
	GtkWidget* vbox=gtk_bin_get_child(GTK_BIN(curwindow));
	GList* list=gtk_container_get_children(GTK_CONTAINER(vbox));
	curtopnb=(GtkWidget*)list->next->data;
	g_list_free(list);
#endif
}
/*Get the current page for notebook*/
static GtkWidget* get_current_page(GtkWidget* notebook){
	int n=gtk_notebook_get_current_page(GTK_NOTEBOOK(notebook));
	/*if(n<0){
	warning_once("get_current_page returns %d\n", n);
	}*/
	return gtk_notebook_get_nth_page(GTK_NOTEBOOK(notebook), n);
}
static drawdata_t* get_current_drawdata(void){
	GtkWidget* topnb=curtopnb;
	GtkWidget* w1=get_current_page(topnb);
	GtkWidget* w2=get_current_page(w1);
	drawdata_t** pdrawdata=(drawdata_t**)g_object_get_data(G_OBJECT(w2), "drawdatawrap");
	if(pdrawdata){
		return *pdrawdata;
	} else{
		return NULL;
	}
}
static GtkWidget* get_topnb(GtkWidget* window){
#if GTK_MAJOR_VERSION>=4 
	GtkWidget* vbox=gtk_window_get_child(GTK_WINDOW(window));
	GtkWidget* toolbar=gtk_widget_get_first_child(vbox);
	GtkWidget* topnb=gtk_widget_get_next_sibling(toolbar);
#else
	GtkWidget* vbox=gtk_bin_get_child(GTK_BIN(window));
	GList* list=gtk_container_get_children(GTK_CONTAINER(vbox));
	GtkWidget* topnb=(GtkWidget*)list->next->data;
	g_list_free(list);
#endif
	return topnb;
}
static GtkWidget* get_toolbar(GtkWidget* window){
#if GTK_MAJOR_VERSION>=4 
	GtkWidget* vbox=gtk_window_get_child(GTK_WINDOW(window));
	GtkWidget* toolbar=gtk_widget_get_first_child(vbox);
#else
	GtkWidget* vbox=gtk_bin_get_child(GTK_BIN(window));
	GList* list=gtk_container_get_children(GTK_CONTAINER(vbox));
	GtkWidget* toolbar=(GtkWidget*)list->data;
	g_list_free(list);
#endif
	return toolbar;
}

static const char* topnb_label_get(GtkWidget* topnb, GtkWidget* subnb){
	GtkWidget* label=gtk_notebook_get_tab_label(GTK_NOTEBOOK(topnb), subnb);
#if GTK_MAJOR_VERSION>=4
	GtkWidget* label2=gtk_widget_get_first_child(label);
#else
	GtkWidget* label2=gtk_bin_get_child(GTK_BIN(label));
#endif
	return gtk_label_get_text(GTK_LABEL(label2));
}
#if GTK_MAJOR_VERSION<4
static void topnb_detach(GtkMenuItem* menu){
	(void)menu;
	GtkWidget* topnb=curtopnb;
	int n=gtk_notebook_get_current_page(GTK_NOTEBOOK(topnb));
	GtkWidget* page=gtk_notebook_get_nth_page(GTK_NOTEBOOK(topnb), n);
	GtkWidget* label=gtk_notebook_get_tab_label(GTK_NOTEBOOK(topnb), page);
	g_object_ref(page);
	g_object_ref(label);
	gtk_notebook_remove_page(GTK_NOTEBOOK(topnb), n);
	GtkWidget* window=create_window(NULL);/*create a new window. */
	topnb=get_topnb(window);/*another topnb */
	gtk_notebook_append_page(GTK_NOTEBOOK(topnb), page, label);
	gtk_notebook_set_tab_detachable(GTK_NOTEBOOK(topnb), page, TRUE);
	gtk_notebook_set_tab_reorderable(GTK_NOTEBOOK(topnb), page, TRUE);
	g_object_unref(page);
	g_object_unref(label);
}
static gboolean tab_button_cb(GtkWidget* widget, GdkEventButton* event, GtkWidget* page){
	/*
	  widget is the event box that is the page label. page is the page in the notebook.
	*/
	GtkWidget* topnb=gtk_widget_get_parent(page);
	int n=gtk_notebook_page_num(GTK_NOTEBOOK(topnb), page);
	gtk_notebook_set_current_page(GTK_NOTEBOOK(topnb), n);
	if(event->button==3){
		if(gtk_menu_get_attach_widget(GTK_MENU(contexmenu))){
			gtk_menu_detach(GTK_MENU(contexmenu));
		}
		gtk_menu_attach_to_widget(GTK_MENU(contexmenu), widget, NULL);
#if GTK_MAJOR_VERSION>=3 		
		gtk_menu_popup_at_pointer(GTK_MENU(contexmenu), NULL);
#else
		gtk_menu_popup(GTK_MENU(contexmenu), NULL, NULL, NULL, NULL, 3, gtk_get_current_event_time());
#endif
	}
	return FALSE;
}
#endif
typedef struct updatetimer_t{
	int pending;
	drawdata_t* drawdata;
}updatetimer_t;
static void update_pixmap(drawdata_t* drawdata){
	/*no more pending updates, do the updating. */
	gint width=drawdata->width;
	gint height=drawdata->height;
	if(drawdata->pixmap){
		if(width!=drawdata->pwidth||height!=drawdata->pheight){
#if GTK_MAJOR_VERSION>=3 
			cairo_surface_destroy(drawdata->pixmap);
#else
			g_object_unref(drawdata->pixmap);
#endif
			drawdata->pixmap=NULL;
		}
	}
	if(!drawdata->pixmap){
	/*Create a new server size pixmap and then draw on it. */
		drawdata->pwidth=width;
		drawdata->pheight=height;
#if GTK_MAJOR_VERSION>=4
		cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
#elif GTK_MAJOR_VERSION>=3 
		drawdata->pixmap=gdk_window_create_similar_surface
		(gtk_widget_get_window(curwindow), CAIRO_CONTENT_COLOR_ALPHA, width, height);
#else
		drawdata->pixmap=gdk_pixmap_new(curwindow->window, width, height, -1);
#endif
	}
	if(drawdata->p0||drawdata->npts){
		/*cairo_t is destroyed in draw */
#if GTK_MAJOR_VERSION>=3 
		cairo_draw(cairo_create(drawdata->pixmap), drawdata, width, height);
#else
		cairo_draw(gdk_cairo_create(drawdata->pixmap), drawdata, width, height);
#endif
	}
	gtk_widget_queue_draw(drawdata->drawarea);
}
static gboolean update_pixmap_timer(updatetimer_t* timer){
	drawdata_t* drawdata=timer->drawdata;
	if(timer->pending==drawdata->pending){
		update_pixmap(drawdata);
	}
	free(timer);
	return FALSE;
}
/**
   Call the update_pixmap after time out. If another call is queue within
   the time, the "pending" counter will mismatch and no action will be taken.
 */
static void delayed_update_pixmap(drawdata_t* drawdata){
	if(!drawdata->pixmap){
		update_pixmap(drawdata);
	} else{
		drawdata->pending++;
		updatetimer_t* tmp=mycalloc(1, updatetimer_t);
		tmp->pending=drawdata->pending;
		tmp->drawdata=drawdata;
		g_timeout_add(20, (GSourceFunc)update_pixmap_timer, tmp);
	}
}

#if  GTK_MAJOR_VERSION>=4
static gboolean
on_resize_event(GtkDrawingArea* widget, int width, int height, gpointer pdata){
	(void)widget;
	drawdata_t* drawdata=*((drawdata_t**)pdata);
	if(width>1&&height>1){
		drawdata->width=width;
		drawdata->height=height;
		delayed_update_pixmap(drawdata);
	}
	return FALSE;
}
#else
/*
   Expose event happens when the widget is first configured and whenever it is
   resized.  Create a new backup pixmap of the appropriate size and draw there.*/
static gboolean
on_configure_event(GtkWidget* widget, GdkEventConfigure* event, gpointer pdata){
	(void)event;
	(void)widget;
	drawdata_t* drawdata=*((drawdata_t**)pdata);
	if(event->width>1&&event->height>1){
		drawdata->width=event->width;
		drawdata->height=event->height;
		delayed_update_pixmap(drawdata);
	}
	return FALSE;
}
#endif
/* Redraw the screen from the backing pixmap */
#if GTK_MAJOR_VERSION>=3
static gboolean
on_draw_event(GtkWidget* widget, cairo_t* cr, gpointer pdata){
	(void)widget;
	drawdata_t* drawdata=*((drawdata_t**)pdata);
	if(drawdata->font_name_version!=font_name_version||!drawdata->drawn){
		update_pixmap(drawdata);
	}
	cairo_set_source_surface(cr, drawdata->pixmap, 0, 0);
	cairo_paint(cr);
	if(drawdata->draw_rect){
		cairo_set_source_rgba(cr, 0, 0, 1, 0.1);
		cairo_set_line_width(cr, 1);
		cairo_rectangle(cr, drawdata->mxdown, drawdata->mydown, drawdata->dxdown, drawdata->dydown);
		cairo_fill_preserve(cr);
		cairo_set_source_rgba(cr, 0, 0, 0, 1);
		cairo_stroke(cr);
		drawdata->draw_rect=0;
	}
	return FALSE;
}
#else
static gboolean
on_expose_event(GtkWidget* widget, GdkEventExpose* event, gpointer pdata){
	drawdata_t* drawdata=*((drawdata_t**)pdata);
	int cumu_effective=cumu&&!drawdata->image&&!drawdata->square;
	if(drawdata->font_name_version!=font_name_version||!drawdata->drawn||drawdata->cumu!=cumu_effective){
		delayed_update_pixmap(drawdata);
	}
	if(drawdata->pixmap){
		gdk_draw_drawable(widget->window,
			widget->style->fg_gc[GTK_WIDGET_STATE(widget)],
			drawdata->pixmap,
			event->area.x, event->area.y,
			event->area.x, event->area.y,
			event->area.width, event->area.height);
		if(drawdata->draw_rect){
			cairo_t* cr=gdk_cairo_create(widget->window);
			cairo_set_source_rgba(cr, 0, 0, 1, 0.1);
			cairo_set_line_width(cr, 1);
			cairo_rectangle(cr, drawdata->mxdown, drawdata->mydown, drawdata->dxdown, drawdata->dydown);
			cairo_fill_preserve(cr);
			cairo_set_source_rgba(cr, 0, 0, 0, 1);
			cairo_stroke(cr);
			cairo_destroy(cr);
			drawdata->draw_rect=0;
		}
	} else{
		warning("pixmap is empty\n");
	}
	return FALSE;
}
#endif
#define FREE(A) free(A); A=NULL;
static void drawdata_free_input(drawdata_t* drawdata){
	/*Only free the input received via fifo from draw.c */
	if(drawdata->image){
		cairo_surface_destroy(drawdata->image);
		drawdata->image=NULL;
	}
	FREE(drawdata->p);
	if(drawdata->npts>0){
		for(int ipts=0; ipts<drawdata->npts; ipts++){
			FREE(drawdata->pts[ipts]);
		}
		FREE(drawdata->pts);
		FREE(drawdata->ptsdim);
	}
	if(drawdata->nstyle>0){
		FREE(drawdata->style);
	}
	if(drawdata->ncir>0){
		FREE(drawdata->cir);
	}
	FREE(drawdata->fig);
	FREE(drawdata->name);
	FREE(drawdata->title);
	FREE(drawdata->xlabel);
	FREE(drawdata->ylabel);
	if(drawdata->legend){
		for(int i=0; i<drawdata->npts; i++){
			FREE(drawdata->legend[i]);
		}
		FREE(drawdata->legend);
	}
}
static void drawdata_free(drawdata_t* drawdata){
	if(drawdata->pixmap){
#if GTK_MAJOR_VERSION>=3 
		cairo_surface_destroy(drawdata->pixmap);
#else
		g_object_unref(drawdata->pixmap);
#endif
		drawdata->pixmap=NULL;
	}
	if(drawdata->cacheplot){
		cairo_surface_destroy(drawdata->cacheplot);
		drawdata->cacheplot=NULL;
	}
	drawdata_free_input(drawdata);
	FREE(drawdata->p0);
	FREE(drawdata);
	LOCK(drawdata_mutex);
	ndrawdata--;
	UNLOCK(drawdata_mutex);
}

/**
   Delete a figure page.
*/
static void delete_page(GtkButton* btn, drawdata_t** drawdatawrap){
	(void)btn;
	info("deleting page\n");
	GtkWidget* root;
	/*First find the root page */
	GtkWidget* topnb=curtopnb;
	root=get_current_page(topnb);
	/*Find the sub page. it may not be the current page */
	int ipage=gtk_notebook_page_num(GTK_NOTEBOOK(root), (*drawdatawrap)->page);
	gtk_notebook_remove_page(GTK_NOTEBOOK(root), ipage);
	drawdata_free(*drawdatawrap);
	free(drawdatawrap);
	if(gtk_notebook_get_n_pages(GTK_NOTEBOOK(root))==0){
		info("delete top page");
		int jpage=gtk_notebook_page_num(GTK_NOTEBOOK(topnb), root);
		gtk_notebook_remove_page(GTK_NOTEBOOK(topnb), jpage);
	}
}

#if ! defined(__linux__)
/**
   Delete a figure page using button-press-event. button works not good in mac.
*/
#if GTK_MAJOR_VERSION>=4
static void delete_page_event(GtkWidget* widget, drawdata_t** drawdatawrap){
	(void)widget;
	delete_page(NULL, drawdatawrap);
}
#else
static void delete_page_event(GtkWidget* widget, GdkEventButton* event, drawdata_t** drawdatawrap){
	(void)widget;
	(void)event;
	if(event->button==1&&event->type==GDK_BUTTON_PRESS){
		delete_page(NULL, drawdatawrap);
	}
}
#endif
#endif
static GtkWidget* subnb_label_new(drawdata_t** drawdatawrap){
	GtkWidget* out;
#if GTK_MAJOR_VERSION>=4
	out=gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
	//GtkWidget* image=gtk_image_new_from_icon_name("window-close");
#else
	out=gtk_hbox_new(FALSE, 0);
	GtkWidget* image=gtk_image_new_from_icon_name("window-close", GTK_ICON_SIZE_MENU);
#endif

#if defined(__linux__)
	GtkWidget* close_btn;
	close_btn=gtk_button_new();
#if GTK_MAJOR_VERSION>=3
	gtk_widget_set_focus_on_click(close_btn, FALSE);
#else
	gtk_button_set_focus_on_click(GTK_BUTTON(close_btn), FALSE);
#endif
#if GTK_MAJOR_VERSION>=4
	gtk_button_set_icon_name(GTK_BUTTON(close_btn), "window-close");
#else
	gtk_container_add(GTK_CONTAINER(close_btn), image);
#endif
	g_signal_connect(close_btn, "clicked", G_CALLBACK(delete_page), drawdatawrap);
	/* make button as small as possible */
#if GTK_MAJOR_VERSION<3
	gtk_button_set_relief(GTK_BUTTON(close_btn), GTK_RELIEF_NONE);
	gtk_widget_modify_style(close_btn, btn_rcstyle);
	gtk_widget_set_size_request(close_btn, 20, 20);
	gtk_container_set_border_width(GTK_CONTAINER(close_btn), 0);
#endif

	gtk_box_pack_end(GTK_BOX(out), close_btn, FALSE, FALSE, 0);
#else//macos


#if GTK_MAJOR_VERSION>=4
	GtkWidget* ebox=gtk_button_new();
	gtk_button_set_icon_name(GTK_BUTTON(ebox), "window-close");
	g_signal_connect(ebox, "clicked", G_CALLBACK(delete_page_event), drawdatawrap);
	gtk_box_append(GTK_BOX(out), ebox);
#else
	GtkWidget* ebox=gtk_event_box_new();
	gtk_container_add(GTK_CONTAINER(ebox), image);
	g_signal_connect(ebox, "button_press_event", G_CALLBACK(delete_page_event), drawdatawrap);
	gtk_box_pack_end(GTK_BOX(out), ebox, FALSE, FALSE, 0);
#endif

#endif

	/* create label for tab */
	drawdata_t* drawdata=*drawdatawrap;
	const gchar* str=drawdata->name;
	GtkWidget* label=gtk_label_new(str);
	//
#if GTK_MAJOR_VERSION >= 3
	gtk_widget_set_halign(label, GTK_ALIGN_START);
	gtk_widget_set_margin_start(label, 0);
	gtk_widget_set_margin_end(label, 0);
	gtk_widget_set_margin_top(label, 0);
	gtk_widget_set_margin_bottom(label, 0);
#else	
	gtk_misc_set_padding(GTK_MISC(label), 0, 0);
	gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
#endif
	gtk_label_set_width_chars(GTK_LABEL(label), 12);
#if GTK_MAJOR_VERSION>=4
	gtk_box_append(GTK_BOX(out), label);
	gtk_widget_show(out);
#else
	gtk_box_pack_start(GTK_BOX(out), label, TRUE, TRUE, 0);
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
	drawdata->offx+=xdiff/drawdata->zoomx;
	drawdata->offy+=ydiff/drawdata->zoomy;
	delayed_update_pixmap(drawdata);/*no need delay since motion notify already did it. */
}

static gboolean motion_notify(GtkWidget* widget, GdkEventMotion* event,
	drawdata_t** drawdatawrap){
	(void)widget;
	drawdata_t* drawdata=*drawdatawrap;
	if(((event->state&GDK_BUTTON1_MASK)||(event->state&GDK_BUTTON3_MASK))
		&&drawdata->valid){
		float x, y;
		x=event->x;
		y=event->y;
		float dx, dy;
		dx=x-drawdata->mxdown;
		dy=y-drawdata->mydown;
		float dt=myclockd()-drawdata->mtdown;
		/*move with left cursor */
		if((fabs(dx)>3 || fabs(dy)>3)&&dt>.16){
			if((event->state&GDK_BUTTON1_MASK)){
				do_move(drawdata, dx, -dy);/*notice the reverse sign. */
				drawdata->mxdown=x;
				drawdata->mydown=y;
			} else if(event->state&GDK_BUTTON3_MASK){/*select and zoom. */
				if(drawdata->square&&!drawdata->image){/*for a square */
					float ratio=1;
					if(drawdata->image){
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
					gtk_widget_queue_draw(widget);
				}
			}
		}
	}
	/*we set the cursor */
	if(event->x>drawdata->xoff&&event->x < drawdata->xoff+drawdata->widthim
		&&event->y > drawdata->yoff&&event->y<drawdata->yoff+drawdata->heightim){
		if(!drawdata->cursorinside){
			drawdata->cursorinside=1;
			gdk_window_set_cursor(gtk_widget_get_window(widget), cursors[0]);
			gtk_widget_set_has_tooltip(widget, 1);
		}
		{
			float x=(event->x-drawdata->xoff)/(float)(drawdata->widthim);
			float y=(event->y-drawdata->yoff)/(float)(drawdata->heightim);
			x=(1.-x)*drawdata->limit0[0]+x*drawdata->limit0[1];
			y=(1.-y)*drawdata->limit0[3]+y*drawdata->limit0[2];
			if(drawdata->xylog[0]!='n') x=pow(10, x);
			if(drawdata->xylog[1]!='n') y=pow(10, y);
			snprintf(drawdata->tooltip, sizeof(drawdata->tooltip), "(%g, %g)", x, y);
			gtk_widget_set_tooltip_text(widget, drawdata->tooltip);
		}
	} else{
		if(drawdata->cursorinside){
			drawdata->cursorinside=0;
			gdk_window_set_cursor(gtk_widget_get_window(widget), NULL);
			gtk_widget_set_has_tooltip(widget, 0);
		}
	}
	return FALSE;
}


static void do_zoom(drawdata_t* drawdata, float xdiff, float ydiff, int mode){
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
	}
	if(drawdata->zoomx<MIN_ZOOM)
		drawdata->zoomx=MIN_ZOOM;
	else if(drawdata->zoomx>MAX_ZOOM){
		drawdata->zoomx=MAX_ZOOM;
	}
	if(drawdata->zoomy<MIN_ZOOM)
		drawdata->zoomy=MIN_ZOOM;
	else if(drawdata->zoomy>MAX_ZOOM){
		drawdata->zoomy=MAX_ZOOM;
	}
	if(mode){/*not zero. */
		float factorx=1/old_zoomx-1/drawdata->zoomx;
		drawdata->offx-=xdiff*factorx;
		float factory=1/old_zoomy-1/drawdata->zoomy;
		drawdata->offy-=ydiff*factory;
	}
	delayed_update_pixmap(drawdata);
}

static gboolean scroll_event(GtkWidget* widget, GdkEventScroll* event,
	drawdata_t** drawdatawrap){
	(void)widget;
	drawdata_t* drawdata=*drawdatawrap;
	float xdiff=event->x-drawdata->centerx;
	float ydiff=-(event->y-drawdata->centery);
	if(event->direction==GDK_SCROLL_UP){
		do_zoom(drawdata, xdiff, ydiff, 1);
	} else if(event->direction==GDK_SCROLL_DOWN){
		do_zoom(drawdata, xdiff, ydiff, -1);
	}
	return FALSE;
}
static gboolean focus_in_handler(GtkWidget* widget, GdkEvent* event, drawdata_t** drawdatawrap){
	(void)event;
	(void)widget;
	drawdata_t* drawdata=*drawdatawrap;
	drawdata->valid=0;
	//dbg_time("focus_in_handler.\n");
	return FALSE;
}
static gboolean button_press(GtkWidget* widget, GdkEventButton* event, drawdata_t** drawdatawrap){
	drawdata_t* drawdata=*drawdatawrap;
	/*Grab focus so the keys work */
	if(!GTK_WIDGET_HAS_FOCUS(widget))
		gtk_widget_grab_focus(widget);

	if(event->x>drawdata->xoff&&event->x < drawdata->xoff+drawdata->widthim
		&&event->y > drawdata->yoff&&event->y<drawdata->yoff+drawdata->heightim){
		drawdata->mxdown=event->x;
		drawdata->mydown=event->y;
		drawdata->mtdown=myclockd();
		drawdata->valid=1;
	} else{
		drawdata->valid=0;
	}
	//dbg_time("button_press %g %g.\n", event->x, event->y);
	return FALSE;
}


static gboolean button_release(GtkWidget* widget, GdkEventButton* event, drawdata_t** drawdatawrap){
	(void)widget;
	drawdata_t* drawdata=*drawdatawrap;
	if(!drawdata->valid) return FALSE;
	float x, y;
	x=event->x;
	y=event->y;
	float dx=x-drawdata->mxdown;
	float dy=y-drawdata->mydown;
	float dt=myclockd()-drawdata->mtdown;
	//dbg_time("button_release %g %g. dx=%g %g. button is %d.dt is %g\n", event->x, event->y, dx, dy, event->button, dt);
	if((fabs(dx)<3&&fabs(dy)<3)||dt<0.16){
		//dbg_time("Ignore accidental click\n");
	} else if(event->button==1){/*move only on left button */
		do_move(drawdata, dx, -dy);
	} else if(event->button==3){/*right button select and zoom. */
		float xx=drawdata->mxdown;
		float yy=drawdata->mydown;
		if(drawdata->square&&!drawdata->image){
			float ratio=1;
			if(drawdata->image){
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
		drawdata->limit_changed=1;
		apply_limit(drawdata);
		update_pixmap(drawdata);
	}

	return FALSE;

}


static void switch_tab(int lr, int ud){
	GtkWidget* topnb=curtopnb;
	if(lr){
		gtk_notebook_set_current_page
		(GTK_NOTEBOOK(topnb),
			gtk_notebook_get_current_page(GTK_NOTEBOOK(topnb))+lr);
	} else if(ud){
		GtkWidget* page=gtk_notebook_get_nth_page
		(GTK_NOTEBOOK(topnb),
			gtk_notebook_get_current_page(GTK_NOTEBOOK(topnb)));
		gtk_notebook_set_current_page
		(GTK_NOTEBOOK(page),
			gtk_notebook_get_current_page(GTK_NOTEBOOK(page))+ud);
	}
}
static gboolean key_press(GtkWidget* widget, GdkEventKey* event,
	drawdata_t** drawdatawrap){
	(void)widget;

	if(event->state&GDK_CONTROL_MASK){
		switch(event->keyval){
		case GDK_Left:
			switch_tab(-1, 0);break;
		case GDK_Right:
			switch_tab(1, 0);break;
		case GDK_Up:
			switch_tab(0, -1);break;
		case GDK_Down:
			switch_tab(0, 1);break;
		}
	} else{
		switch(event->keyval){
		case GDK_plus:
		case GDK_equal:
			do_zoom(*drawdatawrap, 0, 0, 1); break;
		case GDK_minus:
			do_zoom(*drawdatawrap, 0, 0, -1);break;
		case GDK_0:
		case GDK_1:
			do_zoom(*drawdatawrap, 0, 0, 0);break;
		case GDK_Left:
			do_move(*drawdatawrap, -10, 0);break;
		case GDK_Right:
			do_move(*drawdatawrap, 10, 0);break;
		case GDK_Up:
			do_move(*drawdatawrap, 0, 10);break;
		case GDK_Down:
			do_move(*drawdatawrap, 0, -10);break;
		}
	}
	return TRUE;
}
static void page_changed(int topn, int subn){
	//dbg("page changed. sock=%d. topn=%d, subn=%d\n", sock, topn, subn);
	if(sock==-1||sock_idle) return;
	GtkWidget* topnb=curtopnb;
	GtkWidget* subnb, * subpage;
	if(topn>-1){
		subnb=gtk_notebook_get_nth_page(GTK_NOTEBOOK(topnb), topn);
	} else{
		subnb=get_current_page(topnb);
	}
	if(!subnb) return;
	const char* fig=topnb_label_get(topnb, subnb);
	if(subn>-1){
		subpage=gtk_notebook_get_nth_page(GTK_NOTEBOOK(subnb), subn);
	} else{
		subpage=get_current_page(subnb);
	}
	if(!subpage) return;
	const char* fn=subnb_label_get(subnb, subpage);
	//info("send fig=%s, fn=%s", fig, fn);
	if(stwriteint(sock, DRAW_FIGFN)||
		stwritestr(sock, fig)||
		stwritestr(sock, fn)){
		warning("Talk to client failed\n");
		sock_idle=1;
	}
	//info("done\n");
}
/*These signal handlers are called before the notebook page switch is done.*/
static void topnb_page_switch(GtkNotebook* topnb, GtkWidget* page, guint n, GtkWidget* toolbar){
	(void)topnb; (void)page;
	gtk_widget_set_sensitive(toolbar, TRUE);
	page_changed(n, -1);
}
static void subnb_page_switch(GtkNotebook* subnb, GtkWidget* page, guint n, gpointer dummy){
	(void)subnb; (void)page;
	(void)dummy;
	page_changed(-1, n);
}
static void topnb_page_changed(GtkNotebook* topnb, GtkWidget* child, guint n, GtkWidget* toolbar){
	(void)child;
	(void)n;
	int npage=gtk_notebook_get_n_pages(topnb);
	if(npage==0){/*no more pages left. */
		if(g_slist_length(windows)>1){
			GtkWidget* window=gtk_widget_get_parent
			(gtk_widget_get_parent(GTK_WIDGET(topnb)));
#if GTK_MAJOR_VERSION>=3
			gtk_widget_destroy(GTK_WIDGET(window));
#else
			gtk_object_destroy(GTK_OBJECT(window));
#endif
		} else{//this is the last window
			gtk_widget_set_sensitive(toolbar, FALSE);
		}
	} else{
		gtk_widget_set_sensitive(toolbar, TRUE);
		page_changed(n, -1);
	}
	/*gtk_notebook_set_show_tabs(topnb, npage!=1); */
}
#define DO_LPF(T,pold,pnew,n) \
T*p0old=pold; T*p0new=pnew;\
for(long i=0; i<n; i++) p0new[i]=p0new[i]*(lpf)+p0old[i]*(1.-lpf);

/**
   2012-10-27: GTK 3.6 deprecated gdk_threads_enter(). So it is hard to call
   addpage from child threads. Modify the routine so that addpage is called by
   the main thread when gtk is idle.
*/
gboolean addpage(gpointer indata){
	drawdata_t** drawdatawrap=(drawdata_t**)indata;
	drawdata_t* drawdata=*drawdatawrap;
	GtkWidget* drawarea;
	if(!drawdata->fig){
		warning("Must set fig before calling addpage");
	}
	GSList* subnbs=NULL;
	int nsubnb=0;
	GtkWidget* window=0;
	GtkWidget* topnb=0;
	int itab=0;
	for(GSList* p=windows; p&&!nsubnb; p=p->next){
	//scan through all window to find all topnb page that has the same "fig"
		window=(GtkWidget*)p->data;
		topnb=(GtkWidget*)get_topnb(window);
		for(itab=0; itab<gtk_notebook_get_n_pages(GTK_NOTEBOOK(topnb)); itab++){
			GtkWidget* subnb=gtk_notebook_get_nth_page(GTK_NOTEBOOK(topnb), itab);
			const gchar* text=topnb_label_get(topnb, subnb);
			int res=strcmp(text, drawdata->fig);
			if(!res){//found
				subnbs=g_slist_append(subnbs, subnb);
				nsubnb++;/*number of subnbs find. */
				break;
			} else if(res>0){//not found. Insert here.
				break;
			}
		}
	}
	if(!nsubnb){/*subnb not found. create one. */
		GtkWidget* subnb=gtk_notebook_new();
		subnbs=g_slist_append(subnbs, subnb);
		nsubnb++;
		//gtk_container_set_border_width(GTK_CONTAINER(subnb),0);
#if GTK_MAJOR_VERSION>=3 || GTK_MINOR_VERSION >= 24
		gtk_notebook_set_group_name(GTK_NOTEBOOK(subnb), "secondlevel");
#elif GTK_MINOR_VERSION >=12	
		gtk_notebook_set_group(GTK_NOTEBOOK(subnb), "secondlevel");
#endif
		gtk_notebook_set_tab_pos(GTK_NOTEBOOK(subnb), GTK_POS_RIGHT);
		gtk_notebook_set_scrollable(GTK_NOTEBOOK(subnb), TRUE);
		g_signal_connect(subnb, "switch-page", G_CALLBACK(subnb_page_switch), NULL);
		GtkWidget* label=gtk_label_new(drawdata->fig);
		GtkWidget* eventbox=gtk_event_box_new();
		gtk_container_add(GTK_CONTAINER(eventbox), label);
		gtk_widget_add_events(GTK_WIDGET(eventbox), GDK_BUTTON_PRESS);
		gtk_event_box_set_visible_window(GTK_EVENT_BOX(eventbox), FALSE);
		g_signal_connect(eventbox, "button-press-event", G_CALLBACK(tab_button_cb), subnb);
		gtk_notebook_insert_page(GTK_NOTEBOOK(topnb), subnb, eventbox, itab);
		gtk_notebook_set_tab_detachable(GTK_NOTEBOOK(topnb), subnb, FALSE);
		gtk_notebook_set_tab_reorderable(GTK_NOTEBOOK(topnb), subnb, TRUE);
		gtk_widget_show_all(eventbox);
	}
	GtkWidget* page=NULL;
	GtkWidget* subnb=NULL;
	for(GSList* p=subnbs; p&&!page; p=p->next){
	/*scan through all the subnb pages with same label */
		subnb=(GtkWidget*)p->data;
		for(itab=0; itab<gtk_notebook_get_n_pages(GTK_NOTEBOOK(subnb)); itab++){
			GtkWidget* tmp=gtk_notebook_get_nth_page(GTK_NOTEBOOK(subnb), itab);
			const gchar* labeltext=subnb_label_get(subnb, tmp);
			int res=strcmp(labeltext, drawdata->name);
			if(!res){
				page=tmp;
				break;
			} else if(res>0){
				break;
			}
		}
	}
	drawdata_t* drawdata_old=NULL;
	if(page){
		/*
			we use drawdatawrap so that we don't have to modify the data on the g_object.
		*/
		drawdata_old=*(drawdata_t**)g_object_get_data(G_OBJECT(page), "drawdatawrap");
		drawdata_free_input(drawdata_old);
	}
	if(drawdata->p0){/*draw image */
		int nx=drawdata->nx;
		int ny=drawdata->ny;
		if(drawdata_old&&drawdata_old->p0&&lpf<1){
			extern int byte_float;
			if(byte_float==4){
				DO_LPF(float, drawdata_old->p0, drawdata->p0, nx*ny);
			} else if(byte_float==8){
				DO_LPF(double, drawdata_old->p0, drawdata->p0, nx*ny);
			}
		}
		size_t size=0;
		if(nx<=0||ny<=0) error("Please call DRAW_DATA\n");
		if(drawdata->gray){
			drawdata->format=(cairo_format_t)CAIRO_FORMAT_A8;
			size=1;
		} else{
			drawdata->format=(cairo_format_t)CAIRO_FORMAT_ARGB32;
			size=4;
		}
		int stride=cairo_format_stride_for_width(drawdata->format, nx);
		if(!drawdata->limit_manual){
			if(!drawdata->limit_data){
				drawdata->limit_data=mycalloc(4, float);
			}
			drawdata->limit_data[0]=-0.5;
			drawdata->limit_data[1]=drawdata->nx-0.5;
			drawdata->limit_data[2]=-0.5;
			drawdata->limit_data[3]=drawdata->ny-0.5;
		}
		/*convert data from float to int/char. */
		drawdata->p=(unsigned char*)calloc(nx*ny, size);
		flt2pix(nx, ny, !drawdata->gray, drawdata->p0, drawdata->p, drawdata->zlim);
		drawdata->image=cairo_image_surface_create_for_data
		(drawdata->p, drawdata->format, nx, ny, stride);
	}
	if(page){
		/*Instead of freeing drawdata, we replace its content with newdata. this
		  makes dialog continue to work. Do not replace generated data.*/
		drawdata_old->dtime=drawdata->time-drawdata_old->time;//compute frame rate.
		drawdata_old->time=drawdata->time;//record the time
		drawdata_old->image=drawdata->image;
		drawdata_old->p=drawdata->p;
		drawdata_old->p0=drawdata->p0;
		drawdata_old->pts=drawdata->pts;
		drawdata_old->ptsdim=drawdata->ptsdim;
		drawdata_old->npts=drawdata->npts;
		drawdata_old->nstyle=drawdata->nstyle;
		drawdata_old->style=drawdata->style;
		drawdata_old->ncir=drawdata->ncir;
		drawdata_old->cir=drawdata->cir;
		drawdata_old->fig=drawdata->fig;
		drawdata_old->name=drawdata->name;
		drawdata_old->title=drawdata->title;
		drawdata_old->xlabel=drawdata->xlabel;
		drawdata_old->ylabel=drawdata->ylabel;
		drawdata_old->legend=drawdata->legend;
		drawdata_old->grid=drawdata->grid;
		if(drawdata->limit_data){
			if(!drawdata_old->limit_data){
				drawdata_old->limit_data=mycalloc(4, float);
			}
			if(!drawdata_old->cumu){
				for(int i=0; i<4; i++){//data has different size. Reset the zoom/off.
					if(fabs(drawdata_old->limit_data[i]-drawdata->limit_data[i])>0.1){
						do_zoom(drawdata_old, 0, 0, 0);//resets zoom, offset.
					}
				}
			}
			memcpy(drawdata_old->limit_data, drawdata->limit_data, sizeof(float)*4);
		} else{
			/*free(drawdata_old->limit_data);
			  drawdata_old->limit_data=NULL;*/
		}
		drawdata_old->limit_changed=-1;
		//free(drawdata_old->limit_cumu); drawdata_old->limit_cumu=NULL;
		drawdata_old->nx=drawdata->nx;
		drawdata_old->ny=drawdata->ny;
		if(drawdata->zlim[0]||drawdata->zlim[1]){
			if(drawdata->zlim[0]<drawdata_old->zlim[0]
				||(drawdata->zlim[0]-drawdata_old->zlim[0])>0.5*fabs(drawdata_old->zlim[0])){
				//info("update zlim[0] from %g to %g\n", drawdata_old->zlim[0], drawdata->zlim[0]);
				drawdata_old->zlim[0]=drawdata->zlim[0];
			}
			if(drawdata->zlim[1]>drawdata_old->zlim[1]
				||(drawdata_old->zlim[1]-drawdata->zlim[1])>0.5*fabs(drawdata_old->zlim[1])){
					//info("update zlim[1] from %g to %g\n", drawdata_old->zlim[1], drawdata->zlim[1]);
				drawdata_old->zlim[1]=drawdata->zlim[1];
			}
		}
		drawdata_old->format=drawdata->format;
		drawdata_old->gray=drawdata->gray;
		drawdata_old->drawn=0;/*need redraw. */
		if(drawdata_old->square==-1) drawdata_old->square=drawdata->square;//otherwise keep old value.
		/*we preserve the limit instead of off, zoom in case we are drawing curves */
		//if(drawdata_old->npts){
		//    drawdata_old->limit_changed=-1;
		//	}
		{
			free(drawdata);
			LOCK(drawdata_mutex);
			ndrawdata--;
			UNLOCK(drawdata_mutex);
		}
		if(get_current_drawdata()==drawdata_old){/*we are the current page. need to update pixmap */
			update_pixmap(drawdata_old);
		} else{
			/*otherwise, notify client that it is not drawing to active page */
			page_changed(-1, -1);
		}
	} else{
		/*new tab inside the fig to contain the plot. */
		drawdata->page=page=gtk_scrolled_window_new(NULL, NULL);
		drawdata->dtime=INFINITY;
		gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(page),
			GTK_POLICY_AUTOMATIC,
			GTK_POLICY_AUTOMATIC);

		g_object_set_data(G_OBJECT(page), "drawdatawrap", drawdatawrap);
		drawdata->drawarea=drawarea=gtk_drawing_area_new();
		gtk_widget_add_events(drawarea, GDK_BUTTON_PRESS_MASK|
			GDK_BUTTON_RELEASE_MASK|
			GDK_POINTER_MOTION_MASK|GDK_POINTER_MOTION_HINT_MASK|
			GDK_BUTTON1_MOTION_MASK|
			GDK_KEY_PRESS_MASK|
			GDK_KEY_RELEASE_MASK);
#if GTK_MAJOR_VERSION>=3 || GTK_MINOR_VERSION>=22
		gtk_widget_set_can_focus(drawarea, TRUE);
		gtk_widget_set_sensitive(drawarea, TRUE);
#else
		GTK_WIDGET_SET_FLAGS(drawarea, GTK_CAN_FOCUS);
		GTK_WIDGET_SET_FLAGS(drawarea, GTK_SENSITIVE);
#endif
#if GTK_MAJOR_VERSION>=4
		GtkEventController* controller=gtk_event_controller_new(drawarea);
		g_signal_connect(controller, "motion", motion_notify, drawdatawrap);
#else
		g_signal_connect(drawarea, "motion-notify-event",
			G_CALLBACK(motion_notify), drawdatawrap);
#endif
		//notice that there is no button_press event when click through to activate a window.
		//so we also connect focus-in-event
		g_signal_connect(drawarea, "button-press-event",
			G_CALLBACK(button_press), drawdatawrap);
		g_signal_connect(drawarea, "focus-in-event",
		G_CALLBACK(focus_in_handler), drawdatawrap);
		g_signal_connect(drawarea, "button-release-event",
			G_CALLBACK(button_release), drawdatawrap);
		g_signal_connect(drawarea, "scroll-event",
			G_CALLBACK(scroll_event), drawdatawrap);
   /*g_signal_connect(drawarea,"button-release-event",
	 G_CALLBACK(button_release),drawdatawrap);*/
		g_signal_connect(drawarea, "key-press-event",
			G_CALLBACK(key_press), drawdatawrap);
#if GTK_MAJOR_VERSION<3 || GTK_MINOR_VERSION < 6
		gtk_scrolled_window_add_with_viewport
		(GTK_SCROLLED_WINDOW(page), drawarea);
#else
		gtk_container_add(GTK_CONTAINER(page), drawarea);
#endif
		GtkWidget* button=subnb_label_new(drawdatawrap);
		gtk_notebook_insert_page(GTK_NOTEBOOK(subnb), page, button, itab);
		gtk_notebook_set_tab_detachable(GTK_NOTEBOOK(subnb), page, FALSE);
		gtk_notebook_set_tab_reorderable(GTK_NOTEBOOK(subnb), page, TRUE);
#if GTK_MAJOR_VERSION>=3
		g_signal_connect(drawarea, "draw",
			G_CALLBACK(on_draw_event), drawdatawrap);
#else
		g_signal_connect(drawarea, "expose-event",
			G_CALLBACK(on_expose_event), drawdatawrap);
#endif
	/*handles zooming. */
#if GTK_MAJOR_VERSION>=4
		g_signal_connect(drawarea, "resize",
			G_CALLBACK(on_resize_event), drawdatawrap);
#else
		g_signal_connect(drawarea, "configure-event",
			G_CALLBACK(on_configure_event), drawdatawrap);
#endif
		gtk_widget_set_size_request(drawarea, DRAWAREA_MIN_WIDTH,
			DRAWAREA_MIN_HEIGHT);
		gtk_widget_show_all(subnb);
		if(get_current_drawdata()!=drawdata){
			//The new page is not activated. notify client that it is not drawing to active page
			page_changed(-1, -1);
		}
	}
	g_slist_free(subnbs);
	return 0;//return 0 cause the function to be removed frm gdb_threads_idle()
}



static void tool_save(GtkToolButton* button){
	(void)button;
	drawdata_t* drawdata=get_current_drawdata();
	static gchar* folder=NULL;
	static gchar* filename=NULL;
	GtkWidget* dialog;
	cairo_surface_t* surface=NULL;
	int width, height;
retry:
	dialog=gtk_file_chooser_dialog_new
	("Select a file to save",
		GTK_WINDOW(curwindow),
		GTK_FILE_CHOOSER_ACTION_SAVE,
		"_Cancel", GTK_RESPONSE_CANCEL,
		"_Save", GTK_RESPONSE_ACCEPT,
		NULL);
	gtk_file_chooser_set_do_overwrite_confirmation
	(GTK_FILE_CHOOSER(dialog), TRUE);
	if(filename){
		gtk_file_chooser_set_filename(GTK_FILE_CHOOSER(dialog), filename);
	} else if(folder){
		gtk_file_chooser_set_current_folder
		(GTK_FILE_CHOOSER(dialog), folder);
	} else{
		char curpath[PATH_MAX];
		if(!getcwd(curpath, PATH_MAX)){
			strcpy(curpath, HOME);
		}
		gtk_file_chooser_set_current_folder
		(GTK_FILE_CHOOSER(dialog), curpath);
	}

	if(gtk_dialog_run(GTK_DIALOG(dialog))==GTK_RESPONSE_ACCEPT){
		if(folder) g_free(folder);
		folder=gtk_file_chooser_get_current_folder(GTK_FILE_CHOOSER(dialog));
		if(filename) g_free(filename);
		filename=gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
		gtk_widget_destroy(dialog);
	} else{
		gtk_widget_destroy(dialog);
		return;
	}
	char* suffix=rindex(filename, '.');
	if(!suffix){
		char* filename2=stradd(filename, ".png", NULL);
		g_free(filename);
		filename=filename2;
		suffix=rindex(filename, '.');
	}

	if(strcmp(suffix, ".eps")==0){
#if CAIRO_HAS_PS_SURFACE == 1
		width=96*8;
		height=(drawdata)->height*96*8/(drawdata)->width;/*same aspect ratio as widget */
		surface=cairo_ps_surface_create(filename, width, height);
		cairo_ps_surface_set_eps(surface, TRUE);
#else
		error_msg("ps surface is unavailable");
		goto retry;
#endif
	} else if(strcmp(suffix, ".pdf")==0){
#if CAIRO_HAS_PDF_SURFACE == 1
		width=96*8;
		height=(drawdata)->height*96*8/(drawdata)->width;/*same aspect ratio as widget */
		surface=cairo_pdf_surface_create(filename, width, height);
#else
		error_msg("pdf surface is unavailable");
		goto retry;
#endif
	} else if(strcmp(suffix, ".svg")==0){
#if CAIRO_HAS_SVG_SURFACE == 1
		width=96*8;
		height=(drawdata)->height*96*8/(drawdata)->width;/*same aspect ratio as widget */
		surface=cairo_svg_surface_create(filename, width, height);
#else
		error_msg("svg surface is unavailable");
		goto retry;
#endif
	} else if(strcmp(suffix, ".png")==0){
		width=(drawdata)->width;/*same size as the widget. */
		height=(drawdata)->height;
		surface=cairo_image_surface_create
		((cairo_format_t)CAIRO_FORMAT_RGB24, width, height);
	} else if(strcmp(suffix, ".svg")==0){
#if CAIRO_HAS_SVG_SURFACE == 1
		width=(drawdata)->width;/*same size as the widget. */
		height=(drawdata)->height;
		surface=cairo_svg_surface_create(filename, width, height);
#else
		error_msg("svg surface is unavailable");
		goto retry;
#endif
	} else{
		error_msg("%s has unknown suffix\n", filename);
		goto retry;
	}
	if(!surface){
		error_msg("surface is NULL\n");
		return;
	}
	drawdata->dtime=100;//disable FPS
	cairo_draw(cairo_create(surface), drawdata, width, height);
	if(strcmp(suffix, ".png")==0){
		cairo_surface_write_to_png(surface, filename);
	}
	cairo_surface_finish(surface);
	cairo_surface_destroy(surface);
}


static void tool_zoom(GtkToolButton* button, gpointer data){
	(void)button;
	drawdata_t* drawdata=get_current_drawdata();
	int mode=GPOINTER_TO_INT(data);
	do_zoom(drawdata, 0, 0, mode);
}

static void limit_change(GtkSpinButton* spin, gfloat* val){
	*val=gtk_spin_button_get_value(spin);
	drawdata_dialog->limit_changed=1;
	delayed_update_pixmap(drawdata_dialog);
}

static void limit_change2(GtkSpinButton* spin, gfloat* val){
	*val=gtk_spin_button_get_value(spin);
	drawdata_dialog->limit_changed=2;
	delayed_update_pixmap(drawdata_dialog);
}
static void checkbtn_toggle(GtkToggleButton* btn, gint* key){
	*key=gtk_toggle_button_get_active(btn);
	delayed_update_pixmap(drawdata_dialog);
	/*if(key==&cumu){
		gtk_toggle_tool_button_set_active(GTK_TOGGLE_TOOL_BUTTON(toggle_cumu), cumu);
	}*/
}
static void range_changed(GtkRange* range, gfloat* val){
	*val=gtk_range_get_value(range);
	delayed_update_pixmap(drawdata_dialog);
}
static void entry_changed(GtkEditable* entry, char** key){
	free(*key);
	*key=gtk_editable_get_chars(entry, 0, -1);
	drawdata_dialog->drawn=0;
	delayed_update_pixmap(drawdata_dialog);
}
static void spin_changed(GtkSpinButton* spin, gfloat* val){
	*val=gtk_spin_button_get_value(spin);
	drawdata_dialog->drawn=0;
	delayed_update_pixmap(drawdata_dialog);
}
static void toolbutton_cumu_click(GtkToggleToolButton* btn){
	(void)btn;
	drawdata_t* page=get_current_drawdata();
	if(!page->image&&!page->square){
		cumu=gtk_toggle_tool_button_get_active(btn);
		//cumu=1-cumu;
	}
	delayed_update_pixmap(page);
}
static void toolbutton_stop(GtkToolButton* btn){
	(void)btn;
	close(sock); sock=-1; sock_idle=1;
}
static void togglebutton_pause(GtkToggleToolButton* btn){
	if(sock_idle) return;
	if(gtk_toggle_tool_button_get_active(btn)){
		if(stwriteint(sock, DRAW_PAUSE)) sock_idle=1;
	} else{
		if(stwriteint(sock, DRAW_RESUME)) sock_idle=1;
	}
}
static void togglebutton_play(GtkToggleToolButton* btn){
	(void)btn;
	if(stwriteint(sock, DRAW_SINGLE)) sock_idle=1;
}

/**
   Response to the quest to set the zaxis limit (the range of the color bar)
*/
static void tool_property(GtkToolButton* button, gpointer data){
	(void)button;
	(void)data;
	drawdata_t* drawdata=get_current_drawdata();
	if(!drawdata){
		return;
	}
	drawdata_dialog=drawdata;
	GtkWidget* dialog=gtk_dialog_new_with_buttons("Figure Properties",
		GTK_WINDOW(curwindow),
		(GtkDialogFlags)(GTK_DIALOG_DESTROY_WITH_PARENT),
		"_OK", GTK_RESPONSE_ACCEPT,
		NULL);
	GtkWidget* content_area=gtk_dialog_get_content_area(GTK_DIALOG(dialog));
	GtkWidget* vbox=gtk_vbox_new(FALSE, 0);
	GtkWidget* hbox;
	GtkWidget* checkbtn, * label, * entry, * spin;
	int n;
	GtkWidget* spins[6];
	float diff[3];
	diff[0]=(drawdata->limit0[1]-drawdata->limit0[0]);
	diff[1]=(drawdata->limit0[3]-drawdata->limit0[2]);
	if(drawdata->zlim[0]||drawdata->zlim[1]){
		diff[2]=(drawdata->zlim[1]-drawdata->zlim[0]);
		n=6;
	} else{
		n=4;
	}
	for(int i=0; i<n; i++){
	/*divide the separation to 100 steps. */
		if(diff[i/2]<EPS){
			diff[i/2]=1;
		}
		float step=pow(10, floor(log10(fabs(diff[i/2])))-2);
		spins[i]=gtk_spin_button_new_with_range(drawdata->limit0[i]-100*step, drawdata->limit0[i]+100*step, step);

		float* val;
		if(i<4){
			val=&drawdata->limit0[i];
			g_signal_connect(spins[i], "value-changed", G_CALLBACK(limit_change), val);
		} else{
			val=&drawdata->zlim[i-4];
			g_signal_connect(spins[i], "value-changed", G_CALLBACK(limit_change2), val);
		}
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(spins[i]), *val);
	}
	drawdata->spins=spins;

	checkbtn=gtk_check_button_new_with_label("Make image square");
	g_signal_connect(checkbtn, "toggled", G_CALLBACK(checkbtn_toggle), &drawdata->square);
	gtk_box_pack_start(GTK_BOX(vbox), checkbtn, FALSE, FALSE, 0);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbtn), drawdata->square);

	checkbtn=gtk_check_button_new_with_label("Enable grids");
	g_signal_connect(checkbtn, "toggled", G_CALLBACK(checkbtn_toggle), &drawdata->grid);
	gtk_box_pack_start(GTK_BOX(vbox), checkbtn, FALSE, FALSE, 0);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbtn), drawdata->grid);

	checkbtn=gtk_check_button_new_with_label("Put tic inside");
	g_signal_connect(checkbtn, "toggled", G_CALLBACK(checkbtn_toggle), &drawdata->ticinside);
	gtk_box_pack_start(GTK_BOX(vbox), checkbtn, FALSE, FALSE, 0);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbtn), drawdata->ticinside);

	checkbtn=gtk_check_button_new_with_label("Legend on curve");
	g_signal_connect(checkbtn, "toggled", G_CALLBACK(checkbtn_toggle), &drawdata->legendcurve);
	gtk_box_pack_start(GTK_BOX(vbox), checkbtn, FALSE, FALSE, 0);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbtn), drawdata->legendcurve);

	hbox=gtk_hbox_new(FALSE, 0);
	checkbtn=gtk_check_button_new_with_label("Legend.");
	g_signal_connect(checkbtn, "toggled", G_CALLBACK(checkbtn_toggle), &drawdata->legendbox);
	gtk_box_pack_start(GTK_BOX(hbox), checkbtn, FALSE, FALSE, 0);
	gtk_widget_set_sensitive(checkbtn, drawdata->legend!=NULL);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbtn), drawdata->legendbox);
	checkbtn=gtk_hscale_new_with_range(0, 1, 0.01);
	gtk_widget_set_size_request(checkbtn, 80, 20);
	gtk_scale_set_draw_value(GTK_SCALE(checkbtn), 0);
	gtk_range_set_value(GTK_RANGE(checkbtn), drawdata->legendoffx);
	g_signal_connect(checkbtn, "value-changed", G_CALLBACK(range_changed), &drawdata->legendoffx);
	gtk_box_pack_start(GTK_BOX(hbox), gtk_label_new("Position: H"), FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(hbox), checkbtn, TRUE, TRUE, 0);
	checkbtn=gtk_hscale_new_with_range(0, 1, 0.01);
	gtk_widget_set_size_request(checkbtn, 80, 20);
	gtk_scale_set_draw_value(GTK_SCALE(checkbtn), 0);
	gtk_range_set_value(GTK_RANGE(checkbtn), drawdata->legendoffy);
	g_signal_connect(checkbtn, "value-changed", G_CALLBACK(range_changed), &drawdata->legendoffy);
	gtk_box_pack_start(GTK_BOX(hbox), gtk_label_new("V"), FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(hbox), checkbtn, TRUE, TRUE, 0);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

	hbox=gtk_hbox_new(FALSE, 0);
	checkbtn=gtk_check_button_new_with_label("Plot cumulative average (");
	g_signal_connect(checkbtn, "toggled", G_CALLBACK(checkbtn_toggle), &cumu);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbtn), cumu);
	gtk_widget_set_sensitive(checkbtn, (drawdata->npts>0));
	gtk_box_pack_start(GTK_BOX(hbox), checkbtn, FALSE, FALSE, 0);
	checkbtn=gtk_check_button_new_with_label("quadrature ");
	g_signal_connect(checkbtn, "toggled", G_CALLBACK(checkbtn_toggle), &drawdata->cumuquad);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbtn), drawdata->cumuquad);
	gtk_box_pack_start(GTK_BOX(hbox), checkbtn, FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(hbox), gtk_label_new(") from"), FALSE, FALSE, 0);
	spin=gtk_spin_button_new_with_range(drawdata->limit[0], drawdata->limit[1], 1);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin), drawdata->icumu);
	g_signal_connect(spin, "value-changed", G_CALLBACK(spin_changed), &drawdata->icumu);
	gtk_box_pack_start(GTK_BOX(hbox), spin, TRUE, TRUE, 0);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

	//Set update low pass filter 
	hbox=gtk_hbox_new(FALSE, 0);
	spin=gtk_spin_button_new_with_range(0, 1, 0.01);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin), lpf);
	g_signal_connect(spin, "value-changed", G_CALLBACK(spin_changed), &lpf);
	gtk_box_pack_start(GTK_BOX(hbox), spin, TRUE, TRUE, 0);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

	hbox=gtk_hbox_new(FALSE, 0);
	entry=gtk_entry_new();
	label=gtk_label_new("Title"); gtk_label_set_width_chars(GTK_LABEL(label), 6);
	gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(hbox), entry, TRUE, TRUE, 0);
	gtk_entry_set_text(GTK_ENTRY(entry), drawdata->title);
	g_signal_connect(GTK_EDITABLE(entry), "changed", G_CALLBACK(entry_changed), &drawdata->title);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

	hbox=gtk_hbox_new(FALSE, 0);
	entry=gtk_entry_new();
	label=gtk_label_new("X label");gtk_label_set_width_chars(GTK_LABEL(label), 6);
	gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(hbox), entry, TRUE, TRUE, 0);
	gtk_entry_set_text(GTK_ENTRY(entry), drawdata->xlabel);
	g_signal_connect(GTK_EDITABLE(entry), "changed", G_CALLBACK(entry_changed), &drawdata->xlabel);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

	hbox=gtk_hbox_new(FALSE, 0);
	entry=gtk_entry_new();
	label=gtk_label_new("Y label");gtk_label_set_width_chars(GTK_LABEL(label), 6);
	gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(hbox), entry, TRUE, TRUE, 0);
	gtk_entry_set_text(GTK_ENTRY(entry), drawdata->ylabel);
	g_signal_connect(GTK_EDITABLE(entry), "changed", G_CALLBACK(entry_changed), &drawdata->ylabel);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);


	hbox=gtk_hbox_new(FALSE, 0);
	label=gtk_label_new("xmin");gtk_label_set_width_chars(GTK_LABEL(label), 5);
	gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(hbox), spins[0], TRUE, TRUE, 0);
	label=gtk_label_new("xmax");gtk_label_set_width_chars(GTK_LABEL(label), 5);
	gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(hbox), spins[1], TRUE, TRUE, 0);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

	hbox=gtk_hbox_new(FALSE, 0);
	label=gtk_label_new("ymin");gtk_label_set_width_chars(GTK_LABEL(label), 5);
	gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(hbox), spins[2], TRUE, TRUE, 0);
	label=gtk_label_new("ymax");gtk_label_set_width_chars(GTK_LABEL(label), 5);
	gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(hbox), spins[3], TRUE, TRUE, 0);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

	if(n>4){
		hbox=gtk_hbox_new(FALSE, 0);
		label=gtk_label_new("zmin");gtk_label_set_width_chars(GTK_LABEL(label), 5);
		gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
		gtk_box_pack_start(GTK_BOX(hbox), spins[4], TRUE, TRUE, 0);
		label=gtk_label_new("zmax");gtk_label_set_width_chars(GTK_LABEL(label), 5);
		gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
		gtk_box_pack_start(GTK_BOX(hbox), spins[5], TRUE, TRUE, 0);
		gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
	}

	gtk_container_add(GTK_CONTAINER(content_area), vbox);
	gtk_widget_show_all(vbox);
	gtk_dialog_run(GTK_DIALOG(dialog));
	drawdata->spins=NULL;
	gtk_widget_destroy(dialog);
}


static void tool_font_set(GtkFontButton* btn){
#if GTK_MAJOR_VERSION >=3
	const char* font_name_new=gtk_font_chooser_get_font(GTK_FONT_CHOOSER(btn));
#else	
	const char* font_name_new=gtk_font_button_get_font_name(btn);
#endif
	desc=pango_font_description_from_string(font_name_new);
	PangoFontDescription* pfd
		=pango_font_description_from_string(font_name_new);
	font_name_version++;/*tell every expose event to update the figure. */
	if(font_name) free(font_name);
	font_name=strdup(pango_font_description_get_family(pfd));

	PangoStyle style=pango_font_description_get_style(pfd);
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
	}
	PangoWeight weight=pango_font_description_get_weight(pfd);
	if(weight<=500){
		font_weight=CAIRO_FONT_WEIGHT_NORMAL;
	} else{
		font_weight=CAIRO_FONT_WEIGHT_BOLD;
	}
	int size=pango_font_description_get_size(pfd);
	/*get font_size in device unit (dots, pixels); */
	if(pango_font_description_get_size_is_absolute(pfd)){
		font_size=(float)size/(float)PANGO_SCALE;
		fprintf(stderr, "absolute: font_size=%g\n", font_size);
	} else{
		gfloat dpi=gdk_screen_get_resolution(gdk_screen_get_default());
		if(dpi<0) dpi=96;
		font_size=(float)size/(float)PANGO_SCALE*(float)dpi/72.;
	}
	SP_XL=font_size*2.4+8;
	SP_YT=font_size*1.3+12;
	SP_YB=font_size*2.4+12;
	SP_XR=SP_LEG+LEN_LEG+font_size*2;
	pango_font_description_free(pfd);
}
#if GTK_MAJOR_VERSION>=3
static gboolean close_window(GtkWidget* object, GdkEvent* event)
#else
static gboolean close_window(GtkObject* object, GdkEvent* event)
#endif
{
	(void)event;
	windows=g_slist_remove(windows, object);
	info("close_window called\n");

	GtkWidget* topnb=get_topnb(GTK_WIDGET(object));
	GtkWidget* toolbar=get_toolbar(GTK_WIDGET(object));
	g_signal_handlers_disconnect_by_func(topnb, (gpointer)topnb_page_changed, toolbar);
	int npages=gtk_notebook_get_n_pages(GTK_NOTEBOOK(topnb));
	for(int ipage=0; ipage<npages; ipage++){
		GtkWidget* page=gtk_notebook_get_nth_page(GTK_NOTEBOOK(topnb), ipage);
		int ntabs=gtk_notebook_get_n_pages(GTK_NOTEBOOK(page));
		for(int itab=0; itab<ntabs; itab++){
			GtkWidget* tab=gtk_notebook_get_nth_page(GTK_NOTEBOOK(page), itab);
			drawdata_t** drawdatawrap=(drawdata_t**)g_object_get_data(G_OBJECT(tab), "drawdatawrap");
			gtk_widget_hide(tab);
			drawdata_free(*drawdatawrap);
		}
	}
	gtk_widget_destroy(object);
	if(!windows){
		info("sock %d closed\n", sock);
		if(sock!=-1) close(sock);

		gtk_main_quit();
	}
	/*don't call exit here. */
	return FALSE;
}

static gboolean window_state(GtkWidget* window, GdkEvent* event){
	if(event->type==GDK_FOCUS_CHANGE){
		GdkEventFocus* focus=(GdkEventFocus*)event;
		if(focus->in){
			set_cur_window(window);
		}
	}
	return FALSE;
}
//Create a new toolbar item.
#if GTK_MAJOR_VERSION < 4
void new_tool(GtkWidget* toolbar, GtkWidget* child, int toggle, const char* name,
	GCallback callback, gpointer user_data){
	GtkToolItem* item=NULL;

	if(name){
		if(child){
			item=gtk_tool_item_new();
		} else if(toggle){
			item=gtk_toggle_tool_button_new();
		} else{
			item=gtk_tool_button_new(NULL, name);
		}
		if(child){
			gtk_container_add(GTK_CONTAINER(item), child);
		} else{
			gtk_tool_button_set_icon_name(GTK_TOOL_BUTTON(item), name);
		}
	} else{
		item=gtk_separator_tool_item_new();
	}

	gtk_toolbar_insert(GTK_TOOLBAR(toolbar), item, -1);
	if(callback){
		g_signal_connect(item, toggle?"toggled":"clicked", G_CALLBACK(callback), user_data);
	}
}
#else
void new_tool(GtkWidget* toolbar, GtkWidget* item, int toggle, const char* name,
	GCallback callback, gpointer user_data){
	if(!item){
		if(name){
			if(toggle){
				item=gtk_button_new_from_icon_name(name);
			} else{
				item=gtk_toggle_button_new();
				gtk_button_set_icon_name(GTK_BUTTON(item), name);
			}
		} else{
			item=gtk_button_new();
		}
	}
	gtk_button_set_has_frame(GTK_BUTTON(item), false);
	gtk_box_append(GTK_BOX(toolbar), item);

	if(callback){
		g_signal_connect(item, toggle?"toggled":"clicked", G_CALLBACK(callback), user_data);
	}
}
#endif
GtkWidget* create_window(GtkWidget* window){
	if(!cursors[1]){
		GdkDisplay* display=gdk_display_get_default();
		cursors[0]=gdk_cursor_new_from_name(display, "default");
		cursors[1]=gdk_cursor_new_from_name(display, "crosshair");
#if GTK_MAJOR_VERSION < 3
		gtk_rc_parse_string(rc_string_notebook);
#endif
	}
#if GTK_MAJOR_VERSION < 4
	if(!contexmenu){
		contexmenu=gtk_menu_new();
		GtkWidget* item=gtk_menu_item_new_with_label("Detach");
		g_signal_connect(item, "activate", G_CALLBACK(topnb_detach), NULL);
		gtk_menu_shell_append(GTK_MENU_SHELL(contexmenu), item);
		gtk_widget_show_all(contexmenu);
		g_object_ref(contexmenu);
	}
#endif
	if(!window){
		window=gtk_window_new(GTK_WINDOW_TOPLEVEL);
	}
	//gtk_window_stick(GTK_WINDOW(window));
	curwindow=window;
	windows=g_slist_append(windows, window);
	extern GdkPixbuf* icon_main;
	gtk_window_set_icon(GTK_WINDOW(window), icon_main);
	g_object_unref(icon_main);

	char title[80];
	if(iwindow==0){
		snprintf(title, 80, "MAOS Drawdaemon");
	} else{
		snprintf(title, 80, "MAOS Drawdaemon (%d)", iwindow);
	}
	gtk_window_set_title(GTK_WINDOW(window), title);
#if GTK_MAJOR_VERSION>=3
	GtkCssProvider* provider_default=gtk_css_provider_new();
	gtk_css_provider_load_from_data(provider_default, all_style, strlen(all_style), NULL);
	GtkStyleContext* all_context=gtk_widget_get_style_context(window);
	GdkScreen* screen=gtk_style_context_get_screen(all_context);
	gtk_style_context_add_provider_for_screen(screen, GTK_STYLE_PROVIDER(provider_default),
		GTK_STYLE_PROVIDER_PRIORITY_USER);
#endif
#if GTK_MAJOR_VERSION<3
	if(!btn_rcstyle){
		btn_rcstyle=gtk_rc_style_new();
	//This option is used in many places to determine padding between text and border of widgets
		btn_rcstyle->xthickness=btn_rcstyle->ythickness=0;
	}
#endif
#if GTK_MAJOR_VERSION<4
	GtkWidget* toolbar=gtk_toolbar_new();
	gtk_toolbar_set_icon_size(GTK_TOOLBAR(toolbar), GTK_ICON_SIZE_MENU);
	gtk_toolbar_set_style(GTK_TOOLBAR(toolbar), GTK_TOOLBAR_ICONS);
#else
	//Toolbar is deprecated in GTK4. Use box instead
	GtkWidget* toolbar=gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
#endif
	gtk_widget_set_sensitive(toolbar, FALSE);
	new_tool(toolbar, NULL, 0, "document-save-as", G_CALLBACK(tool_save), NULL);
	new_tool(toolbar, NULL, 0, NULL, NULL, NULL); //separator
	new_tool(toolbar, NULL, 0, "zoom-in", G_CALLBACK(tool_zoom), GINT_TO_POINTER(1));
	new_tool(toolbar, NULL, 0, "zoom-fit-best", G_CALLBACK(tool_zoom), GINT_TO_POINTER(0));
	new_tool(toolbar, NULL, 0, "zoom-out", G_CALLBACK(tool_zoom), GINT_TO_POINTER(-1));
	new_tool(toolbar, NULL, 0, NULL, NULL, NULL); //separator
	new_tool(toolbar, NULL, 0, "document-properties", G_CALLBACK(tool_property), NULL);
	new_tool(toolbar, NULL, 1, "edit-copy", G_CALLBACK(toolbutton_cumu_click), NULL);
	new_tool(toolbar, NULL, 0, NULL, NULL, NULL); //separator

	GtkWidget* fontsel=gtk_font_button_new_with_font("Sans 12");
	g_signal_connect(GTK_FONT_BUTTON(fontsel), "font-set", G_CALLBACK(tool_font_set), NULL);
	new_tool(toolbar, fontsel, 0, "font-set", NULL, NULL);
	new_tool(toolbar, NULL, 1, "media-playback-play", G_CALLBACK(togglebutton_play), NULL);
	new_tool(toolbar, NULL, 1, "media-playback-pause", G_CALLBACK(togglebutton_pause), NULL);
	new_tool(toolbar, NULL, 0, "media-playback-stop", G_CALLBACK(toolbutton_stop), NULL);

	GtkWidget* topnb=gtk_notebook_new();
	//gtk_container_set_border_width(GTK_CONTAINER(topnb),2);
#if GTK_MAJOR_VERSION>=3 || GTK_MINOR_VERSION >= 24
	gtk_notebook_set_group_name(GTK_NOTEBOOK(topnb), "toplevel");
#elif GTK_MINOR_VERSION >= 12
	gtk_notebook_set_group(GTK_NOTEBOOK(topnb), "toplevel");
#endif
	gtk_notebook_set_scrollable(GTK_NOTEBOOK(topnb), TRUE);
	g_signal_connect(GTK_NOTEBOOK(topnb), "page-added", G_CALLBACK(topnb_page_changed), toolbar);
	g_signal_connect(GTK_NOTEBOOK(topnb), "page-removed", G_CALLBACK(topnb_page_changed), toolbar);
	g_signal_connect(GTK_NOTEBOOK(topnb), "switch-page", G_CALLBACK(topnb_page_switch), toolbar);
	GtkWidget* vbox=gtk_vbox_new(FALSE, 0);
	gtk_box_pack_start(GTK_BOX(vbox), toolbar, FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(vbox), topnb, TRUE, TRUE, 0);
	gtk_container_add(GTK_CONTAINER(window), vbox);
	g_signal_connect(window, "delete-event", G_CALLBACK(close_window), NULL);
	gtk_window_set_position(GTK_WINDOW(window), GTK_WIN_POS_CENTER);
	gtk_window_set_default_size(GTK_WINDOW(window), 1050, 900);
	gtk_widget_add_events(GTK_WIDGET(window), GDK_FOCUS_CHANGE_MASK);/*in case this is not on. */
	set_cur_window(window);
	g_signal_connect(GTK_WIDGET(window), "event", G_CALLBACK(window_state), NULL);
	gtk_widget_show_all(window);
	gtk_window_present(GTK_WINDOW(window));
	tool_font_set(GTK_FONT_BUTTON(fontsel));/*initialize. */
	iwindow++;
	return window;
}

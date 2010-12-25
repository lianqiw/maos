#include <glib/gprintf.h>
#include <gdk/gdkkeysyms.h>
#include "drawdaemon.h"
#include "icon-draw.h"
#include "mouse_hand.h"
#include "mouse_white.h"

#define DRAWAREA_MIN_WIDTH 440
#define DRAWAREA_MIN_HEIGHT 320
static GSList *windows=NULL;
static int iwindow=0;
static GtkWidget *curwindow=NULL;
static GtkWidget *curtopnb=NULL;
static GtkWidget *topmenu=NULL;
int font_name_version=0;
char *font_name=NULL;
double font_size=9;
cairo_font_slant_t font_style=CAIRO_FONT_SLANT_NORMAL;
cairo_font_weight_t font_weight=CAIRO_FONT_WEIGHT_NORMAL;
static int cursor_type=0;//cursor type of the drawing area.
/*
  Routines in this file are about the GUI.
*/

GdkCursor *cursors[2];
GdkPixbuf *pix_hand=NULL;
GdkPixbuf *pix_arrow=NULL;

static const char *rc_string_notebook={
    "style \"noborder\"{                      \n"
    "GtkNoteBook::draw-border={0 0 0 0}       \n"
    "}                                        \n"
    "class \"GtkNoteBook\" style \"noborder\" \n"
};


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

static void close_window(GtkObject *object);

static void set_cur_window(GtkWidget *window){
    curwindow=window;
    GtkWidget *vbox=gtk_bin_get_child(GTK_BIN(curwindow));
    GList *list=gtk_container_get_children(GTK_CONTAINER(vbox));
    curtopnb=list->next->data;
    g_list_free(list);
}

static drawdata_t *get_current_page(void){
    GtkWidget *topnb=curtopnb;
    int p1=gtk_notebook_get_current_page(GTK_NOTEBOOK(topnb));
    GtkWidget *w1=gtk_notebook_get_nth_page(GTK_NOTEBOOK(topnb),p1);//second notebook
    int p2=gtk_notebook_get_current_page(GTK_NOTEBOOK(w1));
    GtkWidget *w2=gtk_notebook_get_nth_page(GTK_NOTEBOOK(w1),p2);//scrolled window
    drawdata_t **pdrawdata=g_object_get_data(G_OBJECT(w2),"drawdatawrap");
    if(pdrawdata){
	return *pdrawdata;
    }else{
	return NULL;
    }
}
static GtkWidget *get_topnb(GtkWidget *window){
    GtkWidget *vbox=gtk_bin_get_child(GTK_BIN(window));
    GList *list=gtk_container_get_children(GTK_CONTAINER(vbox));
    GtkWidget *topnb=list->next->data;
    g_list_free(list);
    return topnb;
}
static GtkWidget *get_toolbar(GtkWidget *window){
    GtkWidget *vbox=gtk_bin_get_child(GTK_BIN(window));
    GList *list=gtk_container_get_children(GTK_CONTAINER(vbox));
    GtkWidget *toolbar=list->data;
    g_list_free(list);
    return toolbar;
}
/*
static void drag_end(GtkWidget *widget, GdkDragContext *drag_context, gpointer data){
    info("Drag end on %p\n", widget);
}

static gboolean drag_failed(GtkWidget *widget, GdkDragContext *drag_context, 
			    GtkDragResult result,
			    gpointer data){
    info("Drag failed on %p\n", widget);
    return FALSE;
    }*/
static void topnb_page_added(GtkNotebook *topnb, GtkWidget *child, guint n, GtkWidget *toolbar){
    (void)topnb;
    (void)child;
    (void)n;
    gtk_widget_set_sensitive(toolbar, TRUE);
}
static void topnb_page_removed(GtkNotebook *topnb, GtkWidget *child, guint n, GtkWidget *toolbar){
    (void)child;
    (void)n;
    if(gtk_notebook_get_n_pages(GTK_NOTEBOOK(topnb))==0){
	if(g_slist_length(windows)>1){
	    GtkWidget *window=gtk_widget_get_parent
		(gtk_widget_get_parent(GTK_WIDGET(topnb)));
	    gtk_object_destroy(GTK_OBJECT(window));
	}else{
	    gtk_widget_set_sensitive(toolbar, FALSE);
	}
    }
}
static void topnb_detach(GtkMenuItem *menu){
    (void)menu;
    GtkWidget *topnb=curtopnb;
    int n=gtk_notebook_get_current_page(GTK_NOTEBOOK(topnb));
    GtkWidget *page=gtk_notebook_get_nth_page(GTK_NOTEBOOK(topnb),n);
    GtkWidget *label=gtk_notebook_get_tab_label(GTK_NOTEBOOK(topnb),page);
    g_object_ref(page);
    g_object_ref(label);
    gtk_notebook_remove_page(GTK_NOTEBOOK(topnb), n);
    GtkWidget *window=create_window();//create a new window.
    topnb=get_topnb(window);
    gtk_notebook_append_page(GTK_NOTEBOOK(topnb), page, label);
    gtk_notebook_set_tab_detachable(GTK_NOTEBOOK(topnb), page, TRUE);
    gtk_notebook_set_tab_reorderable(GTK_NOTEBOOK(topnb), page, TRUE);
    g_object_unref(page);
    g_object_unref(label);
}
static gboolean tab_button_cb(GtkWidget *widget, GdkEventButton *event, GtkWidget *page){
    /*
      widget is the event box that is the page label. page is the page in the notebook.
     */
    info("button press\n");
    GtkWidget *topnb=gtk_widget_get_parent(page);
    int n=gtk_notebook_page_num(GTK_NOTEBOOK(topnb), page);
    gtk_notebook_set_current_page(GTK_NOTEBOOK(topnb), n);
    if(event->button==3){
	if(gtk_menu_get_attach_widget(GTK_MENU(topmenu))){
	    gtk_menu_detach(GTK_MENU(topmenu));
	}
	gtk_menu_attach_to_widget(GTK_MENU(topmenu), widget, NULL);
	gtk_menu_popup(GTK_MENU(topmenu),NULL,NULL,NULL,NULL,3,gtk_get_current_event_time ());
    }
    return FALSE;
}
typedef struct updatetimer_t{
    int pending;
    drawdata_t *drawdata;
}updatetimer_t;
static void update_pixmap(drawdata_t *drawdata){
    //no more pending updates, do the updating.
    gint width=drawdata->width;
    gint height=drawdata->height;
    gint width2, height2;
    if(drawdata->pixmap){
	gdk_drawable_get_size(drawdata->pixmap, &width2, &height2);
	if(width!=width2 || height !=height2){
	    info("Replacing pixmap\n");
	    g_object_unref(drawdata->pixmap);
	    drawdata->pixmap=NULL;
	}else{
	    info("Keep old pixmap\n");
	}
    }
    if(!drawdata->pixmap){
	//Create a new server size pixmap and then draw on it.
	drawdata->pixmap=gdk_pixmap_new(curwindow->window, width, height, -1);
    }
    cairo_draw(gdk_cairo_create(drawdata->pixmap), drawdata,width,height);
    gtk_widget_queue_draw(drawdata->drawarea);
}
static gboolean update_pixmap_timer(gpointer timer0){
    updatetimer_t *timer=(updatetimer_t*)timer0;
    drawdata_t *drawdata=timer->drawdata;
    if(timer->pending==drawdata->pending){
	update_pixmap(drawdata);
    }
    free(timer);
    return FALSE;
}
static void delayed_update_pixmap(drawdata_t *drawdata){
    drawdata->pending++;
    updatetimer_t *tmp=calloc(1, sizeof(updatetimer_t));
    tmp->pending=drawdata->pending;
    tmp->drawdata=drawdata;
    if(!drawdata->pixmap){
	update_pixmap(drawdata);
    }else{
	g_timeout_add(100, update_pixmap_timer, tmp);
    }
}


/* 
   Expose event happens when the widget is first configured and whenever it is
   resized.  Create a new backup pixmap of the appropriate size and draw there.*/
static gboolean
on_configure_event(GtkWidget *widget, GdkEventConfigure *event, gpointer pdata){
    (void)event; 
    (void)widget;
    drawdata_t* drawdata=*((drawdata_t**)pdata);
    //info("configure event :%dx%d. event:%dx%d\n", 
    //drawdata->width, drawdata->height, event->width, event->height);
    if(event->width>1 && event->height>1){
	drawdata->width=event->width;
	drawdata->height=event->height;
	delayed_update_pixmap(drawdata);
    }
    return FALSE;
}
/* Redraw the screen from the backing pixmap */
static gboolean
on_expose_event(GtkWidget *widget,
		GdkEventExpose *event,
		gpointer pdata){
    drawdata_t* drawdata=*((drawdata_t**)pdata);
    if(drawdata->font_name_version != font_name_version || !drawdata->drawn){
	update_pixmap(drawdata);
    }
    //info("expose event\n");
    if(drawdata->pixmap){
	gdk_draw_drawable(widget->window, 
			  widget->style->fg_gc[GTK_WIDGET_STATE(widget)],
			  drawdata->pixmap,
			  event->area.x, event->area.y,
			  event->area.x, event->area.y,
			  event->area.width, event->area.height);
    }
    return FALSE;
}

static void drawdata_free(drawdata_t *drawdata){
    cairo_surface_destroy(drawdata->image);
    if(drawdata->pixmap){
	g_object_unref(drawdata->pixmap);
	drawdata->pixmap=NULL;
    }

    if(drawdata->p0){
	free(drawdata->p0);
    }
    if(drawdata->p){
	free(drawdata->p);
    }
    if(drawdata->npts>0){
	for(int ipts=0; ipts<drawdata->npts; ipts++){	
	    dfree(drawdata->pts[ipts]);
	}
	free(drawdata->pts);
    }
    if(drawdata->nstyle>0){
	free(drawdata->style);
    }
    if(drawdata->ncir>0){
	free(drawdata->cir);
    }
    if(drawdata->fig!=deffig) free(drawdata->fig);
    if(drawdata->name!=defname) free(drawdata->name);
    if(drawdata->title) free(drawdata->title);
    if(drawdata->xlabel) free(drawdata->xlabel);
    if(drawdata->ylabel) free(drawdata->ylabel);
    free(drawdata->limit);
    free(drawdata);
    pthread_mutex_lock(&mutex_drawdata);
    ndrawdata--;
    info("drawdata deleted, ndrawdata=%d\n", ndrawdata);
    pthread_mutex_unlock(&mutex_drawdata);
}

/**
   Delete a figure page.
*/
static void delete_page(GtkButton *btn, drawdata_t **drawdatawrap){
    (void)btn;
    GtkWidget *root;
    //First find the root page
    GtkWidget *topnb=curtopnb;
    root=gtk_notebook_get_nth_page(GTK_NOTEBOOK(topnb), 
				   gtk_notebook_get_current_page (GTK_NOTEBOOK(topnb)));
    //Find the sub page. it may not be the current page
    int ipage=gtk_notebook_page_num(GTK_NOTEBOOK(root), (*drawdatawrap)->page);
    gtk_notebook_remove_page(GTK_NOTEBOOK(root), ipage);
    drawdata_free(*drawdatawrap);
    free(drawdatawrap);
    if(gtk_notebook_get_n_pages(GTK_NOTEBOOK(root))==0){
	int jpage=gtk_notebook_page_num(GTK_NOTEBOOK(topnb), root);
	gtk_notebook_remove_page(GTK_NOTEBOOK(topnb),jpage);
    }
}

#if ! defined(__linux__)
/**
   Delete a figure page using button-press-event. button works not good in mac.
*/
static void delete_page_event(GtkWidget *widget, GdkEventButton *event, drawdata_t **drawdatawrap){
    (void)widget;
    (void)event;
    if(event->button==1 && event->type==GDK_BUTTON_PRESS){
	delete_page(NULL, drawdatawrap);
    }
}
#endif
static GtkWidget *tab_label_new(drawdata_t **drawdatawrap){
    drawdata_t *drawdata=*drawdatawrap;
    const gchar *str=drawdata->name;
    GtkWidget *label;
    GtkWidget *out;
    out=gtk_hbox_new(FALSE, 0);
#if defined(__linux__)
    GtkWidget *close_btn;
    close_btn=gtk_button_new();
    gtk_button_set_focus_on_click(GTK_BUTTON(close_btn), FALSE);
    gtk_container_add(GTK_CONTAINER(close_btn), 
		      gtk_image_new_from_stock(GTK_STOCK_CLOSE, GTK_ICON_SIZE_MENU));
    g_signal_connect(close_btn, "clicked", G_CALLBACK(delete_page), drawdatawrap);
    /* make button as small as possible */

    GtkRcStyle *rcstyle;
    gtk_button_set_relief(GTK_BUTTON(close_btn), GTK_RELIEF_NONE);
    rcstyle = gtk_rc_style_new();
    rcstyle->xthickness = rcstyle->ythickness = 0;
    gtk_widget_modify_style(close_btn, rcstyle);
    gtk_widget_set_size_request(close_btn,14,14);
    gtk_rc_style_unref(rcstyle);
    gtk_box_pack_end(GTK_BOX(out), close_btn, FALSE, FALSE, 0);
#else
    GtkWidget *ebox=gtk_event_box_new();
    GtkWidget *image=gtk_image_new_from_stock(GTK_STOCK_CLOSE, GTK_ICON_SIZE_MENU);
    gtk_container_add(GTK_CONTAINER(ebox), image);
    g_signal_connect(ebox, "button_press_event", G_CALLBACK(delete_page_event), drawdatawrap);
    gtk_box_pack_end(GTK_BOX(out), ebox, FALSE, FALSE, 0);
#endif

    /* create label for tab */
    label = gtk_label_new(str);
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
    gtk_misc_set_padding(GTK_MISC(label), 0, 0);
    gtk_label_set_width_chars(GTK_LABEL(label), 12);
    gtk_box_pack_start(GTK_BOX(out), label, TRUE, TRUE, 0);
    gtk_widget_show_all(out);
    return out;
}
static void do_move(drawdata_t *drawdata, double xdiff, double ydiff){
    drawdata->offx+=xdiff/drawdata->zoomx;
    drawdata->offy+=ydiff/drawdata->zoomy;
    update_pixmap(drawdata);//no need delay since motion notify already did it.
}

static gboolean motion_notify(GtkWidget *widget, GdkEventMotion *event, 
			      drawdata_t **drawdatawrap){
    (void)widget;
    drawdata_t *drawdata=*drawdatawrap;
    if(event->state & GDK_BUTTON1_MASK && drawdata->valid){//move with left cursor
	double x, y;
	x = event->x;
	y = event->y;
	double dx, dy;
	dx = x - drawdata->mxdown;
	dy = y - drawdata->mydown;

	if(cursor_type==0){//move
	    do_move(drawdata, dx, -dy);//notice the reverse sign.
	    drawdata->mxdown=x;
	    drawdata->mydown=y;
	}else{//select and zoom.
	    if(drawdata->square){//for a square
		if(fabs(dx)<fabs(dy)){
		    dy*=fabs(dx/dy);
		}else{
		    dx*=fabs(dy/dx);
		}
	    }
	    if(drawdata->pixmap){//force a refresh to remove previous rectangule
		gdk_draw_drawable(widget->window, 
				  widget->style->fg_gc[GTK_WIDGET_STATE(widget)],
				  drawdata->pixmap,
				  0,0,0,0,-1,-1);
	    }
	    //do not draw to pixmap, use window
	    cairo_t *cr=gdk_cairo_create(widget->window);
	    cairo_set_antialias(cr,CAIRO_ANTIALIAS_NONE);
	    cairo_set_source_rgba(cr,0,0,1,0.1);
	    cairo_set_line_width(cr, 1);
	    cairo_rectangle(cr, drawdata->mxdown, drawdata->mydown, dx, dy);
	    cairo_fill_preserve(cr);
	    cairo_set_source_rgba(cr, 0,0,0,1);
	    cairo_stroke(cr);
	    cairo_destroy(cr);
	
	}
    }
    //we set the cursor
    if(event->x > drawdata->xoff && event->x < drawdata->xoff + drawdata->widthim 
       && event->y > drawdata->yoff && event->y < drawdata->yoff + drawdata->heightim){
	if(!drawdata->cursorinside){
	    drawdata->cursorinside=1;
	    gdk_window_set_cursor(gtk_widget_get_window(widget), cursors[cursor_type]);
	}
    }else{
	if(drawdata->cursorinside){
	    drawdata->cursorinside=0;
	    gdk_window_set_cursor(gtk_widget_get_window(widget), NULL);
	}
    }
    return FALSE;
}


static void do_zoom(drawdata_t *drawdata, double xdiff, double ydiff, int mode){
    double old_zoomx=drawdata->zoomx;
    double old_zoomy=drawdata->zoomy;
    if(mode==1){//zoom in
	drawdata->zoomx*=1.2;
	drawdata->zoomy*=1.2;
    }else if(mode==-1){
	drawdata->zoomx/=1.2;
	drawdata->zoomy/=1.2;
    }else if(mode==0){//reset everything
	drawdata->zoomx=1;
	drawdata->zoomy=1;
	drawdata->offx=0;
	drawdata->offy=0;
    }
    if(drawdata->zoomx<0.001)
	drawdata->zoomx=0.001;
    else if(drawdata->zoomx>1000){
	drawdata->zoomx=1000;
    }
    if(drawdata->zoomy<0.001)
	drawdata->zoomy=0.001;
    else if(drawdata->zoomy>1000){
	drawdata->zoomy=1000;
    }
    if(mode){//not zero.
	double factorx=1/old_zoomx-1/drawdata->zoomx;
	drawdata->offx-=xdiff*factorx;
	double factory=1/old_zoomy-1/drawdata->zoomy;
	drawdata->offy-=ydiff*factory;
    }
    update_pixmap(drawdata);
}

static gboolean scroll_event(GtkWidget *widget, GdkEventScroll *event, 
			     drawdata_t **drawdatawrap){
    (void)widget;
    drawdata_t *drawdata=*drawdatawrap;
    double xdiff=event->x-drawdata->centerx;
    double ydiff=-(event->y-drawdata->centery);
    if(event->direction == GDK_SCROLL_UP){
	do_zoom(drawdata,xdiff,ydiff,1);
    }else if(event->direction == GDK_SCROLL_DOWN){
	do_zoom(drawdata,xdiff,ydiff,-1);
    }
    return FALSE;
}
static gboolean button_press(GtkWidget *widget, GdkEventButton *event, drawdata_t **drawdatawrap){
    drawdata_t *drawdata=*drawdatawrap;
    //Grab focus so the keys work
    if(!GTK_WIDGET_HAS_FOCUS(widget))
	gtk_widget_grab_focus(widget);

    /*if(event->type==GDK_2BUTTON_PRESS && event->button==1){
    //double click brings the part of the image to center.
    do_move(drawdata,
    (drawdata->centerx-event->x),
    -(drawdata->centery-event->y));//notice the reverse sign.
    }else{*/
    if(event->x > drawdata->xoff && event->x < drawdata->xoff + drawdata->widthim 
       && event->y > drawdata->yoff && event->y < drawdata->yoff + drawdata->heightim){
	drawdata->mxdown=event->x;
	drawdata->mydown=event->y;
	drawdata->valid=1;
    }else{
	drawdata->valid=0;
    }
    //}
    return FALSE;
}


static gboolean button_release(GtkWidget *widget, GdkEventButton *event, drawdata_t **drawdatawrap){
    (void) widget;
    drawdata_t *drawdata=*drawdatawrap;
    if(!drawdata->valid) return FALSE;
    double x, y;
    x = event->x;
    y = event->y;
    if(cursor_type==0 && event->button==1){//move only on left button
	double dx = x - drawdata->mxdown;
	double dy = y - drawdata->mydown;
	do_move(drawdata, dx, -dy);
    }else if(cursor_type==1){//select and zoom.
	double xx = drawdata->mxdown;
	double dx = x - drawdata->mxdown;
	double dy = y - drawdata->mydown;
	if(drawdata->square){
	    if(fabs(dx)<fabs(dy)){
		dy*=fabs(dx/dy);
	    }else{
		dx*=fabs(dy/dx);
	    }
	}
	if(dx<0) xx+=dx;
	double yy=drawdata->mydown;
	if(dy>0) yy+=dy;
	double diffx=(drawdata->limit0[1]-drawdata->limit0[0])/drawdata->widthim;
	double diffy=(drawdata->limit0[3]-drawdata->limit0[2])/drawdata->heightim;
	drawdata->limit0[0]+=diffx*(xx-drawdata->xoff);
	drawdata->limit0[1]=drawdata->limit0[0]+diffx*fabs(dx);
	drawdata->limit0[2]+=diffy*(drawdata->yoff+drawdata->heightim-yy);
	drawdata->limit0[3]=drawdata->limit0[2]+diffy*fabs(dy);
	apply_limit(drawdata);
	update_pixmap(drawdata);
    }
    return FALSE;
}


static void switch_tab(int lr, int ud){
    GtkWidget *topnb=curtopnb;
    if(lr){
	gtk_notebook_set_current_page
	    (GTK_NOTEBOOK(topnb),
	     gtk_notebook_get_current_page(GTK_NOTEBOOK(topnb))+lr);
    }else if(ud){
	GtkWidget *page=gtk_notebook_get_nth_page
	    (GTK_NOTEBOOK(topnb),
	     gtk_notebook_get_current_page(GTK_NOTEBOOK(topnb)));
	gtk_notebook_set_current_page
	    (GTK_NOTEBOOK(page),
	     gtk_notebook_get_current_page(GTK_NOTEBOOK(page))+ud);
    }
}
static gboolean key_press(GtkWidget *widget, GdkEventKey *event, 
			  drawdata_t **drawdatawrap){
    (void)widget;

    if(event->state & GDK_CONTROL_MASK){
	switch(event->keyval){
	case GDK_Left:
	    switch_tab(-1,0);break;
	case GDK_Right:
	    switch_tab(1,0);break;
	case GDK_Up:
	    switch_tab(0,-1);break;
	case GDK_Down:
	    switch_tab(0,1);break;
	}
    }else{
	switch(event->keyval){
	case GDK_plus:
	case GDK_equal:
	    do_zoom(*drawdatawrap,0,0,1); break;
	case GDK_minus:
	    do_zoom(*drawdatawrap,0,0,-1);break;
	case GDK_0:
	case GDK_1:
	    do_zoom(*drawdatawrap,0,0,0);break;
	case GDK_Left:
	    do_move(*drawdatawrap,-10,0);break;
	case GDK_Right:
	    do_move(*drawdatawrap,10,0);break;
	case GDK_Up:
	    do_move(*drawdatawrap,0,10);break;
	case GDK_Down:
	    do_move(*drawdatawrap,0,-10);break;
	}
    }
    return TRUE;
}


void addpage(drawdata_t **drawdatawrap)
{
    drawdata_t *drawdata=*drawdatawrap;
    GtkWidget *drawarea;
    if(!drawdata->fig) error("Must set fig before calling addpage");
    GSList *roots=NULL;
    int nroot=0;
    for(GSList *p=windows; p; p=p->next){
	GtkWidget *window=p->data;
	GtkWidget *topnb=get_topnb(window);
	GtkWidget *root;
	for(int itab=0; itab<gtk_notebook_get_n_pages(GTK_NOTEBOOK(topnb)); itab++){
	    root=gtk_notebook_get_nth_page(GTK_NOTEBOOK(topnb),itab);
	    GtkWidget *eventbox=gtk_notebook_get_tab_label(GTK_NOTEBOOK(topnb),root);
	    GtkWidget *label=gtk_bin_get_child(GTK_BIN(eventbox));
	    const gchar *text=gtk_label_get_text(GTK_LABEL(label));
	    if(!strcmp(text,drawdata->fig)){
		roots=g_slist_append(roots,root);
		nroot++;//number of roots find.
		break;
	    }else{
		root=NULL;
	    }
	}
    }
    if(!nroot){//root not found.
	GtkWidget *root=gtk_notebook_new();
	roots=g_slist_append(roots,root);
	nroot++;
	gtk_container_set_border_width(GTK_CONTAINER(root),2);
#if GTK_MAJOR_VERSION>=3 || GTK_MINOR_VERSION >= 12
	gtk_notebook_set_group(GTK_NOTEBOOK(root), "secondlevel");
#endif
	gtk_notebook_set_tab_pos(GTK_NOTEBOOK(root),GTK_POS_RIGHT);
	gtk_notebook_set_scrollable(GTK_NOTEBOOK(root), TRUE);
	GtkWidget *label=gtk_label_new(drawdata->fig);
	GtkWidget *eventbox=gtk_event_box_new();
	gtk_container_add(GTK_CONTAINER(eventbox), label);
	gtk_widget_add_events (GTK_WIDGET(eventbox), GDK_BUTTON_PRESS);
	gtk_event_box_set_visible_window(GTK_EVENT_BOX(eventbox),FALSE);
	g_signal_connect(eventbox, "button-press-event", G_CALLBACK(tab_button_cb), root);
	gtk_widget_show_all(eventbox);
	GtkWidget *topnb=get_topnb(GTK_WIDGET(windows->data));//chose the first window.
	gtk_notebook_append_page(GTK_NOTEBOOK(topnb),root,eventbox);
	gtk_notebook_set_tab_detachable(GTK_NOTEBOOK(topnb), root, TRUE);
	gtk_notebook_set_tab_reorderable(GTK_NOTEBOOK(topnb), root, TRUE);
    }
    GtkWidget *page=NULL;
    for(GSList *p=roots; p; p=p->next){//scan through all the root pages with same label
	GtkWidget *root=p->data;
	for(int itab=0; itab<gtk_notebook_get_n_pages(GTK_NOTEBOOK(root)); itab++){
	    page=gtk_notebook_get_nth_page(GTK_NOTEBOOK(root),itab);
	    GtkWidget *label=gtk_notebook_get_tab_label(GTK_NOTEBOOK(root),page);
	    GList *list=gtk_container_get_children(GTK_CONTAINER(label));
	    const gchar *labeltext=gtk_label_get_text(GTK_LABEL(g_list_first(list)->data));
	    if(labeltext && !strcmp(labeltext,drawdata->name)){
		break;
	    }else{
		page=NULL;
	    }
	}
    }
    if(page){
	/*
	  we use drawdatawrap so that we don't have to modify the data on the g_object.
	 */
	drawdata_t **drawdata_wrapold=g_object_get_data(G_OBJECT(page),"drawdatawrap");
	drawdata_t *drawdata_old=(*drawdata_wrapold);
	drawdata->drawarea=drawdata_old->drawarea;
	drawdata->page=drawdata_old->page;
	drawdata->zoomx=drawdata_old->zoomx;
	drawdata->zoomy=drawdata_old->zoomy;
	drawdata->offx=drawdata_old->offx;
	drawdata->offy=drawdata_old->offy;
	drawdata->grid=drawdata_old->grid;
	drawdata->width=drawdata_old->width;
	drawdata->height=drawdata_old->height;
	drawdata->square=drawdata_old->square;
	drawdata_free(*drawdata_wrapold);
	*drawdata_wrapold=drawdata;//just replace the data
	if(get_current_page()==drawdata){//we are the current page. need to update pixmap
	    update_pixmap(drawdata);
	}//otherwise, don't have to do anything.
    }else{
	//new tab inside the fig to contain the plot.
	drawdata->page=page=gtk_scrolled_window_new(NULL,NULL);
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(page), 
				       GTK_POLICY_AUTOMATIC,
				       GTK_POLICY_AUTOMATIC);
	
	g_object_set_data(G_OBJECT(page),"drawdatawrap",drawdatawrap);
	drawdata->drawarea=drawarea=gtk_drawing_area_new();
	gtk_widget_add_events (drawarea,GDK_BUTTON_PRESS_MASK|
			       GDK_BUTTON_RELEASE_MASK|
			       GDK_POINTER_MOTION_MASK|GDK_POINTER_MOTION_HINT_MASK|
			       GDK_BUTTON1_MOTION_MASK|
			       GDK_KEY_PRESS_MASK|
			       GDK_KEY_RELEASE_MASK);
	GTK_WIDGET_SET_FLAGS(drawarea,GTK_CAN_FOCUS);
	GTK_WIDGET_SET_FLAGS(drawarea,GTK_SENSITIVE);
	g_signal_connect(drawarea,"motion-notify-event",
			 G_CALLBACK(motion_notify),drawdatawrap);
	g_signal_connect(drawarea,"button-press-event",
			 G_CALLBACK(button_press),drawdatawrap);
	g_signal_connect(drawarea,"button-release-event",
			 G_CALLBACK(button_release),drawdatawrap);
	g_signal_connect(drawarea,"scroll-event",
			 G_CALLBACK(scroll_event), drawdatawrap);
	/*g_signal_connect(drawarea,"button-release-event",
	  G_CALLBACK(button_release),drawdatawrap);*/
	g_signal_connect(drawarea,"key-press-event",
			 G_CALLBACK(key_press),drawdatawrap);
	gtk_scrolled_window_add_with_viewport
	    (GTK_SCROLLED_WINDOW(page), drawarea);
	GtkWidget *button=tab_label_new(drawdatawrap);
	GtkWidget *root=roots->data;//choose first.
	gtk_notebook_append_page(GTK_NOTEBOOK(root), page,button);
	gtk_notebook_set_tab_detachable(GTK_NOTEBOOK(root), page, TRUE);
	gtk_notebook_set_tab_reorderable(GTK_NOTEBOOK(root), page, TRUE);
	g_signal_connect (drawarea, "expose-event", 
	     G_CALLBACK (on_expose_event), drawdatawrap);
	//handles zooming.
	g_signal_connect (drawarea, "configure-event", 
	     G_CALLBACK (on_configure_event), drawdatawrap);
	gtk_widget_set_size_request(drawarea, DRAWAREA_MIN_WIDTH, 
				    DRAWAREA_MIN_HEIGHT);
	gtk_widget_show_all(root);
    }
    //gdk_threads_leave();
}



static void tool_save(GtkToolButton *button){
    (void)button;
    drawdata_t *drawdata=get_current_page();
    static gchar *folder=NULL;
    char *filename;
    GtkWidget *dialog;
    cairo_surface_t *surface;
    int width,height;
 retry:
    dialog=gtk_file_chooser_dialog_new
	("Select a file to save", 
	 GTK_WINDOW(curwindow),
	 GTK_FILE_CHOOSER_ACTION_SAVE,
	 GTK_STOCK_CANCEL,GTK_RESPONSE_CANCEL,
	 GTK_STOCK_OPEN,GTK_RESPONSE_ACCEPT,
	 NULL);
    gtk_file_chooser_set_do_overwrite_confirmation 
	(GTK_FILE_CHOOSER (dialog), TRUE);
    if(folder){
	gtk_file_chooser_set_current_folder
	    (GTK_FILE_CHOOSER(dialog), folder);
    }else{
	gtk_file_chooser_set_current_folder
	    (GTK_FILE_CHOOSER(dialog), HOME);
    }

    if(gtk_dialog_run(GTK_DIALOG(dialog))==GTK_RESPONSE_ACCEPT){
	g_free(folder);
	folder=gtk_file_chooser_get_current_folder(GTK_FILE_CHOOSER(dialog));
	filename=gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
    }else{
	filename=NULL;
    }
    gtk_widget_destroy(dialog);
    if(!filename) return;

    char *suffix=rindex(filename,'.');
    if(!suffix){
	char *filename2=stradd(filename,".png",NULL);
	g_free(filename);
	filename=filename2;
	suffix=rindex(filename,'.');
    }

    if(strcmp(suffix,".eps")==0){
#ifdef CAIRO_HAS_PS_SURFACE
	width=72*8;
	height=(drawdata)->height*72*8/(drawdata)->width;//same aspect ratio as widget
	surface=cairo_ps_surface_create(filename, width,height);
	cairo_ps_surface_set_eps (surface, TRUE);
#else
	error_msg("ps surface is unavailable");
	goto retry;
#endif
    }else if(strcmp(suffix,".png")==0){
	width=(drawdata)->width;//same size as the widget.
	height=(drawdata)->height;
	surface=cairo_image_surface_create
	    ((cairo_format_t)CAIRO_FORMAT_RGB24,width,height);
    }else{
	error_msg("%s has unknown suffix\n",filename);
	goto retry;
    }
    
    cairo_draw(cairo_create(surface), drawdata,width,height);
    if(strcmp(suffix,".png")==0){
	cairo_surface_write_to_png(surface,filename);
    } 
    cairo_surface_finish(surface);
    cairo_surface_destroy(surface);
 
    g_free(filename);
}


static void tool_zoom(GtkToolButton *button, gpointer data){
    (void)button;
    drawdata_t *drawdata=get_current_page();
    int mode=GPOINTER_TO_INT(data);
    do_zoom(drawdata,0,0,mode);
}
static void tool_toggled(GtkToggleToolButton *button, gpointer data){
    int active=gtk_toggle_tool_button_get_active(button);
    if(active){
	cursor_type=(GPOINTER_TO_INT(data));
    }
}


static void limit_change(GtkSpinButton *button, gpointer data){
    (void)button;
    spin_t *spin=data;
    drawdata_t *drawdata=spin->data;
    //update the values
    drawdata->limit0[spin->i]=gtk_spin_button_get_value(GTK_SPIN_BUTTON(spin->w));
    drawdata->limit_changed=1;
    delayed_update_pixmap(drawdata);
}
static void checkbtn_toggle_square(GtkToggleButton *btn, drawdata_t *drawdata){
    drawdata->square=gtk_toggle_button_get_active(btn);
    delayed_update_pixmap(drawdata);
}
static void checkbtn_toggle_grid(GtkToggleButton *btn, drawdata_t *drawdata){
    drawdata->grid=gtk_toggle_button_get_active(btn);
    delayed_update_pixmap(drawdata);
}
static void checkbtn_toggle_ticinside(GtkToggleButton *btn, drawdata_t *drawdata){
    drawdata->ticinside=gtk_toggle_button_get_active(btn);
    delayed_update_pixmap(drawdata);
}
/**
   Response to the quest to set the zaxis limit (the range of the color bar)
*/
static void tool_property(GtkToolButton *button, gpointer data){
    (void)button;
    (void)data;
    drawdata_t *drawdata=get_current_page();
    if(!drawdata){
	return;
    }
    GtkWidget *dialog=gtk_dialog_new_with_buttons("Figure Properties", 
						  GTK_WINDOW(curwindow), 
						  (GtkDialogFlags)(GTK_DIALOG_DESTROY_WITH_PARENT),
						  GTK_STOCK_OK,
						  GTK_RESPONSE_ACCEPT,
						  GTK_STOCK_CANCEL,
						  GTK_RESPONSE_REJECT,
						  NULL);
    GtkWidget *content_area=gtk_dialog_get_content_area(GTK_DIALOG(dialog));
    GtkWidget *table=gtk_table_new(5,4,0);
   
    int n=4;
    spin_t lim[n];
    for(int i=0; i<n; i++){
	double val=drawdata->limit0[i];
	double step=pow(10,floor(log10(fabs(val)))-2);
	if(fabs(val)<1e-20){
	    step=1;
	}
	if(i<2){
	    lim[i].w=gtk_spin_button_new_with_range(drawdata->limit[0],drawdata->limit[1],step);
	}else{
	    lim[i].w=gtk_spin_button_new_with_range(drawdata->limit[2],drawdata->limit[3],step);
	}
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(lim[i].w), val);
	lim[i].val=&drawdata->limit0[i];
	g_signal_connect(lim[i].w, "value-changed", G_CALLBACK(limit_change), &lim[i]);
	lim[i].i=i;
	lim[i].data=drawdata;
    }
    drawdata->spin=lim;
    GtkWidget *checkbtn;
    gint irow=0;
    checkbtn=gtk_check_button_new_with_label("Make image square");
    gtk_table_attach_defaults(GTK_TABLE(table), checkbtn, 0, 4, irow, irow+1);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbtn), drawdata->square);
    g_signal_connect(checkbtn, "toggled", G_CALLBACK(checkbtn_toggle_square), drawdata);

    irow++;
    checkbtn=gtk_check_button_new_with_label("Enable grids");
    gtk_table_attach_defaults(GTK_TABLE(table), checkbtn, 0, 4, irow, irow+1);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbtn), drawdata->grid);
    g_signal_connect(checkbtn, "toggled", G_CALLBACK(checkbtn_toggle_grid), drawdata);
  
    irow++;
    checkbtn=gtk_check_button_new_with_label("Put tic inside");
    gtk_table_attach_defaults(GTK_TABLE(table), checkbtn, 0, 4, irow, irow+1);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbtn), drawdata->ticinside);
    g_signal_connect(checkbtn, "toggled", G_CALLBACK(checkbtn_toggle_ticinside), drawdata);

    irow++;
    gtk_table_attach_defaults(GTK_TABLE(table), gtk_label_new("xmin"), 0, 1, irow, irow+1);
    gtk_table_attach_defaults(GTK_TABLE(table), lim[0].w, 1, 2, irow, irow+1);
    gtk_table_attach_defaults(GTK_TABLE(table), gtk_label_new("xmax"), 2, 3, irow, irow+1);
    gtk_table_attach_defaults(GTK_TABLE(table), lim[1].w, 3, 4, irow, irow+1);

    irow++;
    gtk_table_attach_defaults(GTK_TABLE(table), gtk_label_new("ymin"), 0, 1, irow, irow+1);
    gtk_table_attach_defaults(GTK_TABLE(table), lim[2].w, 1, 2, irow, irow+1);
    gtk_table_attach_defaults(GTK_TABLE(table), gtk_label_new("ymax"), 2, 3, irow, irow+1);
    gtk_table_attach_defaults(GTK_TABLE(table), lim[3].w, 3, 4, irow, irow+1);

    gtk_container_add(GTK_CONTAINER(content_area), table);
    gtk_widget_show_all(table);
    gint result=gtk_dialog_run(GTK_DIALOG(dialog));
    
    switch (result){
    case GTK_RESPONSE_CANCEL://revert all the changes?
	break;
    default:
	break;
    }
    drawdata->spin=NULL;
    gtk_widget_destroy(dialog);
}


static void tool_font_set(GtkFontButton *btn){
    const char *font_name_new=gtk_font_button_get_font_name(btn);
    PangoFontDescription *pfd
	=pango_font_description_from_string(font_name_new);
    font_name_version++;//tell every expose event to update the figure.
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
    }else{
	font_weight=CAIRO_FONT_WEIGHT_BOLD;
    }
    int size=pango_font_description_get_size(pfd);
    //get font_size in device unit (dots, pixels);
    if(pango_font_description_get_size_is_absolute(pfd)){
	font_size=(double)size/(double)PANGO_SCALE;
    }else{
	gdouble dpi=gdk_screen_get_resolution(gtk_widget_get_screen(curwindow));
	font_size=(double)size/(double)PANGO_SCALE*(double)dpi/72.;
    }
    SP_XL=font_size*2.4+8;
    SP_YT=font_size*1.3+12;
    SP_YB=font_size*2.4+12;
    SP_XR=SP_LEG+LEN_LEG+font_size*2;
    pango_font_description_free(pfd);
}

static void close_window(GtkObject *object){
    windows = g_slist_remove(windows, object);
    GtkWidget *topnb=get_topnb(GTK_WIDGET(object));
    GtkWidget *toolbar=get_toolbar(GTK_WIDGET(object));
    g_signal_handlers_disconnect_by_func(topnb, topnb_page_added, toolbar);
    int npages=gtk_notebook_get_n_pages(GTK_NOTEBOOK(topnb));
    for(int ipage=0; ipage<npages; ipage++){
	GtkWidget *page=gtk_notebook_get_nth_page(GTK_NOTEBOOK(topnb),ipage);
	int ntabs=gtk_notebook_get_n_pages(GTK_NOTEBOOK(page));
	for(int itab=0; itab<ntabs; itab++){
	    GtkWidget *tab=gtk_notebook_get_nth_page(GTK_NOTEBOOK(page),itab);
	    drawdata_t **drawdatawrap=g_object_get_data(G_OBJECT(tab),"drawdatawrap");
	    gtk_widget_hide(tab);
	    drawdata_free(*drawdatawrap); 
	}
    }
    if(!windows){
	gtk_main_quit();
    }
    //don't call exit here.
}

static gboolean window_state(GtkWidget *window, GdkEvent *event){
    if(event->type==GDK_FOCUS_CHANGE){
	GdkEventFocus *focus=(GdkEventFocus*)event;
	if(focus->in){
	    set_cur_window(window);
	}
    }
    return FALSE;
}
GtkWidget *create_window(void){
    if(!pix_hand){
	GdkDisplay *display=gdk_display_get_default();
	pix_hand=gdk_pixbuf_new_from_inline(-1, mouse_hand, FALSE, NULL);
	pix_arrow=gdk_pixbuf_new_from_inline(-1, mouse_white, FALSE, NULL);
	cursors[0]=gdk_cursor_new_from_pixbuf(display, pix_hand, 8, 5);
	cursors[1]=gdk_cursor_new_from_pixbuf(display, pix_arrow, 3, 0);
	gtk_rc_parse_string(rc_string_notebook);
    }
    if(!topmenu){
	topmenu = gtk_menu_new();
	GtkWidget *item=gtk_menu_item_new_with_label("Detach");
	g_signal_connect(item, "activate", G_CALLBACK(topnb_detach), NULL);
	gtk_menu_shell_append(GTK_MENU_SHELL(topmenu), item);
	gtk_widget_show_all(topmenu);
	g_object_ref(topmenu);
    }
    GdkPixbuf *icon_main=gdk_pixbuf_new_from_inline(-1,icon_draw,FALSE,NULL);
    GtkWidget *window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    curwindow=window;
    windows= g_slist_append(windows, window);
    gtk_window_set_icon(GTK_WINDOW(window),icon_main);
    g_object_unref(icon_main);
    
    char title[80];
    snprintf(title,80,"MAOS Draw %d(%d)", fifopid, iwindow);
    gtk_window_set_title(GTK_WINDOW(window), title);

    GtkWidget *toolbar=gtk_toolbar_new();
    gtk_toolbar_set_icon_size(GTK_TOOLBAR(toolbar),GTK_ICON_SIZE_MENU);
    gtk_toolbar_set_style(GTK_TOOLBAR(toolbar), GTK_TOOLBAR_ICONS);

    GtkToolItem *item=gtk_tool_button_new_from_stock(GTK_STOCK_SAVE_AS);
    gtk_toolbar_insert(GTK_TOOLBAR(toolbar),item,-1);
    g_signal_connect(item,"clicked",G_CALLBACK(tool_save),NULL);
    item=gtk_separator_tool_item_new();
    gtk_toolbar_insert(GTK_TOOLBAR(toolbar),item,-1);

    item=gtk_radio_tool_button_new(NULL);
    gtk_tool_button_set_icon_widget(GTK_TOOL_BUTTON(item), gtk_image_new_from_pixbuf(pix_hand));
    g_signal_connect(item,"toggled",G_CALLBACK(tool_toggled),GINT_TO_POINTER(0));
    gtk_toolbar_insert(GTK_TOOLBAR(toolbar),item,-1);

    item=gtk_radio_tool_button_new_from_widget(GTK_RADIO_TOOL_BUTTON(item));
    gtk_tool_button_set_icon_widget(GTK_TOOL_BUTTON(item), gtk_image_new_from_pixbuf(pix_arrow));
    g_signal_connect(item,"toggled",G_CALLBACK(tool_toggled),GINT_TO_POINTER(1));
    gtk_toolbar_insert(GTK_TOOLBAR(toolbar),item,-1);
	
	
    item=gtk_separator_tool_item_new();
    gtk_toolbar_insert(GTK_TOOLBAR(toolbar),item,-1);
    item=gtk_tool_button_new_from_stock(GTK_STOCK_ZOOM_IN);
    g_signal_connect(item,"clicked",G_CALLBACK(tool_zoom),GINT_TO_POINTER(1));
    gtk_toolbar_insert(GTK_TOOLBAR(toolbar),item,-1);
    item=gtk_tool_button_new_from_stock(GTK_STOCK_ZOOM_FIT);
    g_signal_connect(item,"clicked",G_CALLBACK(tool_zoom),GINT_TO_POINTER(0));
    gtk_toolbar_insert(GTK_TOOLBAR(toolbar),item,-1);
    item=gtk_tool_button_new_from_stock(GTK_STOCK_ZOOM_OUT);
    g_signal_connect(item,"clicked",G_CALLBACK(tool_zoom),GINT_TO_POINTER(-1));
    gtk_toolbar_insert(GTK_TOOLBAR(toolbar),item,-1);

    item=gtk_separator_tool_item_new();
    gtk_toolbar_insert(GTK_TOOLBAR(toolbar),item,-1);

    item=gtk_tool_button_new_from_stock(GTK_STOCK_PROPERTIES);
    g_signal_connect(item,"clicked",G_CALLBACK(tool_property),NULL);
    gtk_toolbar_insert(GTK_TOOLBAR(toolbar),item,-1);

    item=gtk_separator_tool_item_new();
    gtk_toolbar_insert(GTK_TOOLBAR(toolbar),item,-1);
    gtk_widget_set_sensitive(toolbar, FALSE);
    item=gtk_tool_item_new();
    GtkWidget *fontsel=gtk_font_button_new_with_font("sans 9");
    gtk_container_add(GTK_CONTAINER(item),fontsel);
    g_signal_connect(GTK_FONT_BUTTON(fontsel),"font-set", 
		     G_CALLBACK(tool_font_set),NULL);
    gtk_toolbar_insert(GTK_TOOLBAR(toolbar),item,-1);
    GtkWidget *topnb=gtk_notebook_new();
    gtk_container_set_border_width(GTK_CONTAINER(topnb),2);
#if GTK_MAJOR_VERSION>=3 || GTK_MINOR_VERSION >= 12
    gtk_notebook_set_group(GTK_NOTEBOOK(topnb), "toplevel");
#endif
    gtk_notebook_set_scrollable(GTK_NOTEBOOK(topnb), TRUE);
    g_signal_connect(GTK_NOTEBOOK(topnb), "page-added", G_CALLBACK(topnb_page_added), toolbar);
    g_signal_connect(GTK_NOTEBOOK(topnb), "page-removed", G_CALLBACK(topnb_page_removed), toolbar);
    GtkWidget *vbox=gtk_vbox_new(FALSE,0);
    gtk_box_pack_start(GTK_BOX(vbox),toolbar,FALSE,FALSE,0);
    gtk_box_pack_start(GTK_BOX(vbox),topnb,TRUE,TRUE,0);
    gtk_container_add(GTK_CONTAINER(window), vbox);
    g_signal_connect(window, "destroy", G_CALLBACK (close_window), NULL);
    gtk_window_set_position(GTK_WINDOW(window), 
			    GTK_WIN_POS_CENTER);
    gtk_window_set_default_size(GTK_WINDOW(window), 600, 500);
    gtk_widget_add_events(GTK_WIDGET(window),GDK_FOCUS_CHANGE_MASK);//in case this is not on.
    g_signal_connect(GTK_WIDGET(window),"event", G_CALLBACK(window_state), NULL);
    gtk_widget_show_all(window); 
    gtk_window_present(GTK_WINDOW(window));
    tool_font_set(GTK_FONT_BUTTON(fontsel));//initialize.
    iwindow++;
    return window;
}

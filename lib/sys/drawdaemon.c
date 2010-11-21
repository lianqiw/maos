/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <fcntl.h>
#include <cairo.h>
#include <cairo/cairo-ps.h>
#include <cairo/cairo-svg.h>
#include <pango/pango.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <glib.h>
#include <gtk/gtk.h>
#include <glib/gprintf.h>
#include <gdk/gdkkeysyms.h>
#include "drawdaemon.h"
#include "misc.h"
#include "proc.h"
#include "icon-draw.h"
char *defname="Title";
char *deffig="Figure";
typedef struct drawdata_t{
    //The following are for surfaces
    int nx, ny;   //array size
    cairo_format_t format;
    double *p0;      //original pointer of double
    unsigned char *p;//converted pointer of char or int.
    int gray;       //do we draw in gray scale or in colored
    //The following are for points
    double *ptsx;//points x
    double *ptsy;//points y
    unsigned int npts;//number of points
    int32_t *style;
    unsigned int nstyle;
    double (*cir)[4];
    double *maxmin;
    unsigned int ncir;
    char *name;
    char *title;
    char *xlabel;
    char *ylabel;
    double *limit;//x,y,limit
    cairo_surface_t *image;
    int handlerid1;
    int handlerid2;
    char *fig;
    GtkWidget *page;
    GtkWidget *drawarea;
    GdkPixmap *pixmap;//in server side memory
    int pending;
    int width;
    int height;
    double zoom,offx,offy;
    double mxdown,mydown;
    double scale;
    double centerx, centery;
    int reconfig;
    int font_name_version;
}drawdata_t;
#define DRAWAREA_MIN_WIDTH 440
#define DRAWAREA_MIN_HEIGHT 320
static const char *fifo=NULL;
static GtkWidget *window=NULL;
static GtkWidget *notebook=NULL;
static int font_name_version=0;
static int no_longer_listen=0;
static char *font_name=NULL;
static double font_size=9;
static cairo_font_slant_t font_style=CAIRO_FONT_SLANT_NORMAL;
static cairo_font_weight_t font_weight=CAIRO_FONT_WEIGHT_NORMAL;
//Spaces reserved for title, label, etc
#define SP_LEG 20//space between image and legend
#define LEN_LEG 25 //size of legend
double SP_XL;//space reserved for ylabel
double SP_YT;//space reserved for title
double SP_YB;//space reserved for xlabel
double SP_XR;//space reserved for legend

static const char *rc_string_notebook={
    "style \"noborder\"{                      \n"
    "GtkNoteBook::draw-border={0 0 0 0}       \n"
    "}                                        \n"
    "class \"GtkNoteBook\" style \"noborder\" \n"
};
#define error_msg(A...) {				\
	GtkWidget *dialog0=gtk_message_dialog_new	\
	    (GTK_WINDOW(window),			\
	     GTK_DIALOG_DESTROY_WITH_PARENT,		\
	     GTK_MESSAGE_ERROR,				\
	     GTK_BUTTONS_CLOSE,				\
	     A);					\
	gtk_dialog_run(GTK_DIALOG(dialog0));		\
	gtk_widget_destroy(dialog0);			\
    }
/*
  static void pango_text(cairo_t *cr, double x, double y, double rot, int center,
  const char *text){
  cairo_save(cr);
  cairo_translate(cr, x, y);
  cairo_set_source_rgb(cr,0,0,0);
  if(fabs(rot)>1.e-10)
  cairo_rotate(cr,rot);
  PangoLayout *layout=pango_cairo_create_layout(cr);
  pango_layout_set_font_description(layout, pfd);
  pango_layout_set_text(layout, text, -1);
  int width, height;
  if(center){
  pango_layout_get_pixel_size(layout,&width,&height);
  warning("width=%d, height=%d\n",width,height);
  //cairo_translate(cr,-width*0.5,-height*0.5);
  pango_cairo_update_layout(cr, layout);
  }

  pango_cairo_show_layout(cr,layout);
  g_object_unref(layout);
  cairo_restore(cr);
  }*/
/**
   convert double to int
*/
static unsigned int crp(double x, double x0){
    double res=1.5-4.*fabs(x-x0);
    if(res>1) res=1.;
    else if(res <0) res=0.;
    return (unsigned int)(res*255.);
}
/**
 convert double to char with color map*/
static void 
dbl2pix(long nx, long ny, int color, const double *restrict p,  void *pout, double *info){
    double max,min;
    if(info[0]>info[1]){
	max=info[0]; min=info[1];
    }else{
	maxmindbl(p, nx*ny, &max, &min);
	info[0]=max; info[1]=min;
    }
    if(color){//colored
	int *pi=pout;
	double scale,offset;
	if(fabs(max-min)>1.e-4*fabs(min)){
	    scale=1./(max-min);
	    offset=0;
	}else{
	    scale=0;
	    offset=0.5;
	}
	for(int i=0; i<nx*ny; i++){
	    double x=(p[i]-min)*scale+offset;
	    pi[i]=crp(x,0.75)<<16 | crp(x, 0.5)<<8 | crp(x, 0.25);
	}
    }else{//b/w
	unsigned char *pc=pout;
	double scale=255./(max-min);
	for(int i=0; i<nx*ny; i++){
	    pc[i]=(unsigned char)((p[i]-min)*scale);
	}
    }
}

/**
   correct floor for round off errors.
*/
static double myfloor(double a){
    double b=floor(a);
    if(a-b>1-1.e-5){
	b++;
    }
    return b;
}
/**
   correct ceil for round off errors;*/

static double myceil(double a){
    double b=ceil(a);
    if(b-a>1-1.e-5){
	b--;
    }
    return b;
}
static void cairo_text_left(cairo_t *cr,double x,double y,const char *text){
    cairo_text_extents_t extents;
    cairo_text_extents (cr, text, &extents);
    double x2 = x - (extents.x_bearing);
    double y2 = y - (extents.height/2 + extents.y_bearing);
    cairo_move_to (cr, x2, y2);
    cairo_show_text (cr, text);
}
static void cairo_text_center(cairo_t *cr,double x,double y,const char *text){
    cairo_text_extents_t extents;
    cairo_text_extents (cr, text, &extents);
    double x2 = (x-(extents.width/2 + extents.x_bearing));
    double y2 = (y-(extents.height/2 + extents.y_bearing));
    cairo_move_to (cr, x2, y2);
    cairo_show_text (cr, text);
}
static void cairo_vltext_center(cairo_t *cr, double x, double y, 
				const char *text){
    cairo_text_extents_t extents;
    cairo_text_extents (cr, text, &extents);
    double x2 = (x-(extents.height/2 + extents.y_bearing));
    double y2 = (y+(extents.width/2 + extents.x_bearing));
    cairo_move_to (cr, x2, y2);
    cairo_rotate(cr,-M_PI/2.);
    cairo_show_text (cr, text);
    cairo_rotate(cr,M_PI/2.);
}
static void cairo_text_powindex(cairo_t *cr, double rot,
				double x, double y, int order){
    if(order==0) return;
    char ticval[80];
    double size=font_size*0.6;//half size of cross
    snprintf(ticval,80,"10");
    cairo_text_extents_t extents;
    cairo_text_extents (cr, ticval, &extents); 
    //draw x
    cairo_save(cr);
    cairo_translate(cr,x,y);
    if(fabs(rot)>0){
	cairo_rotate(cr,rot);
    }
    cairo_move_to(cr,0,0);
    cairo_line_to(cr,size,-size);
    cairo_move_to(cr,0,-size);
    cairo_line_to(cr,size,0);
    cairo_stroke(cr);
 
    cairo_move_to(cr,size,0);
    cairo_show_text(cr,ticval);
    //draw index
    cairo_move_to(cr,size+extents.width+3,-extents.height*2/3);
    snprintf(ticval,80,"%d",order);
    cairo_set_font_size(cr, font_size*0.8);
    cairo_show_text(cr,ticval);
    cairo_restore(cr);
}


static void calc_tic(double *tic1, double *dtic, int *ntic, int *order, 
		     double xmax, double xmin){
    double diff=xmax-xmin;
    int order1;
    if(fabs(diff)<fabs(1.e-4*xmax)){
	order1=(int)floor(log10(fabs(xmax)));
	xmax=xmax/pow(10,order1);
	*order=order1;
	*ntic=2;
	*dtic=0;
	*tic1=xmax;
    }else{
	order1=(int)floor(log10(diff));
	diff/=pow(10,order1);
	if(diff<2){
	    order1=order1-1;
	    diff=diff*10;
	}
	double spacing=0;
	if(diff>=12){
	    spacing=2;
	}else if(diff>=6){
	    spacing=1;
	}else{
	    spacing=0.5;
	}
	xmax/=pow(10,order1); 
	xmin/=pow(10,order1);
	*tic1=myceil(xmin/spacing)*spacing;
	*dtic=spacing;
	*ntic=(int)(myfloor(xmax/spacing)-myceil(xmin/spacing)+1);
	*order=order1;
    }
    if(*order<3 && *order>-2){
	*tic1=*tic1*pow(10,*order);
	*dtic=*dtic*pow(10,*order);
	*order=0;
    }
}
static void cairo_draw(cairo_t *cr, drawdata_t *drawdata, int width, int height){
    //Beginning drawing using cairo.
    //fill white background
    drawdata->font_name_version=font_name_version;
    cairo_rectangle(cr,0,0,width,height);
    cairo_set_source_rgb(cr,1,1,1);
    cairo_fill(cr);
    cairo_surface_t *image=drawdata->image;
    cairo_set_antialias(cr,CAIRO_ANTIALIAS_NONE);
    /*
      cairo_font_options_t *fonto= cairo_font_options_create();
      cairo_font_options_set_hint_metrics(fonto,CAIRO_HINT_METRICS_ON);
      cairo_font_options_set_hint_style (fonto, CAIRO_HINT_STYLE_MEDIUM);
      cairo_font_options_set_antialias (fonto,CAIRO_ANTIALIAS_NONE);
      cairo_set_font_options(cr, fonto);
      cairo_font_options_destroy(fonto);
    */
    int widthim, heightim;
    double xmin=drawdata->limit[0];
    double xmax=drawdata->limit[1];
    double ymin=drawdata->limit[2];
    double ymax=drawdata->limit[3];
    double xmin0,ymin0,xmax0,ymax0;
    xmin0=xmin; ymin0=ymin; xmax0=xmax; ymax0=ymax;
    double scale;
    if(drawdata->image){
	double scale1 = (double)(width-SP_XL-SP_XR)/(double)(drawdata->nx);
	double scale2 = (double)(height-SP_YT-SP_YB)/(double)(drawdata->ny);
	scale  = (scale1<scale2?scale1:scale2);
	widthim=(int)(drawdata->nx*scale);
	heightim=(int)(drawdata->ny*scale);
    }else{
	double scale1=(double)(width-SP_XL-SP_XR)/(xmax-xmin);
	double scale2=(double)(height-SP_YT-SP_YB)/(ymax-ymin);
	scale  = (scale1<scale2?scale1:scale2);
	widthim=(int)((xmax-xmin)*scale);
	heightim=(int)((ymax-ymin)*scale);
    }
    drawdata->scale=scale;
    double xoff=((width-widthim-SP_XL-SP_XR)*0.5)+SP_XL;
    double yoff=((height-heightim-SP_YT-SP_YB)*0.5)+SP_YT;
    cairo_translate(cr,xoff, heightim+yoff);
    drawdata->centerx=xoff+widthim*0.5;
    drawdata->centery=yoff+heightim*0.5;
    //flip upside down so that lower left is (0,0);
    cairo_scale(cr,1,-1);
    double zoom=drawdata->zoom;
    if(drawdata->image){
	cairo_save(cr);
	cairo_rectangle(cr,0,0,widthim,heightim);
	cairo_clip(cr);
	cairo_scale(cr,scale*zoom,scale*zoom);
	double ofx=(drawdata->nx/2)*(1/zoom-1)+drawdata->offx/scale;
	double ofy=(drawdata->ny/2)*(1/zoom-1)+drawdata->offy/scale;
	/*The x and y patterns are negated and then set as
	  translation values in the pattern matrix.*/
	cairo_set_source_surface(cr, image, ofx,ofy);
	if(scale*zoom>1)//use nearest filter for up sampling to get clear images
	    cairo_pattern_set_filter
		(cairo_get_source(cr),CAIRO_FILTER_NEAREST);
	cairo_paint(cr);
	cairo_reset_clip(cr);
	double xdiff=xmax-xmin; 
	double ydiff=ymax-ymin;
	xmin0=xmin-(ofx/drawdata->nx)*xdiff;
	xmax0=xmin0+xdiff/zoom;
	ymin0=ymin-(ofy/drawdata->ny)*ydiff;
	ymax0=ymin0+ydiff/zoom;
	cairo_restore(cr);
    }
 
    if(drawdata->npts>0){
	cairo_save(cr);
	cairo_rectangle(cr,0,0,widthim,heightim);
	cairo_clip(cr);
	cairo_set_antialias(cr,CAIRO_ANTIALIAS_GRAY);
	cairo_set_source_rgba(cr,0.2,0.0,1.0,1.0);
	cairo_set_line_width(cr,1);
	int style=3;
	int size=4;
	double centerx=(xmax+xmin)/2-drawdata->offx/scale;
	double centery=(ymax+ymin)/2-drawdata->offy/scale;
	int ncx=widthim/2;
	int ncy=heightim/2;
	//computed from below ix, iy formula by setting ix, iy to 0 and widthim or heightim
	xmax0=(widthim-ncx)/(zoom*scale)+centerx;
	xmin0=(0-ncx)/(zoom*scale)+centerx;
	ymax0=(heightim-ncy)/(zoom*scale)+centery;
	ymin0=(0-ncy)/(zoom*scale)+centery;

	for(unsigned int ips=0; ips<drawdata->npts; ips++){
	    double ix=(drawdata->ptsx[ips]-centerx)*scale*zoom+ncx;
	    double iy=(drawdata->ptsy[ips]-centery)*scale*zoom+ncy;
	    if(ix<0 ||ix>widthim || iy<0 || iy>heightim) continue;
	    if(drawdata->style){
		style=drawdata->style[ips]&0xF;
		size=(drawdata->style[ips]&0xF0)>>4;
		int color=(drawdata->style[ips]&0xFFFFFF00)>>8;
		int r=color/100;
		int g=(color-r*100)/10;
		int b=color-r*100-g*10;
		cairo_set_source_rgba(cr,r*0.11,g*0.11,b*0.11,1.0);
		if(style>4) style=0;
	    }
	    switch(style){
	    case 0:// .
		cairo_arc(cr, ix, iy, 1, 0, 2*M_PI);
		break;
	    case 1:// o
		cairo_arc(cr, ix, iy, size, 0, 2*M_PI);
		break;
	    case 2:// x
		cairo_move_to(cr,ix-size,iy-size);
		cairo_line_to(cr,ix+size,iy+size);
		cairo_move_to(cr,ix+size,iy-size);
		cairo_line_to(cr,ix-size,iy+size);
		break;
	    case 3:// +
		cairo_move_to(cr,ix-size,iy);
		cairo_line_to(cr,ix+size,iy);
		cairo_move_to(cr,ix,iy-size);
		cairo_line_to(cr,ix,iy+size);
		break;
	    case 4:// square []
		cairo_move_to(cr,ix-size,iy-size);
		cairo_line_to(cr,ix+size,iy-size);
		cairo_line_to(cr,ix+size,iy+size);
		cairo_line_to(cr,ix-size,iy+size);
		cairo_line_to(cr,ix-size,iy-size);
		break;
	    default:
		warning("Invalid style\n");
	    }
	    if(drawdata->style)
		cairo_stroke(cr);
	}
	if(!drawdata->style)
	    cairo_stroke(cr);
	cairo_restore(cr);
    }
    if(drawdata->ncir>0){
	cairo_save(cr);
	cairo_rectangle(cr,0,0,widthim,heightim);
	cairo_clip(cr);
	cairo_set_line_width(cr,1);
	cairo_set_antialias(cr,CAIRO_ANTIALIAS_GRAY);

	double centerx=(xmax+xmin)/2-drawdata->offx/scale;
	double centery=(ymax+ymin)/2-drawdata->offy/scale;
	int ncx=widthim/2;
	int ncy=heightim/2;

	for(unsigned int icir=0; icir<drawdata->ncir; icir++){
	    int color=(int)drawdata->cir[icir][3];
	    int r=color/100;
	    int g=(color-r*100)/10;
	    int b=color-r*100-g*10;
	    double ix=(drawdata->cir[icir][0]-centerx)*scale*zoom+ncx;
	    double iy=(drawdata->cir[icir][1]-centery)*scale*zoom+ncy;
	    cairo_set_source_rgba(cr,r*0.11,g*0.11,b*0.11,1.0);
	    cairo_arc(cr,ix,iy, drawdata->cir[icir][2]*scale*zoom,
		      0,M_PI*2);
	    cairo_stroke(cr);
	}
	cairo_restore(cr);
    }
    //Now doing the border, tic, and colobar
    cairo_identity_matrix(cr);
    cairo_translate(cr, xoff, yoff);
    cairo_rectangle(cr, 0, 0, widthim, heightim);
    cairo_set_source_rgba(cr, 0.0, 0.0, 0.0,1.0);
    cairo_set_line_width(cr,1);
    cairo_stroke(cr);//border
    
    cairo_select_font_face(cr, font_name, font_style, font_weight);
    cairo_set_font_size(cr, font_size);
    cairo_identity_matrix(cr);
    
    char ticval[80];
    double tic1, dtic;
    int ntic, order;
    double sep;
    calc_tic(&tic1,&dtic,&ntic,&order,xmax0,xmin0);
    sep=xmax0-xmin0;
    for(int itic=0; itic<ntic; itic++){
	double ticv=tic1+dtic*itic;
	double val=ticv*pow(10,order);
	double frac=(val-xmin0)/sep;
	cairo_move_to(cr,xoff+widthim*frac,yoff+heightim);
	cairo_line_to(cr,xoff+widthim*frac,yoff+heightim+4);
	cairo_stroke(cr);
	snprintf(ticval,80,"%g",ticv);
	cairo_text_center(cr,xoff+widthim*frac,
			  yoff+heightim+font_size*0.6+6,ticval);
	//pango_text(cr,xoff+widthim*frac,yoff+heightim+12,0,1,ticval);
    }
    cairo_text_powindex(cr,0, xoff+widthim-font_size*2,
			yoff+heightim+6+font_size*2.2,order);
    calc_tic(&tic1,&dtic,&ntic,&order,ymax0,ymin0);
    sep=ymax0-ymin0;
    for(int itic=0; itic<ntic; itic++){
	double ticv=tic1+dtic*itic;
	double val=ticv*pow(10,order);
	double frac=(val-ymin0)/sep;
	cairo_move_to(cr,xoff,yoff+heightim*(1-frac));
	cairo_line_to(cr,xoff-4,yoff+heightim*(1-frac));
	cairo_stroke(cr);
	snprintf(ticval,80,"%g",ticv);
	cairo_vltext_center(cr,xoff-font_size*0.6-4,
			    yoff+heightim*(1-frac),ticval);
    }
    cairo_text_powindex(cr,-M_PI/2,xoff-font_size*1.6,
			yoff+font_size*1.8,order);

    if(drawdata->maxmin){//draw colorbar
	cairo_identity_matrix(cr);
	cairo_translate(cr, xoff+widthim+SP_LEG, yoff);
	cairo_rectangle(cr, 0, 0, LEN_LEG,heightim);
	cairo_pattern_t *bar=cairo_pattern_create_linear(0,0,0,heightim);
	cairo_pattern_add_color_stop_rgb(bar,0,0.5625,0,0);
	cairo_pattern_add_color_stop_rgb(bar,0.1111,1,0,0);
	cairo_pattern_add_color_stop_rgb(bar,0.3651,1,1,0);
	cairo_pattern_add_color_stop_rgb(bar,0.6190,0,1,1);
	cairo_pattern_add_color_stop_rgb(bar,0.8730,0,0,1);
	cairo_pattern_add_color_stop_rgb(bar,1,0,0,0.5);
	cairo_set_source(cr,bar);
	cairo_fill(cr);
	cairo_pattern_destroy(bar);
	cairo_rectangle(cr, 0, 0, LEN_LEG,heightim);
	cairo_set_source_rgba(cr, 0.0, 0.0, 0.0,1.0);
	cairo_set_line_width(cr,1);
	cairo_stroke(cr);

	cairo_set_source_rgba(cr, 0.0, 0.0, 0.0,1.0);
	cairo_set_line_width(cr,1);
	calc_tic(&tic1,&dtic,&ntic,&order,
		 drawdata->maxmin[0],drawdata->maxmin[1]);
	sep=drawdata->maxmin[0]-drawdata->maxmin[1];

	for(int itic=0; itic<ntic; itic++){
	    double ticv=tic1+dtic*itic;
	    double val=ticv*pow(10,order);
	    double frac;
	    if(sep>1.e-10*fabs(drawdata->maxmin[0])){
		frac=(val-drawdata->maxmin[1])/sep;
	    }else{
		if(itic==0) frac=0;
		else if(itic==1) frac=1;
		else frac=-1;
	    }
	    cairo_move_to(cr,0,heightim*(1-frac));
	    cairo_line_to(cr,4,heightim*(1-frac));
	    cairo_move_to(cr,LEN_LEG,heightim*(1-frac));
	    cairo_line_to(cr,LEN_LEG-4,heightim*(1-frac));
	    cairo_stroke(cr);
	    snprintf(ticval,80,"%g",ticv);
	    cairo_text_left(cr,LEN_LEG+4,heightim*(1-frac), ticval);
	}
	cairo_text_powindex(cr,0,LEN_LEG/2,-font_size*0.4-2,order);
    }
    cairo_identity_matrix(cr);
    if(drawdata->title){
	cairo_text_center(cr,xoff+widthim/2,yoff-font_size*0.5-4,drawdata->title);
    }
    if(drawdata->xlabel){
	cairo_text_center(cr,xoff+widthim/2,yoff+heightim+8+font_size*1.8,
			  drawdata->xlabel);
    }
    if(drawdata->ylabel){
	cairo_vltext_center(cr,xoff-font_size*1.8-6, 
			    yoff+heightim/2, drawdata->ylabel);
    }
    
    cairo_destroy(cr);
}

typedef struct updatetimer_t{
    int pending;
    gpointer pdata;
}updatetimer_t;
static void update_pixmap(drawdata_t *drawdata){
    //no more pending updates, do the updating.
    int width=drawdata->width;
    int height=drawdata->height;
    if(drawdata->pixmap){
	g_object_unref(drawdata->pixmap);
    }
    //Create a new server size pixmap and then draw on it.
    drawdata->pixmap=gdk_pixmap_new(window->window, width, height, -1);
    cairo_draw(gdk_cairo_create(drawdata->pixmap), drawdata,width,height);
    gtk_widget_queue_draw(drawdata->drawarea);
}
static gboolean update_pixmap_timer(gpointer timer0){
    updatetimer_t *timer=(updatetimer_t*)timer0;
    gpointer pdata=timer->pdata;
    drawdata_t* drawdata=*((drawdata_t**)pdata);
    if(timer->pending==drawdata->pending){
	update_pixmap(drawdata);
    }
    free(timer);
    return FALSE;
}
/* Create a new backing pixmap of the appropriate size and draw there.*/
static gboolean
on_configure_event(GtkWidget *widget, GdkEventConfigure *event, gpointer pdata){
    (void)event; 
    drawdata_t* drawdata=*((drawdata_t**)pdata);
    drawdata->width=widget->allocation.width;
    drawdata->height=widget->allocation.height;
    drawdata->pending++;
    updatetimer_t *tmp=calloc(1, sizeof(updatetimer_t));
    tmp->pending=drawdata->pending;
    tmp->pdata=pdata;
    if(!drawdata->pixmap){
	update_pixmap(drawdata);
    }else{
	g_timeout_add(100, update_pixmap_timer, tmp);
    }
    return FALSE;
}
/* Redraw the screen from the backing pixmap */
static gboolean
on_expose_event(GtkWidget *widget,
		GdkEventExpose *event,
		gpointer pdata){
    drawdata_t* drawdata=*((drawdata_t**)pdata);
    if(drawdata->font_name_version != font_name_version){
	update_pixmap(drawdata);
    }
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
	free(drawdata->ptsx);
	free(drawdata->ptsy);
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
}
static void do_zoom(drawdata_t *drawdata, int mode){
    if(mode==1){
	drawdata->zoom*=1.2;
    }else if(mode==-1){
	drawdata->zoom*=0.8;
    }else if(mode==0){
	drawdata->zoom=1;
	drawdata->offx=0;
	drawdata->offy=0;
    }
    if(drawdata->zoom<0.0001)
	drawdata->zoom=0.0001;
    else if(drawdata->zoom>10000){
	drawdata->zoom=10000;
    }
    update_pixmap(drawdata);
}
/**
   Delete a figure page.
*/
static void delete_page(GtkButton *btn, drawdata_t **drawdatawrap){
    (void)btn;
    GtkWidget *root;
    root=gtk_notebook_get_nth_page(GTK_NOTEBOOK(notebook), 
				   gtk_notebook_get_current_page
				   (GTK_NOTEBOOK(notebook)));
    
    int ipage=gtk_notebook_page_num(GTK_NOTEBOOK(root), 
				    (*drawdatawrap)->page);
    gtk_notebook_remove_page(GTK_NOTEBOOK(root), ipage);
    drawdata_free(*drawdatawrap);
    free(drawdatawrap);
    if(gtk_notebook_get_n_pages(GTK_NOTEBOOK(root))<=0){
	int jpage=gtk_notebook_page_num(GTK_NOTEBOOK(notebook), root);
	gtk_notebook_remove_page(GTK_NOTEBOOK(notebook),jpage);
    }
}

static GtkWidget *tab_label_new(drawdata_t **drawdatawrap){
    drawdata_t *drawdata=*drawdatawrap;
    const gchar *str=drawdata->name;
    GtkWidget *close_btn;
    GtkWidget *label;
    GtkWidget *out;
    GtkRcStyle *rcstyle;
    close_btn=gtk_button_new();
    gtk_button_set_relief(GTK_BUTTON(close_btn), GTK_RELIEF_NONE);
    gtk_button_set_focus_on_click(GTK_BUTTON(close_btn), FALSE);
    gtk_container_add(GTK_CONTAINER(close_btn), 
		      gtk_image_new_from_stock(GTK_STOCK_CLOSE,
					       GTK_ICON_SIZE_MENU));
    g_signal_connect(close_btn, "clicked", 
		     G_CALLBACK(delete_page), drawdatawrap);
    /* make button as small as possible */
    rcstyle = gtk_rc_style_new();
    rcstyle->xthickness = rcstyle->ythickness = 0;
    gtk_widget_modify_style(close_btn, rcstyle);
    gtk_widget_set_size_request(close_btn,14,14);
    /* create label for tab */
    label = gtk_label_new(str);
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
    gtk_misc_set_padding(GTK_MISC(label), 0, 0);
    gtk_label_set_width_chars(GTK_LABEL(label), 12);
    out=gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(out), label, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(out), close_btn, FALSE, FALSE, 0);
    gtk_rc_style_unref(rcstyle);
    gtk_widget_show_all(out);
    return out;
}
static void do_move(drawdata_t *drawdata, double xdiff, double ydiff){
    double ratio=1./(drawdata->zoom);
    drawdata->offx+=xdiff*ratio;
    drawdata->offy+=ydiff*ratio;
    update_pixmap(drawdata);
}
static gboolean motion_notify(GtkWidget *widget, GdkEventMotion *event, 
			      drawdata_t **drawdatawrap){
    (void)widget;
    drawdata_t *drawdata=*drawdatawrap;
    if(event->state & GDK_BUTTON1_MASK){
	double x, y;
	x = event->x;
	y = event->y;
  
	do_move(drawdata, 
		x - drawdata->mxdown, 
		drawdata->mydown - y);//notice the reverse direction on y
	drawdata->mxdown=x;
	drawdata->mydown=y;
    }
    return FALSE;
}
static gboolean scroll_event(GtkWidget *widget, GdkEventScroll *event, 
			     drawdata_t **drawdatawrap){
    (void)widget;
    drawdata_t *drawdata=*drawdatawrap;
    if(event->direction == GDK_SCROLL_UP){
	do_zoom(drawdata,1);
    }else if(event->direction == GDK_SCROLL_DOWN){
	do_zoom(drawdata,-1);
    }
    return FALSE;
}
static gboolean button_press(GtkWidget *widget, GdkEventButton *event,
			     drawdata_t **drawdatawrap){
    drawdata_t *drawdata=*drawdatawrap;
    //Grab focus so the keys work
    if(!GTK_WIDGET_HAS_FOCUS(widget))
	gtk_widget_grab_focus(widget);

    if(event->type==GDK_2BUTTON_PRESS && event->button==1){//double click brings the part of the image to center.
	drawdata->mxdown=-1;
	drawdata->mydown=-1;
	do_move(drawdata,
		(drawdata->centerx-event->x),
		-(drawdata->centery-event->y));//notice the reverse sign.
    }else{
	drawdata->mxdown=event->x;
	drawdata->mydown=event->y;
    }
    return FALSE;
}
/*
  static gboolean button_release(GtkWidget *widget, GdkEventButton *event,
  drawdata_t **drawdatawrap){
  if(event->button!=1) return FALSE;
  drawdata_t *drawdata=*drawdatawrap;
  (void)widget;
  //remember that mouse counts from lower left.
  double xdiff=event->x-drawdata->mxdown;
  double ydiff=-(event->y-drawdata->mydown);
  if(drawdata->mxdown>0 && drawdata->mydown>0 
  && (fabs(xdiff)>1 || fabs(ydiff)>1)){
  do_move(drawdata,xdiff,ydiff);
  }
  
  return FALSE;
  }*/
static void switch_tab(int lr, int ud){
    if(lr){
	gtk_notebook_set_current_page
	    (GTK_NOTEBOOK(notebook),
	     gtk_notebook_get_current_page(GTK_NOTEBOOK(notebook))+lr);
    }else if(ud){
	GtkWidget *page=gtk_notebook_get_nth_page
	    (GTK_NOTEBOOK(notebook),
	     gtk_notebook_get_current_page(GTK_NOTEBOOK(notebook)));
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
	    do_zoom(*drawdatawrap,1); break;
	case GDK_minus:
	    do_zoom(*drawdatawrap,-1);break;
	case GDK_0:
	case GDK_1:
	    do_zoom(*drawdatawrap,0);break;
	case GDK_Left:
	    do_move(*drawdatawrap,10,0);break;
	case GDK_Right:
	    do_move(*drawdatawrap,-10,0);break;
	case GDK_Up:
	    do_move(*drawdatawrap,0,-10);break;
	case GDK_Down:
	    do_move(*drawdatawrap,0,10);break;
	}
    }
    return TRUE;
}
static void addpage(drawdata_t **drawdatawrap)
{
    drawdata_t *drawdata=*drawdatawrap;
    GtkWidget *drawarea;
    GtkWidget *root=NULL;
    //not invoked by signal so need to acquire lock of gtk
    //gdk_threads_enter();
    if(!drawdata->fig) error("Must set fig before calling addpage");
    for(int itab=0; itab<gtk_notebook_get_n_pages(GTK_NOTEBOOK(notebook)); itab++){
	root=gtk_notebook_get_nth_page(GTK_NOTEBOOK(notebook),itab);
	const gchar *label=gtk_notebook_get_tab_label_text(GTK_NOTEBOOK(notebook),root);
	if(!strcmp(label,drawdata->fig)){
	    break;
	}else{
	    root=NULL;
	}
    }
    if(!root){
	root=gtk_notebook_new();
	gtk_container_set_border_width(GTK_CONTAINER(root),0);
	gtk_notebook_set_tab_pos(GTK_NOTEBOOK(root),GTK_POS_RIGHT);
	gtk_notebook_set_scrollable(GTK_NOTEBOOK(root), TRUE);
	gtk_notebook_set_tab_vborder(GTK_NOTEBOOK(root),0);
	gtk_notebook_set_tab_hborder(GTK_NOTEBOOK(root),2);
	GtkWidget *label=gtk_label_new(drawdata->fig);
	gtk_notebook_append_page(GTK_NOTEBOOK(notebook),root,label);
    }
    GtkWidget *page=NULL;
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
    if(page){
	drawdata_t **drawdata_wrapold=g_object_get_data(G_OBJECT(page),"drawdatawrap");
	drawdata_t *drawdata_old=(*drawdata_wrapold);
	drawdata->drawarea=drawdata_old->drawarea;
	drawdata->page=drawdata_old->page;
	drawdata->zoom=drawdata_old->zoom;
	drawdata->offx=drawdata_old->offx;
	drawdata->offy=drawdata_old->offy;
	drawdata_free(*drawdata_wrapold);
	*drawdata_wrapold=drawdata;//just replace the data
	gtk_widget_hide(drawdata->drawarea);
	gtk_widget_show(drawdata->drawarea);
	//g_signal_emit_by_name(drawdata->drawarea,"configure-event");
	//g_signal_emit_by_name(drawdata->drawarea,"expose-event");
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
			       GDK_BUTTON1_MOTION_MASK|
			       GDK_POINTER_MOTION_HINT_MASK|
			       GDK_KEY_PRESS_MASK|
			       GDK_KEY_RELEASE_MASK);
	GTK_WIDGET_SET_FLAGS(drawarea,GTK_CAN_FOCUS);
	GTK_WIDGET_SET_FLAGS(drawarea,GTK_SENSITIVE);
	g_signal_connect(drawarea,"motion-notify-event",
			 G_CALLBACK(motion_notify),drawdatawrap);
	g_signal_connect(drawarea,"button-press-event",
			 G_CALLBACK(button_press),drawdatawrap);
	g_signal_connect(drawarea,"scroll-event",
			 G_CALLBACK(scroll_event), drawdatawrap);
	/*g_signal_connect(drawarea,"button-release-event",
	  G_CALLBACK(button_release),drawdatawrap);*/
	g_signal_connect(drawarea,"key-press-event",
			 G_CALLBACK(key_press),drawdatawrap);
	gtk_scrolled_window_add_with_viewport
	    (GTK_SCROLLED_WINDOW(page), drawarea);
	GtkWidget *button=tab_label_new(drawdatawrap);
	gtk_notebook_append_page(GTK_NOTEBOOK(root), page,button);

	drawdata->handlerid1= g_signal_connect
	    (drawarea, "expose-event", 
	     G_CALLBACK (on_expose_event), drawdatawrap);
	drawdata->handlerid2= g_signal_connect
	    (drawarea, "configure-event", 
	     G_CALLBACK (on_configure_event), drawdatawrap);
	gtk_widget_set_size_request(drawarea, DRAWAREA_MIN_WIDTH, 
				    DRAWAREA_MIN_HEIGHT);
	gtk_widget_show_all(root);
    }
    //gdk_threads_leave();
}
static GtkWidget *get_current_page(void){
    int p1=gtk_notebook_get_current_page(GTK_NOTEBOOK(notebook));
    GtkWidget *w1=gtk_notebook_get_nth_page(GTK_NOTEBOOK(notebook),p1);//second notebook
    int p2=gtk_notebook_get_current_page(GTK_NOTEBOOK(w1));
    GtkWidget *w2=gtk_notebook_get_nth_page(GTK_NOTEBOOK(w1),p2);//scrolled window
    return w2;
}
static void tool_save(GtkToolButton *button){
    (void)button;
    GtkWidget *page=get_current_page();
    drawdata_t **pdrawdata=g_object_get_data(G_OBJECT(page),"drawdatawrap");
    static gchar *folder=NULL;
    char *filename;
    GtkWidget *dialog;
    cairo_surface_t *surface;
    int width,height;
 retry:
    dialog=gtk_file_chooser_dialog_new
	("Select a file to save", 
	 GTK_WINDOW(window),
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
	height=(*pdrawdata)->height*72*8/(*pdrawdata)->width;//same aspect ratio as widget
	surface=cairo_ps_surface_create(filename, width,height);
	cairo_ps_surface_set_eps (surface, TRUE);
#else
	error_msg("ps surface is unavailable");
	goto retry;
#endif
    }else if(strcmp(suffix,".png")==0){
	width=(*pdrawdata)->width;//same size as the widget.
	height=(*pdrawdata)->height;
	surface=cairo_image_surface_create
	    ((cairo_format_t)CAIRO_FORMAT_RGB24,width,height);
    }else{
	error_msg("%s has unknown suffix\n",filename);
	goto retry;
    }
    
    cairo_draw(cairo_create(surface), *pdrawdata,width,height);
    if(strcmp(suffix,".png")==0){
	cairo_surface_write_to_png(surface,filename);
    } 
    cairo_surface_finish(surface);
    cairo_surface_destroy(surface);
 
    g_free(filename);
}

static void tool_zoom(GtkToolButton *button, gpointer data){
    (void)button;
    GtkWidget *page=get_current_page();
    drawdata_t **pdrawdata=g_object_get_data(G_OBJECT(page),"drawdatawrap");
    int mode=GPOINTER_TO_INT(data);
    do_zoom(*pdrawdata,mode);
}
static void tool_font_set(GtkFontButton *btn){
    const char *font_name_new=gtk_font_button_get_font_name(btn);
    PangoFontDescription *pfd
	=pango_font_description_from_string(font_name_new);
    font_name_version++;
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
	gdouble dpi=gdk_screen_get_resolution(gtk_widget_get_screen(window));
	font_size=(double)size/(double)PANGO_SCALE*(double)dpi/72.;
    }
    SP_XL=font_size*2.4+8;
    SP_YT=font_size*1.3+12;
    SP_YB=font_size*2.4+12;
    SP_XR=SP_LEG+LEN_LEG+font_size*2;
    pango_font_description_free(pfd);
}
static void quitdraw(){
    int npages=gtk_notebook_get_n_pages(GTK_NOTEBOOK(notebook));
    for(int ipage=0; ipage<npages; ipage++){
	GtkWidget *page=gtk_notebook_get_nth_page(GTK_NOTEBOOK(notebook),ipage);
	int ntabs=gtk_notebook_get_n_pages(GTK_NOTEBOOK(page));
	for(int itab=0; itab<ntabs; itab++){
	    GtkWidget *tab=gtk_notebook_get_nth_page(GTK_NOTEBOOK(page),itab);
	    drawdata_t **drawdatawrap=g_object_get_data(G_OBJECT(tab),"drawdatawrap");
	    gtk_widget_hide(tab);
	    drawdata_free(*drawdatawrap); 
	}
    }
    remove(fifo);
    gtk_main_quit();
    exit(0);
}

#define FILE_READ(data,size)			\
    nleft=size;start=(gchar*)data;		\
    do{						\
	int nread=fread(start, 1, nleft, fp);	\
	nleft-=nread;			\
	start+=nread;				\
	if(nread < nleft){			\
	    if(feof(fp)){			\
		info("EOF\n");			\
		return 1;			\
	    }else if(ferror(fp)){		\
		info("File error\n");		\
		return 1;			\
	    }else{				\
		info("Unknown error\n");	\
	    }					\
	}					\
    }while(nleft>0)				

#define FILE_READ_INT(cmd)			\
    FILE_READ(&cmd, sizeof(int));

#define FILE_READ_STR(str)			\
    {						\
	int len;				\
	FILE_READ_INT(len);			\
	str=calloc(len, sizeof(char));		\
	FILE_READ(str,len);			\
    }

static int read_fifo(FILE *fp){
    static drawdata_t *drawdata=NULL;
    int cmd=0;
    gchar *start;
    int nleft;
    static int errcount=0;
    while(1){
	cmd=-1;
	FILE_READ_INT(cmd);
	switch (cmd){
	case FIFO_START:
	    if(drawdata){
		warning("FIFO_START: drawdata is not empty\n");
	    }
	    drawdata=calloc(1, sizeof(drawdata_t));
	    drawdata->zoom=1;
	    drawdata->name=defname;
	    drawdata->format=(cairo_format_t)0;
	    drawdata->gray=0;
	    drawdata->fig=NULL;
	    break;
	case FIFO_DATA:
	    {
		int32_t header[2];
		FILE_READ(header, 2*sizeof(int32_t));
		drawdata->nx=header[0];
		drawdata->ny=header[1];
		int nx=drawdata->nx;
		int ny=drawdata->ny;
		drawdata->p0=malloc(sizeof(double)*nx*ny);
		FILE_READ(drawdata->p0, nx*ny*sizeof(double));
	    }
	    break;
	case FIFO_POINTS:
	    FILE_READ_INT(drawdata->npts);
	    drawdata->ptsx=calloc(drawdata->npts, sizeof(double));
	    drawdata->ptsy=calloc(drawdata->npts, sizeof(double));
	    FILE_READ(drawdata->ptsx, sizeof(double)*drawdata->npts);
	    FILE_READ(drawdata->ptsy,sizeof(double)*drawdata->npts);
	    break;
	case FIFO_STYLE:
	    FILE_READ_INT(drawdata->nstyle);
	    drawdata->style=calloc(drawdata->nstyle, sizeof(int32_t));
	    FILE_READ(drawdata->style, sizeof(int32_t)*drawdata->nstyle);
	    break;
	case FIFO_CIRCLE:
	    FILE_READ_INT(drawdata->ncir);
	    drawdata->cir=calloc(4*drawdata->ncir, sizeof(double));
	    FILE_READ(drawdata->cir,sizeof(double)*4*drawdata->ncir);
	    break;
	case FIFO_LIMIT:
	    drawdata->limit=calloc(4, sizeof(double));
	    FILE_READ(drawdata->limit, 4*sizeof(double));
	    break;
	case FIFO_FIG:
	    FILE_READ_STR(drawdata->fig);
	    break;
	case FIFO_NAME:
	    FILE_READ_STR(drawdata->name);
	    break;
	case FIFO_TITLE:
	    FILE_READ_STR(drawdata->title);
	    break;
	case FIFO_XLABEL:
	    FILE_READ_STR(drawdata->xlabel);
	    break;
	case FIFO_YLABEL:
	    FILE_READ_STR(drawdata->ylabel);
	    break;
	case FIFO_MAXMIN:
	    drawdata->maxmin=calloc(2, sizeof(double));
	    FILE_READ(drawdata->maxmin, sizeof(double)*2);
	    break;
	case FIFO_END:
	    {
		if(drawdata->p0){//draw image
		    int nx=drawdata->nx;
		    int ny=drawdata->ny;
		    size_t size=0;
		    if(nx<=0 || ny<=0) error("Please call _DATA\n");
		    if(drawdata->gray){
			drawdata->format = (cairo_format_t)CAIRO_FORMAT_A8;
			size=1;
		    }else{
			drawdata->format = (cairo_format_t)CAIRO_FORMAT_RGB24;
			size=4;
		    }
		    int stride=cairo_format_stride_for_width(drawdata->format, nx);
		    if(!drawdata->limit){
			drawdata->limit=calloc(4, sizeof(double));
			drawdata->limit[0]=0;
			drawdata->limit[1]=drawdata->nx;
			drawdata->limit[2]=0;
			drawdata->limit[3]=drawdata->ny;
		    }
		    //convert data from double to int/char.
		    if(!drawdata->maxmin){
			drawdata->maxmin=calloc(2, sizeof(double));
		    }
		    drawdata->p=calloc(nx*ny, size);
		    dbl2pix(nx, ny, !drawdata->gray, drawdata->p0, drawdata->p, drawdata->maxmin);
		    gdk_threads_enter();//do I need this?
		    drawdata->image= cairo_image_surface_create_for_data 
			(drawdata->p, drawdata->format, nx, ny, stride);
		    gdk_threads_leave();
		}
		if(drawdata->npts>0){
		    if(!drawdata->limit){
			drawdata->limit=calloc(4, sizeof(double));
			maxmindbl(drawdata->ptsx, drawdata->npts,
				  drawdata->limit+1, drawdata->limit);
			maxmindbl(drawdata->ptsy, drawdata->npts,
				  drawdata->limit+3, drawdata->limit+2);
		    }
		    if(drawdata->nstyle>0 && drawdata->nstyle!=drawdata->npts){
			error("nstyle must equal to npts\n");
		    }
		}
		if(!drawdata->fig) drawdata->fig=deffig;
		drawdata_t **drawdatawrap=calloc(1, sizeof(drawdata_t*));
		drawdatawrap[0]=drawdata;
		gdk_threads_enter();
		addpage(drawdatawrap);
		gdk_threads_leave();
		drawdata=NULL;
	    }
	    break;
	case -1:
	    return 1;//read failed.
	    break;
	default:
	    warning("Unknown cmd: %x\n", cmd);
	    if(errcount++>10){
		no_longer_listen=1;
		return 1;
	    }
	    break;
	}//switch
    }//while
}

static void open_fifo(void){
    FILE *fp=NULL;
    if(no_longer_listen){
	return;
    }
    int retrycount=0;
 retry:
    if(fp) fclose(fp);
    info("Try to open fifo\n");
    fp=fopen(fifo,"rb");
    if(!fp){
	perror("open");
	if(!exist(fifo)){
	    warning("fifo %s is removed\n",fifo);
	    if(mkfifo(fifo,0700)){
		error("Error making fifo\n");
	    }
	}
	if(retrycount++<10){
	    sleep(1);
	    goto retry;
	}
    }
    info("Opened\n");
    read_fifo(fp);
    retrycount=0;
    sleep(1);
    goto retry;
}    

int main(int argc, char *argv[])
{
    if(argc!=2)
	error("Wrong number of arguments\n");
    if(!g_thread_supported()){
	g_thread_init(NULL);
	gdk_threads_init();
    }
    gtk_init(&argc, &argv);

    fifo=argv[1];
    int ppid;
    const char *fifo2=strstr(fifo,"drawdaemon");
    if(!fifo2){
	warning("drawdaemon not found in string\n");
	ppid=getpid();
    }else{
	sscanf(fifo2, "drawdaemon_%d.fifo", &ppid);
    }
    {//record the drawdaemon pid. and redirect output
	char fnpid[PATH_MAX];
	snprintf(fnpid, PATH_MAX,"%s/drawdaemon_%d.pid", TEMP, ppid);
	FILE *fp=fopen(fnpid, "w");
	fprintf(fp, "%d", (int)getpid());
	fclose(fp);
    
	char fnlog[PATH_MAX];
	snprintf(fnlog, PATH_MAX,"%s/drawdaemon_%d.log", TEMP, ppid);
	if(!freopen(fnlog, "w", stdout)) {
	    perror("freopen");
	    warning("Error redirect stdout\n");
	}
	if(!freopen(fnlog, "w", stderr)) {
	    perror("freopen");
	    warning("Error redirect stderr\n");
	}
	setbuf(stdout,NULL);//disable buffering.
	setbuf(stderr,NULL);
    }
    if(!window){
	GdkPixbuf *icon_main=gdk_pixbuf_new_from_inline(-1,icon_draw,FALSE,NULL);
	window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
	gtk_window_set_icon(GTK_WINDOW(window),icon_main);
	gtk_window_set_title(GTK_WINDOW(window),"MAOS Draw");
	GtkWidget *toolbar=gtk_toolbar_new();
	gtk_toolbar_set_icon_size(GTK_TOOLBAR(toolbar),GTK_ICON_SIZE_MENU);
	gtk_toolbar_set_style(GTK_TOOLBAR(toolbar), GTK_TOOLBAR_ICONS);

	GtkToolItem *item=gtk_tool_button_new_from_stock(GTK_STOCK_SAVE_AS);
	gtk_toolbar_insert(GTK_TOOLBAR(toolbar),item,-1);
	g_signal_connect(item,"clicked",G_CALLBACK(tool_save),NULL);
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
	item=gtk_tool_item_new();
	GtkWidget *fontsel=gtk_font_button_new_with_font("sans 9");
	gtk_container_add(GTK_CONTAINER(item),fontsel);
	g_signal_connect(GTK_FONT_BUTTON(fontsel),"font-set", 
			 G_CALLBACK(tool_font_set),NULL);
	gtk_toolbar_insert(GTK_TOOLBAR(toolbar),item,-1);
	notebook=gtk_notebook_new();
	gtk_container_set_border_width(GTK_CONTAINER(notebook),0);
	gtk_notebook_set_scrollable(GTK_NOTEBOOK(notebook), TRUE);
	gtk_notebook_set_tab_vborder(GTK_NOTEBOOK(notebook), 0);
	gtk_notebook_set_tab_hborder(GTK_NOTEBOOK(notebook), 2);
	GtkWidget *vbox=gtk_vbox_new(FALSE,0);
	gtk_box_pack_start(GTK_BOX(vbox),toolbar,FALSE,FALSE,0);
	gtk_box_pack_start(GTK_BOX(vbox),notebook,TRUE,TRUE,0);
	gtk_container_add(GTK_CONTAINER(window), vbox);
	g_signal_connect(window, "destroy",
			 G_CALLBACK (quitdraw), NULL);
	gtk_window_set_position(GTK_WINDOW(window), 
				GTK_WIN_POS_CENTER);
	gtk_window_set_default_size(GTK_WINDOW(window), 600, 500);
	gtk_widget_show_all(window); 
	tool_font_set(GTK_FONT_BUTTON(fontsel));//initialize.
    }
    gtk_rc_parse_string(rc_string_notebook);
    if(NULL==g_thread_create((GThreadFunc)open_fifo, NULL, FALSE, NULL)){
	error("Thread create failed.\n");
    }
    gdk_threads_enter();
    gtk_main();
    gdk_threads_leave();
    {
	char fn[PATH_MAX];
	snprintf(fn, PATH_MAX,"%s/drawdaemon_%d.pid", TEMP, ppid);
	remove(fn);
	snprintf(fn, PATH_MAX,"%s/drawdaemon_%d.log", TEMP, ppid);
	remove(fn);
    }
    remove(fifo);
}//main

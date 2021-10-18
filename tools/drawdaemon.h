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
#ifndef AOS_TOOLS_DRAWDAEMON_H
#define AOS_TOOLS_DRAWDAEMON_H

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <fcntl.h>
#include <cairo.h>
#include <cairo-ps.h>
#include <cairo-pdf.h>
#include <cairo-svg.h>
#include <pango/pango.h>
#include <glib.h>
#include <gtk/gtk.h>
#include <stdlib.h>
#include <stdio.h>
#include "../sys/sys.h"
#include "mygtk.h"
#ifndef CAIRO_FORMAT_A8
#define CAIRO_FORMAT_RGB24 0x01
#define CAIRO_FORMAT_A8 0x02
#endif
typedef struct drawdata_t drawdata_t;
extern int sock;
extern int sock_idle;

struct drawdata_t{
	char* fig;
	char* name;
	struct drawdata_t *next;//form an linked list.
	/*First, input data from draw.c */
	/*Draw images. */
	cairo_surface_t* image;
	float* p0;      /*2d array of data. */
	int nx, ny;   /*array size */
	int nmax;     /*allocated size of array*/
	unsigned char* p;/*converted pointer of char or int. */
	/*Draw points */
	float** pts;      /*pts; */
	int(*ptsdim)[2];  /*nx, ny of pts*/
	int npts;        /*number of pts mat, not points. */
	int nptsmax;     /**<allocated size of pts*/
	/*styles*/
	int32_t* style;
	int* style_pts;    /*save pts style for legend */
	int nstyle;
	int nstylemax;     /*memory storeage of style*/
	/*draw circles */
	float(*cir)[4];
	int ncir;
	int ncirmax; /*storage size of cir*/
	
	/*limit */
	float* limit_data;/*x,y,limit of data. might be supplied by user. */
	float* limit_cumu;/*x,y,limit of cumulatively averaged data. */
	float* limit;/*points to either limit_data or limit_cumu */
	float zlim[4];//2 additional elements for i/o in case double is passed in
	int limit_manual; /*limit_data is supplied by user*/
	char xylog[2];
	//misc
	int byte_float; //record the value used
	int ready;      //ready is set to 0 when data is being read and 1 after wards.
	int recycle;    //when set in GUI thread, data to be deleted by drawdaemon_io() 
	int delete;     //when set in io thread, page will be deleted by addpage()
	float io_time;  //time data was received.
	/*The following are for surfaces */
	
	cairo_format_t format;

	int gray;       /*do we draw in gray scale or in colored */
	
	char* title;
	char* xlabel;
	char* ylabel;
	char** legend;
	
	GtkWidget* page;
	GtkWidget* drawarea;
#if GTK_MAJOR_VERSION>=3 
	cairo_surface_t* pixmap;
#else
	GdkPixmap* pixmap;/*server side memory. */
#endif
	char tooltip[64];

	gint pwidth, pheight;/*size of pixmap.*/
	GtkWidget** spins;/*used on the dialog to change limits. */
	cairo_surface_t* cacheplot;/*cache the plot results so that we don't have to redraw during just panning. */
	int pending;/*drawing is pending. */
	int width;/*width of the canvas */
	int height;/*height of the canvas */

	int widthim;/*width of the part of the canvas for drawing */
	int heightim;/*height of the part of the canvas for drawing */
	int widthim_last, heightim_last;/*width,height of last drawing canvas. */

	float zoomx, zoomy;/*zoom level. */
	float zoomxlast, zoomylast;/*last zoom level. */
	float offx, offy;/*off set of the center of the data. */
	float mxdown, mydown, mtdown;/*mouse pointer down. */
	float dxdown, dydown; /*length of rectangular*/
	int draw_rect; /*draw a rectangular with mxdown, mydown, dxdown, dydown*/
	
	float scalex, scaley;/*scale of the data to fit the display. */
	float centerx, centery;
	float xoff, yoff;/*offset of the area to draw figure. */
	int ncxoff, ncyoff;/*offset of ncx, ncy */
	float limit0[4];/*x,y limit of displayed region. */
	
	int square;/*make x/y scaling be the same, for image and coordinate display */
	int valid;/*move is valid. */
	int font_name_version;
	int grid;/*whether we want grid lines. */
	int ticinside;/*put tick inside. */
	int cursorinside;
	int limit_changed;/*limit has changed. */
	int legendbox;/*draw the legend */
	int legendcurve;/*mark each line with legend entry*/
	float legendoffx;/*location of legend along x. */
	float legendoffy;/*location of legend along y */
	int drawn;/*whether we have been drawn.  */
	int cumu;/*plot cumulative mean. */
	int cumuquad;/*make cumulative quadrature */
	int cumuquadlast;
	/*icumu has to be float because it is used by the GtkSpin */
	float icumu;/*plot cumulative mean from this time step if cumu!=0 */
	float icumulast;/*plot cumulative mean from this time step if cumu!=0 */
	int cumulast;/*=0: we are drawing cumu the first time. */
};
extern float io_time1;/*The time this data is received.*/
extern float io_time2;/*The time last data is received.*/
extern char* font_name;
extern float font_size;
extern cairo_font_slant_t font_style;
extern cairo_font_weight_t font_weight;
extern char* fifo;
extern char* defname;
extern char* deffig;
extern int font_name_version;
extern int ndrawdata;
extern pthread_mutex_t mutex_drawdata;
extern int fifopid;
/*Spaces reserved for title, label, etc */
#define SP_LEG 20/*space between image and legend */
#define LEN_LEG 25 /*size of legend */

extern float SP_XL;/*space on x, left */
extern float SP_XR;/*space on x, right */
extern float SP_YT;/*space on y, top */
extern float SP_YB;/*space on y, buttom */
extern PangoFontDescription* desc;
extern pthread_mutex_t drawdata_mutex;
/*from drawdaemon_draw */
void round_limit(float* xmin, float* xmax, int logscale);
void cairo_draw(cairo_t* cr, drawdata_t* drawdata, int width, int height);
void apply_limit(drawdata_t* drawdata);
/*from drawdaemon_gui */
GtkWidget* create_window(GtkWidget *window);
gboolean addpage(gpointer user_data);
int delete_page(gpointer user_data);
/*from drawdaemon_io */
void* listen_draw(void*);
void flt2pix(long nx, long ny, int color, const float* restrict p, void* pout, float* info);
void fmaxmin(const float* p, long n, float* max, float* min);
void round_limit(float* xmin, float* xmax, int logscale);
#endif

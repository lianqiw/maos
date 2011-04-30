/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#include "aos.h"
#ifndef CAIRO_FORMAT_A8
#define CAIRO_FORMAT_RGB24 0x01
#define CAIRO_FORMAT_A8 0x02
#endif
typedef struct drawdata_t drawdata_t;

struct drawdata_t{
    //First, input data from draw.c
    //Draw images.
    cairo_surface_t *image;
    double *p0;      //original pointer of double
    unsigned char *p;//converted pointer of char or int.
    //Draw points
    dmat **pts;      //pts;
    int npts;        //number of pts mat, not points.
    int32_t *style;
    unsigned int nstyle;
    //draw circles
    double (*cir)[4];
    unsigned int ncir;
    //limit
    double *limit_data;//x,y,limit of data. might be suplied by user.
    double *limit_cumu;//x,y,limit of cumulatively averaged data.
    double *limit;//points to either limit_data or limit_cumu
    double *zlim;
    //The following are for surfaces
    int nx, ny;   //array size
    cairo_format_t format;

    int gray;       //do we draw in gray scale or in colored

    char *name;
    char *title;
    char *xlabel;
    char *ylabel;
    char **legend;
    char *fig;


    GtkWidget *page;
    GtkWidget *drawarea;
    GdkPixmap *pixmap;//server side memory.
    GtkWidget **spins;//used on the dialog to change limits.

    int pending;//drawing is pending.
    int width;//width of the canvas
    int height;//height of the canvas

    int widthim;//width of the part of the canvas for drawing
    int heightim;//height of the part of the canvas for drawing
    int widthim_last, heightim_last;//width,height of last drawing canvas.
    
    double zoomx, zoomy;//zoom level.
    double offx,offy;//off set of the center of the data.
    double mxdown,mydown;//mouse pointer down.
    double scalex, scaley;//scale of the data to fit the display.
    double centerx, centery;
    double xoff, yoff;//offset of the area to draw figure.
    
    double limit0[4];//x,y limit of displayed region.

    int square;//make x/y scaling be the same, for image and coordinate display
    int valid;//move is valid.
    int font_name_version;
    int grid;//whether we want grid lines.
    int ticinside;//put tick inside.
    int cursorinside;
    int limit_changed;//limit has changed.
    int legendbox;
    int drawn;//whether we have been drawn. 
    int cumu;//plot cumulative mean.
    int cumuquad;//make cumulative quadrature
    int icumu;//plot cumulative mean from this time step if cumu!=0
    int icumulast;//plot cumulative mean from this time step if cumu!=0
    int cumulast;//=0: we are drawing cumu the first time.
};
extern char *font_name;
extern double font_size;
extern cairo_font_slant_t font_style;
extern cairo_font_weight_t font_weight;
extern char *fifo;
extern char *defname;
extern char *deffig;
extern int font_name_version;
extern int ndrawdata;
extern pthread_mutex_t mutex_drawdata;
extern int fifopid;
//Spaces reserved for title, label, etc
#define SP_LEG 20//space between image and legend
#define LEN_LEG 25 //size of legend

extern double SP_XL;//space on x, left
extern double SP_XR;//space on x, right
extern double SP_YT;//space on y, top
extern double SP_YB;//space on y, buttom
extern PangoFontDescription *desc;

//from drawdaemon_draw
void round_limit(double *xmin, double *xmax);
void cairo_draw(cairo_t *cr, drawdata_t *drawdata, int width, int height);
void apply_limit(drawdata_t *drawdata);
//from drawdaemon_gui
GtkWidget* create_window(void);
void addpage(drawdata_t **drawdatawrap);
//from drawdaemon_io
void open_fifo(void *);
void dbl2pix(long nx, long ny, int color, const double *restrict p,  void *pout, double *info);
#endif

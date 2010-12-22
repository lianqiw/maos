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
typedef struct {
    GtkWidget *w;
    double *val;
    drawdata_t *data;
    int i;
}spin_t;

struct drawdata_t{
    //The following are for surfaces
    int nx, ny;   //array size
    cairo_format_t format;
    double *p0;      //original pointer of double
    unsigned char *p;//converted pointer of char or int.
    int gray;       //do we draw in gray scale or in colored
    //The following are for points
    dmat **pts;      //pts;
    int npts;        //number of pts mat, not points.
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
    char *fig;
    GtkWidget *page;
    GtkWidget *drawarea;
    GdkPixmap *pixmap;//server side memory.
    int pending;
    int width;//width of the canvas
    int height;//height of the canvas

    int widthim;//width of the part of the canvas for drawing
    int heightim;//height of the part of the canvas for drawing
    int widthim_last, heightim_last;//width,height of last drawing canvas.

    double zoomx, zoomy;//zoom level.
    double offx,offy;//off set of the center of the data.
    double mxdown,mydown;
    double mdx, mdy;
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
    int drawn;//whether we have been drawn. 
    spin_t *spin;//used on the dialog to change limits.
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


//from drawdaemon_draw
void round_limit(double *xmin, double *xmax);
void cairo_draw(cairo_t *cr, drawdata_t *drawdata, int width, int height);
void apply_limit(drawdata_t *drawdata);
//from drawdaemon_gui
GtkWidget* create_window(void);
void addpage(drawdata_t **drawdatawrap);
//from drawdaemon_io
void open_fifo(void);
#endif

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
#ifndef AOS_TOOLS_DRAWDAEMON_H
#define AOS_TOOLS_DRAWDAEMON_H

#include "../sys/sys.h"
#include "mygtk.h"
#ifndef CAIRO_FORMAT_A8
#define CAIRO_FORMAT_RGB24 0x01
#define CAIRO_FORMAT_A8 0x02
#endif
typedef struct drawdata_t drawdata_t;
extern int sock;
extern int sock_idle;
/**
 * Canvas coordinate: what is displayed on the window. if scale is sv.
 * Data coordinate:  in internal data. scale sd=sv/scalex
 * Zoomed data coordinate: zoomed data. scale sz=sv/scalex/zoomx
  */
struct drawdata_t{
	struct drawdata_t *next;//form an linked list.
	char* fig; 			//topnb tab
	char* name;			//subnb tab
	char *title;		//figure title
	char *xlabel;		//x axis label
	char *ylabel;		//y axis label
	char **legend;		//legends
	char **legend_ellipsis;	//legends shortened
	char *filename; 	//previous filename to save
	char *filename_gif; //if set, save all new frames and convert to gif
	
	/*First, input data from draw.c */
	/*For 2-d images. data range is referred to as z axis*/
	float* p0;      	/*2d array of data. */
	float *p1;      	/*2d array of data. used when lpf is < 1. */
	int nx, ny;      	/*p0 array size */
	int nx_last, ny_last; /*Previous p0 array size*/
	int nmax;     		/*allocated size of p0 array*/
	float zlim[4];		//zlim. first two elements are for zlog=0; last two elements are for zlog=1.
	int zlim_manual;	//zlim is manually set.
	int zlim_changed;	//zlim has changed.
	int zlog; 	    	/*draw image in log scale*/
	int zlog_last; 		/*zlog status during previous call to cairo_draw()*/
	int gray;       	/*do we draw in gray scale or in colored */
	void* p;			/*converted from p0 of char or int depends on value of gary */
	cairo_surface_t *image;/*image from p.*/

	/*Draw points */
	float** pts;      	/*set of points to draw; */
	int(*ptsdim)[2];  	/*nx, ny of pts*/
	int npts;        	/*number of pts sets. */
	int nptsmax;    	/*allocated size of pts sets*/

	/*styles for points */
	uint32_t* style;
	int* style_pts;    /*save pts style for legend */
	int nstyle;
	int nstylemax;     /*memory storeage of style*/
	/*draw circles */
	float(*cir)[4];
	int ncir;
	int ncirmax; /*storage size of cir*/

	/*dimension of points (x/y) or array size */
	float* limit_data;	/*x,y,limit of data. might be supplied by user. */
	float* limit_cumu;	/*x,y,limit of cumulatively averaged data. */
	float* limit;		/*x,y limit points to either limit_data or limit_cumu */
	float limit0[4];	/*x,y limit of displayed region. equals to limit if no zoom or pan */
	int limit_manual; 	/*limit_data is manually set and should be not changed*/

	int update_zoom; /*1: update to preserve zoomed area. 2: reset zoom.*/
	int update_limit;/*update limit for pts*/
	//drawy x or y axis with log scale
	char xylog[2];  /*draw in log scale x, y axis*/
	
	//misc
	int byte_float; //number of bytes used for float (4 or 8)
	int ready;      //ready is set to 0 when data is being read and 1 after wards.
	unsigned int recycle;//when set in GUI thread, data to be deleted by drawdaemon_io(). operated by atomic operations
	float io_time;  //time data was received.
	/*The following are for surfaces */

#if GTK_VERSION_AFTER(4,10)
	GFile *file;
#endif	
	GtkWidget* drawarea;
	GtkWidget *subnb;
	cairo_surface_t* pixmap;//for expose
	cairo_surface_t* pixmap2;//for cairo_draw in a separate thread
	cairo_surface_t* cacheplot;//cache the plot results so that we don't have to redraw during just panning.
	cairo_surface_t* surface;//for saving to file
	pthread_mutex_t mutex;
	pthread_t thread;
	char tooltip[64];

	gint pwidth, pheight;/*size of pixmap.*/
	gint pwidth2, pheight2;/*size of pixmap2.*/
	GtkWidget** spins;/*used on the dialog to change limits. */

	int cache_width, cache_height; /*width and height of cacheplot.*/
	int iframe;/*used by delayed_update_pixmap to limit rate of drawing. */
	int width, height;/*width and height of the entire drawing area */

	int widthim, heightim;/*width and height of the part of the canvas for drawing */
	int widthim_last, heightim_last;/*width,height of last drawing canvas. */

	float zoomx, zoomy;		/*zoom level. */
	float zoomx_last, zoomy_last;/*last zoom level. */
	float offx, offy;   	/*offset of the plot in visual coordinate.*/
	int ncxoff, ncyoff;    	/*offset of ncx, ncy */
	//to handle drag and drop
	float mxdown, mydown, mtdown;/*mouse pointer down. */
	float dxdown, dydown; 	/*length of rectangular*/
	int draw_rect; 			/*draw a rectangular with mxdown, mydown, dxdown, dydown*/

	float scalex, scaley;  /**<scale of the data size to fit the drawarea size.  scales with window size*/
	float centerx, centery;/**<center of the drawarea data area. changes with window size.*/
	float xoff, yoff;      /*offset of the data area in the drawarea.*/
	
	int square;/*make x/y scaling be the same, for image and coordinate display */
	int region;/*drag region. */
	int font_name_version;
	int grid;/*whether we want grid lines. */
	int ticinside;/*put tick inside. */
	int cursorinside;

	int legendbox;/*draw the legend */
	int legendcurve;/*mark each line with legend entry*/
	float legendoffx;/*location of legend along x. */
	float legendoffy;/*location of legend along y */
	float legbox_width;/*legend box width*/
	float legbox_height;/*legend box height*/
	float legbox_ox; /*legend box origin in x*/
	float legbox_oy; /*legend box origin in y*/
	float legbox_rx; /*mapping movement to legendoffx*/
	float legbox_ry; /*mapping movement to legendoffy*/
	int drawn;/*whether we have been drawn.  */
	int cumu;/*plot cumulative mean. */
	int cumuquad;/*make cumulative quadrature */
	int cumuquadlast;
	int cumuover;/*overlay cumulative plot on top of regular*/
	/*icumu has to be float because it is used by the GtkSpin */
	float icumu;/*plot cumulative mean from this time step if cumu!=0 */
	float icumulast;/*plot cumulative mean from this time step if cumu!=0 */
	int cumulast;/*=0: we are drawing cumu the first time. */
	//Frame counter
	int frame_io;  /**<Frame counter from input.*/
	int frame_draw;/**<Frame counter that is drawn.*/
	int frame_gif; /**<Frame counter saved for gif.*/
	int session;   /**<Increses when a new client connects */
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
extern int fifopid;
extern float lpf;
extern int noellipsis; 	/*do not allow legend ellipsis.*/
/*Spaces reserved for title, label, etc */
extern int desc_font_size;
extern PangoFontDescription* desc;
extern int client_pid;
extern char *client_hostname;
extern char *client_path_full;
extern char *client_path;
extern char *client_exename;
extern int keep_listen;

extern int hide_xlabel;
extern int hide_ylabel;
extern int hide_title;
extern int hide_legend;
extern int hide_colorbar;

extern GdkPixbuf *icon_main;
extern GdkPixbuf *icon_log;
extern GdkPixbuf *icon_avg;

/*from drawdaemon_draw */
void round_limit(float* xmin, float* xmax, int logscale);
void cairo_draw(drawdata_t* drawdata);
gboolean drawarea_refresh(GtkWidget *drawarea);//called by cairo_draw to redraw drawarea
void update_zoom(drawdata_t* drawdata);
/*from drawdaemon_gui */
GtkWidget* create_window(GtkWidget *window);
gboolean addpage(drawdata_t *drawdata);
int delete_page(drawdata_t *drawdata);
/*from drawdaemon_io */
void* listen_draw(void*);
void flt2pix(const float *restrict p, void *pix, long nx, long ny, int gray, float *zlim, int zlim_manual, int zlog);
void fmaxmin(const float* p, long n, float* max, float* min);
void round_limit(float* xmin, float* xmax, int logscale);
gboolean update_title(gpointer data);
gboolean update_fpslabel(gpointer label);
gboolean finalize_gif(void*);
#endif

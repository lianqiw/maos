/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

<<<<<<< Updated upstream
=======

>>>>>>> Stashed changes
#include <math.h>
#include "drawdaemon.h"
/*
  Routines in this file are about drawing in the cairo surface.

  2011-06-02: New idea to improve the effiency of drawdaemon:
  1) Only draw once, to a cached surface with no zoom, or displacement.
  2) Draw the cached surface to the real surface with zoom or displacement when user requests.
*/
const double stroke_dot[2]={1,5};
//const float stroke_dash[2]={10,10};
//const float stroke_solid[2]={10,0};
int maxtic_x=12;/*Maximum number of tics along x */
int maxtic_y=12;/*Maximum number of tics along y */

float SP_XL;/*space reserved for ylabel */
float SP_YT;/*space reserved for title */
float SP_YB;/*space reserved for xlabel */
float SP_XR;/*space reserved for legend */
#define DRAW_NEW 1 //1: draw to a cache surface

int default_color_table[]={0x0000FF,
			   0xFF0000,
			   0x00FF00,
			   0x009999,
			   0x00FFFF,
			   0x9900CC,
			   0xFFCC00,
			   0xFF00FF,
			   0x000000,
			   0x666666,
};
#define default_color(i) default_color_table[i%11]
#define set_color(cr,color)				   \
    cairo_set_source_rgba(cr,				   \
			  ((color>>16)&0xFF)/255.,	   \
			  ((color>>8)&0xFF)/255.,	   \
			  ((color)&0xFF)/255.,1)

/**
   correct floor for round off errors.
*/
static float myfloor(float a){
	float b=floor(a);
	if(a-b>1-1.e-5){
		b++;
	}
	return b;
}
/**
   correct ceil for round off errors;*/

static float myceil(float a){
	float b=ceil(a);
	if(b-a>1-1.e-5){
		b--;
	}
	return b;
}
/**
 * determine text size
 * */
static void pango_size(PangoLayout* layout, const char* text, float* legwidth, float* legheight){
	if(!text){
		warning("text is null, canceled\n");
		*legwidth=0;
		*legheight=0;
		print_backtrace();
	}	
	pango_layout_set_markup(layout, text, -1);
	int width, height;
	pango_layout_get_size(layout, &width, &height);
	*legwidth=(float)width/PANGO_SCALE;
	*legheight=(float)height/PANGO_SCALE;
}
/**
   Drawing text. frac=0: left align. frac=0.5: center. frac=1: right align
*/
static void pango_text(cairo_t* cr, PangoLayout* layout, float x, float y,
	const char* text, float fracx, float fracy, int vertical){
	if(!text){
		warning("text is null, canceled\n");
		print_backtrace();
		return;
	}
	cairo_save(cr);
	pango_layout_set_markup(layout, text, -1);
	if(fracx>0||fracy>0){
		int width, height;
		pango_layout_get_size(layout, &width, &height);
		if(vertical){
			cairo_move_to(cr, x-height*fracx/PANGO_SCALE, y+width*fracy/PANGO_SCALE);
		} else{
			cairo_move_to(cr, x-width*fracx/PANGO_SCALE, y-height*fracy/PANGO_SCALE);
		}
	} else{
		cairo_move_to(cr, x, y);
	}
	if(vertical){
		cairo_rotate(cr, -M_PI/2);
	}
	pango_cairo_update_layout(cr, layout);
	pango_cairo_show_layout(cr, layout);
	cairo_restore(cr);
}

static void pango_text_powindex(cairo_t* cr, PangoLayout* layout, float x, float y, int order, int vertical){
	if(order==0) return;
	char powindex[40];
	snprintf(powindex, 40, "10<sup>%d</sup>", order);
	pango_text(cr, layout, x, y, powindex, 0, 0, vertical);
}


static void calc_tic(float* tic1, float* dtic, int* ntic, int* order,
	float xmax, float xmin, int maxtic, int logscale){
//when logscale is true, the xmin, xmax are already log of data.
//dbg("xmin=%g, xmax=%g, \n", xmin, xmax);
	(void)maxtic;
	float diff=xmax-xmin;
	/*first get the order of magnitude. */
	float rmax=MAX(fabs(xmin), fabs(xmax));
	float order1=floor(log10(rmax));
	if(isnan(order1)||order1<-1000||fabs(order1)<=1){
		order1=0;
	}
	*order=order1;
	const float scale1=pow(10, -order1);
	xmax*=scale1;
	xmin*=scale1;
	diff*=scale1;
	float spacing;
	if(diff<rmax*1e-10){//handle special case
		*tic1=xmin;
		*dtic=xmax-xmin;
		*ntic=2;
	}else{
		if(logscale){
			spacing=1;
		} else{
			spacing=diff*0.1;
			float scale=pow(10., floor(log10(spacing)));
			int ratio=(int)ceil(spacing/scale);//value of the first significant digit rounded up
			//scale ratio up to keep maximum 10 tics
			if(ratio>5){
				ratio=10;
			}else if(ratio>2){
				ratio=5;
			}
			spacing=ratio*scale;//keep a single significant digit
		}

		*tic1=myfloor(xmin/spacing)*spacing;
		*dtic=spacing;
		*ntic=MAX(2, (int)(myceil(xmax/spacing)-myfloor(xmin/spacing)+1));
	}
	
}
/**
   adjust xmin, xmax properly.
*/
void round_limit(float* xmin, float* xmax, int logscale){
	if(*xmin==0 && *xmax==0) return;
	if(logscale){
		if(*xmin<=0) *xmin=0.1;
		if(*xmax<=0) *xmax=1;
		*xmin=log10(*xmin);
		*xmax=log10(*xmax);
	}
	if(fabs(*xmin)<1e-15&&fabs(*xmax)<1e-15){/*both are zero. */
		*xmin=0;
		*xmax=0;
	} else{
		float tic1, dtic;
		int ntic, order;
		calc_tic(&tic1, &dtic, &ntic, &order, *xmax, *xmin, 12, logscale);
		float xmin0=tic1*pow(10, order);
		float xmax0=(tic1+dtic*(ntic-1))*pow(10, order);
#if 1
		*xmin=xmin0;
		*xmax=xmax0;
#else
		if(fabs(xmin0-*xmin)<1e-5*xmin0){
			*xmin=xmin0;
		} else if(*xmin<xmin0){
			*xmin=xmin0-dtic*pow(10, order);
		}
		if(fabs(xmax0-*xmax)<1e-5*xmax0){
			*xmax=xmax0;
		} else if(*xmax>xmax0){
			*xmax=xmax0+dtic*pow(10, order);
		}
#endif
	}
	if(logscale){
		*xmin=pow(10, *xmin);
		*xmax=pow(10, *xmax);
	}
}
/**
 * When limit or limit0 changes, update zoomx,y and offx,y to not alter the limit0.
*/
void update_zoom(drawdata_t* drawdata){
	//info("update_zoom in: zoom=%g %g, off=%g %g\n", drawdata->zoomx, drawdata->zoomy, drawdata->offx, drawdata->offy);
	/*if(drawdata->zoomx==1&&drawdata->zoomy==1&&drawdata->offx==0&&drawdata->offy==0){
		return;
	}*/
	if((drawdata->limit0[0]==0&&drawdata->limit0[1]==0)||(drawdata->limit0[2]==0&&drawdata->limit0[3]==0)){
		return;//never been plotted. do not try to zoom/offset.
	}
	if(drawdata->limit_changed!=1){
		return;
	}
	/*limit0 matches limit in unzoomed state */
	int xlog=drawdata->xylog[0]=='y'?1:0;
	int ylog=drawdata->xylog[1]=='y'?1:0;
	float xmin, xmax, ymin, ymax;
	if(xlog){
		xmin=log10(drawdata->limit[0]);
		xmax=log10(drawdata->limit[1]);
		if(isinf(xmin)) xmin=0;
	} else{
		xmin=drawdata->limit[0];
		xmax=drawdata->limit[1];
	}
	if(ylog){
		ymin=log10(drawdata->limit[2]);
		ymax=log10(drawdata->limit[3]);
		if(isinf(ymin)) ymin=0;
	} else{
		ymin=drawdata->limit[2];
		ymax=drawdata->limit[3];
	}
	float diffx0=(xmax-xmin);
	float diffy0=(ymax-ymin);
	float midx0=(xmin+xmax)*0.5;
	float midy0=(ymin+ymax)*0.5;

	float diffx1=drawdata->limit0[1]-drawdata->limit0[0];
	float diffy1=drawdata->limit0[3]-drawdata->limit0[2];
	float midx1=(drawdata->limit0[1]+drawdata->limit0[0])*0.5;
	float midy1=(drawdata->limit0[3]+drawdata->limit0[2])*0.5;
	if(diffx1<=0) diffx1=1;
	if(diffy1<=0) diffy1=1;
	if(diffx1>diffx0*1e-5&&diffy1>diffy0*1e-5){/*limit allowable range */
	/*the new zoom */
		float ratiox=diffx0/diffx1; if(ratiox==0) ratiox=1;
		float ratioy=diffy0/diffy1; if(ratioy==0) ratioy=1;
		if(drawdata->square&&!drawdata->image){/*make the ratio equal. */
			if(fabs(ratiox-drawdata->zoomx)>1e-2&&fabs(ratioy-drawdata->zoomy)<1e-2){
			/*only x changed */
				ratioy=ratiox;
			} else if(fabs(ratiox-drawdata->zoomx)<1e-2&&fabs(ratioy-drawdata->zoomy)>1e-2){
				ratiox=ratioy;
			} else{
				ratiox=(ratiox+ratioy)*0.5;
				ratioy=ratiox;
			}
		}
		drawdata->zoomx=ratiox;
		drawdata->zoomy=ratioy;
		drawdata->offx=-(midx1-midx0)*drawdata->widthim/(diffx1*drawdata->zoomx);
		drawdata->offy=-(midy1-midy0)*drawdata->heightim/(diffy1*drawdata->zoomy);
	}

	//info("update_zoom: zoom=%g %g, off=%g %g\n", drawdata->zoomx, drawdata->zoomy, drawdata->offx, drawdata->offy);
}
/*
  Definition of style: (bits count from lowest end0
  bits 1-3:the point style.
  bit  4: whether points are connected.
  bits 5-8: size
  bits 9-32:color

  The color follows RGB representation: (Since 2011-02-18)
  bits 32-25: Red. bits 24-17: Green. bits 16-9: Blue.
 */
#define PARSE_STYLE(stylein)					\
    {								\
	style=stylein & 0x7;					\
	connectpts=(stylein&0x8)>>3;				\
	color=(stylein&0xFFFFFF00)>>8;				\
	sym_size=round(((stylein&0xF0)>>4));			\
	if(style>5) style=0;					\
	if(style==0) connectpts=1;				\
    }

static inline void
draw_point(cairo_t* cr, float ix, float iy, long style, float size){
	float size1=size+1;
	/*It is important to round ix, iy to have right symbols. */
	switch(style){
	case 0:/*nothing. just connect lines */
		break;
	case 1:/* o */
		cairo_new_sub_path(cr);
		cairo_arc(cr, ix-0.5, iy-0.5, size, 0, 2*M_PI);
		cairo_move_to(cr, ix, iy);
		break;
	case 2:/* x */
		cairo_move_to(cr, ix-size1, iy-size1);
		cairo_line_to(cr, ix+size, iy+size);
		cairo_move_to(cr, ix+size, iy-size1);
		cairo_line_to(cr, ix-size1, iy+size);
		cairo_move_to(cr, ix, iy);
		break;
	case 3:/* + */
#if DRAW_NEW == 1
		cairo_move_to(cr, ix-size1, iy);
		cairo_line_to(cr, ix+size, iy);
		cairo_move_to(cr, ix, iy-size1);
		cairo_line_to(cr, ix, iy+size);
		cairo_move_to(cr, ix, iy);
#else
	/*iy-1 is because flipping makes effective rounding of 0.5 differently. */
		cairo_move_to(cr, ix-size1, iy-1);
		cairo_line_to(cr, ix+size, iy-1);
		cairo_move_to(cr, ix, iy-size1);
		cairo_line_to(cr, ix, iy+size);
		cairo_move_to(cr, ix, iy);
#endif
		break;
	case 4:/* square [] */
		cairo_move_to(cr, ix-size1, iy-size1);
		cairo_line_to(cr, ix+size, iy-size1);
		cairo_line_to(cr, ix+size, iy+size);
		cairo_move_to(cr, ix+size, iy+size-1);
		cairo_line_to(cr, ix-size1, iy+size-1);
		cairo_move_to(cr, ix-size, iy+size-1);
		cairo_line_to(cr, ix-size, iy-size1);
		cairo_move_to(cr, ix, iy);
		break;
	case 5:/* . */
		cairo_new_sub_path(cr);
		cairo_arc(cr, ix-0.5, iy-0.5, 0., 0, 2*M_PI);
		cairo_arc(cr, ix-0.5, iy-0.5, 1, 0, 2*M_PI);
		cairo_move_to(cr, ix, iy);
		break;
	default:
		warning("Invalid style\n");
	}
}
/**
   Compute the limit from data for line plotting.
 */
static void update_limit(drawdata_t *drawdata){
	/* update max/minimum for both non-cumu and cumulative case. */
	//info("update_limit in\n");
	if((drawdata->limit_manual&&!drawdata->cumu)||drawdata->limit_changed!=-1) return;

	float xmin0=INFINITY, xmax0=-INFINITY, ymin0=INFINITY, ymax0=-INFINITY;
	for(int ipts=0; ipts<drawdata->npts; ipts++){
		const float *ptsx=drawdata->pts[ipts], *ptsy=0;
		if(!ptsx) continue;
		const int ptsnx=drawdata->ptsdim[ipts][0];
		const int ptsny=drawdata->ptsdim[ipts][1];
		float xmin, xmax, ymin=INFINITY, ymax=-INFINITY;
		if(ptsny>1){/*x is supplied */
			fmaxmin(ptsx, ptsnx, &xmax, &xmin);
			ptsy=ptsx+ptsnx;
		} else{/*x is index */
			xmin=0; xmax=(float)(ptsnx-1);
			ptsy=ptsx;
		}
		if(drawdata->cumu){
			int icumu=(int)drawdata->icumu;
			int ips0=0;
			if(icumu<ptsnx){
				ips0=icumu;
			}
			xmin=ips0;
			float y_cumu=0, y=0;

			int first=1;
			if(drawdata->cumuquad){
				for(int ips=ips0; ips<ptsnx; ips++){
					const float tmp=ptsy[ips];
					if(tmp!=0||first){
						y_cumu+=tmp*tmp;
						y=sqrt(y_cumu/(ips-ips0+1));
						if(y>ymax) ymax=y;
						if(y<ymin) ymin=y;
						first=0;
					}
				}
			} else{
				for(int ips=ips0; ips<ptsnx; ips++){
					const float tmp=ptsy[ips];
					if(tmp!=0||first){
						y_cumu+=tmp;
						y=y_cumu/(ips-ips0+1);
						if(y>ymax) ymax=y;
						if(y<ymin) ymin=y;
						first=0;
					}
				}
			}
		} else{
			fmaxmin(ptsy, ptsnx, &ymax, &ymin);
		}
		//info("xmin=%g, xmax=%g, ymin=%g, ymax=%g\n", xmin, xmax, ymin, ymax);
		if(xmin<xmin0) xmin0=xmin;
		if(xmax>xmax0) xmax0=xmax;
		if(ymin<ymin0) ymin0=ymin;
		if(ymax>ymax0) ymax0=ymax;
	}
	if(isinf(ymin0)) ymin0=0;
	if(isinf(ymax0)) ymax0=0;

	int xlog=drawdata->xylog[0]=='y'?1:0;
	int ylog=drawdata->xylog[1]=='y'?1:0;

	
	float xlimit=(xmax0-xmin0)*0.1;
	float ylimit=(ymax0-ymin0)*0.1;
	if(drawdata->cumu!=drawdata->cumulast){
		drawdata->limit_changed=3;//limit is changed. reset zoom.
		drawdata->limit[0]=xmin0;
		drawdata->limit[1]=xmax0;
		drawdata->limit[2]=ymin0;
		drawdata->limit[3]=ymax0;
	}else{
		drawdata->limit_changed=0;
		if(drawdata->limit[0]>xmin0||drawdata->limit[0]+xlimit<xmin0){
			drawdata->limit[0]=xmin0;//only update if the lower limit is below or a threshold above the old result
			drawdata->limit_changed=3;
		}
		if(drawdata->limit[1]<xmax0||drawdata->limit[1]>xlimit+xmax0){
			drawdata->limit[1]=xmax0;//only update if the upper limit is above or a threshold below the old result
			drawdata->limit_changed=3;
		}
		if(drawdata->limit[2]>ymin0||drawdata->limit[2]+ylimit<ymin0){
			drawdata->limit[2]=ymin0;//only update if the limit is below or a threshold above the old result
			drawdata->limit_changed=3;
		}
		if(drawdata->limit[3]<ymax0||drawdata->limit[3]>ylimit+ymax0){
			drawdata->limit[3]=ymax0;//only update if the upper limit is above or a threshold below the old result
			drawdata->limit_changed=3;
		}
	}
	if(drawdata->limit_changed==3){
		round_limit(&drawdata->limit[0], &drawdata->limit[1], xlog);
		round_limit(&drawdata->limit[2], &drawdata->limit[3], ylog);
	}
	//info("update_limit out:%g %g, %g %g, limit changed=%d\n", xmin0, xmax0, ymin0, ymax0, drawdata->limit_changed);
}
/**
   The master routine that draws in the cairo surface.

   meaning of limit_changed
   0: no change from previous plot.
   1: limit0 is changed and update_zoom is needed to update zoom.
   2: z limit is changed
   3: limit is changed and zoom will be reset.
   -1:new data, or switch between cumu and non-cumu. need to run update_limit.
*/
void cairo_draw(cairo_t* cr, drawdata_t* drawdata, int width, int height){
	if(!drawdata->ready) {
		dbg_time("data is not ready, cancelled.\n");
		return;
	}
	/*fill white background */
	//TIC;tic;
	drawdata->font_name_version=font_name_version;
	PangoLayout* layout=pango_cairo_create_layout(cr);
	pango_layout_set_font_description(layout, desc);
	/*
	if(!drawdata->image&&!drawdata->square){
		drawdata->cumu=cumu;
	} else{
		drawdata->cumu=0;
	}*/
	if(drawdata->cumu){
		/*if((int)drawdata->icumulast!=(int)drawdata->icumu){
			free(drawdata->limit_cumu);
			drawdata->limit_cumu=NULL;
		}*/
		if(!drawdata->limit_cumu){
			drawdata->limit_cumu=mycalloc(4, float);
		}
		drawdata->limit=drawdata->limit_cumu;
	} else{
		if(!drawdata->limit_data){
			drawdata->limit_data=mycalloc(4, float);
		}
		drawdata->limit=drawdata->limit_data;
	}
	if(!drawdata->image){
		if(drawdata->cumulast!=drawdata->cumu){
			drawdata->drawn=0;
			drawdata->limit_changed=-1;
		}else if(drawdata->cumu){
			if(drawdata->icumu!=drawdata->icumulast){
				drawdata->drawn=0;//just redraw, do not recompute limit
				//drawdata->limit[0]=drawdata->icumu;
				drawdata->limit_changed=-1;
			}
			if(drawdata->cumuquadlast!=drawdata->cumuquad){
				drawdata->drawn=0;//just redraw, do not recompute limit
			}
		}
		if(drawdata->limit_changed==-1){
			update_limit(drawdata);
		}
		if(drawdata->limit_changed==3){
			drawdata->offx=0;
			drawdata->offy=0;
			drawdata->zoomx=1;
			drawdata->zoomy=1;
			drawdata->drawn=0;
			drawdata->limit_changed=0;
		}
	}
	int xlog=drawdata->xylog[0]=='y'?1:0;
	int ylog=drawdata->xylog[1]=='y'?1:0;
	int widthim, heightim;
	static float xmin, xmax, ymin, ymax;
	if(xlog){
		xmin=log10(drawdata->limit[0]);
		xmax=log10(drawdata->limit[1]);
		if(isinf(xmin)) xmin=0;
	} else{
		xmin=drawdata->limit[0];
		xmax=drawdata->limit[1];
	}
	if(ylog){
		ymin=log10(drawdata->limit[2]);
		ymax=log10(drawdata->limit[3]);
		if(isinf(ymin)) ymin=0;
	} else{
		ymin=drawdata->limit[2];
		ymax=drawdata->limit[3];
	}
	float xmin0, ymin0, xmax0, ymax0;
	xmin0=xmin; ymin0=ymin; xmax0=xmax; ymax0=ymax;
	float xdim, ydim;/*dimension of the data. */
	if(drawdata->image){/*we are drawing an image. */
		xdim=(float)drawdata->nx;
		ydim=(float)drawdata->ny;
	} else{
		xdim=xmax-xmin;
		ydim=ymax-ymin;
	}
	if(xdim==0||!isfinite(xdim)) xdim=1;
	if(ydim==0||!isfinite(ydim)) ydim=1;
	//dbg("xdim=%g, ydim=%g\n", xdim, ydim);
	float sp_xr=20;
	if(drawdata->image){/*there is a colorbar */
		sp_xr=SP_XR;
	}
	float scalex=(float)(width-SP_XL-sp_xr)/xdim;
	float scaley=(float)(height-SP_YT-SP_YB)/ydim;

	if(drawdata->square){
		scaley=(scalex<scaley?scalex:scaley);
		scalex=scaley;
	}
	widthim=(int)(xdim*scalex/2)*2;
	heightim=(int)(ydim*scaley/2)*2;
	maxtic_x=widthim/(3*font_size);
	maxtic_y=heightim/(3*font_size);
	scalex=widthim/xdim;
	scaley=heightim/ydim;
	drawdata->widthim=widthim;
	drawdata->heightim=heightim;
	drawdata->scalex=scalex;
	drawdata->scaley=scaley;
	if(drawdata->limit_changed==1
		||(drawdata->widthim_last!=0
			&&(drawdata->widthim_last!=drawdata->widthim
				||drawdata->heightim_last!=drawdata->heightim))){
		 /*update_zoom: response to canvas resize or zoom/move called by the GUI*/
		update_zoom(drawdata);
		drawdata->limit_changed=0;
		drawdata->drawn=0;
	}
	if(drawdata->limit_changed==2&&drawdata->p){/*zlim changed. */
		const int nx=drawdata->nx;
		const int ny=drawdata->ny;
		const int stride=cairo_format_stride_for_width(drawdata->format, nx);

		flt2pix(nx, ny, !drawdata->gray, drawdata->p0, drawdata->p, drawdata->zlim, drawdata->zlog);
		cairo_surface_destroy(drawdata->image);
		drawdata->image=cairo_image_surface_create_for_data
			(drawdata->p, drawdata->format, nx, ny, stride);
		drawdata->limit_changed=0;
	}
	drawdata->widthim_last=drawdata->widthim;
	drawdata->heightim_last=drawdata->heightim;

	/*Offset in the cairo surface to draw the image. */
	const float xoff=round(((width-widthim-SP_XL-sp_xr)*0.5)+SP_XL);
	const float yoff=round(((height-heightim-SP_YT-SP_YB)*0.5)+SP_YT);
	drawdata->xoff=xoff;/*save for the GUI to use. */
	drawdata->yoff=yoff;
	/*center of the image on the screen. */
	drawdata->centerx=xoff+widthim*0.5;
	drawdata->centery=yoff+heightim*0.5;
	const float zoomx=drawdata->zoomx;/*Zoom of the image when displayed. */
	const float zoomy=drawdata->zoomy;

	cairo_select_font_face(cr, font_name, font_style, font_weight);
	cairo_set_font_size(cr, font_size);
	float linewidth=round(font_size*0.08);
	float ticlength=font_size*0.5;
	const float ticskip=drawdata->ticinside?5:ticlength;
	//float gridskip=drawdata->ticinside?ticlength:0;
	if(drawdata->ticinside){
		ticlength=-ticlength;
	}
	if(linewidth<1) linewidth=1;
	cairo_rectangle(cr, 0, 0, width, height);
	cairo_set_source_rgb(cr, 1, 1, 1);
	cairo_fill(cr);
	cairo_set_antialias(cr, CAIRO_ANTIALIAS_NONE);
	/*{
	  cairo_font_options_t *fonto= cairo_font_options_create();
	  cairo_font_options_set_hint_metrics(fonto,CAIRO_HINT_METRICS_ON);
	  cairo_font_options_set_hint_style (fonto, CAIRO_HINT_STYLE_MEDIUM);
	  cairo_font_options_set_antialias (fonto,CAIRO_ANTIALIAS_SUBPIXEL);
	  cairo_set_font_options(cr, fonto);
	  cairo_font_options_destroy(fonto);
	}*/
	cairo_set_line_width(cr, linewidth);
	/*Save the state of cairo before we drawing the image/points. */
	cairo_save(cr);
	/*flip upside down so that lower left is (0,0); */
	cairo_translate(cr, xoff, heightim+yoff);
	cairo_scale(cr, 1, -1);
	/*clip out an rectangular region to draw. */
	cairo_rectangle(cr, 0, 0, widthim, heightim);
	cairo_clip(cr);
	/*
		2021-09-20:
		Tried to draw tics before points usin cairo_push_group and
		cairo_pop_group. It didn't work because the points were drawing into an
		intermediate buffer to facilitate paning, but it covers all the grids.
	*/
	if(drawdata->nx && drawdata->ny){//draw 2-d image
		cairo_save(cr);
		cairo_scale(cr, scalex*zoomx, scaley*zoomy);
		/*
		  offx, offy are the offset in the cairo window.
		  ofx, ofy are the actual offset in the original data to display.
		 */
		float ofx=(drawdata->nx*0.5)*(1./zoomx-1.)+drawdata->offx/scalex;
		float ofy=(drawdata->ny*0.5)*(1./zoomy-1.)+drawdata->offy/scaley;
		/*The x and y patterns are negated and then set as
		  translation values in the pattern matrix.*/
		cairo_set_source_surface(cr, drawdata->image, ofx, ofy);
		if(scalex*zoomx>1){/*use nearest filter for up sampling to get clear images */
			cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_NEAREST);
		}
		cairo_paint(cr);
		/*cairo_reset_clip(cr); */
		float xdiff=xmax-xmin;
		float ydiff=ymax-ymin;
		/*
		  xmin, xmax is the real min/max of the x axis.
		  ofx is the offset.
		  nx, is the number of elements along x.
		  we can figure out xmin0, xmax0 from ofx and zoom.
		  We can also figure out ofx, zoom, from xmin0, xmax0
		 */
		xmin0=xmin-(ofx/drawdata->nx)*xdiff;
		xmax0=xmin0+xdiff/zoomx;
		ymin0=ymin-(ofy/drawdata->ny)*ydiff;
		ymax0=ymin0+ydiff/zoomy;
		cairo_restore(cr);
		//toc("cairo_draw image");
	}
	if(drawdata->npts>0){//plot points
		cairo_save(cr);
		float centerx=(xmax+xmin)*0.5;
		float centery=(ymax+ymin)*0.5;
		float ncx=widthim*0.5+drawdata->offx*zoomx;
		float ncy=heightim*0.5+drawdata->offy*zoomy;
		cairo_set_antialias(cr, CAIRO_ANTIALIAS_NONE);
		cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_NEAREST);
#if DRAW_NEW == 1
		int new_width0=(int)ceil(widthim*zoomx)+font_size*80;
		int new_height0=(int)ceil(heightim*zoomy)+font_size*80;
		int new_height=new_height0;
		int new_width=new_width0;
		/*Don't use too much buffer. Over flow memory. */
#define ZOOM_MAX 10
		if(zoomx>ZOOM_MAX){
			new_width=widthim*ZOOM_MAX;
		} else{
			drawdata->ncxoff=0;
		}
		if(zoomy>ZOOM_MAX){
			new_height=heightim*ZOOM_MAX;
		} else{
			drawdata->ncyoff=0;
		}

		int new_offx=(int)round((widthim-new_width)*0.5+drawdata->offx*zoomx+drawdata->ncxoff);
		int new_offy=(int)round((heightim-new_height)*0.5+drawdata->offy*zoomy+drawdata->ncyoff);
		if(drawdata->cacheplot&&(cairo_image_surface_get_width(drawdata->cacheplot)!=new_width||
			cairo_image_surface_get_height(drawdata->cacheplot)!=new_height)){
			cairo_surface_destroy(drawdata->cacheplot);
			drawdata->cacheplot=NULL;
		}

		int redraw=0;
		if(!drawdata->cacheplot){
			drawdata->cacheplot=cairo_image_surface_create(CAIRO_FORMAT_ARGB32, new_width, new_height);
			redraw=1;
		}
		if(fabs(drawdata->zoomx-drawdata->zoomxlast)>EPS||fabs(drawdata->zoomy-drawdata->zoomylast)>EPS
			||!drawdata->drawn){
			redraw=1;
		}
		if(new_width<new_width0&&(new_offx>0||-new_offx>(new_width-widthim))){
			drawdata->ncxoff=(int)(-drawdata->offx*zoomx);
			new_offx=(int)round((widthim-new_width)*0.5+drawdata->offx*zoomx+drawdata->ncxoff);
			redraw=1;
		}
		if(new_height<new_height0&&(new_offy>0||-new_offy>(new_height-heightim))){
			drawdata->ncyoff=(int)(-drawdata->offy*zoomy);
			new_offy=(int)round((heightim-new_height)*0.5+drawdata->offy*zoomy+drawdata->ncyoff);
			redraw=1;
		}
		/*center of plot. Important. */
		ncx=(float)new_width*0.5-drawdata->ncxoff;
		ncy=(float)new_height*0.5-drawdata->ncyoff;
		if(redraw){
			cairo_t* cr2=cr;
			cr=cairo_create(drawdata->cacheplot);
			cairo_set_antialias(cr, CAIRO_ANTIALIAS_NONE);
			cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_BILINEAR);
			cairo_set_line_width(cr, linewidth);
			/*Blank it first. */
			cairo_set_source_rgba(cr, 1., 1., 1., 1.);
			cairo_rectangle(cr, 0, 0, new_width, new_height);
			cairo_fill(cr);
#endif
			if(!drawdata->style_pts){
				drawdata->style_pts=mycalloc(drawdata->npts, int);/*save the styles for legend */
			}
			int color=0x0000FF;
			cairo_set_source_rgba(cr, 0.2, 0.0, 1.0, 1.0);/*Blue */
			int style;
			int connectpts;
			if(drawdata->square){/*we are plotting points. */
				style=3;/*default style is + */
				connectpts=0;
			} else{/*we are plotting curves */
				style=0;
				connectpts=1;
			}
			float sym_size=3;//round(3*sqrt(zoomx));
			if(drawdata->nstyle==1){
				PARSE_STYLE(drawdata->style[0]);
			}
			/*computed from below ix, iy formula by setting ix, iy to 0 and widthim or heightim */
			int icumu=(int)drawdata->icumu;
			for(int ipts=0; ipts<drawdata->npts; ipts++){
				const float *pts=drawdata->pts[ipts];
				const int ptsnx=drawdata->ptsdim[ipts][0];
				const int ptsny=drawdata->ptsdim[ipts][1];
				if(!pts||ptsnx==0||ptsny==0) continue;

				cairo_save(cr);
				cairo_set_antialias(cr, CAIRO_ANTIALIAS_NONE);
				if(drawdata->nstyle>1){
					PARSE_STYLE(drawdata->style[ipts]);
				} else if(drawdata->nstyle==0){
					color=default_color(ipts);
				}
				/*save the styles for legend */
				drawdata->style_pts[ipts]=style|connectpts<<3|color<<8|(int)(sym_size)<<4;
				set_color(cr, color);

				const float* ptsx, * ptsy=NULL;
				if(ptsny==2){
					ptsx=pts;
					ptsy=ptsx+ptsnx;
				} else{
					ptsx=0;
					ptsy=pts;
				}
				int ips0=0;
				if(drawdata->cumu&&icumu<ptsnx){
					ips0=icumu;
				}
				int ptstep=1;
				float ix=0, iy=0, y=0, y_cumu=0;
				if(connectpts){/*plot curves. */
					int ips=ips0;
					/*Draw first point */
					if(ptsx){/*don't do round here. */
						ix=(((xlog?log10(ptsx[ips]):ptsx[ips])-centerx)*scalex*zoomx+ncx);
					} else{
						ix=(((xlog?log10(ips):ips)-centerx)*scalex*zoomx+ncx);
					}
					if(isinf(ptsy[ips])){
						iy=0;
					}else{
						iy=(((ylog?log10(ptsy[ips]):ptsy[ips])-centery)*scaley*zoomy+ncy);
					}
					if(drawdata->cumuquad){
						y_cumu=ptsy[ips]*ptsy[ips];
					} else{
						y_cumu=ptsy[ips];
					}
					cairo_move_to(cr, ix, iy);
					/*connect additional points. */
					for(ips++; ips<ptsnx; ips+=ptstep){
						if(ptsx){/*don't do round here. */
							ix=(((xlog?log10(ptsx[ips]):ptsx[ips])-centerx)*scalex*zoomx+ncx);
						} else{
							ix=(((xlog?log10(ips):ips)-centerx)*scalex*zoomx+ncx);
						}
						y=ptsy[ips];

						if(drawdata->cumu){
							if(!isinf(y)){
								if(drawdata->cumuquad){
									y_cumu+=y*y;
									y=sqrt(y_cumu/(ips-ips0+1));
								} else{
									y_cumu+=y;
									y=y_cumu/(ips-ips0+1);
								}
							}
						} 
						if(isinf(y)) y=0;
					
						iy=(((ylog?log10(y):y)-centery)*scaley*zoomy+ncy);

						cairo_line_to(cr, round(ix), round(iy));
						if(ips%100==0){
							cairo_stroke(cr);
							cairo_move_to(cr, round(ix), round(iy));
						}
					}
					cairo_stroke(cr);
				}
				if(!(connectpts&&style==5)){/*plot points. */
					y_cumu=0;
					for(int ips=ips0; ips<ptsnx; ips+=ptstep){
						if(isinf(ptsy[ips])) continue;
						/*Map the coordinate to the image */
						if(ptsx){/*don't do round here. */
							ix=(((xlog?log10(ptsx[ips]):ptsx[ips])-centerx)*scalex*zoomx+ncx);
						} else{
							ix=(((xlog?log10(ips):ips)-centerx)*scalex*zoomx+ncx);
						}

						if(drawdata->cumu){
							if(ptsy[ips]!=0){
								if(drawdata->cumuquad){
									y_cumu+=ptsy[ips]*ptsy[ips];
									y=sqrt(y_cumu/(ips-ips0+1));
								} else{
									y_cumu+=ptsy[ips];
									y=y_cumu/(ips-ips0+1);
								}
							}
						} else{
							y=ptsy[ips];
						}
						iy=(((ylog?log10(y):y)-centery)*scaley*zoomy+ncy);
						draw_point(cr, ix, iy, style, sym_size);
						cairo_stroke(cr);/*stroke each point because color may change. */
					}/*ipts */
				}/*if */
				/*if(!drawdata->style){
				tic;
				dbg("stroke");
				cairo_stroke(cr);//stroke all together.
				toc("stroke");
				}*/
				if(ptsnx>0 && (drawdata->cumu || drawdata->legendcurve)){
					char val[20]={0};
					if(drawdata->legend&&drawdata->legend[ipts]&&connectpts){
						val[0]=drawdata->legend[ipts][0];
					}
					if(drawdata->cumu){
						snprintf(val, sizeof(val), "%c (%.2f)", val[0], y);
					}
					if(val[0]!=0){
						cairo_save(cr);
						cairo_translate(cr, ix+2, round((iy-font_size*0.5)/1)*1);
						cairo_scale(cr, 1, -1);
						pango_text(cr, layout, 0, 0, val, 0, 1, 0);
						cairo_restore(cr);
					}
				}
				cairo_restore(cr);
				//append cumulative average to end of legend.
				if(drawdata->cumu&&drawdata->legend[ipts]){
					char *tmp=drawdata->legend[ipts];
					char val[20]={0};
					snprintf(val, sizeof(val), " (%.2f)", y);
					char *old=strstr(tmp, " (");
					if(old){//remove old
						old[0]=0;
					}

					drawdata->legend[ipts]=stradd(tmp, val, NULL);
					free(tmp);
				}
			}/*iptsy */
#if DRAW_NEW == 1
			cairo_destroy(cr);
			cr=cr2;
		}
		cairo_set_antialias(cr, CAIRO_ANTIALIAS_NONE);
		/*cairo_pattern_set_filter(cairo_get_source(cr),CAIRO_FILTER_NEAREST); */
		cairo_set_source_surface(cr, drawdata->cacheplot, new_offx, new_offy);
		cairo_paint(cr);
#endif

		xmax0=(((widthim)*0.5)/zoomx-drawdata->offx)/scalex+centerx;
		xmin0=((-widthim*0.5)/zoomx-drawdata->offx)/scalex+centerx;
		ymax0=(((heightim)*0.5)/zoomy-drawdata->offy)/scaley+centery;
		ymin0=((-heightim*0.5)/zoomy-drawdata->offy)/scaley+centery;
		
		cairo_restore(cr);
		//toc("cairo_draw pts");
	}
	/*info("im: %d, %d min0: %g, %g, max0 %g, %g. zoom: %g, %g, scale %g, %g, off: %g, %g\n",
		widthim, heightim, xmin0, ymin0, xmax0, ymax0, zoomx, zoomy, scalex, scaley, drawdata->offx, drawdata->offy);*/
	if(drawdata->ncir>0){//plot circles
		cairo_save(cr);
		cairo_set_antialias(cr, CAIRO_ANTIALIAS_GRAY);

		float centerx=(xmax+xmin)/2-drawdata->offx/scalex;
		float centery=(ymax+ymin)/2-drawdata->offy/scaley;
		int ncx=widthim/2;
		int ncy=heightim/2;

		for(int icir=0; icir<drawdata->ncir; icir++){
			int color=(int)drawdata->cir[icir][3];
			set_color(cr, color);
			float ix=(drawdata->cir[icir][0]-centerx)*scalex*zoomx+ncx;
			float iy=(drawdata->cir[icir][1]-centery)*scaley*zoomy+ncy;
			cairo_arc(cr, ix, iy, drawdata->cir[icir][2]*scalex*zoomx, 0, M_PI*2);
			cairo_stroke(cr);
		}
		cairo_restore(cr);
	}
	cairo_restore(cr);//undo the clip and scale

	//now draw the tics 
	cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 1.0);
	/*Now doing the tic, and colobar */
	/*When there is no zoom, panning, limit0 equals to limit. */
	drawdata->limit0[0]=xmin0;
	drawdata->limit0[1]=xmax0;
	drawdata->limit0[2]=ymin0;
	drawdata->limit0[3]=ymax0;
	/*if(drawdata->spins){//dialog is running, update its values 
		for(int i=0; i<4; i++){//update spin button's value. 
			//gtk_spin_button_set_value(GTK_SPIN_BUTTON(drawdata->spins[i]), drawdata->limit0[i]);
		}
		if(drawdata->zlim[0] || drawdata->zlim[1]){
			for(int i=5; i<6; i++){//update spin button's value. 
				//gtk_spin_button_set_value(GTK_SPIN_BUTTON(drawdata->spins[i]), drawdata->zlim[i-4]);
			}
		}
	}*/
	char ticval[80];
	float tic1, dtic;
	int ntic, order;
	float sep;
	calc_tic(&tic1, &dtic, &ntic, &order, xmax0, xmin0, maxtic_x, xlog);
	sep=xmax0-xmin0;

	//draw the x axis
	for(int itic=0; itic<ntic; itic++){
		float ticv=(tic1+dtic*itic); if(fabs(ticv)<1e-6) ticv=0;
		float val=ticv*pow(10, order);
		float frac=(val-xmin0)/sep;
		float xpos=xoff+widthim*frac;

		if(frac>0.0&&frac<1){
			/*draw the vertical grid */
			if(drawdata->grid){
				cairo_save(cr);
				//cairo_set_dash(cr, stroke_dot, 2, gridskip);
				cairo_set_source_rgba(cr, 0.8, 0.8, 0.8, 1);
				cairo_move_to(cr, xpos, yoff);
				cairo_line_to(cr, xpos, yoff+heightim);
				cairo_stroke(cr);
				cairo_restore(cr);
			}
			/*draw the tic */
			cairo_move_to(cr, xpos, yoff+heightim);
			cairo_line_to(cr, xpos, yoff+heightim+ticlength);
			cairo_stroke(cr);
		}
		if(frac>=0&&frac<=1){
			snprintf(ticval, 80, "%g", (xlog?pow(10, ticv):ticv));
			pango_text(cr, layout, xpos, yoff+heightim+font_size*0.6+ticskip+1, ticval, 0.5, 0.5, 0);
		}
	}

	//draw the x axis minor ticks
	for(int itic=-1; itic<ntic; itic++){
		float ticv=(tic1+dtic*itic);
		for(int im=1; im<10; im++){
			float ticvm=ticv+(xlog?log10(1+im):im*0.1)*dtic;
			float valm=ticvm*pow(10, order);
			float fracm=(valm-xmin0)/sep;
			float xposm=xoff+widthim*fracm;
			if((fracm)>0.&&(fracm)<1){
				cairo_move_to(cr, xposm, yoff+heightim);
				cairo_line_to(cr, xposm, yoff+heightim+ticlength/2);
			}
		}
	}
	cairo_stroke(cr);
	//draw the x axis power
	if(order) pango_text_powindex(cr, layout, xoff+widthim-font_size*2, yoff+heightim+6+font_size*1.2, order, 0);
	calc_tic(&tic1, &dtic, &ntic, &order, ymax0, ymin0, maxtic_y, ylog);
	sep=ymax0-ymin0;

	//draw the y axis tic
	for(int itic=0; itic<ntic; itic++){
		float ticv=(tic1+dtic*itic); if(fabs(ticv)<1e-6) ticv=0;
		float val=ticv*pow(10, order);
		float frac=(val-ymin0)/sep;
		float ypos=yoff+heightim*(1-frac);
		if(frac>0.&&frac<1){
			/*draw the horizontal grid */
			if(drawdata->grid){
				cairo_save(cr);
				//cairo_set_dash(cr, stroke_dot, 2, gridskip);
				cairo_set_source_rgba(cr, 0.8, 0.8, 0.8, 1);
				cairo_move_to(cr, xoff, ypos);
				cairo_line_to(cr, xoff+widthim, ypos);
				cairo_stroke(cr);
				cairo_restore(cr);
			}
			/*draw the tic */
			cairo_move_to(cr, xoff, ypos);
			cairo_line_to(cr, xoff-ticlength, ypos);
			cairo_stroke(cr);
		}
		if(frac>=0&&frac<=1){
			snprintf(ticval, 80, "%g", (ylog?pow(10, ticv):ticv));
			pango_text(cr, layout, xoff-font_size*0.6-ticskip+1, ypos, ticval, 0.5, 0.5, 1);
		}
	}

	//draw y axis minor ticks
	for(int itic=-1; itic<ntic; itic++){
		float ticv=(tic1+dtic*itic);
		for(int im=1; im<10; im++){
			float ticvm=ticv+(ylog?log10(im+1):im*0.1)*dtic;
			float valm=ticvm*pow(10, order);
			float fracm=(valm-ymin0)/sep;
			float yposm=yoff+heightim*(1-fracm);
			if((fracm)>0.&&(fracm)<1){
				cairo_move_to(cr, xoff, yposm);
				cairo_line_to(cr, xoff-ticlength/2, yposm);
			}
		}
	}
	cairo_stroke(cr);
	//draws the y axis power
	if(order) pango_text_powindex(cr, layout, xoff-font_size*2.8, yoff+font_size*1.8, order, 1);

	if(drawdata->zlim[0]||drawdata->zlim[1]){/*draw colorbar */
		cairo_save(cr);
		cairo_translate(cr, xoff+widthim+SP_LEG, yoff);
		cairo_rectangle(cr, 0, 0, LEN_LEG, heightim);
		cairo_pattern_t* bar=cairo_pattern_create_linear(0, 0, 0, heightim);
		cairo_pattern_add_color_stop_rgb(bar, 0, 0.5625, 0, 0);
		cairo_pattern_add_color_stop_rgb(bar, 0.1111, 1, 0, 0);
		cairo_pattern_add_color_stop_rgb(bar, 0.3651, 1, 1, 0);
		cairo_pattern_add_color_stop_rgb(bar, 0.6190, 0, 1, 1);
		cairo_pattern_add_color_stop_rgb(bar, 0.8730, 0, 0, 1);
		cairo_pattern_add_color_stop_rgb(bar, 1, 0, 0, 0.5);
		cairo_set_source(cr, bar);
		cairo_fill(cr);
		cairo_pattern_destroy(bar);
		cairo_rectangle(cr, 0, 0, LEN_LEG, heightim);
		cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 1.0);
		cairo_stroke(cr);

		cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 1.0);
		calc_tic(&tic1, &dtic, &ntic, &order, drawdata->zlim[1], drawdata->zlim[0], maxtic_y, 0);
		sep=drawdata->zlim[1]-drawdata->zlim[0];

		for(int itic=0; itic<ntic; itic++){
			float ticv=tic1+dtic*itic; if(fabs(ticv)<1e-6) ticv=0;
			float val=ticv*pow(10, order);
			float frac;
			if(sep>1.e-10*fabs(drawdata->zlim[1])){
				frac=(val-drawdata->zlim[0])/sep;
			} else{
				if(itic==0) frac=0;
				else if(itic==1) frac=1;
				else frac=-1;
			}
			if(frac>0&&frac<1){
				cairo_move_to(cr, 0, heightim*(1-frac));
				cairo_line_to(cr, 4, heightim*(1-frac));
				cairo_move_to(cr, LEN_LEG, heightim*(1-frac));
				cairo_line_to(cr, LEN_LEG-4, heightim*(1-frac));
				cairo_stroke(cr);
			}
			if(frac>=0&&frac<=1){
				snprintf(ticval, 80, "%g", ticv);
				pango_text(cr, layout, LEN_LEG+4, heightim*(1-frac), ticval, 0, 0, 0);
			}
		}
		pango_text_powindex(cr, layout, LEN_LEG/2, -font_size*1.4-2, order, 0);
		cairo_restore(cr);
	}
	if(drawdata->title){
		pango_text(cr, layout, xoff+widthim/2, yoff-font_size*0.5-4, drawdata->title, 0.5, 0.5, 0);
	}
	if(drawdata->xlabel){
		pango_text(cr, layout, xoff+widthim/2, yoff+heightim+8+font_size*1.8, drawdata->xlabel, 0.5, 0.5, 0);
	}
	if(drawdata->ylabel){
		pango_text(cr, layout, xoff-font_size*1.8-6, yoff+heightim/2, drawdata->ylabel, 0.5, 0.5, 1);
	}
	if(drawdata->legend&&drawdata->npts&&drawdata->legendbox){
		int style, color, connectpts;
		float sym_size=0;
		cairo_save(cr);
		cairo_identity_matrix(cr);
		/*draw legend */
		char** legend=drawdata->legend;
		const int npts=drawdata->npts;
		const float linelen=30;/*length of line in legend if exist. */
		float textlen=0;/*maximum legend length. */
		float tall=0;
		float symlen=0;/*length of legend symbol */
		int npts_valid=0;
		/*first figure out the size required of longest legend entry. */
		for(int ipts=0; ipts<npts; ipts++){
			if(!legend[ipts]||!drawdata->pts[ipts]||!drawdata->ptsdim[ipts][0]||!drawdata->ptsdim[ipts][1]) continue;
			float legwidth, legheight;
			npts_valid++;
			pango_size(layout, legend[ipts], &legwidth, &legheight);
			textlen=MAX(textlen, legwidth);/*length of text. */
			tall=MAX(tall, legheight*1.2);/*tall of text. */
			PARSE_STYLE(drawdata->style_pts[ipts]);
			if(connectpts){
				symlen=MAX(symlen, linelen);
			} else{
				symlen=MAX(symlen, sym_size*2);
				tall=MAX(tall, sym_size*2);
			}
		}
		if(textlen){//legend is available
			const float legmarin=3;/*margin inside of box */
			const float legmarout=5;/*margin outside of box */
			const float symmargin=5;/*space before and after symbol */
			float legwidth=textlen+symlen+2*legmarin+symmargin*2;
			float legheight=tall*npts_valid+legmarin*2;
			cairo_translate(cr, xoff+legmarout+drawdata->legendoffx*(widthim-legwidth-2*legmarout),
				yoff+legmarout+drawdata->legendoffy*(heightim-legheight-2*legmarout));
			if(1){//box
				cairo_save(cr);
				cairo_rectangle(cr, 0, 0, legwidth, legheight);
				cairo_set_source_rgba(cr, 1.0, 1.0, 1.0, 1.0);
				cairo_fill_preserve(cr);
				cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 1.0);
				cairo_stroke(cr);
				cairo_restore(cr);
			}
			cairo_translate(cr, legmarin+symmargin, legmarin);
			for(int ipts=0; ipts<npts; ipts++){
				if(!legend[ipts]||!drawdata->pts[ipts]||!drawdata->ptsdim[ipts][0]||!drawdata->ptsdim[ipts][1]) continue;
				PARSE_STYLE(drawdata->style_pts[ipts]);
				set_color(cr, color);
				float ix=symlen*0.5;
				float iy=tall*0.5;
				draw_point(cr, ix, iy, style, sym_size);
				cairo_stroke(cr);
				if(connectpts){
					cairo_move_to(cr, 0, tall*0.5);
					cairo_line_to(cr, symlen, tall*0.5);
					cairo_stroke(cr);
				}
				cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 1.0);
				float toff=tall*1.1;
				cairo_move_to(cr, symlen+symmargin, tall*0.8);
				/*cairo_show_text(cr, legend[ig]); */
				pango_text(cr, layout, symlen+symmargin, toff-tall, legend[ipts], 0, 0, 0);
				cairo_translate(cr, 0, tall);
			}
		}
		cairo_restore(cr);
	}
	/*if(drawdata->dtime<10){
		cairo_save(cr);
		cairo_identity_matrix(cr);
		char fps[10];
		snprintf(fps, 10, "FPS:%5.1f", 1./drawdata->dtime);
		pango_text(cr, layout, font_size*0.3, font_size*0.2, fps, 0, 0, 0);
		cairo_restore(cr);
	}*/
	/*draw the border in the end to avoid covering it.*/
	cairo_rectangle(cr, xoff, yoff, widthim, heightim);
	cairo_stroke(cr);/*border */
	g_object_unref(layout);
	drawdata->icumulast=drawdata->icumu;
	drawdata->cumulast=drawdata->cumu;
	drawdata->cumuquadlast=drawdata->cumuquad;
	drawdata->zoomxlast=drawdata->zoomx;
	drawdata->zoomylast=drawdata->zoomy;
	drawdata->drawn=1;
	//toc("cairo_draw");
}

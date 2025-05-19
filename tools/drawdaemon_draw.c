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
#include <math.h>
#include "drawdaemon.h"
/*
  Routines in this file are about drawing in the cairo surface.

  2011-06-02: New idea to improve the effiency of drawdaemon:
  1) Only draw once, to a cached surface with no zoom, or displacement.
  2) Draw the cached surface to the real surface with zoom or displacement when user requests.
*/
const double stroke_dot[2]={1,5};
const double stroke_dash[2]={10,5};
const double stroke_solid[2]={10,0};

int hide_xlabel=0;
int hide_ylabel=0;
int hide_title=0;
int hide_legend=0;
int hide_colorbar=0;
#define DRAW_NEW 1 //1: draw to a cache surface. speed up draging.

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
static void pango_scale_fontsize(PangoLayout* layout, float fontscale){
	pango_font_description_set_size(desc, desc_font_size*fontscale);//use smaller font size for legends.
	pango_layout_set_font_description(layout, desc);//needed
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
   Drawing text. 
   @param x,y		reference point
   @param fracx		0: left, 0.5: middle, 1: right
   @param fracy		0: top, 0.5: middle, 1: bottom
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
	pango_cairo_show_layout(cr, layout);//top left corner of the layout is drawn at the current point of cairo context.
	cairo_restore(cr);
}

static void pango_text_power(cairo_t *cr, PangoLayout *layout, float x, float y, int power, float fracx, float fracy, int vertical){
	char powindex[40];
	snprintf(powindex, 40, "10<sup>%d</sup>", power);
	pango_text(cr, layout, x, y, powindex, fracx, fracy, vertical);
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
	if(isnan(order1)||order1<-1000||fabs(order1)<=2){
		order1=0;
	}
	*order=order1;
	const float scale1=pow(10, -order1);
	//convert values to between (-10 and 10)
	xmax*=scale1;
	xmin*=scale1;
	diff*=scale1;
	rmax*=scale1;
	float spacing;
	if(diff<rmax*1e-5){//handle special case
		*tic1=round((0.5*(xmax+xmin)-0.1)*10)*0.1;
		*dtic=0.2;
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
	//float oldmin=*xmin, oldmax=*xmax;
	if(logscale){
		if(*xmin<=0) *xmin=0.1;
		if(*xmax<=0) *xmax=*xmin;
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
		//info("xmin0=%g, xmax0=%g, tic1=%g, dtic=%g, ntic=%d\n", *xmin, *xmax, tic1, dtic, ntic);
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
	//info_time("round_limit: [%g, %g] --> [%g, %g]\n", *xmin, *xmax, oldmin, oldmax);
}
/**
 * When limit or limit0 changes, update zoomx,y and offx,y to not alter the limit0.
*/
void update_zoom(drawdata_t* drawdata){
	//dbg_time("update_zoom entr: zoom=%g %g, off=%g %g\n", drawdata->zoomx, drawdata->zoomy, drawdata->offx, drawdata->offy);
	if((drawdata->limit0[0]==0&&drawdata->limit0[1]==0)||(drawdata->limit0[2]==0&&drawdata->limit0[3]==0)){
		return;//never been plotted. do not try to zoom/offset.
	}
	/*limit0 matches limit in unzoomed state */
	int xlog=drawdata->xylog[0]=='y'?1:0;
	int ylog=drawdata->xylog[1]=='y'?1:0;
	float xmin, xmax, ymin, ymax;
	if(xlog){
		xmin=log10(drawdata->limit[0]);
		xmax=log10(drawdata->limit[1]);
		if(isinf(xmin)) xmin=xmax-5;
	} else{
		xmin=drawdata->limit[0];
		xmax=drawdata->limit[1];
	}
	if(ylog){
		ymin=log10(drawdata->limit[2]);
		ymax=log10(drawdata->limit[3]);
		if(isinf(ymin)) ymin=ymax-5;
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
		if(drawdata->square&&!drawdata->p){/*make the ratio equal. */
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
		drawdata->offx=-(midx1-midx0)*drawdata->widthim/(diffx1);
		drawdata->offy=-(midy1-midy0)*drawdata->heightim/(diffy1);
	}
	drawdata->update_zoom=0;
	//dbg_time("update_zoom exit: zoom=%g %g, off=%g %g\n", drawdata->zoomx, drawdata->zoomy, drawdata->offx, drawdata->offy);
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
#define PARSE_STYLE(stylein)				\
{											\
	style=stylein & 0x7;					\
	connectpts=(stylein&0x8)>>3;			\
	color=(stylein&0xFFFFFF00)>>8;			\
	sym_size=round((stylein&0xF0)>>4);		\
	if(style==0||style==7) connectpts=1;	\
	cairo_set_dash(cr, style==7?stroke_dash:stroke_solid, 2, 0);\
}

static inline void
draw_point(cairo_t* cr, float ix, float iy, long style, float size, float zoomx, int rounding){
	if(zoomx>1 && size<20){
		size*=sqrt(zoomx);
	}
	if(rounding){
		ix=round(ix);
		iy=round(iy);
		size=round(size);
	}
	
	/*It is important to round ix, iy to have right symbols. */
	switch(style){
	case 0:/*nothing. just connect lines */
		break;
	case 1:/* o */
		cairo_new_sub_path(cr);
		cairo_arc(cr, ix-0.5, iy-0.5, size, 0, 2*M_PI);
		break;
	case 2:/* x */
		cairo_move_to(cr, ix-size, iy-size);
		cairo_line_to(cr, ix+size, iy+size);
		cairo_move_to(cr, ix+size, iy-size);
		cairo_line_to(cr, ix-size, iy+size);
		break;
	case 3:/* + */
#if DRAW_NEW == 1
		cairo_move_to(cr, ix-size, iy);
		cairo_line_to(cr, ix+size, iy);
		cairo_move_to(cr, ix, iy-size);
		cairo_line_to(cr, ix, iy+size);
#else
	/*iy-1 is because flipping makes effective rounding of 0.5 differently. */
		cairo_move_to(cr, ix-size1, iy-1);
		cairo_line_to(cr, ix+size, iy-1);
		cairo_move_to(cr, ix, iy-size1);
		cairo_line_to(cr, ix, iy+size);
#endif
		break;
	case 4:/* square [] */
		cairo_move_to(cr, ix-size, iy-size);
		cairo_line_to(cr, ix+size, iy-size);
		cairo_line_to(cr, ix+size, iy+size);
		cairo_move_to(cr, ix+size, iy+size);
		cairo_line_to(cr, ix-size, iy+size);
		cairo_move_to(cr, ix-size, iy+size);
		cairo_line_to(cr, ix-size, iy-size);
		break;
	case 5:/* . */
		cairo_new_sub_path(cr);
		cairo_arc(cr, ix-0.5, iy-0.5, 0., 0, 2*M_PI);
		cairo_arc(cr, ix-0.5, iy-0.5, 1, 0, 2*M_PI);
		break;
	default:
		break;//treat as 0
	}
}
/**
   Compute the limit from data for line plotting.
 */
static void update_limit(drawdata_t *drawdata){
	/* update max/minimum for both non-cumu and cumulative case. */
	//if(drawdata->p || (drawdata->limit_manual&&!drawdata->cumu)||drawdata->limit_changed!=-1) return;
	const int xlog=drawdata->xylog[0]=='y'?1:0;
	const int ylog=drawdata->xylog[1]=='y'?1:0;
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
			if (xlog && xmin==0){
				if(ptsx[0]==0 && ptsx[1]>0){//frequency index
					xmin=ptsx[1];
				}else{
					xmin=xmax/1e8;
				}
			}
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

			for(int ips=ips0; ips<ptsnx; ips++){
				const float tmp=ptsy[ips];
				if(isfinite(tmp)){
					if(drawdata->cumuquad){
						y_cumu+=tmp*tmp;
						y=sqrt(y_cumu/(ips-ips0+1));
					}else{
						y_cumu+=tmp;
						y=y_cumu/(ips-ips0+1);
					}
					if(y>ymax) ymax=y;
					if(y<ymin) ymin=y;
				}
			}
		} else{
			//fmaxmin(ptsy, ptsnx, &ymax, &ymin);
			for(int ips=0; ips<ptsnx; ips++){
				const float y=ptsy[ips];
				if(isfinite(y)){
					if(y>ymax) ymax=y;
					if(y<ymin) ymin=y;
				}
			}
		}
		//info("xmin=%g, xmax=%g, ymin=%g, ymax=%g\n", xmin, xmax, ymin, ymax);
		if(xmin<xmin0) xmin0=xmin;
		if(xmax>xmax0) xmax0=xmax;
		if(ymin<ymin0) ymin0=ymin;
		if(ymax>ymax0) ymax0=ymax;
	}
	for(int icir=0; icir<drawdata->ncir; icir++){
		xmax0=MAX(xmax0, drawdata->cir[icir][0]+drawdata->cir[icir][2]);
		xmin0=MIN(xmin0, drawdata->cir[icir][0]-drawdata->cir[icir][2]);
		ymax0=MAX(ymax0, drawdata->cir[icir][1]+drawdata->cir[icir][2]);
		ymin0=MIN(ymin0, drawdata->cir[icir][1]-drawdata->cir[icir][2]);
	}
	if(isinf(ymin0)) ymin0=0;
	if(isinf(ymax0)) ymax0=1;
	int limit_changed=0;
	float ylimit=(ymax0-ymin0)*0.1;
	/*if(drawdata->cumu!=drawdata->cumulast){
		drawdata->update_zoom=2;//limit is changed. reset zoom.
		drawdata->limit[0]=xmin0;
		drawdata->limit[1]=xmax0;
		drawdata->limit[2]=ymin0;
		drawdata->limit[3]=ymax0;
	}else*/
	{
		if(drawdata->limit[0]!=xmin0){
			drawdata->limit[0]=xmin0;
			limit_changed=1;
		}
		if(drawdata->limit[1]!=xmax0){
			drawdata->limit[1]=xmax0;
			limit_changed=1;
		}
		if(!drawdata->limit[2]||drawdata->limit[2]>ymin0||drawdata->limit[2]+ylimit<ymin0){
			drawdata->limit[2]=ymin0;//only update if the limit is below or a threshold above the old result
			limit_changed=1;
		}
		if(!drawdata->limit[3]||drawdata->limit[3]<ymax0||drawdata->limit[3]>ylimit+ymax0){
			drawdata->limit[3]=ymax0;//only update if the upper limit is above or a threshold below the old result
			limit_changed=1;
		}
	}
	round_limit(&drawdata->limit[0], &drawdata->limit[1], xlog);
	round_limit(&drawdata->limit[2], &drawdata->limit[3], ylog);
	if(drawdata->update_zoom!=2 && limit_changed==1&&!(drawdata->zoomx==1&&drawdata->zoomy==1&&drawdata->offx==0&&drawdata->offy==0)){
		drawdata->update_zoom=1;//update zoom to preserve viewing area
	}
	drawdata->update_limit=0;
	//dbg_time("update_limit out:%g %g, %g %g, limit changed=%d\n", xmin0, xmax0, ymin0, ymax0, limit_changed);
}
/*static double get_scale(cairo_t* cr){
	cairo_surface_t *sur=cairo_get_target(cr);
	double xs, ys;
	cairo_surface_get_device_scale(sur, &xs, &ys);
	return (xs+ys)*0.5;
}*/

void thread_cleanup(void *drawdata){
	((drawdata_t*)drawdata)->thread=0;
}
void* cairo_draw_thread(drawdata_t* drawdata){
	pthread_cleanup_push(thread_cleanup, drawdata);

repeat:
	//info_time("cairo_draw entr: limit_changed=%d, width=%d, height=%d\n", drawdata->limit_changed, width ,height);
	if(drawdata->frame_io!=drawdata->frame_draw){
		drawdata->drawn=0;
		//info("%s %s: frame_io=%d, frame_draw=%d\n", drawdata->fig, drawdata->name, drawdata->frame_io, drawdata->frame_draw);
		drawdata->frame_draw=drawdata->frame_io;
	}
	int width=drawdata->width;
	int height=drawdata->height;
	int iframe=drawdata->iframe;
	/*fill white background */
	//TIC;tic;
	drawdata->font_name_version=font_name_version;
	
	cairo_t *crt=NULL;
	cairo_t *cr=NULL;
	pthread_mutex_lock(&drawdata->mutex);
	if(drawdata->surface){//do not use a second pixmap as it maybe vector
		cr=cairo_create(drawdata->surface);
	}else{
		crt=cairo_create(drawdata->pixmap);//maintain a reference to surface.
		if(width!=drawdata->pwidth2 || height!=drawdata->pheight2){
			if(drawdata->pixmap2){
				cairo_surface_destroy(drawdata->pixmap2);
			}
			drawdata->pixmap2=NULL;
			drawdata->drawn=0;
		}
		if(!drawdata->pixmap2){
			drawdata->pwidth2=width;
			drawdata->pheight2=height;
			//int scale=cairo_surface_get_device_scale(target)
			drawdata->pixmap2=cairo_surface_create_similar(drawdata->pixmap, CAIRO_CONTENT_COLOR_ALPHA, width, height);
			//cairo_surface_set_device_scale(drawdata->pixmap2, scale, scale);
		}
		cr=cairo_create(drawdata->pixmap2);
	}
	pthread_mutex_unlock(&drawdata->mutex);

	PangoLayout* layout=pango_cairo_create_layout(cr);
	pango_layout_set_font_description(layout, desc);
	cairo_set_antialias(cr, CAIRO_ANTIALIAS_DEFAULT);
	
	if(drawdata->p||drawdata->square){
		drawdata->cumu=0;
	}
	if(drawdata->cumu){
		/*if((int)drawdata->icumulast!=(int)drawdata->icumu){
			free(drawdata->limit_cumu);
			drawdata->limit_cumu=NULL;
		}*/
		if(!drawdata->limit_cumu){
			drawdata->limit_cumu=mycalloc(4, float);
			drawdata->update_limit=1;
		}
		drawdata->limit=drawdata->limit_cumu;
	} else{
		if(!drawdata->limit_data){
			drawdata->limit_data=mycalloc(4, float);
			drawdata->update_limit=1;
		}
		drawdata->limit=drawdata->limit_data;
	}
	float *zlim=drawdata->zlim+(drawdata->zlog?2:0);
	if(drawdata->p){
		if(drawdata->zlog_last!=drawdata->zlog){
			if(drawdata->zlim_manual){
				if(drawdata->zlog){
					drawdata->zlim[2]=log10(drawdata->zlim[0]);
					drawdata->zlim[3]=log10(drawdata->zlim[1]);
				} else{
					drawdata->zlim[0]=pow(10, drawdata->zlim[2]);
					drawdata->zlim[1]=pow(10, drawdata->zlim[3]);
				}
			}
			drawdata->zlog_last=drawdata->zlog;
			drawdata->zlim_changed=1;
		}
	}else{
		if(drawdata->cumulast!=drawdata->cumu
			||(drawdata->cumu&&(drawdata->icumu!=drawdata->icumulast||drawdata->cumuquadlast!=drawdata->cumuquad))){
			drawdata->drawn=0;
			drawdata->update_limit=1;
			drawdata->update_zoom=2;//reset zoom
		}
		if(drawdata->update_limit){
			update_limit(drawdata);
		}
	}
	if(drawdata->update_zoom==2){
		drawdata->offx=0;
		drawdata->offy=0;
		drawdata->zoomx=1;
		drawdata->zoomy=1;
	}
	
	int xlog=drawdata->xylog[0]=='y'?1:0;
	int ylog=drawdata->xylog[1]=='y'?1:0;
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
	if(drawdata->p){/*we are drawing an image. */
		xdim=(float)drawdata->nx;
		ydim=(float)drawdata->ny;
	} else{
		xdim=xmax-xmin;
		ydim=ymax-ymin;
	}
	if(xdim==0||!isfinite(xdim)) xdim=1;
	if(ydim==0||!isfinite(ydim)) ydim=1;
	//dbg("xdim=%g, ydim=%g\n", xdim, ydim);
	/*int SP_XR=drawdata->p?(!hide_colorbar?SP_XR:10):20;
	int SP_XL=hide_ylabel?0:SP_XL;
	int SP_YT=hide_title?0:SP_YT;
	int SP_YB=hide_xlabel?0:SP_YB;*/
	
	const float linewidth=MAX(1, round(font_size*0.08));
	const float ticlength=(drawdata->ticinside?-1:1)*font_size*0.5;
	const float W_CB=25;//width of colorbar
	const float S_CB=15;//separation between image and colorbar
	const float W_CBL=font_size*3;//width of colorbar label
	const float S_CBL=font_size*0.1; //separation between colorbar tic and text
	const float H_LBL=font_size*1.3; //height of xlabel, ylabel or title
	const float S_LBL=font_size*0.4; //separation at the top and bottom of label
	const float S_TIC=MAX(ticlength, 0)+S_LBL;
	const float SP_XL=hide_ylabel?S_LBL:(S_LBL+H_LBL+S_LBL+H_LBL+S_TIC);
	const float SP_YB=hide_xlabel?S_LBL:(S_LBL+H_LBL+S_LBL+H_LBL+S_TIC);
	const float SP_YT=hide_title?S_LBL:(S_LBL+H_LBL+S_LBL);
	const float SP_XR=(!drawdata->p || hide_colorbar)?S_CB:(S_CB+W_CB+S_TIC+S_CBL+W_CBL);
	/*float SP_XL=font_size*2.6+12;
	float SP_YT=font_size*1.3+12;
	float SP_YB=font_size*2.6+12;
	float SP_XR=SP_LEG+W_LEG+font_size*3;*/
	/*Offset in the cairo surface to draw the image. */

	//compute valid drawing area.
	float scalex0=(float)(width-SP_XL-SP_XR)/xdim;
	float scaley0=(float)(height-SP_YT-SP_YB)/ydim;

	if(drawdata->square){
		scalex0=scaley0=MIN(scalex0,scaley0);
	}
	const int widthim=(int)(xdim*scalex0/2)*2; //data area bounded by axes
	const int heightim=(int)(ydim*scaley0/2)*2;//data area bounded by axes
	
	const float scalex=widthim/xdim;
	const float scaley=heightim/ydim;
	const int maxtic_x=widthim/(3*font_size);
	const int maxtic_y=heightim/(3*font_size);
	
	drawdata->widthim=widthim;
	drawdata->heightim=heightim;
	drawdata->scalex=scalex;
	drawdata->scaley=scaley;

	
	const float xoff=round(((width-widthim-SP_XL-SP_XR)*0.5)+SP_XL);
	const float yoff=round(((height-heightim-SP_YT-SP_YB)*0.5)+SP_YT);

	drawdata->xoff=xoff;/*save for the GUI to use. */
	drawdata->yoff=yoff;
	/*center of the image on the screen. */
	drawdata->centerx=xoff+widthim*0.5;
	drawdata->centery=yoff+heightim*0.5;

	if(drawdata->update_zoom!=2){
		if((drawdata->widthim_last>0&&drawdata->widthim_last!=drawdata->widthim)
			||(drawdata->heightim_last>0 &&drawdata->heightim_last!=drawdata->heightim)){
			drawdata->update_zoom=1;
		}
		if(drawdata->update_zoom==1){
			/*Update zoomx/y, offx/y to preserve limit0 displayed. responses to window resize or zoom by selection.*/
			update_zoom(drawdata);
		}
	}
	drawdata->widthim_last=drawdata->widthim;
	drawdata->heightim_last=drawdata->heightim;

	if(drawdata->zlim_changed&&drawdata->p){/*zlim changed. */
		const int nx=drawdata->nx;
		const int ny=drawdata->ny;
		const cairo_format_t format=(cairo_format_t)(drawdata->gray?CAIRO_FORMAT_A8:CAIRO_FORMAT_ARGB32);
		const int stride=cairo_format_stride_for_width(format, nx);
		flt2pix(drawdata->p0, drawdata->p, nx, ny, drawdata->gray,
			drawdata->zlim+(drawdata->zlog?2:0), drawdata->zlim_manual, drawdata->zlog);
		cairo_surface_destroy(drawdata->image);
		drawdata->image=cairo_image_surface_create_for_data
			(drawdata->p, format, nx, ny, stride);
		drawdata->zlim_changed=0;
	}
	
	
	const float zoomx=drawdata->zoomx;/*Zoom of the image when displayed. */
	const float zoomy=drawdata->zoomy;

	//cairo_select_font_face(cr, font_name, font_style, font_weight);//not used by pango
	//cairo_set_font_size(cr, font_size);//not used by pango
	
	//float gridskip=drawdata->ticinside?ticlength:0;
	
	cairo_rectangle(cr, 0, 0, width, height);
	cairo_set_source_rgb(cr, 1, 1, 1);
	cairo_fill(cr);
	
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
	const float xdiff=xmax-xmin;
	const float ydiff=ymax-ymin;
	float centerx=0, centery=0;
	
	//compute the data range shown
	if(drawdata->nx&&drawdata->ny){//draw 2-d image
		centerx=(drawdata->nx*0.5)*(1./zoomx-1.)+drawdata->offx/scalex/zoomx;
		centery=(drawdata->ny*0.5)*(1./zoomy-1.)+drawdata->offy/scaley/zoomy;
		xmin0=xmin-(centerx/drawdata->nx)*xdiff;
		xmax0=xmin0+xdiff/zoomx;
		ymin0=ymin-(centery/drawdata->ny)*ydiff;
		ymax0=ymin0+ydiff/zoomy;
	}else{
		centerx=(xmax+xmin)*0.5;
		centery=(ymax+ymin)*0.5;
		xmax0=(((widthim)*0.5)-drawdata->offx)/zoomx/scalex+centerx;
		xmin0=((-widthim*0.5)-drawdata->offx)/zoomx/scalex+centerx;
		ymax0=(((heightim)*0.5)-drawdata->offy)/zoomy/scaley+centery;
		ymin0=((-heightim*0.5)-drawdata->offy)/zoomy/scaley+centery;
	}
	/*When there is no zoom, panning, limit0 equals to limit. */
	drawdata->limit0[0]=xmin0;
	drawdata->limit0[1]=xmax0;
	drawdata->limit0[2]=ymin0;
	drawdata->limit0[3]=ymax0;
	char ticval[80];
	float tic1, dtic, sep;
	int ntic, order;
	{//draw the grid and tics.
		cairo_save(cr);
		cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 1.0);

		calc_tic(&tic1, &dtic, &ntic, &order, xmax0, xmin0, maxtic_x, xlog);
		sep=xmax0-xmin0;
		//draw the x axis
		for(int itic=0; itic<ntic; itic++){
			float ticv=(tic1+dtic*itic); if(fabs(ticv)<1e-6) ticv=0;
			float val=ticv*pow(10, order);
			float frac=(val-xmin0)/sep;
			float xpos=xoff+widthim*frac;

			if(frac>0.&&frac<1.){
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
				if(!hide_xlabel){/*draw the tic */
					cairo_move_to(cr, xpos, yoff+heightim);
					cairo_line_to(cr, xpos, yoff+heightim+ticlength);
					cairo_stroke(cr);
				}
			}
			if(frac>=-0.01&&frac<=1.01&&!hide_xlabel){
				if(xlog){
					snprintf(ticval, 80, "10<sup>%g</sup>", ticv);
				} else{
					snprintf(ticval, 80, "%g", ticv);
				}
				float ypos=yoff+heightim+S_TIC;
				//float ypos=yoff+heightim+font_size*0.6+ticskip+1;
				pango_text(cr, layout, xpos, ypos, ticval, 0.5, 0, 0);
			}
		}
		if(!hide_xlabel){
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
			float ypos=yoff+heightim+S_TIC+H_LBL+S_LBL;
			if(order) pango_text_power(cr, layout, xoff+widthim-S_LBL, ypos, order, 1, 0, 0);
			if(drawdata->xlabel){
				pango_text(cr, layout, xoff+widthim/2, ypos, drawdata->xlabel, 0.5, 0, 0);
			}
		}

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
					cairo_set_source_rgba(cr, 0.6, 0.6, 0.6, 1);
					cairo_move_to(cr, xoff, ypos);
					cairo_line_to(cr, xoff+widthim, ypos);
					cairo_stroke(cr);
					cairo_restore(cr);
				}
				if(!hide_ylabel){
					/*draw the tic */
					cairo_move_to(cr, xoff, ypos);
					cairo_line_to(cr, xoff-ticlength, ypos);
					cairo_stroke(cr);
				}
			}
			if(frac>=-0.01&&frac<=1.01&&!hide_ylabel){
				if(ylog){
					snprintf(ticval, 80, "10<sup>%g</sup>", ticv);
				} else{
					snprintf(ticval, 80, "%g", ticv);
				}
				float xpos=xoff-S_TIC;
				pango_text(cr, layout, xpos, ypos, ticval, 1, 0.5, 1);
			}
		}
		if(!hide_ylabel){
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
			float xpos=xoff-S_TIC-H_LBL-S_LBL;
			if(order) pango_text_power(cr, layout, xpos, yoff+S_LBL, order, 1, 1, 1);
			if(drawdata->ylabel){
				pango_text(cr, layout, xpos, yoff+heightim/2, drawdata->ylabel, 1, 0.5, 1);
			}
		}
		if(drawdata->title&&!hide_title){
			pango_text(cr, layout, xoff+widthim/2, yoff-S_LBL, drawdata->title, 0.5, 1, 0);
			if(drawdata->filename_gif){
				char frame[100];
				snprintf(frame, sizeof(frame), "<small>%d</small>", drawdata->frame_draw);
				pango_text(cr, layout, xoff+widthim, yoff-S_LBL, frame, 1, 1, 0);
			}
		}
		cairo_restore(cr);
	}
	char **legend=noellipsis?drawdata->legend:drawdata->legend_ellipsis;

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
		
		/*The x and y patterns are negated and then set as
		  translation values in the pattern matrix.*/
		cairo_set_source_surface(cr, drawdata->image, centerx, centery);
		cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_NEAREST);//use after set_source_surface
		//if(scalex*zoomx>1){/*use nearest filter for up sampling to get clear images */
		//	cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_NEAREST);
		//}
		cairo_paint(cr);
		/*cairo_reset_clip(cr); */
		
		/*
		  xmin, xmax is the real min/max of the x axis.
		  ofx is the offset.
		  nx, is the number of elements along x.
		  we can figure out xmin0, xmax0 from ofx and zoom.
		  We can also figure out ofx, zoom, from xmin0, xmax0
		 */
		
		cairo_restore(cr);
		//toc("cairo_draw image");
	}else if(drawdata->npts>0){//plot points
		cairo_save(cr);
		
		float ncx=widthim*0.5+drawdata->offx;
		float ncy=heightim*0.5+drawdata->offy;
		//cairo_set_antialias(cr, CAIRO_ANTIALIAS_NONE);
		//cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_NEAREST);
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
		int redraw=0;
		int new_offx=(int)round((widthim-new_width)*0.5+drawdata->offx+drawdata->ncxoff);
		int new_offy=(int)round((heightim-new_height)*0.5+drawdata->offy+drawdata->ncyoff);
		if(drawdata->cacheplot&&(drawdata->cache_width!=new_width||drawdata->cache_height!=new_height)){
			redraw+=1;//window resized
		}
		if(fabs(drawdata->zoomx-drawdata->zoomx_last)>1e-3||fabs(drawdata->zoomy-drawdata->zoomy_last)>1e-3){
			redraw+=2;//zoom changed
		}
		if(new_width<new_width0&&(new_offx>0||-new_offx>(new_width-widthim))){
			drawdata->ncxoff=(int)(-drawdata->offx);
			new_offx=(int)round((widthim-new_width)*0.5+drawdata->offx+drawdata->ncxoff);
			//redraw+=4;
		}
		if(new_height<new_height0&&(new_offy>0||-new_offy>(new_height-heightim))){
			drawdata->ncyoff=(int)(-drawdata->offy);
			new_offy=(int)round((heightim-new_height)*0.5+drawdata->offy+drawdata->ncyoff);
			//redraw+=8;
		}
		if(redraw||!drawdata->drawn||!drawdata->cacheplot){
			if(drawdata->cacheplot) cairo_surface_destroy(drawdata->cacheplot);
			drawdata->cacheplot=cairo_surface_create_similar(cairo_get_target(cr), CAIRO_CONTENT_COLOR_ALPHA, new_width, new_height);
			/*double x_scale, y_scale;
			cairo_surface_get_device_scale(cairo_get_target(cr), &x_scale, &y_scale);
			drawdata->cacheplot=cairo_image_surface_create(CAIRO_FORMAT_ARGB32, new_width*x_scale, new_height*y_scale);
			cairo_surface_set_device_scale(drawdata->cacheplot, x_scale, y_scale);*/
			drawdata->cache_width=new_width;
			drawdata->cache_height=new_height;
			if(!redraw) redraw+=16;
		}

		/*center of plot. Important. */
		ncx=(float)new_width*0.5-drawdata->ncxoff;
		ncy=(float)new_height*0.5-drawdata->ncyoff;
		if(redraw){
			//info_time("redraw=%d\n", redraw);
			cairo_t* cr2=cr;
			cr=cairo_create(drawdata->cacheplot);
			cairo_set_antialias(cr, CAIRO_ANTIALIAS_DEFAULT);
			cairo_set_line_width(cr, linewidth);
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
			float sym_size=4;
			if(drawdata->nstyle==1){
				PARSE_STYLE(drawdata->style[0]);
			}
			/*computed from below ix, iy formula by setting ix, iy to 0 and widthim or heightim */
			for(int ipts=0; ipts<drawdata->npts; ipts++){
				const float *pts=drawdata->pts[ipts];
				const int ptsnx=drawdata->ptsdim[ipts][0];
				const int ptsny=drawdata->ptsdim[ipts][1];
				if(!pts||ptsnx==0||ptsny==0) continue;

				cairo_save(cr);
				//cairo_set_antialias(cr, CAIRO_ANTIALIAS_NONE);
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
				int icumu=(int)drawdata->icumu;
				if(icumu>=ptsnx) icumu=0;
				int nlayer=(drawdata->cumu && drawdata->cumuover)?2:1;
				float ix=0, iy=0, y=0, y_cumu=0;
				for(int ilayer=0; ilayer<nlayer; ilayer++){
					int cumu=drawdata->cumu&&(!drawdata->cumuover || ilayer==1);
					if((drawdata->cumuover && ilayer==1)||style==7){
						cairo_set_dash(cr, stroke_dash, 2, 0);
					}
					int ips=cumu?icumu:0;
					int ptstep=1;
					/*Draw first point */
					/*don't do round here. */
					ix=((xlog?log10(ptsx?ptsx[ips]:ips):(ptsx?ptsx[ips]:ips))-centerx)*scalex*zoomx+ncx; if(!isfinite(ix)) ix=0;
					iy=((ylog?log10(ptsy[ips]):ptsy[ips])-centery)*scaley*zoomy+ncy; 
					if(!isfinite(iy)) iy=0;
					if(cumu){
						if(drawdata->cumuquad){
							y_cumu=ptsy[ips]*ptsy[ips];
						} else{
							y_cumu=ptsy[ips];
						}
					}
					cairo_move_to(cr, ix, iy);
					if(style){
						draw_point(cr, ix, iy, style, sym_size, zoomx, 1);
					}
					/*connect additional points. */
					for(ips++; ips<ptsnx; ips+=ptstep){
						ix=((xlog?log10(ptsx?ptsx[ips]:ips):(ptsx?ptsx[ips]:ips))-centerx)*scalex*zoomx+ncx;		
						y=ptsy[ips];
						if(!isfinite(y)) y=0;
						if(cumu){
							if(drawdata->cumuquad){
								y_cumu+=y*y;
								y=sqrt(y_cumu/(ips-icumu+1));
							} else{
								y_cumu+=y;
								y=y_cumu/(ips-icumu+1);
							}
						}

						iy=((ylog?log10(y):y)-centery)*scaley*zoomy+ncy;
						if(connectpts){
							cairo_line_to(cr, ix, iy);
						}else{
							cairo_move_to(cr, ix, iy);
						}
						if(style){
							draw_point(cr, ix, iy, style, sym_size, zoomx, 1);
						}
					}
					cairo_stroke(cr);
				}
				//append cumulative average to end of legend.
				if(legend[ipts]){
					char* old=strstr(legend[ipts], " (");
					if(old){//remove old value.
						old[0]=0;
					}
					if(drawdata->cumu){//add cumulative rms to the end of legend 
						char *tmp=legend[ipts];
						char val[20]={0};
						snprintf(val, sizeof(val), " (%.2f)", y);
						legend[ipts]=stradd(tmp, val, NULL);
						free(tmp);
					}
				}
				//plot legend name on curve
				if(ptsnx>0 && drawdata->legendcurve && connectpts){
					/*char val[20]={0};
					if(legend&&legend[ipts]&&connectpts){
						val[0]=legend[ipts][0];
					}
					if(drawdata->cumu){
						snprintf(val, sizeof(val), "%c (%.2f)", val[0], y);
					}*/
					char *val=legend[ipts];
					if(val && val[0]!=0){
						cairo_save(cr);
						cairo_translate(cr, ix+5, iy); //round((iy-font_size*0.5)/1)*1);
						cairo_scale(cr, 1, -1);
						pango_text(cr, layout, 0, 0, val, 0, 0.5, 0);
						cairo_restore(cr);
					}
				}
				cairo_restore(cr);

			}/*iptsy */
#if DRAW_NEW == 1
			cairo_destroy(cr);
			cr=cr2;
		}
		cairo_set_source_surface(cr, drawdata->cacheplot, new_offx, new_offy);
		cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_NEAREST);
		cairo_paint(cr);
#endif
		cairo_restore(cr);
		//toc("cairo_draw pts");
	}
	/*info("im: %d, %d min0: %g, %g, max0 %g, %g. zoom: %g, %g, scale %g, %g, off: %g, %g\n",
		widthim, heightim, xmin0, ymin0, xmax0, ymax0, zoomx, zoomy, scalex, scaley, drawdata->offx, drawdata->offy);*/
	if(drawdata->ncir>0){//plot circles
		cairo_save(cr);
		//cairo_set_antialias(cr, CAIRO_ANTIALIAS_GRAY);
		centerx=(xmax+xmin)/2-drawdata->offx/zoomx/scalex;
		centery=(ymax+ymin)/2-drawdata->offy/zoomy/scaley;
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

	if((zlim[0]||zlim[1])&&!hide_colorbar){/*draw colorbar */
		cairo_save(cr);
		cairo_translate(cr, xoff+widthim+S_CB, yoff);
		cairo_rectangle(cr, 0, 0, W_CB, heightim);
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
		cairo_rectangle(cr, 0, 0, W_CB, heightim);
		cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 1.0);
		cairo_stroke(cr);

		cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 1.0);
		calc_tic(&tic1, &dtic, &ntic, &order, zlim[1], zlim[0], maxtic_y, drawdata->zlog);
		sep=zlim[1]-zlim[0];

		for(int itic=0; itic<ntic; itic++){
			float ticv=tic1+dtic*itic; if(fabs(ticv)<1e-6) ticv=0;
			float val=ticv*pow(10, order);
			float frac;
			if(sep>1.e-10*fabs(zlim[1])){
				frac=(val-zlim[0])/sep;
			} else{
				if(itic==0) frac=0;
				else if(itic==1) frac=1;
				else frac=-1;
			}
			if((drawdata->ticinside&&frac>0&&frac<1)||(!drawdata->ticinside&&frac>=-0.01&&frac<=1.01)){
				//cairo_move_to(cr, 0, heightim*(1-frac));
				//cairo_line_to(cr, 4, heightim*(1-frac));
				cairo_move_to(cr, W_CB, heightim*(1-frac));
				cairo_line_to(cr, W_CB+ticlength, heightim*(1-frac));
				cairo_stroke(cr);
			}
			if(frac>=-0.01&&frac<=1.01){
				float xpos=W_CB+S_TIC;
				if(drawdata->zlog){
					pango_text_power(cr, layout, xpos, heightim *(1-frac), ticv, 0, 0.5, 0);
				}else{
					snprintf(ticval, 80, "%g", drawdata->zlog?pow(10,ticv):ticv);
					pango_text(cr, layout, xpos, heightim*(1-frac), ticval, 0, 0.5, 0);
				}
			}
		}
		if(order) pango_text_power(cr, layout, W_CB/2, -S_LBL-H_LBL*0.5, order, .5, .5, 0);
		cairo_restore(cr);
	}

	if(legend&&drawdata->npts&&drawdata->legendbox){
		int style, color, connectpts;
		float sym_size=0;
		cairo_save(cr);
		cairo_identity_matrix(cr);
		pango_scale_fontsize(layout, 0.8);
		/*draw legend */
		const int npts=drawdata->npts;
		const float linelen=30;/*length of line in legend if exist. */
		float textlen=0;/*maximum legend length. */
		float tall=0;
		float symlen=0;/*length of legend symbol */
		int npts_valid=0;
		/*first figure out the size required of longest legend entry. */
		cairo_save(cr);
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
		cairo_restore(cr);
		if(textlen){//legend is available
			const float legmarin=3;/*margin inside of box */
			const float legmarout=5;/*margin outside of box */
			const float symmargin=5;/*space before and after symbol */
			//same the legend box information for gesture control of it.
			drawdata->legbox_width=textlen+symlen+2*legmarin+symmargin*2;
			drawdata->legbox_height=tall*npts_valid+legmarin*2;
			drawdata->legbox_ox=xoff+legmarout+drawdata->legendoffx*(widthim-drawdata->legbox_width-2*legmarout);
			drawdata->legbox_oy=yoff+legmarout+drawdata->legendoffy*(heightim-drawdata->legbox_height-2*legmarout);
			drawdata->legbox_rx=1./(widthim-drawdata->legbox_width-2*legmarout);
			drawdata->legbox_ry=1./(heightim-drawdata->legbox_height-2*legmarout);
			cairo_translate(cr, drawdata->legbox_ox, drawdata->legbox_oy);
			if(1){//box
				cairo_save(cr);
				cairo_rectangle(cr, 0, 0, drawdata->legbox_width, drawdata->legbox_height);
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
				draw_point(cr, ix, iy, style, sym_size, 1, 1);
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
		pango_scale_fontsize(layout, 1);
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
	cairo_set_source_rgba(cr, 0, 0, 0, 1);
	cairo_stroke(cr);/*border */
	g_object_unref(layout);
	cairo_destroy(cr);
	drawdata->icumulast=drawdata->icumu;
	drawdata->cumulast=drawdata->cumu;
	drawdata->cumuquadlast=drawdata->cumuquad;
	drawdata->zoomx_last=drawdata->zoomx;
	drawdata->zoomy_last=drawdata->zoomy;
	drawdata->update_zoom=0;
	drawdata->update_limit=0;
	drawdata->drawn=1;
	
	if(crt){
		cairo_set_source_surface(crt, drawdata->pixmap2, 0, 0);
		cairo_paint(crt);
		cairo_destroy(crt);
	}
	if(drawdata->filename_gif && drawdata->frame_draw!=drawdata->frame_gif){
		drawdata->frame_gif=drawdata->frame_draw;
		char filename[PATH_MAX];
		snprintf(filename, PATH_MAX, "%s/%03d.png", drawdata->filename_gif,drawdata->frame_gif);
		//info("filename is %s\n", filename);
		cairo_surface_write_to_png(drawdata->pixmap, filename);
	}
	
	if(drawdata->drawarea){
		g_idle_add((GSourceFunc)drawarea_refresh, drawdata);
	}

	//if(iframe<drawdata->iframe){
	//	info("cairo_draw: iframe has changed from %d to %d, repeat\n", iframe, drawdata->iframe);
	//}
	//toc("cairo_draw");
	if(iframe<drawdata->iframe){
		//info("cairo_draw_thread: finished iframe %d, %d is pending\n", iframe, drawdata->iframe);
		goto repeat;
	}else{
		//info("cairo_draw_thread: finished iframe %d, no more pending\n", iframe);
	}
	pthread_cleanup_pop(1);
	return NULL;
}
/**
   The master routine that draws in the cairo surface.

   meaning of limit_changed
   0: no change from previous plot.
   1: limit0 is changed and update_zoom is needed to update zoom.
   -1:new data, or switch between cumu and non-cumu. need to run update_limit.
*/
void cairo_draw(drawdata_t* drawdata){
	//layout=pango_cairo_create_layout(cr);
	//pango_layout_set_font_description(layout, desc);
	if(!drawdata->ready || drawdata->width==0 || drawdata->width == 0) {
		dbg_time("data is not ready or size is 0, cancelled.\n");
		return;
	}
	if(drawdata->surface){//for file saving, do not thread
		drawdata->drawn=0;
		cairo_draw_thread(drawdata);
	}else if(drawdata->thread){
		//info_time("thread is already running, skip. todo: add a time out.\n");
	}else{
		drawdata->thread=thread_new((thread_fun)cairo_draw_thread, drawdata);
		/*pthread_attr_t attr;
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, 1);
		pthread_create(&drawdata->thread, &attr, cairo_draw, drawdata);*/
	}
	//g_object_unref(layout);
}

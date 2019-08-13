/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "drawdaemon.h"
#include <math.h>
/*
  Routines in this file are about drawing in the cairo surface.

  2011-06-02: New idea to improve the effiency of drawdaemon:
  1) Only draw once, to a cached surface with no zoom, or displacement.
  2) Draw the cached surface to the real surface with zoom or displacement when user requests.
*/
const double stroke_dot[2]={1,5};
//const double stroke_dash[2]={10,10};
//const double stroke_solid[2]={10,0};
int maxtic_x=12;/*Maximum number of tics along x */
int maxtic_y=12;/*Maximum number of tics along y */

double SP_XL;/*space reserved for ylabel */
double SP_YT;/*space reserved for title */
double SP_YB;/*space reserved for xlabel */
double SP_XR;/*space reserved for legend */
#define DRAW_NEW 1

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
/**
   Drawing text. frac=0: left align. frac=0.5: center. frac=1: right align
*/
static void pango_text(cairo_t *cr, PangoLayout *layout, double x, double y, 
		       const char *text, double fracx, double fracy, int vertical){
    cairo_save(cr);
    pango_layout_set_markup(layout, text, -1);
    if(fracx>0 || fracy>0){
	int width, height;
	pango_layout_get_size(layout, &width, &height);
	if(vertical){
	    cairo_move_to(cr, x-height*fracx/PANGO_SCALE, y+width*fracy/PANGO_SCALE);
	}else{
	    cairo_move_to(cr, x-width*fracx/PANGO_SCALE,  y-height*fracy/PANGO_SCALE);
	}
    }else{
	cairo_move_to(cr, x, y);
    }
    if(vertical){
	cairo_rotate(cr, -M_PI/2);
    }
    pango_cairo_update_layout(cr, layout);
    pango_cairo_show_layout(cr, layout);
    cairo_restore(cr);
}

static void pango_text_powindex(cairo_t *cr, PangoLayout *layout, double x, double y, int order, int vertical){
    if(order==0) return;
    char powindex[40];
    snprintf(powindex, 40,"10<sup>%d</sup>", order);
    pango_text(cr, layout, x, y, powindex, 0, 0, vertical);
}


static void calc_tic(double *tic1, double *dtic, int *ntic, int *order, 
		     double xmax, double xmin, int maxtic, int logscale){
    //when logscale is true, the xmin, xmax are already log of data.
    //dbg("xmin=%g, xmax=%g, \n", xmin, xmax);
    (void)maxtic;
    double diff=xmax-xmin;
    /*first get the order of magnitude. */
    double rmax=MAX(fabs(xmin),fabs(xmax));
    double order1=floor(log10(rmax));
    if(isnan(order1) || order1<-1000 || fabs(order1)<=2){
	order1=0;
    }
    const double scale1=pow(10, -order1);
    xmax*=scale1;
    xmin*=scale1;
    diff*=scale1;
    
    double spacing;
    if(logscale || diff<1e-3){
	spacing=1;
    }else{
	spacing=diff*0.1;
	double scale=pow(10., floor(log10(spacing)));
	int ratio=(int)round(spacing/scale);
	const double ratio2[11]={1, 1, 1, 5, 5, 5, 5, 5, 10, 10, 10};
	if(ratio<0 || ratio>10){
	    warning("ratio=%d, spacing=%g\n", ratio, spacing);
	    ratio=1;
	}
	spacing=ratio2[ratio]*scale;
    }
	
    *tic1=myfloor(xmin/spacing)*spacing;
    *dtic=spacing;
    *ntic=(int)(myceil(xmax/spacing)-myfloor(xmin/spacing)+1);
    /*if(*ntic<2 || (*tic1<xmin && *tic1+spacing>xmax)) {
	dbg("*ntic=%d\n", *ntic);
	*ntic=2;
	*dtic=xmax-xmin;
	*tic1=xmin;
	}*/
    *order=(int)order1;
    if(*ntic<2){
	dbg("xmin=%g, xmax=%g, diff=%g, tic1=%g, dtic=%g, ntic=%d, order1=%g\n",
	    xmin, xmax, diff, *tic1, *dtic, *ntic, order1);
    }
}
/**
   adjust xmin, xmax properly.
*/
void round_limit(double *xmin, double *xmax, int logscale){
    if(logscale){
	if(*xmin<=0) *xmin=0.1;
	if(*xmax<=0) *xmax=1;
	*xmin=log10(*xmin);
	*xmax=log10(*xmax);
    }
    if(fabs(*xmin)<1e-15 && fabs(*xmax)<1e-15){/*both are zero. */
	*xmin=-1;
	*xmax=1;
    }else{
	double tic1, dtic;
	int ntic, order;
	calc_tic(&tic1, &dtic, &ntic, &order, *xmax, *xmin, 12, logscale);
	double xmin0=tic1*pow(10,order);
	double xmax0=(tic1+dtic*(ntic-1))*pow(10,order);
#if 1
	*xmin=xmin0;
	*xmax=xmax0;
#else
	if(fabs(xmin0-*xmin)<1e-5*xmin0){
	    *xmin=xmin0;
	}else if(*xmin < xmin0){
	    *xmin=xmin0-dtic*pow(10,order);
	}
	if(fabs(xmax0-*xmax)<1e-5*xmax0){
	    *xmax=xmax0;
	}else if(*xmax > xmax0){
	    *xmax=xmax0+dtic*pow(10,order);
	}
#endif
    }
    if(logscale){
	*xmin=pow(10, *xmin);
	*xmax=pow(10, *xmax);
    }
}
/**
   convert new limit to zoom and off using updated limit0.
*/
void apply_limit(drawdata_t *drawdata){
    /*limit0 matches limit in unzoomed state */
    int xlog=drawdata->xylog[0]=='n'?0:1;
    int ylog=drawdata->xylog[1]=='n'?0:1;
    double xmin, xmax, ymin, ymax;
    if(xlog){
	xmin=log10(drawdata->limit[0]);
	xmax=log10(drawdata->limit[1]);
	if(isinf(xmin)) xmin=0;
    }else{
	xmin=drawdata->limit[0];
	xmax=drawdata->limit[1];
    }
    if(ylog){
	ymin=log10(drawdata->limit[2]);
	ymax=log10(drawdata->limit[3]);
	if(isinf(ymin)) ymin=0;
    }else{
	ymin=drawdata->limit[2];
	ymax=drawdata->limit[3];
    }
    double diffx0=(xmax-xmin);
    double diffy0=(ymax-ymin);
    double midx0=(xmin+xmax)*0.5;
    double midy0=(ymin+ymax)*0.5;

    double diffx1=drawdata->limit0[1]-drawdata->limit0[0];
    double diffy1=drawdata->limit0[3]-drawdata->limit0[2];
    double midx1=(drawdata->limit0[1]+drawdata->limit0[0])*0.5;
    double midy1=(drawdata->limit0[3]+drawdata->limit0[2])*0.5;
    if(diffx1<=0) diffx1=1;
    if(diffy1<=0) diffy1=1;
    if(diffx1 > diffx0*1e-5 && diffy1 > diffy0 *1e-5){/*limit allowable range */
	/*the new zoom */
	double ratiox=diffx0/diffx1; if(ratiox==0) ratiox=1;
	double ratioy=diffy0/diffy1; if(ratioy==0) ratioy=1;
	if(drawdata->square){/*make the ratio equal. */
	    if(fabs(ratiox-drawdata->zoomx)>1e-2 && fabs(ratioy-drawdata->zoomy)<1e-2){
		/*only x changed */
		ratioy=ratiox;
	    }else if(fabs(ratiox-drawdata->zoomx)<1e-2 && fabs(ratioy-drawdata->zoomy)>1e-2){
		ratiox=ratioy;
	    }else{
		ratiox=(ratiox+ratioy)*0.5;
		ratioy=ratiox;
	    }
	}
	drawdata->zoomx=ratiox;
	drawdata->zoomy=ratioy;
	drawdata->offx=-(midx1-midx0)*drawdata->widthim/(diffx1*drawdata->zoomx);
	drawdata->offy=-(midy1-midy0)*drawdata->heightim/(diffy1*drawdata->zoomy);
    }
    
    //dbg("zoom=%g %g, off=%g %g\n", drawdata->zoomx, drawdata->zoomy, drawdata->offx, drawdata->offy);
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
draw_point(cairo_t *cr, double ix, double iy, long style, double size){
    double size1=size+1;
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
	cairo_move_to(cr,ix-size1,iy-size1);
	cairo_line_to(cr,ix+size,iy+size);
	cairo_move_to(cr,ix+size,iy-size1);
	cairo_line_to(cr,ix-size1,iy+size);
	cairo_move_to(cr, ix, iy);
	break;
    case 3:/* + */
#if DRAW_NEW == 1
	cairo_move_to(cr,ix-size1,iy);
	cairo_line_to(cr,ix+size, iy);
	cairo_move_to(cr,ix,iy-size1);
	cairo_line_to(cr,ix,iy+size);
	cairo_move_to(cr, ix, iy);
#else
	/*iy-1 is because flipping makes effective rounding of 0.5 differently. */
	cairo_move_to(cr,ix-size1,iy-1);
	cairo_line_to(cr,ix+size,iy-1);
	cairo_move_to(cr,ix,iy-size1);
	cairo_line_to(cr,ix,iy+size);
	cairo_move_to(cr, ix, iy);
#endif
	break;
    case 4:/* square [] */
	cairo_move_to(cr,ix-size1,iy-size1);
	cairo_line_to(cr,ix+size,iy-size1);
	cairo_line_to(cr,ix+size,iy+size);
	cairo_move_to(cr,ix+size,iy+size-1);
	cairo_line_to(cr,ix-size1,iy+size-1);
	cairo_move_to(cr,ix-size,iy+size-1);
	cairo_line_to(cr,ix-size,iy-size1);
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
void update_limit(drawdata_t *drawdata){
    /*need to update max/minimum. */
    drawdata->limit_changed=0;
    if(drawdata->cumulast!=drawdata->cumu){
	drawdata->offx=0;
	drawdata->offy=0;
	drawdata->zoomx=1;
	drawdata->zoomy=1;
    }
    double gain=1;
    if(drawdata->cumu){
	if(!drawdata->limit_cumu){
	    drawdata->limit_cumu=mycalloc(4,double);
	}
	drawdata->limit=drawdata->limit_cumu;
    }else{
	if(!drawdata->limit_data){
	    drawdata->limit_data=mycalloc(4,double);
	}
	drawdata->limit=drawdata->limit_data;
    }
    /*if(drawdata->dtime<1){
	gain=drawdata->dtime;
	}*/
    double xmin0=INFINITY, xmax0=-INFINITY, ymin0=INFINITY, ymax0=-INFINITY;
    for(int ipts=0; ipts<drawdata->npts; ipts++){
	const double *ptsx=drawdata->pts[ipts], *ptsy=0;
	if(!ptsx) continue;
	const int ptsnx=drawdata->ptsdim[ipts][0];
	const int ptsny=drawdata->ptsdim[ipts][1];
	double xmin, xmax, ymin=INFINITY, ymax=-INFINITY;
	if(ptsny>1){/*x is supplied */
	    dmaxmin(ptsx, ptsnx, &xmax, &xmin);
	    ptsy=ptsx+ptsnx;
	}else{/*x is index */
	    xmin=0; xmax=(double)(ptsnx-1);
	    ptsy=ptsx;
	}
	if(drawdata->cumu){
	    int icumu=(int)drawdata->icumu;
	    int ips0=0;
	    if(icumu<ptsnx){
		ips0=icumu;
	    }
	    xmin=ips0;
	    double y_cumu=0,y=0;
		  
	    int first=1;
	    if(drawdata->cumuquad){
		for(int ips=ips0; ips<ptsnx; ips++){
		    const double tmp=ptsy[ips];
		    if(tmp!=0 || first){
			y_cumu+=tmp*tmp;
			y=sqrt(y_cumu/(ips-ips0+1));
			if(y>ymax) ymax=y;
			if(y<ymin) ymin=y;
			first=0;
		    }
		}
	    }else{
		for(int ips=ips0; ips<ptsnx; ips++){
		    const double tmp=ptsy[ips];
		    if(tmp!=0 || first){
			y_cumu+=tmp;
			y=y_cumu/(ips-ips0+1);
			if(y>ymax) ymax=y;
			if(y<ymin) ymin=y;
			first=0;
		    }
		}
	    } 
	}else{
	    dmaxmin(ptsy, ptsnx, &ymax, &ymin);
	}
	if(xmin<xmin0) xmin0=xmin;
	if(xmax>xmax0) xmax0=xmax;
	if(ymin<ymin0) ymin0=ymin;
	if(ymax>ymax0) ymax0=ymax;
    }
    if(isinf(ymin0)) ymin0=0;
    if(isinf(ymax0)) ymax0=0;

    int xlog=drawdata->xylog[0]=='n'?0:1;
    int ylog=drawdata->xylog[1]=='n'?0:1;
    if(!xlog) round_limit(&xmin0, &xmax0, xlog);
    if(!ylog) round_limit(&ymin0, &ymax0, ylog);
    //dbg("xmin0=%g, xmax0=%g, ymin0=%g, ymax0=%g\n", xmin0, xmax0, ymin0, ymax0);
    drawdata->limit[0]=xmin0;
    drawdata->limit[1]=xmax0;
    drawdata->limit[2]=ymin0*gain+(1-gain)*drawdata->limit[2];
    drawdata->limit[3]=ymax0*gain+(1-gain)*drawdata->limit[3];
}

/**
   The master routine that draws in the cairo surface.
*/
void cairo_draw(cairo_t *cr, drawdata_t *drawdata, int width, int height){
    /*fill white background */
    //TIC;tic;
    drawdata->font_name_version=font_name_version;
    PangoLayout *layout=pango_cairo_create_layout(cr);
    pango_layout_set_font_description(layout, desc);
    cairo_rectangle(cr,0,0,width,height);
    cairo_set_source_rgb(cr,1,1,1);
    cairo_fill(cr);
    cairo_set_antialias(cr,CAIRO_ANTIALIAS_NONE);
    /*{
      cairo_font_options_t *fonto= cairo_font_options_create();
      cairo_font_options_set_hint_metrics(fonto,CAIRO_HINT_METRICS_ON);
      cairo_font_options_set_hint_style (fonto, CAIRO_HINT_STYLE_MEDIUM);
      cairo_font_options_set_antialias (fonto,CAIRO_ANTIALIAS_SUBPIXEL);
      cairo_set_font_options(cr, fonto);
      cairo_font_options_destroy(fonto);
      }*/
    if(!drawdata->image && !drawdata->square){
	drawdata->cumu=cumu;
    }else{
	drawdata->cumu=0;
    }
    if(drawdata->cumu){
	if((int)drawdata->icumulast==(int)drawdata->icumu){
	    drawdata->limit=drawdata->limit_cumu;
	}else{
	    free(drawdata->limit_cumu);
	    drawdata->limit=drawdata->limit_cumu=NULL;
	}
    }else{
	drawdata->limit=drawdata->limit_data;
    }
    if(!drawdata->limit){
	update_limit(drawdata);
    }
    if(drawdata->cumulast!=drawdata->cumu || drawdata->cumuquadlast != drawdata->cumuquad){
	drawdata->drawn=0;
    }
    int xlog=drawdata->xylog[0]=='n'?0:1;
    int ylog=drawdata->xylog[1]=='n'?0:1;

    int widthim, heightim;
    static double xmin, xmax, ymin, ymax;
    if(xlog){
	xmin=log10(drawdata->limit[0]);
	xmax=log10(drawdata->limit[1]);
	if(isinf(xmin)) xmin=0;
    }else{
	xmin=drawdata->limit[0];
	xmax=drawdata->limit[1];
    }
    if(ylog){
	ymin=log10(drawdata->limit[2]);
	ymax=log10(drawdata->limit[3]);
	if(isinf(ymin)) ymin=0;
    }else{
	ymin=drawdata->limit[2];
	ymax=drawdata->limit[3];
    }
    double xmin0,ymin0,xmax0,ymax0;
    xmin0=xmin; ymin0=ymin; xmax0=xmax; ymax0=ymax;
    double xdim, ydim;/*dimension of the data. */
    if(drawdata->image){/*we are drawing an image. */
	xdim=(double)drawdata->nx;
	ydim=(double)drawdata->ny;
    }else{
	xdim=xmax-xmin;
	ydim=ymax-ymin;
    }
    if(xdim==0) xdim=1;
    if(ydim==0) ydim=1;
    //dbg("xdim=%g, ydim=%g\n", xdim, ydim);
    double sp_xr=20;
    if(drawdata->image){/*there is a colorbar */
	sp_xr=SP_XR;
    }
    double scalex = (double)(width-SP_XL-sp_xr)/xdim;
    double scaley = (double)(height-SP_YT-SP_YB)/ydim;
  
    if(drawdata->square){
	double scale  = (scalex<scaley?scalex:scaley);
	scalex = scaley = scale;
    }
    widthim=(int)(xdim*scalex/2)*2;
    heightim=(int)(ydim*scaley/2)*2;
    maxtic_x=widthim/(3*font_size);
    maxtic_y=heightim/(3*font_size);
    scalex = widthim/xdim;
    scaley = heightim/ydim;
    drawdata->widthim=widthim;
    drawdata->heightim=heightim;
    drawdata->scalex=scalex;
    drawdata->scaley=scaley;
    if(drawdata->limit_changed == 1
       || (drawdata->widthim_last !=0 
	   && (drawdata->widthim_last!=drawdata->widthim 
	       || drawdata->heightim_last!=drawdata->heightim))){
	/*canvas is resized, need to adjust zoom/paning */
	apply_limit(drawdata);
	drawdata->limit_changed=0;
	drawdata->drawn=0;
    }
    if(drawdata->limit_changed == 2 && drawdata->p){/*zlim changed. */
	int nx=drawdata->nx;
	int ny=drawdata->ny;
	int stride=cairo_format_stride_for_width(drawdata->format, nx);

	dbl2pix(nx, ny, !drawdata->gray, drawdata->p0, drawdata->p, drawdata->zlim);
	cairo_surface_destroy(drawdata->image);
	drawdata->image= cairo_image_surface_create_for_data 
	    (drawdata->p, drawdata->format, nx, ny, stride);
    }
    drawdata->widthim_last=drawdata->widthim;
    drawdata->heightim_last=drawdata->heightim;

    /*Offset in the cairo surface to draw the image. */
    double xoff=round(((width-widthim-SP_XL-sp_xr)*0.5)+SP_XL);
    double yoff=round(((height-heightim-SP_YT-SP_YB)*0.5)+SP_YT);
    drawdata->xoff=xoff;/*save for the GUI to use. */
    drawdata->yoff=yoff;
    /*center of the image on the screen. */
    drawdata->centerx=xoff+widthim*0.5;
    drawdata->centery=yoff+heightim*0.5;
    double zoomx=drawdata->zoomx;/*Zoom of the image when displayed. */
    double zoomy=drawdata->zoomy;
    cairo_select_font_face(cr, font_name, font_style, font_weight);
    cairo_set_font_size(cr, font_size);
    double linewidth=round(font_size*0.08);
    double ticlength=font_size*0.5;
    double ticskip=drawdata->ticinside?5:ticlength;
    double gridskip=drawdata->ticinside?ticlength:0;
    if(drawdata->ticinside){
	ticlength=-ticlength;
    }
    if(linewidth<1) linewidth=1;
    cairo_set_line_width(cr,linewidth);
    /*Save the state of cairo before we drawing the image/points. */
    cairo_save(cr);
    /*flip upside down so that lower left is (0,0); */
    cairo_translate(cr,xoff, heightim+yoff);
    cairo_scale(cr,1,-1);
    /*clip out an rectangular region to draw. */
    cairo_rectangle(cr,0,0,widthim,heightim);
    cairo_clip(cr);
    if(drawdata->image){
	cairo_save(cr);
	cairo_scale(cr, scalex*zoomx, scaley*zoomy);
	/*
	  offx, offy are the offset in the cairo window.
	  ofx, ofy are the actual offset in the original data to display.
	 */
	double ofx=(drawdata->nx*0.5)*(1./zoomx-1.)+drawdata->offx/scalex;
	double ofy=(drawdata->ny*0.5)*(1./zoomy-1.)+drawdata->offy/scaley;
	/*The x and y patterns are negated and then set as
	  translation values in the pattern matrix.*/
	cairo_set_source_surface(cr, drawdata->image, ofx,ofy);
	if(scalex*zoomx>1){/*use nearest filter for up sampling to get clear images */
	    cairo_pattern_set_filter(cairo_get_source(cr),CAIRO_FILTER_NEAREST);
	}
	cairo_paint(cr);
	/*cairo_reset_clip(cr); */
	double xdiff=xmax-xmin; 
	double ydiff=ymax-ymin;
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
	//toc2("cairo_draw image");
    }
    if(drawdata->npts>0){
	cairo_save(cr);
	double centerx=(xmax+xmin)*0.5;
	double centery=(ymax+ymin)*0.5;
	double ncx=widthim*0.5 + drawdata->offx*zoomx;
	double ncy=heightim*0.5 + drawdata->offy*zoomy;
	cairo_set_antialias(cr,CAIRO_ANTIALIAS_NONE);
	cairo_pattern_set_filter(cairo_get_source(cr),CAIRO_FILTER_NEAREST);
#if DRAW_NEW == 1
	int new_width0=(int)ceil(widthim*zoomx)+font_size*80;
	int new_height0=(int)ceil(heightim*zoomy)+font_size*80;
	int new_height=new_height0;
	int new_width=new_width0;
	/*Don't use too much buffer. Over flow memory. */
#define ZOOM_MAX 10
	if(zoomx>ZOOM_MAX){
	    new_width=widthim*ZOOM_MAX;
	}else{
	    drawdata->ncxoff=0;
	}
	if(zoomy>ZOOM_MAX){
	    new_height=heightim*ZOOM_MAX;
	}else{
	    drawdata->ncyoff=0;
	}

	int new_offx=(int)round((widthim-new_width)*0.5+drawdata->offx*zoomx+ drawdata->ncxoff);
	int new_offy=(int)round((heightim-new_height)*0.5+drawdata->offy*zoomy+drawdata->ncyoff);
	if(drawdata->cacheplot && (cairo_image_surface_get_width(drawdata->cacheplot)!=new_width ||
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
	if(new_width < new_width0 && (new_offx>0 || -new_offx>(new_width-widthim))){
	    drawdata->ncxoff=(int)( - drawdata->offx*zoomx);
	    new_offx=(int)round((widthim-new_width)*0.5+drawdata->offx*zoomx+ drawdata->ncxoff);
	    redraw=1;
	}
	if(new_height<new_height0 && (new_offy>0 || -new_offy>(new_height-heightim))){
	    drawdata->ncyoff=(int)( - drawdata->offy*zoomy);
	    new_offy=(int)round((heightim-new_height)*0.5+drawdata->offy*zoomy+drawdata->ncyoff);
	    redraw=1;
	}
	/*center of plot. Important. */
	ncx=(double)new_width*0.5 - drawdata->ncxoff;
	ncy=(double)new_height*0.5 - drawdata->ncyoff;
	if(redraw){
	    cairo_t *cr2=cr;
	    cr=cairo_create(drawdata->cacheplot);
	    cairo_set_antialias(cr,CAIRO_ANTIALIAS_NONE);
	    cairo_pattern_set_filter(cairo_get_source(cr),CAIRO_FILTER_BILINEAR);
	    cairo_set_line_width(cr,linewidth);
	    /*Blank it first. */
	    cairo_set_source_rgba(cr, 1.,1.,1.,1.);
	    cairo_rectangle(cr, 0, 0, new_width, new_height);
	    cairo_fill(cr);
#endif
	if(!drawdata->style_pts){
	    drawdata->style_pts=mycalloc(drawdata->npts,int);/*save the styles for legend */
	}
	int color=0x0000FF;
	cairo_set_source_rgba(cr,0.2,0.0,1.0,1.0);/*Blue */
	int style;
	int connectpts;
	if(drawdata->square){/*we are plotting points. */
	    style=3;/*default style is + */
	    connectpts=0;
	}else{/*we are plotting curves */
	    style=0;
	    connectpts=1;
	}
	double sym_size=3;//round(3*sqrt(zoomx));
	if(drawdata->nstyle==1){
	    PARSE_STYLE(drawdata->style[0]);
	}
	/*computed from below ix, iy formula by setting ix, iy to 0 and widthim or heightim */
	int icumu=(int)drawdata->icumu;
	for(int ipts=0; ipts<drawdata->npts; ipts++){
	    cairo_save(cr);
	    cairo_set_antialias(cr,CAIRO_ANTIALIAS_NONE);
	    if(drawdata->nstyle>1){
		PARSE_STYLE(drawdata->style[ipts]);
	    }else if(drawdata->nstyle==0){
		color=default_color(ipts);
	    }
	    /*save the styles for legend */
	    drawdata->style_pts[ipts]=style | connectpts<<3 |color << 8 |(int)(sym_size)<<4;
	    set_color(cr, color);

	    const double *pts=drawdata->pts[ipts];
	    const int ptsnx=drawdata->ptsdim[ipts][0];
	    const int ptsny=drawdata->ptsdim[ipts][1];
	    if(!pts || ptsnx==0 || ptsny==0) continue;
	    const double *ptsx, *ptsy=NULL;
	    if(ptsny==2){
		ptsx=pts;
		ptsy=ptsx+ptsnx;
	    }else{
		ptsx=0;
		ptsy=pts;
	    }
	    int ips0=0;
	    if(drawdata->cumu && icumu<ptsnx){
		ips0=icumu;
	    }
	    int ptstep=1;
	    double ix=0, iy=0, y=0, y_cumu=0;
	    if(connectpts){/*plot curves. */
		int ips=ips0;
		/*Draw first point */
		if(ptsx){/*don't do round here. */
		    ix=(((xlog?log10(ptsx[ips]):ptsx[ips])-centerx)*scalex*zoomx+ncx);
		}else{
		    ix=(((xlog?log10(ips):ips)-centerx)*scalex*zoomx+ncx);
		}
		if(!isfinite(ptsy[ips])){
		    continue;
		}
		iy=(((ylog?log10(ptsy[ips]):ptsy[ips])-centery)*scaley*zoomy+ncy);
		if(drawdata->cumuquad){
		    y_cumu=ptsy[ips]*ptsy[ips];
		}else{
		    y_cumu=ptsy[ips];
		}
		cairo_move_to(cr, ix, iy);
		/*connect additional points. */
		for(ips++; ips<ptsnx; ips+=ptstep){
		    if(ptsx){/*don't do round here. */
			ix=(((xlog?log10(ptsx[ips]):ptsx[ips])-centerx)*scalex*zoomx+ncx);
		    }else{
			ix=(((xlog?log10(ips):ips)-centerx)*scalex*zoomx+ncx);
		    }
		    if(!isfinite(ptsy[ips])){
			break;
		    }
		    if(drawdata->cumu){
			if(ptsy[ips]!=0){
			    if(drawdata->cumuquad){
				y_cumu+=ptsy[ips]*ptsy[ips];
				y=sqrt(y_cumu/(ips-ips0+1));
			    }else{
				y_cumu+=ptsy[ips];
				y=y_cumu/(ips-ips0+1);
			    }
			}
		    }else{
			y=ptsy[ips];
		    }
		    iy=(((ylog?log10(y):y)-centery)*scaley*zoomy+ncy);
	
		    cairo_line_to(cr, round(ix), round(iy));
		    if(ips%100==0){
			cairo_stroke(cr);
			cairo_move_to(cr, round(ix), round(iy));
		    }
		}
		cairo_stroke(cr);
	    }
	    if(!(connectpts && style==5)){/*plot points. */
		y_cumu=0;
		for(int ips=ips0; ips<ptsnx; ips+=ptstep){
		    if(!isfinite(ptsy[ips])) break;
		    /*Map the coordinate to the image */
		    if(ptsx){/*don't do round here. */
			ix=(((xlog?log10(ptsx[ips]):ptsx[ips])-centerx)*scalex*zoomx+ncx);
		    }else{
			ix=(((xlog?log10(ips):ips)-centerx)*scalex*zoomx+ncx);
		    }
		    
		    if(drawdata->cumu){
			if(ptsy[ips]!=0){
			    if(drawdata->cumuquad){
				y_cumu+=ptsy[ips]*ptsy[ips];
				y=sqrt(y_cumu/(ips-ips0+1));
			    }else{
				y_cumu+=ptsy[ips];
				y=y_cumu/(ips-ips0+1);
			    }
			}
		    }else{
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
	    if(ptsnx>0 && drawdata->legendcurve){
		cairo_save(cr);
		cairo_translate(cr, ix+2, round((iy-font_size*0.5)/1)*1);
		cairo_scale(cr,1,-1);
		if(drawdata->cumu){
		    char val[80];
		    snprintf(val, 80, "%.2f", y);
		    pango_text(cr, layout, 0, 0, val, 1, 1, 0);
		}
		if(drawdata->legend && connectpts){
		    pango_text(cr, layout, 0, 0, drawdata->legend[ipts], 0, 1, 0);
		}
		cairo_restore(cr);
	    }
	    cairo_restore(cr);
	}/*iptsy */
#if DRAW_NEW == 1
	cairo_destroy(cr);
	cr=cr2;
	}
	cairo_set_antialias(cr,CAIRO_ANTIALIAS_NONE);
	/*cairo_pattern_set_filter(cairo_get_source(cr),CAIRO_FILTER_NEAREST); */
	cairo_set_source_surface(cr, drawdata->cacheplot, new_offx, new_offy);
	cairo_paint(cr);
#endif

	xmax0=(((widthim)*0.5)/zoomx - drawdata->offx)/scalex+centerx;
	xmin0=((-widthim*0.5)/zoomx - drawdata->offx)/scalex+centerx;
	ymax0=(((heightim)*0.5)/zoomy - drawdata->offy)/scaley+centery;
	ymin0=((-heightim*0.5)/zoomy - drawdata->offy)/scaley+centery;

	cairo_restore(cr);
	//toc2("cairo_draw pts");
    }

    if(drawdata->ncir>0){
	cairo_save(cr);
	cairo_set_antialias(cr,CAIRO_ANTIALIAS_GRAY);

	double centerx=(xmax+xmin)/2-drawdata->offx/scalex;
	double centery=(ymax+ymin)/2-drawdata->offy/scaley;
	int ncx=widthim/2;
	int ncy=heightim/2;

	for(int icir=0; icir<drawdata->ncir; icir++){
	    int color=(int)drawdata->cir[icir][3];
	    set_color(cr, color);
	    double ix=(drawdata->cir[icir][0]-centerx)*scalex*zoomx+ncx;
	    double iy=(drawdata->cir[icir][1]-centery)*scaley*zoomy+ncy;
	    cairo_arc(cr,ix,iy, drawdata->cir[icir][2]*scalex*zoomx, 0,M_PI*2);
	    cairo_stroke(cr);
	}
	cairo_restore(cr);
    }
    cairo_restore(cr);
    cairo_set_source_rgba(cr, 0.0, 0.0, 0.0,1.0);
    /*the border */
    cairo_rectangle(cr,xoff,yoff,widthim,heightim);
    cairo_stroke(cr);/*border */
    /*Now doing the tic, and colobar */
    /*When there is no zoom, panning, limit0 equals to limit. */
    drawdata->limit0[0]=xmin0;
    drawdata->limit0[1]=xmax0;
    drawdata->limit0[2]=ymin0;
    drawdata->limit0[3]=ymax0;
    if(drawdata->spins){/*dialog is running, update its values */
	for(int i=0; i<4; i++){/*update spin button's value. */
	    gtk_spin_button_set_value(GTK_SPIN_BUTTON(drawdata->spins[i]), drawdata->limit0[i]);
	}
	if(drawdata->zlim){
	    for(int i=5; i<6; i++){/*update spin button's value. */
		gtk_spin_button_set_value(GTK_SPIN_BUTTON(drawdata->spins[i]), drawdata->zlim[i-4]);
	    }
	}
    }
    char ticval[80];
    double tic1, dtic;
    int ntic, order;
    double sep;
    calc_tic(&tic1,&dtic,&ntic,&order,xmax0,xmin0,maxtic_x, xlog);
    sep=xmax0-xmin0;
    int nXticShown=0;

    for(int itic=0; itic<ntic; itic++){
	double ticv=(tic1+dtic*itic); if(fabs(ticv)<1e-12) ticv=0;
	double val=ticv*pow(10,order);
	double frac=(val-xmin0)/sep;
	double xpos=xoff+widthim*frac;
	/*draw the tic */
	if(frac>0.0 && frac<1){
	    cairo_move_to(cr,xpos,yoff+heightim);
	    cairo_line_to(cr,xpos,yoff+heightim+ticlength);
	    cairo_stroke(cr);
	    /*draw the grid */
	    if(drawdata->grid){
		cairo_save(cr);
		cairo_set_dash(cr, stroke_dot, 2, 0);
		cairo_move_to(cr, xpos, yoff);
		cairo_line_to(cr, xpos, yoff+heightim);
		cairo_stroke(cr);
		cairo_restore(cr);
	    }
	}
	if(frac>=0 && frac<=1){
	    snprintf(ticval,80,"%g",(xlog?pow(10,ticv):ticv));
	    pango_text(cr, layout, xpos, yoff+heightim+font_size*0.6+ticskip+1,ticval, 0.5, 0.5, 0);
	    nXticShown++;
	}
    }

    {//draw the minor ticks

	for(int itic=-1; itic<ntic; itic++){
	    double ticv=(tic1+dtic*itic);
	    for(int im=1; im<10; im++){
		double ticvm=ticv+(xlog?log10(1+im):im*0.1)*dtic;
		double valm=ticvm*pow(10,order);
		double fracm=(valm-xmin0)/sep;
		double xposm=xoff+widthim*fracm;
		if((fracm)>0. && (fracm)<1){
		    cairo_move_to(cr,xposm,yoff+heightim);
		    cairo_line_to(cr,xposm,yoff+heightim+ticlength/2);
		    cairo_stroke(cr);
		    if(drawdata->grid && nXticShown<3){
			cairo_save(cr);
			cairo_set_dash(cr, stroke_dot, 2, 0);
			cairo_move_to(cr, xposm, yoff);
			cairo_line_to(cr, xposm, yoff+heightim);
			cairo_stroke(cr);
			cairo_restore(cr);
		    }
		}
		if((fracm)>=0. && (fracm)<=1 && nXticShown<3){//shown minor tic values
		    snprintf(ticval,80,"%g",(xlog?pow(10,ticvm):ticvm));
		    pango_text(cr, layout, xposm, yoff+heightim+font_size*0.6+ticskip+1,ticval, 0.5, 0.5, 0);
		}
	    }
	}
    }
    pango_text_powindex(cr, layout, xoff+widthim-font_size*2, yoff+heightim+6+font_size*1.2,order, 0);
    calc_tic(&tic1,&dtic,&ntic,&order,ymax0,ymin0,maxtic_y, ylog);
    sep=ymax0-ymin0;
    for(int itic=0; itic<ntic; itic++){
	double ticv=(tic1+dtic*itic); if(fabs(ticv)<1e-12) ticv=0;
	double val=ticv*pow(10,order);
	double frac=(val-ymin0)/sep;
	double ypos=yoff+heightim*(1-frac);
	/*draw the tic */
	if(frac>0. && frac<1){
	    cairo_move_to(cr,xoff,ypos);
	    cairo_line_to(cr,xoff-ticlength,ypos);
	    cairo_stroke(cr);
	    /*draw the grid */
	    if(drawdata->grid){
		cairo_save(cr);
		cairo_set_dash(cr, stroke_dot, 2, gridskip);
		cairo_move_to(cr,xoff,ypos);
		cairo_line_to(cr,xoff+widthim,ypos);
		cairo_stroke(cr);
		cairo_restore(cr);
	    }
	}
	if(frac>=0 && frac<=1){
	    snprintf(ticval,80,"%g",(ylog?pow(10,ticv):ticv));
	    pango_text(cr,layout,xoff-font_size*0.6-ticskip+1, ypos,ticval,0.5,0.5,1);
	}
    }
    {//draw minor ticks
	cairo_save(cr);
	for(int itic=-1; itic<ntic; itic++){
	    double ticv=(tic1+dtic*itic);
	    for(int im=(itic==0?-9:1); im<10; im++){
		double ticvm=ticv+(ylog?log10(im+1):im*0.1)*dtic;
		double valm=ticvm*pow(10,order);
		double fracm=(valm-ymin0)/sep;
		double yposm=yoff+heightim*(1-fracm);
		if((fracm)>0. && (fracm)<1){
		    cairo_move_to(cr,xoff,yposm);
		    cairo_line_to(cr,xoff-ticlength/2,yposm);
		}
	    }
	}
	cairo_stroke(cr);
	cairo_restore(cr);
    }
    pango_text_powindex(cr,layout,xoff-font_size*2.8, yoff+font_size*1.8,order, 1);

    if(drawdata->zlim){/*draw colorbar */
	cairo_save(cr);
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
	cairo_stroke(cr);

	cairo_set_source_rgba(cr, 0.0, 0.0, 0.0,1.0);
	calc_tic(&tic1,&dtic,&ntic,&order, drawdata->zlim[1],drawdata->zlim[0],maxtic_y, 0);
	sep=drawdata->zlim[1]-drawdata->zlim[0];

	for(int itic=0; itic<ntic; itic++){
	    double ticv=tic1+dtic*itic; if(fabs(ticv)<1e-12) ticv=0;
	    double val=ticv*pow(10,order);
	    double frac;
	    if(sep>1.e-10*fabs(drawdata->zlim[1])){
		frac=(val-drawdata->zlim[0])/sep;
	    }else{
		if(itic==0) frac=0;
		else if(itic==1) frac=1;
		else frac=-1;
	    }
	    if(frac>0 && frac<1){
		cairo_move_to(cr,0,heightim*(1-frac));
		cairo_line_to(cr,4,heightim*(1-frac));
		cairo_move_to(cr,LEN_LEG,heightim*(1-frac));
		cairo_line_to(cr,LEN_LEG-4,heightim*(1-frac));
		cairo_stroke(cr);
	    }
	    if(frac>=0 && frac<=1){
		snprintf(ticval,80,"%g",ticv);
		pango_text(cr,layout,LEN_LEG+4,heightim*(1-frac), ticval,0,0,0);
	    }
	}
	pango_text_powindex(cr,layout,LEN_LEG/2,-font_size*1.4-2,order,0);
	cairo_restore(cr);
    }
    if(drawdata->title){
	pango_text(cr,layout,xoff+widthim/2,yoff-font_size*0.5-4,drawdata->title,0.5,0.5,0);
    }
    if(drawdata->xlabel){
	pango_text(cr,layout,xoff+widthim/2,yoff+heightim+8+font_size*1.8, drawdata->xlabel,0.5,0.5,0);
    }
    if(drawdata->ylabel){
	pango_text(cr,layout,xoff-font_size*1.8-6,yoff+heightim/2, drawdata->ylabel,0.5,0.5,1);
    }
    if(drawdata->legend && drawdata->npts && drawdata->legendbox){
	int style, color, connectpts;
	double sym_size=0;
	cairo_save(cr);
	cairo_identity_matrix(cr);
	/*draw legend */
	char **legend=drawdata->legend;
	const int ng=drawdata->npts;
	cairo_text_extents_t extents;
	const double linelen=30;/*length of line in legend if exist. */
	double maxlen=0;/*maximum legend length. */
	double tall=0;
	double leglen=0;/*length of legend symbol */
	/*first figure out the size required of longest legend entry. */
	for(int ig=0; ig<ng; ig++){
	    cairo_text_extents(cr, legend[ig], &extents);
	    maxlen=MAX(maxlen, extents.width+1);/*length of text. */
	    tall=MAX(tall, extents.height+1);/*tall of text. */
	    PARSE_STYLE(drawdata->style_pts[ig]);
	    if(connectpts){
		leglen=MAX(leglen, linelen);
	    }else{
		leglen=MAX(leglen, sym_size*2);
		tall=MAX(tall, sym_size*2);
	    }
	}
	const double legmarin=3;/*margin inside of box */
	const double legmarout=5;/*margin outside of box */
	const double linehead=3;/*space before and after symbol */
	double legwidth=maxlen+leglen+2*legmarin+linehead*2;
	double legheight=tall*ng+legmarin*2;
	cairo_translate(cr, xoff + legmarout + drawdata->legendoffx*(widthim - legwidth - 2*legmarout), 
			yoff + legmarout + drawdata->legendoffy*(heightim -legheight - 2*legmarout));
	if(1){//box
	    cairo_save(cr);
	    cairo_rectangle(cr, 0, 0, legwidth, legheight);
	    cairo_set_source_rgba(cr, 1.0, 1.0, 1.0,1.0);
	    cairo_fill_preserve(cr);
	    cairo_set_source_rgba(cr, 0.0, 0.0, 0.0,1.0);
	    cairo_stroke(cr);
	    cairo_restore(cr);
	}
	cairo_translate(cr, legmarin+linehead, legmarin);
	for(int ig=0; ig<ng; ig++){
	    PARSE_STYLE(drawdata->style_pts[ig]);
	    set_color(cr, color);
	    double ix=leglen*0.5;
	    double iy=tall*0.5;
	    draw_point(cr, ix, iy, style, sym_size);
	    cairo_stroke(cr);
	    if(connectpts){
		cairo_move_to(cr, 0, tall*0.5);
		cairo_line_to(cr, leglen, tall*0.5);
		cairo_stroke(cr);
	    }
	    cairo_set_source_rgba(cr, 0.0, 0.0, 0.0,1.0);
	    cairo_move_to(cr, leglen+linehead, tall*0.8);
	    /*cairo_show_text(cr, legend[ig]); */
	    pango_text(cr, layout, leglen+linehead, -tall*0.2, legend[ig], 0, 0,0);
	    cairo_translate(cr, 0, tall);
	}
	cairo_restore(cr);
    }
    if(drawdata->dtime<10){
	cairo_save(cr);
	cairo_identity_matrix(cr);
	char fps[10];
	snprintf(fps, 10, "FPS:%5.1f", 1./drawdata->dtime);
	pango_text(cr, layout, font_size*0.3, font_size*0.2, fps, 0, 0,0);
	cairo_restore(cr);
    }
    g_object_unref(layout);
    cairo_destroy(cr);
    drawdata->icumulast=drawdata->icumu;
    drawdata->cumulast=drawdata->cumu;
    drawdata->cumuquadlast=drawdata->cumuquad;
    drawdata->zoomxlast=drawdata->zoomx;
    drawdata->zoomylast=drawdata->zoomy;
    drawdata->drawn=1;
    //toc2("cairo_draw");
}

#include "drawdaemon.h"
/*
  Routines in this file are about drawing in the cairo surface.
*/
const double stroke_dot[2]={1,5};
const double stroke_dash[2]={10,10};
const double stroke_solid[2]={10,0};


double SP_XL;//space reserved for ylabel
double SP_YT;//space reserved for title
double SP_YB;//space reserved for xlabel
double SP_XR;//space reserved for legend

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
   Draw some text starting at location (x,y)
*/
static void cairo_text_left(cairo_t *cr,double x,double y,const char *text){
    cairo_text_extents_t extents;
    cairo_text_extents (cr, text, &extents);
    double x2 = x - (extents.x_bearing);
    double y2 = y - (extents.height/2 + extents.y_bearing);
    cairo_move_to (cr, x2, y2);
    cairo_show_text (cr, text);
}
/**
   Draw some text centered at location (x,y)
*/
static void cairo_text_center(cairo_t *cr,double x,double y,const char *text){
    cairo_text_extents_t extents;
    cairo_text_extents (cr, text, &extents);
    double x2 = (x-(extents.width/2 + extents.x_bearing));
    double y2 = (y-(extents.height/2 + extents.y_bearing));
    cairo_move_to (cr, x2, y2);
    cairo_show_text (cr, text);
}
/**
   Draw some text vertically, centered at location (x,y)
*/
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
/**
   Draw a cross from the location (x,y.
 */
static void cairo_text_cross(cairo_t *cr, double x, double y, double size){
    cairo_save(cr);
    size=size*0.5;
    cairo_translate(cr,x,y);
    cairo_move_to(cr,0,0);
    cairo_line_to(cr,-size, -size);
    cairo_move_to(cr,0,0);
    cairo_line_to(cr,-size, size);
    cairo_move_to(cr,0,0);
    cairo_line_to(cr,size, -size);
    cairo_move_to(cr,0,0);
    cairo_line_to(cr,size, size);
    cairo_stroke(cr);
    cairo_restore(cr);
}
static void cairo_text_powindex(cairo_t *cr, double rot,
				double x, double y, int order){
    if(order==0) return;
    cairo_save(cr);
    char ticval[80];
    double size=font_size*0.6;//size of cross
    snprintf(ticval,80,"10");
    cairo_text_extents_t extents;
    cairo_text_extents (cr, ticval, &extents); 
    //draw x
    cairo_translate(cr,x,y);
    if(fabs(rot)>0){
	cairo_rotate(cr,rot);
    }
    cairo_set_line_width(cr, font_size*0.08);
    cairo_set_line_cap(cr,CAIRO_LINE_CAP_SQUARE);
    cairo_text_cross(cr, size*0.5, -size*0.5, size*0.9);

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
    //first get the order of magnitude.
    double rmax=fabs(xmax);
    if(fabs(xmin)>rmax) rmax=fabs(xmin);
    int order1=(int)floor(log10(rmax));

    if(fabs(diff)<fabs(1.e-4*xmax)){//very small separation
	xmax=xmax/pow(10,order1);
	*order=order1;
	*ntic=2;
	*dtic=0;
	*tic1=xmax;
    }else{
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
/**
   adjust xmin, xmax properly.
*/
void round_limit(double *xmin, double *xmax){
    if(fabs(*xmin)<1e-15 && fabs(*xmax)<1e-15){//both are zero.
	*xmin=-1;
	*xmax=1;
    }else{
	double tic1, dtic;
	int ntic, order;
	calc_tic(&tic1, &dtic, &ntic, &order, *xmax, *xmin);
	double xmin0=tic1*pow(10,order);
	double xmax0=(tic1+dtic*(ntic-1))*pow(10,order);
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
    }
}

/**
   Create a color based on index.
*/
static int default_color(int ind){
    //we only have 8 colors
    static int *color=NULL;
    if(!color){
	color=calloc(8, sizeof(int));
	color[0]=0x009;
	color[1]=0x900;
	color[2]=0x090;
	color[3]=0x099;
	color[4]=0x909;
	color[5]=0x990;
	color[6]=0x999;
	color[7]=0x449;
    }
    return color[ind & 7];
}
/**
   convert new limit to zoom and off
*/
void apply_limit(drawdata_t *drawdata){
    //limit0 matches limit in unzoomed state
    double diffx0=(drawdata->limit[1]-drawdata->limit[0]);
    double diffy0=(drawdata->limit[3]-drawdata->limit[2]);
    double midx0=(drawdata->limit[0]+drawdata->limit[1])*0.5;
    double midy0=(drawdata->limit[2]+drawdata->limit[3])*0.5;

    double diffx1=drawdata->limit0[1]-drawdata->limit0[0];
    double diffy1=drawdata->limit0[3]-drawdata->limit0[2];
    double midx1=(drawdata->limit0[1]+drawdata->limit0[0])*0.5;
    double midy1=(drawdata->limit0[3]+drawdata->limit0[2])*0.5;

    if(diffx1 > diffx0*1e-3 && diffy1 > diffy0 *1e-3){//limit allowable range
	//the new zoom
	double ratiox=diffx0/diffx1;
	double ratioy=diffy0/diffy1;
	if(drawdata->square){//make the ratio equal.
	    if(fabs(ratiox-drawdata->zoomx)>1e-2 && fabs(ratioy-drawdata->zoomy)<1e-2){
		//only x changed
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
}
#define PARSE_STYLE(stylein)					\
    {								\
	style=stylein & 0x7;					\
	connectpts=(stylein&0x8)>>3;				\
	color=(stylein&0xFFFFFF00)>>8;				\
	size=round(((stylein&0xF0)>>4) * sqrt(zoomx));		\
	if(style>5) style=0;					\
	if(style==0) connectpts=1;				\
    }
static void inline
draw_point(cairo_t *cr, double ix, double iy, long style, double size){
    double size1=size+1;
    switch(style){
    case 0://nothing. just connect lines
	break;
    case 1:// o
	cairo_new_sub_path(cr);
	cairo_arc(cr, ix-0.5, iy-0.5, size, 0, 2*M_PI);
	cairo_move_to(cr, ix, iy);
	break;
    case 2:// x
	cairo_move_to(cr,ix-size1,iy-size1);
	cairo_line_to(cr,ix+size,iy+size);
	cairo_move_to(cr,ix+size,iy-size1);
	cairo_line_to(cr,ix-size1,iy+size);
	cairo_move_to(cr, ix, iy);
	break;
    case 3:// +
	//-1 is because flipping makes effective rounding of 0.5 differently.
	cairo_move_to(cr,ix-size1,iy-1);
	cairo_line_to(cr,ix+size,iy-1);
	cairo_move_to(cr,ix,iy-size1);
	cairo_line_to(cr,ix,iy+size);
	cairo_move_to(cr, ix, iy);
	break;
    case 4:// square []
	cairo_move_to(cr,ix-size1,iy-size1);
	cairo_line_to(cr,ix+size,iy-size1);
	cairo_line_to(cr,ix+size,iy+size);
	cairo_move_to(cr,ix+size,iy+size-1);
	cairo_line_to(cr,ix-size1,iy+size-1);
	cairo_move_to(cr,ix-size,iy+size-1);
	cairo_line_to(cr,ix-size,iy-size1);
	cairo_move_to(cr, ix, iy);
	break;
    case 5:// .
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
   The master routine that draws in the cairo surface.
*/
void cairo_draw(cairo_t *cr, drawdata_t *drawdata, int width, int height){
    //fill white background
    drawdata->font_name_version=font_name_version;
    cairo_rectangle(cr,0,0,width,height);
    cairo_set_source_rgb(cr,1,1,1);
    cairo_fill(cr);
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
    double xdim, ydim;//dimension of the data.
    if(drawdata->image){//we are drawing an image.
	xdim=(double)drawdata->nx;
	ydim=(double)drawdata->ny;
    }else{
	xdim=xmax-xmin;
	ydim=ymax-ymin;
    }
    double sp_xr=20;
    if(drawdata->image){//there is a colorbar
	sp_xr=SP_XR;
    }
    double scalex = (double)(width-SP_XL-sp_xr)/xdim;
    double scaley = (double)(height-SP_YT-SP_YB)/ydim;
    //info("width=%d, height=%d\n", width, height);
    //info("xmax=%g, xmin=%g, xdim=%g, scalex=%g, scaley=%g\n", xmax, xmin, xdim, scalex, scaley);

    if(drawdata->square){
	double scale  = (scalex<scaley?scalex:scaley);
	scalex = scaley = scale;
    }
    widthim=(int)(xdim*scalex);
    heightim=(int)(ydim*scaley);
    scalex = widthim/xdim;
    scaley = heightim/ydim;
    drawdata->widthim=widthim;
    drawdata->heightim=heightim;
    drawdata->scalex=scalex;
    drawdata->scaley=scaley;
    drawdata->drawn=1;
    if(drawdata->limit_changed
       || (drawdata->widthim_last !=0 
	   && (drawdata->widthim_last!=drawdata->widthim 
	       || drawdata->heightim_last!=drawdata->heightim))){
	//canvas is resized, need to adjust zoom/paning
	apply_limit(drawdata);
	drawdata->limit_changed=0;
    }
    drawdata->widthim_last=drawdata->widthim;
    drawdata->heightim_last=drawdata->heightim;

    //Offset in the cairo surface to draw the image.
    double xoff=round(((width-widthim-SP_XL-sp_xr)*0.5)+SP_XL);
    double yoff=round(((height-heightim-SP_YT-SP_YB)*0.5)+SP_YT);
    drawdata->xoff=xoff;
    drawdata->yoff=yoff;
    //center of the image on the screen.
    drawdata->centerx=xoff+widthim*0.5;
    drawdata->centery=yoff+heightim*0.5;
    double zoomx=drawdata->zoomx;//Zoom of the image when displayed.
    double zoomy=drawdata->zoomy;
    cairo_select_font_face(cr, font_name, font_style, font_weight);
    cairo_set_font_size(cr, font_size);
    double linewidth=round(font_size*0.08);
    double ticlength=font_size*0.5;
    double ticskip=drawdata->ticinside?5:ticlength;
    double gridskip=drawdata->ticinside?ticlength:0;
    if(linewidth<1) linewidth=1;
    cairo_set_line_width(cr,linewidth);
    //Save the state of cairo before we drawing the image/points.
    cairo_save(cr);
    //flip upside down so that lower left is (0,0);
    cairo_translate(cr,xoff, heightim+yoff);
    cairo_scale(cr,1,-1);
    //clip out an rectangular region to draw.
    cairo_rectangle(cr,0,0,widthim,heightim);

    cairo_set_source_rgba(cr, 0.0, 0.0, 0.0,1.0);
 
    cairo_stroke(cr);//border
    cairo_rectangle(cr,0,0,widthim,heightim);
    cairo_clip(cr);
    if(drawdata->image){
	cairo_save(cr);
	cairo_scale(cr,scalex*zoomx,scaley*zoomy);
	/*
	  offx, offy are the offset in the cairo window.
	  ofx, ofy are the actual offset in the original data to display.
	 */
	double ofx=(drawdata->nx*0.5)*(1/zoomx-1)+drawdata->offx/scalex;
	double ofy=(drawdata->ny*0.5)*(1/zoomy-1)+drawdata->offy/scaley;
	/*The x and y patterns are negated and then set as
	  translation values in the pattern matrix.*/
	cairo_set_source_surface(cr, drawdata->image, ofx,ofy);
	if(scalex*zoomx>1){//use nearest filter for up sampling to get clear images
	    cairo_pattern_set_filter(cairo_get_source(cr),CAIRO_FILTER_NEAREST);
	}
	cairo_paint(cr);
	cairo_reset_clip(cr);
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
    }
    int styles[drawdata->npts];//save the styles for legend
    if(drawdata->npts>0){
	cairo_save(cr);
	cairo_set_antialias(cr,CAIRO_ANTIALIAS_NONE);//GRAY
	int color=0x009;
	cairo_set_source_rgba(cr,0.2,0.0,1.0,1.0);
	int style;
	int connectpts;
	if(drawdata->square){//we are plotting points.
	    style=3;//default style is +
	    connectpts=0;
	}else{//we are plotting curves
	    style=0;
	    connectpts=1;
	}
	double size=round(3*sqrt(zoomx));
	if(drawdata->nstyle==1){
	    PARSE_STYLE(drawdata->style[0]);
	}
	double centerx=(xmax+xmin)/2;
	double centery=(ymax+ymin)/2;
	double ncx=widthim*0.5 + drawdata->offx*zoomx;
	double ncy=heightim*0.5 + drawdata->offy*zoomy;
	//computed from below ix, iy formula by setting ix, iy to 0 and widthim or heightim
	xmax0=(((widthim)*0.5)/zoomx - drawdata->offx)/scalex+centerx;
	xmin0=((-widthim*0.5)/zoomx - drawdata->offx)/scalex+centerx;
	ymax0=(((heightim)*0.5)/zoomy - drawdata->offy)/scaley+centery;
	ymin0=((-heightim*0.5)/zoomy - drawdata->offy)/scaley+centery;
	for(int ipts=0; ipts<drawdata->npts; ipts++){
	    dmat *pts=drawdata->pts[ipts];
	    double *ptsx=NULL, *ptsy=NULL;
	    if(pts->ny==2){
		ptsx=pts->p;	
		ptsy=ptsx+pts->nx;
	    }else{
		ptsy=pts->p;
	    }
	    if(drawdata->nstyle>1){
		PARSE_STYLE(drawdata->style[ipts]);
	    }else if(drawdata->nstyle==0){
		color=default_color(ipts);
	    }
	    //save the styles for legend
	    styles[ipts]=style | connectpts<<3 |color << 8 |(int)(size/sqrt(zoomx))<<4;
	    int r=color/100;
	    int g=(color-r*100)/10;
	    int b=color-r*100-g*10;
	    cairo_set_source_rgba(cr,r*0.11,g*0.11,b*0.11,1.0);
	    if(connectpts && style==5){
		unsigned int ips=0;
		double ix, iy;
		if(ptsx){//don't do round here.
		    ix=((ptsx[ips]-centerx)*scalex*zoomx+ncx);
		}else{
		    ix=((ips-centerx)*scalex*zoomx+ncx);
		}
		iy=((ptsy[ips]-centery)*scaley*zoomy+ncy);
		cairo_move_to(cr, ix, iy);
		for(ips=1; ips<pts->nx; ips++){
		    if(ptsx){//don't do round here.
			ix=((ptsx[ips]-centerx)*scalex*zoomx+ncx);
		    }else{
			ix=((ips-centerx)*scalex*zoomx+ncx);
		    }
		    iy=((ptsy[ips]-centery)*scaley*zoomy+ncy);
		    cairo_line_to(cr, ix, iy);
		}
		cairo_stroke(cr);
	    }else{
		double ix, iy;
		for(unsigned int ips=0; ips<pts->nx; ips++){
		    //Mape the coordinate to the image
		    if(ptsx){//don't do round here.
			ix=((ptsx[ips]-centerx)*scalex*zoomx+ncx);
		    }else{
			ix=((ips-centerx)*scalex*zoomx+ncx);
		    }
		    iy=((ptsy[ips]-centery)*scaley*zoomy+ncy);
		    if(connectpts && ips>0){
			cairo_line_to(cr, ix, iy);
		    }
		    draw_point(cr, ix, iy, style, size);
		    if(drawdata->nstyle>0){
			cairo_stroke(cr);//stroke each point because color may change.
		    }
		}//ipts
	    }//if
	    if(!drawdata->style){
		cairo_stroke(cr);//stroke all together.
	    }
	}//iptsy
	cairo_restore(cr);
    }

    if(drawdata->ncir>0){
	cairo_save(cr);
	cairo_set_antialias(cr,CAIRO_ANTIALIAS_GRAY);

	double centerx=(xmax+xmin)/2-drawdata->offx/scalex;
	double centery=(ymax+ymin)/2-drawdata->offy/scaley;
	int ncx=widthim/2;
	int ncy=heightim/2;

	for(unsigned int icir=0; icir<drawdata->ncir; icir++){
	    int color=(int)drawdata->cir[icir][3];
	    int r=color/100;
	    int g=(color-r*100)/10;
	    int b=color-r*100-g*10;
	    double ix=(drawdata->cir[icir][0]-centerx)*scalex*zoomx+ncx;
	    double iy=(drawdata->cir[icir][1]-centery)*scaley*zoomy+ncy;
	    cairo_set_source_rgba(cr,r*0.11,g*0.11,b*0.11,1.0);
	    cairo_arc(cr,ix,iy, drawdata->cir[icir][2]*scalex*zoomx,
		      0,M_PI*2);
	    cairo_stroke(cr);
	}
	cairo_restore(cr);
    }
    cairo_restore(cr);
    //Now doing the border, tic, and colobar
    cairo_set_source_rgba(cr, 0.0, 0.0, 0.0,1.0);
    //When there is no zoom, panning, limit0 equals to limit.
    drawdata->limit0[0]=xmin0;
    drawdata->limit0[1]=xmax0;
    drawdata->limit0[2]=ymin0;
    drawdata->limit0[3]=ymax0;
    if(drawdata->spin){//dialog is running, update its values
	for(int i=0; i<4; i++){//update spin button's value.
	    gtk_spin_button_set_value(GTK_SPIN_BUTTON(drawdata->spin[i].w), drawdata->limit0[i]);
	}
    }
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
	//draw the tic
	cairo_move_to(cr,xoff+widthim*frac,yoff+heightim);
	if(fabs(frac)>0.01 && fabs(frac)<0.99){
	    if(drawdata->ticinside){
		cairo_line_to(cr,xoff+widthim*frac,yoff+heightim-ticlength);
	    }else{
		cairo_line_to(cr,xoff+widthim*frac,yoff+heightim+ticlength);
	    }
	    cairo_stroke(cr);
	    //draw the grid
	    if(drawdata->grid){
		cairo_save(cr);
		cairo_set_dash(cr, stroke_dot, 2, 0);
		cairo_move_to(cr,xoff+widthim*frac,yoff);
		cairo_line_to(cr,xoff+widthim*frac,yoff+heightim);
		cairo_stroke(cr);
		cairo_restore(cr);
	    }
	}
	snprintf(ticval,80,"%g",ticv);
	cairo_text_center(cr,xoff+widthim*frac,
			  yoff+heightim+font_size*0.6+ticskip+1,ticval);
    }
    cairo_text_powindex(cr,0, xoff+widthim-font_size*2,
			yoff+heightim+6+font_size*2.2,order);
    calc_tic(&tic1,&dtic,&ntic,&order,ymax0,ymin0);
    sep=ymax0-ymin0;
    for(int itic=0; itic<ntic; itic++){
	double ticv=tic1+dtic*itic;
	double val=ticv*pow(10,order);
	double frac=(val-ymin0)/sep;
	double yh=yoff+heightim*(1-frac);
	//draw the tic
	if(fabs(frac)>0.01 && fabs(frac)<0.99){
	    cairo_move_to(cr,xoff,yh);
	    if(drawdata->ticinside){
		cairo_line_to(cr,xoff+ticlength,yh);
	    }else{
		cairo_line_to(cr,xoff-ticlength,yh);
	    }
	    cairo_stroke(cr);
	    //draw the grid
	    if(drawdata->grid){
		cairo_save(cr);
		cairo_set_dash(cr, stroke_dot, 2, gridskip);
		cairo_move_to(cr,xoff,yh);
		cairo_line_to(cr,xoff+widthim,yh);
		cairo_stroke(cr);
		cairo_restore(cr);
	    }
	}
	snprintf(ticval,80,"%g",ticv);
	cairo_vltext_center(cr,xoff-font_size*0.6-ticskip+1,
			    yoff+heightim*(1-frac),ticval);
    }
    cairo_text_powindex(cr,-M_PI/2,xoff-font_size*1.6,
			yoff+font_size*1.8,order);

    if(drawdata->maxmin){//draw colorbar
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
	cairo_restore(cr);
    }
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
    if(drawdata->legend && drawdata->npts){
	int style, color, connectpts;
	double size;
	cairo_save(cr);
	cairo_identity_matrix(cr);
	//draw legend
	char **legend=drawdata->legend;
	const int ng=drawdata->npts;
	cairo_text_extents_t extents;
	double linelen=50;//length of line in legend.
	double maxlen=0;//maximum legend length.
	double tall=0;
	double leglen=0;
	for(int ig=0; ig<ng; ig++){
	    cairo_text_extents(cr, legend[ig], &extents);
	    maxlen=MAX(maxlen, extents.width);
	    tall=MAX(tall, extents.height);
	    PARSE_STYLE(styles[ig]);
	    if(connectpts){
		leglen=MAX(leglen, linelen);
	    }else{
		leglen=MAX(leglen, size*2);
		tall=MAX(tall, size*2);
	    }
	}
	const double legmarin=2;
	const double legmarout=5;
	const double linehead=2;
	double legwidth=maxlen+leglen+2*legmarin+linehead*2;
	cairo_translate(cr, xoff+widthim - legwidth - legmarout, yoff + legmarout);
	if(drawdata->legendbox){
	    cairo_save(cr);
	    cairo_rectangle(cr, 0, 0, legwidth, tall*ng+legmarin*2);
	    cairo_set_source_rgba(cr, 1.0, 1.0, 1.0,1.0);
	    cairo_fill_preserve(cr);
	    cairo_set_source_rgba(cr, 0.0, 0.0, 0.0,1.0);
	    cairo_stroke(cr);
	    cairo_restore(cr);
	}
	cairo_translate(cr, legmarin+linehead, legmarin);
	for(int ig=0; ig<ng; ig++){
	    PARSE_STYLE(styles[ig]);
	    int r=color/100;
	    int g=(color-r*100)/10;
	    int b=color-r*100-g*10;
	    cairo_set_source_rgba(cr,r*0.11,g*0.11,b*0.11,1.0);
	    double ix=leglen*0.5;
	    double iy=tall*0.5;
	    draw_point(cr, ix, iy, style, size);
	    cairo_stroke(cr);
	    if(connectpts){
		cairo_move_to(cr, 0, tall*0.5);
		cairo_line_to(cr, leglen, tall*0.5);
		cairo_stroke(cr);
	    }
	    cairo_set_source_rgba(cr, 0.0, 0.0, 0.0,1.0);
	    cairo_move_to(cr, leglen+linehead, tall*0.8);
	    cairo_show_text(cr, legend[ig]);
	    cairo_translate(cr, 0, tall);
	}
	cairo_restore(cr);
    }
    cairo_destroy(cr);
}

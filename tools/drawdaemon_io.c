#include "drawdaemon.h"
/*
  Routines in this file handles I/O
*/

char *defname="Title";
char *deffig="Figure";
static int no_longer_listen=0;
int ndrawdata=0;
pthread_mutex_t mutex_drawdata=PTHREAD_MUTEX_INITIALIZER;
#define FILE_READ(data,size)			\
    nleft=size;start=(gchar*)data;		\
    do{						\
	int nread=fread(start, 1, nleft, fp);	\
	nleft-=nread;				\
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
	    pthread_mutex_lock(&mutex_drawdata);
	    ndrawdata++;
	    info("drawdata created, ndrawdata=%d\n", ndrawdata);
	    pthread_mutex_unlock(&mutex_drawdata);
	    drawdata->zoomx=1;
	    drawdata->zoomy=1;
	    drawdata->square=1;//default to square.
	    drawdata->name=defname;
	    drawdata->format=(cairo_format_t)0;
	    drawdata->gray=0;
	    drawdata->ticinside=1;
	    drawdata->legendbox=1;
	    drawdata->fig=NULL;
	    break;
	case FIFO_DATA://image data.
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
	    {
		int nptsx, nptsy;
		int ipts=drawdata->npts;
		drawdata->npts++;
		FILE_READ_INT(nptsx);
		FILE_READ_INT(nptsy);
		FILE_READ_INT(drawdata->square);
		drawdata->grid=1;
		drawdata->pts=realloc(drawdata->pts, drawdata->npts*sizeof(dmat*));
		drawdata->pts[ipts]=dnew(nptsx, nptsy);
		FILE_READ(drawdata->pts[ipts]->p, sizeof(double)*nptsx*nptsy);
	    }
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
	case FIFO_LEGEND:
	    drawdata->legend=calloc(drawdata->npts, sizeof(char*));
	    for(int i=0; i<drawdata->npts; i++){
		FILE_READ_STR(drawdata->legend[i]);
	    }
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
			double xmin0=INFINITY, xmax0=-INFINITY, ymin0=INFINITY, ymax0=-INFINITY;
			for(int ipts=0; ipts<drawdata->npts; ipts++){
			    dmat *pts=drawdata->pts[ipts];
			    double xmin, xmax, ymin, ymax;
			    if(pts->ny>1){
				maxmindbl(pts->p, pts->nx, &xmax, &xmin);
				maxmindbl(pts->p+pts->nx, pts->nx, &ymax, &ymin);
			    }else{
				xmin=0; xmax=(double)(pts->nx-1);
				maxmindbl(pts->p, pts->nx, &ymax, &ymin);
			    }
			    if(xmin<xmin0) xmin0=xmin;
			    if(ymin<ymin0) ymin0=ymin;
			    if(xmax>xmax0) xmax0=xmax;
			    if(ymax>ymax0) ymax0=ymax;
			}
			drawdata->limit[0]=xmin0;
			drawdata->limit[1]=xmax0;
			drawdata->limit[2]=ymin0;
			drawdata->limit[3]=ymax0;
			round_limit(drawdata->limit, drawdata->limit+1);
			round_limit(drawdata->limit+2, drawdata->limit+3);
		    }
		    if(drawdata->nstyle>1){
			if(drawdata->nstyle!=drawdata->npts){
			    warning("nstyle must equal to npts\n");
			    drawdata->nstyle=0;//disable it.
			    free(drawdata->style);
			}
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

void open_fifo(void){
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




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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include "common.h"
#include "thread.h"
#include "misc.h"
#include "mathmisc.h"
int DRAW_ID=0;
int DRAW_DIRECT=0; //set to 1 to launch drawdaemon directly without going through scheduler.
/*
  draw by calling gtkdraw via daemon
*/
/**
   \file draw.c
   Contains functions for data visualization
*/
#include "draw.h"
#include "bin.h"
#include "sys/scheduler_client.h"
static FILE *pfifo=NULL;
#define DOPPRINT 1
#if DOPPRINT == 0
#undef info
#define info(A...)
#endif
#undef FWRITE
#define FWRITE(ptr,size,nmemb,pfifo)		\
    if(fifo_write(ptr,size,nmemb,pfifo)){	\
	goto done;				\
    }
#define FWRITEINT(pfifo,cmd) {int tmp=cmd; if(fifo_write(&tmp,sizeof(int),1,pfifo)){goto done;}}
#define FWRITESTR(pfifo,arg)				\
    if(arg && strlen(arg)>0){				\
	int slen=strlen(arg)+1;				\
	FWRITEINT(pfifo,slen);				\
	if(fifo_write(arg,sizeof(char),slen,pfifo)){	\
	    goto done;					\
	}						\
    }else{						\
	FWRITEINT(pfifo,0);				\
    }
#define FWRITECMDSTR(pfifo,cmd,arg)			\
    if(arg&&strlen(arg)>0){				\
	FWRITEINT(pfifo, cmd); FWRITESTR(pfifo,arg);	\
    }							\

PNEW(lock);
int disable_draw=0; //if 1, draw will be disabled 
/**
   Open a fifo (created by the drawdaemon) and start writing to it.
*/
static int fifo_open(){
    int fd;
    int retry=0;
    char *fifo_fn;
    if(disable_draw){
	return -1;
    }
    if(pfifo){
	return 0;
    }
    if(!DRAW_ID){
	DRAW_ID=getpid();
    }

 retry:
    //do not free the returned string.
    fifo_fn=scheduler_get_drawdaemon(DRAW_ID, DRAW_DIRECT);
    if(!fifo_fn){
	warning("Unable to find and launch drawdaemon\n"); 
	disable_draw=1;
	return -1;
    }
    retry++;
    fd=open(fifo_fn,O_NONBLOCK|O_WRONLY);
    if(fd==-1){
	perror("fifo_open");
	if(errno==ENXIO && retry>20){
	    warning("Unable to open drawdaemon. Draw is disabled\n");
	    disable_draw=1;
	    return -1;
	}else{
	    sleep(2);
	    if(retry<20){
		goto retry;
	    }else{
		return -1;
	    }
	}
    }else{
	//do not use fdopen, which has NONBLOCK flag and is not what we want.
        pfifo=fopen(fifo_fn,"wb");
        close(fd);
	if(!pfifo){
	    warning("Failed to open fifo\n");
	    return -1;
	}else{
	    return 0;
	}
    }
}
/**
   Write data to the fifo. Handle exceptions.
*/
inline static int fifo_write(const void *ptr, /**<Pointer to the data*/
			     size_t size,     /**<Size of each element*/
			     size_t nmemb,    /**<Number of elements*/
			     FILE *fp         /**<Pointer to the File*/
			     ){
 retry:
    if(fwrite(ptr,size,nmemb,fp)!=nmemb){
	perror("fifo_write");
	if(errno==EAGAIN || errno==EWOULDBLOCK){
	    sleep(1);
	    warning("\nfifowrite: Retrying\n");
	    goto retry;
	}else if(errno==EPIPE){
	    warning("\nBroken pipe.\n");
	    fclose(pfifo); pfifo=NULL;
	    return -1;
	}else{
	    warning("\nUnknown error\n");
	    fclose(pfifo); pfifo=NULL;
	    return -1;
	}
    }
    return 0;
}
/**
   Plot the coordinates ptsx, ptsy using style, and optionally plot ncir circles.
*/
void plot_points(char *fig,          /**<Category of the figure*/
		 long ngroup,        /**<Number of groups to plot*/
		 loc_t **loc,        /**<Plot arrays of loc as grid*/
		 dcell *cell,        /**<If loc ismpety, use cell to plot curves*/
		 const int32_t *style,  /**<Style of each point*/
		 const double *limit,/**<x min, xmax, ymin and ymax*/
		 int ncir,           /**<Number of circles*/
		 double (*pcir)[4],  /**<Data for the circles: x, y origin, radius, and color*/
		 char **legend,/**<ngroup number of char**/
		 const char *title,  /**<title of the plot*/
		 const char *xlabel, /**<x axis label*/
		 const char *ylabel, /**<y axis label*/
		 const char *format, /**<subcategory of the plot.*/
		 ...){
    format2fn;
    LOCK(lock);
    if(fifo_open()){//failed to open.
	goto done;
    }
    //info2("pts");getchar();
    FWRITEINT(pfifo, FIFO_START);
    if(loc){//there are points to plot.
	for(int ig=0; ig<ngroup; ig++){
	    FWRITEINT(pfifo, FIFO_POINTS);
	    FWRITEINT(pfifo, loc[ig]->nloc);
	    FWRITEINT(pfifo, 2);
	    FWRITEINT(pfifo, 1);
	    FWRITE(loc[ig]->locx, sizeof(double), loc[ig]->nloc, pfifo);
	    FWRITE(loc[ig]->locy, sizeof(double), loc[ig]->nloc, pfifo);
	}
	if(cell){
	    warning("both loc and cell are specified\n");
	}
    }else if(cell){
	if(ngroup!=cell->nx*cell->ny){
	    warning("ngroup and dimension of cell mismatch\n");
	    ngroup=cell->nx*cell->ny;
	}
	for(int ig=0; ig<ngroup; ig++){
	    FWRITEINT(pfifo, FIFO_POINTS);
	    FWRITEINT(pfifo, cell->p[ig]->nx);
	    FWRITEINT(pfifo, cell->p[ig]->ny);
	    FWRITEINT(pfifo, 0);
	    FWRITE(cell->p[ig]->p, sizeof(double),cell->p[ig]->nx*cell->p[ig]->ny, pfifo);
	}
    }
    if(style){
	FWRITEINT(pfifo,FIFO_STYLE);
	FWRITEINT(pfifo,ngroup);
	FWRITE(style,sizeof(uint32_t),ngroup,pfifo);
    }
    if(ncir>0){
	FWRITEINT(pfifo, FIFO_CIRCLE);
	FWRITEINT(pfifo, ncir);
	FWRITE(pcir, sizeof(double), ncir*4, pfifo);
    }
    if(limit){//xmin,xmax,ymin,ymax
	FWRITEINT(pfifo, FIFO_LIMIT);
	FWRITE(limit, sizeof(double), 4, pfifo);
    }
    if(format){
	FWRITECMDSTR(pfifo,FIFO_NAME,fn);
    }
    if(legend){
	FWRITEINT(pfifo, FIFO_LEGEND);
	for(int ig=0; ig<ngroup; ig++){
	    FWRITESTR(pfifo, legend[ig]);
	}
    }
    FWRITECMDSTR(pfifo,FIFO_FIG,fig);
    FWRITECMDSTR(pfifo,FIFO_TITLE,title);
    FWRITECMDSTR(pfifo,FIFO_XLABEL,xlabel);
    FWRITECMDSTR(pfifo,FIFO_YLABEL,ylabel);
    FWRITEINT(pfifo,FIFO_END);
    if(fflush(pfifo)){
	//it is important to flush it so that drawdaemon does not stuck waiting.
	warning("Failed to fflush the fifo");
    }
 done:
    UNLOCK(lock); 
}
/**
   Draw an image.
*/
void imagesc(char *fig, /**<Category of the figure*/
	     long nx,   /**<the image is of size nx*ny*/
	     long ny,   /**<the image is of size nx*ny*/
	     const double *limit, /**<x min, xmax, ymin and ymax*/
	     const double *zlim,/**< min,max, the data*/
	     const void *p, /**<The image*/
	     const char *title,  /**<title of the plot*/
	     const char *xlabel, /**<x axis label*/
	     const char *ylabel, /**<y axis label*/
	     const char *format, /**<subcategory of the plot.*/
	     ...){
    format2fn;
    LOCK(lock);
    if(fifo_open()){
	goto done;
    }
    int32_t header[2];
    header[0]=nx;
    header[1]=ny;
    FWRITEINT(pfifo, FIFO_START);
    FWRITEINT(pfifo, FIFO_DATA);
    FWRITE(header, sizeof(int32_t), 2, pfifo);
    FWRITE(p, sizeof(double), nx*ny, pfifo);
    if(zlim){
	FWRITEINT(pfifo,FIFO_ZLIM);
	FWRITE(zlim,sizeof(double),2,pfifo);
    }
    if(limit){//xmin,xmax,ymin,ymax
	FWRITEINT(pfifo, FIFO_LIMIT);
	FWRITE(limit, sizeof(double), 4, pfifo);
    }
    
    if(format){
	FWRITECMDSTR(pfifo,FIFO_NAME,fn);
    }
    FWRITECMDSTR(pfifo,FIFO_FIG,fig);
    FWRITECMDSTR(pfifo,FIFO_TITLE,title);
    FWRITECMDSTR(pfifo,FIFO_XLABEL,xlabel);
    FWRITECMDSTR(pfifo,FIFO_YLABEL,ylabel);
    FWRITEINT(pfifo,FIFO_END);
    if(fflush(pfifo)){
	//it is important to flush it so that drawdaemon does not stuck waiting
	warning("Failed to fflush the fifo");
    }
 done:
    UNLOCK(lock);
}

/**
   Draw the OPD of real and imaginary of complex p defined on nx*ny grid. see imagesc()
*/
void imagesc_cmp_ri(char *fig, long nx, long ny, const double *limit, const double *zlim,
		    const dcomplex *p, const char *title, const char *xlabel, const char *ylabel,
		    const char *format,...){
    if(disable_draw) return;
    format2fn;

    double *pr,*pi;
    pr=malloc(nx*ny*sizeof(double));
    pi=malloc(nx*ny*sizeof(double));
    for(int i=0; i<nx*ny; i++){
	pr[i]=creal(p[i]);
	pi[i]=cimag(p[i]);
    }
    imagesc(fig, nx, ny, limit, zlim,pr, title, xlabel, ylabel, "%s real",fn);
    free(pr);
    imagesc(fig, nx, ny, limit, zlim,pi, title, xlabel, ylabel, "%s imag",fn);
    free(pi);
}
/**
   Draw the OPD of abs and phase of complex p defined on nx*ny grid. see imagesc()
*/
void imagesc_cmp_ap(char *fig, long nx, long ny, const double *limit,const double *zlim,
		    const dcomplex *p, const char *title, const char *xlabel, const char *ylabel,
		    const char *format,...){
    if(disable_draw) return;
    format2fn;
    double *pr,*pi;
    int isreal=1;
    pr=malloc(nx*ny*sizeof(double));
    pi=calloc(nx*ny,sizeof(double));
    for(int i=0; i<nx*ny; i++){
	pr[i]=cabs(p[i]);
	if(pr[i]>1.e-10){
	    pi[i]=atan2(cimag(p[i]),creal(p[i]));
	    if(isreal&&fabs(pi[i])>1.e-10) isreal=0;
	}
    }
    imagesc(fig, nx, ny, limit, zlim, pr, title, xlabel, ylabel,
	    "%s abs",fn);
    free(pr);
    if(!isreal){
	imagesc(fig, nx,  ny, limit, zlim, pi,title, xlabel, ylabel,
		"%s phi",fn);
    }
    free(pi);
}
/**
   Draw the OPD of abs of complex p defined on nx*ny grid. see imagesc()
*/
void imagesc_cmp_abs(char *fig, long nx, long ny, const double *limit,const double *zlim,
		     const dcomplex *p, const char *title, const char *xlabel, const char *ylabel,
		     const char *format,...){
    if(disable_draw) return;
    format2fn;
    double *pr;
    pr=malloc(nx*ny*sizeof(double));
    for(int i=0; i<nx*ny; i++){
	pr[i]=cabs(p[i]);
    }
    imagesc(fig, nx, ny, limit, zlim, pr, title, xlabel, ylabel,
	    "%s abs",fn);
    free(pr);
}
#include "dmat.h"
#include "cmat.h"
#include "loc.h"
/*
  The following routines applies the imagesc_* functions onto
  dmat,cmat,loc,map,etc, data types. 
   
  The first argument is the type of the plot. Plots with the
  same types are grouped into a single tab.
   
  The last argument is the title of the plot. Plots with the
  same type and title as existing plots will replace the already
  existing plots.
*/
/**
   Mapping the floating point numbers onto screen with scaling similar to matlab
   imagesc.  . see imagesc()
*/

void ddraw(char *fig, const dmat *A, double *xylim, double *zlim,
	   const char *title, const char *xlabel, const char *ylabel,
	   const char *format,...){
    format2fn;
    imagesc(fig,A->nx,A->ny,xylim,zlim,A->p,title, xlabel, ylabel,"%s",fn);
}

/**
   Mapping the absolution value of complex array. see imagesc()
*/
void cdrawabs(char *fig, const cmat *A, double *xylim, double *zlim,
	      const char *title, const char *xlabel, const char *ylabel,
	      const char *format,...){
  
    format2fn;
    imagesc_cmp_abs(fig,A->nx,A->ny,xylim,zlim,A->p,title,xlabel,ylabel,"%s abs",fn);
}

/**
   Mapping the real/imaginary part of complex array. see imagesc()
*/
void cdrawri(char *fig, const cmat *A, double *xylim, double *zlim,
	     const char *title, const char *xlabel, const char *ylabel,
	     const char *format,...){
    
    format2fn;
    imagesc_cmp_ri(fig,A->nx,A->ny,xylim,zlim,A->p,title, xlabel, ylabel,"%s",fn);
}

/**
   Mapping the absolute and phase of complex array. see imagesc()
*/
void cdraw(char *fig, const cmat *A, double *xylim, double *zlim,
	   const char *title, const char *xlabel, const char *ylabel,
	   const char *format,...){
    format2fn;
    imagesc_cmp_ap(fig,A->nx,A->ny,xylim,zlim,A->p,title, xlabel, ylabel,"%s",fn);
}

/**
   like ddraw, acting on map object. see imagesc()
*/
void drawmap(char *fig, const map_t *map,  double *zlim,
	     const char *title, const char *xlabel, const char *ylabel,
	     const char *format,...){
    format2fn;
    double limit[4];
    limit[0]=map->ox-map->dx/2;
    limit[1]=map->ox+(map->nx-0.5)*map->dx;
    limit[2]=map->oy-map->dx/2;
    limit[3]=map->oy+(map->ny-0.5)*map->dx;
    imagesc(fig, map->nx, map->ny, limit, zlim, map->p,  title, xlabel, ylabel,"%s",fn);
}
/**
   Plot the loc on the screen. see imagesc()
*/
void drawloc(char *fig, loc_t *loc, double *zlim,
	     const char *title, const char *xlabel, const char *ylabel,
	     const char* format,...){
    format2fn;
    loc_create_map_npad(loc,0);
    int npad=loc->map->npad;
    int nxm=loc->map->nx;
    int nx=loc->map->nx-npad*2;
    int ny=loc->map->ny-npad*2;
    double *opd0=calloc(nx*ny, sizeof(double));
    for(int iy=0; iy<ny; iy++){
	for(int ix=0; ix<nx; ix++){
	    opd0[ix+iy*nx]=(loc->map->p[(ix+npad)+(iy+npad)*nxm]>0);
	}
    }
    double limit[4];
    limit[0]=loc->map->ox+loc->dx*(npad-1/2);
    limit[1]=limit[0]+loc->dx*nx;
    limit[2]=loc->map->oy+loc->dx*(npad-1/2);
    limit[3]=limit[2]+loc->dx*ny;
    imagesc(fig, nx, ny,limit,zlim,opd0,  title, xlabel, ylabel,"%s",fn);
    free(opd0);
}

/**
   Plot the opd using coordinate loc. see imagesc()
*/
void drawopd(char *fig, loc_t *loc, const double *opd,  double *zlim,
	     const char *title, const char *xlabel, const char *ylabel,
	     const char* format,...){

    format2fn;
    loc_create_map_npad(loc,0);

    int npad=loc->map->npad;
    int nxm=loc->map->nx;
    int nx=loc->map->nx-npad*2;
    int ny=loc->map->ny-npad*2;
    double *opd0=calloc(nx*ny, sizeof(double));
    for(int iy=0; iy<ny; iy++){
	for(int ix=0; ix<nx; ix++){
	    int ii=loc->map->p[(ix+npad)+(iy+npad)*nxm];
	    if(ii){
		opd0[ix+iy*nx]=opd[ii-1];
	    }else{
		opd0[ix+iy*nx]=NAN;
	    }
	}
    }
    double limit[4];
    limit[0]=loc->map->ox+loc->dx*(npad-1/2);
    limit[1]=loc->map->ox+loc->dx*(nx+npad-1/2);
    limit[2]=loc->map->oy+loc->dx*(npad-1/2);
    limit[3]=loc->map->oy+loc->dx*(ny+npad-1/2);
    imagesc(fig,nx,ny, limit,zlim,opd0,  title, xlabel, ylabel,"%s",fn);
    free(opd0);
}
/**
   Plot opd*amp with coordinate loc. see imagesc()
*/
void drawopdamp(char *fig, loc_t *loc, const double *opd, const double *amp, double *zlim,
		const char *title, const char *xlabel, const char *ylabel,
		const char* format,...){
    format2fn;
    (void)fig;
    loc_create_map_npad(loc,0);

    int npad=loc->map->npad;
    int nxm=loc->map->nx;
    int nx=loc->map->nx-npad*2;
    int ny=loc->map->ny-npad*2;
    double *opd0=calloc(nx*ny, sizeof(double));
    for(int iy=0; iy<ny; iy++){
	for(int ix=0; ix<nx; ix++){
	    int ii=loc->map->p[(ix+npad)+(iy+npad)*nxm];
	    if(ii){
		ii--;
		if(fabs(amp[ii])<EPS){
		    opd0[ix+iy*nx]=NAN;
		}else{
		    opd0[ix+iy*nx]=opd[ii]*amp[ii];
		}
	    }else{
		opd0[ix+iy*nx]=NAN;
	    }
	}
    }
    double limit[4];
    limit[0]=loc->map->ox+loc->dx*(npad-1);
    limit[1]=limit[0]+loc->dx*nx;
    limit[2]=loc->map->oy+loc->dx*(npad-1);
    limit[3]=limit[2]+loc->dx*ny;
    imagesc(fig, nx,ny,limit,zlim, opd0,  title, xlabel, ylabel,"%s",fn);
    free(opd0);
}

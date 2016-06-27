/*
  Copyright 2009-2016 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/socket.h>
#include <fcntl.h>
#include <errno.h>
#include <search.h>
#include "../math/mathdef.h"
#include "draw.h"
int DRAW_ID=0;
int DRAW_DIRECT=0;
int disable_draw=0; /*if 1, draw will be disabled  */
PNEW(lock);
#define MAXDRAW 1024
static int sock_helper=-1;
int listening=0;
int draw_single=0;//1: Only draw active frame. 0: draw all frames.
/*If not null, only draw those that match draw_fig and draw_fn*/
/**
   Contains functions for data visualization. 

   2013-02-20 

   A new scheme:
 
   The draw() routines pass sock to drawdaemon if need to launch it. 
   
   There are two scenarios where drawdaemon is used to display graphics.

   1) monitor (remote machine) opens a socket to connect to scheduler to request
   showing progress of a maos run. The scheduler passes the socket to maos. The
   monitor passes the other end of the socket to drawdaemon by fork()+exec(). In
   the end maos has direct connect to drawdaemon using sockets. draw() will
   write to the sockets the information to display.

   2) maos wants to launch drawdaemon. draw() will now talk to draw_helper()
   (running in a separate lightweight process) with pid through a socket pair
   (AF_UNIX type). The draw_helper() will then create another socket_pair, pass
   one end of it to draw(), and the other end to drawdaemon() by fork()+exec().


   todo: 
   1) code drawopdmap to call routines to copy data from gpu to cpu when needed to avoid unnecessary copy.
*/
/* List of list*/
typedef struct list_t{
    char *key;
    struct list_t *next;
    struct list_t *child;//child list
}list_t;
int sock_ndraw=0;//number of displays
typedef struct{
    int fd;
    int pause;
    list_t *list;
    char *figfn[2];
}sockinfo_t;
sockinfo_t sock_draws[MAXDRAW];

#define CATCH(A) if(A){\
	info("stwrite to %d failed with %s\n", \
	     sock_draw, strerror(errno));				\
	warning("\n\n\nwrite to sock_draw=%d failed\n\n\n",sock_draw);	\
	if(sock_helper<0&&!DRAW_DIRECT){				\
	    disable_draw=1;						\
	    warning("disable draw\n");					\
	}								\
	draw_remove(sock_draw,0);					\
	continue;							\
    }
#define STWRITESTR(A) CATCH(stwritestr(sock_draw,A))
#define STWRITEINT(A) CATCH(stwriteint(sock_draw,A))
#define STWRITECMDSTR(cmd,str) CATCH(stwriteint(sock_draw,cmd) || stwritestr(sock_draw,str))
#define STWRITE(p,len) CATCH(stwrite(sock_draw,p,len));

/**
   Listen to drawdaemon for update of fig, fn. The hold values are stored in figfn.
*/
static void listen_drawdaemon(sockinfo_t *sock_data){
    listening=1;
    int sock_draw=sock_data->fd;
    char **figfn=sock_data->figfn;
    //info2("draw is listening to drawdaemon at %d\n", sock_draw);
    int cmd;
    while(!streadint(sock_draw, &cmd)){
	switch(cmd){
	case DRAW_FIGFN:
	    {
		char *fig=0, *fn=0;
		streadstr(sock_draw, &fig);
		streadstr(sock_draw, &fn);
		if(figfn[0] && figfn[1] && (strcmp(figfn[0], fig) || strcmp(figfn[1], fn))){
		    info2("draw %d switch to fig=%s, fn=%s\n", sock_draw, fig, fn);
		}
		free(figfn[0]);
		free(figfn[1]);
		figfn[0]=fig;
		figfn[1]=fn;
	    }
	    break;
	case DRAW_PAUSE:
	    sock_data->pause=1;
	    info2("draw %d paused\n", sock_draw);
	    break;
	case DRAW_RESUME:
	    sock_data->pause=0;
	    info2("draw %d resumed\n", sock_draw);
	    break;
	default:
	    warning("cmd=%d is not understood\n", cmd);
	}
    }
    //info2("draw stop lisening to drawdaemon at %d\n", sock_draw);
    listening=0;
}

static int list_search(list_t **head, list_t **node, const char *key, int add){
    list_t *p=0;
    for(p=*head; p; p=p->next){
	if(!strcmp(p->key, key)){
	    break;
	}
    }
    int ans=p?1:0;
    if(add){
	if(!p){
	    p=mycalloc(1,list_t);
	    p->key=strdup(key);
	    p->next=*head;
	    *head=p;
	}
    }
    if(node) *node=p;
    return ans;
}
static void list_destroy(list_t **head){
    list_t *p;
    for(p=*head; p; p=*head){
	*head=p->next;
	if(p->child){
	    list_destroy(&p->child);
	}
	free(p);
    }
}
/*Add fd to list of drawing socks*/
int draw_add(int fd){
    if(fd==-1) return -1;
    if(sock_helper==-1){//externally added
	sock_helper=-2;
    }
    for(int ifd=0; ifd<sock_ndraw; ifd++){
	if(sock_draws[ifd].fd==fd){//already found
	    return 0;
	}
    }
    if(sock_ndraw<MAXDRAW){
	memset(&sock_draws[sock_ndraw], 0, sizeof(sockinfo_t));
	sock_draws[sock_ndraw].fd=fd;
	thread_new((thread_fun)listen_drawdaemon, &sock_draws[sock_ndraw]);
	sock_ndraw++;
	return 0;
    }else{
	return -1;
    }
}
static void draw_remove(int fd, int reuse){
    if(sock_ndraw<=0 || fd<0) return;
    int found=0;
    for(int ifd=0; ifd<sock_ndraw; ifd++){
	if(sock_draws[ifd].fd==fd){
	    found=1;
	    list_destroy(&sock_draws[ifd].list);
	    free(sock_draws[ifd].figfn[0]);
	    free(sock_draws[ifd].figfn[1]);
	}else if(found){//shift left
	    memcpy(&sock_draws[ifd-1], &sock_draws[ifd], sizeof(sockinfo_t));
	}
    }
    if(reuse){
	scheduler_send_socket(fd, DRAW_ID);
    }
    close(fd);
    if(found){
	sock_ndraw--;
    }else{
	warning("draw_remove: fd=%d is not found\n", fd);
    }
}
static int launch_drawdaemon(){
    int sv2[2];
    /*one end of sv2 will be passed back to draw_helper, the other end of sv2
      will be passed to drawdaemon.*/
    if(!socketpair(AF_UNIX, SOCK_STREAM, 0, sv2)){
	if(spawn_drawdaemon(sv2[1])){
	    warning("spawn drawdaemon failed\n");
	    close(sv2[0]);
	    sv2[0]=-1;
	}
	close(sv2[1]);
    }else{
	perror("socketpair");
	warning("socket pair failed, cannot launch drawdaemon\n");
	disable_draw=1;
	sv2[0]=-1;
    }
    return sv2[0];
}

/**
   A helper routine that forks in the early stage to launch drawdaemon.  We
   don't use the main routine to launch drawdaemon because it takes very long to
   fork a process that consumes Gig's of memory.

   We don't use drawdaemon to launch it because drawdaemon may not have the
   right DISPLAY related environment and also hard to work during ssh session
   even if environment is passed.
*/
void draw_helper(void){
    int sv[2];
    if(socketpair(AF_UNIX, SOCK_STREAM, 0, sv)){
	perror("socketpair");
	warning("socketpair failed, disable drawing.\n"); 
	disable_draw=1;
	sock_helper=-1;
    }
    pid_t pid=fork();
    if(pid<0){
	close(sv[0]); 
	close(sv[1]);
	warning("Fork failed. Return.\n");
	sock_helper=-1;
    }
    if(pid){/*Parent */
	sock_helper=sv[0];
	close(sv[1]);
    }else{/*Child */
	close(sv[0]);
	int cmd;
	while(!streadint(sv[1], &cmd)){
	    int sd=launch_drawdaemon();
	    if(stwritefd(sv[1], sd)){
		warning("write %d to sv[1] failed\n", sd);
	    }
	}
	close(sv[1]);
	_exit(0);//shall never return.
    }
}
/**
   Open a connection to drawdaemon. sock_draw may be set externally, in which case helper=-1.
*/
static int open_drawdaemon(){
    signal(SIGPIPE, SIG_IGN);
    if(disable_draw){
	return -1;
    }
    if(!sock_ndraw){
	if(DRAW_ID<=0){
	    DRAW_ID=getsid(0);
	    if(DRAW_ID<0){
		DRAW_ID=1;
	    }
	}
	int sock=-1;
	if(scheduler_recv_socket(&sock, DRAW_ID)){
	    sock=-1;
	}
	info("sock=%d, DRAW_ID=%d\n", sock, DRAW_ID);
	if(sock==-1){
	    if(DRAW_DIRECT){//directly fork and launch
		sock=launch_drawdaemon();
	    }else if(sock_helper!=-2){//use help to launch
		if(sock_helper==-1){
		    draw_helper();
		}
		if(stwriteint(sock_helper, DRAW_ID) || streadfd(sock_helper, &sock)){
		    sock=-1;
		    disable_draw=1;
		    close(sock_helper);
		    sock_helper=-1;
		    warning("Unable to talk to the helper to launch drawdaemon\n");
		}
	    }
	}
	if(sock!=-1){
	    draw_add(sock);
	}
    }
    return sock_ndraw==0?-1:0;
}

/* 
   Check whether we need to draw current page. True if
   1) not yet exist in the list. it is inserted into the list if add is valid and set as current page.
   2) is current page.
   3) when draw_single is not set.
   4) pause must not be set set.
*/
static int check_figfn(int ifd,  const char *fig, const char *fn, int add){
    if(disable_draw) return 0;
    list_t *child=0;
    list_search(&sock_draws[ifd].list, &child, fig, add);
    int found=0;
    if(child){
	found=list_search(&child->child, NULL, fn, add);
    }
    char **figfn=sock_draws[ifd].figfn;
    return !sock_draws[ifd].pause && (!draw_single || !found || (!mystrcmp(figfn[0], fig) && !mystrcmp(figfn[1], fn)));
}
/**
   Tell drawdaemon that this client will no long use the socket. Send the socket to scheduler for future reuse.
*/
void draw_final(int reuse){
    LOCK(lock);
    for(int ifd=0; ifd<sock_ndraw; ifd++){
	int sock_draw=sock_draws[ifd].fd;
	STWRITEINT(DRAW_FINAL);
	draw_remove(sock_draw, reuse);
    }
    UNLOCK(lock);
}

/**
   Check whether what we are drawing is current page.
*/
int draw_current(const char *fig, const char *fn){
    if(disable_draw) return 0;
    int current=0;
    for(int ifd=0; ifd<sock_ndraw; ifd++){
	/*Draw only if 1) first time (check with check_figfn), 2) is current active*/
	if(check_figfn(ifd, fig, fn, 0)){
	    current=1;
	}
    }
    return current;
}
/**
   Plot the coordinates ptsx, ptsy using style, and optionally plot ncir circles.
*/
void plot_points(const char *fig,    /**<Category of the figure*/
		 long ngroup,        /**<Number of groups to plot*/
		 loc_t **loc,        /**<Plot arrays of loc as grid*/
		 const dcell *dc,    /**<If loc isempty, use cell to plot curves*/
		 const int32_t *style,/**<Style of each point*/
		 const double *limit,/**<x min, xmax, ymin and ymax*/
		 const char *xylog,  /**<Whether use logscale for x, y*/
		 const dmat *cir,    /**<Data for the circles: x, y origin, radius, and color*/
		 const char *const* legend, /**<ngroup number of char**/
		 const char *title,  /**<title of the plot*/
		 const char *xlabel, /**<x axis label*/
		 const char *ylabel, /**<y axis label*/
		 const char *format, /**<subcategory of the plot.*/
		 ...){
    if(disable_draw) return;
    format2fn;
    LOCK(lock);
    if(open_drawdaemon()){/*failed to open. */
	warning("Failed to open\n");
	goto end;
    }
    for(int ifd=0; ifd<sock_ndraw; ifd++){
	/*Draw only if 1) first time (check with check_figfn), 2) is current active*/
	int sock_draw=sock_draws[ifd].fd;
	if(!check_figfn(ifd, fig, fn, 1)) continue;
	STWRITEINT(DRAW_START);
	if(loc){/*there are points to plot. */
	    for(int ig=0; ig<ngroup; ig++){
		STWRITEINT(DRAW_POINTS);
		STWRITEINT(loc[ig]->nloc);
		STWRITEINT(2);
		STWRITEINT(1);
		STWRITE(loc[ig]->locx, sizeof(double)*loc[ig]->nloc);
		STWRITE(loc[ig]->locy, sizeof(double)*loc[ig]->nloc);
	    }
	    if(dc){
		warning("both loc and dc are specified\n");
	    }
	}else if(dc){
	    if(ngroup!=dc->nx*dc->ny){
		warning("ngroup and dimension of dc mismatch\n");
		ngroup=dc->nx*dc->ny;
	    }
	    for(int ig=0; ig<ngroup; ig++){
		int nx=0, ny=0;
		double *p=NULL;
		if(dc->p[ig]){
		    nx=dc->p[ig]->nx;
		    ny=dc->p[ig]->ny;
		    p=dc->p[ig]->p;
		}
		STWRITEINT(DRAW_POINTS);
		STWRITEINT(nx);
		STWRITEINT(ny);
		STWRITEINT(0);
		if(p){
		    STWRITE(p, sizeof(double)*dc->p[ig]->nx*dc->p[ig]->ny);
		}
	    }
	}
	if(style){
	    STWRITEINT(DRAW_STYLE);
	    STWRITEINT(ngroup);
	    STWRITE(style,sizeof(uint32_t)*ngroup);
	}
	if(cir){
	    if(cir->nx!=4){
		error("Cir should have 4 rows\n");
	    }
	    STWRITEINT(DRAW_CIRCLE);
	    STWRITEINT(cir->ny);
	    STWRITE(cir->p, sizeof(double)*cir->nx*cir->ny);
	}
	if(limit){/*xmin,xmax,ymin,ymax */
	    STWRITEINT(DRAW_LIMIT);
	    STWRITE(limit, sizeof(double)*4);
	}
	if(xylog){
	    STWRITEINT(DRAW_XYLOG);
	    STWRITE(xylog, sizeof(const char)*2);
	}
	if(format){
	    STWRITEINT(DRAW_NAME);
	    STWRITESTR(fn);
	}
	if(legend){
	    STWRITEINT(DRAW_LEGEND);
	    for(int ig=0; ig<ngroup; ig++){
		STWRITESTR(legend[ig]);
	    }
	}
	STWRITECMDSTR(DRAW_FIG,fig);
	STWRITECMDSTR(DRAW_TITLE,title);
	STWRITECMDSTR(DRAW_XLABEL,xlabel);
	STWRITECMDSTR(DRAW_YLABEL,ylabel);
	STWRITEINT(DRAW_END);
    }
 end:
    UNLOCK(lock); 
}

/**
   Draw an image.
*/
//Data structure for imagesc
struct imagesc_t{
    char *fig; /**<Category of the figure*/
    long nx;   /**<the image is of size nx*ny*/
    long ny;   /**<the image is of size nx*ny*/
    double *limit; /**<x min; xmax; ymin and ymax*/
    double *zlim;/**< min;max; the data*/
    double *p; /**<The image*/
    char *title;  /**<title of the plot*/
    char *xlabel; /**<x axis label*/
    char *ylabel; /**<y axis label*/
    char *fn;
};
static void imagesc_do(struct imagesc_t *data){
    LOCK(lock);
    if(!open_drawdaemon()){
	char *fig=data->fig;
	long nx=data->nx;
	long ny=data->ny;
	double *limit=data->limit;
	double *zlim=data->zlim;
	double *p=data->p;
	char *title=data->title;
	char *xlabel=data->xlabel;
	char *ylabel=data->ylabel;
	char *fn=data->fn;
	for(int ifd=0; ifd<sock_ndraw; ifd++){
	    /*Draw only if 1) first time (check with check_figfn), 2) is current active*/
	    int sock_draw=sock_draws[ifd].fd;
	    if(!check_figfn(ifd, fig, fn, 1)) continue;
	    int32_t header[2];
	    header[0]=nx;
	    header[1]=ny;
	    STWRITEINT(DRAW_START);
	    STWRITEINT(DRAW_DATA);
	    STWRITE(header, sizeof(int32_t)*2);
	    STWRITE(p, sizeof(double)*nx*ny);
	    if(zlim){
		STWRITEINT(DRAW_ZLIM);
		STWRITE(zlim,sizeof(double)*2);
	    }
	    if(limit){/*xmin,xmax,ymin,ymax */
		STWRITEINT(DRAW_LIMIT);
		STWRITE(limit, sizeof(double)*4);
	    }
	    if(fn){
		STWRITECMDSTR(DRAW_NAME,fn);
	    }
	    STWRITECMDSTR(DRAW_FIG,fig);
	    STWRITECMDSTR(DRAW_TITLE,title);
	    STWRITECMDSTR(DRAW_XLABEL,xlabel);
	    STWRITECMDSTR(DRAW_YLABEL,ylabel);
	    STWRITEINT(DRAW_END);
	}
    }
    UNLOCK(lock);
    free(data->fig);
    free(data->limit);
    free(data->zlim);
    free(data->p);
    free(data->title);
    free(data->xlabel);
    free(data->ylabel);
    free(data->fn);
}
void imagesc(const char *fig, /**<Category of the figure*/
	     long nx,   /**<the image is of size nx*ny*/
	     long ny,   /**<the image is of size nx*ny*/
	     const double *limit, /**<x min, xmax, ymin and ymax*/
	     const double *zlim,/**< min,max, the data*/
	     const double *p, /**<The image*/
	     const char *title,  /**<title of the plot*/
	     const char *xlabel, /**<x axis label*/
	     const char *ylabel, /**<y axis label*/
	     const char *format, /**<subcategory of the plot.*/
	     ...){
    format2fn;
    if(disable_draw || !draw_current(fig, fn)) return;
    //Skip this drawing if line is busy.
    if((TRYLOCK(lock))){
	return;
    }else{
	UNLOCK(lock);
	//We copy all the data and put the imagesc job into a task
	//The task will free the data after it finishes.
	struct imagesc_t data;
	data.nx=nx;
	data.ny=ny;
#define datastrdup(x) data.x=(x)?strdup(x):0
#define datamemdup(x, size,type)			\
	if(x){						\
	    data.x=mymalloc(size,type);			\
	    memcpy(data.x, x,size*sizeof(type));	\
	}else{						\
	    data.x=0;					\
	}
	datastrdup(fig);
	datamemdup(limit, 4, double);
	datamemdup(zlim, 2, double);
	datamemdup(p, nx*ny, double);
	datastrdup(title);
	datastrdup(xlabel);
	datastrdup(ylabel);
	data.fn=format?strdup(fn):0;
#undef datastrdup
#undef datamemdup
	long group;
	QUEUE(group, imagesc_do, &data, 1, 0);
    }
}

/**
   Draw the OPD of real and imaginary of complex p defined on nx*ny grid. see imagesc()
*/
void imagesc_cmp_ri(const char *fig, long nx, long ny, const double *limit, const double *zlim,
		    const dcomplex *p, const char *title, const char *xlabel, const char *ylabel,
		    const char *format,...){
    format2fn;
    if(disable_draw || !draw_current(fig, fn)) return;

    double *pr,*pi;
    pr=mymalloc(nx*ny,double);
    pi=mymalloc(nx*ny,double);
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
void imagesc_cmp_ap(const char *fig, long nx, long ny, const double *limit,const double *zlim,
		    const dcomplex *p, const char *title, const char *xlabel, const char *ylabel,
		    const char *format,...){
    format2fn;
    if(disable_draw || !draw_current(fig, fn)) return;
    double *pr,*pi;
    int isreal=1;
    pr=mymalloc(nx*ny,double);
    pi=mycalloc(nx*ny,double);
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
void imagesc_cmp_abs(const char *fig, long nx, long ny, const double *limit,const double *zlim,
		     const dcomplex *p, const char *title, const char *xlabel, const char *ylabel,
		     const char *format,...){
    format2fn;
    if(disable_draw|| !draw_current(fig, fn)) return;
    double *pr;
    pr=mymalloc(nx*ny,double);
    for(int i=0; i<nx*ny; i++){
	pr[i]=cabs(p[i]);
    }
    imagesc(fig, nx, ny, limit, zlim, pr, title, xlabel, ylabel,
	    "%s abs",fn);
    free(pr);
}
#include "../math/mathdef.h"
#include "../math/mathdef.h"
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

void ddraw(const char *fig, const dmat *A, double *xylim, double *zlim,
	   const char *title, const char *xlabel, const char *ylabel,
	   const char *format,...){
    format2fn;
    imagesc(fig,A->nx,A->ny,xylim,zlim,A->p,title, xlabel, ylabel,"%s",fn);
}

/**
   Mapping the absolution value of complex array. see imagesc()
*/
void cdrawabs(const char *fig, const cmat *A, double *xylim, double *zlim,
	      const char *title, const char *xlabel, const char *ylabel,
	      const char *format,...){
  
    format2fn;
    imagesc_cmp_abs(fig,A->nx,A->ny,xylim,zlim,A->p,title,xlabel,ylabel,"%s abs",fn);
}

/**
   Mapping the real/imaginary part of complex array. see imagesc()
*/
void cdrawri(const char *fig, const cmat *A, double *xylim, double *zlim,
	     const char *title, const char *xlabel, const char *ylabel,
	     const char *format,...){
    
    format2fn;
    imagesc_cmp_ri(fig,A->nx,A->ny,xylim,zlim,A->p,title, xlabel, ylabel,"%s",fn);
}

/**
   Mapping the absolute and phase of complex array. see imagesc()
*/
void cdraw(const char *fig, const cmat *A, double *xylim, double *zlim,
	   const char *title, const char *xlabel, const char *ylabel,
	   const char *format,...){
    format2fn;
    imagesc_cmp_ap(fig,A->nx,A->ny,xylim,zlim,A->p,title, xlabel, ylabel,"%s",fn);
}

/**
   like ddraw, acting on map object. see imagesc()
*/
void drawmap(const char *fig, const map_t *map,  double *zlim,
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
void drawloc(const char *fig, loc_t *loc, double *zlim,
	     const char *title, const char *xlabel, const char *ylabel,
	     const char* format,...){
    format2fn;
    loc_create_map(loc);
    int npad=loc->npad;
    int nxm=loc->map->nx;
    int nx=loc->map->nx-npad*2;
    int ny=loc->map->ny-npad*2;
    double *opd0=mycalloc(nx*ny,double);
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
void drawopd(const char *fig, loc_t *loc, const double *opd,  double *zlim,
	     const char *title, const char *xlabel, const char *ylabel,
	     const char* format,...){

    format2fn;
    loc_create_map(loc);

    int npad=loc->npad;
    int nxm=loc->map->nx;
    int nx=loc->map->nx-npad*2;
    int ny=loc->map->ny-npad*2;
    dmat *opd0=dnew(nx,ny);
    for(int iy=0; iy<ny; iy++){
	for(int ix=0; ix<nx; ix++){
	    long ii=loc->map->p[(ix+npad)+(iy+npad)*nxm];
	    if(ii>0){
		opd0->p[ix+iy*nx]=opd[ii-1];
	    }else{
		opd0->p[ix+iy*nx]=NAN;
	    }
	}
    }
    double limit[4];
    limit[0]=loc->map->ox+fabs(loc->dx)*(npad-1/2);
    limit[1]=loc->map->ox+fabs(loc->dx)*(nx+npad-1/2);
    limit[2]=loc->map->oy+fabs(loc->dy)*(npad-1/2);
    limit[3]=loc->map->oy+fabs(loc->dy)*(ny+npad-1/2);
    imagesc(fig,nx,ny, limit,zlim,opd0->p,  title, xlabel, ylabel,"%s",fn);
    dfree(opd0);
}
/**
   Plot opd*amp with coordinate loc. see imagesc()
*/
void drawopdamp(const char *fig, loc_t *loc, const double *opd, const double *amp, double *zlim,
		const char *title, const char *xlabel, const char *ylabel,
		const char* format,...){
    format2fn;
    (void)fig;
    loc_create_map(loc);

    int npad=loc->npad;
    int nxm=loc->map->nx;
    int nx=loc->map->nx-npad*2;
    int ny=loc->map->ny-npad*2;
    double ampthres;
    dmaxmin(amp, loc->nloc, &ampthres, 0);
    ampthres*=0.5;
    double *opd0=mycalloc(nx*ny,double);
    for(int iy=0; iy<ny; iy++){
	for(int ix=0; ix<nx; ix++){
	    long ii=loc->map->p[(ix+npad)+(iy+npad)*nxm]-1;
	    if(ii>-1 && amp[ii]>ampthres){
		opd0[ix+iy*nx]=opd[ii];
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
